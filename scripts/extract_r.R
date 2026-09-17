#!/usr/bin/env Rscript
# R extractor for SoilGrids (soilDB) + CHELSA (terra) + DEM (CopernicusDEM when available)
# Reads data/points.csv (columns: id,lat,lon, optional depth_cm)

# Required R packages (must be installed in the conda env before running)
req <- c("soilDB","dplyr","terra","httr","jsonlite","sf")
missing_pkgs <- req[!req %in% installed.packages()[,"Package"]]
if(length(missing_pkgs) > 0){
  message("Missing required R packages: ", paste(missing_pkgs, collapse=", "))
  message("Install them in the 'globalalpine-r' conda env before running. Example:")
  message("  conda activate globalalpine-r")
  message("  mamba install -n globalalpine-r -c conda-forge r-soildb r-terra r-sf r-httr r-jsonlite")
  message("If CopernicusDEM is desired and not available on conda, install from CRAN inside the env:")
  message("  R -e \"install.packages('CopernicusDEM', repos='https://cloud.r-project.org')\"")
  stop("Missing R packages: ", paste(missing_pkgs, collapse=", "))
}

suppressPackageStartupMessages({
  library(soilDB)
  library(dplyr)
  library(terra)
  library(httr)
  library(jsonlite)
  library(sf)
})

# optional CopernicusDEM
have_copdem <- requireNamespace("CopernicusDEM", quietly=TRUE)
if(have_copdem) library(CopernicusDEM)

# Read points
pts_file <- "data/points.csv"
if(!file.exists(pts_file)) stop("points.csv not found at data/points.csv")
pts <- read.csv(pts_file, stringsAsFactors=FALSE)
# Ensure columns
if(!all(c("lon","lat") %in% names(pts))) stop("points.csv must contain lon and lat columns")

# Soil variables and depth intervals to fetch
vars <- c("soc","phh2o")
depths <- c("0-5","5-15","15-30","30-60","60-100","100-200")

message("Querying SoilGrids for variables: ", paste(vars, collapse=","))
soil_res <- fetchSoilGrids(pts %>% mutate(id = ifelse(is.null(id), as.character(seq_len(n())), as.character(id))),
                           variables = vars,
                           depth_intervals = depths,
                           loc.names = c("id","lat","lon"),
                           verbose = TRUE)

# Replace ISRIC nodata sentinel with NA
soil_res[soil_res == -32768] <- NA

# Helper: compute 0-30cm thickness-weighted average using first three intervals
thicknesses <- c(5,10,15,30,40,100) # cm
compute_0_30 <- function(row, var){
  cols <- paste0(var, "_", depths)
  vals <- as.numeric(row[cols[1:3]])
  weights <- thicknesses[1:3] / 30
  if(all(is.na(vals))) return(NA_real_)
  sum(vals * weights, na.rm=TRUE) / sum(weights[!is.na(vals)])
}

# Linear interpolation helper using midpoints
midpoints <- c((0+5)/2, (5+15)/2, (15+30)/2, (30+60)/2, (60+100)/2, (100+200)/2)
linear_interp <- function(row, var, depth_cm){
  cols <- paste0(var, "_", depths)
  vals <- as.numeric(row[cols])
  if(is.na(depth_cm)) return(NA_real_)
  # if depth <=30, use first three midpoints
  if(depth_cm <= 30){
    use_idx <- which(!is.na(vals[1:3]))
    if(length(use_idx)==0) return(NA_real_)
    x <- midpoints[use_idx]
    y <- vals[use_idx]
    return(as.numeric(approx(x,y,xout=depth_cm,rule=2)$y))
  } else {
    # use all non-NA bands up to band containing depth
    valid_idx <- which(!is.na(vals))
    if(length(valid_idx)==0) return(NA_real_)
    # ensure we have bands up to depth
    band_idx <- which(cumsum(thicknesses) >= depth_cm)[1]
    if(is.na(band_idx)) band_idx <- length(depths)
    use_idx <- valid_idx[valid_idx <= band_idx]
    if(length(use_idx)==0) return(NA_real_)
    x <- midpoints[use_idx]
    y <- vals[use_idx]
    return(as.numeric(approx(x,y,xout=depth_cm,rule=2)$y))
  }
}

# combine results with input points order
out <- pts
for(v in vars){
  out[[paste0('soilgrid_nodepth_',v)]] <- apply(soil_res, 1, compute_0_30, var=v)
}
# if depth column exists, compute depth-interpolated values per rule
if('depth_cm' %in% names(pts)){
  out$depth_cm <- pts$depth_cm
  for(v in vars){
    out[[paste0('soilgrid_depth_',v)]] <- mapply(function(i, d) linear_interp(soil_res[i,], v, d), seq_len(nrow(soil_res)), out$depth_cm)
  }
} else {
  for(v in vars) out[[paste0('soilgrid_depth_',v)]] <- NA_real_
}

# Fetch CHELSA SCD raster (climatology) and sample using terra
chelsa_url <- 'https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_scd_1981-2010_V.2.1.tif'
message('Downloading CHELSA SCD to temp and sampling...')
chelsa_tmp <- tempfile(fileext='.tif')
GET(chelsa_url, write_disk(chelsa_tmp, overwrite=TRUE), progress())
chelsa_r <- rast(chelsa_tmp)
pts_sp <- vect(pts[,c('lon','lat')], geom=c('lon','lat'), crs='EPSG:4326')
chelsa_vals <- extract(chelsa_r, pts_sp)
out$chelsa_scd_1981_2010 <- chelsa_vals[,2]

# DEM: use CopernicusDEM if available, else skip
if(have_copdem){
  message('Using CopernicusDEM package to fetch DEM and derivatives')
  # CopernicusDEM provides access functions; attempt to get values per point
  # Here we use CopernicusDEM::get_elevation_point if available (fallback generic)
  if('get_elevation_point' %in% ls('package:CopernicusDEM')){
    elevs <- sapply(seq_len(nrow(pts)), function(i) CopernicusDEM::get_elevation_point(pts$lon[i], pts$lat[i]))
    out$dem_elevation <- elevs
    # For slope/aspect use terra::terrain if CopernicusDEM returns raster; skipped here for simplicity
  } else {
    message('CopernicusDEM installed but no get_elevation_point; skipping DEM derivs')
    out$dem_elevation <- NA_real_
  }
} else {
  message('CopernicusDEM not installed; skipping DEM derivs. To enable, run install.packages("CopernicusDEM") in R')
  out$dem_elevation <- NA_real_
}

# Print resulting table
print(out)

# Save to results
if(!dir.exists('results')) dir.create('results')
write.csv(out, 'results/sampled_soilgrids_r.csv', row.names=FALSE)
message('Wrote results/sampled_soilgrids_r.csv')
