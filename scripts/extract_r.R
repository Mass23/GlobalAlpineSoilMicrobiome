#!/usr/bin/env Rscript
# Single-file R extractor: SoilGrids, CHELSA (rchelsa or direct), and DEM.
# No external config: all logic here. Fail loudly for debugging.

# Required packages
req <- c('soilDB','dplyr','terra','httr','jsonlite','sf','arrow')
missing_pkgs <- req[!req %in% installed.packages()[,'Package']]
if(length(missing_pkgs)>0){
  stop('Missing required R packages: ', paste(missing_pkgs, collapse=', '),
       '\nInstall via the conda env file envs/conda_env_r_geodata-download.yml and then run this script inside that env.')
}

suppressPackageStartupMessages({
  library(soilDB)
  library(dplyr)
  library(terra)
  library(httr)
  library(jsonlite)
  library(sf)
  library(arrow)
})

# Optional DEM package
have_copdem <- requireNamespace('CopernicusDEM', quietly=TRUE)
if(have_copdem) library(CopernicusDEM)

# Read input points
pts_file <- 'data/points.csv'
if(!file.exists(pts_file)) stop('data/points.csv not found')
pts <- read.csv(pts_file, stringsAsFactors=FALSE)
if(!all(c('lon','lat') %in% names(pts))) stop('points.csv must contain lon and lat columns')

# Output files (fixed)
csv_path <- 'results/sampled_soilgrids_r.csv'
parquet_path <- 'results/sampled_soilgrids.parquet'
if(!dir.exists(dirname(csv_path))) dir.create(dirname(csv_path), recursive=TRUE)
if(!dir.exists(dirname(parquet_path))) dir.create(dirname(parquet_path), recursive=TRUE)

# SoilGrids sampling
vars <- c('soc','phh2o')
depths <- c('0-5','5-15','15-30','30-60','60-100','100-200')
message('Fetching SoilGrids...')
soil_res <- fetchSoilGrids(dplyr::mutate(pts, id = ifelse(is.null(id), as.character(seq_len(n())), as.character(id))),
                           variables = vars, depth_intervals = depths,
                           loc.names = c('id','lat','lon'), verbose=TRUE)
soil_res[soil_res == -32768] <- NA

thicknesses <- c(5,10,15,30,40,100)
compute_0_30 <- function(row,var){
  cols <- paste0(var,'_',depths)
  vals <- as.numeric(row[cols[1:3]])
  weights <- thicknesses[1:3]/30
  if(all(is.na(vals))) return(NA_real_)
  sum(vals*weights, na.rm=TRUE)/sum(weights[!is.na(vals)])
}
midpoints <- c((0+5)/2,(5+15)/2,(15+30)/2,(30+60)/2,(60+100)/2,(100+200)/2)
linear_interp <- function(row,var,depth_cm){
  cols <- paste0(var,'_',depths)
  vals <- as.numeric(row[cols])
  if(is.na(depth_cm)) return(NA_real_)
  if(depth_cm <= 30){
    use_idx <- which(!is.na(vals[1:3])); if(length(use_idx)==0) return(NA_real_)
    x <- midpoints[use_idx]; y <- vals[use_idx]; return(as.numeric(approx(x,y,xout=depth_cm,rule=2)$y))
  } else {
    valid_idx <- which(!is.na(vals)); if(length(valid_idx)==0) return(NA_real_)
    band_idx <- which(cumsum(thicknesses) >= depth_cm)[1]; if(is.na(band_idx)) band_idx <- length(depths)
    use_idx <- valid_idx[valid_idx <= band_idx]; if(length(use_idx)==0) return(NA_real_)
    x <- midpoints[use_idx]; y <- vals[use_idx]; return(as.numeric(approx(x,y,xout=depth_cm,rule=2)$y))
  }
}

out <- pts
for(v in vars) out[[paste0('soilgrid_nodepth_',v)]] <- apply(soil_res,1,compute_0_30,var=v)
if('depth_cm' %in% names(pts)){
  out$depth_cm <- pts$depth_cm
  for(v in vars) out[[paste0('soilgrid_depth_',v)]] <- mapply(function(i,d) linear_interp(soil_res[i,],v,d), seq_len(nrow(soil_res)), out$depth_cm)
} else {
  for(v in vars) out[[paste0('soilgrid_depth_',v)]] <- NA_real_
}

# CHELSA: try rchelsa first
if(requireNamespace('rchelsa', quietly=TRUE)){
  message('Using rchelsa to extract CHELSA variables')
  if('chelsa_bioclim' %in% getNamespaceExports('rchelsa')){
    res_chelsa <- rchelsa::chelsa_bioclim(points=pts[,c('lon','lat')])
    if(is.data.frame(res_chelsa)) for(nm in names(res_chelsa)) out[[nm]] <- res_chelsa[[nm]]
    else stop('rchelsa returned unexpected result; aborting')
  } else stop('rchelsa installed but expected function chelsa_bioclim not found')
} else {
  message('rchelsa not available: downloading CHELSA bioclim TIFFs sequentially')
  base <- 'https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio'
  bios <- sprintf('CHELSA_bio%d_1981-2010_V.2.1.tif', 1:19)
  others <- c('CHELSA_scd_1981-2010_V.2.1.tif')
  all_files <- c(bios, others)
  pts_sp <- terra::vect(pts[,c('lon','lat')], geom=c('lon','lat'), crs='EPSG:4326')
  for(fn in all_files){
    url <- file.path(base, fn)
    message('Downloading ', url)
    tmpf <- tempfile(fileext='.tif')
    res <- httr::GET(url, httr::write_disk(tmpf, overwrite=TRUE), httr::progress())
    httr::stop_for_status(res)
    r <- terra::rast(tmpf)
    vals <- terra::extract(r, pts_sp)
    colname <- gsub('CHELSA_|_1981-2010_V.2.1.tif','', fn)
    out[[colname]] <- vals[,2]
    unlink(tmpf)
  }
}

# DEM
if(have_copdem){
  message('Using CopernicusDEM to fetch elevations')
  if('get_elevation_point' %in% ls('package:CopernicusDEM')){
    out$dem_elevation <- sapply(seq_len(nrow(pts)), function(i) CopernicusDEM::get_elevation_point(pts$lon[i], pts$lat[i]))
  } else stop('CopernicusDEM present but helper missing; abort')
} else {
  message('CopernicusDEM not installed; setting dem_elevation to NA')
  out$dem_elevation <- NA_real_
}

# Write outputs; fail if parquet cannot be written
message('Writing outputs...')
write.csv(out, csv_path, row.names=FALSE)
if(!is.character(parquet_path) || parquet_path=='') stop('Invalid parquet_path')
arrow::write_parquet(out, parquet_path)
message('Wrote ', csv_path, ' and ', parquet_path)
