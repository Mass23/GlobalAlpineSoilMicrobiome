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

# CHELSA: prefer rchelsa for bioclim + monthly; fail loudly if monthly extraction not available
if(requireNamespace('rchelsa', quietly=TRUE)){
  message('Using rchelsa to extract CHELSA variables (bioclim + monthly)')
  exports <- getNamespaceExports('rchelsa')
  # bioclim
  if('chelsa_bioclim' %in% exports){
    res_bio <- tryCatch(rchelsa::chelsa_bioclim(points = pts[,c('lon','lat')]), error = function(e) stop('rchelsa::chelsa_bioclim failed: ', e$message))
    if(is.data.frame(res_bio)){
      for(nm in names(res_bio)) out[[nm]] <- res_bio[[nm]]
    } else stop('rchelsa::chelsa_bioclim returned unexpected result; aborting')
  } else stop('rchelsa installed but function chelsa_bioclim not found; abort')

  # monthly variables: try a set of likely function names; require one to exist
  monthly_fns <- c('chelsa_monthly','chelsa_get_monthly','chelsa_monthly_extract','chelsa_monthly_values','chelsa_monthly_ts')
  monthly_fn <- NULL
  for(fn in monthly_fns) if(fn %in% exports){ monthly_fn <- fn; break }
  if(is.null(monthly_fn)) stop('rchelsa installed but no known monthly extraction function found; aborting. Please check rchelsa documentation or update the package.')

  monthly_vars <- c('tas','tasmin','tasmax','prec')
  for(var in monthly_vars){
    message('Extracting monthly variable: ', var, ' using rchelsa::', monthly_fn)
    res_month <- tryCatch(do.call(getFromNamespace(monthly_fn, 'rchelsa'), list(points = pts[,c('lon','lat')], variable = var)),
                          error = function(e) stop('rchelsa monthly extract failed for ', var, ': ', e$message))
    # Expect res_month to be a data.frame or matrix with 12 columns (months)
    if(is.data.frame(res_month) || is.matrix(res_month)){
      # Ensure columns correspond to months; create column names var_01..var_12
      for(m in seq_len(ncol(res_month))){
        colname <- sprintf('%s_month%02d', var, m)
        out[[colname]] <- res_month[,m]
      }
    } else stop('rchelsa monthly extract returned unexpected result for ', var)
  }
} else {
  stop('rchelsa not installed: installer should have installed it. Aborting because monthly CHELSA variables are required.')
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
