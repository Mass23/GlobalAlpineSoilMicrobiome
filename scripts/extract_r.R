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

# CHELSA: direct streaming from envicloud WSL bucket (monthly per-year + climatology bioclim SSP370 2011-2040)
library(xml2)

# Helper: list links in an HTML directory page and return absolute URLs
list_links <- function(url){
  res <- httr::GET(url)
  if(httr::status_code(res) != 200) return(character(0))
  doc <- tryCatch(xml2::read_html(httr::content(res, as='text', encoding='UTF-8')),
                  error = function(e) return(character(0)))
  nodes <- xml2::xml_find_all(doc, './/a')
  hrefs <- xml2::xml_attr(nodes, 'href')
  hrefs <- hrefs[!is.na(hrefs)]
  # make absolute
  hrefs <- sapply(hrefs, function(h){
    if(grepl('^https?://', h)) return(h)
    # handle relative links
    paste0(sub('/+$','', url), '/', sub('^/+', '', h))
  }, USE.NAMES = FALSE)
  hrefs
}

# Base URLs
monthly_base <- 'https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/monthly/'
bioclim_base <- 'https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/'

# Monthly variables list provided
monthly_vars <- c('clt','cmi','hurs','pet','47','pr','prec','ps','rsds','sfcWind','spei12','spi12','tas','tasmax','tasmin','tz','vpd')

# Helper to get available years for a given variable and month by parsing directory listings
get_available_years_for_var_month <- function(var, month){
  years <- integer(0)
  # list top-level months directory; expect year subdirs
  top_links <- list_links(monthly_base)
  # filter links that look like year directories, e.g., '1981/'
  year_dirs <- unique(gsub('.*/', '', top_links[grepl('/$', top_links)]))
  year_dirs <- year_dirs[grepl('^\\d{4}/?$', year_dirs)]
  year_dirs <- as.integer(gsub('/','',year_dirs))
  year_dirs <- sort(year_dirs)
  for(y in year_dirs){
    # list files in year dir
    year_url <- paste0(monthly_base, y, '/')
    files <- list_links(year_url)
    # match files containing var and month (e.g., 'pr_01' or 'tas_07')
    patt <- paste0(var, '_', sprintf('%02d', month))
    matches <- files[grepl(patt, files, ignore.case=TRUE) & grepl('\\.tif$', files, ignore.case=TRUE)]
    if(length(matches)>0) years <- c(years, y)
  }
  years
}

# Helper to download and sample one TIFF URL for all points
download_and_sample <- function(url){
  tmpf <- tempfile(fileext='.tif')
  message('Downloading ', url)
  res <- httr::GET(url, httr::write_disk(tmpf, overwrite=TRUE))
  httr::stop_for_status(res)
  r <- terra::rast(tmpf)
  pts_sp <- terra::vect(pts[,c('lon','lat')], geom=c('lon','lat'), crs='EPSG:4326')
  vals <- terra::extract(r, pts_sp)
  unlink(tmpf)
  vals[,2]
}

# For each monthly var, build a matrix values[point, year_index]
for(var in monthly_vars){
  message('Processing monthly variable: ', var)
  # for each point, get sample year and month
  sample_years <- as.integer(format(as.Date(pts$sample_date), '%Y'))
  sample_months <- as.integer(format(as.Date(pts$sample_date), '%m'))
  # Find available years for the variable/month by scanning top-level yearly dirs
  unique_months <- sort(unique(sample_months))
  available_years <- integer(0)
  for(m in unique_months){
    yrs <- get_available_years_for_var_month(var, m)
    available_years <- sort(unique(c(available_years, yrs)))
  }
  if(length(available_years)==0) stop('No available CHELSA monthly TIFFs found for variable ', var, ' on the envicloud host')

  # For each available year, download the file for each month present among points and sample
  year_vals <- list()
  year_list <- available_years
  for(y in year_list){
    # create storage for months in this year
    year_vals[[as.character(y)]] <- matrix(NA_real_, nrow=nrow(pts), ncol=length(unique_months))
    colnames(year_vals[[as.character(y)]]) <- as.character(unique_months)
    # list year dir files once
    year_url <- paste0(monthly_base, y, '/')
    files <- list_links(year_url)
    for(idx_m in seq_along(unique_months)){
      m <- unique_months[idx_m]
      patt <- paste0(var, '_', sprintf('%02d', m))
      match_files <- files[grepl(patt, files, ignore.case=TRUE) & grepl('\\.tif$', files, ignore.case=TRUE)]
      if(length(match_files)==0) next
      url <- match_files[1]
      vals <- tryCatch(download_and_sample(url), error = function(e) rep(NA_real_, nrow(pts)))
      year_vals[[as.character(y)]][, idx_m] <- vals
    }
  }

  # Build years vector and matrix
  years_vec <- as.integer(names(year_vals))
  if(length(years_vec)==0) stop('No sampled yearly values for var ', var)
  years_vec <- sort(years_vec)
  vals_mat <- matrix(NA_real_, nrow=nrow(pts), ncol=length(years_vec))
  colnames(vals_mat) <- as.character(years_vec)
  for(i in seq_along(years_vec)){
    y <- as.character(years_vec[i])
    for(pi in seq_len(nrow(pts))){
      m <- sample_months[pi]
      col_idx <- which(unique_months==m)
      if(length(col_idx)==0) next
      vals_mat[pi,i] <- year_vals[[y]][, col_idx]
    }
  }

  # For each point, decide value for its sample year using 10-year regression if needed
  result_vec <- numeric(nrow(pts))
  for(pi in seq_len(nrow(pts))){
    sy <- sample_years[pi]
    if(is.na(sy)){
      result_vec[pi] <- NA_real_; next
    }
    if(sy %in% years_vec){
      result_vec[pi] <- vals_mat[pi, which(years_vec==sy)]
    } else if(sy < min(years_vec)){
      result_vec[pi] <- vals_mat[pi,1]
    } else {
      recent_idx <- which(years_vec >= (max(years_vec)-9))
      x <- years_vec[recent_idx]
      yvals <- vals_mat[pi, recent_idx]
      valid <- !is.na(yvals)
      if(sum(valid) >= 2){
        # linear regression over last up to 10 years
        fit <- lm(yvals[valid] ~ x[valid])
        pred <- predict(fit, newdata=data.frame(x=sy))
        result_vec[pi] <- as.numeric(pred)
      } else {
        result_vec[pi] <- vals_mat[pi, length(years_vec)]
      }
    }
  }
  # attach result column named var_sampled
  colname <- paste0(var, '_sampled')
  out[[colname]] <- result_vec
}

# Bioclim SSP370 2011-2040: iterate GCM models and compute median across models
models <- c('GFDL-ESM4', 'IPSL-CM6A-LR', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'UKESM1-0-LL')
start_period <- '2011-2040'

# For each bio variable, collect per-model sampled values and compute median
for(b in 1:19){
  bioname <- paste0('bio', b)
  per_model_vals <- matrix(NA_real_, nrow=nrow(pts), ncol=length(models))
  colnames(per_model_vals) <- models
  for(mi in seq_along(models)){
    mdl <- models[mi]
    # construct expected URL pattern: bioclim_base/bioXX/2011-2040/{MODEL}/ssp370/FILE
    dir_url <- paste0(bioclim_base, bioname, '/', start_period, '/', mdl, '/ssp370/')
    files <- list_links(dir_url)
    tif_files <- files[grepl('\\.tif$', files, ignore.case=TRUE)]
    if(length(tif_files)==0){
      stop('No bioclim TIFFs found in ', dir_url, '; check path or model naming')
    }
    # pick first matching tif (should be one)
    f <- tif_files[1]
    per_model_vals[,mi] <- download_and_sample(f)
  }
  # compute median across models per point (na.rm=TRUE)
  med <- apply(per_model_vals, 1, function(x) median(x, na.rm=TRUE))
  out[[bioname]] <- med
  # also keep per-model columns
  for(mi in seq_along(models)) out[[paste0(bioname,'_',models[mi])]] <- per_model_vals[,mi]
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
