#!/usr/bin/env Rscript
# Deterministic CHELSA extractor: construct exact URLs (V.2.1) and sample rasters.
# Strict mode: fail if any expected URL is not readable.

req <- c('terra','sf','dplyr','arrow','httr')
missing_pkgs <- req[!req %in% installed.packages()[,'Package']]
if(length(missing_pkgs)>0) stop('Missing R packages: ', paste(missing_pkgs, collapse=', '))

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(arrow)
  library(httr)
})

pts <- read.csv('data/points.csv', stringsAsFactors = FALSE)
if(!all(c('site','lon','lat','sample_date') %in% names(pts))) stop("data/points.csv must contain columns: site, lon, lat, sample_date")
pts$sample_date <- as.Date(pts$sample_date)

monthly_vars <- c('pet','pr','tas','tasmax','tasmin','rsds','sfcWind','hurs','clt','cmi','ps','spei12','spi12','prec','tz','vpd')
bios <- sprintf('bio%02d', 1:19)
models <- c('GFDL-ESM4','IPSL-CM6A-LR','MPI-ESM1-2-HR','MRI-ESM2-0','UKESM1-0-LL')

bucket_base <- 'https://os.unil.cloud.switch.ch/chelsa02/'

# helpers to construct /vsicurl/ paths
chelsa_monthly_url <- function(var, year, month){
  mm <- sprintf('%02d', as.integer(month))
  sprintf('/vsicurl/https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/monthly/%s/%d/CHELSA_%s_%s_%d_V.2.1.tif',
          var, as.integer(year), var, mm, as.integer(year))
}

chelsa_bioclim_future_url <- function(bio, period, gcm, ssp){
  bio_num <- sprintf('%02d', as.integer(bio))
  gcm_lower <- tolower(gcm)
  sprintf('/vsicurl/https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/bio%s/%s/%s/%s/CHELSA_%s_%s_bio%s_%s_V.2.1.tif',
          bio_num, period, gcm, ssp, gcm_lower, ssp, bio_num, period)
}

chelsa_bioclim_hist_url <- function(bio){
  bio_num <- sprintf('%02d', as.integer(bio))
  sprintf('/vsicurl/https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/bio%s/1981-2010/CHELSA_bio%s_1981-2010_V.2.1.tif',
          bio_num, bio_num)
}

# sample helper: accept /vsicurl/... path or https; ensure /vsicurl/ for terra
sample_one_raster <- function(url, pts_df){
  if(grepl('^/vsicurl/', url)){
    fetch_url <- url
    # check underlying https exists
    https_url <- sub('^/vsicurl/', '', url)
  } else if(grepl('^https?://', url)){
    fetch_url <- paste0('/vsicurl/', url)
    https_url <- url
  } else {
    stop('sample_one_raster expects /vsicurl/ or https URL')
  }
  # HEAD-check to avoid HTML error pages
  h <- try(httr::HEAD(https_url, httr::timeout(30)), silent = TRUE)
  if(inherits(h, 'try-error')) stop(paste('HEAD failed for', https_url, h))
  if(httr::status_code(h) != 200) stop(paste('Non-200 for', https_url, httr::status_code(h)))
  ct <- tolower(httr::headers(h)[['content-type']])
  if(!is.null(ct) && grepl('html', ct)) stop(paste('URL returned HTML:', https_url))
  r <- try(terra::rast(fetch_url), silent = TRUE)
  if(inherits(r, 'try-error')) stop(paste('terra::rast failed for', fetch_url, r))
  v <- terra::vect(pts_df[,c('lon','lat')], geom = c('lon','lat'), crs = 'EPSG:4326')
  res <- terra::extract(r, v)
  return(res[,2])
}

out <- pts

# Monthly variables: construct exact /vsicurl/ URLs per variable/year/month and sample
for(var in monthly_vars){
  vals <- numeric(nrow(pts))
  for(i in seq_len(nrow(pts))){
    yr <- as.integer(format(pts$sample_date[i], '%Y'))
    mm <- as.integer(format(pts$sample_date[i], '%m'))
    url <- chelsa_monthly_url(var, yr, mm)
    # sample_one_raster accepts /vsicurl/ URL
    vals[i] <- sample_one_raster(url, pts[i, , drop = FALSE])
  }
  out[[paste0(var, '_sampled')]] <- vals
}

# Bioclim climatologies: use historical 1981-2010 climatologies (constructed path) and sample
for(b in bios){
  # b is like 'bio01'.. convert to number
  bio_num <- as.integer(sub('bio', '', b))
  url <- chelsa_bioclim_hist_url(bio_num)
  vals <- sample_one_raster(url, pts)
  out[[b]] <- vals
}

# Write outputs (CSV + Parquet)
write.csv(out, 'results/chelsa_climate.csv', row.names = FALSE)
arrow::write_parquet(out, 'results/chelsa_climate.parquet')
cat('Wrote results/chelsa_climate.csv and results/chelsa_climate.parquet\n')
