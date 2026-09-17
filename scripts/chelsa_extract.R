#!/usr/bin/env Rscript
# Deterministic CHELSA extractor: construct exact URLs (V.2.1) and sample rasters.
# Strict mode: fail if any expected URL is not readable.

req <- c('terra','sf','dplyr','arrow')
missing_pkgs <- req[!req %in% installed.packages()[,'Package']]
if(length(missing_pkgs)>0) stop('Missing R packages: ', paste(missing_pkgs, collapse=', '))

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(arrow)
})

pts <- read.csv('data/points.csv', stringsAsFactors = FALSE)
if(!all(c('site','lon','lat','sample_date') %in% names(pts))) stop("data/points.csv must contain columns: site, lon, lat, sample_date")
pts$sample_date <- as.Date(pts$sample_date)

monthly_vars <- c('pet','pr','tas','tasmax','tasmin','rsds','sfcWind','hurs','clt','cmi','ps','spei12','spi12','prec','tz','vpd')
bios <- sprintf('bio%02d', 1:19)
models <- c('GFDL-ESM4','IPSL-CM6A-LR','MPI-ESM1-2-HR','MRI-ESM2-0','UKESM1-0-LL')

monthly_base <- 'os.unil.cloud.switch.ch/chelsa02/chelsa/global/monthly/'
bioclim_base <- 'os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/'

sample_one_raster <- function(url, pts_df){
  # Ensure the URL has a protocol; default to https if missing
  if(!grepl('^https?://', url)) url <- paste0('https://', url)
  r <- try(terra::rast(url), silent=TRUE)
  if(inherits(r,'try-error')) stop(paste('Failed to open raster URL:', url))
  v <- terra::vect(pts_df[,c('lon','lat')], geom = c('lon','lat'), crs = 'EPSG:4326')
  res <- terra::extract(r, v)
  return(res[,2])
}

out <- pts

# Monthly variables: sample file for point's sample_date year/month (strict - fail if missing)
for(var in monthly_vars){
  vals <- numeric(nrow(pts))
  for(i in seq_len(nrow(pts))){
    yr <- format(pts$sample_date[i], '%Y')
    mm <- format(pts$sample_date[i], '%m')
    url <- paste0(monthly_base, var, '/', yr, '/CHELSA_', var, '_', mm, '_', yr, '_V.2.1.tif')
    vals[i] <- sample_one_raster(url, pts[i, , drop = FALSE])
  }
  out[[paste0(var, '_sampled')]] <- vals
}

# Bioclim climatologies: sample each model file and compute median per point (strict; fail if file missing)
for(b in bios){
  mat <- matrix(NA_real_, nrow = nrow(pts), ncol = length(models))
  for(mi in seq_along(models)){
    model <- models[mi]
    file_name <- paste0('CHELSA_', tolower(model), '_ssp370_', b, '_2011-2040_V.2.1.tif')
    url <- paste0(bioclim_base, b, '/2011-2040/', model, '/ssp370/', file_name)
    mat[,mi] <- sample_one_raster(url, pts)
  }
  out[[b]] <- apply(mat, 1, function(x) median(x, na.rm = TRUE))
}

# Write outputs (CSV + Parquet)
write.csv(out, 'results/chelsa_climate.csv', row.names = FALSE)
arrow::write_parquet(out, 'results/chelsa_climate.parquet')
cat('Wrote results/chelsa_climate.csv and results/chelsa_climate.parquet\n')
