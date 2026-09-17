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

variables <- c(sprintf('bio%02d', 1:19), 'fcf','fgd','scd',
               'clt','cmi','hurs','pet','pr','rsds','sfcWind','spei12','spi12','tas','tasmax','tasmin','vpd')
# candidate monthly variables (user-supplied) — exact list from your message
monthly_candidates <- c('clt','cmi','hurs','pet','pr','rsds','sfcWind','spei12','spi12','tas','tasmax','tasmin','vpd')

# metadata for monthly variables (short label, unit, description)
monthly_metadata <- list(
  clt = list(label='Total Cloud Cover Percentage', unit='percent'),
  cmi = list(label='Climate Moisture Index', unit='kg m-2 month-1'),
  hurs = list(label='Near-Surface Relative Humidity', unit='percent'),
  pet = list(label='Potential Evapotranspiration', unit='kg m-2 month-1'),
  pr = list(label='Precipitation', unit='kg m-2 month-1'),
  rsds = list(label='Surface Downwelling Shortwave Flux', unit='MJ m-2'),
  sfcWind = list(label='Near-Surface Wind Speed', unit='m s-1'),
  spei12 = list(label='Standardized Precipitation Evapotranspiration Index', unit='unitless'),
  spi12 = list(label='Standardized Precipitation Index', unit='unitless'),
  tas = list(label='Daily Mean Near-Surface Air Temperature', unit='K'),
  tasmax = list(label='Daily Max Near-Surface Air Temperature', unit='K'),
  tasmin = list(label='Daily Min Near-Surface Air Temperature', unit='K'),
  vpd = list(label='Vapor Pressure Deficit', unit='Pa')
)

# monthly_vars will be determined after helper functions are declared and will default to candidates if discovery fails
monthly_vars <- NULL


# bios variables are a subset of variables
bios <- sprintf('bio%02d', 1:19)
# models/institutions for projections (kept for future-proofing)
institutions <- c('UKESM1-0-LL','GFDL-ESM4','IPSL-CM6A-LR','MPI-ESM1-2-HR','MRI-ESM2-0')

# scaling and offsets (from user's previous code)
vars_scales <- setNames(sapply(variables, function(v) if(grepl('^bio', v) || v %in% monthly_vars) 0.1 else 1), variables)
offsets_temp <- c(sprintf('bio%02d', c(10,11,1,5,6,8,9)), 'tas','tasmin','tasmax')
vars_offsets <- setNames(sapply(variables, function(v) if(v %in% offsets_temp) -273.15 else 0), variables)

bucket_base <- 'https://os.unil.cloud.switch.ch/chelsa02/'

# helpers to construct /vsicurl/ paths
chelsa_monthly_url <- function(var, year, month){
  mm <- sprintf('%02d', as.integer(month))
  sprintf('https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/monthly/%s/%d/CHELSA_%s_%s_%d_V.2.1.tif',
          var, as.integer(year), var, mm, as.integer(year))
}

chelsa_bioclim_future_url <- function(bio, period, gcm, ssp){
  bio_num <- sprintf('%02d', as.integer(bio))
  gcm_lower <- tolower(gcm)
  sprintf('https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/bio%s/%s/%s/%s/CHELSA_%s_%s_bio%s_%s_V.2.1.tif',
          bio_num, period, gcm, ssp, gcm_lower, ssp, bio_num, period)
}

chelsa_bioclim_hist_url <- function(bio){
  bio_num <- sprintf('%02d', as.integer(bio))
  sprintf('https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/bio%s/1981-2010/CHELSA_bio%s_1981-2010_V.2.1.tif',
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

# discover available monthly vars by probing several representative years (avoid depending on sample_date range)
probe_years <- c(1979, 1984, 1990, 2000, 2010, 2018, 2020, 2022)
probe_month <- 1

included_monthly <- c()
for(var in monthly_candidates){
  ok <- FALSE
  for(yr in probe_years){
    test_url <- chelsa_monthly_url(var, yr, probe_month)
    https_url <- sub('^/vsicurl/', '', test_url)
    h <- try(httr::HEAD(https_url, httr::timeout(10)), silent = TRUE)
    if(!inherits(h, 'try-error') && httr::status_code(h) == 200){
      ct <- tolower(httr::headers(h)[['content-type']])
      if(is.null(ct) || !grepl('html', ct)){
        ok <- TRUE
        cat('CHELSA: monthly variable available (found year', yr, ') and will be used:', var, '\n')
        break
      }
    }
  }
  if(!ok){
    cat('CHELSA: monthly variable NOT available (no probe years matched), skipping:', var, '\n')
  } else {
    included_monthly <- c(included_monthly, var)
  }
}

monthly_vars <- included_monthly



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
  # apply scale & offset
  scale_val <- as.numeric(vars_scales[[var]])
  offset_val <- as.numeric(vars_offsets[[var]])
  out[[paste0(var, '_sampled')]] <- as.numeric(vals) * scale_val + offset_val
}

# Bioclim climatologies: use historical 1981-2010 climatologies (constructed path) and sample
for(b in bios){
  # b is like 'bio01'.. convert to number
  bio_num <- as.integer(sub('bio', '', b))
  url <- chelsa_bioclim_hist_url(bio_num)
  vals <- sample_one_raster(url, pts)
  # apply scale & offset
  scale_val <- as.numeric(vars_scales[[b]])
  offset_val <- as.numeric(vars_offsets[[b]])
  out[[b]] <- as.numeric(vals) * scale_val + offset_val
}

# Write outputs (CSV + Parquet)
write.csv(out, 'results/chelsa_climate.csv', row.names = FALSE)
arrow::write_parquet(out, 'results/chelsa_climate.parquet')
cat('Wrote results/chelsa_climate.csv and results/chelsa_climate.parquet\n')
