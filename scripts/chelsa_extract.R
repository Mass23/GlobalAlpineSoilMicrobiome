#!/usr/bin/env Rscript
# Deterministic CHELSA extractor: construct exact URLs (V.2.1) and sample rasters.
# Strict mode: fail if any expected URL is not readable.

req <- c('terra','sf','dplyr','arrow','httr','xml2')
missing_pkgs <- req[!req %in% installed.packages()[,'Package']]
if(length(missing_pkgs)>0) stop('Missing R packages: ', paste(missing_pkgs, collapse=', '))

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(arrow)
  library(httr)
  library(xml2)
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

# scaling and offsets will be initialized after monthly_vars is set (so monthly vars get the correct scale)
offsets_temp <- c(sprintf('bio%02d', c(10,11,1,5,6,8,9)), 'tas','tasmin','tasmax')

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

# list bucket prefix helper: returns folders and files for a given prefix
chelsa_list <- function(prefix){
  # prefix like 'chelsa/global/monthly/tasmax/'
  base <- 'https://os.unil.cloud.switch.ch/chelsa02/'
  resp <- try(httr::GET(base, query = list(prefix = prefix, delimiter = '/'), timeout(30)), silent = TRUE)
  if(inherits(resp, 'try-error')){
    return(list(error = paste('Failed to list bucket prefix', prefix, resp)))
  }
  if(httr::status_code(resp) != 200) return(list(error = paste('Non-200 from bucket list for', prefix, httr::status_code(resp))))
  xml <- try(xml2::read_xml(httr::content(resp, as = 'text', encoding = 'UTF-8')), silent = TRUE)
  if(inherits(xml, 'try-error')) return(list(error = paste('Failed to parse XML listing for', prefix)))
  ns <- xml_ns(xml)
  folders <- xml_text(xml_find_all(xml, './/d1:CommonPrefixes/d1:Prefix', ns))
  files <- xml_text(xml_find_all(xml, './/d1:Contents/d1:Key', ns))
  return(list(folders = folders, files = files))
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
  # HEAD-check to avoid HTML error pages; on failure provide detailed bucket listing and stop
  h <- try(httr::HEAD(https_url, httr::timeout(30)), silent = TRUE)
  if(inherits(h, 'try-error')){
    # attempt to discover available prefixes for helpful error
    prefix_dir <- NULL
    m <- regexec('/chelsa02/chelsa/global/monthly/([^/]+)/', https_url)
    mm <- regmatches(https_url, m)
    if(length(mm) && length(mm[[1]])>=2) prefix_dir <- paste0('chelsa/global/monthly/', mm[[1]][2], '/')
    if(is.null(prefix_dir)){
      stop(paste('HEAD failed for', https_url, h))
    } else {
      listing <- chelsa_list(prefix_dir)
      stop(sprintf('HEAD failed for %s; bucket listing for %s: folders=%s files=%s', https_url, prefix_dir,
                   paste(head(listing$folders,10), collapse=', '), paste(head(listing$files,10), collapse=', ')))
    }
  }
  if(httr::status_code(h) != 200){
    # provide helpful listing and then stop (strict error)
    prefix_dir <- NULL
    m <- regexec('/chelsa02/chelsa/global/monthly/([^/]+)/', https_url)
    mm <- regmatches(https_url, m)
    if(length(mm) && length(mm[[1]])>=2) prefix_dir <- paste0('chelsa/global/monthly/', mm[[1]][2], '/')
    if(!is.null(prefix_dir)){
      listing <- chelsa_list(prefix_dir)
      stop(sprintf('Non-200 (%s) for %s. Bucket listing for %s: folders=%s files=%s', httr::status_code(h), https_url, prefix_dir,
                   paste(head(listing$folders,10), collapse=', '), paste(head(listing$files,10), collapse=', ')))
    } else {
      stop(paste('Non-200 (', httr::status_code(h), ') for', https_url))
    }
  }
  ct <- tolower(httr::headers(h)[['content-type']])
  if(!is.null(ct) && grepl('html', ct)){
    # try listing and fail strictly
    prefix_dir <- NULL
    m <- regexec('/chelsa02/chelsa/global/monthly/([^/]+)/', https_url)
    mm <- regmatches(https_url, m)
    if(length(mm) && length(mm[[1]])>=2) prefix_dir <- paste0('chelsa/global/monthly/', mm[[1]][2], '/')
    listing <- if(!is.null(prefix_dir)) chelsa_list(prefix_dir) else list()
    stop(paste('URL returned HTML:', https_url, 'bucket listing (if available):', paste(head(listing$files,10), collapse=', ')))
  }
  r <- try(terra::rast(fetch_url), silent = TRUE)
  if(inherits(r, 'try-error')){
    stop(paste('terra::rast failed for', fetch_url, r))
  }
  v <- terra::vect(pts_df[,c('lon','lat')], geom = c('lon','lat'), crs = 'EPSG:4326')
  res <- terra::extract(r, v)
  return(res[,2])
}

out <- pts

# Use the exact monthly variables list provided (no discovery)
monthly_vars <- c('tasmax','tasmin','vpd','tas','spi12','spei12','sfcWind','rsds','pr','pet','hurs','cmi','clt')
cat('CHELSA: monthly variables to be sampled:', paste(monthly_vars, collapse=', '), '\n')

# now initialize scaling and offsets (monthly_vars must be present first)
vars_scales <- setNames(sapply(variables, function(v) if(grepl('^bio', v) || v %in% monthly_vars) 0.1 else 1), variables)
vars_offsets <- setNames(sapply(variables, function(v) if(v %in% offsets_temp) -273.15 else 0), variables)

# Helper to find exact S3 key for monthly var/year/month — fails if not exactly one match
find_monthly_key <- function(var, year, month){
  mm <- sprintf('%02d', as.integer(month))
  prefix <- paste0('chelsa/global/monthly/', var, '/', as.integer(year), '/')
  listing <- chelsa_list(prefix)
  if(!is.null(listing$error)) stop(listing$error)
  files <- listing$files
  expected <- paste0(prefix, sprintf('CHELSA_%s_%s_%d_V.2.1.tif', var, mm, as.integer(year)))
  exact_matches <- files[files == expected]
  if(length(exact_matches) == 1) return(paste0('https://os.unil.cloud.switch.ch/chelsa02/', exact_matches))
  if(length(exact_matches) > 1) stop(paste('Multiple exact matches for', expected, ':', paste(exact_matches, collapse=', ')))
  # try suffix match if listing contains longer keys
  suffix_pattern <- paste0('CHELSA_', var, '_', mm, '_', as.integer(year), '_V.2.1.tif$')
  partial_matches <- files[grepl(suffix_pattern, files)]
  if(length(partial_matches) == 1) return(paste0('https://os.unil.cloud.switch.ch/chelsa02/', partial_matches))
  if(length(partial_matches) > 1) stop(paste('Multiple partial matches for', expected, ':', paste(partial_matches, collapse=', ')))
  stop(paste('No exact key found for', expected, '\navailable (head):', paste(head(files,20), collapse='; ')))
}

# relaxed finder for historic trend sampling — returns NULL when missing/ambiguous (no stop)
find_monthly_key_relaxed <- function(var, year, month){
  mm <- sprintf('%02d', as.integer(month))
  prefix <- paste0('chelsa/global/monthly/', var, '/', as.integer(year), '/')
  listing <- chelsa_list(prefix)
  if(!is.null(listing$error)) return(NULL)
  files <- listing$files
  expected <- paste0(prefix, sprintf('CHELSA_%s_%s_%d_V.2.1.tif', var, mm, as.integer(year)))
  exact_matches <- files[files == expected]
  if(length(exact_matches) == 1) return(paste0('https://os.unil.cloud.switch.ch/chelsa02/', exact_matches))
  if(length(exact_matches) > 1) return(NULL)
  suffix_pattern <- paste0('CHELSA_', var, '_', mm, '_', as.integer(year), '_V.2.1.tif$')
  partial_matches <- files[grepl(suffix_pattern, files)]
  if(length(partial_matches) == 1) return(paste0('https://os.unil.cloud.switch.ch/chelsa02/', partial_matches))
  return(NULL)
}

# relaxed sampler for trend construction: returns NA vector on missing/failed access
sample_var_year_month_relaxed <- function(var, year, month){
  key <- find_monthly_key_relaxed(var, year, month)
  if(is.null(key)) return(rep(NA_real_, nrow(pts)))
  # try to read raster but don't stop on failure
  res <- try(sample_one_raster(key, pts), silent = TRUE)
  if(inherits(res, 'try-error')) return(rep(NA_real_, nrow(pts)))
  return(as.numeric(res))
}

# Helper to find bioclim historical key — strict
find_bioclim_key <- function(bio_num){
  bn <- sprintf('%02d', as.integer(bio_num))
  prefix <- paste0('chelsa/global/bioclim/bio', bn, '/1981-2010/')
  listing <- chelsa_list(prefix)
  if(!is.null(listing$error)) stop(listing$error)
  files <- listing$files
  expected <- paste0(prefix, sprintf('CHELSA_bio%s_1981-2010_V.2.1.tif', bn))
  exact_matches <- files[files == expected]
  if(length(exact_matches) == 1) return(paste0('https://os.unil.cloud.switch.ch/chelsa02/', exact_matches))
  if(length(exact_matches) > 1) stop(paste('Multiple exact matches for', expected, ':', paste(exact_matches, collapse=', ')))
  partial_matches <- files[grepl(paste0('CHELSA_bio', bn, '_1981-2010_V.2.1.tif$'), files)]
  if(length(partial_matches) == 1) return(paste0('https://os.unil.cloud.switch.ch/chelsa02/', partial_matches))
  if(length(partial_matches) > 1) stop(paste('Multiple partial matches for', expected, ':', paste(partial_matches, collapse=', ')))
  stop(paste('No exact key found for', expected, '\navailable (head):', paste(head(files,20), collapse='; ')))
}

# Monthly variables: strict resolution with preflight and extrapolation for years beyond available data
# Build available years map per var
available_years_map <- list()
for(var in monthly_vars){
  prefix <- paste0('chelsa/global/monthly/', var, '/')
  listing <- chelsa_list(prefix)
  if(!is.null(listing$error)) stop(listing$error)
  # listing$folders contain prefixes like 'chelsa/global/monthly/tasmax/1984/'
  yrs <- integer(0)
  if(length(listing$folders) > 0){
    yrs <- as.integer(gsub('.*/([0-9]{4})/', '\\1', listing$folders))
    yrs <- yrs[!is.na(yrs)]
  }
  # fallback: extract years from file keys if no folders parsed
  if(length(yrs) == 0 && length(listing$files) > 0){
    # file keys like 'chelsa/global/monthly/tasmax/1984/CHELSA_tasmax_01_1984_V.2.1.tif'
    yr_matches <- regmatches(listing$files, regexec('/([0-9]{4})/', listing$files))
    yrs2 <- unique(na.omit(sapply(yr_matches, function(x) if(length(x)>=2) as.integer(x[2]) else NA_integer_)))
    yrs <- sort(yrs2)
  }
  available_years_map[[var]] <- sort(unique(yrs))
}

# cache for resolved keys per (var,year,month)
key_cache <- new.env(parent = emptyenv())
get_monthly_key_cached <- function(var, year, month){
  k <- paste(var, year, sprintf('%02d', as.integer(month)), sep='|')
  if(exists(k, envir = key_cache)) return(get(k, envir = key_cache))
  # attempt to find exact key; if not found, this will stop (strict) in find_monthly_key
  url <- find_monthly_key(var, year, month)
  assign(k, url, envir = key_cache)
  return(url)
}

# function to sample a specific (var, year, month) for all points; returns raw values (no scaling)
sample_var_year_month <- function(var, year, month){
  url <- get_monthly_key_cached(var, year, month)
  sample_one_raster(url, pts)
}

for(var in monthly_vars){
  yrs_needed <- unique(as.integer(format(pts$sample_date, '%Y')))
  months_needed <- unique(as.integer(format(pts$sample_date, '%m')))
  avail <- available_years_map[[var]]
  if(length(avail) == 0) stop(paste('No available years found for monthly variable', var))
  max_avail <- max(avail)

  # prepare result vector (raw)
  vals_raw <- rep(NA_real_, nrow(pts))

  # first handle years that are available directly
  for(i in seq_len(nrow(pts))){
    yr <- as.integer(format(pts$sample_date[i], '%Y'))
    mm <- as.integer(format(pts$sample_date[i], '%m'))
    if(yr <= max_avail){
      # strict: must exist exactly
      url <- get_monthly_key_cached(var, yr, mm)
      vals_raw[i] <- sample_one_raster(url, pts[i, , drop = FALSE])
    }
  }

  # Now handle years beyond max_avail: extrapolate using last up to 10 years
  exceed_idx <- which(as.integer(format(pts$sample_date, '%Y')) > max_avail)
  if(length(exceed_idx) > 0){
    years_for_trend <- sort(avail)
    if(length(years_for_trend) < 2) stop(paste('Not enough years to extrapolate for', var))
    # sample trend years for the specific months present among exceed_idx — use relaxed sampler to take as much data as possible
    months_to_process <- unique(as.integer(format(pts$sample_date[exceed_idx], '%m')))
    # build matrix: rows=years_for_trend, cols=points
    trend_mat <- matrix(NA_real_, nrow = length(years_for_trend), ncol = nrow(pts))
    for(ri in seq_along(years_for_trend)){
      y <- years_for_trend[ri]
      for(mm in months_to_process){
        # sample that year's month with relaxed sampler (may return NAs)
        vals_y <- sample_var_year_month_relaxed(var, y, mm)
        cols <- which(as.integer(format(pts$sample_date, '%m')) == mm)
        trend_mat[ri, cols] <- vals_y[cols]
      }
    }
    # for each point that needs extrapolation, fit a linear trend on available years
    for(ii in exceed_idx){
      mm <- as.integer(format(pts$sample_date[ii], '%m'))
      # construct x (years) and y (values) from rows where trend_mat has non-NA for this column
      valid_rows <- which(!is.na(trend_mat[, ii]))
      if(length(valid_rows) < 2) stop(paste('Insufficient historical samples to extrapolate for site', pts$site[ii], 'variable', var))
      x <- years_for_trend[valid_rows]
      yvals <- trend_mat[valid_rows, ii]
      # fit simple linear regression (slope/intercept)
      m <- lm(yvals ~ x)
      target_year <- as.integer(format(pts$sample_date[ii], '%Y'))
      pred <- predict(m, newdata = data.frame(x = target_year))
      vals_raw[ii] <- as.numeric(pred)
    }
  }

  # apply scale & offset
  scale_val <- as.numeric(vars_scales[[var]])
  offset_val <- as.numeric(vars_offsets[[var]])
  out[[paste0(var, '_sampled')]] <- as.numeric(vals_raw) * scale_val + offset_val
}

# Bioclim climatologies: discover exact keys (strict) and sample
for(b in bios){
  bio_num <- as.integer(sub('bio', '', b))
  url <- find_bioclim_key(bio_num)
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
