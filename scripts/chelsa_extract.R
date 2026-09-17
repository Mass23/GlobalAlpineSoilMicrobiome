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
})

pts <- read.csv('data/points.csv', stringsAsFactors = FALSE)
if(!all(c('site','lon','lat','sample_date') %in% names(pts))) stop("data/points.csv must contain columns: site, lon, lat, sample_date")
pts$sample_date <- as.Date(pts$sample_date)

monthly_vars <- c('pet','pr','tas','tasmax','tasmin','rsds','sfcWind','hurs','clt','cmi','ps','spei12','spi12','prec','tz','vpd')
bios <- sprintf('bio%02d', 1:19)
models <- c('GFDL-ESM4','IPSL-CM6A-LR','MPI-ESM1-2-HR','MRI-ESM2-0','UKESM1-0-LL')

bucket_base <- 'https://os.unil.cloud.switch.ch/chelsa02/'

# helper: list bucket prefix using S3 XML listing
list_chelsa <- function(prefix){
  # returns list(folders=..., files=...)
  url <- bucket_base
  resp <- try(httr::GET(url, query = list(prefix = prefix, delimiter = '/'), timeout(30)), silent = TRUE)
  if(inherits(resp, 'try-error')) stop(paste('Failed to list bucket prefix:', prefix, resp))
  if(httr::status_code(resp) != 200) stop(paste('Non-200 listing for prefix', prefix, 'status', httr::status_code(resp)))
  txt <- httr::content(resp, as='text', encoding='UTF-8')
  xml <- xml2::read_xml(txt)
  ns <- xml2::xml_ns(xml)
  prefixes <- xml2::xml_text(xml2::xml_find_all(xml, './/d1:CommonPrefixes/d1:Prefix', ns))
  keys <- xml2::xml_text(xml2::xml_find_all(xml, './/d1:Contents/d1:Key', ns))
  return(list(folders = prefixes, files = keys))
}

# construct /vsicurl/ path for a given key
vsicurl_for_key <- function(key){
  paste0('/vsicurl/', bucket_base, key)
}

sample_one_raster <- function(url, pts_df){
  # url is expected to be a /vsicurl/... path or full https
  # if it starts with /vsicurl/, leave as is; else ensure https
  if(grepl('^/vsicurl/', url)) {
    fetch_url <- url
  } else if(!grepl('^https?://', url)) {
    fetch_url <- paste0('/vsicurl/', bucket_base, url)
  } else {
    fetch_url <- paste0('/vsicurl/', url)
  }
  # Try to open with terra using vsicurl streaming
  r <- try(terra::rast(fetch_url), silent=TRUE)
  if(inherits(r,'try-error')) stop(paste('Failed to open raster URL with terra::rast:', fetch_url, '-', r))
  v <- terra::vect(pts_df[,c('lon','lat')], geom = c('lon','lat'), crs = 'EPSG:4326')
  res <- terra::extract(r, v)
  return(res[,2])
}

out <- pts

# Monthly variables: sample file for point's sample_date year/month (strict - discover via listing)
for(var in monthly_vars){
  vals <- numeric(nrow(pts))
  # cache per year-month
  cache_urls <- list()
  for(i in seq_len(nrow(pts))){
    yr <- format(pts$sample_date[i], '%Y')
    mm <- format(pts$sample_date[i], '%m')
    key_expected <- paste0('chelsa/global/monthly/', var, '/', yr, '/CHELSA_', var, '_', mm, '_', yr, '_V.2.1.tif')
    if(!is.null(cache_urls[[key_expected]])){
      url <- cache_urls[[key_expected]]
    } else {
      # list the prefix to find exact key
      lst <- list_chelsa(paste0('chelsa/global/monthly/', var, '/', yr, '/'))
      files <- lst$files
      match_key <- files[basename(files) == basename(key_expected)]
      if(length(match_key) == 0){
        stop(paste('Expected monthly file not found in bucket for', var, yr, mm, 'expected key', key_expected))
      }
      url <- vsicurl_for_key(match_key[1])
      cache_urls[[key_expected]] <- url
    }
    vals[i] <- sample_one_raster(url, pts[i, , drop = FALSE])
  }
  out[[paste0(var, '_sampled')]] <- vals
}

# Bioclim climatologies: discover available CHELSA_bio files and sample (prefer 1981-2010)
for(b in bios){
  # search for keys under chelsa/global/bioclim/ that contain the bio name
  lst_top <- list_chelsa('chelsa/global/bioclim/')
  all_files <- lst_top$files
  # if nothing at top, try to walk folders
  if(length(all_files) == 0){
    # drill into folders
    for(pref in lst_top$folders){
      sub <- list_chelsa(pref)
      all_files <- c(all_files, sub$files)
    }
  }
  pattern <- paste0('CHELSA_', b, '_')
  candidates <- all_files[grepl(pattern, basename(all_files), ignore.case = FALSE)]
  if(length(candidates) == 0){
    stop(paste('No bioclim files found for', b))
  }
  # prefer 1981-2010 if present
  prefer_idx <- grep('1981-2010', candidates)
  if(length(prefer_idx)>0) chosen <- candidates[prefer_idx[1]] else chosen <- candidates[1]
  url <- vsicurl_for_key(chosen)
  # sample the single-band climatology raster for all points
  vals <- sample_one_raster(url, pts)
  out[[b]] <- vals
}

# Write outputs (CSV + Parquet)
write.csv(out, 'results/chelsa_climate.csv', row.names = FALSE)
arrow::write_parquet(out, 'results/chelsa_climate.parquet')
cat('Wrote results/chelsa_climate.csv and results/chelsa_climate.parquet\n')
