#!/usr/bin/env Rscript
# Merge chelsa_climate, soilgrids_sampled, and other_data by 'site' and write final CSV+Parquet.
req <- c('dplyr','arrow')
missing_pkgs <- req[!req %in% installed.packages()[,'Package']]
if(length(missing_pkgs)>0) stop('Missing R packages: ', paste(missing_pkgs, collapse=', '))

suppressPackageStartupMessages({
  library(dplyr)
  library(arrow)
})

chelsa_csv <- 'data/chelsa_climate.csv'
soil_csv <- 'data/soilgrids_soil.csv'
other_csv <- 'data/other_data.csv'

final_csv <- 'results/all_downloaded_data.csv'

if(!file.exists(chelsa_csv)) stop('Missing ', chelsa_csv)
if(!file.exists(soil_csv)) stop('Missing ', soil_csv)
if(!file.exists(other_csv)) stop('Missing ', other_csv)

chelsa <- read.csv(chelsa_csv, stringsAsFactors = FALSE)
soil <- read.csv(soil_csv, stringsAsFactors = FALSE)
other <- read.csv(other_csv, stringsAsFactors = FALSE)

# left join in order: chelsa <- soil <- other by site, keep chelsa order
merged <- chelsa %>%
  left_join(soil, by = 'site') %>%
  left_join(other, by = 'site')

write.csv(merged, final_csv, row.names = FALSE)
cat('Wrote', final_csv, '\n')
