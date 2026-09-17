#!/usr/bin/env Rscript
# Merge chelsa_climate, soilgrids_sampled, and other_data by 'site' and write final CSV+Parquet.
req <- c('dplyr','arrow')
missing_pkgs <- req[!req %in% installed.packages()[,'Package']]
if(length(missing_pkgs)>0) stop('Missing R packages: ', paste(missing_pkgs, collapse=', '))

suppressPackageStartupMessages({
  library(dplyr)
  library(arrow)
})

chelsa_parq <- 'results/chelsa_climate.parquet'
soil_csv <- 'results/soilgrids_sampled.csv'
other_parq <- 'results/other_data.parquet'
final_csv <- 'results/sampled_soilgrids_r.csv'
final_parq <- 'results/sampled_soilgrids.parquet'

if(!file.exists(chelsa_parq)) stop('Missing ', chelsa_parq)
if(!file.exists(soil_csv)) stop('Missing ', soil_csv)
if(!file.exists(other_parq)) stop('Missing ', other_parq)

chelsa <- arrow::read_parquet(chelsa_parq)
soil <- read.csv(soil_csv, stringsAsFactors = FALSE)
other <- arrow::read_parquet(other_parq)

# left join in order: chelsa <- soil <- other by site, keep chelsa order
merged <- chelsa %>%
  left_join(soil, by = 'site') %>%
  left_join(other, by = 'site')

write.csv(merged, final_csv, row.names = FALSE)
arrow::write_parquet(merged, final_parq)
cat('Wrote', final_csv, 'and', final_parq, '\n')
