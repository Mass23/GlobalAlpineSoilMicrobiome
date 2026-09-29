library(dplyr)
library(ggplot2)
library(purrr)
library(maps)
library(tidyr)
library(raster)
library(terra)
library(data.table)
library(sf)
library(rbiom)
library(vegan)
library(phyloseq)
library(geosphere)

add_alpha_metrics_batched <- function(biom_file, metadata,
                                             id_col = "Accessions",
                                             batch_size = 1000) {
  biom_samples <- biom_file$samples
  meta_samples <- metadata[[id_col]]
  
  # keep only valid shared samples
  common <- intersect(biom_samples, meta_samples)
  
  biom_samples = biom_samples[biom_samples %in% common]
  metadata = metadata[match(biom_samples, metadata[[id_col]]), ]
  
  n = length(biom_samples)
  batches = split(seq_len(n), ceiling(seq_len(n) / batch_size))
  
  total_n_reads = numeric(n)
  shannon = numeric(n)
  n_otus_1k = rep(NA_real_, n)
  n_otus_5k = rep(NA_real_, n)
  n_otus_10k = rep(NA_real_, n)
  
  for (b in batches) {
    
    to_keep <- biom_samples[b]
    if (length(to_keep) == 0) next
    
    biom_b = biom_file[to_keep]
    biom_b_no_unmapped = rbiom::subset_taxa(biom_b, !grepl("Unmapped", biom_b$otus))
    
    mat = as.matrix(biom_b$counts)    
    mat = mat[rownames(mat) != 'Unmapped',]

    mat_t = t(mat)
    
    # Shannon
    shannon[b] = diversity(mat_t, index = "shannon")

    # total reads per sample (vector)
    total_n_reads[b] = rowSums(mat_t)
    
    # alpha diversity (includes shannon + richness at multiple depths)
    biom_1k = rbiom::rarefy(biom_b_no_unmapped, 1000)
    biom_5k = rbiom::rarefy(biom_b_no_unmapped, 5000)
    biom_10k = rbiom::rarefy(biom_b_no_unmapped, 10000)
    
    mat1k = as.matrix(biom_1k$counts)    
    mat1k_t = t(mat1k)
    
    mat5k = as.matrix(biom_5k$counts)    
    mat5k_t = t(mat5k)
    
    mat10k = as.matrix(biom_10k$counts)    
    mat10k_t = t(mat10k)
    
    # Observed OTUs (richness) at rarefaction depths
    n_otus_1k[b]  = specnumber(mat1k_t)
    n_otus_5k[b]  = specnumber(mat5k_t)
    n_otus_10k[b] = specnumber(mat10k_t)
    
  }
  
  metadata$alpha_total_n_reads = total_n_reads
  metadata$alpha_shannon = shannon
  metadata$alpha_n_otus_1k = n_otus_1k
  metadata$alpha_n_otus_5k = n_otus_5k
  metadata$alpha_n_otus_10k = n_otus_10k

  return(metadata)
}


setwd('~/Documents/MACE/AlpineSoilMicrobiome')

sample_env = read.csv('data/microbeatlas/metadata_soilCurated_OTU97_nextMAPrelease.tsv', sep='\t', header = T)
sample_env = sample_env %>% mutate(Accessions = MAP_SID) %>%
  separate(MAP_SID, into = c("Reads_acc", "Sample_acc"), sep = "\\.")

sample_env_filtered = sample_env %>% filter(is.numeric(Latitude) & is.numeric(Longitude))

# Filter only with date available
sample_env_date <- sample_env_filtered[sample_env_filtered$Date_collected != '',]
sample_env_geo <- sample_env_date[!is.na(sample_env_date$Date_collected),]
# 508'484

# Filter Latitude and Longitude that make sense
sample_env_geo$Latitude = as.numeric(sample_env_geo$Latitude)
sample_env_geo$Longitude = as.numeric(sample_env_geo$Longitude)

sample_env_geo = sample_env_geo %>% filter(Longitude > -180)
sample_env_geo = sample_env_geo %>% filter(Longitude < 180)
sample_env_geo = sample_env_geo %>% filter(Latitude > -90)
sample_env_geo = sample_env_geo %>% filter(Latitude < 90)
# 419'952

# Keep only samples that fall on Land
# Based on this: https://www.naturalearthdata.com/downloads/10m-physical-vectors/
land <- vect("data/ne_10m_land/ne_10m_land.shp")
pts <- vect(sample_env_geo, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
on_land <- extract(land, pts)
to_keep = na.omit(on_land[on_land[,2] == 'Land',1])
sample_env_land <- sample_env_geo[on_land[,2] == 'Land',]
# 419'952

# Categorise Alpine samples
# Based on this: https://figshare.com/articles/dataset/Global_distribution_and_bioclimatic_characterization_of_alpine_biomes/11710002?file=33157427
v <- vect("data/global_alpine_30m_v1_1/global_alpine_30m_v1_1.shp")
alpine_union <- aggregate(v)
pts <- vect(sample_env_land %>% dplyr::select("Longitude", "Latitude") %>% rename(lon=Longitude, lat=Latitude), crs = "EPSG:4326")
alpine <- extract(v, pts)
sample_env_land$Alpine <- ifelse(is.na(alpine[,4]), 0L, 1L)
sample_env_alpine = sample_env_land[!is.na(sample_env_land$Latitude),] # to be sure
sample_env_alpine = sample_env_land[!is.na(sample_env_land$Longitude),]
# 401'861

# How many alpine/non-alpine samples?
sample_env_alpine$Alpine = as.factor(sample_env_alpine$Alpine)  
table(sample_env_alpine$Alpine)
#      0      1 
# 396677   5184 --> Alpine samples!




biom_file <- read_biom("data/microbeatlas/otuTable_soilCurated_OTU97_nextMAPrelease.biom")
table(sample_env_alpine$Accessions %in% biom_file$samples)

biom_filtered <- biom_file[sample_env_alpine$Accessions]
table(sample_env_alpine$Accessions %in% biom_filtered$samples)
table(biom_filtered$samples %in% sample_env_alpine$Accessions)


# add unmapped stats to sample_env_filtered
library(slam)
biom_counts = col_sums(biom_filtered$counts)

biom_unmapped = rbiom::subset_taxa(biom_filtered, 'Unmapped')
biom_unmapped_counts = col_sums(biom_unmapped$counts)

biom_filtered = rbiom::subset_taxa(biom_filtered, !grepl("Unmapped", biom_filtered$otus))
biom_filtered_counts = col_sums(biom_filtered$counts)

sample_env_alpine$total_raw_counts = map_dbl(sample_env_alpine$Accessions, function(x) biom_counts[x])
sample_env_alpine$unmapped_raw_counts = map_dbl(sample_env_alpine$Accessions, function(x) biom_unmapped_counts[x])
sample_env_alpine$total_filtered_counts = map_dbl(sample_env_alpine$Accessions, function(x) biom_filtered_counts[x])

sample_env_alpine = sample_env_alpine %>% filter(total_filtered_counts >= 1000)
table(sample_env_alpine$Alpine)

sample_env_filtered_alpha = add_alpha_metrics_batched(biom_filtered, sample_env_alpine)

sample_env_filtered_alpha$Alpine = as.character(sample_env_filtered_alpha$Alpine)
sample_env_filtered_alpha$Alpine[sample_env_filtered_alpha$Alpine == "1"] = 'Alpine'
sample_env_filtered_alpha$Alpine[sample_env_filtered_alpha$Alpine == "0"] = 'Lowland'

sample_env_filtered_alpha$Alpine_regions = ''
sample_env_filtered_alpha$Alpine_regions[(sample_env_filtered_alpha$Longitude < -20)] = 'North America'
sample_env_filtered_alpha$Alpine_regions[(sample_env_filtered_alpha$Longitude < 0) & (sample_env_filtered_alpha$Latitude < 0)] = 'South America'
sample_env_filtered_alpha$Alpine_regions[(sample_env_filtered_alpha$Longitude > 0) & (sample_env_filtered_alpha$Latitude < 0)] = 'Oceania'
sample_env_filtered_alpha$Alpine_regions[(sample_env_filtered_alpha$Longitude > -20) & (sample_env_filtered_alpha$Longitude < 30) & (sample_env_filtered_alpha$Latitude > 0)] = 'Europe'
sample_env_filtered_alpha$Alpine_regions[(sample_env_filtered_alpha$Longitude > 30) & (sample_env_filtered_alpha$Latitude > 0)] = 'Asia'
sample_env_filtered_alpha$Alpine_regions[sample_env_filtered_alpha$Alpine == 'Lowland'] = 'Lowland'
table(sample_env_filtered_alpha$Alpine_regions, sample_env_filtered_alpha$Alpine)
table(sample_env_filtered_alpha$Alpine)

sample_env_filtered_alpha$Alpine_regions <- factor(sample_env_filtered_alpha$Alpine_regions, levels = c("Lowland", "Oceania", "Europe", "Asia", "North America", "South America"))
sample_env_filtered_alpha$Alpine <- factor(sample_env_filtered_alpha$Alpine, levels = c("Lowland", "Alpine"))

# Add new column - distance to nearest alpine region (the biome polygon), not to alpine sample points
#sample_pts <- vect(sample_env_filtered_alpha, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
#sample_env_filtered_alpha_dist <- sample_env_filtered_alpha
#sample_env_filtered_alpha_dist$dist_to_alpine <- as.numeric(terra::distance(sample_pts, alpine_union, unit = "km"))
#sample_env_filtered_alpha_dist$dist_to_alpine[sample_env_filtered_alpha_dist$Alpine == "Alpine"] <- 0

#sample_env_filtered_alpha_dist$Alpine <- factor(sample_env_filtered_alpha_dist$Alpine, levels = c("Lowland", "Perialpine", "Alpine"))
#sample_env_filtered_alpha_dist$Alpine[(sample_env_filtered_alpha_dist$Alpine == 'Lowland') &
#                                       (sample_env_filtered_alpha_dist$dist_to_alpine < 20)] = 'Perialpine'

biom_filtered = biom_filtered[sample_env_filtered_alpha$Accessions]

write.csv(sample_env_filtered_alpha, file = 'data/sample_data_filtered.csv', quote = F, row.names = F)
write.biom(biom_filtered, file = 'data/microbeatlas/samples-otus.97.mapped.metag.minfilter.refilt_filtered.biom')


# ---- Stratified random points (Lowland / Perialpine / Alpine), land only ----

sample_category_points <- function(category, n = 1000, alpine_vect, land_vect, peri_km = 20) {
  out <- data.frame(Longitude = numeric(0), Latitude = numeric(0), dist_to_alpine = numeric(0))

  while (nrow(out) < n) {
    lon <- runif(5000, -180, 180)
    lat <- asin(runif(5000, -1, 1)) * 180 / pi  # area-weighted, avoids pole oversampling

    candidates <- vect(data.frame(lon, lat), geom = c("lon", "lat"), crs = "EPSG:4326")

    on_land <- !is.na(extract(land_vect, candidates)[, 2])
    candidates <- candidates[on_land]; lon <- lon[on_land]; lat <- lat[on_land]
    if (length(candidates) == 0) next

    in_alpine <- !is.na(extract(alpine_vect, candidates)[, 2])
    d_km <- apply(distance(candidates, alpine_vect, unit = "km"), 1, min)

    keep <- switch(category,
      "Alpine"     = in_alpine,
      "Perialpine" = !in_alpine & d_km < peri_km,
      "Lowland"    = !in_alpine & d_km >= peri_km
    )

    out <- rbind(out, data.frame(Longitude = lon[keep], Latitude = lat[keep], dist_to_alpine = d_km[keep]))
  }

  out <- out[seq_len(n), ]
  out$Category <- category
  out
}

set.seed(123)

stratified_points <- bind_rows(
  sample_category_points("Alpine", alpine_vect = alpine_union, land_vect = land),
  sample_category_points("Perialpine", alpine_vect = alpine_union, land_vect = land),
  sample_category_points("Lowland", alpine_vect = alpine_union, land_vect = land)
)

write.csv(stratified_points, file = 'data/stratified_points.csv', quote = F, row.names = F)
