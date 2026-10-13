# If GDAL warns 'PROJ: proj_create_from_database ... share/proj failed', point PROJ at the
# folder containing proj.db BEFORE loading sf/terra (find it with: find $CONDA_PREFIX -name proj.db):
# Sys.setenv(PROJ_DATA = "/path/to/folder/with/proj.db", PROJ_LIB = "/path/to/folder/with/proj.db")

# The input BIOM is HDF5 and is only read: switch HDF5 file locking off, so a second R process (or a
# crashed one) holding a lock cannot make read_biom() fail with "unable to lock file, errno = 11".
# Set BEFORE the packages are loaded. Do NOT run two copies of this script at the same time.
Sys.setenv(HDF5_USE_FILE_LOCKING = "FALSE")

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
library(slam)

message(">>> microbeatlas_filter_data.R v2026-10-13b (resume: CSV read with the brace-aware loader) <<<")
args <- commandArgs(trailingOnly = TRUE)
# Snakemake passes the data folder as the first argument; fall back to the default path when run by hand
data_folder = if (length(args) >= 1) as.character(args[1]) else "/home/renku/work/epfl-altshuler/GlobalAlpineSoilMicrobiome/data"

# ---- Test switch ----
# TRUE: keep only the first N_TEST samples of the BIOM (and the matching metadata rows)
# and write outputs with a _TEST suffix so the real files are never overwritten.
# Note: the geographic / land / depth filters still apply afterwards, so fewer than
# N_TEST samples will remain, and the first samples may contain no alpine ones.
TEST_MODE <- FALSE   # set TRUE to test
N_TEST    <- 100

# ---- Resume switch ----
# FALSE (default): if the filtered CSV already exists, skip the metadata/land/alpine/alpha
#                  steps and only load + filter + write the BIOM (or do nothing if both exist).
# TRUE: ignore existing outputs and redo everything.
FORCE_RERUN <- FALSE

# ---- Output paths ----
suffix   <- if (TEST_MODE) "_TEST" else ""
csv_out  <- paste0('/home/renku/work/epfl-altshuler/GlobalAlpineSoilMicrobiome/data/sample_data_filtered', suffix, '.csv')
biom_out <- paste0('/home/renku/work/epfl-altshuler/GlobalAlpineSoilMicrobiome/data/microbeatlas/samples-otus.97.mapped.metag.minfilter.refilt_filtered', suffix, '.biom')

# ---- Helper functions ----

# Reads sample_data_filtered.csv even when it was written WITHOUT quotes (quote = F).
#  * commas inside {...} are protected: everything within braces is ONE field;
#  * rows that still have extra fields (unbraced commas in free text such as
#    "Washington state, USA") get the extras merged back into the geo_col column;
#  * rows that still do not fit (non-numeric coordinates, Alpine not Alpine/Lowland, too few
#    fields) are dropped and counted.
# A properly quoted CSV (first data line starts with a quote) is read with read.csv as usual.
read_filtered_table <- function(path, geo_col = "GPT_Geotext") {
  first <- readLines(path, n = 2)
  if (length(first) >= 2 && substr(first[2], 1, 1) == "\"")
    return(read.csv(path, stringsAsFactors = FALSE))
  
  lines <- readLines(path)
  hdr <- strsplit(lines[1], ",", fixed = TRUE)[[1]]; n <- length(hdr); lines <- lines[-1]
  
  m <- gregexpr("\\{[^{}]*\\}", lines)                       # protect commas inside {...}
  regmatches(lines, m) <- lapply(regmatches(lines, m), function(z) gsub(",", "\001", z, fixed = TRUE))
  parts <- strsplit(lines, ",", fixed = TRUE); rm(lines, m)
  
  g <- match(geo_col, hdr)
  extra_rows <- which(lengths(parts) > n)
  if (length(extra_rows)) {
    if (is.na(g)) stop("Rows with extra fields, but column '", geo_col, "' not found in the header")
    parts[extra_rows] <- lapply(parts[extra_rows], function(p) {
      e <- length(p) - n
      rest <- if (g + e < length(p)) p[(g + e + 1):length(p)] else character(0)
      c(p[seq_len(g - 1)], paste(p[g:(g + e)], collapse = ","), rest)
    })
  }
  ok_len <- lengths(parts) == n
  mat <- matrix(unlist(parts[ok_len], use.names = FALSE), ncol = n, byrow = TRUE); rm(parts)
  colnames(mat) <- hdr
  
  valid <- rep(TRUE, nrow(mat))
  if (all(c("Latitude", "Longitude", "Alpine") %in% hdr))
    valid <- is.finite(suppressWarnings(as.numeric(mat[, "Latitude"]))) &
    is.finite(suppressWarnings(as.numeric(mat[, "Longitude"]))) &
    mat[, "Alpine"] %in% c("Alpine", "Lowland")
  message(sprintf("%s: %d lines | %d with extra fields merged into %s | %d with too few fields | %d invalid rows dropped",
                  basename(path), length(ok_len), length(extra_rows), geo_col, sum(!ok_len), sum(!valid)))
  mat <- mat[valid, , drop = FALSE]
  
  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  df[] <- lapply(df, function(x) type.convert(gsub("\001", ",", x, fixed = TRUE), as.is = TRUE))
  df
}

# Write the BIOM to a temp file, then rename, so a crash mid-write can't leave a
# truncated .biom that would later be mistaken for a finished one.
write_biom_safe <- function(biom, path) {
  tmp <- paste0(path, ".tmp")
  write_biom(biom, file = tmp)
  stopifnot(file.rename(tmp, path))
}

# Column totals that work whatever rbiom stores counts as
# (slam simple_triplet_matrix in older versions, Matrix sparse matrix in newer ones).
col_totals <- function(m, rows = NULL) {
  if (!is.null(rows)) m <- m[rows, , drop = FALSE]
  if (inherits(m, "simple_triplet_matrix")) slam::col_sums(m) else Matrix::colSums(m)
}

# Alpha diversity, all metrics from vegan.
#   - alpha_total_n_reads : total reads (Unmapped excluded)
#   - alpha_shannon       : vegan::diversity(index = "shannon"), natural log
#   - alpha_n_otus        : vegan::specnumber() at full depth (depth-dependent!)
#   - alpha_pielou        : Pielou J = Shannon / ln(S), derived from the two vegan outputs
#                           above (vegan has no dedicated function); NA when S < 2
#   - alpha_n_otus_<d>k   : vegan::rarefy() = expected richness at depth d
#                           (Hurlbert; deterministic, no random subsampling);
#                           NA for samples with fewer than d reads
# Counts are densified one batch at a time only (lower batch_size if memory is tight).
add_alpha_metrics_batched <- function(biom_file, metadata,
                                      id_col = "Accessions",
                                      batch_size = 1000,
                                      depths = c(1000, 5000, 10000)) {
  
  # keep only valid shared samples, in BIOM order
  ids      <- intersect(biom_file$samples, metadata[[id_col]])
  metadata <- metadata[match(ids, metadata[[id_col]]), ]
  batches  <- split(ids, ceiling(seq_along(ids) / batch_size))
  
  alpha <- lapply(batches, function(b) {
    mat_t <- t(as.matrix(biom_file[b]$counts))                        # samples x OTUs
    mat_t <- mat_t[, colnames(mat_t) != "Unmapped", drop = FALSE]
    mat_t <- mat_t[, colSums(mat_t) > 0, drop = FALSE]                # drop OTUs absent from the batch
    
    out <- data.frame(
      id                  = rownames(mat_t),
      alpha_total_n_reads = rowSums(mat_t),
      alpha_shannon       = vegan::diversity(mat_t, index = "shannon"),
      alpha_n_otus        = vegan::specnumber(mat_t),
      row.names = NULL
    )
    out$alpha_pielou <- ifelse(out$alpha_n_otus > 1,
                               out$alpha_shannon / log(out$alpha_n_otus), NA_real_)
    
    # expected richness at fixed depths; only for samples that reach the depth
    for (d in depths) {
      ok <- out$alpha_total_n_reads >= d
      r  <- rep(NA_real_, nrow(out))
      if (any(ok)) r[ok] <- as.numeric(vegan::rarefy(mat_t[ok, , drop = FALSE], sample = d))
      out[[paste0("alpha_n_otus_", d / 1000, "k")]] <- r
    }
    out
  })
  
  alpha <- dplyr::bind_rows(alpha)
  names(alpha)[names(alpha) == "id"] <- id_col
  
  # join by accession (not by position)
  dplyr::left_join(metadata, alpha, by = id_col)
}

biom_path <- file.path(data_folder, "microbeatlas/otuTable_soilCurated_OTU97_nextMAPrelease.biom")

# =====================================================================================
# Three cases:
#   1. CSV and BIOM outputs both exist  -> nothing to do
#   2. CSV exists, BIOM output missing  -> load BIOM, filter to CSV accessions, write BIOM
#   3. otherwise (or FORCE_RERUN)       -> full pipeline
# =====================================================================================

if (!FORCE_RERUN && file.exists(csv_out) && file.exists(biom_out)) {
  
  message("Both ", csv_out, " and ", biom_out, " already exist - nothing to do. ",
          "(Set FORCE_RERUN <- TRUE to redo.)")
  
} else if (!FORCE_RERUN && file.exists(csv_out)) {
  
  # ---------------------------------------------------------------- FAST PATH
  message(csv_out, " already exists: skipping metadata/land/alpine/alpha steps, ",
          "only filtering the BIOM.")
  
  # The CSV was written with quote = F (commas inside text fields): read it with the same brace-aware
  # parser as scripts 1 and 2, so the filtered BIOM keeps exactly the samples those scripts will use.
  csv_tab <- read_filtered_table(csv_out)
  stopifnot("column 'Accessions' not found in the CSV" = "Accessions" %in% names(csv_tab))
  acc <- unique(csv_tab$Accessions)
  message(nrow(csv_tab), " CSV rows, ", length(acc), " unique accessions")
  stopifnot("no accession could be read from the CSV" = length(acc) > 0)
  rm(csv_tab)
  
  biom_file <- read_biom(biom_path)
  stopifnot("some accessions in the CSV are missing from the BIOM" = all(acc %in% biom_file$samples))
  
  biom_filtered <- biom_file[acc]
  rm(biom_file)
  # same order of operations as the full pipeline: subset samples, then drop the Unmapped row
  biom_filtered <- rbiom::subset_taxa(biom_filtered, !grepl("Unmapped", biom_filtered$otus))
  
  write_biom_safe(biom_filtered, biom_out)
  message("Wrote ", biom_out, " (", length(acc), " samples)")
  
} else {
  
  # ---------------------------------------------------------------- FULL PIPELINE
  sample_env = read.csv(file.path(data_folder, 'microbeatlas/metadata_soilCurated_OTU97_nextMAPrelease.tsv'), sep='\t', header = T)
  sample_env = sample_env %>% mutate(Accessions = MAP_SID) %>%
    separate(MAP_SID, into = c("Reads_acc", "Sample_acc"), sep = "\\.")
  
  # Read the BIOM once, up front (needed for the test subset)
  biom_file <- read_biom(biom_path)
  
  if (TEST_MODE) {
    test_ids  <- head(biom_file$samples, N_TEST)
    biom_file <- biom_file[test_ids]
    sample_env <- sample_env[sample_env$Accessions %in% test_ids, ]
    message("TEST_MODE: ", length(test_ids), " BIOM samples, ", nrow(sample_env),
            " of them found in the metadata")
    stopifnot(nrow(sample_env) > 0)
  }
  
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
  land <- vect(file.path(data_folder, "ne_10m_land/ne_10m_land.shp"))
  pts <- vect(sample_env_geo, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
  on_land <- extract(land, pts)
  stopifnot(nrow(on_land) == nrow(pts))                      # one row per point, so the mask below is aligned
  sample_env_land <- sample_env_geo[on_land[,2] %in% "Land", ]   # sea points (NA) are dropped, not turned into NA rows
  message(nrow(sample_env_land), " samples on land")
  
  # Categorise Alpine samples
  # Based on this: https://figshare.com/articles/dataset/Global_distribution_and_bioclimatic_characterization_of_alpine_biomes/11710002?file=33157427
  v <- vect(file.path(data_folder, "global_alpine_30m_v1_1/global_alpine_30m_v1_1.shp"))
  # No dissolve (aggregate) needed: sf builds a spatial index (STRtree) over the polygons,
  # so testing ~400k points against the undissolved layer is fast.
  sf::sf_use_s2(FALSE)                                   # planar GEOS predicates, as in terra
  alpine_sf <- sf::st_as_sf(v)
  # If st_intersects complains about invalid geometries, uncomment:
  # alpine_sf <- sf::st_make_valid(alpine_sf)
  # The layer declares geographic WGS84 (lon/lat); the sample coordinates are WGS84 lon/lat too.
  # Reuse the layer's own CRS object for the points instead of looking up EPSG:4326, so this
  # works even when PROJ's database is not found (see the PROJ_DATA note at the top).
  print(sf::st_crs(alpine_sf)$input)                    # what the shapefile declares
  stopifnot(!is.na(sf::st_crs(alpine_sf)), isTRUE(sf::st_is_longlat(alpine_sf)))
  pts_sf <- sf::st_as_sf(sample_env_land, coords = c("Longitude", "Latitude"),
                         crs = sf::st_crs(alpine_sf), remove = FALSE)
  # one logical per point (TRUE if it intersects any alpine polygon): aligned by construction
  sample_env_land$Alpine <- as.integer(lengths(sf::st_intersects(pts_sf, alpine_sf)) > 0)
  rm(pts_sf, alpine_sf)
  sample_env_alpine <- sample_env_land
  # 401'861
  
  # How many alpine/non-alpine samples?
  sample_env_alpine$Alpine = as.factor(sample_env_alpine$Alpine)
  print(table(sample_env_alpine$Alpine))
  #      0      1
  # 396677   5184 --> Alpine samples!
  
  print(table(sample_env_alpine$Accessions %in% biom_file$samples))
  
  biom_filtered <- biom_file[sample_env_alpine$Accessions]
  print(table(sample_env_alpine$Accessions %in% biom_filtered$samples))
  print(table(biom_filtered$samples %in% sample_env_alpine$Accessions))
  
  # add unmapped stats to sample_env_filtered
  cnt <- biom_filtered$counts
  is_unmapped <- rownames(cnt) == "Unmapped"
  stopifnot(any(is_unmapped))
  biom_counts          = col_totals(cnt)                  # all reads
  biom_unmapped_counts = col_totals(cnt, is_unmapped)     # reads in the "Unmapped" row
  rm(cnt)
  
  biom_filtered = rbiom::subset_taxa(biom_filtered, !grepl("Unmapped", biom_filtered$otus))
  biom_filtered_counts = col_totals(biom_filtered$counts) # reads without "Unmapped"
  
  # Vectorised lookups by name (replaces one map_dbl() call per sample).
  # Fail loudly if an accession is missing rather than silently getting NA.
  stopifnot(all(sample_env_alpine$Accessions %in% names(biom_counts)),
            all(sample_env_alpine$Accessions %in% names(biom_unmapped_counts)),
            all(sample_env_alpine$Accessions %in% names(biom_filtered_counts)))
  sample_env_alpine$total_raw_counts      = as.numeric(biom_counts[sample_env_alpine$Accessions])
  sample_env_alpine$unmapped_raw_counts   = as.numeric(biom_unmapped_counts[sample_env_alpine$Accessions])
  sample_env_alpine$total_filtered_counts = as.numeric(biom_filtered_counts[sample_env_alpine$Accessions])
  
  sample_env_alpine = sample_env_alpine %>% filter(total_filtered_counts >= 1000)
  print(table(sample_env_alpine$Alpine))
  
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
  print(table(sample_env_filtered_alpha$Alpine_regions, sample_env_filtered_alpha$Alpine))
  print(table(sample_env_filtered_alpha$Alpine))
  
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
  
  # CSV first, then BIOM: if the BIOM step crashes, the next run takes the fast path
  # (quote = F kept as in your version: scripts 1 and 2 read this file with read_filtered_table())
  write.csv(sample_env_filtered_alpha, file = csv_out, quote = F, row.names = F)
  write_biom_safe(biom_filtered, biom_out)
  message("Wrote ", csv_out, " and ", biom_out)
}