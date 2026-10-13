# If GDAL warns 'PROJ: proj_create_from_database ... share/proj failed', point PROJ at the
# folder containing proj.db BEFORE loading sf/terra (find it with: find $CONDA_PREFIX -name proj.db):
# Sys.setenv(PROJ_DATA = "/path/to/folder/with/proj.db", PROJ_LIB = "/path/to/folder/with/proj.db")

# The input BIOM is HDF5 and is only read: switch HDF5 file locking off, so a second R process (or a
# crashed one) holding a lock cannot make read_biom() fail with "unable to lock file, errno = 11".
# Set BEFORE the packages are loaded. Do NOT run two copies of this script at the same time.
Sys.setenv(HDF5_USE_FILE_LOCKING = "FALSE")

library(dplyr)
library(tidyr)
library(Matrix)
library(terra)
library(sf)
library(rbiom)
library(vegan)

# =============================================================================
# microbeatlas_filter_data.R
# Builds the filtered SAMPLE table (coordinates, alpine flag, read counts, alpha metrics), written as a TSV
# (the text fields contain commas but no tabs), and the list of sample accessions to keep. It does NOT write a BIOM any more: the BIOM is subset by the
# Snakemake rule `filter_biom` (biom subset-table, envs/biom_filter/environment.yml), using
#   <data_folder>/sample_accessions.txt     (one accession per line)
# Usage  : Rscript scripts/microbeatlas_filter_data.R <data_folder> [tsv_out] [accessions_out]
#   data_folder : smk_config.yaml data_folder: the raw inputs (MicrobeAtlas files, shapefiles) AND the derived
#                 data live there
# Outputs: <data_folder>/sample_data_filtered.tsv (tab-separated, no quoting) and
#          <data_folder>/sample_accessions.txt, unless tsv_out /
#          accessions_out are given.
# =============================================================================
message(">>> microbeatlas_filter_data.R v2026-10-13d (sample table as TSV + accession list; BIOM subset by Snakemake) <<<")
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
# FALSE (default): if the filtered TSV already exists, do not recompute anything: only (re)write the
#                  accession list from it.
# TRUE: ignore the existing TSV and redo everything.
FORCE_RERUN <- FALSE

# ---- Output paths ----
# data_folder (from smk_config.yaml) holds the raw inputs and receives the outputs, unless other paths are
# given as 2nd / 3rd argument.
suffix  <- if (TEST_MODE) "_TEST" else ""
add_suffix <- function(path) sub("(\\.[^.]+)$", paste0(suffix, "\\1"), path)
tsv_out <- add_suffix(if (length(args) >= 2) args[2] else file.path(data_folder, "sample_data_filtered.tsv"))
acc_out <- add_suffix(if (length(args) >= 3) args[3] else file.path(data_folder, "sample_accessions.txt"))
dir.create(dirname(tsv_out), recursive = TRUE, showWarnings = FALSE)

# ---- Helper functions ----

# The sample table is a TSV without quoting: its text fields contain commas (and braces) but no tabs.
# Tabs / line breaks inside a text field are replaced by a space when writing; the file is written to a temp
# file and renamed, so a crash cannot leave a truncated table. Always READ it with quote = "" (a stray double
# quote in a text field must not swallow the following lines).
write_tsv_safe <- function(df, path) {
  chr <- vapply(df, is.character, logical(1))
  df[chr] <- lapply(df[chr], function(x) gsub("[\t\r\n]+", " ", x))
  tmp <- paste0(path, ".tmp")
  write.table(df, file = tmp, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  stopifnot(file.rename(tmp, path))
}
read_tsv_table <- function(path)
  read.delim(path, quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE, fill = FALSE)

# BIOM counts as a sparse matrix (OTUs x samples), whatever rbiom stores them as
# (slam simple_triplet_matrix in older versions, Matrix sparse matrix in newer ones)
to_sparse <- function(m) {
  if (inherits(m, "simple_triplet_matrix"))
    Matrix::sparseMatrix(i = m$i, j = m$j, x = m$v, dims = c(m$nrow, m$ncol),
                         dimnames = list(m$dimnames[[1]], m$dimnames[[2]]))
  else as(m, "dgCMatrix")
}

# Alpha diversity, all metrics from vegan, from a SPARSE count matrix x (OTUs x samples, no "Unmapped"
# row). Only one batch of samples is densified at a time, after dropping the OTUs absent from the batch.
#   - alpha_total_n_reads : total reads (Unmapped excluded)
#   - alpha_shannon       : vegan::diversity(index = "shannon"), natural log
#   - alpha_n_otus        : vegan::specnumber() at full depth (depth-dependent!)
#   - alpha_pielou        : Pielou J = Shannon / ln(S), derived from the two vegan outputs
#                           above (vegan has no dedicated function); NA when S < 2
#   - alpha_n_otus_<d>k   : vegan::rarefy() = expected richness at depth d
#                           (Hurlbert; deterministic, no random subsampling);
#                           NA for samples with fewer than d reads
add_alpha_metrics_batched <- function(x, metadata,
                                      id_col = "Accessions",
                                      batch_size = 1000,
                                      depths = c(1000, 5000, 10000)) {
  
  # keep only valid shared samples
  ids      <- intersect(colnames(x), metadata[[id_col]])
  metadata <- metadata[match(ids, metadata[[id_col]]), ]
  batches  <- split(ids, ceiling(seq_along(ids) / batch_size))
  
  alpha <- lapply(batches, function(b) {
    xb <- x[, b, drop = FALSE]
    xb <- xb[Matrix::rowSums(xb) > 0, , drop = FALSE]                 # drop OTUs absent from the batch
    mat_t <- t(as.matrix(xb))                                         # samples x OTUs
    
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
# Two cases:
#   1. TSV exists (and not FORCE_RERUN) -> only write the accession list from it
#   2. otherwise                        -> full pipeline: sample table + accession list
# =====================================================================================

if (!FORCE_RERUN && file.exists(tsv_out)) {
  
  # ---------------------------------------------------------------- ACCESSION LIST ONLY
  message(tsv_out, " already exists: only writing the accession list.")
  tab <- read_tsv_table(tsv_out)
  stopifnot("column 'Accessions' not found in the table" = "Accessions" %in% names(tab))
  acc <- unique(tab$Accessions)
  stopifnot("no accession could be read from the table" = length(acc) > 0)
  writeLines(acc, acc_out)
  message("Wrote ", acc_out, " (", length(acc), " accessions; ", nrow(tab), " table rows)")
  
} else {
  
  # ---------------------------------------------------------------- FULL PIPELINE
  sample_env = read.csv(file.path(data_folder, 'microbeatlas/metadata_soilCurated_OTU97_nextMAPrelease.tsv'), sep='\t', header = T)
  sample_env = sample_env %>% mutate(Accessions = MAP_SID) %>%
    separate(MAP_SID, into = c("Reads_acc", "Sample_acc"), sep = "\\.")
  
  # Read the BIOM counts once, as a sparse matrix. NO rbiom subsetting (biom_file[...]) anywhere:
  # column selection on the sparse matrix is fast.
  message("Reading ", biom_path)
  biom_file <- read_biom(biom_path)
  x <- to_sparse(biom_file$counts)                                  # OTUs x samples
  rm(biom_file); invisible(gc())
  stopifnot("the BIOM has no 'Unmapped' row" = "Unmapped" %in% rownames(x))
  message(sprintf("BIOM: %d OTUs x %d samples", nrow(x), ncol(x)))
  
  if (TEST_MODE) {
    test_ids  <- head(colnames(x), N_TEST)
    x <- x[, test_ids, drop = FALSE]
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
  land <- terra::vect(file.path(data_folder, "ne_10m_land/ne_10m_land.shp"))
  pts <- terra::vect(sample_env_geo, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
  on_land <- terra::extract(land, pts)
  stopifnot(nrow(on_land) == nrow(pts))                      # one row per point, so the mask below is aligned
  sample_env_land <- sample_env_geo[on_land[,2] %in% "Land", ]   # sea points (NA) are dropped, not turned into NA rows
  message(nrow(sample_env_land), " samples on land")
  
  # Categorise Alpine samples
  # Based on this: https://figshare.com/articles/dataset/Global_distribution_and_bioclimatic_characterization_of_alpine_biomes/11710002?file=33157427
  v <- terra::vect(file.path(data_folder, "global_alpine_30m_v1_1/global_alpine_30m_v1_1.shp"))
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
  
  # keep the samples that are in the BIOM
  in_biom <- sample_env_alpine$Accessions %in% colnames(x)
  print(table(in_biom))
  sample_env_alpine <- sample_env_alpine[in_biom, ]
  
  # Read counts per sample from the sparse matrix (column selection only, no rbiom subsetting)
  xs <- x[, sample_env_alpine$Accessions, drop = FALSE]
  rm(x); invisible(gc())
  biom_counts          <- as.numeric(Matrix::colSums(xs))               # all reads
  biom_unmapped_counts <- as.numeric(xs["Unmapped", ])                  # reads in the "Unmapped" row
  sample_env_alpine$total_raw_counts      = biom_counts
  sample_env_alpine$unmapped_raw_counts   = biom_unmapped_counts
  sample_env_alpine$total_filtered_counts = biom_counts - biom_unmapped_counts   # reads without "Unmapped"
  
  sample_env_alpine = sample_env_alpine %>% filter(total_filtered_counts >= 1000)
  print(table(sample_env_alpine$Alpine))
  
  # sparse count matrix without the Unmapped row, for the samples that pass the depth filter
  xs <- xs[rownames(xs) != "Unmapped", sample_env_alpine$Accessions, drop = FALSE]
  sample_env_filtered_alpha = add_alpha_metrics_batched(xs, sample_env_alpine)
  rm(xs); invisible(gc())
  
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
  
  # Sample table: TSV, no quoting (text fields contain commas but no tabs)
  write_tsv_safe(sample_env_filtered_alpha, tsv_out)
  # Accession list for the Snakemake rule filter_biom (biom subset-table -s): one accession per line
  writeLines(unique(sample_env_filtered_alpha$Accessions), acc_out)
  message("Wrote ", tsv_out, " and ", acc_out, " (", nrow(sample_env_filtered_alpha), " samples)")
}
