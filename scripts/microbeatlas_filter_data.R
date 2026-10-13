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
# Builds the filtered SAMPLE table (coordinates, alpine flag, read counts, alpha metrics) and the
# list of sample accessions to keep. It does NOT write a BIOM any more: the BIOM is subset by the
# Snakemake rule `filter_biom` (biom subset-table, envs/biom_filter/environment.yml), using
#   <data_folder>/sample_accessions.txt     (one accession per line)
# Usage  : Rscript scripts/microbeatlas_filter_data.R <data_folder> [csv_out] [accessions_out]
#   data_folder : smk_config.yaml data_folder: the raw inputs (MicrobeAtlas files, shapefiles) AND the derived
#                 data live there
# Outputs: <data_folder>/sample_data_filtered.csv and <data_folder>/sample_accessions.txt, unless csv_out /
#          accessions_out are given.
# =============================================================================
message(">>> microbeatlas_filter_data.R v2026-10-13c (sample table + accession list; BIOM subset by Snakemake) <<<")
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
# FALSE (default): if the filtered CSV already exists, do not recompute anything: only (re)write the
#                  accession list from it.
# TRUE: ignore the existing CSV and redo everything.
FORCE_RERUN <- FALSE

# ---- One-time repair of an UNQUOTED csv ----
# An existing sample_data_filtered.csv written with quote = F cannot be read by pandas (prepare_points).
# Set TRUE and run the script ONCE by hand: the CSV is parsed with the brace-aware loader and rewritten
# with quoted fields (via a temp file). Then set it back to FALSE. It is not done automatically, so a
# Snakemake job can never modify a file that another job is reading.
REWRITE_QUOTED <- FALSE

# ---- Output paths ----
# data_folder (from smk_config.yaml) holds the raw inputs and receives the outputs, unless other paths are
# given as 2nd / 3rd argument.
suffix  <- if (TEST_MODE) "_TEST" else ""
add_suffix <- function(path) sub("(\\.[^.]+)$", paste0(suffix, "\\1"), path)
csv_out <- add_suffix(if (length(args) >= 2) args[2] else file.path(data_folder, "sample_data_filtered.csv"))
acc_out <- add_suffix(if (length(args) >= 3) args[3] else file.path(data_folder, "sample_accessions.txt"))
dir.create(dirname(csv_out), recursive = TRUE, showWarnings = FALSE)

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
#   1. CSV exists (and not FORCE_RERUN) -> only write the accession list from it
#   2. otherwise                        -> full pipeline: sample table + accession list
# =====================================================================================

if (!FORCE_RERUN && file.exists(csv_out)) {
  
  # ---------------------------------------------------------------- ACCESSION LIST ONLY
  message(csv_out, " already exists: only writing the accession list.")
  # The CSV may have been written with quote = F (commas inside text fields): the loader copes with both
  csv_tab <- read_filtered_table(csv_out)
  stopifnot("column 'Accessions' not found in the CSV" = "Accessions" %in% names(csv_tab))
  acc <- unique(csv_tab$Accessions)
  stopifnot("no accession could be read from the CSV" = length(acc) > 0)
  
  first_line <- readLines(csv_out, n = 2)[2]
  if (!startsWith(first_line, "\"")) {
    if (REWRITE_QUOTED) {
      tmp <- paste0(csv_out, ".tmp")
      write.csv(csv_tab, tmp, row.names = FALSE)
      stopifnot(file.rename(tmp, csv_out))
      message("Rewrote ", csv_out, " with quoted fields (", nrow(csv_tab), " rows): pandas can read it now")
    } else {
      warning(csv_out, " is NOT quoted: pandas (rule prepare_points) cannot read it. ",
              "Set REWRITE_QUOTED <- TRUE and run this script once by hand.")
    }
  }
  writeLines(acc, acc_out)                                  # written last: newer than the CSV
  message("Wrote ", acc_out, " (", length(acc), " accessions; ", nrow(csv_tab), " CSV rows)")
  
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
  
  # Sample table: quoted CSV (text fields can contain commas; scripts 1-2 read both quoted and unquoted)
  write.csv(sample_env_filtered_alpha, file = csv_out, row.names = F)
  # Accession list for the Snakemake rule filter_biom (biom subset-table -s): one accession per line
  writeLines(unique(sample_env_filtered_alpha$Accessions), acc_out)
  message("Wrote ", csv_out, " and ", acc_out, " (", nrow(sample_env_filtered_alpha), " samples)")
}