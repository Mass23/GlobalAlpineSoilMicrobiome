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
# Builds the filtered SAMPLE table (coordinates, alpine flag, perialpine_alpine, BioProject, read counts, alpha metrics), written as a TSV
# (the text fields contain commas but no tabs), and the list of sample accessions to keep. It does NOT write a BIOM any more: the BIOM is subset by the
# Snakemake rule `filter_biom` (biom subset-table, envs/biom_filter/environment.yml), using
#   <data_folder>/sample_accessions.txt     (one accession per line)
# New columns (v2026-10-19a):
#   bioproject         BioProject accession of the sample's run (PRJNA..., PRJEB..., PRJDB...), from a metadata column if the metadata has
#                      one, else looked up in the ENA portal API from Reads_acc (cached in <data_folder>/run_to_bioproject.tsv; needs internet,
#                      once). Samples of one BioProject usually come from one study: they share protocols and are not independent.
#   perialpine_alpine  Lowland / Perialpine / Alpine: Perialpine = a LOWLAND sample (column Alpine) within PERIALPINE_KM (20) km of an
#                      alpine polygon. The column Alpine itself is unchanged.
#   cluster_id         one id per CLUSTER of samples: samples are linked if they share a BioProject or lie within CLUSTER_RADIUS_M (100) m of each other, and
#                      a cluster is a connected group (A and B share a project, B and C a place: A, B and C are one cluster). It is the single grouping
#                      column for the random effects of the analyses (samples of one cluster are not independent). A location shared by more than
#                      MAX_PROJECTS_PER_LOCATION BioProjects (a city / country centroid) does not link them, and a cluster holding more than
#                      MAX_CLUSTER_SHARE (5%) of all samples (a chain of links) is split back into BioProject clusters.
#   If the filtered TSV already exists (FORCE_RERUN = FALSE) and lacks these columns, only they are computed and added to it.
# Also dropped: samples whose GPT_Geotext contains one of EXCLUDE_GEOTEXT (agriculture / agricultural / crop; case-insensitive, anywhere
#   in the text), in the full pipeline and, if the TSV already exists, in that table too (its accession list is rewritten).
# Usage  : Rscript scripts/microbeatlas_filter_data.R <data_folder> [tsv_out] [accessions_out]
#   data_folder : smk_config.yaml data_folder: the raw inputs (MicrobeAtlas files, shapefiles) AND the derived
#                 data live there
# Outputs: <data_folder>/sample_data_filtered.tsv (tab-separated, no quoting) and
#          <data_folder>/sample_accessions.txt, unless tsv_out /
#          accessions_out are given.
# =============================================================================
message(">>> microbeatlas_filter_data.R v2026-10-19c (sample table as TSV + accession list; + bioproject, perialpine_alpine, cluster_id; GPT_Geotext filter) <<<")
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

# ---- Perialpine / BioProject ----
PERIALPINE_KM <- 20               # lowland samples within this distance (km) of an alpine polygon are "Perialpine" in perialpine_alpine
BIOPROJECT_BATCH   <- 80          # run accessions per ENA request
BIOPROJECT_WORKERS <- 4           # parallel ENA requests (ENA allows up to 50 per second; be polite)
RETRY_MISSING_BIOPROJECT <- FALSE # TRUE: ask ENA again for the runs it did not know the last time

# ---- Clusters (column cluster_id) ----
CLUSTER_RADIUS_M <- 100           # samples within this distance (metres) are at the same location
MAX_PROJECTS_PER_LOCATION <- 5    # a location with more BioProjects than this is a generic coordinate and does not link them (Inf: always link)
MAX_CLUSTER_SHARE <- 0.05         # a cluster with more than this share of ALL samples is split back into BioProject clusters (a chain of links)

# ---- Land-use words ----
# Samples whose GPT_Geotext contains one of these words (case-insensitive, anywhere in the text: "Crop" also drops "cropland", "crops")
# are dropped. Set to character(0) to switch the filter off.
EXCLUDE_GEOTEXT <- c("agriculture", "agricultural", "crop")

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

# ---- Drop samples by the words in their GPT_Geotext -----------------------------------------------------------------------------
drop_geotext <- function(tab, words = EXCLUDE_GEOTEXT) {
  if (!length(words)) return(tab)
  if (!"GPT_Geotext" %in% names(tab)) stop("column GPT_Geotext not found: it is needed to drop samples by their geographic text (EXCLUDE_GEOTEXT)")
  txt <- tolower(as.character(tab$GPT_Geotext)); txt[is.na(txt)] <- ""
  hits <- lapply(words, function(w) grepl(tolower(w), txt, fixed = TRUE))
  any_hit <- Reduce(`|`, hits)
  message(sprintf("GPT_Geotext filter: %d of %d samples dropped (%s; a sample can contain several words)", sum(any_hit), nrow(tab),
                  paste(sprintf("'%s': %d", words, vapply(hits, sum, numeric(1))), collapse = "; ")))
  if (any(any_hit)) message("  e.g. ", paste(substr(head(unique(trimws(as.character(tab$GPT_Geotext)[any_hit])), 4), 1, 60), collapse = " | "))
  tab[!any_hit, , drop = FALSE]
}

# ---- BioProject of every sample ----------------------------------------------------------------------------------------
# Taken from a metadata column that already holds BioProject accessions (PRJNA / PRJEB / PRJDB), if there is one. Otherwise it is
# looked up from the RUN accession (Reads_acc: SRR / ERR / DRR) in the ENA portal API, which mirrors the INSDC / NCBI SRA:
# study_accession = the BioProject (e.g. PRJNA123456), secondary_study_accession = the SRA study (SRP...). The answers are cached
# in <data_folder>/run_to_bioproject.tsv, so the lookup runs only once, an interrupted run resumes, and a rerun is instant.
# Runs that ENA does not know are cached as NA (RETRY_MISSING_BIOPROJECT = TRUE asks for them again).
find_project_col <- function(tab) {
  for (nm in names(tab)) {
    v <- as.character(tab[[nm]]); v <- head(v[!is.na(v) & nzchar(v)], 2000)
    if (length(v) && mean(grepl("^PRJ[EDN][A-Z]?[0-9]+$", v)) > 0.5) return(nm)
  }
  NULL
}

# one ENA request for a batch of run accessions: a data.frame (run_accession, study_accession, secondary_study_accession),
# or NULL if the service could not be reached (the batch is then NOT cached and is asked again next time)
ena_fetch <- function(runs, tries = 3) {
  options(timeout = max(300, getOption("timeout")))
  q <- paste0("run_accession=", runs, collapse = " OR ")
  u <- paste0("https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&format=tsv&limit=0",
              "&fields=run_accession,study_accession,secondary_study_accession&query=", utils::URLencode(q, reserved = TRUE))
  for (k in seq_len(tries)) {
    out <- tryCatch({
      txt <- readLines(url(u), warn = FALSE)
      if (length(txt) < 1) data.frame(run_accession = character(0), study_accession = character(0), secondary_study_accession = character(0))
      else utils::read.delim(text = txt, sep = "\t", quote = "", colClasses = "character", na.strings = "")
    }, error = function(e) NULL)
    if (!is.null(out)) return(out)
    Sys.sleep(5 * k)
  }
  NULL
}

bioproject_lookup <- function(runs, cache_file, batch = BIOPROJECT_BATCH, workers = BIOPROJECT_WORKERS, fetch = ena_fetch) {
  runs <- unique(runs[!is.na(runs) & nzchar(runs)])
  empty <- data.frame(run_accession = character(0), bioproject = character(0), sra_study = character(0))
  cache <- if (file.exists(cache_file)) read_tsv_table(cache_file) else empty
  if (RETRY_MISSING_BIOPROJECT) cache <- cache[!is.na(cache$bioproject), , drop = FALSE]
  todo <- setdiff(runs, cache$run_accession)
  message(sprintf("BioProject lookup: %d runs, %d already in %s, %d to ask ENA for", length(runs), length(runs) - length(todo),
                  basename(cache_file), length(todo)))
  if (length(todo)) {
    groups <- split(todo, ceiling(seq_along(todo) / batch))
    waves  <- split(seq_along(groups), ceiling(seq_along(groups) / (max(workers, 1) * 5)))
    n_fail <- 0L; n_done <- 0L; t0 <- Sys.time()
    for (w in waves) {
      one_batch <- function(g) {
        r <- fetch(g)
        if (is.null(r)) return(NULL)
        found <- data.frame(run_accession = r$run_accession, bioproject = r$study_accession, sra_study = r$secondary_study_accession)
        miss <- setdiff(g, r$run_accession)                         # asked for but not returned: cached as "not found"
        if (length(miss)) found <- rbind(found, data.frame(run_accession = miss, bioproject = NA_character_, sra_study = NA_character_))
        found
      }
      res <- if (workers > 1) parallel::mclapply(groups[w], one_batch, mc.cores = workers) else lapply(groups[w], one_batch)
      ok  <- vapply(res, function(x) is.data.frame(x), logical(1))
      n_fail <- n_fail + sum(!ok)
      if (any(ok)) {
        new <- do.call(rbind, res[ok]); new <- new[!new$run_accession %in% cache$run_accession, , drop = FALSE]
        write.table(new, cache_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", append = file.exists(cache_file),
                    col.names = !file.exists(cache_file))                  # appended after every wave: an interruption loses nothing
        cache <- rbind(cache, new)
      }
      n_done <- n_done + sum(lengths(groups[w]))
      message(sprintf("  %d / %d runs asked (%.0f min elapsed)", n_done, length(todo), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    }
    if (n_fail) warning(n_fail, " batch(es) could not be fetched from ENA (network?): rerun to ask for them again")
  }
  cache
}

add_bioproject <- function(tab) {
  if ("bioproject" %in% names(tab)) return(tab)
  col <- find_project_col(tab)
  if (!is.null(col)) {
    message("BioProject taken from the metadata column '", col, "'")
    tab$bioproject <- as.character(tab[[col]])
  } else {
    if (!"Reads_acc" %in% names(tab)) stop("the table has no Reads_acc (run accession) column to look the BioProject up from")
    lk <- bioproject_lookup(tab$Reads_acc, file.path(data_folder, "run_to_bioproject.tsv"))
    tab$bioproject <- lk$bioproject[match(tab$Reads_acc, lk$run_accession)]
  }
  tab$bioproject[!is.na(tab$bioproject) & !nzchar(tab$bioproject)] <- NA_character_
  n_ok <- sum(!is.na(tab$bioproject)); top <- sort(table(tab$bioproject), decreasing = TRUE)
  message(sprintf("bioproject: %d of %d samples have one, %d distinct projects; largest: %s", n_ok, nrow(tab), length(top),
                  paste(sprintf("%s (%d)", names(head(top, 6)), as.integer(head(top, 6))), collapse = ", ")))
  if (n_ok < 0.9 * nrow(tab)) warning("more than 10% of the samples have no BioProject (runs unknown to ENA, or the service was not reachable)")
  tab
}

# ---- Clusters: samples of one BioProject or at one location -----------------------------------------------------------------------
# Location clusters: samples within `radius_m` metres of an earlier one share its cluster (exact great-circle distances; 3D coordinates and a
# hash of cells of the radius size, so it also works for hundreds of thousands of samples). Returns an integer id per sample.
cluster_sites <- function(lat, lon, radius_m) {
  n <- length(lat); r <- max(radius_m, 0); cell <- max(r, 1); R <- 6371008.8
  la <- lat * pi / 180; lo <- lon * pi / 180
  X <- R * cos(la) * cos(lo); Y <- R * cos(la) * sin(lo); Z <- R * sin(la)
  cx <- floor(X / cell); cy <- floor(Y / cell); cz <- floor(Z / cell)
  off <- 2^17; key <- function(a, b, c) sprintf("%.0f", (a + off) + 2^19 * (b + off) + 2^38 * (c + off))
  nb <- as.matrix(expand.grid(-1:1, -1:1, -1:1))
  seen <- new.env(hash = TRUE, parent = emptyenv()); id <- integer(n); k <- 0L
  for (i in seq_len(n)) {
    js <- unlist(mget(key(cx[i] + nb[, 1], cy[i] + nb[, 2], cz[i] + nb[, 3]), envir = seen, ifnotfound = list(NULL)), use.names = FALSE)
    if (length(js)) {
      d2 <- (X[js] - X[i])^2 + (Y[js] - Y[i])^2 + (Z[js] - Z[i])^2; j <- which.min(d2)
      if (d2[j] <= r^2) { id[i] <- id[js[j]]; next }
    }
    k <- k + 1L; id[i] <- k; kk <- key(cx[i], cy[i], cz[i]); seen[[kk]] <- c(seen[[kk]], i)
  }
  id
}

# Connected groups: two samples are in the same cluster if they share a (non-NA) value in ANY of the label vectors, directly or through a
# chain of such samples (A and B share a BioProject, B and C a location: A, B and C are one cluster). Returns an integer id per sample.
connected_clusters <- function(labels) {
  codes <- lapply(labels, function(g) as.integer(factor(g)))             # NA stays NA: no link
  comp <- seq_along(codes[[1]])
  repeat {
    old <- comp
    for (g in codes) {
      ok <- !is.na(g)
      m <- tapply(comp[ok], g[ok], min)                                   # the smallest cluster id in every group
      comp[ok] <- as.integer(m[g[ok]])
    }
    if (identical(comp, old)) break
  }
  match(comp, unique(comp))
}

add_cluster_id <- function(tab, radius_m = CLUSTER_RADIUS_M, max_projects = MAX_PROJECTS_PER_LOCATION, max_share = MAX_CLUSTER_SHARE) {
  if ("cluster_id" %in% names(tab)) return(tab)
  ok  <- is.finite(tab$Longitude) & is.finite(tab$Latitude)
  loc <- rep(NA_integer_, nrow(tab)); loc[ok] <- cluster_sites(tab$Latitude[ok], tab$Longitude[ok], radius_m)
  bp  <- if ("bioproject" %in% names(tab)) as.character(tab$bioproject) else rep(NA_character_, nrow(tab))
  bp[!is.na(bp) & !nzchar(bp)] <- NA_character_
  if (all(is.na(bp))) warning("no BioProject known: the clusters are the locations only")
  # a location shared by many BioProjects is a generic coordinate (a city or country centroid, not a site): it must not glue unrelated
  # projects together, so it does not link them (samples of one project are still linked by the project)
  has <- !is.na(bp) & !is.na(loc)
  n_proj <- tapply(bp[has], loc[has], function(x) length(unique(x)))
  generic <- as.integer(names(n_proj)[n_proj > max_projects])
  loc_link <- loc; loc_link[loc %in% generic] <- NA_integer_
  cl <- connected_clusters(list(bp, loc_link))
  # a cluster that holds a large share of ALL samples is not a group of related samples but a chain of links that has glued unrelated
  # projects together: inside it only the BioProject links count (samples without a BioProject stay linked by location among themselves)
  big <- as.integer(names(which(table(cl) > max_share * nrow(tab))))
  if (length(big)) {
    sub <- which(cl %in% big)
    cl_sub <- connected_clusters(list(bp[sub], ifelse(is.na(bp[sub]), loc_link[sub], NA_integer_)))
    cl2 <- cl; cl2[sub] <- max(cl) + cl_sub; cl <- match(cl2, unique(cl2))
    message(sprintf("  %d cluster(s) held more than %.0f%% of the samples (%d samples): split back into BioProject clusters", length(big), 100 * max_share, length(sub)))
  }
  tab$cluster_id <- sprintf("C%06d", cl)
  sz <- table(tab$cluster_id)
  message(sprintf("cluster_id: %d clusters for %d samples (samples of one BioProject or within %g m of each other are linked); %d singletons, %d clusters of >= 5 samples; largest cluster %d samples (%.1f%%)",
                  length(sz), nrow(tab), radius_m, sum(sz == 1), sum(sz >= 5), max(sz), 100 * max(sz) / nrow(tab)))
  if (length(generic)) message(sprintf("  %d locations hold more than %d BioProjects (generic coordinates): they do not link projects", length(generic), max_projects))
  if (max(sz) > 0.25 * nrow(tab)) warning("one cluster holds more than 25% of the samples: unrelated projects may have been linked through shared coordinates (lower MAX_PROJECTS_PER_LOCATION)")
  tab
}

# ---- Perialpine: lowland samples within PERIALPINE_KM of an alpine polygon -------------------------------------------------------
# New column perialpine_alpine (factor: Lowland < Perialpine < Alpine); the column Alpine (Lowland / Alpine) is NOT changed. Distances
# are exact great-circle distances on the sphere (the s2 library: an index over the polygons, so 400k points take minutes, not hours).
# Points already flagged Alpine stay Alpine; the other points are "Perialpine" if an alpine polygon lies within PERIALPINE_KM.
add_perialpine <- function(tab, alpine_sf = NULL, km = PERIALPINE_KM) {
  if ("perialpine_alpine" %in% names(tab)) return(tab)
  if (is.null(alpine_sf)) alpine_sf <- sf::st_as_sf(terra::vect(alpine_shp))
  is_alp <- as.character(tab$Alpine) == "Alpine"
  idx <- which(!is_alp & is.finite(tab$Longitude) & is.finite(tab$Latitude))
  message(sprintf("perialpine_alpine: distance (s2, great-circle) of %d lowland samples to %d alpine polygons, threshold %g km", length(idx), nrow(alpine_sf), km))
  alp_s2 <- tryCatch(sf::st_as_s2(sf::st_geometry(alpine_sf), rebuild = TRUE),
                     error = function(e) stop("s2 could not use the alpine polygons (", conditionMessage(e), "): try sf::st_make_valid() on them first"))
  near <- logical(nrow(tab))
  for (ch in split(idx, ceiling(seq_along(idx) / 50000))) {
    pts <- s2::s2_geog_point(tab$Longitude[ch], tab$Latitude[ch])
    near[ch] <- lengths(s2::s2_dwithin_matrix(pts, alp_s2, distance = km * 1000)) > 0
    message(sprintf("  %d / %d points done", max(match(ch, idx)), length(idx)))
  }
  tab$perialpine_alpine <- factor(ifelse(is_alp, "Alpine", ifelse(near, "Perialpine", "Lowland")), levels = c("Lowland", "Perialpine", "Alpine"))
  print(table(tab$perialpine_alpine))
  tab
}

biom_path <- file.path(data_folder, "microbeatlas/otuTable_soilCurated_OTU97_nextMAPrelease.biom")
alpine_shp <- file.path(data_folder, "global_alpine_30m_v1_1/global_alpine_30m_v1_1.shp")

# =====================================================================================
# Two cases:
#   1. TSV exists (and not FORCE_RERUN) -> only write the accession list from it
#   2. otherwise                        -> full pipeline: sample table + accession list
# =====================================================================================

if (!FORCE_RERUN && file.exists(tsv_out)) {
  
  # ---------------------------------------------------------------- ACCESSION LIST ONLY
  message(tsv_out, " already exists: not recomputing it; only dropping samples by GPT_Geotext, adding missing columns (bioproject, cluster_id, perialpine_alpine) and writing the accession list.")
  tab <- read_tsv_table(tsv_out)
  stopifnot("column 'Accessions' not found in the table" = "Accessions" %in% names(tab))
  cols_before <- names(tab); n_before <- nrow(tab)
  tab <- drop_geotext(tab)
  tab <- add_bioproject(tab)
  tab <- add_cluster_id(tab)
  tab <- add_perialpine(tab)
  if (nrow(tab) != n_before || !identical(names(tab), cols_before)) {
    write_tsv_safe(tab, tsv_out)
    message("Updated ", tsv_out, ": ", n_before - nrow(tab), " samples dropped, columns added: ",
            if (length(setdiff(names(tab), cols_before))) paste(setdiff(names(tab), cols_before), collapse = ", ") else "none")
  }
  acc <- unique(tab$Accessions)
  stopifnot("no accession could be read from the table" = length(acc) > 0)
  writeLines(acc, acc_out)
  message("Wrote ", acc_out, " (", length(acc), " accessions; ", nrow(tab), " table rows)")
  
} else {
  
  # ---------------------------------------------------------------- FULL PIPELINE
  sample_env = read.csv(file.path(data_folder, 'microbeatlas/metadata_soilCurated_OTU97_nextMAPrelease.tsv'), sep='\t', header = T)
  sample_env = sample_env %>% mutate(Accessions = MAP_SID) %>%
    separate(MAP_SID, into = c("Reads_acc", "Sample_acc"), sep = "\\.")
  sample_env <- drop_geotext(sample_env)                    # no agriculture / crop samples
  
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
  v <- terra::vect(alpine_shp)
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
  rm(pts_sf)                                             # alpine_sf is kept for the perialpine step below
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
  
  # perialpine_alpine: Lowland / Perialpine (a lowland sample within PERIALPINE_KM km of an alpine polygon) / Alpine. The column Alpine
  # (Lowland / Alpine) is left as it is.
  sample_env_filtered_alpha <- add_perialpine(sample_env_filtered_alpha, alpine_sf)
  rm(alpine_sf)

  # BioProject of every sample (metadata column or ENA lookup of the run accession; cached)
  sample_env_filtered_alpha <- add_bioproject(sample_env_filtered_alpha)

  # one cluster id per group of samples that share a BioProject or a location (the grouping column of the random effects)
  sample_env_filtered_alpha <- add_cluster_id(sample_env_filtered_alpha)

  # Sample table: TSV, no quoting (text fields contain commas but no tabs)
  write_tsv_safe(sample_env_filtered_alpha, tsv_out)
  # Accession list for the Snakemake rule filter_biom (biom subset-table -s): one accession per line
  writeLines(unique(sample_env_filtered_alpha$Accessions), acc_out)
  message("Wrote ", tsv_out, " and ", acc_out, " (", nrow(sample_env_filtered_alpha), " samples)")
}