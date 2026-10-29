# =============================================================================
# Script 6 - Alpine specialist OTUs. Who are they (taxonomy), where are they (space), and how do they vary with the season and the soil depth (CLR abundance)?
#
# ALPINE SPECIALIST OTU = an OTU (not an aggregated taxon) that is (1) a niche SPECIALIST in the pooled alpine + lowland environment and (2) ALPINE: 
#   (1) niche breadth as in script 7: the abundance-weighted variance of the environment (z-scored principal components of climate, soil and terrain) of the units where the OTU lives, against a null model
#       that permutes the environment among the units: z <= -z_cut AND a breadth at least min_effect below the null mean;
#   (2) at least alpine_share_min of its abundance is in alpine units (the draws have as many lowland as alpine units, so 0.5 = no preference) and this share is above chance (z of the share against the
#       permutation of the group labels >= z_cut), and it is present in at least min_occ_alpine of the alpine units.
#   Units: random draws as in script 4 (one sample per UNIT, Alpine / matched Lowland, n_draws); a unit is a cluster, split by mountain range (alpine) or by 5-degree cell (lowland) when it spans several regions (split_units). The calls must hold in at least min_draw_frac of the draws. The mirror image (narrow and lowland) is
#   counted as "Lowland specialist" for contrast. OTUs rather than genera: they are the finest unit, but the most noisy: check the numbers of the log.
# FIGURE: A the taxonomic tree (domain > phylum > class > order > family) of the families with the most alpine specialist OTUs, with the number of specialist OTUs (and their share of the OTUs of the family;
#   * : the family is enriched in specialists, Fisher test, FDR); B the spatial repartition: ONE DOT PER MOUNTAIN RANGE (GMBA ranges as in script 3; continents without the GMBA file) coloured by the
#   estimated number of alpine specialist OTUs per sample, and how widespread the specialists are (the number of continents where they occur, against the other alpine OTUs); C the seasonal and the depth
#   variation of their CLR abundance (GAMs on the alpine samples: cyclic day of the year, and log soil depth), by phylum of the specialists and for all of them, and the estimated number of specialist OTUs
#   per sample of every mountain range.
# DEFINITION vs DISPLAY: the specialists are DEFINED on the draws (one sample per unit: a multi-region study contributes one sample per region); the shares of reads and the estimates of panel B / C are
#   SHOWN with all the samples (at most max_per_cluster per cluster and range in the models).
# GEOGRAPHY: every sample is placed (and assigned to its mountain range and continent) by its OWN coordinates. A cluster is a set of samples of one BioProject or within 100 m of each other and can span
#   several regions, so the mean position of a cluster is not a place and is never used. tables/alpine_cluster_spread.csv gives the spread of every cluster and
#   tables/alpine_samples_outside_ranges.csv lists the alpine samples that lie in no mountain range (to check their coordinates).
# ESTIMATE PER MOUNTAIN RANGE: the number of specialist OTUs DETECTED in a sample depends on how deeply it was sequenced and on how many reads did not map. A negative binomial GAM of the number of specialist
#   OTUs per sample (alpine samples of the ranges with at least min_range_clusters clusters, at most max_per_cluster per cluster, and lowland samples as a reference) on the mountain range, log10 sequencing
#   depth, the share of unmapped reads and a random cluster effect gives, for every range, the expected number of specialist OTUs of a typical sample at the SAME depth (the median of the alpine samples) and
#   the same share of unmapped reads (95% interval), next to the raw mean. It is a per-sample number: the total of a range grows with the number of samples.
# SOIL DEPTH: Depth_m is read as text: plain numbers, ranges ("0-10" = 5), values with units (cm, mm, m); values without a unit are metres or centimetres (depth_unit = "auto" decides from the values:
#   up to 3 = m). The axis is log10(depth in cm + 1), so that a depth of 0 is kept. The log says how many values were read and shows those that were not understood.
# CLR: the parts are the specialists of each of the n_parts richest phyla, the other specialists and all non-specialist reads; CLR = log2((count + pseudocount) / geometric mean of the parts). "All specialists" is
#   the 2-part version (specialists vs the rest, CLR = half the log2 ratio): its scale differs from the phylum parts, compare the shapes. Partial effects are centred.
# Output: fig6_alpine_specialist_otus.png and tables/range_specialist_estimates.csv, otu_niche_breadth.csv, alpine_specialist_otus.csv, family_enrichment.csv, cluster_specialist_share.csv, continent_occupancy.csv, specialist_clr_models.csv
# Input : sample_data_filtered.tsv (Date_collected, Depth_m for panel C), results/all_sampled_data.csv, the filtered BIOM, the OTU taxonomy file (needed)
# Usage : Rscript 6_alpine_specialist_otus.R [sample_table] [all_sampled_data.csv] [filtered_biom] [out_folder]
# Needs : R packages mgcv, rbiom, Matrix, dplyr, tidyr, ggplot2, maps
# NOT YET TESTED ON REAL DATA. Caveats: compositional data; niche breadth relative to the sampled environment; an OTU that is found in few clusters (one study) can look specialised: see the continents
# and the number of clusters; the CLR of the phylum parts depends on the parts chosen.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(Matrix)
  library(rbiom)
  library(mgcv)
  library(grid)
})

TEST_MODE <- FALSE   # TRUE: a tiny run
REUSE_FITS <- TRUE   # reuse the saved GMBA ranges if they match

args <- commandArgs(trailingOnly = TRUE)
data_folder     <- "/home/renku/work/epfl-altshuler/GlobalAlpineSoilMicrobiome/data"
results_default <- file.path(dirname(data_folder), "results")
sample_table  <- if (length(args) >= 1) as.character(args[1]) else file.path(data_folder, "sample_data_filtered.tsv")
env_table     <- if (length(args) >= 2) as.character(args[2]) else file.path(results_default, "all_sampled_data.csv")
filtered_biom <- if (length(args) >= 3) as.character(args[3]) else
  file.path(data_folder, "microbeatlas", "otuTable_soilCurated_OTU97_nextMAPrelease_filtered.biom")
out_folder    <- if (length(args) >= 4) as.character(args[4]) else file.path(results_default, "6_alpine_specialist_otus")

CFG <- list(
  sample_file = sample_table, env_file = env_table, biom_file = filtered_biom,
  id_table = "Sample_acc", id_env = "site", biom_id = "Accessions",
  group_col = "perialpine_alpine", alpine_value = "Alpine", lowland_value = "Lowland", cluster_col = "cluster_id",
  rare_depth = 1000, lowland_reps = 20, n_draws = 10, nperm = 99,
  env_sets = list(
    Climate = c("^bio[0-9]{2}$", "^(fcf|fgd|scd)$", "^(clt|cmi|hurs|pet|pr|rsds|sfcWind|spei12|spi12|tas|tasmax|tasmin|vpd)$"),
    Soil    = c("_30cm_avg$"),
    Terrain = c("^elevation_m$", "^slope_deg$", "^aspect_northness$", "^aspect_eastness$", "^tpi_m$", "^roughness_m$")),
  use_technical = FALSE, max_na = 0.2, skew_log = 2, pca_var = 0.9, pcs_max = 6, k_space = 30,
  breadth_vars = character(0),
  taxonomy_file = file.path(data_folder, "microbeatlas", "taxonomy_soilCurated_OTU97_nextMAPrelease.tsv"),
  tax_rank = NULL, tax_fill = TRUE, tax_exclude = "chloroplast|mitochond|Eukaryota|Eukarya", tax_id_col = NULL, tax_lineage_col = NULL, tax_skip = 0,   # tax_rank = NULL: the OTUs themselves
  min_samples = 5, min_clusters = 5,       # an OTU must be present in at least this many samples and clusters (of the analysed ones) to be tested
  z_cut = 1.64, min_effect = 0.10,         # narrow: z <= -z_cut and a breadth at least min_effect below the null mean
  alpine_share_min = 0.80, min_occ_alpine = 0.05, min_draw_frac = 0.5,   # alpine: at least this share of the abundance in alpine units, present in this share of the alpine units, in this share of the draws
  min_otus_family = 3, max_families = 60, n_parts = 6, pseudocount = 0.5,
  depth_unit = "auto",                      # the unit of Depth_m values that carry none: "m", "cm" or "auto" (up to 3 = m)
  gmba_path = file.path(data_folder, "GMBA"), range_level = NA, range_max_km = 50, min_range_clusters = 3,   # the mountain ranges (as in script 3); a range needs at least this many alpine clusters
  split_units = TRUE,                      # the units of the draws: an alpine cluster is split by mountain range and a lowland cluster by 5-degree cell (FALSE: one unit per cluster)
  max_per_cluster = 10, min_n = 150, min_depth_n = 100, max_lowland_ref = 2000, n_range_labels = 6,
  cols = c(Alpine = "#2C5F9E", Lowland = "#7A9A7E")
)
if (TEST_MODE) { out_folder <- paste0(out_folder, "_TEST"); CFG$n_draws <- 2; CFG$nperm <- 29; CFG$lowland_reps <- 2; CFG$k_space <- 15 }
dir.create(file.path(out_folder, "tables"), recursive = TRUE, showWarnings = FALSE)
out_t <- function(...) file.path(out_folder, "tables", ...)
set.seed(1)

theme_s <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), plot.tag = element_text(face = "bold", size = 14), plot.title = element_text(size = 9.5, face = "bold"))

theme_fig <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), plot.tag = element_text(face = "bold", size = 16), strip.background = element_rect(fill = "grey95"), strip.text = element_text(face = "bold"))
to_sparse <- function(m) {
  if (inherits(m, "simple_triplet_matrix")) Matrix::sparseMatrix(i = m$i, j = m$j, x = m$v, dims = c(m$nrow, m$ncol), dimnames = list(m$dimnames[[1]], m$dimnames[[2]]))
  else as(m, "dgCMatrix")
}
cluster_ids <- function(t_g) { cl <- as.character(t_g[[CFG$cluster_col]]); no <- is.na(cl) | !nzchar(cl); cl[no] <- paste0("none", which(no)); cl }
yes <- function(x) tolower(trimws(as.character(x))) %in% c("yes", "true", "1", "y")
skew <- function(x) mean(((x - mean(x)) / sd(x))^3)
clean_env <- function(M) {
  M <- as.matrix(M); M <- M[, colMeans(is.na(M)) <= CFG$max_na, drop = FALSE]
  for (j in seq_len(ncol(M))) { x <- M[, j]; x[is.na(x)] <- median(x, na.rm = TRUE); if (sd(x) > 0 && min(x) >= 0 && skew(x) > CFG$skew_log) x <- log1p(x); M[, j] <- x }
  M <- M[, apply(M, 2, sd) > 1e-8, drop = FALSE]; if (ncol(M)) scale(M) else M
}
env_pcs <- function(M, cap) {                                                                           # a set reduced to its first principal components (no use of the response)
  M <- M[, apply(M, 2, sd) > 1e-8, drop = FALSE]; if (ncol(M) == 0 || nrow(M) < 4) return(NULL)
  pc <- prcomp(M, scale. = TRUE); cum <- cumsum(pc$sdev^2) / sum(pc$sdev^2)
  k <- min(cap, max(2, which(cum >= CFG$pca_var)[1], na.rm = TRUE), ncol(pc$x)); scale(pc$x[, seq_len(k), drop = FALSE], scale = FALSE)
}

# The taxonomy table of the OTUs (taxonomy_soilCurated_OTU97_nextMAPrelease.tsv): TAB-separated, a header (OID, OTU, Tax, ..., OrigTax, ..., RemapSeqs), long fields (lists of species /
# accessions) and the sequence. Read from a sample of its lines:
#   * the OTU ID is the column OID (or, without header, the column whose values match the OTUs of the BIOM, exactly or by the number at the end of the id: "97_123" ~ "123");
#   * the LINEAGE is the column Tax (or, without header, the short column with the most ";" separators); the ranks are its first six elements after tax_skip, with [organelle]
#     tags, (confidence) values and prefixes (d__) removed. The raw Tax and OrigTax strings are kept to recognise chloroplasts, mitochondria and eukaryotes.
# tax_id_col / tax_lineage_col (column numbers) override the detection. Only these columns are read from the whole file.
read_taxonomy <- function(file, otu_ids) {
  smp <- readLines(file, n = 400, warn = FALSE); fl <- strsplit(smp, "\t", fixed = TRUE); nc <- max(lengths(fl))
  if (nc < 2) { message("The taxonomy file has no tab-separated columns: ", file); return(NULL) }
  M <- t(vapply(fl, function(z) { length(z) <- nc; z[is.na(z)] <- ""; z }, character(nc)))
  hdr <- tolower(trimws(M[1, ])); has_hdr <- sum(hdr %in% c("oid", "otu", "tax", "taxonomy", "lineage", "origtax", "repsequenceid", "seqcount", "remapseqs")) >= 2
  message(sprintf("Taxonomy file %s: %d columns%s", basename(file), nc, if (has_hdr) paste0("; header: ", paste(substr(M[1, seq_len(min(nc, 21))], 1, 14), collapse = " ")) else paste0("; the last four fields of the first line: ", paste(substr(M[1, max(1, nc - 3):nc], 1, 40), collapse = " | "))))
  if (has_hdr) M <- M[-1, , drop = FALSE]
  nkey <- function(x) sub("^.*?([0-9]+)$", "\\1", x, perl = TRUE); okey <- nkey(otu_ids)
  uniq <- vapply(seq_len(nc), function(j) length(unique(M[, j])) / nrow(M), 1)                          # an id column has a different value in every row (a length column does not)
  exact <- vapply(seq_len(nc), function(j) mean(M[, j] %in% otu_ids), 1); numk <- vapply(seq_len(nc), function(j) mean(grepl("[0-9]+$", M[, j]) & nkey(M[, j]) %in% okey), 1)
  last_max <- function(v) { w <- which(v >= max(v) - 1e-9); w[length(w)] }
  h_id <- if (has_hdr) match(intersect(c("oid", "otu_id", "otuid", "id"), hdr)[1], hdr) else NA_integer_
  if (!is.null(CFG$tax_id_col)) { idc <- CFG$tax_id_col; use_exact <- exact[idc] >= 0.3
  } else if (!is.na(h_id) && max(exact[h_id], numk[h_id]) >= 0.3) { idc <- h_id; use_exact <- exact[h_id] >= 0.3
  } else { use_exact <- max(exact * uniq) >= 0.3; idc <- if (use_exact) last_max(exact * uniq) else last_max(numk * uniq)
    if (max(exact * uniq, numk * uniq) < 0.3) { message(sprintf("  no column of the taxonomy file matches the OTU ids of the BIOM (best: %.0f%%); examples of the BIOM ids: %s; set tax_id_col: no taxonomic analysis", 100 * max(exact * uniq, numk * uniq), paste(head(otu_ids, 3), collapse = ", "))); return(NULL) } }
  cand <- setdiff(seq_len(nc), idc)
  sc <- vapply(cand, function(j) mean(grepl(";", M[, j], fixed = TRUE)), 1); len <- vapply(cand, function(j) median(nchar(M[, j])), 1)
  nsep <- vapply(cand, function(j) mean(nchar(M[, j]) - nchar(gsub(";", "", M[, j], fixed = TRUE))), 1); okc <- sc > 0.5 & len < 600
  h_lin <- if (has_hdr) match(intersect(c("tax", "taxonomy", "lineage"), hdr)[1], hdr) else NA_integer_
  lc <- if (!is.null(CFG$tax_lineage_col)) CFG$tax_lineage_col else if (!is.na(h_lin)) h_lin else if (any(okc)) cand[okc][which.max(nsep[okc])] else NA_integer_
  if (is.na(lc)) { message("  no column of the taxonomy file looks like a lineage (short strings with ';'); set tax_lineage_col: no taxonomic analysis"); return(NULL) }
  oc <- if (has_hdr) match("origtax", hdr) else NA_integer_
  cc <- rep("NULL", nc); cc[na.omit(unique(c(idc, lc, oc)))] <- "character"
  tx <- read.delim(file, header = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE, colClasses = cc, fill = TRUE, na.strings = character(0), col.names = paste0("V", seq_len(nc)))
  id_f <- tx[[paste0("V", idc)]]; raw_f <- tx[[paste0("V", lc)]]; orig_f <- if (!is.na(oc)) tx[[paste0("V", oc)]] else rep("", nrow(tx))
  bi <- if (use_exact) match(id_f, otu_ids) else match(nkey(id_f), okey)
  ok <- !is.na(bi) & !duplicated(bi); id_b <- otu_ids[bi[ok]]; raw <- raw_f[ok]; excl <- paste(raw_f[ok], orig_f[ok])
  message(sprintf("  OTU id = column %d%s (%s), lineage = column %d%s%s; %d of the %d OTUs of the BIOM have a taxonomy", idc, if (has_hdr) paste0(" '", hdr[idc], "'") else "", if (use_exact) "same ids as the BIOM" else "matched by the number at the end of the id",
                  lc, if (has_hdr) paste0(" '", hdr[lc], "'") else "", if (!is.na(oc)) paste0(", original taxonomy = column ", oc) else "", length(id_b), length(otu_ids)))
  clean <- function(z) { z <- trimws(z); z <- gsub("\\[[^]]*\\]", "", z); z <- sub("\\s*\\([0-9.]+\\)$", "", z); z <- sub("^[A-Za-z]__", "", z)
    z[!nzchar(z) | grepl("^(uncultured|unidentified|unclassified|incertae sedis|unassigned|NA|environmental sample|metagenome)", z, ignore.case = TRUE)] <- NA_character_; z }
  parts <- strsplit(raw, ";", fixed = TRUE); if (CFG$tax_skip > 0) parts <- lapply(parts, function(z) if (length(z) > CFG$tax_skip) z[-seq_len(CFG$tax_skip)] else character(0))
  ranks <- matrix(t(vapply(parts, function(z) { z <- clean(z); length(z) <- 6; z }, character(6))), ncol = 6, dimnames = list(NULL, c("domain", "phylum", "class", "order", "family", "genus")))
  set.seed(2); ex <- sample(length(raw), min(3, length(raw)))
  for (i in ex) message(sprintf("  example: OTU %s | %s  ->  %s", id_b[i], substr(raw[i], 1, 80), paste(ifelse(is.na(ranks[i, ]), "-", ranks[i, ]), collapse = " > ")))
  message(sprintf("  OTUs with a phylum / family / genus: %.0f%% / %.0f%% / %.0f%%; domains: %s", 100 * mean(!is.na(ranks[, 2])), 100 * mean(!is.na(ranks[, 5])), 100 * mean(!is.na(ranks[, 6])),
                  paste(head(paste0(names(sort(table(ranks[, 1]), decreasing = TRUE)), " ", sort(table(ranks[, 1]), decreasing = TRUE)), 4), collapse = ", ")))
  list(ids = id_b, ranks = ranks, raw = raw, excl = excl)
}
tax_agg <- function(tax, rank, ids, fill = CFG$tax_fill) {                                               # an indicator matrix taxa x OTUs for a rank
  rk <- if (is.numeric(rank)) rank else match(tolower(rank), c("domain", "phylum", "class", "order", "family", "genus")); if (is.na(rk)) stop("tax_rank must be domain, phylum, class, order, family or genus")
  R <- tax$ranks[match(ids, tax$ids), , drop = FALSE]; has <- !is.na(R[, rk]); lab <- rep(NA_character_, length(ids)); filled <- rep(FALSE, length(ids))
  lab[has] <- apply(R[has, seq_len(rk), drop = FALSE], 1, function(z) paste(z[!is.na(z)], collapse = ";"))
  if (isTRUE(fill)) {                                                                                  # lineage stops above the rank: an 'unclassified' taxon of the lowest known rank
    rest <- which(!has & !is.na(R[, 1]))
    if (length(rest)) { low <- vapply(rest, function(i) max(which(!is.na(R[i, seq_len(rk)]))), 1L)
      lab[rest] <- vapply(seq_along(rest), function(k) paste0(paste(R[rest[k], seq_len(low[k])][!is.na(R[rest[k], seq_len(low[k])])], collapse = ";"), ";(unclassified ", colnames(R)[rk], ")"), ""); filled[rest] <- TRUE }
  }
  assigned <- !is.na(lab); taxa <- sort(unique(lab[assigned]))
  list(Agg = Matrix::sparseMatrix(i = match(lab[assigned], taxa), j = which(assigned), x = 1, dims = c(length(taxa), length(ids)), dimnames = list(taxa, ids)), assigned = assigned, filled = filled, n_taxa = length(taxa), lab = lab, R = R)
}


# Continents by boxes of latitude / longitude, evaluated in this order (a rule of thumb: check tables/umap_samples.csv and edit the boxes if
# an alpine sample is on the wrong continent; samples that fit nowhere are "Other")
continent_of <- function(lat, lon) {
  out <- rep("Other", length(lat))
  set <- function(cond, name) out[cond & out == "Other"] <<- name
  set(lat < -60, "Antarctica")
  set(lon >= -161 & lon <= -154 & lat >= 18 & lat <= 23, "Oceania")                                  # Hawaii
  set(lon >= -170 & lon <= -50 & lat >= 12, "North America")                                         # incl. Mexico, Central America north of 12 N, Greenland west of 50 W
  set(lon >= -85 & lon <= -30 & lat < 12, "South America")
  set((lon >= 110 & lat < -9) | (lon >= 140 & lat < 0), "Oceania")
  set((lon >= -20 & lon <= 34 & lat >= -36 & lat <= 37) | (lon > 34 & lon <= 52 & lat >= -36 & lat <= 12), "Africa")
  set(lat >= 35 & lon >= -25 & lon <= 45 & !(lon >= 26 & lat < 43), "Europe")                        # the Caucasus is counted to Europe, Anatolia to Asia
  set(lat > -10 & lon > 26, "Asia")                                                                  # incl. the Middle East, the Himalaya, Siberia, Japan
  out
}
# The regions of the grid: the ranges of the GMBA Mountain Inventory v2 (polygons of one hierarchical selection, in an equal-area projection) or, without the
# file, the continents
load_gmba <- function() {
  p <- CFG$gmba_path
  if (is.null(p) || !nzchar(p) || !(file.exists(p) || dir.exists(p))) { message("GMBA Mountain Inventory not found (", p, "): the ranges are the continents"); return(NULL) }
  if (!requireNamespace("sf", quietly = TRUE)) { message("the package sf is not available (conda: r-sf): the ranges are the continents"); return(NULL) }
  shp <- p
  if (dir.exists(p)) {
    f <- list.files(p, pattern = "\\.shp$", recursive = TRUE, full.names = TRUE)
    f_std <- f[grepl("standard", basename(f), ignore.case = TRUE) & !grepl("basic|broad", basename(f), ignore.case = TRUE)]
    shp <- if (length(f_std)) f_std[1] else if (length(f)) f[1] else NA_character_
  }
  if (is.na(shp)) { message("no .shp file in ", p, ": the ranges are the continents"); return(NULL) }
  cache <- file.path(out_folder, "gmba_ranges.rds"); key <- list(shp = shp, size = file.size(shp), level = CFG$range_level)
  if (REUSE_FITS && file.exists(cache)) { cc <- readRDS(cache); if (identical(cc$key, key)) { message("Reusing the saved GMBA ranges: ", cache); return(cc$polys) } }
  message("Reading the GMBA Mountain Inventory (slow, once): ", shp)
  g <- sf::st_read(shp, quiet = TRUE)
  if (!"DBaseName" %in% names(g)) { message("the GMBA file lacks the column DBaseName: the ranges are the continents"); return(NULL) }
  sel <- if (is.na(CFG$range_level)) {
    if ("Select_300" %in% names(g)) { sc <- trimws(as.character(g$Select_300)); !is.na(sc) & nzchar(sc) & !(toupper(sc) %in% c("0", "FALSE", "F", "N", "NO", "NA")) }
    else if (grepl("300", basename(shp))) rep(TRUE, nrow(g))                                        # the layer of the 292 selected ranges: all of its polygons
    else { message("the GMBA file lacks the column Select_300: set range_level to a hierarchical level (1-10)"); return(NULL) }
  } else { if (!"Hier_Lvl" %in% names(g)) { message("the GMBA file lacks the column Hier_Lvl"); return(NULL) }; !is.na(g$Hier_Lvl) & g$Hier_Lvl == CFG$range_level }
  g <- g[sel, ]
  if (!nrow(g)) { message("no GMBA polygon selected (range_level = ", CFG$range_level, "): the ranges are the continents"); return(NULL) }
  old <- suppressMessages(sf::sf_use_s2(FALSE)); on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  g <- sf::st_make_valid(g)
  g <- sf::st_transform(g, "+proj=moll +lon_0=0 +datum=WGS84 +units=m +no_defs")                     # equal-area and planar: containment, nearest polygon, area
  polys <- sf::st_sf(name = as.character(g$DBaseName), area = as.numeric(sf::st_area(g)), geometry = sf::st_geometry(g))
  message(sprintf("GMBA: %d polygons selected (%s)", nrow(polys), if (is.na(CFG$range_level)) "Select_300" else paste("hierarchical level", CFG$range_level)))
  saveRDS(list(key = key, polys = polys), cache)
  polys
}
gmba_polys <- load_gmba()
region_kind <- if (!is.null(gmba_polys)) "mountain ranges (GMBA Mountain Inventory v2; Snethlage et al. 2022)" else "continents"
region_of <- if (!is.null(gmba_polys)) function(lon, lat) {                                          # the range of every point
  pts <- sf::st_transform(sf::st_as_sf(data.frame(lon = lon, lat = lat), coords = c("lon", "lat"), crs = 4326), sf::st_crs(gmba_polys))
  hit <- sf::st_intersects(pts, gmba_polys)
  out <- vapply(hit, function(h) if (length(h)) gmba_polys$name[h[which.min(gmba_polys$area[h])]] else NA_character_, "")   # the smallest polygon that contains the point
  miss <- which(is.na(out))
  if (length(miss) && is.finite(CFG$range_max_km)) {
    nn <- sf::st_nearest_feature(pts[miss, ], gmba_polys)
    d <- as.numeric(sf::st_distance(pts[miss, ], gmba_polys[nn, ], by_element = TRUE)) / 1000
    out[miss[d <= CFG$range_max_km]] <- gmba_polys$name[nn[d <= CFG$range_max_km]]
  }
  out[is.na(out)] <- "Other"; out
} else function(lon, lat) continent_of(lat, lon)

## ---- 1. Data ------------------------------------------------------------------------------------------------------------------------------
raw <- read.delim(CFG$sample_file, quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE, fill = FALSE)
need <- c(CFG$id_table, CFG$biom_id, CFG$group_col, CFG$cluster_col, "Latitude", "Longitude", "alpha_total_n_reads")
miss <- setdiff(need, names(raw)); if (length(miss)) stop("Missing columns in ", CFG$sample_file, ": ", paste(miss, collapse = ", "))
tab <- raw[raw[[CFG$group_col]] %in% c(CFG$alpine_value, CFG$lowland_value), ]
tab <- tab[is.finite(tab$Latitude) & is.finite(tab$Longitude) & is.finite(tab$alpha_total_n_reads) & tab$alpha_total_n_reads >= CFG$rare_depth, ]
env <- read.csv(CFG$env_file, check.names = FALSE, stringsAsFactors = FALSE); env <- env[!duplicated(env[[CFG$id_env]]), ]
mm <- match(tab[[CFG$id_table]], env[[CFG$id_env]]); tab <- tab[!is.na(mm), ]; env <- env[mm[!is.na(mm)], , drop = FALSE]
is_alp <- tab[[CFG$group_col]] == CFG$alpine_value
rep_cols <- head(sort(grep("^rep_[0-9]+$", names(env), value = TRUE)), CFG$lowland_reps)
alp_ix <- which(is_alp); low_ix <- sort(unique(unlist(lapply(rep_cols, function(rc) which(!is_alp & yes(env[[rc]]))))))
message(sprintf("%d alpine samples and %d lowland samples (the union of %d replicates) have environmental data", length(alp_ix), length(low_ix), length(rep_cols)))
if (length(alp_ix) < 50) stop("fewer than 50 alpine samples with environmental data: check CFG$id_table / CFG$id_env")
has_low <- length(low_ix) >= 50; if (!has_low) { low_ix <- integer(0); message("No lowland replicate with environmental data (columns rep_XX): alpine only") }

message("Reading the BIOM: ", CFG$biom_file)
use <- c(alp_ix, low_ix); biom <- rbiom::read_biom(CFG$biom_file); ids <- intersect(tab[[CFG$biom_id]][use], biom$samples); biom <- biom[ids]
X <- to_sparse(biom$counts); X <- X[rownames(X) != "Unmapped", , drop = FALSE]; rm(biom); invisible(gc())
col_of <- match(tab[[CFG$biom_id]], colnames(X))
tax <- NULL
if (!is.null(CFG$taxonomy_file)) { if (file.exists(CFG$taxonomy_file)) tax <- read_taxonomy(CFG$taxonomy_file, rownames(X)) else message("Taxonomy file not found (", CFG$taxonomy_file, "): the OTUs are used") }
if (!is.null(tax) && !is.null(CFG$tax_exclude) && nzchar(CFG$tax_exclude)) {                              # plant / organelle DNA is not part of the soil bacterial community
  ex <- rownames(X) %in% tax$ids[grepl(CFG$tax_exclude, tax$excl, ignore.case = TRUE)]
  for (gp in c("Alpine", "Lowland")) { cols <- unique(col_of[if (gp == "Alpine") alp_ix else low_ix]); cols <- cols[!is.na(cols)]; if (!length(cols)) next
    v <- Matrix::colSums(X[ex, cols, drop = FALSE]) / pmax(Matrix::colSums(X[, cols, drop = FALSE]), 1)
    message(sprintf("  chloroplast / mitochondria / eukaryote OTUs (%d): %s samples: median %.1f%% of the reads; %.0f%% of the samples have more than 20%%, %.0f%% more than 50%%", sum(ex), gp, 100 * median(v), 100 * mean(v > 0.2), 100 * mean(v > 0.5))) }
  X <- X[!ex, , drop = FALSE]; message(sprintf("  these OTUs are removed: %d OTUs remain", nrow(X)))
}
tot <- Matrix::colSums(X); good <- !is.na(col_of) & tot[pmax(col_of, 1)] >= CFG$rare_depth
alp_ix <- alp_ix[good[alp_ix]]; low_ix <- low_ix[good[low_ix]]; if (length(low_ix) < 50) { has_low <- FALSE; low_ix <- integer(0) }
if (length(alp_ix) < 50) stop("fewer than 50 alpine samples with reads in the BIOM")
Prel <- X %*% Matrix::Diagonal(x = 1 / tot); Xcnt <- X; rm(X)                                                    # (the counts are kept for the CLR of panel C)
if (!is.null(tax) && !is.null(CFG$tax_rank)) {                                                           # the OTUs aggregated to a taxonomic rank
  ag <- tax_agg(tax, CFG$tax_rank, rownames(Prel)); mass <- Matrix::rowSums(Prel)
  message(sprintf("Aggregated to %s: %d taxa from %d OTUs; %.1f%% of the reads have a %s, %.1f%% %s; %.1f%% have no lineage and are left out", CFG$tax_rank, ag$n_taxa, sum(ag$assigned), 100 * sum(mass[ag$assigned & !ag$filled]) / sum(mass), CFG$tax_rank,
                  100 * sum(mass[ag$filled]) / sum(mass), if (isTRUE(CFG$tax_fill)) "are in 'unclassified' taxa of a higher rank" else "are left out", 100 * sum(mass[!ag$assigned]) / sum(mass)))
  Prel <- ag$Agg %*% Prel; ag_lab <- ag$lab; ag_R <- ag$R; rm(ag)
}
unit_lab <- if (!is.null(tax) && !is.null(CFG$tax_rank)) paste0(CFG$tax_rank, " taxa") else "OTUs"
an_ix <- sort(unique(c(alp_ix, low_ix))); cl_all <- cluster_ids(tab)
Wind <- Matrix::sparseMatrix(i = seq_along(an_ix), j = match(cl_all[an_ix], sort(unique(cl_all[an_ix]))), x = 1, dims = c(length(an_ix), length(unique(cl_all[an_ix]))))
Pan <- Prel[, col_of[an_ix], drop = FALSE] > 0; n_s <- Matrix::rowSums(Pan); n_c <- Matrix::rowSums((Pan %*% Wind) > 0)
taxa <- which(n_s >= CFG$min_samples & n_c >= CFG$min_clusters)
message(sprintf("Prevalence filter: %d of %d %s are present in at least %d samples and %d clusters (of %d analysed samples)", length(taxa), nrow(Prel), unit_lab, CFG$min_samples, CFG$min_clusters, length(an_ix)))
if (length(taxa) < 20) stop("fewer than 20 taxa pass the prevalence filter")

# the variables (all analysed samples; cleaned once)
Zs <- list()
for (nm in names(CFG$env_sets)) { cols <- unique(unlist(lapply(CFG$env_sets[[nm]], function(p) grep(p, names(env), value = TRUE)))); cols <- cols[vapply(env[cols], is.numeric, TRUE)]
  M <- clean_env(env[an_ix, cols, drop = FALSE]); if (ncol(M)) { Zs[[nm]] <- M; message(sprintf("  %s: %d of %d variables kept", nm, ncol(M), length(cols))) } }
if (!length(Zs)) stop("no environmental variables found: check CFG$env_sets and the columns of ", CFG$env_file)
tech <- NULL
if (isTRUE(CFG$use_technical) && all(c("unmapped_raw_counts", "total_raw_counts") %in% names(tab))) {
  tv <- cbind(log10_reads = log10(tab$alpha_total_n_reads[an_ix]), unmapped_share = tab$unmapped_raw_counts[an_ix] / tab$total_raw_counts[an_ix])
  for (j in 1:2) tv[!is.finite(tv[, j]), j] <- median(tv[is.finite(tv[, j]), j]); tech <- scale(tv)
}
Zv <- list(); for (v in intersect(CFG$breadth_vars, names(env))) { M <- clean_env(env[an_ix, v, drop = FALSE]); if (ncol(M)) Zv[[v]] <- as.numeric(M[, 1]) }          # single variables for the niche breadth
set_of <- c(unlist(lapply(names(Zs), function(nm) setNames(rep(nm, ncol(Zs[[nm]])), colnames(Zs[[nm]])))), if (!is.null(tech)) setNames(rep("Technical", 2), colnames(tech)), Latitude = "Geography", Longitude = "Geography")

## ---- 2. The random draws: one sample per cluster ------------------------------------------------------------------------------------------
# UNITS of the draws. A cluster (samples of one BioProject or within 100 m of each other) can span several regions: with split_units = TRUE an alpine cluster is split by MOUNTAIN RANGE (a sample is given to the
# range it lies in, as in panel B) and a lowland cluster by 5-degree cell, so that a multi-region study contributes one sample PER REGION to every draw, each at its own place. split_units = FALSE: one unit per cluster.
gmba_polys <- load_gmba(); region_kind <- if (!is.null(gmba_polys)) "mountain ranges (GMBA)" else "continents"
rng_tab <- rep(NA_character_, nrow(tab)); rng_tab[alp_ix] <- region_of(tab$Longitude[alp_ix], tab$Latitude[alp_ix])
cell_of <- function(ix) paste(floor(tab$Latitude[ix] / 5), floor(tab$Longitude[ix] / 5))
alp_cl <- split(alp_ix, if (CFG$split_units) paste(cl_all[alp_ix], rng_tab[alp_ix], sep = " | ") else cl_all[alp_ix]); n_u <- length(alp_cl); pick1 <- function(i) i[sample.int(length(i), 1)]
draw_alpine <- function() unname(vapply(alp_cl, pick1, 1L))
if (has_low) { low_cl <- split(low_ix, if (CFG$split_units) paste(cl_all[low_ix], cell_of(low_ix), sep = " | ") else cl_all[low_ix])
  cen <- data.frame(lat = vapply(low_cl, function(i) mean(tab$Latitude[i]), 1), lon = vapply(low_cl, function(i) mean(tab$Longitude[i]), 1)); cell <- paste(floor(cen$lat / 5), floor(cen$lon / 5))
  draw_lowland <- function() { o <- sample(nrow(cen)); rk <- ave(seq_along(o), cell[o], FUN = seq_along); ord <- o[order(rk, match(cell[o], sample(unique(cell[o]))))][seq_len(min(n_u, nrow(cen)))]; unname(vapply(low_cl[ord], pick1, 1L)) } }
message(sprintf("Draws: %d; every draw: %d alpine units (%d clusters%s)%s", CFG$n_draws, n_u, length(unique(cl_all[alp_ix])), if (CFG$split_units) ", split by mountain range" else "", if (has_low) sprintf(" + %d lowland units (of %d from %d clusters%s), one random sample each", min(n_u, length(low_cl)), length(low_cl), length(unique(cl_all[low_ix])), if (CFG$split_units) ", split by 5-degree cell" else "") else ""))
unit_data <- function(idx) {
  r <- match(idx, an_ix); P <- Prel[taxa, col_of[idx], drop = FALSE]; s <- Matrix::colSums(P); ok <- s > 0; idx <- idx[ok]; r <- r[ok]
  P <- as.matrix(Matrix::t(P[, ok, drop = FALSE] %*% Matrix::Diagonal(x = 1 / s[ok])))                  # units x taxa: relative abundances among the kept taxa
  list(idx = idx, r = r, P = P, Z = lapply(Zs, function(M) M[r, , drop = FALSE]), tech = if (!is.null(tech)) tech[r, , drop = FALSE], ll = cbind(tab$Longitude[idx], tab$Latitude[idx]))
}




## ---- 3. Niche breadth and alpine affinity of every OTU in every draw ---------------------------------------------------------------------------------------
niche_otu <- function(P, X, isA, nperm) {                         # P: units x OTUs (relative abundances), X: units x axes (z-scored), isA: alpine units
  W <- t(P); W <- W / pmax(rowSums(W), 1e-300); n <- ncol(W); X2 <- X^2; a <- as.numeric(isA)
  obs <- rowMeans(W %*% X2 - (W %*% X)^2); sh <- as.numeric(W %*% a)
  s1 <- s2 <- h1 <- h2 <- numeric(nrow(W))
  for (b in seq_len(nperm)) { idx <- sample.int(n); v <- rowMeans(W %*% X2[idx, , drop = FALSE] - (W %*% X[idx, , drop = FALSE])^2); s1 <- s1 + v; s2 <- s2 + v^2; h <- as.numeric(W %*% a[idx]); h1 <- h1 + h; h2 <- h2 + h^2 }
  mu <- s1 / nperm; hm <- h1 / nperm
  list(z = (obs - mu) / sqrt(pmax(s2 / nperm - mu^2, 1e-12)), ratio = obs / pmax(mu, 1e-12), share = sh, z_share = (sh - hm) / sqrt(pmax(h2 / nperm - hm^2, 1e-12)),
       marginality = sqrt(rowSums((W %*% X)^2) / ncol(X)), occ_alp = colMeans(P[isA, , drop = FALSE] > 0), occ_low = colMeans(P[!isA, , drop = FALSE] > 0))
}
res <- list(); t0 <- Sys.time(); message(sprintf("Niche breadth of %d OTUs in %d draws (%d permutations) ...", length(taxa), CFG$n_draws, CFG$nperm))
for (r in seq_len(CFG$n_draws)) {
  a <- draw_alpine(); l <- draw_lowland(); u <- unit_data(c(a, l)); sets <- lapply(u$Z, env_pcs, cap = CFG$pcs_max); sets <- sets[!vapply(sets, is.null, TRUE)]; X <- do.call(cbind, lapply(sets, scale))
  nz <- niche_otu(u$P, X, u$idx %in% a, CFG$nperm)
  res[[r]] <- data.frame(draw = r, otu = colnames(u$P), z = nz$z, ratio = nz$ratio, share = nz$share, z_share = nz$z_share, marg = nz$marginality, occ_alp = nz$occ_alp, occ_low = nz$occ_low, stringsAsFactors = FALSE)
  if (r == 1 || r %% 5 == 0 || r == CFG$n_draws) message(sprintf("  draw %d of %d (%.1f min)", r, CFG$n_draws, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
res <- bind_rows(res); res <- res[is.finite(res$z) & is.finite(res$share), ]
res$narrow <- res$z <= -CFG$z_cut & res$ratio <= 1 - CFG$min_effect
res$alp_spec <- res$narrow & res$share >= CFG$alpine_share_min & res$z_share >= CFG$z_cut; res$low_spec <- res$narrow & res$share <= 1 - CFG$alpine_share_min & res$z_share <= -CFG$z_cut
otu <- res %>% group_by(otu) %>% summarise(draws = n(), z_breadth = mean(z), breadth_ratio = mean(ratio), share_alpine = mean(share), z_share = mean(z_share), marginality = mean(marg), occ_alpine = mean(occ_alp), occ_lowland = mean(occ_low),
                                           frac_alpine_specialist = mean(alp_spec), frac_lowland_specialist = mean(low_spec), frac_narrow = mean(narrow), .groups = "drop") %>% as.data.frame()
enough <- otu$draws >= CFG$n_draws / 2
otu$category <- ifelse(enough & otu$frac_alpine_specialist >= CFG$min_draw_frac & otu$occ_alpine >= CFG$min_occ_alpine, "Alpine specialist",
                ifelse(enough & otu$frac_lowland_specialist >= CFG$min_draw_frac & otu$occ_lowland >= CFG$min_occ_alpine, "Lowland specialist", ifelse(enough & otu$frac_narrow >= CFG$min_draw_frac, "Other specialist", "Not specialist")))
R6 <- tax$ranks[match(otu$otu, tax$ids), seq_len(6), drop = FALSE]; R6[is.na(R6)] <- ""; rk_names <- c("domain", "phylum", "class", "order", "family", "genus"); for (k in 1:6) otu[[rk_names[k]]] <- R6[, k]
write.csv(otu, out_t("otu_niche_breadth.csv"), row.names = FALSE)
ct <- table(otu$category); message(sprintf("%d OTUs tested: %s", nrow(otu), paste(sprintf("%s %d (%.1f%%)", names(ct), as.integer(ct), 100 * as.integer(ct) / nrow(otu)), collapse = ", ")))
spec <- otu[otu$category == "Alpine specialist", ]; if (nrow(spec) < 5) stop("fewer than 5 alpine specialist OTUs: the figure needs more (check the criteria and the numbers above)")
message(sprintf("  alpine specialists: median share of the abundance in alpine units %.2f, present in a median of %.0f%% of the alpine units; phyla: %s", median(spec$share_alpine), 100 * median(spec$occ_alpine),
                paste(head(sprintf("%s %d", names(sort(table(ifelse(nzchar(spec$phylum), spec$phylum, "(no phylum)")), decreasing = TRUE)), as.integer(sort(table(ifelse(nzchar(spec$phylum), spec$phylum, "(no phylum)")), decreasing = TRUE))), 6), collapse = ", ")))

## ---- 4. Panel A: the taxonomy of the specialists ----------------------------------------------------------------------------------------------------------
tested <- otu[otu$category %in% c("Alpine specialist", "Lowland specialist", "Other specialist", "Not specialist"), ]; tested$is_spec <- tested$category == "Alpine specialist"
fam_lab <- ifelse(nzchar(tested$family), paste(tested$domain, tested$phylum, tested$class, tested$order, tested$family, sep = ";"), NA_character_)
fe <- tested[!is.na(fam_lab), ] %>% mutate(fam = fam_lab[!is.na(fam_lab)]) %>% group_by(fam, domain, phylum, class, order, family) %>% summarise(n_spec = sum(is_spec), n_otus = n(), .groups = "drop") %>% as.data.frame()
tot_s <- sum(fe$n_spec); tot_o <- sum(fe$n_otus); fe$share <- fe$n_spec / fe$n_otus
fe$p_enrich <- mapply(function(a, n) fisher.test(matrix(c(a, n - a, tot_s - a, (tot_o - n) - (tot_s - a)), 2), alternative = "greater")$p.value, fe$n_spec, fe$n_otus); fe$fdr <- p.adjust(fe$p_enrich, "BH")
fe <- fe[order(-fe$n_spec, fe$p_enrich), ]; write.csv(fe, out_t("family_enrichment.csv"), row.names = FALSE)
fs <- fe[fe$n_spec >= CFG$min_otus_family, ]; fs <- head(fs, CFG$max_families)
message(sprintf("Panel A: %d families have at least %d alpine specialist OTUs (%d shown); enriched in specialists (FDR < 0.05): %d", sum(fe$n_spec >= CFG$min_otus_family), CFG$min_otus_family, nrow(fs), sum(fs$fdr < 0.05)))
if (nrow(fs) < 3) { message("  too few families for the tree: the tree is drawn from the families with at least 1 specialist OTU"); fs <- head(fe[fe$n_spec >= 1, ], CFG$max_families) }
fs <- fs[order(fs$domain, fs$phylum, fs$class, fs$order, fs$family), ]; fs$phy <- ifelse(nzchar(fs$phylum), fs$phylum, "(no phylum)")
phy_all <- unique(ifelse(nzchar(spec$phylum), spec$phylum, "(no phylum)")); pal <- setNames(grDevices::hcl.colors(max(3, length(phy_all)), "Dark 3")[seq_along(phy_all)], phy_all)
tree_plots <- function(fm) {
  n <- nrow(fm); fm$y <- n:1; rn <- c("domain", "phylum", "class", "order")
  for (k in 1:4) { pk <- fm[[rn[k]]]; fm[[paste0("p", k)]] <- ifelse(nzchar(pk), pk, paste0("(no ", rn[k], ")")) }; for (k in 1:4) fm[[paste0("k", k)]] <- do.call(paste, c(lapply(1:k, function(j) fm[[paste0("p", j)]]), sep = "|"))
  vseg <- list(); hseg <- list(); child <- data.frame(key = fm$fam, y = fm$y, parent = fm$k4, phy = fm$phy, stringsAsFactors = FALSE)
  for (d in 4:1) { nd <- child %>% group_by(parent) %>% summarise(ymin = min(y), ymax = max(y), yn = mean(y), phy = first(phy), .groups = "drop")
    vseg[[d]] <- data.frame(x = d, ymin = nd$ymin, ymax = nd$ymax, phy = nd$phy); hseg[[d]] <- data.frame(x = d, xend = d + 1, y = child$y, phy = child$phy)
    if (d > 1) child <- data.frame(key = nd$parent, y = nd$yn, parent = fm[[paste0("k", d - 1)]][match(nd$parent, fm[[paste0("k", d)]])], phy = nd$phy, stringsAsFactors = FALSE) else root <- nd }
  vseg <- bind_rows(vseg); hseg <- bind_rows(hseg); stub <- data.frame(x = 0.4, xend = 1, y = root$yn)
  blocks <- fm %>% group_by(phy) %>% summarise(ymin = min(y) - 0.5, ymax = max(y) + 0.5, tips = n(), .groups = "drop") %>% arrange(desc(ymax)); blocks$band <- rep(c(TRUE, FALSE), length.out = nrow(blocks))
  band <- geom_rect(data = blocks[blocks$band, ], aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = "grey94")
  ylim_ <- scale_y_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)); no_y <- theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank())
  pm <- theme(plot.margin = margin(6, 4, 6, 4), panel.grid.minor = element_blank())
  pT <- ggplot() + band + geom_segment(data = vseg, aes(x = x, xend = x, y = ymin, yend = ymax, colour = phy), linewidth = 0.35) + geom_segment(data = hseg, aes(x = x, xend = xend, y = y, yend = y, colour = phy), linewidth = 0.35) +
    geom_segment(data = stub, aes(x = x, xend = xend, y = y, yend = y), colour = "grey40", linewidth = 0.35) + geom_point(data = fm, aes(5, y, colour = phy), size = 1.1) + geom_text(data = fm, aes(5.2, y, label = family), hjust = 0, size = 2.3) +
    geom_rect(data = blocks, aes(xmin = -0.95, xmax = -0.2, ymin = ymin + 0.1, ymax = ymax - 0.1, fill = phy), inherit.aes = FALSE) + geom_text(data = blocks[blocks$tips >= 4, ], aes(-0.575, (ymin + ymax) / 2, label = phy), angle = 90, size = 2.6, colour = "white", fontface = "bold") +
    scale_colour_manual(values = pal, guide = "none") + scale_fill_manual(values = pal, guide = "none") + scale_x_continuous(limits = c(-0.95, 9), expand = c(0, 0), breaks = 0, labels = "0") + ylim_ +
    labs(x = " ", y = NULL, title = "Taxonomy of the alpine specialist OTUs") + theme_s + pm + no_y + theme(plot.margin = margin(6, 4, 6, 24), axis.text.x = element_text(colour = "white"), axis.ticks.x = element_line(colour = "white"), axis.title.x = element_text(colour = "white"), panel.grid.major.x = element_blank(), panel.border = element_blank())
  pBar <- ggplot(fm) + band + geom_col(aes(n_spec, y, fill = phy), orientation = "y", width = 0.7) + geom_text(aes(n_spec, y, label = sprintf("%d (%.0f%%)%s", n_spec, 100 * share, ifelse(fdr < 0.05, " *", ""))), hjust = -0.1, size = 2.3) +
    scale_fill_manual(values = pal, guide = "none") + scale_x_continuous(expand = expansion(mult = c(0, 0.45))) + ylim_ + labs(x = "Specialist OTUs (% of the family; * enriched)", y = NULL, title = "Alpine specialist OTUs") + theme_s + pm + no_y
  list(tree = pT, bars = pBar, n = n)
}
tp <- tree_plots(fs)

## ---- 5. Panel B: space ----------------------------------------------------------------------------------------------------------------------------------
sp_rows <- match(spec$otu, rownames(Prel)); cols_all <- col_of[an_ix]
Psp <- Prel[sp_rows, cols_all, drop = FALSE]
dat_s <- data.frame(cluster = cl_all[an_ix], group = ifelse(is_alp[an_ix], "Alpine", "Lowland"), lat = tab$Latitude[an_ix], lon = tab$Longitude[an_ix], share = as.numeric(Matrix::colSums(Psp)), n_spec = as.numeric(Matrix::colSums(Psp > 0)), stringsAsFactors = FALSE)
# every sample is placed and assigned to its mountain range by ITS OWN coordinates: a cluster (samples of one BioProject or within 100 m of each other) can span several regions, so the mean position
# of a cluster is not a place (it can fall in the sea or in the lowlands) and is never used for geography.
hav <- function(lat1, lon1, lat2, lon2) { r <- pi / 180; a <- sin((lat2 - lat1) * r / 2)^2 + cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2; 2 * 6371 * asin(pmin(1, sqrt(a))) }
dat_s$bioproject <- if ("bioproject" %in% names(tab)) as.character(tab$bioproject[an_ix]) else NA_character_; dat_s$sample <- as.character(tab[[CFG$id_table]][an_ix])
ia <- dat_s$group == "Alpine"; dat_s$range <- "Lowland"; dat_s$range[ia] <- rng_tab[an_ix][ia]; dat_s$continent <- NA_character_; dat_s$continent[ia] <- continent_of(dat_s$lat[ia], dat_s$lon[ia])
spread <- dat_s[ia, ] %>% group_by(cluster) %>% summarise(samples = n(), bioproject = first(bioproject), spread_km = if (n() > 1) max(hav(median(lat), median(lon), lat, lon)) else 0, ranges = n_distinct(range[range != "Other"]), .groups = "drop") %>% arrange(desc(spread_km)) %>% as.data.frame()
write.csv(spread, out_t("alpine_cluster_spread.csv"), row.names = FALSE); wide <- spread[spread$spread_km > 100, ]
message(sprintf("Alpine clusters: %d, of which %d have samples more than 100 km from the median position of the cluster (largest: %s); a cluster = samples of one BioProject or within 100 m of each other, so its mean position is not a place: every sample is placed by its own coordinates",
                nrow(spread), nrow(wide), if (nrow(wide)) paste(sprintf("%s spans %.0f km, %d samples, %d ranges, BioProject %s", head(wide$cluster, 2), head(wide$spread_km, 2), head(wide$samples, 2), head(wide$ranges, 2), head(wide$bioproject, 2)), collapse = "; ") else "none"))
out_s <- dat_s[ia & dat_s$range == "Other", c("sample", "cluster", "bioproject", "lat", "lon")]
if ("GPT_Geotext" %in% names(tab)) out_s$geotext <- as.character(tab$GPT_Geotext[an_ix][ia & dat_s$range == "Other"])
write.csv(out_s, out_t("alpine_samples_outside_ranges.csv"), row.names = FALSE)
message(sprintf("  alpine samples that lie in no %s (none within %g km): %d of %d (%.1f%%), %d clusters; listed in tables/alpine_samples_outside_ranges.csv (coordinates, BioProject%s) to check their coordinates", if (!is.null(gmba_polys)) "mountain range" else "region", CFG$range_max_km, nrow(out_s), sum(ia), 100 * nrow(out_s) / sum(ia), length(unique(out_s$cluster)), if ("geotext" %in% names(out_s)) ", the GPT geographic text" else ""))
cl_share <- dat_s %>% group_by(cluster, group) %>% summarise(share = mean(share), n_spec = mean(n_spec), samples = n(), .groups = "drop") %>% as.data.frame()
rc_share <- dat_s[ia & dat_s$range != "Other", ] %>% group_by(range, cluster) %>% summarise(samples = n(), share = mean(share), n_spec = mean(n_spec), .groups = "drop") %>% as.data.frame(); write.csv(rc_share, out_t("range_cluster_specialist_share.csv"), row.names = FALSE)
message(sprintf("Panel B: mean share of the reads from alpine specialists: alpine clusters %.2f%% (median %.2f%%), lowland clusters %.3f%%", 100 * mean(cl_share$share[cl_share$group == "Alpine"]), 100 * median(cl_share$share[cl_share$group == "Alpine"]), 100 * mean(cl_share$share[cl_share$group == "Lowland"])))
alp_c <- sort(unique(cl_all[alp_ix])); isa_s <- which(is_alp[an_ix]); Wa <- Matrix::sparseMatrix(i = isa_s, j = match(cl_all[an_ix][isa_s], alp_c), x = 1, dims = c(length(an_ix), length(alp_c)))      # samples x alpine clusters
cont_s <- factor(dat_s$continent[isa_s]); cm <- Matrix::sparseMatrix(i = isa_s, j = as.integer(cont_s), x = 1, dims = c(length(an_ix), nlevels(cont_s)))      # samples x continents (the continent of the sample itself)
n_cont <- function(rows) { Pb <- Prel[rows, cols_all, drop = FALSE] > 0; list(clusters = Matrix::rowSums((Pb %*% Wa) > 0), continents = Matrix::rowSums((Pb %*% cm) > 0)) }
ns <- n_cont(sp_rows); oth <- otu[otu$category == "Not specialist" & otu$occ_alpine >= CFG$min_occ_alpine, ]; no <- n_cont(match(oth$otu, rownames(Prel)))
cont_tab <- bind_rows(data.frame(set = "Alpine specialists", continents = ns$continents), data.frame(set = "Other alpine OTUs", continents = no$continents)); cont_tab$continents <- pmin(cont_tab$continents, 4)
cont_sum <- cont_tab %>% count(set, continents) %>% group_by(set) %>% mutate(share = n / sum(n)) %>% ungroup() %>% as.data.frame(); write.csv(cont_sum, out_t("continent_occupancy.csv"), row.names = FALSE)
spec$alpine_clusters <- as.integer(ns$clusters); spec$continents <- as.integer(ns$continents); write.csv(spec, out_t("alpine_specialist_otus.csv"), row.names = FALSE)
message(sprintf("  continents occupied: alpine specialists median %d (%.0f%% in a single continent) vs other alpine OTUs median %d (%.0f%%); clusters occupied: specialists median %d, others median %d", median(ns$continents), 100 * mean(ns$continents == 1), median(no$continents), 100 * mean(no$continents == 1), median(ns$clusters), median(no$clusters)))
# --- one dot per mountain range, and the estimated number of specialist OTUs per sample at the same sequencing depth and share of unmapped reads
us_all <- if (all(c("unmapped_raw_counts", "total_raw_counts") %in% names(tab))) tab$unmapped_raw_counts / tab$total_raw_counts else rep(NA_real_, nrow(tab)); us_all[!is.finite(us_all)] <- median(us_all[is.finite(us_all)], na.rm = TRUE); us_all <- pmin(pmax(us_all, 0), 1)
dat_s$log10_reads <- log10(tab$alpha_total_n_reads[an_ix]); dat_s$unmapped_share <- us_all[an_ix]
lon_circ <- function(l) atan2(mean(sin(l * pi / 180)), mean(cos(l * pi / 180))) * 180 / pi
rng_all <- dat_s[ia & dat_s$range != "Other", ] %>% group_by(range) %>% summarise(clusters = n_distinct(cluster), bioprojects = n_distinct(bioproject), samples = n(), lat = median(lat), lon = lon_circ(lon), .groups = "drop") %>% left_join(rc_share %>% group_by(range) %>% summarise(mean_share_of_reads = mean(share), .groups = "drop"), by = "range") %>% as.data.frame()
rng_keep <- rng_all$range[rng_all$clusters >= CFG$min_range_clusters]
message(sprintf("Spatial units: %s: %d with alpine samples, %d with at least %d clusters (analysed): %s", region_kind, nrow(rng_all), length(rng_keep), CFG$min_range_clusters, paste(sprintf("%s (%d clusters, %d samples)", head(rng_all$range[order(-rng_all$samples)], 6), head(rng_all$clusters[order(-rng_all$samples)], 6), head(sort(rng_all$samples, decreasing = TRUE), 6)), collapse = ", ")))
cap_s <- function(ix) unlist(lapply(split(ix, paste(dat_s$cluster[ix], dat_s$range[ix])), function(i) if (length(i) > CFG$max_per_cluster) sample(i, CFG$max_per_cluster) else i), use.names = FALSE)
ia_s <- which(dat_s$group == "Alpine" & dat_s$range %in% rng_keep); il_s <- cap_s(which(dat_s$group == "Lowland")); if (length(il_s) > CFG$max_lowland_ref) il_s <- sample(il_s, CFG$max_lowland_ref)
ms <- dat_s[c(cap_s(ia_s), il_s), ]; ms$range <- factor(ifelse(ms$group == "Lowland", "Lowland", ms$range), levels = c("Lowland", rng_keep)); ms$cluster <- factor(ms$cluster)
std_reads <- median(ms$log10_reads[ms$group == "Alpine"]); std_un <- median(ms$unmapped_share[ms$group == "Alpine"])
f_rng <- as.formula(paste("n_spec ~ range + s(log10_reads, k = 4) +", if (sd(ms$unmapped_share) > 1e-6) "s(unmapped_share, k = 4) +", "s(cluster, bs = \"re\")"))
mR <- if (length(rng_keep) >= 3) tryCatch(mgcv::bam(f_rng, family = mgcv::nb(), data = ms, method = "fREML", discrete = TRUE), error = function(e) { message("  range model failed: ", conditionMessage(e)); NULL }) else NULL
rng_est <- NULL
if (!is.null(mR)) { nd <- data.frame(range = factor(levels(ms$range), levels = levels(ms$range)), log10_reads = std_reads, unmapped_share = std_un, cluster = ms$cluster[1])
  pr <- predict(mR, nd, type = "link", se.fit = TRUE, exclude = "s(cluster)", newdata.guaranteed = TRUE)
  rng_est <- data.frame(range = levels(ms$range), est = exp(pr$fit), lo = exp(pr$fit - 1.96 * pr$se.fit), hi = exp(pr$fit + 1.96 * pr$se.fit), stringsAsFactors = FALSE)
  rng_est$raw_mean <- as.numeric(tapply(ms$n_spec, ms$range, mean))[match(rng_est$range, levels(ms$range))]
  rng_est$distinct_specialists_observed <- vapply(rng_est$range, function(r) { w <- which(dat_s$range == r & dat_s$group == "Alpine"); if (r == "Lowland" || !length(w)) NA_integer_ else sum(Matrix::rowSums(Psp[, w, drop = FALSE] > 0) > 0) }, 1L)
  rng_est <- left_join(rng_est, rng_all, by = "range"); write.csv(rng_est, out_t("range_specialist_estimates.csv"), row.names = FALSE)
  ra <- rng_est[rng_est$range != "Lowland", ]; ra <- ra[order(-ra$est), ]
  message(sprintf("  estimated alpine specialist OTUs per sample at the same depth (median %.0f reads) and unmapped share (%.1f%%): %s; lowland reference %.2f", 10^std_reads, 100 * std_un, paste(sprintf("%s %.1f", head(ra$range, 4), head(ra$est, 4)), collapse = ", "), rng_est$est[rng_est$range == "Lowland"]))
  message(sprintf("  model of the number of specialist OTUs per sample: depth smooth edf %.1f, unmapped-share smooth edf %.1f", summary(mR)$s.table["s(log10_reads)", "edf"], if ("s(unmapped_share)" %in% rownames(summary(mR)$s.table)) summary(mR)$s.table["s(unmapped_share)", "edf"] else NA)) }
wm <- tryCatch(ggplot2::map_data("world"), error = function(e) NULL)
rd <- if (!is.null(rng_est)) rng_est[rng_est$range != "Lowland", ] else rng_all[rng_all$range %in% rng_keep, ]; if (is.null(rng_est)) rd$est <- rd$mean_share_of_reads
rd_few <- rng_all[!rng_all$range %in% rng_keep, ]                                                                    # ranges with too few clusters for an estimate: a grey open circle
pMap <- ggplot() + (if (!is.null(wm)) geom_polygon(data = wm, aes(long, lat, group = group), fill = "grey92", colour = NA)) + geom_point(data = rd_few, aes(lon, lat), shape = 1, colour = "grey35", size = 2.2, stroke = 0.7) + geom_point(data = rd[order(rd$est), ], aes(lon, lat, colour = est, size = clusters), alpha = 0.9) +
  geom_text(data = transform(head(rd[order(-rd$est), ], CFG$n_range_labels), hj = ifelse(lon > 140, 1, ifelse(lon < -140, 0, 0.5))), aes(lon, lat, label = sub("\\s*\\(.*$", "", range), hjust = hj), vjust = -1.4, size = 2.3, check_overlap = TRUE, show.legend = FALSE, colour = "grey20", inherit.aes = FALSE) +
  scale_colour_viridis_c(option = "magma", direction = -1, name = if (is.null(rng_est)) "Share of reads" else "Estimated specialist\nOTUs per sample") + scale_size_continuous(range = c(2.2, 7), name = "Clusters") +
  coord_fixed(1.25, xlim = c(-180, 180), ylim = c(-60, 80), expand = FALSE) + labs(x = NULL, y = NULL, title = sprintf("Where: one dot per mountain range (%d ranges; grey circle: %d with too few clusters)", nrow(rd) + nrow(rd_few), nrow(rd_few))) + theme_s + theme(legend.position = "bottom", legend.key.width = unit(0.9, "cm"), legend.title = element_text(size = 8), legend.text = element_text(size = 7.5), legend.box = "horizontal")
pRng <- ggplot() + theme_void() + labs(title = "Estimated number per mountain range - not available")
if (!is.null(rng_est)) { rp <- rng_est[order(rng_est$est), ]; rp <- tail(rp, 31); if (!"Lowland" %in% rp$range) rp <- rbind(rng_est[rng_est$range == "Lowland", ], rp); rp <- rp[order(rp$est), ]; rp$lab <- ifelse(is.na(rp$clusters), as.character(rp$range), sprintf("%s (%d; %.0f%% of reads)", rp$range, rp$clusters, 100 * rp$mean_share_of_reads)); rp$lab <- factor(rp$lab, levels = rp$lab); rp$kind <- ifelse(rp$range == "Lowland", "Lowland", "Alpine")
  pRng <- ggplot(rp, aes(y = lab)) + geom_linerange(aes(xmin = lo, xmax = hi, colour = kind), linewidth = 0.5) + geom_point(aes(x = raw_mean), shape = 1, size = 1.8, colour = "grey35") + geom_point(aes(x = est, colour = kind), size = 2.2) +
    scale_colour_manual(values = c(Alpine = "#C0392B", Lowland = "#7A9A7E"), guide = "none") + scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
    labs(x = "Specialist OTUs per sample\n(filled: same depth and unmapped\nshare, 95% interval; open: raw)", y = NULL, title = "How many per mountain range?") + theme_s + theme(axis.text.y = element_text(size = 6.5), axis.title.x = element_text(size = 8.5), plot.title.position = "plot") }
pCont <- ggplot(cont_sum, aes(factor(continents, levels = 1:4, labels = c("1", "2", "3", "4+")), share, fill = set)) + geom_col(position = position_dodge(width = 0.8), width = 0.75, colour = "grey30", linewidth = 0.2) +
  scale_fill_manual(values = c("Alpine specialists" = "#C0392B", "Other alpine OTUs" = "grey70"), name = NULL) + scale_y_continuous(labels = scales::percent) + labs(x = "Continents occupied", y = "Share of the OTUs", title = "How widespread?") + theme_s + theme(legend.position = "bottom", legend.text = element_text(size = 7), plot.title.position = "plot") + guides(fill = guide_legend(nrow = 2))

## ---- 6. Panel C: CLR abundance of the specialists along the season and the soil depth -----------------------------------------------------------------------
ds <- trimws(as.character(if ("Date_collected" %in% names(tab)) tab$Date_collected else env$sample_date))
prec <- ifelse(grepl("^\\d{4}-\\d{2}-\\d{2}$", ds), "day", ifelse(grepl("^\\d{4}-\\d{2}$", ds), "month", "other")); dte <- suppressWarnings(as.Date(ifelse(prec == "day", ds, ifelse(prec == "month", paste0(ds, "-15"), NA)), format = "%Y-%m-%d"))
doy <- as.integer(format(dte, "%j")); doy_al <- ifelse(tab$Latitude < 0, (doy + 182) %% 365, doy)
parse_depth <- function(x) {
  x0 <- tolower(trimws(as.character(x))); x0[x0 %in% c("", "na", "n/a", "nan", "none", "null", "unknown", "not applicable", "missing", "not collected", "not provided")] <- NA
  val <- rep(NA_real_, length(x0)); un <- rep(NA_character_, length(x0)); ok <- which(!is.na(x0)); x1 <- gsub("(?<=[0-9]),(?=[0-9])", ".", x0[ok], perl = TRUE)
  nums <- regmatches(x1, gregexpr("[0-9]*\\.?[0-9]+", x1)); val[ok] <- vapply(nums, function(z) if (length(z) >= 2) mean(as.numeric(z[1:2])) else if (length(z) == 1) as.numeric(z) else NA_real_, 1)
  un[ok] <- ifelse(grepl("mm", x1), "mm", ifelse(grepl("cm", x1), "cm", ifelse(grepl("(^|[^a-z])m([^a-z]|$)|meter|metre", x1), "m", NA_character_)))
  list(raw = x0, value = val, unit = un) }
if ("Depth_m" %in% names(tab)) { pd <- parse_depth(tab$Depth_m); no_u <- is.na(pd$unit) & is.finite(pd$value)
  u_def <- if (CFG$depth_unit == "auto") { q <- suppressWarnings(quantile(pd$value[no_u], 0.95, na.rm = TRUE)); if (is.finite(q) && q > 3) "cm" else "m" } else CFG$depth_unit
  depth_cm <- ifelse(is.na(pd$unit), pd$value * (if (u_def == "m") 100 else 1), ifelse(pd$unit == "cm", pd$value, ifelse(pd$unit == "mm", pd$value / 10, pd$value * 100))); depth_cm[!is.finite(depth_cm) | depth_cm < 0] <- NA
  bad <- unique(pd$raw[is_alp & !is.na(pd$raw) & !is.finite(depth_cm)])
  message(sprintf("Soil depth (Depth_m): %d of %d alpine samples have a non-empty value, %d were read as a depth (values without a unit taken as %s); median %.0f cm, range %.0f-%.0f cm, %d distinct values%s",
                  sum(is_alp & !is.na(pd$raw)), sum(is_alp), sum(is_alp & is.finite(depth_cm)), u_def, median(depth_cm[is_alp], na.rm = TRUE), suppressWarnings(min(depth_cm[is_alp], na.rm = TRUE)), suppressWarnings(max(depth_cm[is_alp], na.rm = TRUE)), length(unique(depth_cm[is_alp & is.finite(depth_cm)])),
                  if (length(bad)) paste0("; NOT understood, e.g.: ", paste(head(bad, 6), collapse = " | ")) else ""))
} else { depth_cm <- rep(NA_real_, nrow(tab)); message("Soil depth: the table has no Depth_m column") }
us <- if (all(c("unmapped_raw_counts", "total_raw_counts") %in% names(tab))) tab$unmapped_raw_counts / tab$total_raw_counts else rep(NA_real_, nrow(tab)); us[!is.finite(us)] <- median(us[is.finite(us)], na.rm = TRUE); us <- pmin(pmax(us, 0), 1)
message(sprintf("Panel C: %d alpine samples with a month / day, %d with a soil depth (median %.0f cm)", sum(is_alp & !is.na(doy_al)), sum(is_alp & !is.na(depth_cm)), median(depth_cm[is_alp], na.rm = TRUE)))
al <- data.frame(idx = alp_ix, cluster = factor(cl_all[alp_ix]), lat = tab$Latitude[alp_ix], lon = tab$Longitude[alp_ix], log10_reads = log10(tab$alpha_total_n_reads[alp_ix]), unmapped_share = us[alp_ix], doy_al = doy_al[alp_ix], log10_depth = log10(depth_cm[alp_ix] + 1), stringsAsFactors = FALSE)
cap_ix <- unlist(lapply(split(seq_len(nrow(al)), al$cluster), function(i) if (length(i) > CFG$max_per_cluster) sample(i, CFG$max_per_cluster) else i), use.names = FALSE); al <- al[sort(cap_ix), ]
# the richest phyla of the specialists, ties in alphabetical order
sp_ph <- ifelse(nzchar(spec$phylum), spec$phylum, "(no phylum)"); tb_ <- table(sp_ph); top_ph <- head(names(tb_)[order(-as.integer(tb_), names(tb_))], CFG$n_parts); part_of <- ifelse(sp_ph %in% top_ph, sp_ph, "Other phyla")
parts <- c(top_ph, if (any(part_of == "Other phyla")) "Other phyla"); cnt <- Xcnt[, col_of[al$idx], drop = FALSE]; totc <- Matrix::colSums(cnt)
Cs <- do.call(rbind, lapply(parts, function(p) Matrix::colSums(cnt[match(spec$otu[part_of == p], rownames(cnt)), , drop = FALSE]))); rownames(Cs) <- parts; rest <- totc - Matrix::colSums(Cs)
L <- log(rbind(Cs, `Non-specialists` = rest) + CFG$pseudocount); clr <- sweep(L, 2, colMeans(L)) / log(2); for (p in parts) al[[paste0("clr_", p)]] <- clr[p, ]
sall <- Matrix::colSums(Cs); al$`clr_All specialists` <- 0.5 * log2((sall + CFG$pseudocount) / (rest + CFG$pseudocount))
k_sp <- max(5, min(CFG$k_space, nrow(unique(al[, c("lat", "lon")])) - 2)); base_terms <- c("s(log10_reads, k = 4)", if (sd(al$unmapped_share) > 1e-6) "s(unmapped_share, k = 4)", sprintf("s(lat, lon, bs = \"sos\", k = %d)", k_sp), "s(cluster, bs = \"re\")")
all_parts <- c("All specialists", parts); mods <- list(); curves <- list()
fit_one <- function(part, xvar, term, d, grid, knots = NULL) { d$cluster <- droplevels(d$cluster); d$y <- d[[paste0("clr_", part)]]
  m <- tryCatch(mgcv::bam(as.formula(paste("y ~", paste(c(term, base_terms), collapse = " + "))), data = d, method = "fREML", discrete = TRUE, knots = knots), error = function(e) NULL); if (is.null(m)) return(NULL)
  nd <- d[rep(1, length(grid)), ]; for (v in names(d)) if (is.numeric(d[[v]])) nd[[v]] <- median(d[[v]], na.rm = TRUE); nd[[xvar]] <- grid; nd$cluster <- d$cluster[1]
  pr <- tryCatch(predict(m, nd, type = "terms", se.fit = TRUE, discrete = FALSE), error = function(e) NULL); tn <- paste0("s(", xvar, ")"); st <- summary(m)$s.table; if (is.null(pr) || !tn %in% colnames(pr$fit)) return(NULL)
  list(row = data.frame(part = part, variable = xvar, edf = st[tn, "edf"], F = st[tn, "F"], p = st[tn, "p-value"], n = nrow(d), stringsAsFactors = FALSE), curve = data.frame(part = part, variable = xvar, x = grid, effect = pr$fit[, tn], se = pr$se.fit[, tn]))
}
ds_ <- al[!is.na(al$doy_al), ]; dd_ <- al[is.finite(al$log10_depth), ]
for (p in all_parts) {
  if (nrow(ds_) >= CFG$min_n) { o <- fit_one(p, "doy_al", "s(doy_al, bs = \"cc\", k = 6)", ds_, seq(0, 364, length.out = 80), list(doy_al = c(0, 365))); if (!is.null(o)) { mods[[length(mods) + 1]] <- o$row; curves[[length(curves) + 1]] <- o$curve } }
  if (nrow(dd_) >= CFG$min_depth_n && length(unique(dd_$log10_depth)) >= 3) { o <- fit_one(p, "log10_depth", "s(log10_depth, k = 4)", dd_, seq(min(dd_$log10_depth), max(dd_$log10_depth), length.out = 60)); if (!is.null(o)) { mods[[length(mods) + 1]] <- o$row; curves[[length(curves) + 1]] <- o$curve } } }
mods <- bind_rows(mods); curves <- bind_rows(curves); if (nrow(mods)) { write.csv(mods, out_t("specialist_clr_models.csv"), row.names = FALSE)
  for (i in seq_len(nrow(mods))) message(sprintf("  %-18s %-12s edf %.1f, p %s (n = %d)", mods$part[i], mods$variable[i], mods$edf[i], format.pval(mods$p[i], digits = 2), mods$n[i])) } else message("  panel C: too few samples with a date / depth")
pcol <- c(pal, `Other phyla` = "grey55", `All specialists` = "black"); leg_levels <- c("All specialists", parts)       # the same colours and order in the season and the depth panels (and the phylum colours of panel A)
cplot <- function(vr, xlab, title, xlabs = NULL, n_have = NA, n_need = NA) { cc <- curves[curves$variable == vr, ]; if (!nrow(cc)) return(ggplot() + theme_void() + labs(title = sprintf("%s: %d usable alpine samples (at least %d needed)", title, as.integer(n_have), as.integer(n_need)))); cc$part <- factor(cc$part, levels = leg_levels)
  g <- ggplot(cc, aes(x, effect, colour = part, fill = part)) + geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) + geom_ribbon(aes(ymin = effect - 1.96 * se, ymax = effect + 1.96 * se), alpha = 0.12, colour = NA) + geom_line(aes(linewidth = part == "All specialists")) +
    scale_linewidth_manual(values = c(`TRUE` = 1.3, `FALSE` = 0.7), guide = "none") + scale_colour_manual(values = pcol, limits = leg_levels, guide = "none") + scale_fill_manual(values = pcol, limits = leg_levels, guide = "none") + labs(x = xlab, y = "Partial effect on CLR abundance", title = title) + theme_s + theme(legend.position = "none")
  if (!is.null(xlabs)) g <- g + scale_x_continuous(breaks = log10(xlabs + 1), labels = xlabs); g }
pSeas <- cplot("doy_al", "Day of the year (southern hemisphere shifted by half a year)", "Season", n_have = nrow(ds_), n_need = CFG$min_n); pDep <- cplot("log10_depth", "Soil depth (cm; log10(x + 1) scale)", "Soil depth", c(0, 2, 5, 10, 20, 30, 50, 100), n_have = nrow(dd_), n_need = CFG$min_depth_n)

# ONE legend for the season and the depth panels: built from the fixed palette, independent of the data of either panel
leg_plot <- ggplot(data.frame(part = factor(leg_levels, levels = leg_levels), x = 1, y = 1), aes(x, y, colour = part, linewidth = part == "All specialists")) + geom_line() + scale_colour_manual(values = pcol, limits = leg_levels, name = NULL) +
  scale_linewidth_manual(values = c(`TRUE` = 1.3, `FALSE` = 0.7), guide = "none") + theme_s + theme(legend.position = "bottom", legend.text = element_text(size = 7.5), legend.margin = margin(0, 0, 0, 0)) + guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = c(1.3, rep(0.7, length(leg_levels) - 1)))))
leg_grob <- tryCatch({ g <- ggplotGrob(leg_plot); i <- grep("^guide-box", vapply(g$grobs, function(x) x$name, "")); i <- i[vapply(g$grobs[i], function(x) !inherits(x, "zeroGrob"), TRUE)]; if (length(i)) g$grobs[[i[1]]] else NULL }, error = function(e) NULL)

## ---- 7. Figure -----------------------------------------------------------------------------------------------------------------------------------------------
theme_s <- theme_s
cap <- sprintf(paste0("Figure 6. Alpine specialist OTUs (%d of %d tested OTUs). Specialist: the niche breadth in the pooled alpine + lowland environment (abundance-weighted variance of z-scored climate, soil and terrain principal components of the units) is narrower than a permutation null (z <= -%.2f and at least %.0f%% below the null mean) and at least %.0f%% of the abundance is in alpine units, above chance (random draws of one sample per unit, a cluster split by mountain range or 5-degree cell, %d draws, calls in at least %.0f%% of them). ",
                      "A: taxonomic tree of the %d families with the most alpine specialist OTUs (bars: number of specialist OTUs, in brackets their share of the OTUs of the family; * enriched, Fisher test, FDR < 0.05; colours: phyla). ",
                      "B: one dot per mountain range (colour: the estimated number of alpine specialist OTUs per sample at the same sequencing depth and share of unmapped reads; size: number of clusters; grey open circle: a range with too few clusters, i.e. independent studies, for an estimate) and the number of continents where the OTUs occur (alpine specialists vs the other alpine OTUs). ",
                      "C: partial effect of the day of the year (cyclic) and of the soil depth on the CLR abundance of the specialists of each phylum (parts of the composition: specialists of the phyla shown, other specialists, all other reads; the same colours as the tree of A, one legend for both panels) and of all specialists together (2-part log ratio; another scale), GAMs on alpine samples with log10 depth, unmapped share, coordinates and a random cluster effect. B, right: the estimated number of alpine specialist OTUs per sample for every mountain range with at least a few clusters (number of clusters and the mean share of the reads from specialists over all samples in brackets; filled: expected number at the median alpine sequencing depth and share of unmapped reads from a negative binomial GAM with a random cluster effect, 95%% interval; open: raw mean; green: lowland reference)."),
               nrow(spec), nrow(otu), CFG$z_cut, 100 * CFG$min_effect, 100 * CFG$alpine_share_min, CFG$n_draws, 100 * CFG$min_draw_frac, nrow(fs))
cl_ <- strwrap(cap, 250); H <- max(10, 0.15 * tp$n + 3.2)
# the right block: 12 equal columns; top row: the map (6), the continents (2) and the mountain ranges (4); middle row: season (6) and soil depth (6); bottom row: their shared legend
NEST <- function() grid.layout(3, 12, widths = unit(rep(1, 12), "null"), heights = unit.c(unit(1.12, "null"), unit(1, "null"), unit(0.42, "in")))
png(file.path(out_folder, "fig6_alpine_specialist_otus.png"), width = 22, height = H + 0.17 * length(cl_), units = "in", res = 300)
grid.newpage(); pushViewport(viewport(layout = grid.layout(2, 3, widths = unit(c(5.0, 2.4, 11.5), "null"), heights = unit.c(unit(1, "null"), unit(0.17 * length(cl_) + 0.2, "in")))))
gs <- tryCatch({ g <- lapply(list(tp$tree, tp$bars), ggplotGrob); mh <- do.call(grid::unit.pmax, lapply(g, function(z) z$heights)); lapply(g, function(z) { z$heights <- mh; z }) }, error = function(e) NULL)
for (i in 1:2) { if (!is.null(gs)) { pushViewport(viewport(layout.pos.row = 1, layout.pos.col = i)); grid.draw(gs[[i]]); popViewport() } else print(list(tp$tree, tp$bars)[[i]], vp = viewport(layout.pos.row = 1, layout.pos.col = i)) }
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3)); pushViewport(viewport(layout = NEST()))
print(pMap, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:6)); print(pCont, vp = viewport(layout.pos.row = 1, layout.pos.col = 7:8)); print(pRng, vp = viewport(layout.pos.row = 1, layout.pos.col = 9:12)); print(pSeas, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:6)); print(pDep, vp = viewport(layout.pos.row = 2, layout.pos.col = 7:12)); if (!is.null(leg_grob)) { pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1:12)); grid.draw(leg_grob); popViewport() }
popViewport(2)
tag <- function(txt, vp) grid.text(txt, x = unit(0.05, "in"), y = unit(1, "npc") - unit(0.12, "in"), just = c("left", "top"), gp = gpar(fontsize = 17, fontface = "bold"), vp = vp)
tag("A", viewport(layout.pos.row = 1, layout.pos.col = 1)); pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3)); pushViewport(viewport(layout = NEST()))
tag("B", viewport(layout.pos.row = 1, layout.pos.col = 1)); tag("C", viewport(layout.pos.row = 2, layout.pos.col = 1)); popViewport(2)
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1:3)); grid.text(paste(cl_, collapse = "\n"), x = 0.01, just = "left", gp = gpar(fontsize = 9)); popViewport(2); dev.off()
message("Done. Outputs in ", out_folder)
writeLines(capture.output(sessionInfo()), file.path(out_folder, "sessionInfo.txt"))
file.create(file.path(out_folder, ".done"))   # flag file for Snakemake