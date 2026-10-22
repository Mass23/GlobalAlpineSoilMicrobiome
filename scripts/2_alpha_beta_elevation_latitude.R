#!/usr/bin/env Rscript
# =============================================================================
# 2_alpha_beta_elevation_latitude.R  --  Figure 2
# Alpha AND beta diversity along ELEVATION and LATITUDE, alpine vs lowland.
#
# PART A  alpha diversity (Observed OTUs, Shannon, Pielou J) and the UNMAPPED READ PROPORTION
#   vs elevation and latitude (the unmapped proportion is modelled as a response too: Supplementary
#   Figure 2; there the unmapped terms are not used as covariates), and Supplementary Figure 3:
#   spatial surfaces of the unmapped proportion (where the reference is less complete), lowland vs alpine
#   Three GAM specifications per metric:
#     net : alpine + s(log10 reads) + s(unmapped proportion) + s(unmapped, by = alpine-ordered)
#           + s(elevation) + s(elevation, by = alpine-ordered)   [lowland = reference,
#           + s(latitude)  + s(latitude,  by = alpine-ordered)    by-smooths = alpine - lowland]
#           + sos spatial spline   -> MAIN: gradients net of space, sequencing depth and the
#                                     unmapped proportion (script 1: it differs between groups
#                                     and its relation to diversity is non-linear)
#     raw : the same WITHOUT the spatial spline -> raw gradients
#     net_nounm : like net but WITHOUT the unmapped terms (tables only, to see what they change)
#   The ordered-factor difference smooths test whether each gradient differs in shape
#   between alpine and lowland (alpha/group_difference_tests.csv, printed on the figure).
#   NB: the sos spline absorbs much of the latitudinal signal (concurvity).
#
# PART B  beta diversity (community dissimilarity) between pairs of samples
#   * SAMPLES. The ALPHA models (Figure 2A, B) use all samples: they condition on elevation, latitude and space, so an unbalanced
#     design does not bias the curves. Samples of one BioProject, or at one location, are not independent: the column cluster_id
#     (microbeatlas_filter_data.R) links such samples (a connected group = one cluster), and the alpha models use ONE random intercept per
#     cluster (cluster_re), plus a spherical spline for the regional trend; the lowland group is first thinned to one sample per location
#     (thin_lowland: a random effect with ~100k levels is not feasible in mgcv).
#   * BETA DIVERSITY compares the alpine samples with the same number of evenly spread non-alpine samples, in each of n_rep replicates:
#     every replicate takes ONE random alpine sample per cluster (a different one each time, so no study or place is over-represented
#     and no random effects are needed) and, from the draw rep_r of define_lowland_sets.py (random draws evenly spaced over latitude /
#     longitude, flags rep_01 ... in the earth-data table), the same number of non-alpine samples, at most one per thin_m metres and
#     cut down evenly over latitude / longitude cells (even_cell_deg). Up to n_pairs pairs per pair type (within lowland, within alpine,
#     alpine vs lowland) are drawn at random within geographic distance classes (dist_breaks_km) crossed with elevation-difference
#     classes (elev_breaks_km: short distances and large elevation differences would hardly be represented by random pairs), shared as
#     evenly as possible over the classes, from all samples together (no latitude zones). The same replicates are meant for the analyses
#     that need soil data (Figures 3-5): analyse replicate by replicate and combine.
#   * Bray-Curtis (relative abundances) and Jaccard (presence/absence after rarefaction to rare_depth reads), computed
#     only for the drawn pairs from sparse matrices (checked against vegan::vegdist in TEST_MODE).
#   * ONE ordinary LINEAR REGRESSION per metric and model on the dissimilarity itself (effects in units of dissimilarity) and on the pairs of
#     ALL replicates: an intercept and one slope per pair type, linear nuisance terms for the differences in depth and unmapped
#     proportion. (A beta regression does not converge with the pile of dissimilarities at exactly 1.)
#       elev: elevation difference + log10 geographic distance      -> elevation ACCOUNTING FOR SPACE (distance is a covariate, not drawn)
#       lat : latitude difference + elevation difference            -> latitude ACCOUNTING FOR ELEVATION
#     The lines of Figure 2C are this regression, the ribbons its 95% interval, the open circles the observed dissimilarity (adjusted to
#     the typical values of the other variables) in 8 bins. The stars in the panels test the slope of every pair type and the differences
#     between the slopes of two pair types.
#   * STATISTICS. Pairs that share a sample, and samples of one cluster, are not independent: an ordinary regression treats thousands of pairs as
#     thousands of observations and gives far too small p-values (simulation: 63-67% false positives at p < 0.05 for a true null). The
#     estimates are the ordinary least-squares ones; their standard errors are cluster-robust for DYADIC data (Aronow, Samii & Assenova
#     2015: the scores of all pairs that share a CLUSTER are summed), so the number of independent units is the number of clusters
#     (BioProject / location groups, column cluster_id), however many replicates there are: in a simulation this had the correct 5% false
#     positives and about twice the power of testing each replicate separately and averaging; more than ~20 replicates add nothing (they redraw
#     from the same clusters). t-tests with (clusters - 1) degrees of freedom; no adjustment for multiple tests.
#     Tables: beta/model_coefficients.csv (all coefficients, with the p-values of the ordinary regression in the column p_naive_lm_NOT_VALID,
#     for comparison only), beta/linear_effects.csv (slopes per pair type), beta/linear_effect_type_contrasts.csv (differences between pair
#     types), samples_per_replicate.csv, pair_counts.csv, share_of_pairs_at_dissimilarity_1.csv and model_checks.csv (does every fit reproduce
#     the observed mean dissimilarity? a warning is printed if not).
#   CAVEAT: the tests cannot see spatial autocorrelation between DIFFERENT clusters (neighbouring regions resemble each other): judge
#   robustness by the size of the effects.
#
# ALPHA GRADIENT AXES (gradient_axes): "raw" (default) = elevation and latitude; the sensitivity analysis "normalised" (elevation and latitude as
#   residuals against a lat/lon surface / an elevation relation) reruns the alpha part only, into <out_folder>_axes_normalised.
#   alpha/concurvity.csv says how well elevation, latitude and the spatial spline can be told apart.
#
# FIGURES (in <out_folder>)
#   fig2_alpha_beta_gradients.png : A = alpha diversity vs elevation (3 metrics), B = vs latitude (signed),
#                                   C = Bray-Curtis vs elevation difference (accounting for space) and vs latitude
#                                   difference (accounting for elevation): effects per pair type
#   figS1_beta_jaccard.png        : Supplementary Figure 1 = same as C with Jaccard
#   alpha/ and beta/              : full grids, tables, diagnostics
#
# Input : <data_folder>/sample_data_filtered.tsv, elevation from env_file and, for part B,
#         the filtered BIOM written by microbeatlas_filter_data.R
# Output: <out_folder>/alpha/ and <out_folder>/beta/
# Usage : Rscript 2_alpha_beta_elevation_latitude.R [sample_table] [env_table] [filtered_biom] [out_folder]
#         Snakemake passes the exact paths (all file names are defined in the Snakefile only)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mgcv)
  library(Matrix)
  library(rbiom)
  library(geosphere)
  library(vegan)
  library(grid)
})

## ---- 0. Switches, folders, configuration ------------------------------------------
RUN_ALPHA <- TRUE    # part A (alpha diversity GAMs)
RUN_BETA  <- TRUE    # part B (beta diversity; slow: reads the BIOM)
TEST_MODE <- FALSE   # TRUE: tiny beta run (few samples, 1 replicate) to check that it works
REUSE_FITS <- TRUE   # reuse the saved alpha models (alpha/fits.rds, alpha/figS3_model.rds) and the saved beta pairs (beta/pairs_all.rds) if they match the data and settings

args <- commandArgs(trailingOnly = TRUE)
# Snakemake passes: args[1] = sample table (TSV), args[2] = all_sampled_data.csv, args[3] = filtered BIOM,
# args[4] = output folder. The defaults are for running by hand.
data_folder     <- "/home/renku/work/epfl-altshuler/GlobalAlpineSoilMicrobiome/data"
results_default <- file.path(dirname(data_folder), "results")
sample_table  <- if (length(args) >= 1) as.character(args[1]) else file.path(data_folder, "sample_data_filtered.tsv")
env_table     <- if (length(args) >= 2) as.character(args[2]) else file.path(results_default, "all_sampled_data.csv")
filtered_biom <- if (length(args) >= 3) as.character(args[3]) else
  file.path(data_folder, "microbeatlas", "otuTable_soilCurated_OTU97_nextMAPrelease_filtered.biom")
out_folder    <- if (length(args) >= 4) as.character(args[4]) else file.path(results_default, "2_alpha_beta_elevation_latitude")

CFG <- list(
  # --- input ---
  alpha_file  = sample_table,
  env_file    = env_table,      # all_sampled_data.csv: elevation (DEM) etc., written by download_earth_data
  alpha_id    = "Sample_acc",   # sample ID column in sample_data_filtered.tsv   (check!)
  env_id      = "site",         # matching ID column in env_file                  (check!)
  env_elev    = "elevation_m",  # elevation column (metres) in env_file: written by download_others.py
  biom_file   = filtered_biom,  # written by the Snakemake rule filter_biom
  biom_id     = "Accessions",   # column of alpha_file that matches the BIOM sample IDs
  # --- part A ---
  abs_lat     = FALSE,          # alpha panels B: FALSE = signed latitude, TRUE = |latitude|
  dist_breaks_km = c(0, 10, 100, 1000, 5000, Inf),   # pairs are drawn at random WITHIN these geographic distance classes (km) ...
  elev_breaks_km = c(0, 0.1, 0.5, 1, 2, Inf),        # ... crossed with these classes of the elevation difference (km): every elevation
                                                     #   difference is represented, also among the lowland pairs (equal numbers per class)
  thin_m      = 100,            # a LOCATION = samples within this radius (metres). Beta subsamples: at most ONE sample per location and replicate
                                #   (metres, chosen at random); 0 = only distinct coordinates (exact duplicates removed)
  cluster_re    = TRUE,           # ONE random intercept per cluster (column cluster_id from microbeatlas_filter_data.R: samples that share a
                                  #   BioProject or lie within thin_m metres of each other, directly or through a chain of such samples)
  min_cluster_n = 5,              # clusters with fewer samples (in the model data) share ONE level, as do samples without a cluster_id
  thin_alpine   = FALSE,          # thin a group to one random sample per location BEFORE modelling (cheaper, but discards samples): alpine ...
  thin_lowland  = TRUE,           # ... and lowland (data-rich: thinning costs little)
  max_re_levels = 3000,           # safety: the random effect may have at most this many levels (the cost grows steeply with the levels)
  gradient_axes = "raw",        # x axes of the ALPHA gradients (Figure 2A, B): "raw" = elevation and latitude; "normalised" = elevation as the residual
                                #   of a lat/lon surface, latitude as the residual of an elevation relation (a SENSITIVITY analysis: alpha part
                                #   only, results in <out_folder>_axes_normalised)
  k_norm        = 30,           # basis dimension of the lat/lon surface (elevation baseline of "normalised")
  concurvity    = TRUE,         # table of how well elevation, latitude and the spatial spline can be predicted from each other (alpha/concurvity.csv)
  concurvity_n  = 15000,        # rows of the random sample on which the concurvity is computed
  k_space     = 60,             # basis dimension of the sos spatial spline
  k_space_alpine = 30,          # Supplementary Figure 3: basis dimension of the alpine-specific surface
  grid_res    = 2,              # Supplementary Figure 3: grid resolution (degrees)
  min_support = 5,              # Supplementary Figure 3: min. samples within the 3x3 cell neighbourhood to show a cell
  map_ylim    = c(-60, 83),     # Supplementary Figure 3: latitude window of the maps
  n_rug       = 3000,           # rug marks drawn per group (plot only)
  moran_n_max = 1500, moran_nperm = 199,
  # --- part B ---
  n_alpine_max = Inf,           # alpine samples (= clusters) per replicate: one per cluster, ALL clusters (Inf); fewer only in tests
  even_cell_deg = 5,            # the non-alpine samples are cut down evenly over latitude / longitude cells of this size (degrees)
  n_rep       = 20,             # replicates: each draws one alpine sample per cluster and the non-alpine samples again from scratch
  use_replicate_sets = TRUE,    # non-alpine samples of replicate r = those flagged rep_<r> by define_lowland_sets.py (evenly spaced over lat/lon);
                                #   FALSE (or no rep_ columns in the earth-data table): a random draw from all non-alpine samples
  rare_depth  = 1000,           # rarefaction depth for Jaccard (samples below are dropped)
  n_pairs     = 1000,           # pairs per pair type and replicate (shared as evenly as possible over the distance x elevation classes)
  min_pairs   = 100,            # pair types with fewer pairs are dropped from the models
  k_smooth    = 5,              # basis dimension of every 1-D smooth (alpha and beta models): small = smooth curves
  bs_smooth   = "ts",           # thin plate spline WITH shrinkage (the null space is penalised too)
  select      = FALSE,          # TRUE would add a second shrinkage penalty on top of "ts"
  beta_points_max = 60000,      # beta panels: number of pairs drawn as points (random subset; 0 = none)
  beta_point_alpha = 0.06,      # transparency of those points
  beta_zoom   = TRUE,           # beta panels: zoom the y axis on the bulk of the pairs (the few far below are not drawn; FALSE: 0-1)
  beta_zoom_q = 0.05,           # lower edge of that window = this quantile of the dissimilarities (or the lowest ribbon edge, if lower)
  # --- main figure ---
  main_alpha_metrics = c("richness", "shannon", "pielou"),   # alpha metrics shown in panels A, B
  # --- common ---
  nthreads    = 4,
  cols        = c(Lowland = "#7C9885", Alpine = "#2B6CB0"),
  cols_beta   = c("Within lowland" = "#7C9885", "Within alpine" = "#2B6CB0",
                  "Alpine vs lowland" = "#D98E04")
)
if (TEST_MODE) {
  CFG$n_alpine_max <- 150; CFG$n_rep <- 1; CFG$n_pairs <- 300
  out_folder <- paste0(out_folder, "_TEST")
}
if (!CFG$gradient_axes %in% c("raw", "normalised")) stop('CFG$gradient_axes must be "raw" or "normalised"')
if (CFG$gradient_axes != "raw") {                      # a sensitivity analysis of the alpha gradients: its own folder, no beta part
  out_folder <- paste0(out_folder, "_axes_", CFG$gradient_axes)
  message("gradient_axes = '", CFG$gradient_axes, "': sensitivity analysis of the alpha gradients only (the beta part is not run); results in ", out_folder)
  RUN_BETA <- FALSE
}
dir.create(file.path(out_folder, "alpha"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_folder, "beta"),  recursive = TRUE, showWarnings = FALSE)
out_a <- function(...) file.path(out_folder, "alpha", ...)
out_b <- function(...) file.path(out_folder, "beta", ...)
set.seed(1)

# term of a 1-D smooth (used by the alpha AND the beta models): thin plate spline with shrinkage, small k
# Spatial thinning: keep one randomly chosen sample per `radius_m` metres (greedy, in random order; radius 0 = distinct
# coordinates only). Returns the indices kept (at most n_max, in the order drawn). Exact great-circle neighbourhoods (3D
# coordinates on the globe, a hash of cells of the radius size), so it also works for hundreds of thousands of samples.
thin_sites <- function(lat, lon, radius_m, n_max = Inf) {
  n <- length(lat); ord <- sample.int(n); r <- max(radius_m, 0); cell <- max(r, 1)
  R <- 6371008.8
  la <- lat * pi / 180; lo <- lon * pi / 180
  X <- R * cos(la) * cos(lo); Y <- R * cos(la) * sin(lo); Z <- R * sin(la)
  cx <- floor(X / cell); cy <- floor(Y / cell); cz <- floor(Z / cell)
  off <- 2^17; key <- function(a, b, c) sprintf("%.0f", (a + off) + 2^19 * (b + off) + 2^38 * (c + off))
  nb <- as.matrix(expand.grid(-1:1, -1:1, -1:1))          # the 27 neighbouring cells
  seen <- new.env(hash = TRUE, parent = emptyenv()); keep <- integer(0)
  for (i in ord) {
    js <- unlist(mget(key(cx[i] + nb[, 1], cy[i] + nb[, 2], cz[i] + nb[, 3]), envir = seen, ifnotfound = list(NULL)), use.names = FALSE)
    if (length(js) && any((X[js] - X[i])^2 + (Y[js] - Y[i])^2 + (Z[js] - Z[i])^2 <= r^2)) next
    k <- key(cx[i], cy[i], cz[i]); seen[[k]] <- c(seen[[k]], i); keep[length(keep) + 1L] <- i
    if (length(keep) >= n_max) break
  }
  keep
}

# Location clusters: samples within `radius_m` metres of an earlier one share its cluster (same 3D-coordinate hash as thin_sites).
# Returns an integer cluster id per sample.
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

beta_types  <- c("Within lowland", "Within alpine", "Alpine vs lowland")     # pair types

sm1 <- function(v, by = NULL) sprintf("s(%s, bs = '%s', k = %d%s)", v, CFG$bs_smooth, CFG$k_smooth,
                                      if (is.null(by)) "" else paste0(", by = ", by))

# Plot helpers (shared by the figures)
theme_fig <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "none",
        plot.tag = element_text(face = "bold", size = 16),
        strip.background = element_rect(fill = "grey95"), strip.text = element_text(face = "bold"))

# Significance stars: *** p < 0.001, ** p < 0.01, * p < 0.05, ns otherwise
stars <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

# Extract the legend of a ggplot (robust to ggplot2 versions)
get_legend <- function(p) {
  g <- ggplotGrob(p + theme(legend.position = "bottom"))
  for (i in grep("guide-box", g$layout$name))
    if (!inherits(g$grobs[[i]], "zeroGrob")) return(g$grobs[[i]])
  NULL
}

# Assemble rows of 2 panels + a shared legend + a caption into one PNG
assemble_fig <- function(file, rows, leg, caption, width, height, row_heights = rep(1, length(rows))) {
  nr <- length(rows)
  cap_lines <- strwrap(caption, 160)
  png(file, width = width, height = height, units = "in", res = 300)
  grid.newpage()
  # the plot rows share what is left; the legend and caption rows have an absolute height (the caption's depends
  # on its number of lines), so they can neither overlap each other nor be cut off
  pushViewport(viewport(layout = grid.layout(nr + 2, 2,
                        heights = unit.c(unit(row_heights, "null"), unit(0.5, "in"),
                                         unit(0.17 * length(cap_lines) + 0.2, "in")))))
  for (i in seq_len(nr)) for (j in seq_along(rows[[i]])) {
    vp <- viewport(layout.pos.row = i, layout.pos.col = if (length(rows[[i]]) == 1) 1:2 else j)
    el <- rows[[i]][[j]]
    if (inherits(el, "ggplot")) print(el, vp = vp) else { pushViewport(vp); grid.draw(el); popViewport() }   # a plot or a legend
  }
  if (!is.null(leg)) {
    legs <- if (inherits(leg, "grob")) list(leg) else leg
    for (k in seq_along(legs)) {
      pushViewport(viewport(layout.pos.row = nr + 1, layout.pos.col = if (length(legs) == 1) 1:2 else k))
      grid.draw(legs[[k]]); popViewport()
    }
  }
  pushViewport(viewport(layout.pos.row = nr + 2, layout.pos.col = 1:2))
  grid.text(paste(cap_lines, collapse = "\n"), x = 0.01, just = "left", gp = gpar(fontsize = 9))
  popViewport(); dev.off()
}

## ---- 1. Data (shared by both parts) -----------------------------------------------------------------------
# The sample table is a TSV without quoting (its text fields contain commas): read it with quote = "" so that
# a stray double quote cannot swallow the following lines; fill = FALSE makes a misaligned row an error.
read_sample_table <- function(path)
  read.delim(path, quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE, fill = FALSE)

raw <- read_sample_table(CFG$alpha_file)
required <- c("Longitude", "Latitude", "Alpine", "alpha_total_n_reads",
              "alpha_n_otus", "alpha_shannon", "alpha_pielou", CFG$alpha_id,
              "unmapped_raw_counts", "total_raw_counts",
              if (RUN_BETA) CFG$biom_id)
miss <- setdiff(required, names(raw))
if (length(miss)) stop("Missing columns in ", CFG$alpha_file, ": ", paste(miss, collapse = ", "))

# Elevation from the earth-data table, joined by sample ID (to be replaced by the merged table)
if (!file.exists(CFG$env_file)) stop("env_file not found: ", CFG$env_file)
env <- read.csv(CFG$env_file, stringsAsFactors = FALSE)
if (!all(c(CFG$env_id, CFG$env_elev) %in% names(env)))
  stop("env_file lacks column(s): check CFG$env_id / CFG$env_elev")
rep_cols_env <- grep("^rep_[0-9]+$", names(env), value = TRUE)          # lowland replicates written by define_lowland_sets.py
env <- data.frame(.env_id = env[[CFG$env_id]], .elev_m = env[[CFG$env_elev]], env[rep_cols_env], check.names = FALSE)
env <- env[!duplicated(env$.env_id), ]
raw <- dplyr::left_join(raw, env, by = setNames(".env_id", CFG$alpha_id))
message(sprintf("Elevation joined for %d / %d rows", sum(is.finite(raw$.elev_m)), nrow(raw)))



## =============================== PART A: ALPHA DIVERSITY ===============================
if (RUN_ALPHA) {

d <- raw %>%
  mutate(
    Alpine      = factor(Alpine, levels = c("Lowland", "Alpine")),
    lon = Longitude, lat = Latitude,
    glat        = if (CFG$abs_lat) abs(Latitude) else Latitude,
    elevation   = .elev_m / 1000,                         # km
    reads       = alpha_total_n_reads, log10_reads = log10(alpha_total_n_reads),
    richness    = round(alpha_n_otus),
    shannon     = alpha_shannon, pielou = alpha_pielou,
    unmapped_prop = unmapped_raw_counts / total_raw_counts
  ) %>%
  filter(is.finite(lon), is.finite(lat), abs(lon) <= 180, abs(lat) <= 90,
         is.finite(elevation), is.finite(reads), reads > 0,
         is.finite(richness), is.finite(shannon), shannon > 0, is.finite(pielou), is.finite(unmapped_prop),
         !is.na(Alpine))

# ---- x axes of the alpha gradients ------------------------------------------------------------------------------------------------------
# "raw": elevation (km) and latitude. Elevation and latitude are not independent in the data (high samples are at some latitudes), and
# both are tied to space, so the effect of one can carry the other. A sensitivity analysis (own output folder):
#  * "normalised": elevation minus the elevation expected from a smooth lat/lon surface of the sampled elevations (local relief
#    relative to the region), and latitude minus the latitude expected at that elevation. The baselines are fitted to these samples,
#    so alpine samples help define them: a sensitivity analysis, not a replacement of the raw axes.
d$elevation_raw <- d$elevation; d$glat_raw <- d$glat
if (CFG$gradient_axes == "normalised") {
  fe <- bam(elevation_raw ~ s(lat, lon, bs = "sos", k = CFG$k_norm), data = d, method = "fREML", discrete = TRUE, nthreads = CFG$nthreads)
  d$elevation <- as.numeric(d$elevation_raw - fitted(fe))
  fl <- bam(glat_raw ~ s(elevation_raw, bs = "ts", k = 8), data = d, method = "fREML", discrete = TRUE, nthreads = CFG$nthreads)
  d$glat <- as.numeric(d$glat_raw - fitted(fl))
  message(sprintf("normalised axes: the lat/lon surface explains %.0f%% of the variance of elevation, elevation explains %.0f%% of that of latitude; correlation of the normalised with the raw axis: elevation %.2f, latitude %.2f",
                  100 * summary(fe)$r.sq, 100 * summary(fl)$r.sq, cor(d$elevation, d$elevation_raw), cor(d$glat, d$glat_raw)))
}
axes_txt <- switch(CFG$gradient_axes, raw = "elevation and latitude",
  normalised = "elevation normalised for the lat/lon surface (residual) and latitude normalised for elevation (residual)")

# Samples of one BioProject, or at one location (within CFG$thin_m metres), are not independent: left alone they make the curves and
# p-values overconfident. The column cluster_id (microbeatlas_filter_data.R) links such samples and gives every connected group one id; ONE
# random intercept per cluster is used (the effective number of independent units is then the number of clusters). A group can first be
# thinned to one random sample per location (thin_alpine / thin_lowland), which makes the random effect cheaper.
set.seed(7); n0 <- nrow(d); keep_rows <- rep(TRUE, nrow(d))
for (g in levels(d$Alpine)) {
  if (isTRUE(if (g == "Alpine") CFG$thin_alpine else CFG$thin_lowland)) {
    ix <- which(d$Alpine == g); keep_rows[ix] <- FALSE
    keep_rows[ix[thin_sites(d$lat[ix], d$lon[ix], CFG$thin_m)]] <- TRUE
  }
}
d <- d[keep_rows, ]
use_re <- isTRUE(CFG$cluster_re)
re_rhs <- ""; re_labels <- character(0)
if (use_re) {
  if ("cluster_id" %in% names(d)) cl <- as.character(d$cluster_id) else {
    message("  no column 'cluster_id' in the sample table (run microbeatlas_filter_data.R): the clusters are the locations only")
    cl <- paste0("L", cluster_sites(d$lat, d$lon, CFG$thin_m))
  }
  cl[is.na(cl) | !nzchar(cl)] <- "pooled_small"
  tb <- table(cl); cl[cl %in% names(tb)[tb < CFG$min_cluster_n]] <- "pooled_small"
  tb <- sort(table(cl[cl != "pooled_small"]), decreasing = TRUE)
  if (length(tb) > CFG$max_re_levels - 1) {
    cl[cl %in% names(tb)[seq(CFG$max_re_levels, length(tb))]] <- "pooled_small"
    message(sprintf("  cluster random effect: the %d smallest clusters are pooled as well (max_re_levels = %d)", length(tb) - (CFG$max_re_levels - 1), CFG$max_re_levels))
  }
  d$re_cluster <- factor(cl)
  re_rhs <- "+ s(re_cluster, bs = 're')"; re_labels <- "s(re_cluster)"
  big <- sort(table(cl[cl != "pooled_small"]), decreasing = TRUE)
  message(sprintf("Cluster random effect: %d levels (%d clusters with >= %d samples + one pooled level holding %.1f%% of the samples); largest clusters: %s",
                  nlevels(d$re_cluster), length(big), CFG$min_cluster_n, 100 * mean(cl == "pooled_small"),
                  paste(sprintf("%s (%.1f%%)", names(head(big, 3)), 100 * as.integer(head(big, 3)) / nrow(d)), collapse = ", ")))
} else d$re_cluster <- factor(rep("none", nrow(d)))
excl_for <- function(m) intersect(c("s(lat,lon)", re_labels), vapply(m$smooth, function(s) s$label, ""))   # no spatial term, no cluster effect
message(sprintf("%d samples in the models (of %d): alpine %d, lowland %d (thinning: alpine %s, lowland %s)",
                nrow(d), n0, sum(d$Alpine == "Alpine"), sum(d$Alpine == "Lowland"), isTRUE(CFG$thin_alpine), isTRUE(CFG$thin_lowland)))
how_txt <- function(g) if (isTRUE(if (g == "Alpine") CFG$thin_alpine else CFG$thin_lowland)) sprintf("%s samples thinned to at most one per %g m", tolower(g), CFG$thin_m) else sprintf("all %s samples", tolower(g))
alpha_txt <- sprintf("A, B use %s and %s, %s, and a spherical spline for the regional trend.", how_txt("Alpine"), how_txt("Lowland"),
                     if (use_re) sprintf("one random intercept per cluster (samples of one BioProject or within %g m of each other are linked)", CFG$thin_m) else "no cluster random effect")

# Pielou J must be strictly in (0, 1) for beta regression (squeeze only if needed)
if (any(d$pielou <= 0 | d$pielou >= 1)) {
  n <- nrow(d); d$pielou <- (d$pielou * (n - 1) + 0.5) / n
  message("Pielou J contained 0/1 values: applied the Smithson-Verkuilen squeeze.")
}
# Unmapped proportion as a RESPONSE (beta regression): squeezed copy; unmapped_prop stays as the covariate
d$unmapped <- d$unmapped_prop
if (any(d$unmapped <= 0 | d$unmapped >= 1)) {
  n <- nrow(d); d$unmapped <- (d$unmapped * (n - 1) + 0.5) / n
  message("Unmapped proportion contained 0/1 values: applied the Smithson-Verkuilen squeeze.")
}
d$AlpineO <- factor(d$Alpine, levels = c("Lowland", "Alpine"), ordered = TRUE)   # for difference smooths
print(table(alpine = d$Alpine))
write.csv(select(d, any_of(c(CFG$alpha_id, "Alpine", "lon", "lat", "glat", "elevation", "reads",
                              "richness", "shannon", "pielou", "unmapped_prop", "unmapped"))),
          out_a("analysis_table.csv"), row.names = FALSE)

## ---- A1. Models -----------------------------------------------------------------------
families <- list(richness = nb(), shannon = Gamma(link = "log"), pielou = betar(link = "logit"),
                 unmapped = betar(link = "logit"))
resp_lab <- c(richness = "Observed OTUs", shannon = "Shannon index", pielou = "Pielou J",
              unmapped = "Unmapped reads (proportion)")

base_rhs <- paste("Alpine +", sm1("log10_reads"), "+", sm1("elevation"), "+", sm1("elevation", "AlpineO"),
                  "+", sm1("glat"), "+", sm1("glat", "AlpineO"), re_rhs)
unm_rhs   <- paste("+", sm1("unmapped_prop"), "+", sm1("unmapped_prop", "AlpineO"))   # curves differ by group (script 1)
space_rhs <- sprintf("+ s(lat, lon, bs = 'sos', k = %d)", CFG$k_space)
specs <- c(net       = paste(base_rhs, unm_rhs, space_rhs),   # MAIN
           raw       = paste(base_rhs, unm_rhs),              # no spatial term
           net_nounm = paste(base_rhs, space_rhs))            # sensitivity: without the unmapped terms
spec_lab <- c(net = "Net of space (sos)", raw = "No spatial term",
              net_nounm = "Net of space, not adjusted for unmapped")

# the unmapped proportion as a response: same structure, without the unmapped terms (net / raw only)
specs_unm <- c(net = paste(base_rhs, space_rhs), raw = base_rhs)
spec_formula <- function(resp, sp) if (resp == "unmapped") specs_unm[[sp]] else specs[[sp]]

# a saved model is reused only if it has the same number of observations and the same formula as now
fits_file <- out_a("fits.rds")
same_model <- function(m, resp, sp) {
  strip <- function(f) gsub("\\s+", "", paste(deparse(f), collapse = ""))
  !is.null(m) && length(m$y) == nrow(d) && identical(strip(formula(m)), strip(as.formula(paste(resp, "~", spec_formula(resp, sp)))))
}
old_fits <- if (REUSE_FITS && file.exists(fits_file)) readRDS(fits_file) else list()
fits <- list()
for (resp in names(families)) {
  sp_names <- if (resp == "unmapped") names(specs_unm) else names(specs)
  for (sp in sp_names) {
    nm <- paste(resp, sp, sep = "__")
    if (same_model(old_fits[[nm]], resp, sp)) {
      message("Reusing saved model ", nm); fits[[nm]] <- old_fits[[nm]]
    } else {
      message("Fitting ", nm)
      fits[[nm]] <- bam(as.formula(paste(resp, "~", spec_formula(resp, sp))), data = d, family = families[[resp]],
                        method = "fREML", discrete = TRUE, select = CFG$select, nthreads = CFG$nthreads)
    }
  }
}
saveRDS(fits, fits_file)
# how strong is the clustering of samples? standard deviation of the cluster random intercept (link scale)
if (use_re) {
  re_sd <- bind_rows(lapply(names(fits), function(nm) {
    sd_loc <- tryCatch({ vc <- suppressWarnings(gam.vcomp(fits[[nm]], rescale = FALSE))
                         if (is.list(vc) && !is.matrix(vc)) vc <- vc$vc
                         as.numeric(if (is.matrix(vc)) vc[re_labels, 1] else vc[re_labels]) }, error = function(e) NA_real_)
    data.frame(model = nm, cluster_sd = if (length(sd_loc) == 1) sd_loc else NA_real_)
  }))
  write.csv(re_sd, out_a("cluster_random_effect_sd.csv"), row.names = FALSE)
  message(sprintf("Cluster random effect: SD %s (link scale; per model in alpha/cluster_random_effect_sd.csv)",
                  paste(sprintf("%s %.2f", re_sd$model[grepl("__net$", re_sd$model)], re_sd$cluster_sd[grepl("__net$", re_sd$model)]), collapse = ", ")))
}

## ---- A1b. Concurvity of elevation, latitude and space --------------------------------------------------------------
# Concurvity = how well a smooth term can be predicted from the others (the non-linear analogue of collinearity): near 1, the effect of
# elevation cannot be told from that of latitude or of the spatial spline. It is computed for the three smooths that matter here (elevation,
# latitude, spherical spline) in a model WITHOUT the group-difference smooths (a main smooth and its alpine difference of the same variable
# are collinear by construction), on a random sample (all alpine samples up to a third). Two models, because the Alpine group term matters:
# alpine means high elevation by definition, so with the group in the model elevation is "explained" by it. Rows:
#   "the other two smooths"           a term against the other two smooths together (model without the Alpine group)
#   "the other two smooths + Alpine"  the same with the Alpine group term in the model
#   "<term>"                          pairwise, against one other smooth (model without the Alpine group)
# "observed" = for the fitted terms, "worst" = the worst case over all possible fits, "estimate" = a rough average. Above ~0.8 = warning.
if (isTRUE(CFG$concurvity)) {
  set.seed(3)
  ia <- which(d$Alpine == "Alpine"); il <- which(d$Alpine == "Lowland")
  n_a <- min(length(ia), CFG$concurvity_n %/% 3); n_l <- min(length(il), CFG$concurvity_n - n_a)
  take <- c(ia[sample.int(length(ia), n_a)], il[sample.int(length(il), n_l)])
  cc_tab <- tryCatch({
    sm3 <- sprintf("%s + %s + s(lat, lon, bs = 'sos', k = %d)", sm1("elevation"), sm1("glat"), CFG$k_space)
    f0 <- gam(as.formula(paste("log1p(richness) ~", sm3)), data = d[take, ], method = "REML")             # the three smooths only
    f1 <- gam(as.formula(paste("log1p(richness) ~ Alpine +", sm3)), data = d[take, ], method = "REML")    # ... and the Alpine group
    nm <- c("s(elevation)", "s(glat)", "s(lat,lon)")
    c0 <- mgcv::concurvity(f0, full = TRUE); p0 <- mgcv::concurvity(f0, full = FALSE)
    c1 <- mgcv::concurvity(f1, full = TRUE); p1 <- mgcv::concurvity(f1, full = FALSE)
    row <- function(term, versus, m, i, j = term) data.frame(term = term, versus = versus, observed = m$observed[i, j], worst = m$worst[i, j], estimate = m$estimate[i, j])
    out <- rbind(
      data.frame(term = nm, versus = "the other two smooths", observed = c0["observed", nm], worst = c0["worst", nm], estimate = c0["estimate", nm]),
      data.frame(term = nm, versus = "the other two smooths + Alpine", observed = c1["observed", nm], worst = c1["worst", nm], estimate = c1["estimate", nm]))
    for (i in 1:2) for (j in (i + 1):3) out <- rbind(out, row(nm[i], nm[j], p0, nm[i], nm[j]))
    rownames(out) <- NULL; out
  }, error = function(e) { message("  concurvity could not be computed: ", conditionMessage(e)); NULL })
  if (!is.null(cc_tab)) {
    write.csv(cc_tab, out_a("concurvity.csv"), row.names = FALSE)
    show <- function(v) paste(sprintf("%s %.2f / %.2f", cc_tab$term[cc_tab$versus == v], cc_tab$observed[cc_tab$versus == v], cc_tab$worst[cc_tab$versus == v]), collapse = "; ")
    message(sprintf("Concurvity (%d samples; observed / worst), against the other two smooths: %s", length(take), show("the other two smooths")))
    message("  with the Alpine group term as well: ", show("the other two smooths + Alpine"))
    message("  pairwise: ", paste(sprintf("%s vs %s %.2f / %.2f", cc_tab$term, cc_tab$versus, cc_tab$observed, cc_tab$worst)[!cc_tab$versus %in% c("the other two smooths", "the other two smooths + Alpine")], collapse = "; "))
    if (any(cc_tab$observed[cc_tab$versus == "the other two smooths"] > 0.8))
      message("  -> above 0.8 against the other two smooths (see alpha/concurvity.csv): those terms can hardly be separated; their separate effects are weakly determined")
  }
}

## ---- A2. Diagnostics -------------------------------------------------------------------
# Moran's I of deviance residuals (inverse-distance weights, permutation p-value)
moran_perm <- function(m, d, n_max = CFG$moran_n_max, nperm = CFG$moran_nperm) {
  idx <- sample(nrow(d), min(n_max, nrow(d)))
  r <- residuals(m, type = "deviance")[idx]; z <- r - mean(r)
  D <- geosphere::distm(cbind(d$lon[idx], d$lat[idx]), fun = geosphere::distHaversine) / 1000
  D[D < 1] <- 1; W <- 1 / D; diag(W) <- 0
  I <- function(z) (length(z) / sum(W)) * drop(crossprod(z, W %*% z)) / sum(z^2)
  obs <- I(z); perm <- replicate(nperm, I(sample(z)))
  c(moran_I = obs, moran_p = (1 + sum(perm >= obs)) / (nperm + 1))
}

comp <- list(); smooths <- list(); conc <- list()
for (nm in names(fits)) {
  m <- fits[[nm]]; parts <- strsplit(nm, "__")[[1]]; message("Checking ", nm)
  sm <- summary(m)
  smooths[[nm]] <- data.frame(model = nm, response = parts[1], spec = parts[2],
                              term = rownames(sm$s.table), sm$s.table, row.names = NULL)
  writeLines(c(paste(parts[1], "~", spec_formula(parts[1], parts[2])), "", capture.output(print(sm)),
               "", capture.output(print(k.check(m)))), out_a(paste0("report_", nm, ".txt")))
  tryCatch({
    png(out_a(paste0("resid_", nm, ".png")), 1200, 500, res = 120)
    tryCatch({
      par(mfrow = c(1, 2))
      qq.gam(m, rep = 0, main = "QQ (deviance residuals)")
      plot(fitted(m), residuals(m, type = "deviance"), pch = 16, cex = 0.3,
           col = rgb(0, 0, 0, 0.2), xlab = "Fitted", ylab = "Deviance residual")
    }, finally = dev.off())
  }, error = function(e) message("  residual plot failed for ", nm))
  if (parts[2] == "net")
    conc[[nm]] <- tryCatch(round(concurvity(m, full = TRUE), 2), error = function(e) NULL)
  mo <- tryCatch(moran_perm(m, d), error = function(e) c(moran_I = NA, moran_p = NA))
  comp[[nm]] <- data.frame(response = parts[1], spec = parts[2], n = nrow(d), AIC = AIC(m),
                           dev_expl = sm$dev.expl, moran_I = mo[["moran_I"]], moran_p = mo[["moran_p"]])
}
comp <- bind_rows(comp)
write.csv(comp, out_a("model_comparison.csv"), row.names = FALSE)
write.csv(bind_rows(smooths), out_a("smooth_terms.csv"), row.names = FALSE)   # edf, F/Chi.sq, p per smooth
capture.output(print(conc), file = out_a("concurvity_net_models.txt"))
print(comp)

## ---- A3. Predicted gradients -----------------------------------------------------------------
# Vary one gradient within each group's 1-99% range; the other gradient is held at the
# group median, depth at the overall median, unmapped proportion at the group median,
# spatial term excluded.
predict_gradient <- function(m, var) {
  excl <- intersect(c("s(lat,lon)", re_labels), vapply(m$smooth, function(s) s$label, ""))   # no spatial term, no location random effect
  bind_rows(lapply(levels(d$Alpine), function(a) {
    sub <- d[d$Alpine == a, ]
    nd <- data.frame(Alpine = factor(a, levels = levels(d$Alpine)),
                     re_cluster = factor(levels(d$re_cluster)[1], levels = levels(d$re_cluster)),
                     AlpineO = factor(a, levels = levels(d$Alpine), ordered = TRUE),
                     elevation = median(sub$elevation), glat = median(sub$glat),
                     log10_reads = median(d$log10_reads), unmapped_prop = median(sub$unmapped_prop),
                     lat = median(sub$lat), lon = median(sub$lon))
    nd <- nd[rep(1, 100), ]
    nd[[var]] <- seq(quantile(sub[[var]], .01), quantile(sub[[var]], .99), length.out = 100)
    p <- predict(m, nd, type = "link", se.fit = TRUE, exclude = excl, newdata.guaranteed = TRUE, discrete = FALSE)
    li <- family(m)$linkinv
    data.frame(Alpine = nd$Alpine, x = nd[[var]], fit = li(p$fit),
               lo = li(p$fit - 1.96 * p$se.fit), hi = li(p$fit + 1.96 * p$se.fit))
  }))
}

grad_lab <- c(elevation = switch(CFG$gradient_axes, raw = "Elevation (km)", normalised = "Normalised elevation (km)"),
              glat = if (CFG$gradient_axes == "normalised") "Normalised latitude (\u00B0)" else if (CFG$abs_lat) "Absolute latitude (\u00B0)" else "Latitude (\u00B0)")
pred <- bind_rows(lapply(names(fits), function(nm) {
  parts <- strsplit(nm, "__")[[1]]
  if (!parts[2] %in% c("net", "raw")) return(NULL)        # net_nounm: tables only
  bind_rows(lapply(names(grad_lab), function(v)
    cbind(response = parts[1], spec = parts[2], gradient = grad_lab[[v]],
          predict_gradient(fits[[nm]], v))))
}))
pred$response <- factor(resp_lab[pred$response], levels = resp_lab)
pred$gradient <- factor(pred$gradient, levels = grad_lab)
pred$spec     <- factor(spec_lab[pred$spec], levels = spec_lab[c("net", "raw")])
write.csv(pred, out_a("predicted_gradients.csv"), row.names = FALSE)

# Alpine-vs-lowland difference smooths: does the gradient shape differ between groups?
# difference smooths of the unmapped proportion are kept in a separate table (not part of the figure)
write.csv(bind_rows(smooths) %>% filter(grepl("AlpineO", term), grepl("unmapped_prop", term)),
          out_a("group_difference_tests_unmapped.csv"), row.names = FALSE)
diff_tests <- bind_rows(smooths) %>%
  filter(grepl("AlpineO", term), grepl("elevation|glat", term)) %>%
  mutate(gradient = factor(ifelse(grepl("elevation", term), grad_lab[["elevation"]], grad_lab[["glat"]]),
                           levels = grad_lab),
         response = factor(resp_lab[response], levels = resp_lab))
write.csv(diff_tests, out_a("group_difference_tests.csv"), row.names = FALSE)
ann <- diff_tests %>% filter(spec == "net") %>%
  mutate(label = sprintf("Alpine vs lowland shape: %s (p %s)", stars(p.value),
                         ifelse(p.value < 0.001, "< 0.001", sprintf("= %.3f", p.value))))

# rug marks (a subsample, for drawing only)
rug <- bind_rows(lapply(names(grad_lab), function(v)
  d %>% group_by(Alpine) %>% slice_sample(n = CFG$n_rug) %>% ungroup() %>%
    transmute(Alpine, gradient = factor(grad_lab[[v]], levels = grad_lab), x = .data[[v]])))

## ---- A4. Alpha-diversity figure -------------------------------------------------------------------------------
n_grp <- table(d$Alpine)
grp_labels <- setNames(sprintf("%s (n = %s)", names(n_grp), format(as.integer(n_grp), big.mark = ",")),
                       names(n_grp))

g <- ggplot(pred, aes(x, fit, colour = Alpine)) +
  geom_ribbon(data = filter(pred, spec == spec_lab[["net"]]),
              aes(ymin = lo, ymax = hi, fill = Alpine), alpha = 0.25, colour = NA) +
  geom_line(aes(linetype = spec), linewidth = 1) +
  geom_rug(data = rug, aes(x = x, colour = Alpine), inherit.aes = FALSE,
           sides = "b", alpha = 0.05, length = unit(0.03, "npc")) +
  geom_text(data = ann, aes(label = label), x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5,
            size = 3.2, inherit.aes = FALSE) +
  facet_grid(response ~ gradient, scales = "free", switch = "both") +
  scale_colour_manual(values = CFG$cols, labels = grp_labels, name = NULL) +
  scale_fill_manual(values = CFG$cols, guide = "none") +
  scale_linetype_manual(values = c("solid", "dashed"), name = NULL) +
  labs(x = NULL, y = NULL,
       subtitle = paste("GAM fits: solid = net of space (95% CI), dashed = no spatial term;",
                        "p = test that the alpine gradient differs in shape from the lowland one (net model)")) +
  theme_bw(base_size = 12) +
  theme(strip.placement = "outside", strip.background = element_blank(),
        strip.text = element_text(face = "bold"), panel.grid.minor = element_blank(),
        legend.position = "bottom")
ggsave(out_a("alpha_full_grid.png"), g, width = 10, height = 12.5, dpi = 300)
saveRDS(g, out_a("alpha_full_grid_plot.rds"))

# Panel builder for the main figure: one gradient, the chosen alpha metrics as facets
main_alpha_labels <- resp_lab[CFG$main_alpha_metrics]
ann_short <- ann %>% mutate(label = sub("Alpine vs lowland shape: ", "", label))
alpha_panel <- function(gr_label, tag, resp_keep = main_alpha_labels) {
  pr <- pred[pred$gradient == gr_label & pred$response %in% resp_keep, ]
  pr$response <- droplevels(pr$response)
  an <- ann_short[ann_short$gradient == gr_label & ann_short$response %in% resp_keep, ]
  rg <- rug[rug$gradient == gr_label, ]
  ggplot(pr, aes(x, fit, colour = Alpine)) +
    geom_ribbon(data = pr[pr$spec == spec_lab[["net"]], ],
                aes(ymin = lo, ymax = hi, fill = Alpine), alpha = 0.25, colour = NA) +
    geom_line(aes(linetype = spec), linewidth = 1) +
    geom_rug(data = rg, aes(x = x, colour = Alpine), inherit.aes = FALSE,
             sides = "b", alpha = 0.05, length = unit(0.03, "npc")) +
    geom_text(data = an, aes(label = label), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
              size = 3, inherit.aes = FALSE) +
    facet_wrap(~response, nrow = 1, scales = "free_y", labeller = label_wrap_gen(12)) +
    scale_colour_manual(values = CFG$cols, name = NULL) +
    scale_fill_manual(values = CFG$cols, guide = "none") +
    scale_linetype_manual(values = c("solid", "dashed"), name = NULL) +
    labs(tag = tag, x = gr_label, y = NULL) + theme_fig
}

# Supplementary Figure 2: unmapped read proportion along elevation (A) and latitude (B)
pS2A <- alpha_panel(grad_lab[["elevation"]], "A", resp_lab[["unmapped"]])
pS2B <- alpha_panel(grad_lab[["glat"]], "B", resp_lab[["unmapped"]])
assemble_fig(file.path(out_folder, "figS2_unmapped_gradients.png"), list(list(pS2A, pS2B)),
             get_legend(pS2A),
             paste("Supplementary Figure 2. Proportion of reads that could not be mapped to the reference OTUs along",
                   "elevation (A) and latitude (B), alpine vs lowland. Beta-regression GAM partial effects adjusted for",
                   "sequencing depth; solid = net of space (95% CI), dashed = no spatial term; stars = test that the alpine",
                   "gradient differs in shape from the lowland one (net model; *** p < 0.001, ** p < 0.01, * p < 0.05, ns)."),
             10, 4.8)

## ---- A5. Supplementary Figure 3: spatial surfaces of the unmapped proportion ----------------------------
# Where is the reference database less complete? One spatial spline (sos) for the lowland samples and
# an alpine-specific one (difference surface), for the unmapped proportion at the median sequencing depth.
# No elevation / latitude terms here: the surface is the regional pattern itself.
m_sp_formula <- as.formula(sprintf(paste("unmapped ~ Alpine + %s + s(lat, lon, bs = 'sos', k = %d)",
                                           "+ s(lat, lon, bs = 'sos', by = AlpineO, k = %d) %s"),
                                     sm1("log10_reads"), CFG$k_space, CFG$k_space_alpine, re_rhs))
strip_f <- function(f) gsub("\\s+", "", paste(deparse(f), collapse = ""))
m_sp_file <- out_a("figS3_model.rds")
m_sp <- if (REUSE_FITS && file.exists(m_sp_file)) readRDS(m_sp_file) else NULL
if (is.null(m_sp) || length(m_sp$y) != nrow(d) || !identical(strip_f(formula(m_sp)), strip_f(m_sp_formula))) {
  m_sp <- bam(m_sp_formula, data = d, family = betar(link = "logit"), method = "fREML", discrete = TRUE,
              select = CFG$select, nthreads = CFG$nthreads)
  saveRDS(m_sp, m_sp_file)
} else message("Reusing saved model for Supplementary Figure 3")
writeLines(capture.output(print(summary(m_sp))), out_a("figS3_spatial_model_summary.txt"))

res <- CFG$grid_res; nx <- 360 / res; ny <- 180 / res
grid <- expand.grid(lon = seq(-180 + res / 2, 180 - res / 2, by = res),      # lon varies fastest
                    lat = seq(-90 + res / 2, 90 - res / 2, by = res))

# number of samples within the 3x3 cell neighbourhood of every cell (support of the data)
neighbourhood_n <- function(sub) {
  cx <- pmin(floor((sub$lon + 180) / res) + 1, nx); cy <- pmin(floor((sub$lat + 90) / res) + 1, ny)
  m <- as.matrix(table(factor(cy, levels = seq_len(ny)), factor(cx, levels = seq_len(nx))))
  out <- matrix(0, ny, nx)
  for (di in -1:1) for (dj in -1:1) {
    ri <- seq_len(ny) + di; ci <- seq_len(nx) + dj
    okr <- ri >= 1 & ri <= ny; okc <- ci >= 1 & ci <= nx
    out[okr, okc] <- out[okr, okc] + m[ri[okr], ci[okc]]
  }
  as.vector(t(out))                                                          # same order as `grid`
}

grp_fit <- bind_rows(lapply(levels(d$Alpine), function(g) {
  nd <- grid
  nd$Alpine  <- factor(g, levels = levels(d$Alpine))
  nd$AlpineO <- factor(g, levels = levels(d$Alpine), ordered = TRUE)
  nd$log10_reads <- median(d$log10_reads)
  nd$re_cluster <- factor(levels(d$re_cluster)[1], levels = levels(d$re_cluster))             # no cluster random effect on the map
  nd$fit <- as.numeric(predict(m_sp, nd, type = "response", newdata.guaranteed = TRUE, discrete = FALSE,
                               exclude = intersect(re_labels, vapply(m_sp$smooth, function(s) s$label, ""))))
  nd$n_near <- neighbourhood_n(d[d$Alpine == g, ])
  nd$fit_shown <- ifelse(nd$n_near >= CFG$min_support, nd$fit, NA_real_)
  nd$group <- g
  nd
}))
write.csv(grp_fit[, c("group", "lon", "lat", "fit", "n_near")], out_a("figS3_spatial_surfaces.csv"), row.names = FALSE)

lims <- range(grp_fit$fit_shown, na.rm = TRUE)
world <- ggplot2::map_data("world")
surface_panel <- function(g, tag, title) {
  gg <- grp_fit[grp_fit$group == g, ]
  p <- ggplot() +
    geom_raster(data = gg, aes(lon, lat, fill = fit_shown)) +
    geom_polygon(data = world, aes(long, lat, group = group), fill = NA, colour = "grey40", linewidth = 0.15) +
    scale_fill_viridis_c(option = "C", limits = lims, na.value = "grey93",
                         name = "Predicted unmapped\nproportion") +
    coord_cartesian(xlim = c(-180, 180), ylim = CFG$map_ylim, expand = FALSE) +
    labs(tag = tag, title = title, x = NULL, y = NULL) +
    theme_fig + theme(panel.background = element_rect(fill = "grey93"))
  if (g == "Alpine")                                                         # show where the alpine samples are
    p <- p + geom_point(data = d[d$Alpine == "Alpine", ], aes(lon, lat), size = 0.25, alpha = 0.5, colour = "white")
  p
}
pS3A <- surface_panel("Lowland", "A", sprintf("Lowland samples (n = %s)", format(sum(d$Alpine == "Lowland"), big.mark = ",")))
pS3B <- surface_panel("Alpine",  "B", sprintf("Alpine samples (n = %s; white dots)", format(sum(d$Alpine == "Alpine"), big.mark = ",")))
assemble_fig(file.path(out_folder, "figS3_spatial_unmapped.png"), list(list(pS3A), list(pS3B)),
             get_legend(pS3A + theme(legend.key.width = unit(1.4, "cm"))),
             paste("Supplementary Figure 3. Spatial surface (spherical spline) of the proportion of reads that could not be",
                   "mapped to the reference OTUs, at the median sequencing depth, for lowland (A) and alpine (B) samples.",
                   sprintf("Higher values (lighter colours) mark regions that are less complete in the reference. Cells with fewer than %d samples", CFG$min_support),
                   sprintf("within their %d x %d degree neighbourhood are not shown (grey), as the surface there is extrapolated.", 3 * res, 3 * res)),
             12, 9)

}  # end RUN_ALPHA

## =============================== PART B: BETA DIVERSITY ================================
if (RUN_BETA) {

meta <- raw %>%
  transmute(
    id        = .data[[CFG$biom_id]],
    Alpine    = factor(Alpine, levels = c("Lowland", "Alpine")),
    lon       = Longitude, lat = Latitude,
    elevation = .elev_m / 1000,                                     # km
    log10_reads = log10(alpha_total_n_reads),
    unmapped  = unmapped_raw_counts / total_raw_counts,
    across(all_of(rep_cols_env)),                                   # replicate flags (yes / no) of the lowland samples
    across(any_of("cluster_id"))                                    # cluster of the sample (BioProject / location)
  ) %>%
  filter(is.finite(lon), is.finite(lat), abs(lon) <= 180, abs(lat) <= 90,
         is.finite(elevation), is.finite(log10_reads), is.finite(unmapped), !is.na(Alpine), !is.na(id))
print(table(group = meta$Alpine))

## ---- B1. Replicate subsamples ------------------------------------------------------------------
# Every replicate compares ONE alpine sample per cluster with the SAME NUMBER of non-alpine samples. cluster_id (microbeatlas_filter_data.R)
# links the samples of one BioProject or one location (within 100 m); a random sample of each cluster is taken, a different one in every
# replicate, so the alpine samples of a replicate come from different studies / places and no cluster is over-represented. The non-alpine
# samples of replicate r are drawn from rep_r (define_lowland_sets.py: random draws evenly spread over latitude / longitude) and
# thinned to at most one per thin_m metres; if there are more than alpine samples they are cut down EVENLY (one per lat/lon cell in turn).
rep_cols <- grep("^rep_[0-9]+$", names(meta), value = TRUE)
use_sets <- isTRUE(CFG$use_replicate_sets) && length(rep_cols) >= CFG$n_rep
message(if (use_sets) sprintf("Replicate pools: the non-alpine samples flagged in %s ... %s (random draws evenly spaced over latitude / longitude)", rep_cols[1], rep_cols[CFG$n_rep])
        else "No replicate pools (use_replicate_sets = FALSE, or fewer than n_rep columns rep_01, ... in the earth-data table): the non-alpine samples are a plain random draw from all of them, cut down evenly over latitude / longitude")
is_alp <- meta$Alpine == "Alpine"
if ("cluster_id" %in% names(meta)) cl <- as.character(meta$cluster_id) else {
  message("  no column 'cluster_id' in the sample table (run microbeatlas_filter_data.R): the clusters are the locations only")
  cl <- rep(NA_character_, nrow(meta)); cl[is_alp] <- paste0("L", cluster_sites(meta$lat[is_alp], meta$lon[is_alp], CFG$thin_m))
}
no_cl <- is.na(cl) | !nzchar(cl); cl[no_cl] <- paste0("none", seq_len(sum(no_cl)))      # a sample without a cluster is its own cluster
meta$cl <- cl
alp <- meta[is_alp, ]
message(sprintf("Alpine samples: %s in %s clusters (BioProject / location); every replicate takes one random sample per cluster",
                format(nrow(alp), big.mark = ","), format(length(unique(alp$cl)), big.mark = ",")))

# order in which to take samples so that any first n are spread as evenly as possible over cell_deg-degree latitude / longitude cells:
# one sample per cell in turn (random sample within a cell, random cell order in a partial round)
even_order <- function(lat, lon, cell_deg = CFG$even_cell_deg) {
  cell <- paste(floor(lat / cell_deg), floor(lon / cell_deg))
  perm <- sample.int(length(lat)); rk <- integer(length(lat))
  rk[perm] <- ave(seq_along(perm), cell[perm], FUN = seq_along)                # rank of a sample within its cell
  order(rk, match(cell, sample(unique(cell))))
}
sel <- lapply(seq_len(CFG$n_rep), function(r) {
  set.seed(100 + r)
  a <- alp[sample.int(nrow(alp)), ]; a <- a[!duplicated(a$cl), ]                 # one random sample per cluster
  if (nrow(a) > CFG$n_alpine_max) a <- a[sample.int(nrow(a), CFG$n_alpine_max), ]
  pool <- if (use_sets) meta[meta$Alpine == "Lowland" & meta[[rep_cols[r]]] %in% "yes", ] else meta[meta$Alpine == "Lowland", ]
  if (!use_sets) pool <- pool[sample.int(nrow(pool), min(nrow(pool), 20L * nrow(a))), ]     # a manageable random draw before thinning
  pool <- pool[thin_sites(pool$lat, pool$lon, CFG$thin_m), ]
  n_use <- min(nrow(a), nrow(pool))                                               # the same number in both groups
  if (n_use < nrow(a)) a <- a[sample.int(nrow(a), n_use), ]
  bind_rows(a, pool[even_order(pool$lat, pool$lon)[seq_len(n_use)], ])
})
cells_of <- function(x) length(unique(paste(floor(x$lat / CFG$even_cell_deg), floor(x$lon / CFG$even_cell_deg))))
samples_per_rep <- bind_rows(lapply(seq_along(sel), function(r) {
  x <- sel[[r]]; data.frame(replicate = r, alpine = sum(x$Alpine == "Alpine"), non_alpine = sum(x$Alpine == "Lowland"),
                            alpine_clusters = length(unique(x$cl[x$Alpine == "Alpine"])),
                            cells_alpine = cells_of(x[x$Alpine == "Alpine", ]), cells_non_alpine = cells_of(x[x$Alpine == "Lowland", ]))
}))
write.csv(samples_per_rep, out_b("samples_per_replicate.csv"), row.names = FALSE)
message(sprintf("Every replicate: %s alpine samples (one per cluster) and %s non-alpine samples; the non-alpine ones cover %s latitude / longitude cells (%g degrees), the alpine ones %s",
                if (min(samples_per_rep$alpine) == max(samples_per_rep$alpine)) format(samples_per_rep$alpine[1], big.mark = ",") else sprintf("%d-%d", min(samples_per_rep$alpine), max(samples_per_rep$alpine)),
                if (min(samples_per_rep$non_alpine) == max(samples_per_rep$non_alpine)) format(samples_per_rep$non_alpine[1], big.mark = ",") else sprintf("%d-%d", min(samples_per_rep$non_alpine), max(samples_per_rep$non_alpine)),
                round(median(samples_per_rep$cells_non_alpine)), CFG$even_cell_deg, round(median(samples_per_rep$cells_alpine))))
stopifnot(all(samples_per_rep$alpine == samples_per_rep$non_alpine))              # the groups really are the same size
stopifnot(all(vapply(sel, function(x) !any(duplicated(x$cl[x$Alpine == "Alpine"])), TRUE)))   # one alpine sample per cluster
all_ids <- unique(unlist(lapply(sel, `[[`, "id")))
message(length(all_ids), " distinct samples across ", CFG$n_rep, " replicate(s)")

to_sparse <- function(m) {
  if (inherits(m, "simple_triplet_matrix"))
    Matrix::sparseMatrix(i = m$i, j = m$j, x = m$v, dims = c(m$nrow, m$ncol),
                         dimnames = list(m$dimnames[[1]], m$dimnames[[2]]))
  else as(m, "dgCMatrix")
}

cache_file <- out_b("subsample_counts.rds")
if (file.exists(cache_file) && all(all_ids %in% colnames(X_cached <- readRDS(cache_file)))) {
  X_all <- X_cached; message("Reusing cached counts: ", cache_file)
} else {
  message("Reading BIOM (slow step): ", CFG$biom_file)
  biom <- rbiom::read_biom(CFG$biom_file)
  ids_ok <- intersect(all_ids, biom$samples)
  if (length(ids_ok) < length(all_ids))
    warning(length(all_ids) - length(ids_ok), " selected samples are not in the BIOM")
  biom <- biom[ids_ok]
  X_all <- to_sparse(biom$counts)                                  # OTUs x samples
  X_all <- X_all[rownames(X_all) != "Unmapped", , drop = FALSE]
  X_all <- X_all[Matrix::rowSums(X_all) > 0, , drop = FALSE]
  rm(biom); invisible(gc())
  saveRDS(X_all, cache_file)
}
message(sprintf("Count matrix: %d OTUs x %d samples", nrow(X_all), ncol(X_all)))

# Rarefy every sample (column) of a sparse count matrix to `depth` reads, without
# replacement. All samples must have >= depth reads. Works on non-zero entries only.
rarefy_sparse <- function(x, depth) {
  newv <- x@x
  for (j in seq_len(ncol(x))) {
    idx  <- seq.int(x@p[j] + 1L, x@p[j + 1L])
    cts  <- x@x[idx]
    pool <- rep.int(seq_along(cts), cts)
    newv[idx] <- tabulate(pool[sample.int(length(pool), depth)], nbins = length(cts))
  }
  x@x <- newv
  Matrix::drop0(x)
}

## ---- B2. Random pairs and their dissimilarities for one replicate -----------------------------------------------
make_pairs <- function(r) {
  set.seed(200 + r)
  S <- sel[[r]] %>% filter(id %in% colnames(X_all))
  X <- X_all[, S$id, drop = FALSE]
  tot <- Matrix::colSums(X)
  ok <- tot >= CFG$rare_depth
  S <- S[ok, ]; X <- X[, ok, drop = FALSE]; tot <- tot[ok]

  # up to n_pairs pairs per pair type (indices into S), drawn at random WITHIN strata = geographic distance class x elevation-difference
  # class (CFG$dist_breaks_km, CFG$elev_breaks_km), the n_pairs shared as evenly as possible over the strata that exist. Random pairs are
  # mostly thousands of km apart and, among lowland samples, hardly differ in elevation: without stratification short distances and
  # large elevation differences would hardly be represented. Stratifying on predictors of the models does not bias the fitted curves; it
  # makes them well determined over the whole range. Candidates are random index pairs (not all pairs: there are millions). The natural
  # (unstratified) median distance and elevation difference of the pair type are kept (nat_geo, nat_delev): they are the "typical" values
  # at which other variables are held in the figures.
  allocate_even <- function(capacity, n) {                 # n over the strata as evenly as possible, never above a stratum's size
    alloc <- integer(length(capacity)); active <- which(capacity > 0); remaining <- min(n, sum(capacity))
    while (remaining > 0 && length(active)) {
      share <- max(remaining %/% length(active), 1L)
      for (k in active[sample.int(length(active))]) {
        give <- min(share, capacity[k] - alloc[k], remaining)
        alloc[k] <- alloc[k] + give; remaining <- remaining - give
        if (remaining == 0) break
      }
      active <- active[alloc[active] < capacity[active]]
    }
    alloc
  }
  draw_pairs <- function(i1, i2, within) {
    M <- 300L * CFG$n_pairs                                                    # random candidate pairs
    a <- i1[sample.int(length(i1), M, TRUE)]; b <- i2[sample.int(length(i2), M, TRUE)]
    key <- if (within) pmin(a, b) * 1e6 + pmax(a, b) else a * 1e6 + b
    ok <- a != b & !duplicated(key); a <- a[ok]; b <- b[ok]
    d <- geosphere::distHaversine(cbind(S$lon[a], S$lat[a]), cbind(S$lon[b], S$lat[b])) / 1000          # km
    e <- abs(S$elevation[a] - S$elevation[b])                                                           # km
    strata <- 100L * cut(d, CFG$dist_breaks_km, labels = FALSE, include.lowest = TRUE) +
                     cut(e, CFG$elev_breaks_km, labels = FALSE, include.lowest = TRUE)
    members <- split(seq_along(strata), strata)
    alloc <- allocate_even(lengths(members), CFG$n_pairs)
    pick <- unlist(Map(function(ix, k) if (k >= length(ix)) ix else ix[sample.int(length(ix), k)], members, alloc), use.names = FALSE)
    data.frame(i = a[pick], j = b[pick], nat_delev = median(e), nat_geo = median(d))
  }
  # n_pairs for EACH pair type: within non-alpine (lowland), within alpine, between (the first sample of a between pair is the alpine one)
  iL <- which(S$Alpine == "Lowland"); iA <- which(S$Alpine == "Alpine")
  one <- function(d, type) if (is.null(d) || !nrow(d)) NULL else data.frame(type = type, i = d$i, j = d$j, nat_delev = d$nat_delev, nat_geo = d$nat_geo)
  pairs <- bind_rows(one(if (length(iL) > 1) draw_pairs(iL, iL, TRUE), beta_types[1]),
                     one(if (length(iA) > 1) draw_pairs(iA, iA, TRUE), beta_types[2]),
                     one(if (length(iA) > 0 && length(iL) > 0) draw_pairs(iA, iL, FALSE), beta_types[3]))
  if (nrow(pairs) == 0) return(NULL)

  # Distances for the drawn pairs only (sparse, vectorised over pairs):
  #   Bray-Curtis on relative abundances p, q:  sum|p - q| / 2  (= 1 - sum min(p, q))
  #   Jaccard on rarefied presence/absence:     1 - |A & B| / |A | B|
  Prel <- X %*% Matrix::Diagonal(x = 1 / tot)
  Br <- rarefy_sparse(X, CFG$rare_depth); Br@x[] <- 1
  Pa <- Prel[, pairs$i, drop = FALSE]; Pb <- Prel[, pairs$j, drop = FALSE]
  Ba <- Br[, pairs$i, drop = FALSE];   Bb <- Br[, pairs$j, drop = FALSE]
  pairs$bray <- as.numeric(Matrix::colSums(abs(Pa - Pb))) / 2
  inter <- as.numeric(Matrix::colSums(Ba * Bb))
  pairs$jaccard <- as.numeric(1 - inter / (Matrix::colSums(Ba) + Matrix::colSums(Bb) - inter))

  if (TEST_MODE) {                       # check against vegan::vegdist on a few pairs
    err <- 0
    for (k in head(seq_len(nrow(pairs)), 25)) {
      a <- as.numeric(Prel[, pairs$i[k]]); b <- as.numeric(Prel[, pairs$j[k]])
      err <- max(err, abs(as.numeric(vegan::vegdist(rbind(a, b), method = "bray")) - pairs$bray[k]))
      ja <- as.numeric(Br[, pairs$i[k]]); jb <- as.numeric(Br[, pairs$j[k]])
      err <- max(err, abs(as.numeric(vegan::vegdist(rbind(ja, jb), method = "jaccard", binary = TRUE)) - pairs$jaccard[k]))
    }
    message(sprintf("  TEST_MODE: max |difference| vs vegan::vegdist = %.2e", err))
    stopifnot(err < 1e-8)
  }

  pairs$rep_id <- r
  pairs$d_elev <- abs(S$elevation[pairs$i] - S$elevation[pairs$j])
  pairs$d_lat  <- abs(S$lat[pairs$i] - S$lat[pairs$j])
  pairs$geo_km <- geosphere::distHaversine(cbind(S$lon[pairs$i], S$lat[pairs$i]),
                                           cbind(S$lon[pairs$j], S$lat[pairs$j])) / 1000
  pairs$d_depth <- abs(S$log10_reads[pairs$i] - S$log10_reads[pairs$j])   # difference in sequencing depth
  pairs$d_unm   <- abs(S$unmapped[pairs$i] - S$unmapped[pairs$j])          # difference in unmapped proportion
  pairs$ca <- S$cl[pairs$i]; pairs$cb <- S$cl[pairs$j]                       # the clusters of the two samples (the independent units of the standard errors)
  pairs[, c("rep_id", "type", "bray", "jaccard", "d_elev", "d_lat", "geo_km", "d_depth", "d_unm", "nat_delev", "nat_geo", "ca", "cb")]
}

## ---- B3. Models ---------------------------------------------------------------------------------------
# ONE ordinary LINEAR REGRESSION per metric and model on the pairs of all replicates, with standard errors clustered by cluster (see the header).
# In every replicate the alpine samples (one per cluster) are compared with the same number of evenly spread non-alpine samples. Three pair
# types (within lowland = non-alpine, within alpine, alpine vs lowland) each have an intercept and their own slopes. The dissimilarity itself
# (0-1) is modelled, so effects are in units of dissimilarity (a beta regression does not converge here: many dissimilarities are exactly 1).
#
# Cluster-robust variance for dyadic data (Aronow, Samii & Assenova 2015): pairs that share a node (here: a cluster) are dependent, so the
# scores of all pairs of a node are summed per node
dyadic_vcov <- function(X, e, node_a, node_b) {
  s <- X * e                                                      # score of every pair (n x k)
  bread <- tryCatch(solve(crossprod(X)), error = function(err) NULL)
  if (is.null(bread)) return(NULL)
  S <- rowsum(rbind(s, s), c(node_a, node_b))                     # sum of the scores of all pairs of each node
  meat <- crossprod(S) - crossprod(s)                             # (a pair is in the sums of both of its nodes: count it once)
  G <- nrow(S)
  V <- bread %*% meat %*% bread * (G / (G - 1))
  attr(V, "G") <- G                                               # number of nodes = independent units
  V
}
# Design-matrix rows of a linear model at new data
design_at <- function(m, newdata) {
  tt <- delete.response(terms(m))
  model.matrix(tt, model.frame(tt, newdata, xlev = m$xlevels), contrasts.arg = m$contrasts)
}
beta_metrics <- c(bray = "Bray\u2013Curtis", jaccard = "Jaccard (rarefied)")
beta_vars    <- c(d_elev = "Elevation difference (km)", d_lat = "Latitude difference (\u00B0)")    # shown in the figures; geographic distance is in the
                                                                                                    # elevation model as a covariate, but is not drawn
curve_model  <- c(d_elev = "elev", d_lat = "lat")      # which model each curve comes from

prep_pairs <- function(P) {
  P <- as.data.frame(P)
  tab <- P %>% group_by(type) %>%
    summarise(n = n(), u = min(n_distinct(round(d_lat, 3)), n_distinct(round(d_elev, 3)), n_distinct(round(geo_km, 3))), .groups = "drop")
  keep <- tab$type[tab$n >= CFG$min_pairs & tab$u >= 20]
  dropped <- setdiff(unique(P$type), keep)
  if (length(dropped)) message("    dropped (fewer than ", CFG$min_pairs, " pairs or fewer than 20 distinct values): ", paste(dropped, collapse = "; "))
  P <- P[P$type %in% keep, ]
  P$type <- factor(P$type, levels = intersect(beta_types, keep))
  P$log10_geo <- log10(pmax(P$geo_km, 0.1))
  P
}

# (elev) elevation difference accounting for geographic distance; (lat) latitude difference accounting for the elevation difference.
# One intercept and one slope per pair TYPE; depth and unmapped-proportion differences as linear nuisance terms.
beta_formula <- function(metric, kind)
  as.formula(paste(metric, "~ type + d_depth + d_unm +", if (kind == "elev") "type:d_elev + type:log10_geo" else "type:d_lat + type:d_elev"))
clip01 <- function(x) pmin(pmax(x, 0), 1)

# typical values per pair type (the natural median distance and elevation difference of the type, before stratification), and the
# medians of the nuisance differences: the values at which the other variables are held when one curve is drawn
typ <- function(P) {
  tp <- P %>% group_by(type) %>% summarise(d_elev = median(nat_delev), d_lat = median(d_lat), log10_geo = log10(pmax(median(nat_geo), 0.1)), .groups = "drop")
  tp$d_depth <- median(P$d_depth); tp$d_unm <- median(P$d_unm); tp
}
curve_nd <- function(tp, P, v) bind_rows(lapply(levels(P$type), function(ty) {
  sub <- P[P$type == ty, ]; x <- seq(quantile(sub[[v]], .01), quantile(sub[[v]], .99), length.out = 80)
  nd <- as.data.frame(tp[tp$type == ty, ])[rep(1, length(x)), ]; nd$type <- factor(ty, levels = levels(P$type)); nd[[v]] <- x
  nd$.x <- x; nd
}))

# difference in slope between pair types: within alpine minus within lowland, alpine vs lowland minus within lowland, and minus within alpine
type_compare <- list(c(2, 1), c(3, 1), c(3, 2))                 # indices into beta_types: first minus second

# ---- all replicates: the random pairs and their dissimilarities. Saved (beta/pairs_all.rds) and reused on a rerun if the
# settings and the drawn samples are the same, so a problem in the models does not cost the distance computation again.
pairs_key <- list(dist_breaks = CFG$dist_breaks_km, elev_breaks = CFG$elev_breaks_km, n_rep = CFG$n_rep, n_alpine_max = CFG$n_alpine_max,
                  n_pairs = CFG$n_pairs, rare_depth = CFG$rare_depth, n_ids = length(all_ids), ids_chars = sum(nchar(all_ids)), use_sets = use_sets,
                  thin_m = CFG$thin_m, even_cell_deg = CFG$even_cell_deg, test = TEST_MODE, design = "one-per-cluster-v2")
pairs_file <- out_b("pairs_all.rds")
P_all <- NULL
if (REUSE_FITS && file.exists(pairs_file)) {
  cached <- readRDS(pairs_file)
  if (identical(cached$key, pairs_key)) { P_all <- cached$pairs; message("Reusing saved pairs: ", pairs_file) }
}
if (is.null(P_all)) {
  beta_pairs <- list()
  for (r in seq_len(CFG$n_rep)) {
    message("Replicate ", r, ": computing dissimilarities")
    pr_r <- make_pairs(r)
    if (is.null(pr_r)) { message("  no pairs drawn: skipping replicate"); next }
    beta_pairs[[length(beta_pairs) + 1]] <- pr_r
  }
  P_all <- bind_rows(beta_pairs); rm(beta_pairs)
  saveRDS(list(key = pairs_key, pairs = P_all), pairs_file)
}

P <- prep_pairs(P_all)
message(sprintf("%d pairs from %d replicates in %d pair types", nrow(P), length(unique(P$rep_id)), nlevels(P$type)))
print(table(type = P$type))
write.csv(P %>% count(rep_id, type), out_b("pair_counts.csv"), row.names = FALSE)
# share of pairs with dissimilarity exactly 1 (no shared OTU): the pile that makes a beta regression fail
at_one <- bind_rows(lapply(names(beta_metrics), function(mt)
  P %>% group_by(type) %>% summarise(share_at_1 = mean(.data[[mt]] >= 1 - 1e-12), n_pairs = n(), .groups = "drop") %>% mutate(metric = mt)))
write.csv(at_one, out_b("share_of_pairs_at_dissimilarity_1.csv"), row.names = FALSE)
for (i in which(at_one$metric == "bray")) message(sprintf("  %s: %.1f%% of the Bray-Curtis pairs are exactly 1", at_one$type[i], 100 * at_one$share_at_1[i]))
beta_raw <- P[, c("type", "d_elev", "d_lat", "log10_geo", names(beta_metrics))]     # every pair (for the plots)

# ---- the pooled regressions: lines with their 95% interval, the adjusted observed values and a check of every fit
tp <- typ(P)
curve_cache <- list(); for (v in names(beta_vars)) curve_cache[[v]] <- curve_nd(tp, P, v)
curves <- list(); adj <- list(); checks <- list(); fits <- list()
for (metric in names(beta_metrics)) for (kind in c("elev", "lat")) {
  f <- beta_formula(metric, kind); message("  fitting ", metric, " / ", kind, ": ", deparse(f), " (", format(nrow(P), big.mark = ","), " pairs)")
  m <- tryCatch(lm(f, data = P), error = function(e) { message("    fit failed: ", conditionMessage(e)); NULL }); if (is.null(m)) next
  V <- dyadic_vcov(model.matrix(m), resid(m), P$ca, P$cb); if (is.null(V)) { message("    the covariance matrix could not be computed"); next }
  key <- paste(metric, kind); fits[[key]] <- list(m = m, V = V, metric = metric, kind = kind)
  for (v in names(curve_model)[curve_model == kind]) {
    nd <- curve_cache[[v]]; Xg <- design_at(m, nd); fi <- drop(Xg %*% coef(m)); se <- sqrt(rowSums((Xg %*% V) * Xg))
    curves[[paste(key, v)]] <- data.frame(metric = metric, variable = v, type = as.character(nd$type), x = nd$.x, fit = clip01(fi), lo = clip01(fi - 1.96 * se), hi = clip01(fi + 1.96 * se))
    # the observed dissimilarity of every pair adjusted to the typical values of the other variables (residual + the line at the pair's x):
    # binned, these are the open circles that the line has to run through
    nd_i <- as.data.frame(tp[match(as.character(P$type), as.character(tp$type)), ]); nd_i$type <- P$type; nd_i[[v]] <- P[[v]]
    adj[[paste(key, v)]] <- data.frame(metric = metric, variable = v, type = as.character(P$type), x = P[[v]], y = P[[metric]] - fitted(m) + predict(m, nd_i))
  }
  fv <- tapply(fitted(m), P$type, mean); ob <- tapply(P[[metric]], P$type, mean)       # does the fit reproduce the observed means?
  checks[[key]] <- data.frame(metric = metric, model = kind, type = names(ob), observed_mean = as.numeric(ob), fitted_mean = as.numeric(fv[names(ob)]), n_pairs = as.numeric(table(P$type)[names(ob)]))
}
if (!length(fits)) stop("no beta model could be fitted: see the messages above")
checks <- bind_rows(checks); checks$difference <- checks$fitted_mean - checks$observed_mean
write.csv(checks, out_b("model_checks.csv"), row.names = FALSE)
worst <- max(abs(checks$difference), na.rm = TRUE)
message(sprintf("Model check: the fitted mean dissimilarity differs from the observed one by at most %.3f (largest over %d fits x pair types)%s", worst, nrow(checks),
                if (worst > 0.02) "  -> LARGER THAN 0.02: look at beta/model_checks.csv, a fit may be wrong" else ""))
preds_o <- bind_rows(curves); adj_all <- bind_rows(adj)
write.csv(preds_o, out_b("predicted_curves_overall.csv"), row.names = FALSE)

# ---- the tests: coefficients, slopes per pair type and differences between pair types (t-tests, clusters - 1 degrees of freedom)
G_cl <- attr(fits[[1]]$V, "G"); dfree <- max(G_cl - 1, 1)
message(sprintf("Tests: %s pairs from %s clusters (the independent units of the standard errors)", format(nrow(P), big.mark = ","), format(G_cl, big.mark = ",")))
tfun <- function(est, se) data.frame(estimate = est, se = se, lo = est - qt(0.975, dfree) * se, hi = est + qt(0.975, dfree) * se, p = 2 * pt(-abs(est / se), dfree))
coef_all <- bind_rows(lapply(fits, function(ft) { m <- ft$m; b <- coef(m)
  cbind(metric = ft$metric, model = ft$kind, term = names(b), tfun(unname(b), sqrt(diag(ft$V))), p_naive_lm_NOT_VALID = unname(summary(m)$coefficients[names(b), 4]), row.names = NULL) }))
write.csv(coef_all, out_b("model_coefficients.csv"), row.names = FALSE)
slope_tab <- bind_rows(lapply(fits, function(ft) { b <- coef(ft$m); idx <- grep(":(d_elev|d_lat|log10_geo)$", names(b)); nm <- names(b)[idx]
  cbind(metric = ft$metric, model = ft$kind, group = sub(":(d_elev|d_lat|log10_geo)$", "", sub("^type", "", nm)), variable = sub("^.*:", "", nm),
        tfun(unname(b[idx]), sqrt(diag(ft$V)[idx])), p_naive_lm_NOT_VALID = unname(summary(ft$m)$coefficients[nm, 4]), row.names = NULL) }))
write.csv(slope_tab, out_b("linear_effects.csv"), row.names = FALSE)               # slopes (dissimilarity per unit of the variable) per pair type
slope_con <- bind_rows(lapply(fits, function(ft) {
  b <- coef(ft$m); nm <- names(b); V <- ft$V
  vars <- unique(sub("^.*:", "", nm[grepl(":", nm)]))
  bind_rows(lapply(type_compare, function(cp) bind_rows(lapply(vars, function(vr) {
    i1 <- paste0("type", beta_types[cp[1]], ":", vr); i2 <- paste0("type", beta_types[cp[2]], ":", vr)
    if (!all(c(i1, i2) %in% nm)) return(NULL)
    v <- setNames(numeric(length(nm)), nm); v[i1] <- 1; v[i2] <- -1
    cbind(metric = ft$metric, model = ft$kind, variable = vr, g1 = beta_types[cp[1]], g2 = beta_types[cp[2]], tfun(sum(v * b), sqrt(drop(t(v) %*% V %*% v))), row.names = NULL)
  }))))
}))
write.csv(slope_con, out_b("linear_effect_type_contrasts.csv"), row.names = FALSE)    # differences between the slopes of two pair types
pfmt <- function(p) format.pval(p, digits = 2, eps = 1e-16)
for (mt in names(beta_metrics)) for (kind in c("elev", "lat")) {
  d <- slope_tab[slope_tab$metric == mt & slope_tab$model == kind & slope_tab$variable == c(elev = "d_elev", lat = "d_lat")[[kind]], ]
  if (nrow(d)) message(sprintf("  %s, slope of %s per %s: %s", beta_metrics[[mt]], c(elev = "the elevation difference", lat = "the latitude difference")[[kind]], c(elev = "km", lat = "degree")[[kind]],
                               paste(sprintf("%s %+.4f (p %s)", d$group, d$estimate, pfmt(d$p)), collapse = "; ")))
}

## ---- B4. Beta-diversity figures ----------------------------------------------------------------------------------------
# Zoom of the y axis: dissimilarities pile up near 1 with a long tail towards 0, which flattens the curves. The window starts
# at the lower of the beta_zoom_q quantile of the panel's pairs and the lowest fitted value (minus a margin): pairs below it are
# not drawn (they are in the models) and a ribbon that runs below it is clipped there, so one badly determined ribbon cannot
# undo the zoom. `by` = the columns that define one panel. Returns the points and the curves to draw.
zoom_panel <- function(ob, pm, by) {
  if (!isTRUE(CFG$beta_zoom) || !nrow(ob)) return(list(ob = ob, pm = pm))
  lim <- ob %>% group_by(across(all_of(by))) %>% summarise(q = quantile(y, CFG$beta_zoom_q), .groups = "drop") %>%
    left_join(pm %>% group_by(across(all_of(by))) %>% summarise(f = min(fit), .groups = "drop"), by = by) %>%
    mutate(edge = pmin(q, f, na.rm = TRUE) - 0.01)
  out <- ob %>% left_join(lim[, c(by, "edge")], by = by) %>% filter(y >= edge) %>% select(-edge)
  pm  <- pm %>% left_join(lim[, c(by, "edge")], by = by) %>% mutate(lo = pmax(lo, edge)) %>% select(-edge)
  message(sprintf("    y-axis zoom: %d of %d points below the window are not drawn", nrow(ob) - nrow(out), nrow(ob)))
  list(ob = out, pm = pm)
}
y_breaks <- function(l) { b <- scales::extended_breaks()(l); b[b <= 1 + 1e-9] }                 # dissimilarities cannot exceed 1
y_scale  <- scale_y_continuous(expand = expansion(mult = c(0.03, 0.06)), breaks = y_breaks)
thin_txt <- sprintf("Each replicate: one random alpine sample per cluster (BioProject / location) and the same number of non-alpine samples, evenly spread over latitude / longitude and at most one per %g m", CFG$thin_m)
pts_note <- if (CFG$beta_points_max > 0) sprintf("points = a random subset of %s pairs (the y axis is zoomed on the bulk of the pairs: the few with a much lower dissimilarity are not drawn and ribbons are clipped at the lower edge, but all pairs are in the models)", format(CFG$beta_points_max, big.mark = ",")) else "pairs not drawn"
rob_txt <- "ordinary linear regression on the pairs of all replicates (lines = regression, ribbons = 95% interval, standard errors cluster-robust for pair data: pairs sharing a cluster are dependent), open circles = observed dissimilarity adjusted to the typical values of the other variables (mean in 8 bins of the x variable); stars = tests of the slope of each pair type and of the differences between two slopes (L = within lowland, A = within alpine, X = alpine vs lowland; *** p < 0.001, ** p < 0.01, * p < 0.05, ns; not adjusted for multiple tests)"
lab_beta <- function(d) {
  d$metric   <- factor(beta_metrics[d$metric], levels = beta_metrics)
  d$variable <- factor(beta_vars[d$variable], levels = beta_vars)
  d
}
type_scales <- list(scale_colour_manual(values = CFG$cols_beta, name = NULL), scale_fill_manual(values = CFG$cols_beta, guide = "none"),
                    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1, linewidth = 1.2, linetype = "solid"))))
fac_type <- function(d) { d$type <- factor(d$type, levels = beta_types); d }
# random subset of the pairs for drawing (all pairs are in the models)
pair_points <- function(vars, metrics, max_n) {
  ob <- bind_rows(lapply(metrics, function(mt) bind_rows(lapply(vars, function(v)
    data.frame(metric = mt, variable = v, type = beta_raw$type, x = beta_raw[[v]], y = beta_raw[[mt]])))))
  set.seed(5)
  if (max_n <= 0) ob[0, ] else if (nrow(ob) > max_n) ob[sample.int(nrow(ob), max_n), ] else ob
}

# the slopes of the three pair types and the differences between them as two lines of text for a panel (elevation: per km; latitude: per 10 degrees)
beta_annot <- function(metric, v) {
  kind <- curve_model[[v]]
  s <- slope_tab[slope_tab$metric == metric & slope_tab$model == kind & slope_tab$variable == v, ]
  d <- slope_con[slope_con$metric == metric & slope_con$model == kind & slope_con$variable == v, ]
  if (!nrow(s)) return("")
  fac <- if (v == "d_lat") 10 else 1; ts <- c("Within lowland" = "L", "Within alpine" = "A", "Alpine vs lowland" = "X")
  paste(paste0("Slope per ", if (v == "d_lat") "10 degrees" else "km", ":  ", paste(sprintf("%s %+.3f %s", ts[s$group], fac * s$estimate, stars(s$p)), collapse = "   ")),
        if (nrow(d)) paste0("Differences:  ", paste(sprintf("%s-%s %+.3f %s", ts[d$g1], ts[d$g2], fac * d$estimate, stars(d$p)), collapse = "   ")), sep = "\n")
}

# One panel of the effects: one metric against one pair variable, a line per pair type
beta_overall_panel <- function(metric, v, tag, subtitle) {
  pm <- preds_o[preds_o$metric == metric & preds_o$variable == v, ]
  if (!nrow(pm)) return(NULL)
  ob <- pair_points(v, metric, CFG$beta_points_max); pm <- fac_type(pm); ob <- fac_type(ob)
  # observed dissimilarity adjusted to the typical values of the other variables, mean in 8 bins of x per pair type: the curves must run through these circles
  bn <- adj_all[adj_all$metric == metric & adj_all$variable == v, ] %>% group_by(type) %>% mutate(bin = ntile(x, 8)) %>%
    group_by(type, bin) %>% summarise(x = mean(x), y = mean(y), .groups = "drop") %>% as.data.frame() %>% fac_type()
  pm$mid <- pm$fit; ob$.p <- 1; pm$.p <- 1
  if (nrow(ob)) { z <- zoom_panel(ob, pm, ".p"); ob <- z$ob; pm <- z$pm }
  ggplot() +
    geom_ribbon(data = pm, aes(x, ymin = lo, ymax = hi, fill = type), alpha = 0.18) +
    geom_point(data = ob, aes(x, y, colour = type), size = 0.35, alpha = CFG$beta_point_alpha, shape = 16) +
    geom_line(data = pm, aes(x, mid, colour = type), linewidth = 1.1) +
    geom_point(data = bn, aes(x, y, colour = type), shape = 21, fill = "white", size = 2.6, stroke = 1.1) +
    annotate("label", x = -Inf, y = -Inf, label = beta_annot(metric, v), hjust = -0.03, vjust = -0.15, size = 2.7, lineheight = 1.1, fill = "white", alpha = 0.8, linewidth = 0) +
    type_scales + y_scale +
    labs(tag = tag, x = beta_vars[[v]], y = paste(beta_metrics[[metric]], "dissimilarity"), subtitle = subtitle) +
    theme_fig + theme(plot.subtitle = element_text(size = 9, colour = "grey30"))
}

# Full grid (all metrics x all pair variables): beta/
{
  pm <- fac_type(lab_beta(preds_o)); ob <- fac_type(lab_beta(pair_points(names(beta_vars), names(beta_metrics), 3 * CFG$beta_points_max)))
  pm$mid <- pm$fit
  if (nrow(ob)) { z <- zoom_panel(ob, pm, c("metric", "variable")); ob <- z$ob; pm <- z$pm }
  g <- ggplot() +
    geom_ribbon(data = pm, aes(x, ymin = lo, ymax = hi, fill = type), alpha = 0.18) +
    geom_point(data = ob, aes(x, y, colour = type), size = 0.3, alpha = CFG$beta_point_alpha, shape = 16) +
    geom_line(data = pm, aes(x, mid, colour = type), linewidth = 1.1) +
    facet_grid(metric ~ variable, scales = "free", switch = "both") +
    type_scales + y_scale + labs(x = NULL, y = "Dissimilarity") +
    theme_bw(base_size = 12) +
    theme(strip.placement = "outside", strip.background = element_blank(), strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank(), legend.position = "bottom")
  ggsave(out_b("beta_overall_grid.png"), g, width = 9.5, height = 7.5, dpi = 300)
}
# Supplementary Figure 1: Jaccard (rarefied presence/absence)
pS1A <- beta_overall_panel("jaccard", "d_elev", "A", "accounting for geographic distance")
pS1B <- beta_overall_panel("jaccard", "d_lat", "B", "accounting for the elevation difference")
if (!is.null(pS1A) && !is.null(pS1B)) {
  s1_caption <- paste("Supplementary Figure 1. Jaccard (rarefied presence/absence) dissimilarity between pairs of samples against the difference in elevation (A; accounting for geographic distance)",
                      "and in latitude (B; accounting for the elevation difference); colour = pair type (within non-alpine, within alpine, between).", paste0(rob_txt, "."),
                      paste0(thin_txt, "."), paste0(pts_note, "."))
  assemble_fig(file.path(out_folder, "figS1_beta_jaccard.png"), list(list(pS1A, pS1B)), get_legend(pS1A), s1_caption, 14, 6.5)
}

}  # end RUN_BETA

## ======================= MAIN FIGURE (needs both parts) =====================================
# A: alpha diversity along elevation (row 1); B: along latitude (row 2); C: Bray-Curtis dissimilarity against the elevation difference
# (accounting for geographic distance) and the latitude difference (accounting for elevation), per pair type: one alpine sample per
# cluster against the same number of evenly spread non-alpine samples, in each of the replicates
if (RUN_ALPHA && RUN_BETA) {
  pA <- alpha_panel(grad_lab[["elevation"]], "A")
  pB <- alpha_panel(grad_lab[["glat"]], "B")
  pC <- beta_overall_panel("bray", "d_elev", "C", "accounting for geographic distance")
  pD <- beta_overall_panel("bray", "d_lat", "", "accounting for the elevation difference")
  if (is.null(pC) || is.null(pD)) stop("main figure: the Bray-Curtis models are missing")
  main_caption <- paste(
    "A, B: alpha diversity (observed OTUs, Shannon index, Pielou J) along elevation (A) and latitude (B); GAM partial effects",
    "adjusted for sequencing depth and unmapped proportion, solid = net of space (95% CI), dashed = no spatial term; stars = test that the alpine gradient differs in shape from the",
    "lowland one (net model; *** p < 0.001, ** p < 0.01, * p < 0.05, ns).",
    "C: Bray-Curtis dissimilarity between pairs of samples (colour = pair type: within non-alpine, within alpine, between alpine and non-alpine)",
    "against the difference in elevation (left; accounting for geographic distance) and the difference in latitude (right; accounting for the elevation difference): one effect per pair type.",
    sprintf("%s; up to %d pairs per pair type and replicate (random within classes of geographic distance and elevation difference), %d replicates.", thin_txt, CFG$n_pairs, CFG$n_rep),
    paste0("Effects: ", rob_txt, "."),
    paste0(pts_note, ". "), alpha_txt,
    "Jaccard: Supplementary Figure 1.")
  assemble_fig(file.path(out_folder, "fig2_alpha_beta_gradients.png"),
               list(list(pA), list(pB), list(get_legend(pA)), list(pC, pD)), get_legend(pC), main_caption,
               width = 14, height = 15, row_heights = c(1, 1, 0.17, 1.3))
  saveRDS(list(A = pA, B = pB, C = pC, D = pD), file.path(out_folder, "fig2_panels.rds"))
} else if (RUN_ALPHA) {                                  # no beta part (sensitivity analyses): panels A and B only
  pA <- alpha_panel(grad_lab[["elevation"]], "A")
  pB <- alpha_panel(grad_lab[["glat"]], "B")
  ab_caption <- paste(sprintf("Alpha diversity (observed OTUs, Shannon index, Pielou J) along %s (axes: gradient_axes = '%s').", axes_txt, CFG$gradient_axes),
                      "A, B: GAM partial effects adjusted for sequencing depth and unmapped proportion, solid = net of space (95% CI), dashed = no spatial term; stars = test that the alpine",
                      "gradient differs in shape from the lowland one (net model; *** p < 0.001, ** p < 0.01, * p < 0.05, ns).", alpha_txt)
  assemble_fig(file.path(out_folder, "fig2AB_alpha_gradients.png"), list(list(pA), list(pB)), get_legend(pA), ab_caption,
               width = 14, height = 8.5, row_heights = c(1, 1))
  saveRDS(list(A = pA, B = pB), file.path(out_folder, "fig2AB_panels.rds"))
  message("The beta part was not run: only panels A and B were drawn (fig2AB_alpha_gradients.png)")
} else message("Main figure skipped: it needs RUN_ALPHA and RUN_BETA both TRUE.")

message("Done. Outputs in ", out_folder)
writeLines(capture.output(sessionInfo()), file.path(out_folder, "sessionInfo.txt"))
file.create(file.path(out_folder, ".done"))   # flag file for Snakemake