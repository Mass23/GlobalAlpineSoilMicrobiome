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
#   * SAMPLE POOLS. Figure 2 uses ALL samples: the gradient models condition on elevation, latitude and distance, so an
#     unbalanced design does not bias the curves. Samples that share a location are not independent, so ALPHA models keep every
#     alpine sample and give each location (samples within thin_m metres) a random intercept, plus a spherical spline for the
#     regional trend; the lowland group (data-rich) is thinned to one sample per location (a random effect with ~100k levels is
#     not feasible in mgcv; see alpha_cluster_*). BETA pairs are drawn from one random sample per location in each replicate:
#     different replicates draw different samples of a location, so all samples are used across the replicates. The three
#     fixed lowland sets of
#     define_lowland_sets.py (random / spatially even / elevation representation, as many lowland samples as alpine ones, kept
#     for the soil data) would discard most of the lowland data. The even coverage of the gradients comes from the way the
#     pairs are drawn instead (next point). The three sets are for the analyses that need soil data or compare alpine and
#     lowland samples as groups: Figure 3 and 5 (elevation set; headline), Figure 4 and 5 (spatially even set), all
#     with the random set as the baseline check (flags: subset_* columns of all_sampled_data.csv).
#   * Pairs are drawn at random within geographic distance classes (dist_breaks_km: short distances would hardly be
#     represented by random pairs) crossed with elevation-difference classes (elev_breaks_km: lowland pairs hardly differ in
#     elevation otherwise), shared as evenly as possible over the classes, and WITHIN latitude zones (lat_breaks / lat_zone_names; default: southern hemisphere, northern
#     mid-latitudes 0-60, northern high latitudes > 60) and within three pair types: within lowland, within alpine,
#     alpine vs lowland. A "cell" = pair type x zone. Samples are spatially thinned (at most one per thin_m metres, and
#     no shared coordinates) and every replicate draws n_per_group samples per group and zone and n_pairs pairs per cell.
#   * Bray-Curtis (relative abundances) and Jaccard (presence/absence after rarefaction to rare_depth reads), computed
#     only for the drawn pairs from sparse matrices (checked against vegan::vegdist in TEST_MODE).
#   * Beta-regression GAMs per metric, each on the pooled pairs of all n_rep replicates with a random intercept for the
#     replicate, an intercept per cell and LINEAR gradient effects (beta_linear; logit scale; elevation difference in km,
#     latitude difference in degrees, distance in log10 km), smooth terms for the differences in depth and unmapped proportion:
#       elev: elevation difference + log10 geographic distance      -> elevation ACCOUNTING FOR SPACE (and the distance decay)
#       lat : latitude difference + elevation difference            -> latitude ACCOUNTING FOR ELEVATION
#     at two levels: OVERALL = one effect per pair type, common to the latitude zones (an intercept per cell, so the effect is
#     estimated within zones): the main figure; BY ZONE = one effect per pair type and zone (a "cell"): Supplementary Figure 4,
#     the zone comparisons and the standardised dissimilarity. Slopes with intervals: beta/linear_effects.csv, zone differences
#     of the slopes: beta/linear_effect_zone_contrasts.csv.
#   * STANDARDISED DISSIMILARITY (beta/standardized_dissimilarity.csv, standardized_contrasts.csv, Figure 2b): from the
#     elev model, the dissimilarity two samples of a cell are predicted to have at the SAME geographic distance (ref_km),
#     the typical elevation difference of the pair type and median depth / unmapped differences. This removes what makes
#     regions incomparable otherwise (spatial extent, elevation range, sequencing), so that zones and pair types can be
#     compared: e.g. alpine pairs in the northern high latitudes vs the northern mid-latitudes, or the southern vs the
#     northern hemisphere. The contrasts come with 95% intervals.
#   CAVEAT: pairs reuse samples (within a replicate each sample is in ~2 of the n_pairs pairs, and a sample can be drawn
#   in several replicates). The random effect for the replicate does not remove that dependence, so the intervals are
#   descriptive and optimistic; judge robustness by the size of the effects and their consistency.
#
# FIGURES (in <out_folder>)
#   fig2_alpha_beta_gradients.png : A = alpha diversity vs elevation (3 metrics), B = vs latitude (signed),
#                                   C = Bray-Curtis vs elevation difference (accounting for space) and vs latitude
#                                   difference (accounting for elevation): overall effects per pair type
#   figS4_beta_by_zone.png        : Supplementary Figure 4 = the same by latitude zone (and the distance decay)
#   fig2b_standardized_dissimilarity.png : standardised dissimilarity per zone and pair type
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
  lat_breaks  = c(0, 60),       # latitude zones (degrees; beta pairs are only drawn WITHIN a zone): -90 to 0, 0 to 60, 60 to 90
  lat_zone_names = c("Southern hemisphere", "Northern mid-latitudes", "Northern high latitudes"),   # one name per zone
  dist_breaks_km = c(0, 10, 100, 1000, 5000, Inf),   # pairs are drawn at random WITHIN these geographic distance classes (km) ...
  elev_breaks_km = c(0, 0.1, 0.5, 1, 2, Inf),        # ... crossed with these classes of the elevation difference (km): every elevation
                                                     #   difference is represented, also among the lowland pairs (equal numbers per class)
  ref_km      = c(100, 1000),   # standardised dissimilarity: reference geographic distances (km)
  ref_zone    = "Northern mid-latitudes",   # zone the other zones are compared with in standardized_contrasts.csv (else: the zone with most pairs)
  thin_m      = 100,            # a LOCATION = samples within this radius (metres). Beta subsamples: at most ONE sample per location and replicate
                                #   (metres, chosen at random); 0 = only distinct coordinates (exact duplicates removed)
  alpha_cluster_alpine  = "re",   # alpha models (part A): samples within thin_m metres of each other share a location. ALPINE group:
                                  #   "re" = keep every sample, random intercept per location; "thin" = one sample per location; "none"
  alpha_cluster_lowland = "thin", # LOWLAND group: "thin" = one sample per location (a data-rich group: a random effect with ~100k levels
                                  #   is not feasible in mgcv); "re" = random intercept per re_cell_deg cell; "none"
  re_cell_deg   = 1,              # size (degrees) of the cells of the lowland random effect (alpha_cluster_lowland = "re")
  max_re_levels = 3000,           # safety: a random effect may have at most this many levels (the cost grows steeply with the levels)
  k_space     = 60,             # basis dimension of the sos spatial spline
  k_space_alpine = 30,          # Supplementary Figure 3: basis dimension of the alpine-specific surface
  grid_res    = 2,              # Supplementary Figure 3: grid resolution (degrees)
  min_support = 5,              # Supplementary Figure 3: min. samples within the 3x3 cell neighbourhood to show a cell
  map_ylim    = c(-60, 83),     # Supplementary Figure 3: latitude window of the maps
  n_rug       = 3000,           # rug marks drawn per group (plot only)
  moran_n_max = 1500, moran_nperm = 199,
  # --- part B ---
  n_per_group = 1000,           # random samples per group AND latitude zone per replicate
  n_rep       = 20,             # replicate random subsamples (each re-drawn from scratch), pooled with a random effect
  rare_depth  = 1000,           # rarefaction depth for Jaccard (samples below are dropped)
  n_pairs     = 1000,           # pairs per pair type, latitude zone and replicate (shared as evenly as possible over the distance x elevation classes)
  min_pairs   = 100,            # cells (pair type x zone) with fewer pairs are dropped from the models
  k_smooth    = 5,              # basis dimension of every 1-D smooth (alpha and beta models): small = smooth curves
  bs_smooth   = "ts",           # thin plate spline WITH shrinkage (the null space is penalised too)
  select      = FALSE,          # TRUE would add a second shrinkage penalty on top of "ts"
  beta_linear = TRUE,           # beta gradient effects (elevation difference, latitude difference, log10 distance) are LINEAR (FALSE: smooth, bs = "ts")
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
  CFG$n_per_group <- 200; CFG$n_rep <- 1; CFG$n_pairs <- 300
  out_folder <- paste0(out_folder, "_TEST")
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

# Latitude zones of the beta analysis: pairs are only drawn within a zone (never across zones)
zone_edges  <- c(-90, CFG$lat_breaks, 90)
zone_range  <- sprintf("%g to %g", head(zone_edges, -1), tail(zone_edges, -1))
if (length(CFG$lat_zone_names) != length(zone_range))
  stop("CFG$lat_zone_names needs ", length(zone_range), " names (one per zone: ", paste(zone_range, collapse = ", "), ")")
zone_names  <- CFG$lat_zone_names
zone_desc   <- sprintf("%s: %s degrees", zone_names, zone_range)
lat_zone    <- function(lat) as.character(cut(lat, breaks = c(-Inf, CFG$lat_breaks, Inf), labels = zone_names, right = FALSE))
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
env <- data.frame(.env_id = env[[CFG$env_id]], .elev_m = env[[CFG$env_elev]])
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

# Samples that share a location (within CFG$thin_m metres: one study / site) are not independent: left alone they make the curves and
# p-values overconfident. Per group: "re" = keep every sample and give each location a random intercept (the effective number of
# independent units is then the number of locations); "thin" = keep one random sample per location; "none" = nothing.
alpha_how <- c(Alpine = CFG$alpha_cluster_alpine, Lowland = CFG$alpha_cluster_lowland)
if (!all(alpha_how %in% c("re", "thin", "none"))) stop('CFG$alpha_cluster_alpine / _lowland must be "re", "thin" or "none"')
set.seed(7); n0 <- nrow(d); d$re_loc <- NA_character_; keep_rows <- rep(TRUE, nrow(d))
for (g in levels(d$Alpine)) {
  ix <- which(d$Alpine == g)
  if (alpha_how[[g]] == "thin") {
    keep_rows[ix] <- FALSE; keep_rows[ix[thin_sites(d$lat[ix], d$lon[ix], CFG$thin_m)]] <- TRUE
  } else if (alpha_how[[g]] == "re") {
    id <- if (g == "Alpine") cluster_sites(d$lat[ix], d$lon[ix], CFG$thin_m) else {
      cs <- CFG$re_cell_deg
      repeat {                                      # coarsen the cells until the random effect is affordable
        key <- paste(floor(d$lat[ix] / cs), floor(d$lon[ix] / cs))
        if (length(unique(key)) <= CFG$max_re_levels || cs > 40) break
        cs <- cs * 2
      }
      if (cs != CFG$re_cell_deg) message(sprintf("  lowland random effect: cells coarsened to %g degrees (max_re_levels = %d)", cs, CFG$max_re_levels))
      match(key, unique(key))
    }
    d$re_loc[ix] <- paste0(substr(g, 1, 1), id)
  }
}
d <- d[keep_rows, ]
d$re_ind <- as.numeric(!is.na(d$re_loc))             # 1 = this row has a location random effect, 0 = it has none (thinned or "none")
use_re <- sum(d$re_ind) > 0
d$re_loc <- if (use_re) { f <- factor(d$re_loc); f[is.na(f)] <- levels(f)[1]; f } else factor(rep("none", nrow(d)))
if (use_re && nlevels(d$re_loc) > CFG$max_re_levels)
  stop(sprintf("the location random effect has %d levels (> max_re_levels = %d): thin that group (alpha_cluster_* = \"thin\") or raise max_re_levels",
               nlevels(d$re_loc), CFG$max_re_levels))
re_rhs <- if (use_re) "+ s(re_loc, bs = 're', by = re_ind)" else ""         # label: s(re_loc):re_ind
re_labels <- "s(re_loc):re_ind"
message(sprintf("Alpha data: %d -> %d samples (alpine %d, lowland %d); location random effect with %d levels (alpine: %s, lowland: %s)",
                n0, nrow(d), sum(d$Alpine == "Alpine"), sum(d$Alpine == "Lowland"), if (use_re) nlevels(d$re_loc) else 0L,
                alpha_how[["Alpine"]], alpha_how[["Lowland"]]))
how_txt <- function(g, noun) switch(alpha_how[[g]],
  re   = sprintf("all %s samples with a random intercept per %s", tolower(g), noun),
  thin = sprintf("%s samples thinned to at most one per %g m", tolower(g), CFG$thin_m),
  none = sprintf("all %s samples", tolower(g)))
alpha_txt <- sprintf("A, B use %s and %s, and a spherical spline for the regional trend.", how_txt("Alpine", sprintf("location (samples within %g m)", CFG$thin_m)),
                     how_txt("Lowland", sprintf("%g-degree cell", CFG$re_cell_deg)))

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
# how strong is the clustering of samples within locations? standard deviation of the location random intercept (link scale)
if (use_re) {
  re_sd <- bind_rows(lapply(names(fits), function(nm) {
    sd_loc <- tryCatch({ vc <- suppressWarnings(gam.vcomp(fits[[nm]], rescale = FALSE))
                         if (is.list(vc) && !is.matrix(vc)) vc <- vc$vc
                         as.numeric(if (is.matrix(vc)) vc[re_labels, 1] else vc[re_labels]) }, error = function(e) NA_real_)
    data.frame(model = nm, location_sd = if (length(sd_loc) == 1) sd_loc else NA_real_)
  }))
  write.csv(re_sd, out_a("location_random_effect_sd.csv"), row.names = FALSE)
  message(sprintf("Location random effect: SD %s (link scale; per model in alpha/location_random_effect_sd.csv)",
                  paste(sprintf("%s %.2f", re_sd$model[grepl("__net$", re_sd$model)], re_sd$location_sd[grepl("__net$", re_sd$model)]), collapse = ", ")))
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
                     re_loc = factor(levels(d$re_loc)[1], levels = levels(d$re_loc)), re_ind = 0,
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

grad_lab <- c(elevation = "Elevation (km)",
              glat = if (CFG$abs_lat) "Absolute latitude (\u00B0)" else "Latitude (\u00B0)")
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
  show_zones <- gr_label == grad_lab[["glat"]] && !isTRUE(CFG$abs_lat) && length(CFG$lat_breaks) > 0   # zone boundaries of the beta part
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
    (if (show_zones) geom_vline(xintercept = CFG$lat_breaks, linetype = 3, colour = "grey50", linewidth = 0.35)) +
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
  nd$re_loc <- factor(levels(d$re_loc)[1], levels = levels(d$re_loc)); nd$re_ind <- 0        # no location random effect on the map
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
    unmapped  = unmapped_raw_counts / total_raw_counts
  ) %>%
  filter(is.finite(lon), is.finite(lat), abs(lon) <= 180, abs(lat) <= 90,
         is.finite(elevation), is.finite(log10_reads), is.finite(unmapped), !is.na(Alpine), !is.na(id))
print(table(group = meta$Alpine))
meta$zone <- lat_zone(meta$lat)
print(table(zone = meta$zone, group = meta$Alpine))

## ---- B1. Replicate subsamples and counts -------------------------------------------------------
# Each replicate draws, within each group AND latitude zone, up to n_per_group samples in random order such that no two of
# them are closer than CFG$thin_m metres (and none share coordinates): the pairs then come from different places, not from
# one site or study, and every zone is represented by a comparable number of samples
set.seed(99)
for (g in levels(meta$Alpine)) for (z in zone_names) {
  ix <- which(meta$Alpine == g & meta$zone == z)
  if (!length(ix)) { message(sprintf("  %s / %s: no samples", g, z)); next }
  capped <- length(ix) > 20000                        # counting all the thinned samples is only done for smaller sets
  n_av <- length(thin_sites(meta$lat[ix], meta$lon[ix], CFG$thin_m, n_max = if (capped) CFG$n_per_group else Inf))
  message(sprintf("  spatial thinning (one sample per %g m), %s / %s: %s%d of %d samples can be drawn%s", CFG$thin_m, g, z,
                  if (capped) "at least " else "", n_av, length(ix),
                  if (n_av < CFG$n_per_group) sprintf(" (fewer than n_per_group = %d: each replicate takes all of them)", CFG$n_per_group) else ""))
}
sel <- lapply(seq_len(CFG$n_rep), function(r) {
  set.seed(100 + r)
  bind_rows(lapply(levels(meta$Alpine), function(g) bind_rows(lapply(zone_names, function(z) {
    mg <- meta[meta$Alpine == g & meta$zone == z, ]
    if (nrow(mg) == 0) return(NULL)
    mg[thin_sites(mg$lat, mg$lon, CFG$thin_m, n_max = CFG$n_per_group), ]
  }))))
})
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

  # up to n_pairs pairs (indices into S), drawn at random WITHIN strata = geographic distance class x elevation-difference class
  # (CFG$dist_breaks_km, CFG$elev_breaks_km), the n_pairs shared as evenly as possible over the strata that exist. Random pairs are
  # mostly thousands of km apart and, among lowland samples, hardly differ in elevation: without stratification short distances
  # and large elevation differences would hardly be represented. Stratifying on predictors of the models does not bias the fitted
  # curves; it makes them well determined over the whole range. The natural (unstratified) median distance and elevation difference
  # of the cell are kept (nat_geo, nat_delev): they are the "typical" values at which other variables are held in the figures.
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
    D <- geosphere::distm(cbind(S$lon[i1], S$lat[i1]), cbind(S$lon[i2], S$lat[i2]), fun = geosphere::distHaversine) / 1000   # km
    E <- abs(outer(S$elevation[i1], S$elevation[i2], "-"))                                                                     # km
    idx <- if (within) which(upper.tri(D), arr.ind = TRUE) else which(is.finite(D), arr.ind = TRUE)
    d <- D[idx]; e <- E[idx]
    strata <- 100L * cut(d, CFG$dist_breaks_km, labels = FALSE, include.lowest = TRUE) +
                     cut(e, CFG$elev_breaks_km, labels = FALSE, include.lowest = TRUE)
    members <- split(seq_along(strata), strata)
    alloc <- allocate_even(lengths(members), CFG$n_pairs)
    pick <- unlist(Map(function(ix, k) if (k >= length(ix)) ix else ix[sample.int(length(ix), k)], members, alloc), use.names = FALSE)
    data.frame(i = i1[idx[pick, 1]], j = i2[idx[pick, 2]], nat_delev = median(e), nat_geo = median(d))
  }
  # n_pairs for EACH pair type in EACH latitude zone: both samples of a pair always come from the same zone
  pl <- list()
  for (z in zone_names) {
    iL <- which(S$Alpine == "Lowland" & S$zone == z); iA <- which(S$Alpine == "Alpine" & S$zone == z)
    one <- function(d, type) if (is.null(d) || !nrow(d)) NULL else data.frame(type = type, zone = z, i = d$i, j = d$j,
                                                                                   nat_delev = d$nat_delev, nat_geo = d$nat_geo)
    pl[[length(pl) + 1]] <- one(if (length(iL) > 1) draw_pairs(iL, iL, TRUE), beta_types[1])
    pl[[length(pl) + 1]] <- one(if (length(iA) > 1) draw_pairs(iA, iA, TRUE), beta_types[2])
    pl[[length(pl) + 1]] <- one(if (length(iA) > 0 && length(iL) > 0) draw_pairs(iA, iL, FALSE), beta_types[3])
  }
  pairs <- bind_rows(pl)
  if (nrow(pairs) == 0) return(NULL)
  pairs$cell <- paste(pairs$type, pairs$zone, sep = " | ")

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
  pairs$d_lat  <- abs(S$lat[pairs$i] - S$lat[pairs$j])            # both samples are in the same zone: at most the zone's width
  pairs$geo_km <- geosphere::distHaversine(cbind(S$lon[pairs$i], S$lat[pairs$i]),
                                           cbind(S$lon[pairs$j], S$lat[pairs$j])) / 1000
  pairs$d_depth <- abs(S$log10_reads[pairs$i] - S$log10_reads[pairs$j])   # difference in sequencing depth
  pairs$d_unm   <- abs(S$unmapped[pairs$i] - S$unmapped[pairs$j])          # difference in unmapped proportion
  pairs[, c("rep_id", "cell", "type", "zone", "bray", "jaccard", "d_elev", "d_lat", "geo_km", "d_depth", "d_unm", "nat_delev", "nat_geo")]
}

## ---- B3. Models ---------------------------------------------------------------------------------------
# A "cell" = pair type x latitude zone (e.g. "Within alpine | Northern high latitudes"). Every cell gets its own curves, so
# zones and pair types can be compared at the same elevation difference, distance, depth and unmapped proportion.
beta_metrics <- c(bray = "Bray\u2013Curtis", jaccard = "Jaccard (rarefied)")
beta_vars    <- c(d_elev = "Elevation difference (km)", d_lat = "Latitude difference within zone (\u00B0)",
                  log10_geo = "Geographic distance (log10 km)")
curve_model  <- c(d_elev = "elev", log10_geo = "elev", d_lat = "lat")      # which model each curve comes from
cell_levels  <- unlist(lapply(beta_types, function(t) paste(t, zone_names, sep = " | ")))

prep_pairs <- function(P) {
  P <- as.data.frame(P)
  tab <- P %>% group_by(cell) %>%
    summarise(n = n(), u = min(n_distinct(round(d_lat, 3)), n_distinct(round(d_elev, 3)), n_distinct(round(geo_km, 3))), .groups = "drop")
  keep <- tab$cell[tab$n >= CFG$min_pairs & tab$u >= 20]
  dropped <- setdiff(unique(P$cell), keep)
  if (length(dropped)) message("    dropped (fewer than ", CFG$min_pairs, " pairs or fewer than 20 distinct values): ", paste(dropped, collapse = "; "))
  P <- P[P$cell %in% keep, ]
  P$cell <- factor(P$cell, levels = intersect(cell_levels, keep))
  P$type <- factor(P$type, levels = beta_types); P$zone <- factor(P$zone, levels = zone_names)
  P$rep_id <- factor(P$rep_id)                                        # replicate: random effect
  P$log10_geo <- log10(pmax(P$geo_km, 0.1))
  P
}

# Models per metric: (elev) elevation difference accounting for geographic distance; (lat) latitude difference accounting for the
# elevation difference. Each at two levels: "overall" = one effect (slope) per pair TYPE, common to the latitude zones, with a
# separate intercept per cell (pair type x zone) so that the slope is estimated WITHIN zones and zone differences, which are
# not the same for every pair type, cannot leak into it = the main figure; "cells" = one effect per cell (pair type x zone) =
# the supplementary figure, the zone comparisons and the standardised dissimilarity. The gradient effects are LINEAR
# (CFG$beta_linear; on the logit scale, in elevation km, latitude degrees and log10 km of distance), the nuisance differences in
# depth and unmapped proportion are smooth.
gterm <- function(v, by) if (isTRUE(CFG$beta_linear)) paste0(by, ":", v) else sm1(v, by)
fit_beta <- function(P, metric, kind, level) {
  n <- nrow(P); P$y <- (P[[metric]] * (n - 1) + 0.5) / n          # squeeze into (0, 1)
  by <- if (level == "cells") "cell" else "type"
  gradient <- if (kind == "elev") c(gterm("d_elev", by), gterm("log10_geo", by)) else c(gterm("d_lat", by), gterm("d_elev", by))
  terms <- c("cell", sm1("d_depth"), sm1("d_unm"), gradient)             # an intercept per cell in both models
  if (nlevels(P$rep_id) > 1) terms <- c(terms, "s(rep_id, bs = 're')")                  # random intercept: replicate
  bam(as.formula(paste("y ~", paste(terms, collapse = " + "))), data = P,
      family = betar(link = "logit"), method = "fREML", discrete = TRUE,
      select = CFG$select, nthreads = CFG$nthreads)
}

# Curve of one pair variable per CELL (the other variables at the cell's typical values, before stratification); the replicate
# effect is left out (a typical replicate); the ribbon is the model's 95% confidence interval
predict_curve <- function(m, P, v) {
  excl <- intersect("s(rep_id)", vapply(m$smooth, function(s) s$label, ""))
  bind_rows(lapply(levels(P$cell), function(cl) {
    sub <- P[P$cell == cl, ]
    x <- seq(quantile(sub[[v]], .01), quantile(sub[[v]], .99), length.out = 80)
    nd <- data.frame(cell = factor(cl, levels = levels(P$cell)), d_elev = median(sub$nat_delev), d_lat = median(sub$d_lat),
                     log10_geo = log10(pmax(median(sub$nat_geo), 0.1)), d_depth = median(sub$d_depth), d_unm = median(sub$d_unm),
                     rep_id = factor(levels(P$rep_id)[1], levels = levels(P$rep_id)))
    nd <- nd[rep(1, length(x)), ]; nd[[v]] <- x
    p <- predict(m, nd, type = "link", se.fit = TRUE, exclude = excl, newdata.guaranteed = TRUE, discrete = FALSE)
    li <- family(m)$linkinv
    data.frame(cell = cl, type = as.character(sub$type[1]), zone = as.character(sub$zone[1]), x = x,
               fit = li(p$fit), lo = li(p$fit - 1.96 * p$se.fit), hi = li(p$fit + 1.96 * p$se.fit))
  }))
}

# Curve of one pair variable per pair TYPE from the "overall" model: the cell intercepts are averaged over the type's zones
# (weights = share of the type's pairs in each zone), so the curve is the effect for the typical mix of zones
predict_overall <- function(m, P, v) {
  excl <- intersect("s(rep_id)", vapply(m$smooth, function(s) s$label, ""))
  li <- family(m)$linkinv; V <- vcov(m); b <- coef(m); zs <- levels(P$zone)
  bind_rows(lapply(levels(P$type), function(ty) {
    sub <- P[P$type == ty, ]
    x <- seq(quantile(sub[[v]], .01), quantile(sub[[v]], .99), length.out = 80)
    w <- as.numeric(table(factor(sub$zone, levels = zs)))
    w[!paste(ty, zs, sep = " | ") %in% levels(P$cell)] <- 0; w <- w / sum(w)
    Xw <- 0
    for (k in which(w > 0)) {
      cl <- paste(ty, zs[k], sep = " | ")
      nd <- data.frame(type = factor(ty, levels = levels(P$type)), cell = factor(cl, levels = levels(P$cell)), zone = factor(zs[k], levels = zs),
                       d_elev = median(sub$nat_delev), d_lat = median(sub$d_lat), log10_geo = log10(pmax(median(sub$nat_geo), 0.1)),
                       d_depth = median(sub$d_depth), d_unm = median(sub$d_unm), rep_id = factor(levels(P$rep_id)[1], levels = levels(P$rep_id)))
      nd <- nd[rep(1, length(x)), ]; nd[[v]] <- x
      Xw <- Xw + w[k] * predict(m, nd, type = "lpmatrix", exclude = excl, newdata.guaranteed = TRUE, discrete = FALSE)
    }
    eta <- as.numeric(Xw %*% b); se <- sqrt(rowSums((Xw %*% V) * Xw))
    data.frame(type = ty, x = x, fit = li(eta), lo = li(eta - 1.96 * se), hi = li(eta + 1.96 * se))
  }))
}

# The linear effects themselves (slopes on the logit scale per unit of the variable, with 95% intervals; p-values descriptive)
slope_table <- function(m, level) {
  b <- coef(m); se <- sqrt(diag(vcov(m)))
  idx <- grep(":(d_elev|d_lat|log10_geo)$", names(b))
  if (!length(idx)) return(NULL)
  nm <- names(b)[idx]
  data.frame(level = level, group = sub(":(d_elev|d_lat|log10_geo)$", "", sub("^(type|cell)", "", nm)), variable = sub("^.*:", "", nm),
             slope = unname(b[idx]), se = unname(se[idx]), lo = unname(b[idx] - 1.96 * se[idx]), hi = unname(b[idx] + 1.96 * se[idx]),
             p = unname(2 * pnorm(-abs(b[idx] / se[idx]))), row.names = NULL)
}
# difference in slope between each zone and the reference zone, within a pair type (cells model)
slope_contrasts <- function(m, ref_zone) {
  b <- coef(m); V <- vcov(m)
  idx <- grep("^cell.*:(d_elev|d_lat|log10_geo)$", names(b)); if (!length(idx)) return(NULL)
  lab <- sub(":(d_elev|d_lat|log10_geo)$", "", sub("^cell", "", names(b)[idx])); var <- sub("^.*:", "", names(b)[idx])
  ty <- sub(" \\| .*$", "", lab); zn <- sub("^.* \\| ", "", lab)
  out <- list()
  for (i in seq_along(idx)) if (zn[i] != ref_zone) {
    j <- which(ty == ty[i] & zn == ref_zone & var == var[i]); if (length(j) != 1) next
    est <- b[idx[i]] - b[idx[j]]; se <- sqrt(V[idx[i], idx[i]] + V[idx[j], idx[j]] - 2 * V[idx[i], idx[j]])
    out[[length(out) + 1]] <- data.frame(pair_type = ty[i], variable = var[i], contrast = sprintf("%s minus %s", zn[i], ref_zone),
                                         slope_difference = unname(est), lo = unname(est - 1.96 * se), hi = unname(est + 1.96 * se),
                                         p = unname(2 * pnorm(-abs(est / se))), row.names = NULL)
  }
  if (length(out)) bind_rows(out) else NULL
}

# STANDARDISED DISSIMILARITY: the dissimilarity two samples of a cell are predicted to have at the SAME geographic distance
# (CFG$ref_km), the same elevation difference (the typical one of the pair type, equal in all zones), and the same differences
# in sequencing depth and unmapped proportion (the medians of all pairs). Distance, elevation and depth are what differ between
# regions otherwise, so this is what lets zones and pair types be compared ("are alpine samples at high latitude more similar
# to each other than alpine samples in the temperate zone?"). Only where the reference distance lies within the cell's
# observed range (5-95th percentile). 95% intervals from 2000 draws of the model coefficients; contrasts are differences of
# those draws (each zone minus the reference zone within a pair type; alpine minus lowland within a zone).
std_dissim <- function(m, P) {
  excl <- intersect("s(rep_id)", vapply(m$smooth, function(s) s$label, ""))
  ref_delev <- tapply(P$nat_delev, P$type, median)      # the typical (natural) elevation difference of each pair type
  rows <- list()
  for (cl in levels(P$cell)) {
    sub <- P[P$cell == cl, ]; rg <- quantile(sub$geo_km, c(0.05, 0.95))
    for (rk in CFG$ref_km) if (rk >= rg[1] && rk <= rg[2])
      rows[[length(rows) + 1]] <- data.frame(cell = cl, type = as.character(sub$type[1]), zone = as.character(sub$zone[1]),
                                             ref_km = rk, n_pairs = nrow(sub), d_elev = unname(ref_delev[as.character(sub$type[1])]))
  }
  if (!length(rows)) return(NULL)
  info <- bind_rows(rows)
  ref_zone <- if (CFG$ref_zone %in% info$zone) CFG$ref_zone else names(which.max(table(P$zone)))     # the zone the others are compared with
  nd <- data.frame(cell = factor(info$cell, levels = levels(P$cell)), d_elev = info$d_elev, log10_geo = log10(info$ref_km),
                   d_depth = median(P$d_depth), d_unm = median(P$d_unm), d_lat = 0,
                   rep_id = factor(levels(P$rep_id)[1], levels = levels(P$rep_id)))
  Xp <- predict(m, nd, type = "lpmatrix", exclude = excl, newdata.guaranteed = TRUE, discrete = FALSE)
  li <- family(m)$linkinv
  set.seed(1)
  draws <- li(Xp %*% t(mgcv::rmvn(2000, coef(m), vcov(m))))
  info$estimate <- as.numeric(li(Xp %*% coef(m)))
  info$lo <- apply(draws, 1, quantile, 0.025); info$hi <- apply(draws, 1, quantile, 0.975)

  cons <- list()
  add_con <- function(i1, i2, comparison, label, rk) {
    dd <- draws[i2, ] - draws[i1, ]
    cons[[length(cons) + 1]] <<- data.frame(comparison = comparison, contrast = label, ref_km = rk,
                                            estimate = info$estimate[i2] - info$estimate[i1],
                                            lo = unname(quantile(dd, 0.025)), hi = unname(quantile(dd, 0.975)),
                                            share_of_draws_above_0 = mean(dd > 0))
  }
  for (rk in unique(info$ref_km)) {
    ir <- which(info$ref_km == rk)
    for (ty in beta_types) {                                       # each zone vs the reference zone, within a pair type
      iz <- ir[info$type[ir] == ty]; i0 <- iz[info$zone[iz] == ref_zone]
      if (length(i0) == 1) for (i1 in iz[info$zone[iz] != ref_zone])
        add_con(i0, i1, "zone", sprintf("%s: %s minus %s", ty, info$zone[i1], ref_zone), rk)
    }
    for (zz in zone_names) {                                       # alpine vs lowland, within a zone
      i1 <- ir[info$type[ir] == beta_types[1] & info$zone[ir] == zz]; i2 <- ir[info$type[ir] == beta_types[2] & info$zone[ir] == zz]
      if (length(i1) == 1 && length(i2) == 1) add_con(i1, i2, "type", sprintf("%s: %s minus %s", zz, beta_types[2], beta_types[1]), rk)
    }
  }
  list(std = info, con = if (length(cons)) bind_rows(cons) else NULL)
}

# ---- all replicates: the random pairs and their dissimilarities. Saved (beta/pairs_all.rds) and reused on a rerun if the
# settings and the drawn samples are the same, so a problem in the models does not cost the distance computation again.
pairs_key <- list(dist_breaks = CFG$dist_breaks_km, elev_breaks = CFG$elev_breaks_km, n_rep = CFG$n_rep, n_per_group = CFG$n_per_group, n_pairs = CFG$n_pairs, rare_depth = CFG$rare_depth,
                  n_ids = length(all_ids), ids_chars = sum(nchar(all_ids)), lat_breaks = CFG$lat_breaks,
                  zones = zone_names, thin_m = CFG$thin_m, test = TEST_MODE, design = "cells-v3")
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
message(sprintf("%d pairs from %d replicates in %d cells (pair type x latitude zone)", nrow(P), nlevels(P$rep_id), nlevels(P$cell)))
print(table(cell = P$cell))
pair_counts <- P %>% count(rep_id, cell, type, zone)
# share of pairs with dissimilarity exactly 1 (no shared OTU): a pile at the boundary of a beta regression
at_one <- bind_rows(lapply(names(beta_metrics), function(mt)
  P %>% group_by(cell) %>% summarise(share_at_1 = mean(.data[[mt]] >= 1 - 1e-12), n_pairs = n(), .groups = "drop") %>% mutate(metric = mt)))
write.csv(at_one, out_b("share_of_pairs_at_dissimilarity_1.csv"), row.names = FALSE)
for (i in which(at_one$metric == "bray")) message(sprintf("  %s: %.1f%% of the Bray-Curtis pairs are exactly 1", at_one$cell[i], 100 * at_one$share_at_1[i]))
write.csv(pair_counts, out_b("pair_counts.csv"), row.names = FALSE)
beta_raw <- P[, c("cell", "type", "zone", "d_elev", "d_lat", "log10_geo", names(beta_metrics))]   # every pair (for the plots)

one <- function(x) if (length(x) == 1 && is.numeric(x)) as.numeric(x) else NA_real_     # always exactly one number
ref_zone_g <- if (CFG$ref_zone %in% levels(P$zone)) CFG$ref_zone else names(which.max(table(P$zone)))
preds_o <- list(); preds_c <- list(); summ <- list(); std_list <- list(); con_list <- list(); slope_list <- list(); slopecon_list <- list()
for (metric in names(beta_metrics)) for (kind in c("elev", "lat")) for (level in c("overall", "cells")) {
  message("  fitting ", metric, " / ", kind, " / ", level, " model (pooled, random effect: replicate)")
  m <- tryCatch(fit_beta(P, metric, kind, level), error = function(e) { message("  FAILED: ", conditionMessage(e)); NULL })
  if (is.null(m)) next
  key <- paste(metric, kind, level)
  vars_k <- names(curve_model)[curve_model == kind]
  if (level == "overall") {
    preds_o[[key]] <- bind_rows(lapply(vars_k, function(v) cbind(metric = metric, variable = v, predict_overall(m, P, v))))
  } else {
    preds_c[[key]] <- bind_rows(lapply(vars_k, function(v) cbind(metric = metric, variable = v, predict_curve(m, P, v))))
  }
  sm_m <- summary(m)
  st <- as.data.frame(sm_m$s.table); if (nrow(st)) st$term <- rownames(st) else st <- data.frame(term = character(0))
  rep_sd <- tryCatch({                  # standard deviation of the replicate random effect (NA if returned in another shape)
    vc <- suppressWarnings(gam.vcomp(m, rescale = FALSE))
    if (is.list(vc) && !is.matrix(vc)) vc <- vc$vc
    one(if (is.matrix(vc)) vc["s(rep_id)", 1] else vc["s(rep_id)"])
  }, error = function(e) NA_real_)
  summ[[key]] <- data.frame(metric = metric, model = kind, level = level, dev_expl = one(sm_m$dev.expl), n_pairs = nrow(P),
                            replicate_sd = rep_sd, st, check.names = FALSE, stringsAsFactors = FALSE)
  if (isTRUE(CFG$beta_linear)) {
    sl <- slope_table(m, level); if (!is.null(sl)) slope_list[[key]] <- cbind(metric = metric, model = kind, sl)
    if (level == "cells") { sc <- slope_contrasts(m, ref_zone_g); if (!is.null(sc)) slopecon_list[[key]] <- cbind(metric = metric, model = kind, sc) }
  }
  if (kind == "elev" && level == "cells") {
    sdz <- tryCatch(std_dissim(m, P), error = function(e) { message("  standardised dissimilarity FAILED: ", conditionMessage(e)); NULL })
    if (!is.null(sdz)) { std_list[[metric]] <- cbind(metric = metric, sdz$std)
                         if (!is.null(sdz$con)) con_list[[metric]] <- cbind(metric = metric, sdz$con) }
  }
  saveRDS(m, out_b(paste0("model_", metric, "_", kind, "_", level, ".rds")))
}
if (!length(summ)) stop("no beta model could be fitted: see the FAILED messages above")
preds_o <- bind_rows(preds_o); preds_c <- bind_rows(preds_c)
write.csv(bind_rows(summ), out_b("models_summary.csv"), row.names = FALSE)      # p-values: descriptive only
write.csv(preds_o, out_b("predicted_curves_overall.csv"), row.names = FALSE)
write.csv(preds_c, out_b("predicted_curves_by_zone.csv"), row.names = FALSE)
if (length(slope_list)) {
  write.csv(bind_rows(slope_list), out_b("linear_effects.csv"), row.names = FALSE)              # slopes (logit scale) per pair type / cell
  write.csv(bind_rows(slopecon_list), out_b("linear_effect_zone_contrasts.csv"), row.names = FALSE)
  message("Linear effects (logit scale per unit; overall = per pair type, common to the zones) are in beta/linear_effects.csv; zone differences in the slopes in beta/linear_effect_zone_contrasts.csv")
}
std_all <- bind_rows(std_list); con_all <- bind_rows(con_list)
write.csv(std_all, out_b("standardized_dissimilarity.csv"), row.names = FALSE)
write.csv(con_all, out_b("standardized_contrasts.csv"), row.names = FALSE)
if (nrow(con_all)) {
  message("Standardised Bray-Curtis dissimilarity, contrasts (95% interval; > 0 = the first-named group is MORE dissimilar):")
  cb <- con_all[con_all$metric == "bray", ]
  for (i in seq_len(nrow(cb))) message(sprintf("  at %g km | %s: %+.3f [%+.3f, %+.3f]", cb$ref_km[i], cb$contrast[i], cb$estimate[i], cb$lo[i], cb$hi[i]))
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
thin_txt <- if (CFG$thin_m > 0) sprintf("Samples are spatially thinned: at most one per %g m", CFG$thin_m) else "Samples with identical coordinates are removed (one is kept)"
pts_note <- if (CFG$beta_points_max > 0) sprintf("points = a random subset of %s pairs (the y axis is zoomed on the bulk of the pairs: the few with a much lower dissimilarity are not drawn and ribbons are clipped at the lower edge, but all pairs are in the models)", format(CFG$beta_points_max, big.mark = ",")) else "pairs not drawn"
eff_txt  <- if (isTRUE(CFG$beta_linear)) "linear effects (logit scale; distance as log10 km)" else "smooth effects"
lab_beta <- function(d) {
  d$metric   <- factor(beta_metrics[d$metric], levels = beta_metrics)
  d$variable <- factor(beta_vars[d$variable], levels = beta_vars)
  d
}
type_scales <- list(scale_colour_manual(values = CFG$cols_beta, name = NULL), scale_fill_manual(values = CFG$cols_beta, guide = "none"),
                    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1, linewidth = 1.2, linetype = "solid"))))
fac_type_zone <- function(d) {
  d$type <- factor(d$type, levels = beta_types)
  if ("zone" %in% names(d)) d$zone <- factor(d$zone, levels = zone_names)
  d
}
# random subset of the pairs for drawing (all pairs are in the models)
pair_points <- function(vars, metrics, max_n) {
  ob <- bind_rows(lapply(metrics, function(mt) bind_rows(lapply(vars, function(v)
    data.frame(metric = mt, variable = v, type = beta_raw$type, zone = beta_raw$zone, x = beta_raw[[v]], y = beta_raw[[mt]])))))
  set.seed(5)
  if (max_n <= 0) ob[0, ] else if (nrow(ob) > max_n) ob[sample.int(nrow(ob), max_n), ] else ob
}

# One panel of the OVERALL effects: one metric against one pair variable, a curve per pair type (common to the zones)
beta_overall_panel <- function(metric, v, tag, subtitle) {
  pm <- preds_o[preds_o$metric == metric & preds_o$variable == v, ]
  if (!nrow(pm)) return(NULL)
  ob <- pair_points(v, metric, CFG$beta_points_max); pm <- fac_type_zone(pm); ob <- fac_type_zone(ob)
  pm$mid <- pm$fit; ob$.p <- 1; pm$.p <- 1
  if (nrow(ob)) { z <- zoom_panel(ob, pm, ".p"); ob <- z$ob; pm <- z$pm }
  ggplot() +
    geom_ribbon(data = pm, aes(x, ymin = lo, ymax = hi, fill = type), alpha = 0.18) +
    geom_point(data = ob, aes(x, y, colour = type), size = 0.35, alpha = CFG$beta_point_alpha, shape = 16) +
    geom_line(data = pm, aes(x, mid, colour = type), linewidth = 1.1) +
    type_scales + y_scale +
    labs(tag = tag, x = beta_vars[[v]], y = paste(beta_metrics[[metric]], "dissimilarity"), subtitle = subtitle) +
    theme_fig + theme(plot.subtitle = element_text(size = 9, colour = "grey30"))
}

# One metric by latitude ZONE: rows = zones, columns = pair variables, a curve per pair type
beta_zone_grid <- function(metric) {
  pm <- preds_c[preds_c$metric == metric, ]
  if (!nrow(pm)) return(NULL)
  ob <- pair_points(names(beta_vars), metric, 3 * CFG$beta_points_max)
  pm <- fac_type_zone(lab_beta(pm)); ob <- fac_type_zone(lab_beta(ob)); pm$mid <- pm$fit
  if (nrow(ob)) { z <- zoom_panel(ob, pm, c("variable", "zone")); ob <- z$ob; pm <- z$pm }
  ggplot() +
    geom_ribbon(data = pm, aes(x, ymin = lo, ymax = hi, fill = type, group = cell), alpha = 0.18) +
    geom_point(data = ob, aes(x, y, colour = type), size = 0.3, alpha = CFG$beta_point_alpha, shape = 16) +
    geom_line(data = pm, aes(x, mid, colour = type, group = cell), linewidth = 1) +
    facet_grid(zone ~ variable, scales = "free", switch = "x") +            # zones in rows, the x variables as column titles under the axes
    type_scales + y_scale +
    labs(x = NULL, y = paste(beta_metrics[[metric]], "dissimilarity")) +
    theme_bw(base_size = 12) +
    theme(strip.placement = "outside", strip.background = element_blank(), strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank(), legend.position = "bottom")
}

# Overall full grid (all metrics x all pair variables): beta/
{
  pm <- fac_type_zone(lab_beta(preds_o)); ob <- fac_type_zone(lab_beta(pair_points(names(beta_vars), names(beta_metrics), 3 * CFG$beta_points_max)))
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
  ggsave(out_b("beta_overall_grid.png"), g, width = 13, height = 7.5, dpi = 300)
}
gz <- beta_zone_grid("jaccard"); if (!is.null(gz)) ggsave(out_b("beta_by_zone_jaccard.png"), gz, width = 13, height = 9, dpi = 300)

# Supplementary Figure 4: Bray-Curtis by latitude zone
gz <- beta_zone_grid("bray")
if (!is.null(gz)) {
  cap4 <- paste0("Supplementary Figure 4. Bray-Curtis dissimilarity between pairs of samples by latitude zone (rows; ", paste(zone_desc, collapse = "; "),
                 ") against the elevation difference (accounting for geographic distance), the latitude difference within the zone (accounting for the elevation difference) and the geographic distance (accounting for the elevation difference); colour = pair type. ",
                 "One beta-regression GAM per column group (elevation and distance share a model) on the pooled pairs of ", CFG$n_rep, " replicates with a random intercept for the replicate and ", eff_txt,
                 " per pair type and zone; lines = fitted effects, ribbons = 95% confidence intervals (descriptive and optimistic: pairs reuse samples), ", pts_note, ". ", thin_txt,
                 ". Slopes: beta/linear_effects.csv; zone differences in the slopes: beta/linear_effect_zone_contrasts.csv.")
  gz <- gz + labs(caption = paste(strwrap(cap4, 150), collapse = "\n")) + theme(plot.caption = element_text(hjust = 0, size = 8.5))
  ggsave(file.path(out_folder, "figS4_beta_by_zone.png"), gz, width = 13, height = 10.5, dpi = 300)
}

# Standardised dissimilarity (Figure 2b): the dissimilarity predicted at the SAME geographic distance (columns), the typical
# elevation difference of the pair type and median depth / unmapped differences, per latitude zone and pair type
if (nrow(std_all)) {
  sdp <- std_all
  sdp$type <- factor(sdp$type, levels = beta_types); sdp$zone <- factor(sdp$zone, levels = zone_names)
  sdp$metric <- factor(beta_metrics[sdp$metric], levels = beta_metrics)
  ref_levels <- sprintf("at %s km apart", formatC(sort(unique(sdp$ref_km)), format = "d", big.mark = ","))
  sdp$ref <- factor(sprintf("at %s km apart", formatC(sdp$ref_km, format = "d", big.mark = ",")), levels = ref_levels)
  cap2 <- paste("Figure 2b. Dissimilarity between two samples of the same latitude zone, standardised: predicted at the same geographic distance (columns),",
                "the typical elevation difference of the pair type (the same in all zones) and the median differences in sequencing depth and unmapped",
                "proportion, from the by-zone model of Supplementary Figure 4 (elevation and distance). Points = model estimate, bars = 95% interval (descriptive and optimistic: pairs reuse",
                paste0("samples); shown only where the distance lies within the 5-95th percentile of the cell's pairs. Zones: ", paste(zone_desc, collapse = "; "),
                       ". ", thin_txt, ". Differences between zones (each minus the reference zone) and pair types with intervals: beta/standardized_contrasts.csv."))
  g2 <- ggplot(sdp, aes(zone, estimate, colour = type)) +
    geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.6), size = 0.5, linewidth = 0.8) +
    facet_grid(metric ~ ref, scales = "free_y", switch = "y") +
    scale_colour_manual(values = CFG$cols_beta, name = NULL) +
    scale_x_discrete(labels = scales::label_wrap(14), drop = FALSE) +
    labs(x = NULL, y = "Standardised dissimilarity", caption = paste(strwrap(cap2, 130), collapse = "\n")) +
    theme_bw(base_size = 12) +
    theme(strip.placement = "outside", strip.background = element_blank(), strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank(), legend.position = "bottom", plot.caption = element_text(hjust = 0, size = 8.5))
  ggsave(file.path(out_folder, "fig2b_standardized_dissimilarity.png"), g2, width = 11, height = 7.5, dpi = 300)
}

# Supplementary Figure 1: Jaccard (rarefied presence/absence), overall effects
pS1A <- beta_overall_panel("jaccard", "d_elev", "A", "accounting for geographic distance")
pS1B <- beta_overall_panel("jaccard", "d_lat", "B", "accounting for the elevation difference")
if (!is.null(pS1A) && !is.null(pS1B)) {
  s1_caption <- paste("Supplementary Figure 1. Jaccard (rarefied presence/absence) dissimilarity between pairs of samples against the difference in elevation (A; accounting for geographic distance)",
                      "and in latitude (B; accounting for the elevation difference); colour = pair type; pairs are always from the same latitude zone and the effects are common to the zones (an intercept per pair type and zone); one beta-regression GAM per panel on the pooled pairs of",
                      CFG$n_rep, paste0("replicates with a random intercept for the replicate and ", eff_txt, "; lines = fitted effects, ribbons = 95% confidence intervals (descriptive and optimistic), "),
                      paste0(pts_note, "."), paste0(thin_txt, ". By zone: beta/beta_by_zone_jaccard.png."))
  assemble_fig(file.path(out_folder, "figS1_beta_jaccard.png"), list(list(pS1A, pS1B)), get_legend(pS1A), s1_caption, 14, 6.5)
}

}  # end RUN_BETA

## ======================= MAIN FIGURE (needs both parts) =====================================
# A: alpha diversity along elevation (row 1); B: along latitude (row 2); C: Bray-Curtis dissimilarity against the elevation
# difference (accounting for geographic distance) and the latitude difference (accounting for elevation): the OVERALL effects
# per pair type (common to the latitude zones); by zone: Supplementary Figure 4
if (RUN_ALPHA && RUN_BETA) {
  pA <- alpha_panel(grad_lab[["elevation"]], "A")
  pB <- alpha_panel(grad_lab[["glat"]], "B")
  pC <- beta_overall_panel("bray", "d_elev", "C", "accounting for geographic distance")
  pD <- beta_overall_panel("bray", "d_lat", "", "accounting for the elevation difference")
  if (is.null(pC) || is.null(pD)) stop("main figure: the Bray-Curtis models are missing")
  main_caption <- paste(
    "A, B: alpha diversity (observed OTUs, Shannon index, Pielou J) along elevation (A) and latitude (B; dotted vertical lines = the zone boundaries used for the beta analysis); GAM partial effects",
    "adjusted for sequencing depth and unmapped proportion, solid = net of space (95% CI), dashed = no spatial term; stars = test that the alpine gradient differs in shape from the",
    "lowland one (net model; *** p < 0.001, ** p < 0.01, * p < 0.05, ns).",
    "C: Bray-Curtis dissimilarity between pairs of samples (colour = pair type: within lowland, within alpine, between alpine and lowland) against the difference in elevation",
    "(left; accounting for geographic distance) and the difference in latitude (right; accounting for the elevation difference): the overall effects, one per pair type and common to the latitude zones",
    paste0("(pairs always come from the same zone and every pair type and zone has its own intercept; by zone: Supplementary Figure 4; zones: ", paste(zone_desc, collapse = "; "), ")."),
    sprintf("Each replicate subsample draws up to %d samples per group and zone and up to %d pairs per pair type and zone (random within classes of geographic distance and elevation difference); the pairs of %d replicates are pooled in one beta-regression GAM per panel with a random intercept for the replicate and %s; lines = fitted effects, ribbons = 95%% confidence intervals (descriptive and optimistic: pairs reuse samples),",
            CFG$n_per_group, CFG$n_pairs, CFG$n_rep, eff_txt), paste0(pts_note, "."),
    paste0(thin_txt, " in every replicate subsample of C (each replicate draws a different sample per location). ", alpha_txt),
    "Dissimilarity standardised for distance and compared between zones: Figure 2b. Jaccard: Supplementary Figure 1.")
  assemble_fig(file.path(out_folder, "fig2_alpha_beta_gradients.png"),
               list(list(pA), list(pB), list(get_legend(pA)), list(pC, pD)), get_legend(pC), main_caption,
               width = 14, height = 15, row_heights = c(1, 1, 0.17, 1.3))
  saveRDS(list(A = pA, B = pB, C = pC, D = pD), file.path(out_folder, "fig2_panels.rds"))
} else message("Main figure skipped: it needs RUN_ALPHA and RUN_BETA both TRUE.")

message("Done. Outputs in ", out_folder)
writeLines(capture.output(sessionInfo()), file.path(out_folder, "sessionInfo.txt"))
file.create(file.path(out_folder, ".done"))   # flag file for Snakemake