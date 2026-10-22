# =============================================================================
# Script 3 - Distance decay of community dissimilarity: are mountains "islands of diversity"?
#
# Question: do alpine soils differ MORE from each other with distance than lowland soils do (alpine samples more dissimilar
#   at the same geographic distance, and/or dissimilarity rising faster with distance)? Peri-alpine soils (lowland samples
#   within PERIALPINE_KM of an alpine polygon, column perialpine_alpine from microbeatlas_filter_data.R) are the control for
#   "being next to mountains": if peri-alpine soils behave like alpine ones, it is the region, not the alpine habitat.
#
# DESIGN (n_rep replicates; the same logic as script 2)
#   * Three groups: alpine, peri-alpine, lowland (= neither). Pairs of samples are WITHIN a group: the distance decay of each group.
#   * Every replicate draws, in each group, ONE random sample per cluster (cluster_id: samples of one BioProject or one location
#     are linked), at most one per thin_m metres (no pairs at 0 km from different studies), and cuts every group down to the SAME
#     number of samples (the smallest group decides, at most n_per_group), EVENLY over latitude / longitude cells (even_cell_deg):
#     one sample per cell in turn. So the groups have the same size, no study or place dominates, and no region (Europe, North
#     America) dominates a group.
#   * Up to n_pairs pairs per group, drawn at random WITHIN classes of geographic distance (dist_breaks_km), shared as evenly as
#     possible over the classes: random pairs are mostly thousands of km apart, so short distances would hardly be represented.
#     Stratifying on the predictor does not bias the fitted curve; it makes it well determined over the whole range.
#   * Bray-Curtis (relative abundances) and Jaccard (presence/absence after rarefaction to rare_depth reads), computed only for the
#     drawn pairs from sparse matrices (checked against vegan::vegdist in TEST_MODE).
#   * One LINEAR model with INTERACTIONS per replicate and metric on the dissimilarity itself (0-1; a beta regression does not
#     converge, many dissimilarities are exactly 1): dissimilarity ~ group * (log10 distance - log10 of the reference distance) + linear
#     nuisance terms for the differences in sequencing depth and unmapped proportion (df_dist > 1: a natural spline instead of the line).
#     Distance is centred at the reference distance (1,000 km), so the coefficients read directly:
#       (Intercept)               dissimilarity of the reference group (lowland) at 1,000 km
#       group<G>                  <G> minus lowland: the difference in LEVEL at 1,000 km
#       distance                  rise of the lowland dissimilarity per 10-fold distance
#       group<G>:distance         <G> minus lowland: the difference in the RISE per 10-fold distance (the interaction)
#     (tables/model_coefficients.csv: pooled over the replicates).
#
# STATISTICAL TESTS (tables/tests.csv; stars in the figures)
#   First two JOINT tests per metric (Wald tests on the pooled coefficients): do the groups differ in level at the reference distance (the group
#   terms), and do they differ in the rise with distance (the group x distance interactions)? Then the contrasts between groups (alpine - lowland, peri-alpine - lowland, alpine - peri-alpine) of
#     - the dissimilarity at the SAME distance (ref_km)                      -> "more dissimilar across distance"
#     - the TURNOVER: the rise of the dissimilarity per 10-fold increase of distance, over the distance range that all groups have in
#       common (common_q)                                                    -> "dissimilarity rises faster with distance"
#   Pairs that share a sample are not independent, so ordinary standard errors would be far too small (in a simulation: 67% false
#   positives at p < 0.05). The standard errors are therefore cluster-robust for DYADIC data (Aronow, Samii & Assenova 2015: the scores
#   of all pairs that share a sample are summed). The 20 replicates are pooled: the estimate is their mean, its standard error the root
#   mean square of the within-replicate standard errors (the averaging is not credited: the replicates redraw from the same samples). In
#   a simulation this gave about 2% false positives under a true null and 72% power for a true difference of half a standard error unit
#   (see the header of the stats block). p-values are adjusted for the number of contrasts (Benjamini-Hochberg) over each analysis.
#   Also reported: the range of the estimates between replicates and the share of replicates in which the test alone has p < 0.05.
#   LIMITS: the tests cannot see spatial autocorrelation between DIFFERENT samples (neighbouring regions resemble each other), and the
#   clusters are the sampling units: judge robustness by the size of the effects and their consistency between replicates.
#
# MAP (panels C, D): how distinct is the alpine community of a region from the alpine communities of the OTHER regions? Regions =
#   region_deg-degree latitude / longitude cells with at least min_region_n alpine samples. In every replicate (one random alpine sample per
#   cluster, ALL clusters) the mean Bray-Curtis dissimilarity to the alpine samples of other regions is computed (C), and the same after
#   removing the distance decay (D): the residual of the dissimilarity from one alpine distance decay over all between-region pairs. D > 0:
#   more distinct than the distance to the other regions explains (an "island"), D < 0: more similar. Median over the replicates.
#
# SUPPLEMENTARY (figS5_latitude_zones.png): the same analysis within latitude zones (lat_breaks: southern hemisphere, temperate 0-60
#   degrees N, polar > 60 degrees N): every zone has its own matched groups, curves and tests. A group needs min_zone_clusters clusters
#   in a zone to be analysed there; zones without enough samples are left out.
#
# FIGURE (in <out_folder>)  fig3_distance_decay.png : A, B = Bray-Curtis / Jaccard dissimilarity against distance per group (lines =
#   median over the replicates, ribbons = range between replicates, open circles = observed dissimilarity adjusted to the median depth /
#   unmapped differences, mean in 8 bins of distance; stars = tests); C, D = map of the distinctiveness of alpine regions.
# TABLES (in <out_folder>/tables)  tests.csv (all contrasts and group values), curves.csv, map_regions.csv, model_checks.csv,
#   samples_per_replicate.csv, pair_counts.csv, distance_range_per_group.csv, share_of_pairs_at_dissimilarity_1.csv
#
# Input : <data_folder>/sample_data_filtered.tsv (needs perialpine_alpine and cluster_id: rerun prepare_microbeatlas_data if absent)
#         and the filtered BIOM written by microbeatlas_filter_data.R
# Usage : Rscript 3_distance_decay.R [sample_table] [filtered_biom] [out_folder]
#         Snakemake passes the exact paths (all file names are defined in the Snakefile only)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(Matrix)
  library(rbiom)
  library(geosphere)
  library(vegan)
  library(grid)
  library(splines)
})

## ---- 0. Switches, folders, configuration ------------------------------------------
TEST_MODE  <- FALSE   # TRUE: tiny run (few samples, 1 replicate) to check that it works
REUSE_FITS <- TRUE    # reuse the saved pairs (pairs_all.rds) and counts if they match the data and settings

args <- commandArgs(trailingOnly = TRUE)
data_folder     <- "/home/renku/work/epfl-altshuler/GlobalAlpineSoilMicrobiome/data"
results_default <- file.path(dirname(data_folder), "results")
sample_table  <- if (length(args) >= 1) as.character(args[1]) else file.path(data_folder, "sample_data_filtered.tsv")
filtered_biom <- if (length(args) >= 2) as.character(args[2]) else
  file.path(data_folder, "microbeatlas", "otuTable_soilCurated_OTU97_nextMAPrelease_filtered.biom")
out_folder    <- if (length(args) >= 3) as.character(args[3]) else file.path(results_default, "3_distance_decay")

CFG <- list(
  sample_file = sample_table,
  biom_file   = filtered_biom,  # written by the Snakemake rule filter_biom
  biom_id     = "Accessions",   # column of the sample table that matches the BIOM sample IDs
  group_col   = "perialpine_alpine",                       # Lowland / Perialpine / Alpine (microbeatlas_filter_data.R)
  groups      = c("Lowland", "Perialpine", "Alpine"),      # order of the groups in tables and figures
  n_rep       = 20,             # replicates: each draws the samples (one per cluster, even cut) and the pairs again from scratch
  n_per_group = 1000,           # at most this many samples per group and replicate (the same number in every group: the smallest decides)
  even_cell_deg = 5,            # groups are cut down evenly over latitude / longitude cells of this size (degrees)
  thin_m      = 100,            # at most one sample per this many metres within a group (0 = only distinct coordinates)
  rare_depth  = 1000,           # rarefaction depth for Jaccard (samples below are dropped)
  n_pairs     = 2000,           # pairs per group and replicate (shared as evenly as possible over the distance classes)
  dist_breaks_km = c(0, 3, 10, 30, 100, 300, 1000, 3000, 10000, Inf),   # pairs are drawn WITHIN these geographic distance classes (km)
  df_dist     = 1,              # 1 = a straight line in log10 distance per group; 2, 3 = a natural spline (a curve)
  ref_km      = c(100, 1000),   # reference distances (km) at which the groups are compared (only where all groups have pairs)
  common_q    = c(0.05, 0.95),  # the distance range common to all groups = from the highest 5th to the lowest 95th percentile
  lat_breaks  = c(0, 60),       # latitude zones of the supplementary figure (degrees): below 0, 0 to 60, above 60
  zone_names  = c("Southern hemisphere", "Northern 0-60\u00B0 (temperate)", "Northern >60\u00B0 (polar)"),
  min_zone_clusters = 10,       # a group needs at least this many clusters in a zone to be analysed there
  region_deg  = 5,              # alpine regions of the map: latitude / longitude cells of this size (degrees)
  min_region_n = 3,             # a region needs at least this many alpine samples (clusters) to appear on the map
  map_max_samples = 800,        # at most this many alpine samples in the all-pairs dissimilarity of the map (a random subset above that)
  points_max  = 40000,          # pairs drawn as points (random subset; the models use all pairs)
  point_alpha = 0.05,
  cols        = c(Lowland = "#7C9885", Perialpine = "#CC79A7", Alpine = "#2B6CB0")
)
if (TEST_MODE) {
  out_folder <- paste0(out_folder, "_TEST")
  CFG$n_per_group <- 60; CFG$n_rep <- 1; CFG$n_pairs <- 300; CFG$min_zone_clusters <- 8
}
dir.create(file.path(out_folder, "tables"), recursive = TRUE, showWarnings = FALSE)
out_t <- function(...) file.path(out_folder, "tables", ...)
set.seed(1)

# ---- helpers ------------------------------------------------------------------------------------------------------
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

# Order in which to take samples so that any first n are spread as evenly as possible over cell_deg-degree latitude / longitude
# cells: one sample per cell in turn (random sample within a cell, random cell order in a partial round)
even_order <- function(lat, lon, cell_deg = CFG$even_cell_deg) {
  cell <- paste(floor(lat / cell_deg), floor(lon / cell_deg))
  perm <- sample.int(length(lat)); rk <- integer(length(lat))
  rk[perm] <- ave(seq_along(perm), cell[perm], FUN = seq_along)                # rank of a sample within its cell
  order(rk, match(cell, sample(unique(cell))))
}

fmt_km <- function(x) format(signif(x, 2), big.mark = ",", scientific = FALSE, trim = TRUE)                    # distances for labels and messages
qlo <- function(x) unname(quantile(x, 0.025, na.rm = TRUE)); qhi <- function(x) unname(quantile(x, 0.975, na.rm = TRUE))
clip01 <- function(x) pmin(pmax(x, 0), 1)
metrics <- c(bray = "Bray\u2013Curtis", jaccard = "Jaccard (rarefied)")
metrics_ascii <- c(bray = "Bray-Curtis", jaccard = "Jaccard (rarefied)")

theme_fig <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "none",
        plot.tag = element_text(face = "bold", size = 16),
        strip.background = element_rect(fill = "grey95"), strip.text = element_text(face = "bold"))

# Extract the legend of a ggplot (robust to ggplot2 versions)
get_legend <- function(p) {
  g <- ggplotGrob(p + theme(legend.position = "bottom"))
  for (i in grep("guide-box", g$layout$name))
    if (!inherits(g$grobs[[i]], "zeroGrob")) return(g$grobs[[i]])
  NULL
}

# Assemble rows of panels (1 to 3 per row) + a shared legend + a caption into one PNG
assemble_fig <- function(file, rows, leg, caption, width, height, row_heights = rep(1, length(rows))) {
  nr <- length(rows); nc <- max(lengths(rows))
  cap_lines <- strwrap(caption, 165)
  png(file, width = width, height = height, units = "in", res = 300)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nr + 2, nc,
                        heights = unit.c(unit(row_heights, "null"), unit(0.5, "in"),
                                         unit(0.17 * length(cap_lines) + 0.2, "in")))))
  for (i in seq_len(nr)) for (j in seq_along(rows[[i]])) {
    vp <- viewport(layout.pos.row = i, layout.pos.col = if (length(rows[[i]]) == 1) seq_len(nc) else j)
    print(rows[[i]][[j]], vp = vp)
  }
  if (!is.null(leg)) {
    pushViewport(viewport(layout.pos.row = nr + 1, layout.pos.col = seq_len(nc))); grid.draw(leg); popViewport()
  }
  pushViewport(viewport(layout.pos.row = nr + 2, layout.pos.col = seq_len(nc)))
  grid.text(paste(cap_lines, collapse = "\n"), x = 0.01, just = "left", gp = gpar(fontsize = 9))
  popViewport(); dev.off()
}

# BEGIN STATS
# ---- Statistical tests for pair data ---------------------------------------------------------------------------------------------
# Pairs that share a sample are NOT independent (a sample with an unusual community is unusual in all its pairs), so ordinary
# standard errors of a regression on pairs are far too small. The variance used here is the cluster-robust variance for dyadic data
# (Aronow, Samii & Assenova 2015): the scores of all pairs that share a sample are summed per sample.
dyadic_vcov <- function(X, e, node_a, node_b) {
  s <- X * e                                                      # score of every pair (n x k)
  bread <- tryCatch(solve(crossprod(X)), error = function(err) NULL)
  if (is.null(bread)) return(NULL)
  S <- rowsum(rbind(s, s), c(node_a, node_b))                     # sum of the scores of all pairs of each sample
  meat <- crossprod(S) - crossprod(s)                             # (a pair is in the sums of both of its samples: count it once)
  G <- nrow(S)
  bread %*% meat %*% bread * (G / (G - 1))
}
# Design-matrix rows of a linear model at new data (also for natural splines: the knots of the fit are used)
design_at <- function(m, newdata) {
  tt <- delete.response(terms(m))
  model.matrix(tt, model.frame(tt, newdata, xlev = m$xlevels), contrasts.arg = m$contrasts)
}
# Combine the estimates of the replicates: the pooled estimate is the MEAN over the replicates; its standard error is the root mean square
# of the within-replicate (dyadic-robust) standard errors, i.e. the averaging over replicates is NOT credited (the replicates redraw
# from the same samples, so they are not independent). Validated by simulation: about 2% false positives at p < 0.05 under a true
# null (conservative), 5% for a single replicate. Also given: the range of the estimates between replicates and the share of
# replicates in which the test alone gives p < 0.05.
pool_replicates <- function(est, se) {
  ok <- is.finite(est) & is.finite(se); est <- est[ok]; se <- se[ok]; m <- length(est)
  if (m < 1) return(data.frame(estimate = NA_real_, se = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, range_lo = NA_real_, range_hi = NA_real_, share_replicates_p_below_0.05 = NA_real_, n_replicates = 0))
  q <- mean(est); s_w <- sqrt(mean(se^2))
  data.frame(estimate = q, se = s_w, lo = q - qnorm(0.975) * s_w, hi = q + qnorm(0.975) * s_w, p = 2 * pnorm(-abs(q / s_w)),
             range_lo = unname(quantile(est, 0.025)), range_hi = unname(quantile(est, 0.975)),
             share_replicates_p_below_0.05 = mean(2 * pnorm(-abs(est / se)) < 0.05), n_replicates = m)
}
# Joint (omnibus) Wald test that several coefficients are all zero, from the replicates: mean coefficient vector and mean covariance matrix
# (the same pooling as for single contrasts: the averaging is not credited, so it is conservative)
pool_wald <- function(beta_list, V_list, idx) {
  b <- rowMeans(do.call(cbind, lapply(beta_list, function(x) x[idx]))); V <- Reduce(`+`, lapply(V_list, function(x) x[idx, idx, drop = FALSE])) / length(V_list)
  w <- tryCatch(drop(t(b) %*% solve(V, b)), error = function(e) NA_real_)
  p_r <- vapply(seq_along(beta_list), function(r) tryCatch(pchisq(drop(t(beta_list[[r]][idx]) %*% solve(V_list[[r]][idx, idx, drop = FALSE], beta_list[[r]][idx])), length(idx), lower.tail = FALSE),
                                                           error = function(e) NA_real_), 1)
  data.frame(wald_chisq = w, df = length(idx), p = pchisq(w, length(idx), lower.tail = FALSE), share_replicates_p_below_0.05 = mean(p_r < 0.05, na.rm = TRUE), n_replicates = length(beta_list))
}
# END STATS


## ---- 1. Data ------------------------------------------------------------------------------------------------------
# The sample table is a TSV without quoting (its text fields contain commas): read it with quote = "" so that a stray double
# quote cannot swallow the following lines; fill = FALSE makes a misaligned row an error.
raw <- read.delim(CFG$sample_file, quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE, fill = FALSE)
required <- c("Longitude", "Latitude", CFG$group_col, "cluster_id", "alpha_total_n_reads", "unmapped_raw_counts", "total_raw_counts", CFG$biom_id)
miss <- setdiff(required, names(raw))
if (length(miss)) stop("Missing columns in ", CFG$sample_file, ": ", paste(miss, collapse = ", "),
                       "  (perialpine_alpine and cluster_id are written by microbeatlas_filter_data.R: rerun prepare_microbeatlas_data)")
# the group column must be plausible: most non-alpine samples are Lowland, some Perialpine (a table in which every non-alpine sample is
# Perialpine, or none is, comes from wrong distances to the alpine polygons: rerun prepare_microbeatlas_data, it recomputes the column)
local({
  f <- as.character(raw[[CFG$group_col]]); f <- f[f %in% c("Lowland", "Perialpine")]
  if (!(length(f) > 0 && any(f == "Lowland") && mean(f == "Perialpine") < 0.5))
    stop("the column ", CFG$group_col, " is not plausible (non-alpine samples: ", sum(f == "Lowland"), " Lowland, ", sum(f == "Perialpine"),
         " Perialpine): rerun prepare_microbeatlas_data so that it is recomputed")
})
meta <- raw %>%
  transmute(id = .data[[CFG$biom_id]], group = factor(.data[[CFG$group_col]], levels = CFG$groups),
            lon = Longitude, lat = Latitude, cl = as.character(cluster_id),
            reads = alpha_total_n_reads, log10_reads = log10(alpha_total_n_reads),
            unmapped = unmapped_raw_counts / total_raw_counts) %>%
  filter(is.finite(lon), is.finite(lat), abs(lon) <= 180, abs(lat) <= 90, is.finite(log10_reads), reads >= CFG$rare_depth,
         is.finite(unmapped), !is.na(group), !is.na(id))
no_cl <- is.na(meta$cl) | !nzchar(meta$cl); meta$cl[no_cl] <- paste0("none", seq_len(sum(no_cl)))      # a sample without a cluster is its own cluster
meta$zone <- cut(meta$lat, c(-Inf, CFG$lat_breaks, Inf), labels = CFG$zone_names, right = FALSE)
cat("\nSamples and clusters per group:\n")
print(meta %>% group_by(group) %>% summarise(samples = n(), clusters = n_distinct(cl), .groups = "drop") %>% as.data.frame())
cat("\nClusters per group and latitude zone:\n")
print(as.data.frame.matrix(tapply(meta$cl, list(zone = meta$zone, group = meta$group), function(x) length(unique(x)))))
if (any(table(meta$group) == 0)) stop("a group has no samples: ", paste(names(table(meta$group))[table(meta$group) == 0], collapse = ", "))

## ---- 2. Replicate subsamples: the same number of samples in every group ----------------------------------------------
# One random sample per cluster (a different one in every replicate), at most one per thin_m metres, then every group is cut down to
# the same n, evenly over latitude / longitude cells: n = the smallest group (the group that limits the comparison is reported).
# A large group (lowland: thousands of clusters) is first cut evenly to 3 x n_per_group: thinning is slow on many samples and only
# n_per_group of them are needed (the cut is even, so it does not favour any region).
draw_reps <- function(meta_a, label, seed0) {
  cl_n <- tapply(meta_a$cl, meta_a$group, function(x) length(unique(x))); cl_n[is.na(cl_n)] <- 0
  min_cl <- if (label == "all") 2 else CFG$min_zone_clusters
  use <- CFG$groups[cl_n[CFG$groups] >= min_cl]
  if (label == "all" && length(use) < length(CFG$groups)) stop("a group has fewer than 2 clusters")
  if (!"Alpine" %in% use || length(use) < 2) {
    message(sprintf("  %s: not analysed (groups with at least %d clusters: %s; the alpine group and one more are needed)", label, min_cl, if (length(use)) paste(use, collapse = ", ") else "none"))
    return(NULL)
  }
  sel <- lapply(seq_len(CFG$n_rep), function(r) {
    set.seed(seed0 + r)
    pools <- lapply(use, function(g) {
      x <- meta_a[meta_a$group == g, ]; x <- x[sample.int(nrow(x)), ]; x <- x[!duplicated(x$cl), ]
      if (nrow(x) > 3 * CFG$n_per_group) x <- x[even_order(x$lat, x$lon)[seq_len(3 * CFG$n_per_group)], ]
      x[thin_sites(x$lat, x$lon, CFG$thin_m), ]
    })
    names(pools) <- use
    n <- min(CFG$n_per_group, vapply(pools, nrow, 1))
    out <- bind_rows(lapply(pools, function(x) x[even_order(x$lat, x$lon)[seq_len(n)], ]))
    attr(out, "pool_sizes") <- vapply(pools, nrow, 1)
    out
  })
  list(sel = sel, groups = use)
}
analyses <- list(all = draw_reps(meta, "all", 100))
for (z in seq_along(CFG$zone_names)) analyses[[CFG$zone_names[z]]] <- draw_reps(meta[meta$zone == CFG$zone_names[z], ], CFG$zone_names[z], 1000 * z + 100)
analyses <- analyses[!vapply(analyses, is.null, TRUE)]

cells_of <- function(x) length(unique(paste(floor(x$lat / CFG$even_cell_deg), floor(x$lon / CFG$even_cell_deg))))
samples_per_rep <- bind_rows(lapply(names(analyses), function(an) bind_rows(lapply(seq_along(analyses[[an]]$sel), function(r) {
  x <- analyses[[an]]$sel[[r]]; ps <- attr(x, "pool_sizes")
  bind_rows(lapply(analyses[[an]]$groups, function(g) data.frame(analysis = an, replicate = r, group = g, pool_one_per_cluster = unname(ps[g]), samples = sum(x$group == g),
                                                                 clusters = length(unique(x$cl[x$group == g])), cells = cells_of(x[x$group == g, ]))))
}))))
write.csv(samples_per_rep, out_t("samples_per_replicate.csv"), row.names = FALSE)
for (an in names(analyses)) {
  x1 <- analyses[[an]]$sel[[1]]; ps1 <- attr(x1, "pool_sizes"); n1 <- sum(x1$group == analyses[[an]]$groups[1])
  message(sprintf("%s: %d samples per group in every replicate (one per cluster, at least %g m apart, evenly cut over %g-degree cells); groups: %s; pool of replicate 1: %s; limiting group: %s",
                  an, n1, CFG$thin_m, CFG$even_cell_deg, paste(analyses[[an]]$groups, collapse = ", "), paste(sprintf("%s %d", names(ps1), ps1), collapse = ", "),
                  paste(names(ps1)[ps1 == min(ps1)], collapse = "/")))
  stopifnot(all(vapply(analyses[[an]]$sel, function(x) length(unique(as.integer(table(droplevels(x$group))))) == 1, TRUE)))                       # the groups really are the same size
  stopifnot(all(vapply(analyses[[an]]$sel, function(x) !any(duplicated(x$cl[x$group == "Alpine"])), TRUE)))             # one alpine sample per cluster
}
n_main <- sum(analyses$all$sel[[1]]$group == analyses$all$groups[1])
if (n_main < 100) warning("fewer than 100 samples per group: the distance decays will be poorly determined")

# Samples of the map: ALL alpine clusters, one random sample each (a different one in every replicate), at most one per thin_m metres
reg_of <- function(lat, lon) paste(floor(lat / CFG$region_deg), floor(lon / CFG$region_deg))
map_sets <- lapply(seq_len(CFG$n_rep), function(r) {
  set.seed(5000 + r)
  a <- meta[meta$group == "Alpine", ]; a <- a[sample.int(nrow(a)), ]; a <- a[!duplicated(a$cl), ]; a <- a[thin_sites(a$lat, a$lon, CFG$thin_m), ]
  if (nrow(a) > CFG$map_max_samples) a <- a[sample.int(nrow(a), CFG$map_max_samples), ]
  a$region <- reg_of(a$lat, a$lon); a
})

all_ids <- unique(c(unlist(lapply(analyses, function(a) unlist(lapply(a$sel, `[[`, "id")))), unlist(lapply(map_sets, `[[`, "id"))))
message(length(all_ids), " distinct samples across all analyses and replicates")

to_sparse <- function(m) {
  if (inherits(m, "simple_triplet_matrix"))
    Matrix::sparseMatrix(i = m$i, j = m$j, x = m$v, dims = c(m$nrow, m$ncol),
                         dimnames = list(m$dimnames[[1]], m$dimnames[[2]]))
  else as(m, "dgCMatrix")
}
cache_file <- file.path(out_folder, "subsample_counts.rds")
if (file.exists(cache_file) && all(all_ids %in% colnames(X_cached <- readRDS(cache_file)))) {
  X_all <- X_cached; message("Reusing cached counts: ", cache_file)
} else {
  message("Reading BIOM (slow step): ", CFG$biom_file)
  biom <- rbiom::read_biom(CFG$biom_file)
  ids_ok <- intersect(all_ids, biom$samples)
  if (length(ids_ok) < length(all_ids)) warning(length(all_ids) - length(ids_ok), " selected samples are not in the BIOM")
  biom <- biom[ids_ok]
  X_all <- to_sparse(biom$counts)                                  # OTUs x samples
  X_all <- X_all[rownames(X_all) != "Unmapped", , drop = FALSE]
  X_all <- X_all[Matrix::rowSums(X_all) > 0, , drop = FALSE]
  rm(biom); invisible(gc())
  saveRDS(X_all, cache_file)
}
message(sprintf("Count matrix: %d OTUs x %d samples", nrow(X_all), ncol(X_all)))

# Rarefy every sample (column) of a sparse count matrix to `depth` reads, without replacement. All samples must have >= depth reads.
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

## ---- 3. Random pairs and their dissimilarities for one replicate -------------------------------------------------------
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
make_pairs <- function(S0, seed) {
  set.seed(seed)
  S <- S0 %>% filter(id %in% colnames(X_all))
  X <- X_all[, S$id, drop = FALSE]
  tot <- Matrix::colSums(X)
  ok <- tot >= CFG$rare_depth
  S <- S[ok, ]; X <- X[, ok, drop = FALSE]; tot <- tot[ok]

  # up to n_pairs pairs per group (indices into S), drawn at random WITHIN geographic distance classes (random index pairs, not all pairs)
  draw_pairs <- function(i1) {
    M <- 300L * CFG$n_pairs                                                    # random candidate pairs
    a <- i1[sample.int(length(i1), M, TRUE)]; b <- i1[sample.int(length(i1), M, TRUE)]
    key <- pmin(a, b) * 1e6 + pmax(a, b)
    keep <- a != b & !duplicated(key); a <- a[keep]; b <- b[keep]
    d <- geosphere::distHaversine(cbind(S$lon[a], S$lat[a]), cbind(S$lon[b], S$lat[b])) / 1000          # km
    strata <- cut(d, CFG$dist_breaks_km, labels = FALSE, include.lowest = TRUE)
    members <- split(seq_along(strata), strata)
    alloc <- allocate_even(lengths(members), CFG$n_pairs)
    pick <- unlist(Map(function(ix, k) if (k >= length(ix)) ix else ix[sample.int(length(ix), k)], members, alloc), use.names = FALSE)
    data.frame(i = a[pick], j = b[pick])
  }
  pairs <- bind_rows(lapply(unique(as.character(S$group)), function(g) {
    ix <- which(S$group == g)
    if (length(ix) < 2) return(NULL)
    cbind(group = g, draw_pairs(ix), stringsAsFactors = FALSE)
  }))
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

  pairs$geo_km <- geosphere::distHaversine(cbind(S$lon[pairs$i], S$lat[pairs$i]), cbind(S$lon[pairs$j], S$lat[pairs$j])) / 1000
  pairs$d_depth <- abs(S$log10_reads[pairs$i] - S$log10_reads[pairs$j])   # difference in sequencing depth
  pairs$d_unm   <- abs(S$unmapped[pairs$i] - S$unmapped[pairs$j])          # difference in unmapped proportion
  pairs$a <- S$id[pairs$i]; pairs$b <- S$id[pairs$j]                       # the two samples of a pair (for the dyadic standard errors)
  pairs[, c("group", "bray", "jaccard", "geo_km", "d_depth", "d_unm", "a", "b")]
}

# Map: for one replicate, the mean Bray-Curtis dissimilarity of the alpine samples of every region to the alpine samples of OTHER regions,
# raw and after removing the alpine distance decay (residual from one line over all between-region pairs)
map_rep <- function(S0, r) {
  set.seed(7000 + r)
  S <- S0 %>% filter(id %in% colnames(X_all)); X <- X_all[, S$id, drop = FALSE]; tot <- Matrix::colSums(X)
  S <- S[tot >= CFG$rare_depth, ]; X <- X[, tot >= CFG$rare_depth, drop = FALSE]; tot <- tot[tot >= CFG$rare_depth]
  D <- as.matrix(vegan::vegdist(t(as.matrix(X %*% Matrix::Diagonal(x = 1 / tot))), method = "bray"))
  G <- geosphere::distm(cbind(S$lon, S$lat), fun = geosphere::distHaversine) / 1000
  other <- outer(S$region, S$region, "!=")
  fit <- lm(D[other] ~ log10(pmax(G[other], 0.1)))
  Rm <- matrix(NA_real_, nrow(D), ncol(D)); Rm[other] <- resid(fit)
  regs <- names(which(table(S$region) >= CFG$min_region_n))
  bind_rows(lapply(regs, function(g) { i <- which(S$region == g); j <- which(S$region != g)
    data.frame(replicate = r, region = g, n_samples = length(i), raw = mean(D[i, j]), residual = mean(Rm[i, j]), mean_km = mean(G[i, j])) }))
}

# All replicates of all analyses. Saved (pairs_all.rds) and reused on a rerun if the settings and the drawn samples are the same, so a
# problem in the models does not cost the distance computation again.
pairs_key <- list(dist_breaks = CFG$dist_breaks_km, n_rep = CFG$n_rep, n_per_group = CFG$n_per_group, n_pairs = CFG$n_pairs,
                  rare_depth = CFG$rare_depth, n_ids = length(all_ids), ids_chars = sum(nchar(all_ids)), thin_m = CFG$thin_m,
                  even_cell_deg = CFG$even_cell_deg, lat_breaks = CFG$lat_breaks, min_zone_clusters = CFG$min_zone_clusters,
                  region_deg = CFG$region_deg, min_region_n = CFG$min_region_n, analyses = names(analyses), test = TEST_MODE, design = "decay-v2")
pairs_file <- file.path(out_folder, "pairs_all.rds")
PP <- NULL
if (REUSE_FITS && file.exists(pairs_file)) {
  cached <- readRDS(pairs_file)
  if (identical(cached$key, pairs_key)) { PP <- cached$data; message("Reusing saved pairs: ", pairs_file) }
}
if (is.null(PP)) {
  PP <- list(pairs = list(), map = NULL)
  for (an in names(analyses)) {
    pr <- list()
    for (r in seq_len(CFG$n_rep)) {
      message(an, ", replicate ", r, ": computing dissimilarities")
      x <- make_pairs(analyses[[an]]$sel[[r]], 200 + r + 100 * match(an, names(analyses)))
      if (is.null(x)) { message("  no pairs drawn: skipping replicate"); next }
      x$rep_id <- r; pr[[length(pr) + 1]] <- x
    }
    PP$pairs[[an]] <- bind_rows(pr)
  }
  message("Map: dissimilarities among all alpine clusters (", CFG$n_rep, " replicates)")
  PP$map <- bind_rows(lapply(seq_len(CFG$n_rep), function(r) map_rep(map_sets[[r]], r)))
  saveRDS(list(key = pairs_key, data = PP), pairs_file)
}

## ---- 4. One analysis: models, tests, curves ---------------------------------------------------------------------------------
# dissimilarity ~ group * f(log10 distance) + depth + unmapped differences, one model per replicate and metric. Predictions at median
# depth / unmapped differences, the same for all groups. Tests: dyadic-robust contrasts between groups (see the header).
cmp_all <- list(c("Alpine", "Lowland"), c("Perialpine", "Lowland"), c("Alpine", "Perialpine"))      # first minus second
short <- c(Alpine = "A", Perialpine = "P", Lowland = "L")
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", " ns"))))
analyse <- function(P, groups_a, label) {
  P <- as.data.frame(P); P$group <- factor(P$group, levels = groups_a); P$log10_geo <- log10(pmax(P$geo_km, 0.1))
  reps <- sort(unique(P$rep_id))
  rng <- P %>% group_by(group) %>% summarise(min_km = min(geo_km), q05_km = quantile(geo_km, CFG$common_q[1]), median_km = median(geo_km),
                                             q95_km = quantile(geo_km, CFG$common_q[2]), max_km = max(geo_km), n_pairs = n(), .groups = "drop") %>% as.data.frame()
  common <- c(max(rng$q05_km), min(rng$q95_km))
  refs <- CFG$ref_km[CFG$ref_km >= common[1] & CFG$ref_km <= common[2]]
  if (!length(refs) && is.finite(common[1]) && common[2] > common[1]) refs <- signif(10^(log10(common[1]) + c(1, 2) / 3 * (log10(common[2]) - log10(common[1]))), 2)
  ref_main <- if (length(refs)) refs[which.min(abs(log10(refs) - 3))] else NA_real_
  decades <- log10(common[2]) - log10(common[1]); use_turnover <- is.finite(decades) && decades >= 0.5
  message(sprintf("%s: distance range common to all groups (%g-%g percentile): %s to %s km; reference distances: %s", label, 100 * CFG$common_q[1], 100 * CFG$common_q[2],
                  fmt_km(common[1]), fmt_km(common[2]), if (length(refs)) paste(fmt_km(refs), collapse = ", ") else "none"))
  med_depth <- median(P$d_depth); med_unm <- median(P$d_unm)
  c0 <- if (is.finite(ref_main)) log10(ref_main) else median(P$log10_geo)                      # distance is centred at the reference distance
  fml <- function(metric) as.formula(sprintf("%s ~ group * %s + d_depth + d_unm", metric, if (CFG$df_dist > 1) sprintf("ns(log10_geo, df = %d)", CFG$df_dist) else sprintf("I(log10_geo - %s)", format(c0, digits = 10))))
  nd_at <- function(m, g, x) data.frame(group = factor(g, levels = levels(m$model$group)), log10_geo = x, d_depth = med_depth, d_unm = med_unm)
  dsg <- function(m, g, x) design_at(m, nd_at(m, g, x))                                      # design rows (for the contrasts)
  pred_at <- function(m, g, x) as.numeric(predict(m, nd_at(m, g, x)))                        # predictions, vectorised (curves, adjusted values)
  grid_g <- bind_rows(lapply(groups_a, function(g) { x <- P$log10_geo[P$group == g]; data.frame(group = g, x = seq(quantile(x, 0.01), quantile(x, 0.99), length.out = 80)) }))
  cmp <- Filter(function(cp) all(cp %in% groups_a), cmp_all)
  curves <- list(); tst <- list(); adj <- list(); checks <- list(); bl <- list(); vl <- list()
  for (metric in names(metrics)) for (r in reps) {
    Pr <- P[P$rep_id == r, ]; Pr$group <- droplevels(Pr$group)
    if (nlevels(Pr$group) < length(groups_a) || any(table(Pr$group) < 50)) next
    m <- tryCatch(lm(fml(metric), data = Pr), error = function(e) NULL); if (is.null(m)) next
    V <- dyadic_vcov(model.matrix(m), resid(m), Pr$a, Pr$b); if (is.null(V)) next
    key <- paste(metric, r)
    bl[[metric]][[length(bl[[metric]]) + 1]] <- coef(m); vl[[metric]][[length(vl[[metric]]) + 1]] <- V
    curves[[key]] <- data.frame(replicate = r, metric = metric, group = grid_g$group, x = grid_g$x, fit = clip01(pred_at(m, grid_g$group, grid_g$x)))
    # observed dissimilarity of every pair adjusted to the median depth / unmapped differences: residual + the group's curve at the pair's
    # distance; binned, these are the open circles that the curve has to run through
    adj[[key]] <- data.frame(metric = metric, group = as.character(Pr$group), x = Pr$log10_geo, y = Pr[[metric]] - fitted(m) + pred_at(m, as.character(Pr$group), Pr$log10_geo))
    fv <- tapply(fitted(m), Pr$group, mean); ob <- tapply(Pr[[metric]], Pr$group, mean)
    checks[[key]] <- data.frame(analysis = label, replicate = r, metric = metric, group = names(ob), observed_mean = as.numeric(ob), fitted_mean = as.numeric(fv[names(ob)]), n_pairs = as.numeric(table(Pr$group)[names(ob)]))
    one <- function(C, ...) { est <- drop(C %*% coef(m)); se <- sqrt(drop(C %*% V %*% t(C))); data.frame(replicate = r, metric = metric, estimate = est, se = se, ...) }
    for (g in groups_a) {
      for (rk in refs) tst[[length(tst) + 1]] <- one(dsg(m, g, log10(rk)), kind = "group", quantity = "level", ref_km = rk, g1 = g, g2 = NA_character_)
      if (use_turnover) tst[[length(tst) + 1]] <- one((dsg(m, g, log10(common[2])) - dsg(m, g, log10(common[1]))) / decades, kind = "group", quantity = "turnover", ref_km = NA_real_, g1 = g, g2 = NA_character_)
    }
    for (cp in cmp) {
      for (rk in refs) tst[[length(tst) + 1]] <- one(dsg(m, cp[1], log10(rk)) - dsg(m, cp[2], log10(rk)), kind = "contrast", quantity = "level", ref_km = rk, g1 = cp[1], g2 = cp[2])
      if (use_turnover) tst[[length(tst) + 1]] <- one(((dsg(m, cp[1], log10(common[2])) - dsg(m, cp[1], log10(common[1]))) - (dsg(m, cp[2], log10(common[2])) - dsg(m, cp[2], log10(common[1])))) / decades,
                                                      kind = "contrast", quantity = "turnover", ref_km = NA_real_, g1 = cp[1], g2 = cp[2])
    }
  }
  if (!length(curves)) { message("  ", label, ": no model could be fitted (needs at least 50 pairs per group in a replicate: the smallest group has ", min(table(P$group)) / length(reps), " pairs per replicate on average; raise min_zone_clusters or n_pairs, or lower the requirement)"); return(NULL) }
  tst <- bind_rows(tst)
  tst_key <- paste(tst$metric, tst$kind, tst$quantity, ifelse(is.na(tst$ref_km), "-", tst$ref_km), tst$g1, ifelse(is.na(tst$g2), "-", tst$g2))      # (split() would drop NA keys)
  pooled <- bind_rows(lapply(split(tst, tst_key), function(d)
    cbind(d[1, c("metric", "kind", "quantity", "ref_km", "g1", "g2")], pool_replicates(d$estimate, d$se))))
  rownames(pooled) <- NULL
  pooled$p_adj <- NA_real_; ic <- pooled$kind == "contrast"; pooled$p_adj[ic] <- p.adjust(pooled$p[ic], "BH")          # Benjamini-Hochberg over all contrasts of the analysis
  pooled$analysis <- label
  # the interaction model's own coefficients (pooled over the replicates) and the joint tests (only for the straight-line model: the
  # coefficients of a spline are not interpretable one by one)
  coefs <- NULL; omni <- NULL
  if (CFG$df_dist == 1) {
    ref_g <- levels(P$group)[1]; c0_km <- 10^c0
    lab <- function(nm) {
      g <- sub("^group", "", sub(":.*$", "", nm)); is_i <- grepl(":", nm); is_g <- grepl("^group", nm)
      ifelse(nm == "(Intercept)", sprintf("%s: dissimilarity at %s km", ref_g, fmt_km(c0_km)),
        ifelse(is_g & !is_i, sprintf("%s minus %s: difference at %s km (level)", g, ref_g, fmt_km(c0_km)),
          ifelse(is_g & is_i, sprintf("%s minus %s: difference in the rise per 10-fold distance (interaction)", g, ref_g),
            ifelse(grepl("^I\\(", nm), sprintf("%s: rise per 10-fold distance", ref_g), nm))))
    }
    coefs <- bind_rows(lapply(names(bl), function(mt) {
      B <- do.call(rbind, bl[[mt]]); SE <- do.call(rbind, lapply(vl[[mt]], function(V) sqrt(diag(V))))
      bind_rows(lapply(colnames(B), function(nm) cbind(analysis = label, metric = mt, term = nm, label = lab(nm), pool_replicates(B[, nm], SE[, nm]))))
    }))
    omni <- bind_rows(lapply(names(bl), function(mt) {
      nm <- names(bl[[mt]][[1]]); lev <- grep("^group[^:]*$", nm); slo <- grep("^group.*:", nm)
      bind_rows(if (length(lev)) cbind(analysis = label, metric = mt, quantity = "level", ref_km = c0_km, pool_wald(bl[[mt]], vl[[mt]], lev)),
                if (length(slo)) cbind(analysis = label, metric = mt, quantity = "slope", ref_km = NA_real_, pool_wald(bl[[mt]], vl[[mt]], slo)))
    }))
    if (nrow(omni)) omni$p_adj <- p.adjust(omni$p, "BH")
  }
  curves <- bind_rows(curves)
  curves_sum <- curves %>% group_by(metric, group, x) %>% summarise(lo = qlo(fit), hi = qhi(fit), fit = median(fit), .groups = "drop") %>% as.data.frame()
  checks <- bind_rows(checks); checks$difference <- checks$fitted_mean - checks$observed_mean
  list(label = label, groups = groups_a, P = P, curves_sum = curves_sum, adj = bind_rows(adj), tests = pooled, coefs = coefs, omni = omni, common = common, refs = refs, ref_main = ref_main,
       use_turnover = use_turnover, rng = cbind(analysis = label, rng), checks = checks)
}
res <- list()
for (an in names(analyses)) {
  P_an <- PP$pairs[[an]]
  if (is.null(P_an) || !nrow(P_an)) next
  message(sprintf("\n=== Analysis: %s (%s) ===", an, paste(analyses[[an]]$groups, collapse = ", ")))
  print(table(group = P_an$group))
  res[[an]] <- analyse(P_an, analyses[[an]]$groups, an)
}
if (is.null(res$all)) stop("the main analysis could not be fitted")
main <- res$all

# ---- tables
tests_all <- bind_rows(lapply(res, `[[`, "tests"))
omni_all <- bind_rows(lapply(res, `[[`, "omni"))
if (nrow(omni_all)) write.csv(omni_all[, c("analysis", "metric", "quantity", "ref_km", "wald_chisq", "df", "p", "p_adj", "share_replicates_p_below_0.05", "n_replicates")], out_t("tests_joint.csv"), row.names = FALSE)
coefs_all <- bind_rows(lapply(res, `[[`, "coefs"))
if (nrow(coefs_all)) write.csv(coefs_all[, c("analysis", "metric", "term", "label", "estimate", "se", "lo", "hi", "p", "range_lo", "range_hi", "share_replicates_p_below_0.05", "n_replicates")], out_t("model_coefficients.csv"), row.names = FALSE)
write.csv(tests_all[, c("analysis", "metric", "kind", "quantity", "ref_km", "g1", "g2", "estimate", "se", "lo", "hi", "p", "p_adj", "range_lo", "range_hi",
                        "share_replicates_p_below_0.05", "n_replicates")], out_t("tests.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(names(res), function(an) cbind(analysis = an, res[[an]]$curves_sum) %>% mutate(distance_km = 10^x))), out_t("curves.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(res, `[[`, "rng")), out_t("distance_range_per_group.csv"), row.names = FALSE)
checks_all <- bind_rows(lapply(res, `[[`, "checks")); write.csv(checks_all, out_t("model_checks.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(names(res), function(an) res[[an]]$P %>% count(rep_id, group) %>% mutate(analysis = an))), out_t("pair_counts.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(names(res), function(an) bind_rows(lapply(names(metrics), function(mt)
  res[[an]]$P %>% group_by(group) %>% summarise(share_at_1 = mean(.data[[mt]] >= 1 - 1e-12), n_pairs = n(), .groups = "drop") %>% mutate(metric = mt, analysis = an))))),
  out_t("share_of_pairs_at_dissimilarity_1.csv"), row.names = FALSE)
worst <- max(abs(checks_all$difference), na.rm = TRUE)
message(sprintf("Model check: the fitted mean dissimilarity differs from the observed one by at most %.3f%s", worst, if (worst > 0.02) "  -> LARGER THAN 0.02: look at tables/model_checks.csv" else ""))

# ---- what the tests say (main analysis)
fmt <- function(x) sprintf("%+.3f", x)
if (nrow(omni_all)) {
  message("\nJoint tests, all samples (do the groups differ at all?):")
  d <- omni_all[omni_all$analysis == "all", ]
  for (i in seq_len(nrow(d))) message(sprintf("  %s | groups differ in %s: Wald chi2(%d) = %.2f, p = %.3g (adj. %.3g), significant alone in %.0f%% of the replicates", metrics_ascii[[d$metric[i]]],
                                              if (d$quantity[i] == "level") sprintf("level at %s km", fmt_km(d$ref_km[i])) else "the rise with distance (interaction)", d$df[i], d$wald_chisq[i], d$p[i], d$p_adj[i], 100 * d$share_replicates_p_below_0.05[i]))
}
if (nrow(coefs_all)) {
  message("\nCoefficients of the interaction model, all samples (tables/model_coefficients.csv), Bray-Curtis:")
  d <- coefs_all[coefs_all$analysis == "all" & coefs_all$metric == "bray", ]
  for (i in seq_len(nrow(d))) message(sprintf("  %-95s %s [%s, %s], p %s", d$label[i], fmt(d$estimate[i]), fmt(d$lo[i]), fmt(d$hi[i]), format.pval(d$p[i], digits = 2, eps = 1e-16)))
}
message("\nDo mountains look like islands of diversity? Contrasts between groups, all samples (pooled estimate [95% CI], adjusted p, share of replicates with p < 0.05):")
for (mt in names(metrics)) {
  message("  ", metrics_ascii[[mt]], ":")
  d <- main$tests[main$tests$metric == mt & main$tests$kind == "contrast", ]
  for (i in seq_len(nrow(d))) message(sprintf("    %s | %s minus %s: %s [%s, %s], p = %.3g (adj. %.3g), %.0f%% of replicates",
                                              if (d$quantity[i] == "level") sprintf("dissimilarity at %s km", fmt_km(d$ref_km[i])) else "turnover per 10-fold distance", d$g1[i], d$g2[i],
                                              fmt(d$estimate[i]), fmt(d$lo[i]), fmt(d$hi[i]), d$p[i], d$p_adj[i], 100 * d$share_replicates_p_below_0.05[i]))
}
message("(> 0 for alpine minus lowland = alpine samples MORE dissimilar at the same distance / dissimilarity rising faster with distance = the island pattern)")

## ---- 5. Map summary ----------------------------------------------------------------------------------------------------------
alp_all <- meta[meta$group == "Alpine", ]; alp_all$region <- reg_of(alp_all$lat, alp_all$lon)
cent <- alp_all %>% group_by(region) %>% summarise(lon = mean(lon), lat = mean(lat), n_clusters = n_distinct(cl), .groups = "drop")
map_sum <- PP$map %>% group_by(region) %>% summarise(raw = median(raw), raw_lo = qlo(raw), raw_hi = qhi(raw), residual = median(residual), residual_lo = qlo(residual), residual_hi = qhi(residual),
                                                      mean_km = median(mean_km), n_samples = median(n_samples), n_replicates = n(), .groups = "drop") %>%
  left_join(cent, by = "region") %>% as.data.frame()
write.csv(map_sum, out_t("map_regions.csv"), row.names = FALSE)
message(sprintf("Map: %d alpine regions (%g-degree cells with at least %d alpine samples); distinctiveness (distance-adjusted residual) from %+.3f to %+.3f",
                nrow(map_sum), CFG$region_deg, CFG$min_region_n, min(map_sum$residual), max(map_sum$residual)))

## ---- 6. Figures -----------------------------------------------------------------------------------------------------------------
x_scale <- scale_x_continuous(breaks = log10(c(1, 10, 100, 1000, 10000)), labels = c("1", "10", "100", "1,000", "10,000"))
col_scales <- list(scale_colour_manual(values = CFG$cols, name = NULL, breaks = CFG$groups, labels = c(Lowland = "Lowland", Perialpine = "Peri-alpine", Alpine = "Alpine")),
                   scale_fill_manual(values = CFG$cols, guide = "none"),
                   guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1, linewidth = 1.2, linetype = "solid"))))
set.seed(5)

# the test results of one analysis and metric as a few lines of text for a panel: A = alpine, P = peri-alpine, L = lowland; stars = adjusted p
annot_text <- function(rs, metric) {
  d <- rs$tests[rs$tests$metric == metric & rs$tests$kind == "contrast", ]
  line <- function(dd) paste(sprintf("%s-%s %s %s", short[dd$g1], short[dd$g2], fmt(dd$estimate), stars(dd$p_adj)), collapse = "   ")
  out <- character(0)
  om <- rs$omni[rs$omni$metric == metric, ]
  if (!is.null(om) && nrow(om)) out <- c(out, paste0("Groups differ in level at ", fmt_km(om$ref_km[om$quantity == "level"][1]), " km: p = ", signif(om$p[om$quantity == "level"][1], 2), "   in rise with distance: p = ", signif(om$p[om$quantity == "slope"][1], 2)))
  lv <- d[d$quantity == "level" & d$ref_km %in% rs$ref_main, ]; if (nrow(lv)) out <- c(out, paste0("At ", fmt_km(rs$ref_main), " km:  ", line(lv)))
  tv <- d[d$quantity == "turnover", ];                          if (nrow(tv)) out <- c(out, paste0("Rise per 10-fold distance:  ", line(tv)))
  paste(out, collapse = "\n")
}

decay_panel <- function(rs, metric, tag, title = NULL, show_x = TRUE) {
  cv <- rs$curves_sum[rs$curves_sum$metric == metric, ]; cv$group <- factor(cv$group, levels = CFG$groups)
  bn <- rs$adj[rs$adj$metric == metric, ] %>% group_by(group) %>% mutate(bin = ntile(x, 8)) %>% group_by(group, bin) %>%
    summarise(x = mean(x), y = mean(y), .groups = "drop") %>% as.data.frame(); bn$group <- factor(bn$group, levels = CFG$groups)
  pts <- rs$P[sample.int(nrow(rs$P), min(nrow(rs$P), CFG$points_max)), ]; pts$group <- factor(pts$group, levels = CFG$groups)
  ylo <- max(0, floor((min(bn$y, cv$lo, na.rm = TRUE) - 0.06) * 20) / 20)               # zoom: pairs below the window are not drawn (they are in the models)
  ggplot() +
    geom_point(data = pts, aes(log10_geo, .data[[metric]], colour = group), size = 0.3, alpha = CFG$point_alpha, shape = 16) +
    geom_ribbon(data = cv, aes(x, ymin = lo, ymax = hi, fill = group), alpha = 0.2) +
    geom_line(data = cv, aes(x, fit, colour = group), linewidth = 1.1) +
    geom_point(data = bn, aes(x, y, colour = group), shape = 21, fill = "white", size = 2.6, stroke = 1.1) +
    annotate("text", x = -Inf, y = -Inf, label = annot_text(rs, metric), hjust = -0.03, vjust = -0.5, size = 2.9, lineheight = 1.1) +
    x_scale + scale_y_continuous(breaks = function(l) { b <- scales::extended_breaks()(l); b[b <= 1 + 1e-9] }) + coord_cartesian(ylim = c(ylo, 1.005)) + col_scales +
    labs(tag = tag, title = title, x = if (show_x) "Geographic distance between samples (km)" else NULL, y = paste(metrics[[metric]], "dissimilarity")) + theme_fig +
    theme(plot.title = element_text(size = 10, face = "bold"))
}
pA <- decay_panel(main, "bray", "A"); pB <- decay_panel(main, "jaccard", "B")

# map panels: alpine regions coloured by their distinctiveness from the other alpine regions
world <- map_data("world")
map_panel <- function(value, tag, title, lab, diverging, show_size = TRUE) {
  d <- map_sum[order(abs(map_sum[[value]])), ]; d$v <- d[[value]]
  lim <- max(abs(d$v), na.rm = TRUE)
  ggplot() +
    geom_polygon(data = world, aes(long, lat, group = group), fill = "grey93", colour = "grey78", linewidth = 0.15) +
    geom_point(data = d, aes(lon, lat, colour = v, size = n_clusters), alpha = 0.9) +
    (if (diverging) scale_colour_gradient2(low = "#2166AC", mid = "grey92", high = "#B2182B", midpoint = 0, limits = c(-lim, lim), name = lab)
     else scale_colour_viridis_c(option = "magma", direction = -1, name = lab)) +
    scale_size_continuous(range = c(1.8, 6), name = "Alpine clusters") +
    guides(colour = guide_colourbar(title.position = "top", barwidth = unit(2.4, "in"), barheight = unit(0.13, "in")),
           size = if (show_size) guide_legend(title.position = "top", nrow = 1) else "none") +
    coord_quickmap(xlim = c(-180, 180), ylim = c(-58, 78), expand = FALSE) +
    labs(tag = tag, title = title, x = NULL, y = NULL) + theme_fig +
    theme(legend.position = "bottom", legend.box = "horizontal", legend.title = element_text(size = 9), legend.text = element_text(size = 8),
          plot.title = element_text(size = 10, face = "bold"), panel.grid = element_blank(), panel.background = element_rect(fill = "white"))
}
pC <- map_panel("raw", "C", "How different is the alpine community of a region from the other alpine regions?", "Bray-Curtis dissimilarity", FALSE, show_size = TRUE)
pD <- map_panel("residual", "D", "... compared with what their distance explains (> 0: more distinct than expected)", "Residual dissimilarity (observed - expected from distance)", TRUE, show_size = FALSE)

stat_txt <- paste("Tests: a linear model with a group x distance interaction per replicate (distance centred at the reference distance). Joint tests first (do the groups differ in level, and in the rise with distance?), then the contrasts between groups (stars; A = alpine, P = peri-alpine, L = lowland) with standard errors that are",
                  "cluster-robust for pair data (pairs sharing a sample are dependent; Aronow et al. 2015), pooled over the replicates (mean estimate; root mean square of the standard errors), p-values adjusted over",
                  "all contrasts (Benjamini-Hochberg); *** p < 0.001, ** p < 0.01, * p < 0.05, ns. The share of replicates in which each test alone is significant is in tables/tests.csv.")
caption <- paste(
  "Distance decay of community dissimilarity in alpine, peri-alpine (lowland within 20 km of an alpine area) and lowland soils: do mountains behave like islands of diversity (alpine samples more dissimilar at the same distance, or dissimilarity rising faster with distance)?",
  sprintf("A, B: Bray-Curtis and Jaccard dissimilarity between pairs of samples of the same group against geographic distance (linear models%s). Lines = median over %d replicates, ribbons = range between replicates (2.5-97.5%%), open circles = observed dissimilarity adjusted to the median differences in sequencing depth and unmapped proportion (mean in 8 bins of distance).",
          if (CFG$df_dist > 1) sprintf(", natural spline with %d df", CFG$df_dist) else ", one straight line per group", length(PP$map$replicate[!duplicated(PP$map$replicate)])),
  sprintf("Each replicate: %d samples per group (one random sample per BioProject / location cluster, at least %g m apart, cut down evenly over %g-degree latitude / longitude cells) and up to %d pairs per group (random within classes of geographic distance).", n_main, CFG$thin_m, CFG$even_cell_deg, CFG$n_pairs),
  stat_txt,
  sprintf("C, D: alpine regions (%g-degree cells with at least %d alpine samples; point size = number of alpine clusters): mean Bray-Curtis dissimilarity of the alpine samples of a region to the alpine samples of all other regions (C) and the same after removing the alpine distance decay (D: residual from one line over all between-region pairs; red = more distinct than the distance to the other regions explains). Median over replicates (one random alpine sample per cluster).",
          CFG$region_deg, CFG$min_region_n))
assemble_fig(file.path(out_folder, "fig3_distance_decay.png"), list(list(pA, pB), list(pC, pD)), get_legend(pA), caption, width = 14, height = 12.5, row_heights = c(1, 0.85))
saveRDS(list(A = pA, B = pB, C = pC, D = pD), file.path(out_folder, "fig3_panels.rds"))

# Supplementary: the same analysis within latitude zones
zone_res <- res[names(res) != "all"]
if (length(zone_res)) {
  zpanels <- lapply(c("bray", "jaccard"), function(mt) lapply(seq_along(CFG$zone_names), function(z) {
    rs <- zone_res[[CFG$zone_names[z]]]
    if (is.null(rs)) ggplot() + annotate("text", x = 0, y = 0, label = "not enough samples", colour = "grey40") + theme_void() +
      labs(title = if (mt == "bray") CFG$zone_names[z] else NULL) + theme(plot.title = element_text(size = 10, face = "bold", hjust = 0))
    else decay_panel(rs, mt, "", if (mt == "bray") sprintf("%s  (%d samples per group)", CFG$zone_names[z], sum(analyses[[CFG$zone_names[z]]]$sel[[1]]$group == rs$groups[1])) else NULL)
  }))
  leg_src <- zone_res[[which.max(vapply(zone_res, function(x) length(x$groups), 1))]]          # the legend of the zone with the most groups
  zleg <- get_legend(decay_panel(leg_src, "bray", ""))
  zcap <- paste("Supplementary Figure 5. The distance decay of Figure 3 within latitude zones (columns; the Southern hemisphere, the Northern hemisphere between 0 and 60 degrees = temperate, above 60 degrees = polar): Bray-Curtis (top) and Jaccard (bottom) dissimilarity between pairs of samples of the same group against geographic distance.",
                "Every zone has its own matched groups (a group needs at least", CFG$min_zone_clusters, "clusters in the zone; groups and sample sizes are given above the panels), replicates, curves and tests. Lines = median over replicates, ribbons = range between replicates, open circles = observed dissimilarity adjusted to the median depth / unmapped differences.",
                stat_txt, "Zones in which the alpine group and one more group do not have enough clusters are left out (grey 'not enough samples').")
  assemble_fig(file.path(out_folder, "figS5_latitude_zones.png"), c(zpanels[1], zpanels[2]), zleg, zcap, width = 17, height = 9.5)
}

message("Done. Outputs in ", out_folder)
writeLines(capture.output(sessionInfo()), file.path(out_folder, "sessionInfo.txt"))
file.create(file.path(out_folder, ".done"))   # flag file for Snakemake