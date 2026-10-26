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
#     America) dominates a group. The replicates only redraw which sample represents a cluster and which pairs are used.
#   * Up to n_pairs pairs per group and replicate, drawn at random WITHIN classes of geographic distance (dist_breaks_km), shared as
#     evenly as possible over the classes: random pairs are mostly thousands of km apart, so short distances would hardly be
#     represented. Stratifying on the predictor does not bias the fitted line; it makes it well determined over the whole range.
#   * Bray-Curtis (relative abundances) and Jaccard (presence/absence after rarefaction to rare_depth reads), computed only for the
#     drawn pairs from sparse matrices (checked against vegan::vegdist in TEST_MODE).
#
# MODEL AND STATISTICS (panels A, B; tables/model_coefficients.csv, tables/tests.csv)
#   ONE ordinary linear regression per metric on the pairs of ALL replicates: dissimilarity ~ group * (log10 distance - log10 of the
#   reference distance) + linear nuisance terms for the differences in sequencing depth and unmapped proportion (a straight line per group;
#   df_dist > 1 gives a natural spline instead). Distance is centred at the reference distance (1,000 km), so the coefficients read directly:
#       (Intercept)               dissimilarity of the reference group (lowland) at 1,000 km
#       group<G>                  <G> minus lowland: the difference in LEVEL at 1,000 km
#       distance                  rise of the lowland dissimilarity per 10-fold distance
#       group<G>:distance         <G> minus lowland: the difference in the RISE per 10-fold distance (the interaction)
#   Alpine - peri-alpine is the difference of the two group terms (tables/tests.csv, with the other two). The lines of the figure are this
#   regression (exactly straight), the ribbons its 95% interval, the open circles the observed dissimilarity (adjusted to the median depth /
#   unmapped differences) in 8 bins of distance.
#   Pairs that share a sample, and samples of one cluster, are not independent: an ordinary regression treats thousands of pairs as thousands
#   of observations and gives far too small p-values (simulation: 63-67% false positives at p < 0.05 for a true null). The estimates are the
#   ordinary least-squares ones; their standard errors are cluster-robust for DYADIC data (Aronow, Samii & Assenova 2015: the scores of all
#   pairs that share a CLUSTER are summed), so the number of independent units is the number of clusters, however many replicates there are.
#   In a simulation this had the correct 5% false positives and about twice the power of testing every replicate separately and averaging;
#   more replicates than ~20 add nothing (the replicates redraw from the same clusters). t-tests with (clusters - 1) degrees of freedom. The p-values of the ordinary
#   regression are also in the coefficient table (column p_naive_lm_NOT_VALID), for comparison only. No adjustment for multiple tests.
#   LIMITS: the tests cannot see spatial autocorrelation between DIFFERENT clusters (neighbouring regions resemble each other): judge
#   robustness by the size of the effects.
#
# UMAP (panel C): a two-dimensional map of the community composition: ALL alpine samples (umap_n = Inf) and, from each of the peri-alpine and the
#   lowland group, as many samples (at most; spread evenly over latitude / longitude, at most one per thin_m metres), coloured by group with 95%
#   ellipses. The relative abundances of the umap_top_otus most abundant OTUs are Hellinger-transformed (square root), reduced to umap_pcs principal
#   components and embedded with UMAP (package uwot, Euclidean distance on the components: a few minutes for ~9,000 samples). UMAP keeps the
#   NEIGHBOURHOODS of samples, not the distances between groups: the position and size of a cluster and the gaps between clusters are not
#   quantitative, the axes have no meaning, and the picture depends somewhat on umap_neighbors and umap_min_dist (fixed seed: reproducible).
#   Needs the R packages uwot and RSpectra (conda: r-uwot, r-rspectra).
#
# REGION GRID (panel D): how different are alpine communities within and between regions? The regions are the continents of the alpine samples (rule-based
#   boxes of latitude / longitude in continent_of(): check tables/region_pair_dissimilarity.csv; continents with fewer than region_min_clusters alpine
#   clusters are left out; regions with fewer than region_warn_clusters clusters are marked with a dagger: their distributions rest on few cluster pairs, and
#   a region with 2 clusters on ONE within-region pair, whose violin then shows only the variation between the samples drawn for those two clusters).
#   In each of n_rep replicates ONE random alpine sample per cluster is taken (a different one every time, at most one per thin_m metres) and the
#   Bray-Curtis dissimilarity of ALL pairs of these samples is computed. Every cell of the grid (lower triangle: region x region; the diagonal = pairs within
#   a region) is a violin of the dissimilarities of the pairs between the two regions, on ONE vertical scale for all cells, coloured by the median (blue =
#   similar, red = dissimilar). The dissimilarities pile up near 1, so the vertical scale is stretched there (region_axis = "log1m": -log10(1 - d + region_eps),
#   the same in every cell; "linear" for a linear scale); the values are the ticks at the left of every row, with a light guide line in every cell. Rows and
#   columns are ordered by hierarchical clustering (average linkage) of the median dissimilarities, so that regions with similar communities are next to
#   each other. Pairs of clusters far apart are more dissimilar anyway (the distance decay of panels A, B): the median distance of the pairs of every cell
#   is in the table.
#
# SUPPLEMENTARY (figS6_overall_dissimilarity.png): the OVERALL dissimilarity of the three groups, crude and adjusted. A, B: the pairwise
#   dissimilarities (and the geographic distances of the pairs) of every group; C, D: the differences between the groups (A - L, P - L, A - P) estimated
#   with four regressions that add what could explain them: (1) group only, (2) + the differences in sequencing depth and unmapped proportion,
#   (3) + geographic distance (one slope for all groups), (4) + a distance slope per group (the main model), evaluated at the mean distance of the
#   pairs in the range common to all groups. All with the cluster-robust standard errors of the main analysis. If a crude difference shrinks
#   when distance is added, it was a difference in where the pairs are, not in the communities; a difference that survives is real.
#   tables/overall_dissimilarity_by_group.csv, tables/overall_dissimilarity_contrasts.csv.
#
# SUPPLEMENTARY (figS5_latitude_zones.png): the same analysis within the two hemispheres (lat_breaks = 0: southern and northern
#   hemisphere): every hemisphere has its own matched groups, lines and tests. A group needs min_zone_clusters clusters in it; hemispheres
#   without enough samples are left out.
#
# FIGURE (in <out_folder>)  fig3_distance_decay.png : A, B = Bray-Curtis / Jaccard dissimilarity against distance per group (with the tests);
#   C = UMAP (alpine / peri-alpine / lowland); D = violin grid of the dissimilarities between alpine regions (continents).
# TABLES (in <out_folder>/tables)  model_coefficients.csv, tests.csv, curves.csv, umap_samples.csv, umap_summary.csv, region_pair_dissimilarity.csv, overall_dissimilarity_by_group.csv, overall_dissimilarity_contrasts.csv, model_checks.csv,
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
  lat_breaks  = 0,              # latitude zones of the supplementary figure (degrees): below 0 and above 0 = the two hemispheres
  zone_names  = c("Southern hemisphere", "Northern hemisphere"),
  min_zone_clusters = 15,       # a group needs at least this many clusters in a zone to be analysed there
  umap_n      = Inf,            # UMAP: this many alpine samples (Inf = all of them) and at most as many peri-alpine and lowland samples (evenly spaced)
  umap_top_otus = 10000,        # UMAP: the most abundant OTUs used (a dense samples x OTUs matrix: ~14,000 x 10,000 = 1.1 GB)
  umap_pcs    = 100,            # UMAP: principal components of the Hellinger-transformed data that UMAP works on
  umap_neighbors = 30,          # UMAP: size of the local neighbourhood (larger = more of the global structure, smaller = more detail)
  umap_min_dist = 0.3,          # UMAP: how tightly points may be packed (0 = tight clumps, 1 = spread out)
  umap_threads = 4,             # UMAP: threads
  umap_seed   = 1,              # UMAP: random seed (the embedding is reproducible)
  region_min_clusters = 2,      # panel D: a continent needs at least this many alpine clusters to be a region of the grid
  region_warn_clusters = 5,     # panel D: regions with fewer clusters are marked with a dagger (interpret with care)
  region_axis = "log1m",        # panel D: vertical scale of the violins: "log1m" = stretched near 1 (-log10(1 - d + region_eps)), or "linear"
  region_eps  = 0.01,           # panel D: the smaller, the more the scale is stretched near 1 (a dissimilarity of exactly 1 is shown at -log10(eps))
  points_max  = 40000,          # pairs drawn as points (random subset; the models use all pairs)
  point_alpha = 0.05,
  cols        = c(Lowland = "#7C9885", Perialpine = "#CC79A7", Alpine = "#2B6CB0"),
  cols_continent = c(Europe = "#0072B2", Asia = "#D55E00", "North America" = "#E69F00", "South America" = "#CC79A7", Africa = "#B8860B", Oceania = "#56B4E9", Antarctica = "#555555")
)
if (TEST_MODE) {
  out_folder <- paste0(out_folder, "_TEST")
  CFG$n_per_group <- 60; CFG$n_rep <- 1; CFG$n_pairs <- 300; CFG$min_zone_clusters <- 8; CFG$umap_n <- 150; CFG$umap_top_otus <- 300; CFG$umap_pcs <- 20; CFG$umap_neighbors <- 15
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
  g <- ggplotGrob(p + theme(legend.position = "top"))
  for (i in grep("guide-box", g$layout$name))
    if (!inherits(g$grobs[[i]], "zeroGrob")) return(g$grobs[[i]])
  NULL
}

# Assemble rows of panels (1 to 3 per row), a shared legend ON TOP, and a caption below into one PNG
assemble_fig <- function(file, rows, leg, caption, width, height, row_heights = rep(1, length(rows))) {
  nr <- length(rows); nc <- max(lengths(rows))
  cap_lines <- strwrap(caption, 165)
  png(file, width = width, height = height, units = "in", res = 300)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nr + 2, nc,
                        heights = unit.c(unit(0.45, "in"), unit(row_heights, "null"), unit(0.17 * length(cap_lines) + 0.2, "in")))))
  if (!is.null(leg)) { pushViewport(viewport(layout.pos.row = 1, layout.pos.col = seq_len(nc))); grid.draw(leg); popViewport() }
  for (i in seq_len(nr)) for (j in seq_along(rows[[i]])) {
    vp <- viewport(layout.pos.row = i + 1, layout.pos.col = if (length(rows[[i]]) == 1) seq_len(nc) else j)
    print(rows[[i]][[j]], vp = vp)
  }
  pushViewport(viewport(layout.pos.row = nr + 2, layout.pos.col = seq_len(nc)))
  grid.text(paste(cap_lines, collapse = "\n"), x = 0.01, just = "left", gp = gpar(fontsize = 9))
  popViewport(); dev.off()
}

# BEGIN STATS
# ---- Statistical tests for pair data ---------------------------------------------------------------------------------------------
# Pairs that share a sample or a cluster are NOT independent (a cluster with an unusual community is unusual in all its pairs), so ordinary
# standard errors of a regression on pairs are far too small. The variance used here is the cluster-robust variance for dyadic data
# (Aronow, Samii & Assenova 2015): the scores of all pairs that share a node (here: a cluster) are summed per node.
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
# Design-matrix rows of a linear model at new data (also for natural splines: the knots of the fit are used)
design_at <- function(m, newdata) {
  tt <- delete.response(terms(m))
  model.matrix(tt, model.frame(tt, newdata, xlev = m$xlevels), contrasts.arg = m$contrasts)
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

# UMAP samples: ALL alpine samples (umap_n = Inf) and the same number of lowland samples, spread evenly over latitude / longitude cells
# (at most one per thin_m metres); a smaller umap_n cuts the alpine samples down evenly as well
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
set.seed(9000)
# the samples of one group for the UMAP: an even cut over latitude / longitude cells (one per cell in turn); the non-alpine groups are first cut
# to 3 x n (thinning is slow on many samples) and thinned to at most one per thin_m metres
umap_pick <- function(pool, n, thin = TRUE) {
  if (!nrow(pool)) return(pool)
  if (thin) {
    if (nrow(pool) > 3 * n) pool <- pool[even_order(pool$lat, pool$lon)[seq_len(3 * n)], ]
    pool <- pool[thin_sites(pool$lat, pool$lon, CFG$thin_m), ]
  }
  pool[even_order(pool$lat, pool$lon)[seq_len(min(n, nrow(pool)))], ]
}
n_umap    <- min(CFG$umap_n, sum(meta$group == "Alpine"))
alp_umap  <- umap_pick(meta[meta$group == "Alpine", ], n_umap, thin = FALSE)
peri_umap <- umap_pick(meta[meta$group == "Perialpine", ], n_umap)
low_umap  <- umap_pick(meta[meta$group == "Lowland", ], n_umap)
alp_umap$continent <- continent_of(alp_umap$lat, alp_umap$lon); peri_umap$continent <- rep("Peri-alpine", nrow(peri_umap)); low_umap$continent <- rep("Lowland", nrow(low_umap))
message(sprintf("UMAP samples: %d alpine (%s), %d peri-alpine and %d lowland (evenly spaced over latitude / longitude, at most one per %g m)", nrow(alp_umap),
                paste(sprintf("%s %d", names(table(alp_umap$continent)), table(alp_umap$continent)), collapse = ", "), nrow(peri_umap), nrow(low_umap), CFG$thin_m))

# Region grid: in every replicate ONE random alpine sample per cluster (ALL clusters), at most one per thin_m metres, with its continent
region_sets <- lapply(seq_len(CFG$n_rep), function(r) {
  set.seed(7000 + r)
  a <- meta[meta$group == "Alpine", ]; a <- a[sample.int(nrow(a)), ]; a <- a[!duplicated(a$cl), ]; a <- a[thin_sites(a$lat, a$lon, CFG$thin_m), ]
  a$region <- continent_of(a$lat, a$lon); a
})
reg_n <- table(region_sets[[1]]$region); reg_keep <- names(reg_n)[reg_n >= CFG$region_min_clusters]; reg_warn <- reg_keep[reg_n[reg_keep] < CFG$region_warn_clusters]
message(sprintf("Region grid: %d alpine clusters; continents: %s; regions with at least %d clusters: %s", nrow(region_sets[[1]]), paste(sprintf("%s %d", names(reg_n), reg_n), collapse = ", "),
                CFG$region_min_clusters, if (length(reg_keep)) paste(reg_keep, collapse = ", ") else "none"))
all_ids <- unique(c(unlist(lapply(analyses, function(a) unlist(lapply(a$sel, `[[`, "id")))), alp_umap$id, peri_umap$id, low_umap$id, unlist(lapply(region_sets, `[[`, "id"))))
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
  pairs$a <- S$id[pairs$i]; pairs$b <- S$id[pairs$j]                       # the two samples of a pair
  pairs$ca <- S$cl[pairs$i]; pairs$cb <- S$cl[pairs$j]                     # and their clusters (the independent units of the standard errors)
  pairs[, c("group", "bray", "jaccard", "geo_km", "d_depth", "d_unm", "a", "b", "ca", "cb")]
}

# All replicates of all analyses. Saved (pairs_all.rds) and reused on a rerun if the settings and the drawn samples are the same, so a
# problem in the models does not cost the distance computation again.
pairs_key <- list(dist_breaks = CFG$dist_breaks_km, n_rep = CFG$n_rep, n_per_group = CFG$n_per_group, n_pairs = CFG$n_pairs,
                  rare_depth = CFG$rare_depth, n_ids = length(all_ids), ids_chars = sum(nchar(all_ids)), thin_m = CFG$thin_m,
                  even_cell_deg = CFG$even_cell_deg, lat_breaks = CFG$lat_breaks, min_zone_clusters = CFG$min_zone_clusters,
                  analyses = names(analyses), test = TEST_MODE, design = "decay-v4")
pairs_file <- file.path(out_folder, "pairs_all.rds")
PP <- NULL
if (REUSE_FITS && file.exists(pairs_file)) {
  cached <- readRDS(pairs_file)
  if (identical(cached$key, pairs_key)) { PP <- cached$data; message("Reusing saved pairs: ", pairs_file) }
}
if (is.null(PP)) {
  PP <- list(pairs = list())
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
  saveRDS(list(key = pairs_key, data = PP), pairs_file)
}

## ---- 4. One analysis: the pooled regression and its tests ------------------------------------------------------------------------
# dissimilarity ~ group * (log10 distance - log10 reference distance) + depth + unmapped differences: ONE ordinary linear regression per
# metric on the pairs of all replicates, with standard errors clustered by cluster (see the header). Predictions at median depth / unmapped
# differences, the same for all groups.
cmp_all <- list(c("Alpine", "Lowland"), c("Perialpine", "Lowland"), c("Alpine", "Perialpine"))      # first minus second
short <- c(Alpine = "A", Perialpine = "P", Lowland = "L")
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", " ns"))))
fmt <- function(x) sprintf("%+.3f", x)
analyse <- function(P, groups_a, label) {
  P <- as.data.frame(P); P$group <- factor(P$group, levels = groups_a); P$log10_geo <- log10(pmax(P$geo_km, 0.1))
  rng <- P %>% group_by(group) %>% summarise(min_km = min(geo_km), q05_km = quantile(geo_km, CFG$common_q[1]), median_km = median(geo_km),
                                             q95_km = quantile(geo_km, CFG$common_q[2]), max_km = max(geo_km), n_pairs = n(), .groups = "drop") %>% as.data.frame()
  common <- c(max(rng$q05_km), min(rng$q95_km))
  refs <- CFG$ref_km[CFG$ref_km >= common[1] & CFG$ref_km <= common[2]]
  if (!length(refs) && is.finite(common[1]) && common[2] > common[1]) refs <- signif(10^(log10(common[1]) + c(1, 2) / 3 * (log10(common[2]) - log10(common[1]))), 2)
  ref_main <- if (length(refs)) refs[which.min(abs(log10(refs) - 3))] else NA_real_              # the reference distance: the one nearest to 1,000 km inside the common range
  message(sprintf("%s: distance range common to all groups (%g-%g percentile): %s to %s km; reference distance: %s km", label, 100 * CFG$common_q[1], 100 * CFG$common_q[2],
                  fmt_km(common[1]), fmt_km(common[2]), if (is.finite(ref_main)) fmt_km(ref_main) else "none (the median distance)"))
  if (any(table(P$group) < 50)) { message("  ", label, ": fewer than 50 pairs in a group: not analysed"); return(NULL) }
  med_depth <- median(P$d_depth); med_unm <- median(P$d_unm)
  c0 <- if (is.finite(ref_main)) log10(ref_main) else median(P$log10_geo)                      # distance is centred at the reference distance
  fml <- function(metric) as.formula(sprintf("%s ~ group * %s + d_depth + d_unm", metric, if (CFG$df_dist > 1) sprintf("ns(log10_geo, df = %d)", CFG$df_dist) else sprintf("I(log10_geo - %s)", format(c0, digits = 10))))
  nd_at <- function(m, g, x) data.frame(group = factor(g, levels = levels(m$model$group)), log10_geo = x, d_depth = med_depth, d_unm = med_unm)
  pred_at <- function(m, g, x) as.numeric(predict(m, nd_at(m, g, x)))                        # predictions, vectorised (lines, adjusted values)
  grid_g <- bind_rows(lapply(groups_a, function(g) { x <- P$log10_geo[P$group == g]; data.frame(group = g, x = seq(quantile(x, 0.01), quantile(x, 0.99), length.out = 80)) }))
  curves <- list(); adj <- list(); checks <- list(); fits <- list()
  for (metric in names(metrics)) {
    m <- tryCatch(lm(fml(metric), data = P), error = function(e) NULL); if (is.null(m)) next
    V <- dyadic_vcov(model.matrix(m), resid(m), P$ca, P$cb); if (is.null(V)) next
    fits[[metric]] <- list(m = m, V = V)
    Xg <- design_at(m, nd_at(m, grid_g$group, grid_g$x)); f <- drop(Xg %*% coef(m)); se <- sqrt(rowSums((Xg %*% V) * Xg))
    curves[[metric]] <- data.frame(metric = metric, group = grid_g$group, x = grid_g$x, fit = clip01(f), lo = clip01(f - 1.96 * se), hi = clip01(f + 1.96 * se))
    # observed dissimilarity of every pair adjusted to the median depth / unmapped differences: residual + the group's line at the pair's
    # distance; binned, these are the open circles that the line has to run through
    adj[[metric]] <- data.frame(metric = metric, group = as.character(P$group), x = P$log10_geo, y = P[[metric]] - fitted(m) + pred_at(m, as.character(P$group), P$log10_geo))
    fv <- tapply(fitted(m), P$group, mean); ob <- tapply(P[[metric]], P$group, mean)
    checks[[metric]] <- data.frame(analysis = label, metric = metric, group = names(ob), observed_mean = as.numeric(ob), fitted_mean = as.numeric(fv[names(ob)]), n_pairs = as.numeric(table(P$group)[names(ob)]))
  }
  if (!length(fits)) { message("  ", label, ": no model could be fitted"); return(NULL) }
  G <- attr(fits[[1]]$V, "G"); dfree <- max(G - 1, 1)
  message(sprintf("  %s: %s pairs from %s clusters (the independent units of the tests)", label, format(nrow(P), big.mark = ","), format(G, big.mark = ",")))
  # the coefficients of the interaction model and the differences between groups (only for the straight-line model: the coefficients of a spline are
  # not interpretable one by one)
  coefs <- NULL; tests <- NULL
  if (CFG$df_dist == 1) {
    ref_g <- levels(P$group)[1]; c0_km <- 10^c0
    lab <- function(nm) {
      g <- sub("^group", "", sub(":.*$", "", nm)); is_i <- grepl(":", nm); is_g <- grepl("^group", nm)
      ifelse(nm == "(Intercept)", sprintf("%s: dissimilarity at %s km", ref_g, fmt_km(c0_km)),
        ifelse(is_g & !is_i, sprintf("%s minus %s: difference in level at %s km", g, ref_g, fmt_km(c0_km)),
          ifelse(is_g & is_i, sprintf("%s minus %s: difference in the rise per 10-fold distance (interaction)", g, ref_g),
            ifelse(grepl("^I\\(", nm), sprintf("%s: rise per 10-fold distance", ref_g), nm))))
    }
    tfun <- function(est, se) { t <- est / se; data.frame(estimate = est, se = se, lo = est - qt(0.975, dfree) * se, hi = est + qt(0.975, dfree) * se, p = 2 * pt(-abs(t), dfree)) }
    coefs <- bind_rows(lapply(names(fits), function(mt) { m <- fits[[mt]]$m; V <- fits[[mt]]$V; b <- coef(m)
      cbind(analysis = label, metric = mt, term = names(b), label = lab(names(b)), tfun(unname(b), sqrt(diag(V))), p_naive_lm_NOT_VALID = unname(summary(m)$coefficients[names(b), 4]), row.names = NULL) }))
    # differences between two groups = differences of their coefficients (the reference group has none): level at the reference distance, and rise
    term_of <- function(nm, g, qt) if (qt == "level") nm[nm == paste0("group", g)] else nm[startsWith(nm, paste0("group", g, ":"))]
    cvec <- function(nm, g1, g2, qt) { v <- setNames(numeric(length(nm)), nm); for (t in term_of(nm, g1, qt)) v[t] <- v[t] + 1; for (t in term_of(nm, g2, qt)) v[t] <- v[t] - 1; v }
    cmp <- Filter(function(cp) all(cp %in% groups_a), cmp_all)
    tests <- bind_rows(lapply(names(fits), function(mt) bind_rows(lapply(cmp, function(cp) bind_rows(lapply(c("level", "slope"), function(qt) {
      m <- fits[[mt]]$m; V <- fits[[mt]]$V; v <- cvec(names(coef(m)), cp[1], cp[2], qt)
      cbind(analysis = label, metric = mt, quantity = qt, ref_km = if (qt == "level") c0_km else NA_real_, g1 = cp[1], g2 = cp[2], tfun(sum(v * coef(m)), sqrt(drop(t(v) %*% V %*% v))))
    }))))))
  } else message("  tests need df_dist = 1 (a straight line per group)")
  curves_sum <- bind_rows(curves); checks <- bind_rows(checks); checks$difference <- checks$fitted_mean - checks$observed_mean
  list(label = label, groups = groups_a, P = P, curves_sum = curves_sum, adj = bind_rows(adj), tests = tests, coefs = coefs, common = common, ref_main = ref_main,
       n_clusters = G, rng = cbind(analysis = label, rng), checks = checks)
}
res <- list()
for (an in names(analyses)) {
  P_an <- PP$pairs[[an]]
  if (is.null(P_an) || !nrow(P_an)) next
  message(sprintf("\n=== Analysis: %s (%s) ===", an, paste(analyses[[an]]$groups, collapse = ", ")))
  print(table(group = P_an$group))
  res[[an]] <- analyse(P_an, analyses[[an]]$groups, an)
}
res <- res[!vapply(res, is.null, TRUE)]
if (is.null(res$all)) stop("the main analysis could not be fitted")
main <- res$all

# ---- tables
coefs_all <- bind_rows(lapply(res, `[[`, "coefs")); tests_all <- bind_rows(lapply(res, `[[`, "tests"))
if (nrow(coefs_all)) write.csv(coefs_all[, c("analysis", "metric", "term", "label", "estimate", "se", "lo", "hi", "p", "p_naive_lm_NOT_VALID")], out_t("model_coefficients.csv"), row.names = FALSE)
if (nrow(tests_all)) write.csv(tests_all[, c("analysis", "metric", "quantity", "ref_km", "g1", "g2", "estimate", "se", "lo", "hi", "p")], out_t("tests.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(names(res), function(an) cbind(analysis = an, res[[an]]$curves_sum) %>% mutate(distance_km = 10^x))), out_t("curves.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(res, `[[`, "rng")), out_t("distance_range_per_group.csv"), row.names = FALSE)
checks_all <- bind_rows(lapply(res, `[[`, "checks")); write.csv(checks_all, out_t("model_checks.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(names(res), function(an) res[[an]]$P %>% count(rep_id, group) %>% mutate(analysis = an))), out_t("pair_counts.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(names(res), function(an) bind_rows(lapply(names(metrics), function(mt)
  res[[an]]$P %>% group_by(group) %>% summarise(share_at_1 = mean(.data[[mt]] >= 1 - 1e-12), n_pairs = n(), .groups = "drop") %>% mutate(metric = mt, analysis = an))))),
  out_t("share_of_pairs_at_dissimilarity_1.csv"), row.names = FALSE)
worst <- max(abs(checks_all$difference), na.rm = TRUE)
message(sprintf("Model check: the fitted mean dissimilarity differs from the observed one by at most %.3f%s", worst, if (worst > 0.02) "  -> LARGER THAN 0.02: look at tables/model_checks.csv" else ""))

# ---- what the models say (main analysis)
pfmt <- function(p) format.pval(p, digits = 2, eps = 1e-16)
if (nrow(coefs_all)) for (mt in names(metrics)) {
  message(sprintf("\nCoefficients of the interaction model, all samples, %s (estimate [95%% CI], p with cluster-robust standard errors | naive lm p, NOT valid):", metrics_ascii[[mt]]))
  d <- coefs_all[coefs_all$analysis == "all" & coefs_all$metric == mt, ]
  for (i in seq_len(nrow(d))) message(sprintf("  %-92s %s [%s, %s], p %s | %s", d$label[i], fmt(d$estimate[i]), fmt(d$lo[i]), fmt(d$hi[i]), pfmt(d$p[i]), pfmt(d$p_naive_lm_NOT_VALID[i])))
}
if (nrow(tests_all)) {
  d <- tests_all[tests_all$analysis == "all" & tests_all$metric == "bray", ]
  message("\nDifferences between groups, Bray-Curtis, all samples (alpine minus peri-alpine is not a coefficient of the model):")
  for (i in seq_len(nrow(d))) message(sprintf("  %s | %s minus %s: %s [%s, %s], p %s", if (d$quantity[i] == "level") sprintf("level at %s km", fmt_km(d$ref_km[i])) else "rise per 10-fold distance", d$g1[i], d$g2[i], fmt(d$estimate[i]), fmt(d$lo[i]), fmt(d$hi[i]), pfmt(d$p[i])))
}
message("(> 0 for alpine minus lowland = alpine samples MORE dissimilar at the same distance (level) / dissimilarity rising faster with distance (rise))")

## ---- 5. UMAP ---------------------------------------------------------------------------------------------------------------------------
# The relative abundances of the umap_top_otus most abundant OTUs (as in the full community), Hellinger-transformed (square root), centred and
# reduced to umap_pcs principal components (a partial SVD, RSpectra), then embedded with UMAP (uwot; Euclidean distance on the components; the first
# two components are the starting layout). No distance matrix between all samples is computed, so it takes a few minutes.
for (pkg in c("uwot", "RSpectra")) if (!requireNamespace(pkg, quietly = TRUE)) stop("the UMAP needs the R package ", pkg, " (conda: r-", tolower(pkg), "; add it to envs/r_analyses/environment.yml)")
umap_file <- file.path(out_folder, "umap_result.rds")
umap_ids <- c(alp_umap$id, peri_umap$id, low_umap$id)
umap_key <- list(n = length(umap_ids), chars = sum(nchar(umap_ids)), top = CFG$umap_top_otus, pcs = CFG$umap_pcs, nn = CFG$umap_neighbors, md = CFG$umap_min_dist,
                 seed = CFG$umap_seed, rare = CFG$rare_depth, test = TEST_MODE)
nm <- NULL
if (REUSE_FITS && file.exists(umap_file)) { umap_cached <- readRDS(umap_file); if (identical(umap_cached$key, umap_key)) { nm <- umap_cached$nm; message("Reusing the saved UMAP: ", umap_file) } }
if (is.null(nm)) {
  S <- bind_rows(alp_umap, peri_umap, low_umap) %>% filter(id %in% colnames(X_all))
  X <- X_all[, S$id, drop = FALSE]; tot <- Matrix::colSums(X); ok <- tot >= CFG$rare_depth
  S <- S[ok, ]; X <- X[, ok, drop = FALSE]; tot <- tot[ok]
  Prel <- X %*% Matrix::Diagonal(x = 1 / tot)
  rs <- Matrix::rowSums(Prel); top <- order(rs, decreasing = TRUE)[seq_len(min(CFG$umap_top_otus, nrow(Prel)))]
  Hraw <- sqrt(as.matrix(Matrix::t(Prel[top, , drop = FALSE])))                           # samples x the most abundant OTUs, Hellinger-transformed
  kpc <- min(CFG$umap_pcs, ncol(Hraw) - 1, nrow(Hraw) - 1)
  message(sprintf("UMAP: %d samples (%d alpine, %d peri-alpine, %d lowland) x %d OTUs (%.1f%% of all reads), %d principal components",
                  nrow(Hraw), sum(S$group == "Alpine"), sum(S$group == "Perialpine"), sum(S$group == "Lowland"), ncol(Hraw), 100 * sum(rs[top]) / sum(rs), kpc))
  Hc <- sweep(Hraw, 2, colMeans(Hraw)); rm(Hraw)
  sv <- RSpectra::svds(Hc, k = kpc); pcs <- sv$u %*% diag(sv$d, nrow = kpc); rm(Hc, sv)
  emb <- uwot::umap(pcs, n_neighbors = min(CFG$umap_neighbors, nrow(pcs) - 1), min_dist = CFG$umap_min_dist, metric = "euclidean",
                    init = pcs[, 1:2], init_sdev = 1e-4, n_threads = CFG$umap_threads, seed = CFG$umap_seed, verbose = FALSE)
  nm <- data.frame(id = S$id, group = as.character(S$group), continent = S$continent, lon = S$lon, lat = S$lat, UMAP1 = emb[, 1], UMAP2 = emb[, 2])
  saveRDS(list(key = umap_key, nm = nm), umap_file)
}
nm$group <- factor(nm$group, levels = CFG$groups)
write.csv(nm, out_t("umap_samples.csv"), row.names = FALSE)
write.csv(data.frame(group = levels(nm$group), samples = as.integer(table(nm$group))), out_t("umap_summary.csv"), row.names = FALSE)

## ---- 5b. Region grid: the dissimilarities of ALL pairs of alpine clusters within and between regions -----------------------------------
# Per replicate: the Bray-Curtis dissimilarity of every pair of the alpine samples (one per cluster) of the regions kept, from sparse matrices
region_pairs <- function(S0) {
  S <- S0 %>% filter(id %in% colnames(X_all), region %in% reg_keep)
  X <- X_all[, S$id, drop = FALSE]; tot <- Matrix::colSums(X); ok <- tot >= CFG$rare_depth
  S <- S[ok, ]; X <- X[, ok, drop = FALSE]; tot <- tot[ok]
  Prel <- X %*% Matrix::Diagonal(x = 1 / tot)
  ij <- which(upper.tri(matrix(0, nrow(S), nrow(S))), arr.ind = TRUE)
  bray <- as.numeric(Matrix::colSums(abs(Prel[, ij[, 1], drop = FALSE] - Prel[, ij[, 2], drop = FALSE]))) / 2
  data.frame(ra = S$region[ij[, 1]], rb = S$region[ij[, 2]], ca = S$cl[ij[, 1]], cb = S$cl[ij[, 2]], bray = bray,
             geo_km = geosphere::distHaversine(cbind(S$lon[ij[, 1]], S$lat[ij[, 1]]), cbind(S$lon[ij[, 2]], S$lat[ij[, 2]])) / 1000, stringsAsFactors = FALSE)
}
rp <- NULL
if (length(reg_keep) >= 2) {
  rp_file <- file.path(out_folder, "region_pairs.rds")
  rp_key <- list(ids = sum(nchar(unlist(lapply(region_sets, `[[`, "id")))), n = CFG$n_rep, keep = reg_keep, rare = CFG$rare_depth, test = TEST_MODE)
  if (REUSE_FITS && file.exists(rp_file)) { rp_c <- readRDS(rp_file); if (identical(rp_c$key, rp_key)) { rp <- rp_c$rp; message("Reusing the saved region pairs: ", rp_file) } }
  if (is.null(rp)) {
    message("Region grid: dissimilarities of all pairs of alpine clusters, ", CFG$n_rep, " replicates")
    rp <- bind_rows(lapply(seq_along(region_sets), function(r) cbind(region_pairs(region_sets[[r]]), rep_id = r)))
    saveRDS(list(key = rp_key, rp = rp), rp_file)
  }
} else message("Region grid: fewer than 2 regions with at least ", CFG$region_min_clusters, " alpine clusters: panel D is left out")
grid_cells <- NULL; reg_order <- NULL
if (!is.null(rp)) {
  # every pair once, ordered so that region a comes before region b alphabetically (the grid is symmetric)
  sw <- rp$ra > rp$rb; tmp <- rp$ra[sw]; rp$ra[sw] <- rp$rb[sw]; rp$rb[sw] <- tmp
  cell_stats <- rp %>% group_by(ra, rb) %>% summarise(n_pairs = n(), n_cluster_pairs = n_distinct(paste(pmin(ca, cb), pmax(ca, cb))), median = median(bray), q25 = quantile(bray, 0.25),
                                                      q75 = quantile(bray, 0.75), median_km = median(geo_km), .groups = "drop") %>% as.data.frame()
  cell_stats$few_clusters <- cell_stats$ra %in% reg_warn | cell_stats$rb %in% reg_warn                       # a region with fewer than region_warn_clusters clusters
  write.csv(cell_stats, out_t("region_pair_dissimilarity.csv"), row.names = FALSE)
  # rows and columns in the order of a hierarchical clustering of the median dissimilarities
  Mmed <- matrix(NA_real_, length(reg_keep), length(reg_keep), dimnames = list(reg_keep, reg_keep))
  for (i in seq_len(nrow(cell_stats))) { Mmed[cell_stats$ra[i], cell_stats$rb[i]] <- cell_stats$median[i]; Mmed[cell_stats$rb[i], cell_stats$ra[i]] <- cell_stats$median[i] }
  Mfill <- Mmed; Mfill[is.na(Mfill)] <- max(Mmed, na.rm = TRUE)
  reg_order <- reg_keep[if (length(reg_keep) > 2) hclust(as.dist(Mfill), method = "average")$order else seq_along(reg_keep)]
  message(sprintf("Region grid: %s pairs of %d regions; order after clustering: %s; median dissimilarity from %.3f to %.3f", format(nrow(rp), big.mark = ","), length(reg_keep),
                  paste(reg_order, collapse = ", "), min(cell_stats$median), max(cell_stats$median)))
}

## ---- 6. Figures -----------------------------------------------------------------------------------------------------------------
x_scale <- scale_x_continuous(breaks = log10(c(1, 10, 100, 1000, 10000)), labels = c("1", "10", "100", "1,000", "10,000"))
col_scales <- list(scale_colour_manual(values = CFG$cols, name = NULL, breaks = CFG$groups, labels = c(Lowland = "Lowland", Perialpine = "Peri-alpine", Alpine = "Alpine")),
                   scale_fill_manual(values = CFG$cols, guide = "none"),
                   guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1, linewidth = 1.2, linetype = "solid"))))
set.seed(5)

# the differences between groups of one analysis and metric as two lines of text for a panel (A = alpine, P = peri-alpine, L = lowland):
# the estimates and the stars of the interaction model's terms (level at the reference distance; rise per 10-fold distance)
annot_text <- function(rs, metric) {
  d <- rs$tests[rs$tests$metric == metric, ]
  if (is.null(d) || !nrow(d)) return("")
  line <- function(dd) paste(sprintf("%s-%s %s %s", short[dd$g1], short[dd$g2], fmt(dd$estimate), stars(dd$p)), collapse = "   ")
  paste(c(paste0("Level at ", fmt_km(rs$ref_main), " km:  ", line(d[d$quantity == "level", ])), paste0("Rise per 10-fold distance:  ", line(d[d$quantity == "slope", ]))), collapse = "\n")
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

# C: all samples coloured by group (lowland first = at the bottom, alpine on top), D: the alpine samples coloured by continent; 95% ellipses of the
# categories with at least 10 samples
pt_size  <- function(n) if (n < 1500) 1.1 else 0.5
pt_alpha <- function(n) min(0.7, max(0.25, 3000 / n))                                       # visible with few points, not a blob with many
nm <- nm[order(nm$group), ]
umap_guides <- guides(colour = guide_legend(nrow = 1, override.aes = list(size = 3.5, alpha = 1, linewidth = 1)))
pC <- ggplot(nm, aes(UMAP1, UMAP2, colour = group)) +
  geom_point(size = pt_size(nrow(nm)), alpha = pt_alpha(nrow(nm)), shape = 16) +
  stat_ellipse(data = nm[nm$group %in% names(which(table(nm$group) >= 10)), ], aes(group = group), level = 0.95, type = "norm", linewidth = 0.9) +
  scale_colour_manual(values = CFG$cols, name = NULL, breaks = CFG$groups, labels = c(Lowland = "Lowland", Perialpine = "Peri-alpine", Alpine = "Alpine")) +
  umap_guides + coord_equal() + labs(tag = "C", x = "UMAP 1", y = "UMAP 2") + theme_fig +
  theme(legend.position = "top", legend.text = element_text(size = 9))
# D: the grid of violins (lower triangle, diagonal = within a region): one vertical scale for all cells, fill = median dissimilarity (blue = similar,
# red = dissimilar), rows and columns ordered by the clustering above. Drawn as polygons on a common coordinate system: cell (row i, column j) has its centre
# at x = j, y = n - i + 1. The dissimilarities pile up near 1, so the vertical scale is stretched there (region_axis = "log1m": -log10(1 - d + eps), the same
# in every cell); the dissimilarity values are the ticks at the left of every row, with a light guide line in every cell.
pD <- NULL
if (!is.null(rp)) {
  nR <- length(reg_order); lo <- floor(quantile(rp$bray, 0.005) * 20) / 20; lo <- min(lo, floor(min(cell_stats$q25) * 20) / 20 - 0.05)
  use_log <- identical(CFG$region_axis, "log1m")
  tr <- if (use_log) function(d) -log10(pmax(1 - d, 0) + CFG$region_eps) else function(d) d
  u_lo <- tr(lo); u_hi <- tr(1)
  ymap <- function(u, cy) cy + ((u - u_lo) / (u_hi - u_lo) - 0.5) * 0.9                       # a value on the (stretched) scale -> the height in the cell
  # tick values: from 1 downwards, a tick is kept if it is far enough from the previous one on the stretched scale
  cand <- if (use_log) c(1, 0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0) else sort(unique(c(1, pretty(c(lo, 1), n = 4))), decreasing = TRUE)
  cand <- cand[cand >= lo - 1e-9 & cand <= 1 + 1e-9]
  tick_vals <- numeric(0)
  for (cv in cand) if (!length(tick_vals) || (tr(tail(tick_vals, 1)) - tr(cv)) / (u_hi - u_lo) >= 0.17) tick_vals <- c(tick_vals, cv)
  tick_lab <- ifelse(tick_vals == 1, "1", ifelse(abs(tick_vals * 10 - round(tick_vals * 10)) < 1e-9, sprintf("%.1f", tick_vals), sprintf("%.2f", tick_vals)))
  polys <- list(); frames <- list(); meds <- list(); labs_n <- list(); guides_l <- list()
  for (i in seq_len(nR)) for (j in seq_len(i)) {
    a <- reg_order[j]; b <- reg_order[i]; ra <- min(a, b); rb <- max(a, b)                      # (a < b alphabetically, as stored)
    sel_c <- rp$ra == ra & rp$rb == rb; v <- rp$bray[sel_c]; cx <- j; cy <- nR - i + 1
    frames[[length(frames) + 1]] <- data.frame(xmin = cx - 0.47, xmax = cx + 0.47, ymin = cy - 0.47, ymax = cy + 0.47)
    guides_l[[length(guides_l) + 1]] <- data.frame(x = cx - 0.47, xend = cx + 0.47, y = ymap(tr(tick_vals), cy))
    if (length(v) < 5) next
    dn <- density(tr(v), from = u_lo, to = u_hi, n = 100, adjust = 1.3)                         # the density on the stretched scale
    hw <- 0.45 * dn$y / max(dn$y)                                                               # equal maximum width: the shape, not the number of pairs
    mv <- median(v)
    polys[[length(polys) + 1]] <- data.frame(cell = paste(i, j), x = c(cx - hw, rev(cx + hw)), y = c(ymap(dn$x, cy), rev(ymap(dn$x, cy))), median = mv)
    meds[[length(meds) + 1]] <- data.frame(x = cx - 0.38, xend = cx + 0.38, y = ymap(tr(mv), cy), yend = ymap(tr(mv), cy))
    labs_n[[length(labs_n) + 1]] <- data.frame(x = cx - 0.44, y = cy - 0.43, label = format(n_distinct(paste(pmin(rp$ca[sel_c], rp$cb[sel_c]), pmax(rp$ca[sel_c], rp$cb[sel_c]))), big.mark = ","))
  }
  polys <- bind_rows(polys); frames <- bind_rows(frames); meds <- bind_rows(meds); labs_n <- bind_rows(labs_n); guides_l <- bind_rows(guides_l)
  ticks <- bind_rows(lapply(seq_len(nR), function(i) data.frame(y = ymap(tr(tick_vals), nR - i + 1), label = tick_lab)))     # the values at the left of every row
  reg_lab <- sprintf("%s (%d)%s", reg_order, as.integer(reg_n[reg_order]), ifelse(reg_order %in% reg_warn, " \u2020", ""))
  note <- paste0("Bar: median\nSmall number: pairs of clusters\nBrackets: alpine clusters of the region",
                 if (any(reg_order %in% reg_warn)) sprintf("\n\u2020 fewer than %d clusters: interpret with care", CFG$region_warn_clusters) else "")      # the rest is in the caption
  pD <- ggplot() +
    geom_rect(data = frames, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "white", colour = "grey80", linewidth = 0.3) +
    geom_segment(data = guides_l, aes(x = x, xend = xend, y = y, yend = y), colour = "grey88", linewidth = 0.25, linetype = "dotted") +
    geom_polygon(data = polys, aes(x, y, group = cell, fill = median), colour = "grey25", linewidth = 0.25) +
    geom_segment(data = meds, aes(x = x, xend = xend, y = y, yend = yend), colour = "black", linewidth = 0.5) +
    geom_segment(data = ticks, aes(x = 0.49, xend = 0.53, y = y, yend = y), colour = "grey30", linewidth = 0.3) +
    geom_text(data = ticks, aes(x = 0.47, y = y, label = label), hjust = 1, size = 2.5, colour = "grey30") +
    geom_text(data = labs_n, aes(x, y, label = label), hjust = 0, size = 2.2, colour = "grey35") +
    annotate("text", x = nR + 0.45, y = nR + 0.45, hjust = 1, vjust = 1, size = 2.6, lineheight = 1.1, colour = "grey25", label = note) +
    scale_fill_gradientn(colours = c("#2166AC", "#F7F7F7", "#B2182B"), name = "Median Bray-Curtis dissimilarity", limits = range(cell_stats$median),
                         values = scales::rescale(c(min(cell_stats$median), mean(range(cell_stats$median)), max(cell_stats$median)), from = range(cell_stats$median))) +
    scale_x_continuous(breaks = seq_len(nR), labels = reg_lab, limits = c(0, nR + 0.5), expand = c(0, 0)) +
    scale_y_continuous(breaks = nR - seq_len(nR) + 1, labels = reg_lab, limits = c(0.5, nR + 0.5), expand = c(0, 0)) +
    coord_fixed() + labs(tag = "D", x = NULL, y = NULL) + theme_fig +
    guides(fill = guide_colourbar(title.position = "top", barwidth = unit(3, "in"), barheight = unit(0.13, "in"))) +
    theme(legend.position = "top", legend.title = element_text(size = 9), legend.text = element_text(size = 8), panel.grid = element_blank(), panel.border = element_blank(),
          axis.text.x = element_text(angle = 40, hjust = 1, size = 9), axis.text.y = element_text(size = 9), axis.ticks = element_blank())
}

stat_txt <- paste("Statistics: one linear regression per metric on the pairs of all replicates with a group x distance interaction (distance centred at the reference distance); the stars are the tests of its terms (level at the reference distance and rise per 10-fold distance of A = alpine, P = peri-alpine, L = lowland; A-L = alpine minus lowland), with standard errors that are cluster-robust for pair data (pairs sharing a cluster are dependent; Aronow et al. 2015; the clusters are the independent units); *** p < 0.001, ** p < 0.01, * p < 0.05, ns; not adjusted for multiple tests. All coefficients: tables/model_coefficients.csv.")
caption <- paste(
  "Distance decay of community dissimilarity in alpine, peri-alpine (lowland within 20 km of an alpine area) and lowland soils, and the ordination of alpine and lowland soils.",
  sprintf("A, B: Bray-Curtis and Jaccard dissimilarity between pairs of samples of the same group against geographic distance (one straight line per group%s). Lines = regression, ribbons = 95%% interval, open circles = observed dissimilarity adjusted to the median differences in sequencing depth and unmapped proportion (mean in 8 bins of distance).",
          if (CFG$df_dist > 1) sprintf("; natural spline with %d df", CFG$df_dist) else ""),
  sprintf("Each of %d replicates: %d samples per group (one random sample per BioProject / location cluster, at least %g m apart, cut down evenly over %g-degree latitude / longitude cells) and up to %d pairs per group (random within classes of geographic distance); the pairs of all replicates are analysed together.", CFG$n_rep, n_main, CFG$thin_m, CFG$even_cell_deg, CFG$n_pairs),
  stat_txt,
  sprintf("C: UMAP of the community composition (Hellinger-transformed relative abundances of the %d most abundant OTUs, %d principal components) of %s alpine, %s peri-alpine and %s lowland samples (the non-alpine ones spread evenly over latitude / longitude); ellipses = 95%% of the samples of a group. UMAP keeps neighbourhoods, not distances: the axes, the size of the clusters and the gaps between them are not quantitative.",
          CFG$umap_top_otus, CFG$umap_pcs, format(sum(nm$group == "Alpine"), big.mark = ","), format(sum(nm$group == "Perialpine"), big.mark = ","), format(sum(nm$group == "Lowland"), big.mark = ",")),
  if (!is.null(pD)) sprintf("D: distributions of the Bray-Curtis dissimilarities of all pairs of alpine clusters within and between regions (continents; %d replicates of one random alpine sample per cluster, pairs pooled): one violin per cell, the same vertical scale in all cells%s (values at the left of every row), coloured by the median (red = more dissimilar); rows and columns ordered by hierarchical clustering of the medians. Far-apart regions are more dissimilar anyway (A, B): the median distance of every cell is in tables/region_pair_dissimilarity.csv.%s", CFG$n_rep, if (identical(CFG$region_axis, "log1m")) ", stretched near 1" else "",
                         if (any(reg_order %in% reg_warn)) sprintf(" A dagger marks regions with fewer than %d alpine clusters (%s): their distributions rest on few pairs of clusters, and a region with 2 clusters on a single within-region pair, whose violin only shows the variation between the samples drawn for those clusters; interpret with care.", CFG$region_warn_clusters, paste(reg_order[reg_order %in% reg_warn], collapse = ", ")) else ""))
assemble_fig(file.path(out_folder, "fig3_distance_decay.png"), list(list(pA, pB), if (!is.null(pD)) list(pC, pD) else list(pC)), get_legend(pA), caption, width = 14, height = 15.5, row_heights = c(1, 1.2))
saveRDS(list(A = pA, B = pB, C = pC, D = pD), file.path(out_folder, "fig3_panels.rds"))

# Supplementary: the same analysis within the two hemispheres
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
  zcap <- paste("Supplementary Figure 5. The distance decay of Figure 3 within the two hemispheres (columns: Southern and Northern hemisphere): Bray-Curtis (top) and Jaccard (bottom) dissimilarity between pairs of samples of the same group against geographic distance.",
                "Every hemisphere has its own matched groups (a group needs at least", CFG$min_zone_clusters, "clusters in it; sample sizes are given above the panels), regression and tests. Lines = regression, ribbons = 95% interval, open circles = observed dissimilarity adjusted to the median depth / unmapped differences.",
                stat_txt, "A hemisphere in which the alpine group and one more group do not have enough clusters is left out ('not enough samples').")
  assemble_fig(file.path(out_folder, "figS5_latitude_zones.png"), c(zpanels[1], zpanels[2]), zleg, zcap, width = 12, height = 9.5)
}

## ---- 7. Supplementary: the overall dissimilarity of the groups, crude and adjusted ----------------------------------------------------
# Does the picture of the distance decay agree with a plain comparison of the dissimilarity of the groups? The pairs of the main analysis (all replicates,
# all distances); the differences between the groups under four regressions that add, one after the other, what could explain a difference.
P0 <- main$P; groups0 <- main$groups
in_common <- P0$geo_km >= main$common[1] & P0$geo_km <= main$common[2]
m0 <- mean(P0$log10_geo[in_common])                                         # the mean distance of the pairs in the range common to all groups (log10 km)
model_labels <- c("1  Crude: group only", "2  + sequencing depth, unmapped proportion", "3  + geographic distance (one slope for all)",
                  sprintf("4  + a distance slope per group (main model), at %s km", fmt_km(10^m0)))
fmls <- list(~ group, ~ group + d_depth + d_unm, ~ group + d_depth + d_unm + log10_geo, as.formula(sprintf("~ group * I(log10_geo - %s) + d_depth + d_unm", format(m0, digits = 10))))
contr <- bind_rows(lapply(names(metrics), function(mt) bind_rows(lapply(seq_along(fmls), function(k) {
  m <- lm(update.formula(fmls[[k]], as.formula(paste(mt, "~ ."))), data = P0)
  V <- dyadic_vcov(model.matrix(m), resid(m), P0$ca, P0$cb); if (is.null(V)) return(NULL)
  dfree <- max(attr(V, "G") - 1, 1); b <- coef(m); nm <- names(b)
  bind_rows(lapply(Filter(function(cp) all(cp %in% groups0), cmp_all), function(cp) {
    v <- setNames(numeric(length(nm)), nm)
    if (paste0("group", cp[1]) %in% nm) v[paste0("group", cp[1])] <- 1; if (paste0("group", cp[2]) %in% nm) v[paste0("group", cp[2])] <- -1   # the reference group has no coefficient
    est <- sum(v * b); se <- sqrt(drop(t(v) %*% V %*% v))
    data.frame(metric = mt, model_n = k, model = model_labels[k], contrast = paste(short[cp[1]], short[cp[2]], sep = " - "), estimate = est, se = se,
               lo = est - qt(0.975, dfree) * se, hi = est + qt(0.975, dfree) * se, p = 2 * pt(-abs(est / se), dfree), stringsAsFactors = FALSE)
  }))
}))))
write.csv(contr, out_t("overall_dissimilarity_contrasts.csv"), row.names = FALSE)
desc <- bind_rows(lapply(names(metrics), function(mt) P0 %>% group_by(group) %>%
  summarise(metric = mt, n_pairs = n(), n_clusters = n_distinct(c(ca, cb)), mean = mean(.data[[mt]]), median = median(.data[[mt]]), q25 = quantile(.data[[mt]], 0.25), q75 = quantile(.data[[mt]], 0.75),
            share_at_1 = mean(.data[[mt]] >= 1 - 1e-12), median_km = median(geo_km), .groups = "drop")))
write.csv(as.data.frame(desc), out_t("overall_dissimilarity_by_group.csv"), row.names = FALSE)
if (nrow(contr)) {
  message("\nOverall dissimilarity of the groups, Bray-Curtis: difference between groups (estimate [95% CI]) as the regressions add depth / unmapped, distance, and a slope per group:")
  for (cn in unique(contr$contrast)) {
    d <- contr[contr$metric == "bray" & contr$contrast == cn, ]
    message(sprintf("  %s : %s", cn, paste(sprintf("(%d) %s [%s, %s]%s", d$model_n, fmt(d$estimate), fmt(d$lo), fmt(d$hi), stars(d$p)), collapse = "   ")))
  }
}

long_d <- bind_rows(lapply(names(metrics), function(mt) data.frame(group = P0$group, metric = metrics[[mt]], y = P0[[mt]])))
ylo_d <- max(0, floor(quantile(long_d$y, 0.01) * 20) / 20)
grp_lab <- c(Lowland = "Lowland", Perialpine = "Peri-alpine", Alpine = "Alpine")
violin_box <- list(geom_violin(scale = "width", colour = "grey30", linewidth = 0.3, alpha = 0.85), geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.3),
                   scale_fill_manual(values = CFG$cols, guide = "none"), scale_x_discrete(labels = grp_lab))
pS1 <- ggplot(long_d, aes(group, y, fill = group)) + violin_box +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.4, fill = "white", colour = "black") +
  facet_wrap(~ metric, nrow = 1) + coord_cartesian(ylim = c(ylo_d, 1)) + labs(tag = "A", x = NULL, y = "Dissimilarity of pairs of samples") + theme_fig
pS2 <- ggplot(P0, aes(group, log10_geo, fill = group)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = log10(main$common[1]), ymax = log10(main$common[2]), fill = "black", alpha = 0.07) + violin_box +
  annotate("label", x = 0.55, y = log10(main$common[2]), label = "range common to all groups", hjust = 0, vjust = -0.2, size = 2.8, colour = "grey25", fill = "white", alpha = 0.8, linewidth = 0) +
  scale_y_continuous(breaks = log10(c(1, 10, 100, 1000, 10000)), labels = c("1", "10", "100", "1,000", "10,000")) +
  labs(tag = "B", x = NULL, y = "Geographic distance of the pairs (km)") + theme_fig
cols_contrast <- c("A - L" = "#2B6CB0", "P - L" = "#CC79A7", "A - P" = "#6A51A3")
forest <- function(mt, tag, ttl) {
  d <- contr[contr$metric == mt, ]; d$model <- factor(d$model, levels = rev(model_labels)); d$contrast <- factor(d$contrast, levels = names(cols_contrast))
  ggplot(d, aes(estimate, model, colour = contrast)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_linerange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.65), linewidth = 0.7) +
    geom_point(position = position_dodge(width = 0.65), size = 2.2) +
    geom_text(aes(x = hi, label = stars(p)), position = position_dodge(width = 0.65), hjust = -0.25, vjust = 0.4, size = 2.6, show.legend = FALSE) +
    scale_colour_manual(values = cols_contrast, name = NULL, labels = c("A - L" = "Alpine - lowland", "P - L" = "Peri-alpine - lowland", "A - P" = "Alpine - peri-alpine")) +
    scale_y_discrete(labels = scales::label_wrap(30)) + scale_x_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    labs(tag = tag, title = ttl, x = "Difference in dissimilarity (estimate and 95% interval)", y = NULL) + theme_fig +
    theme(plot.title = element_text(size = 10, face = "bold"), legend.position = "none")        # one legend for the figure, on top
}
pS3 <- forest("bray", "C", "Bray-Curtis"); pS4 <- forest("jaccard", "D", "Jaccard (rarefied)")
s6cap <- paste("Supplementary Figure 6. The overall dissimilarity of alpine, peri-alpine and lowland soils, crude and adjusted.",
               sprintf("A: dissimilarity of the pairs of samples of each group (all replicates and distances; violin, box = quartiles, diamond = mean; %s pairs). B: the geographic distance of the same pairs (shaded: the range common to all groups): the groups do not cover the same distances, so the crude comparison in A mixes differences in the communities with differences in where the pairs are.", format(nrow(P0), big.mark = ",")),
               "C, D: differences between the groups (A = alpine, P = peri-alpine, L = lowland) estimated by four linear regressions that add, one after the other, the differences in sequencing depth and unmapped proportion, the geographic distance (log10), and a distance slope per group (the main model of Figure 3, evaluated at the mean distance of the pairs in the common range); points = estimates, bars = 95% intervals, standard errors cluster-robust for pair data (pairs sharing a cluster are dependent; the clusters are the independent units); *** p < 0.001, ** p < 0.01, * p < 0.05, ns; not adjusted for multiple tests. A difference that shrinks towards 0 when distance is added reflects where the pairs are, not the communities.")
assemble_fig(file.path(out_folder, "figS6_overall_dissimilarity.png"), list(list(pS1, pS2), list(pS3, pS4)), get_legend(pS3 + theme(legend.position = "top")), s6cap, width = 13, height = 11.5)

message("Done. Outputs in ", out_folder)
writeLines(capture.output(sessionInfo()), file.path(out_folder, "sessionInfo.txt"))
file.create(file.path(out_folder, ".done"))   # flag file for Snakemake