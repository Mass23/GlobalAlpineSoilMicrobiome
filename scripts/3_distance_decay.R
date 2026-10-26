# =============================================================================
# Script 3 - Distance decay of community dissimilarity: are mountains "islands of diversity"?
#
# Question: do alpine soils differ MORE from each other with distance than lowland soils do (alpine samples more dissimilar
#   at the same geographic distance, and/or dissimilarity rising faster with distance)? And how different are the alpine communities of
#   different mountain ranges? Two groups: alpine, and lowland = non-alpine samples more than PERIALPINE_KM from an alpine polygon (column
#   perialpine_alpine from microbeatlas_filter_data.R; the "peri-alpine" samples closer than that are not used here).
#
# DESIGN (n_rep replicates; the same logic as script 2)
#   * Pairs of samples are WITHIN a group: the distance decay of each group.
#   * Every replicate draws, in each group, ONE random sample per cluster (cluster_id: samples of one BioProject or one location
#     are linked), at most one per thin_m metres (no pairs at 0 km from different studies), and cuts both groups down to the SAME
#     number of samples (the smaller group decides, at most n_per_group), EVENLY over latitude / longitude cells (even_cell_deg):
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
#   The lines of the figure are this regression (exactly straight), the ribbons its 95% interval, the open circles the observed dissimilarity
#   (adjusted to the median depth / unmapped differences) in n_bins bins of distance (equal numbers of pairs), pooled over all replicates; the small dots are the same
#   bin means of each single replicate (rep_dots), so that the n_rep replicates are visible (their spread is the noise of the random draws).
#   Pairs that share a sample, and samples of one cluster, are not independent: an ordinary regression treats thousands of pairs as thousands
#   of observations and gives far too small p-values (simulation: 63-67% false positives at p < 0.05 for a true null). The estimates are the
#   ordinary least-squares ones; their standard errors are cluster-robust for DYADIC data (Aronow, Samii & Assenova 2015: the scores of all
#   pairs that share a CLUSTER are summed), so the number of independent units is the number of clusters, however many replicates there are.
#   In a simulation this had the correct 5% false positives and about twice the power of testing every replicate separately and averaging.
#   t-tests with (clusters - 1) degrees of freedom. The p-values of the ordinary regression are also in the coefficient table (column
#   p_naive_lm_NOT_VALID), for comparison only. No adjustment for multiple tests.
#   LIMITS: the tests cannot see spatial autocorrelation between DIFFERENT clusters (neighbouring regions resemble each other): judge
#   robustness by the size of the effects.
#
# HOW MANY REPLICATES? The replicates add no independent information (the clusters are the independent units); they average out WHICH sample
#   represents a cluster and WHICH pairs were drawn. So the number is enough when the result stops moving: tables/replicate_convergence.csv gives
#   the alpine - lowland level and rise (Bray-Curtis) from the first k replicates and from the odd and the even replicates separately; the two halves
#   should differ by much less than the standard error (the log says by how much; the noise of the full set is about half of that difference).
#
# UMAP (panel C): a two-dimensional map of the community composition of ALL alpine samples (umap_n = Inf) and as many lowland samples (spread
#   evenly over latitude / longitude, at most one per thin_m metres). Lowland samples are grey, alpine samples are coloured by mountain range (the
#   n_top_ranges ranges with the most samples and the largest range of every continent; the others together), each range with an ellipse (the central half of its samples, robust: umap_ellipse_*), shown as small multiples (all alpine against lowland, then one panel per continent, umap_facet); the colours are one colour family per continent (Europe reds, Asia blues, North America greens, South America oranges, Oceania purples, Africa golds), the ranges of a continent being shades of it (the largest the darkest), and the legend is grouped by continent. The relative abundances of the umap_top_otus most abundant OTUs are
#   Hellinger-transformed (square root), reduced to umap_pcs principal components and embedded with UMAP (package uwot, Euclidean distance on the
#   components: a few minutes for ~9,000 samples). UMAP keeps the NEIGHBOURHOODS of samples, not the distances between groups: the position and
#   size of a cluster and the gaps between clusters are not quantitative, the axes have no meaning, and the picture depends somewhat on
#   umap_neighbors and umap_min_dist (fixed seed: reproducible). Needs the R packages uwot and RSpectra (conda: r-uwot, r-rspectra).
#
# UMAP TEST AND SPREAD (panel C text, tables/umap_similarity_test.csv, tables/umap_spread.csv, Supplementary Figure 9): are alpine samples more similar to each other
#   than lowland samples in the UMAP, and how do the spreads compare? Pairs of samples of one group from DIFFERENT clusters (a cluster is one study / place: its samples
#   resemble each other, and an alpine cluster has dozens of samples while lowland samples are spread over many clusters, so same-cluster pairs would make alpine look
#   similar trivially): in each of umap_test_rep replicates one random sample per cluster, umap_test_pairs random pairs per group. The difference alpine - lowland in the
#   mean distance of the pairs (negative: alpine samples are closer, more similar) is estimated by a linear regression with cluster-robust standard errors for pair data (the
#   clusters are the independent units), crude and at the same geographic distance (alpine regions are also closer together in space). The spread is measured on the
#   CLUSTER CENTRES (the mean position of the samples of a cluster, the clusters counting equally): their mean distance to the group centre, their total variance and the
#   area of their 95% ellipse, plus the distance of the samples to their own cluster centre and the share of the total variance that lies between clusters; the alpine / lowland
#   ratios have a 95% interval and p-value from a bootstrap that resamples whole clusters (umap_boot resamples). Everything is done in the 2D UMAP (as drawn) and, as the
#   check that matters, in the space the UMAP was made from (the umap_pcs principal components of the Hellinger-transformed data; the Euclidean distance there is the
#   Hellinger distance between communities and is not distorted like a UMAP). The UMAP axes are not metric: its distances and areas are only a description of the picture.
#
# MOUNTAIN RANGES: the GMBA Mountain Inventory v2 (Snethlage et al. 2022, Scientific Data 9:149, doi:10.1038/s41597-022-01256-y; free download from
#   earthenv.org/mountains, doi:10.48601/earthenv-t9k2-1407; put the unzipped shapefile in gmba_path). The inventory is hierarchical: its polygons of
#   all levels overlap in one file. range_level = NA uses the ranges and systems of the column Select_300 (292, made for global comparisons); a
#   number 1-10 uses all polygons of that hierarchical level (Hier_Lvl). An alpine sample is given to the smallest selected polygon that contains it,
#   or, if it lies outside all of them, to the nearest one closer than range_max_km; the rest is left out. Without the GMBA file the "ranges" are the
#   CONTINENTS of the alpine samples (rule-based boxes of latitude / longitude in continent_of()).
#
# RANGE NETWORK (panel D): a network of the mountain ranges. Nodes: all ranges with at least net_min_clusters alpine clusters (and the coloured ones); colour as in C
#   (the other ranges dark grey), size = number of alpine CLUSTERS; only the coloured ranges are named. The distance between two ranges is the MEDIAN
#   dissimilarity of all pairs of alpine clusters of the two ranges: in each of n_rep replicates ONE random alpine sample per cluster is taken (a different
#   one every time) and the Bray-Curtis dissimilarity of ALL pairs of these samples is computed. Every pair of clusters then gets ONE value, the mean over the
#   replicates (so that every pair of clusters, and with it every cluster, counts equally, however often it was drawn and however many samples it has), and
#   everything below works on these pairs of clusters. The dissimilarities are adjusted to the median differences
#   in sequencing depth and unmapped proportion (the coefficients are estimated within pairs of ranges, from the pairs of clusters) but NOT for geographic distance, so the network
#   shows how different the ranges are, distance included. The nodes are placed by classical multidimensional scaling of these medians (ranges with similar
#   communities lie close together; the log gives the share of the variation that the two axes show); the edges are, for every range, the links to its net_k
#   most similar ranges plus the minimum spanning tree (so that the network is connected), coloured by the median dissimilarity. The supplementary figure S7
#   shows the same network next to a version adjusted to the same geographic distance as well (a common slope), the other view; its layout is rotated onto
#   that of Figure 3D to be comparable. tables/range_pair_dissimilarity.csv has the raw and the adjusted medians, the median distance and the number of pairs
#   of clusters of every pair of ranges, the diagonal being the pairs within a range. Supplementary figure S8: the distributions behind the medians of the
#   coloured ranges as a grid of violins of the values of the pairs of clusters (cells with fewer than 5 pairs show the values themselves; one scale for all cells, stretched near 1: region_axis; rows and columns ordered by hierarchical clustering of the
#   medians). Ranges with fewer than region_warn_clusters clusters are marked with a dagger there: their distributions rest on few cluster pairs, and a
#   range with 2 clusters on ONE within-range pair: interpret with care.
#
# SUPPLEMENTARY (figS5_latitude_zones.png): the same analysis within the two hemispheres (lat_breaks = 0: southern and northern
#   hemisphere): every hemisphere has its own matched groups, lines and tests. A group needs min_zone_clusters clusters in it; hemispheres
#   without enough samples are left out.
# SUPPLEMENTARY (figS6_overall_dissimilarity.png): the OVERALL dissimilarity of the two groups, crude and adjusted. A, B: the pairwise
#   dissimilarities (and the geographic distances of the pairs) of every group; C, D: the difference alpine - lowland estimated with four
#   regressions that add what could explain it: (1) group only, (2) + the differences in sequencing depth and unmapped proportion, (3) + geographic
#   distance (one slope for all groups), (4) + a distance slope per group (the main model), evaluated at the mean distance of the pairs in the range
#   common to the groups. All with the cluster-robust standard errors of the main analysis. If a crude difference shrinks when distance is added,
#   it was a difference in where the pairs are, not in the communities; a difference that survives is real.
# SUPPLEMENTARY (figS7_range_dissimilarity_space.png, figS8_range_violins.png): see RANGE NETWORK.
#
# FIGURE (in <out_folder>)  fig3_distance_decay.png : A, B = Bray-Curtis / Jaccard dissimilarity against distance per group (with the tests);
#   C = UMAP (lowland grey, alpine by mountain range); D = network of the mountain ranges (median dissimilarity between them).
# TABLES (in <out_folder>/tables)  model_coefficients.csv, tests.csv, curves.csv, replicate_convergence.csv, umap_samples.csv, umap_summary.csv,
#   range_pair_dissimilarity.csv, overall_dissimilarity_by_group.csv, overall_dissimilarity_contrasts.csv, model_checks.csv, samples_per_replicate.csv,
#   pair_counts.csv, distance_range_per_group.csv, share_of_pairs_at_dissimilarity_1.csv
#
# Input : <data_folder>/sample_data_filtered.tsv (needs perialpine_alpine and cluster_id: rerun prepare_microbeatlas_data if absent),
#         the filtered BIOM written by microbeatlas_filter_data.R and (optional) the GMBA shapefile in gmba_path
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
  groups      = c("Lowland", "Alpine"),                    # the groups used, in this order (the peri-alpine samples are not used)
  n_rep       = 40,             # replicates: each draws the samples (one per cluster, even cut) and the pairs again from scratch
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
  umap_n      = Inf,            # UMAP: this many alpine samples (Inf = all of them) and as many lowland samples (evenly spaced)
  umap_top_otus = 10000,        # UMAP: the most abundant OTUs used (a dense samples x OTUs matrix: ~9,500 x 10,000 = 0.8 GB)
  umap_pcs    = 100,            # UMAP: principal components of the Hellinger-transformed data that UMAP works on
  umap_neighbors = 30,          # UMAP: size of the local neighbourhood (larger = more of the global structure, smaller = more detail)
  umap_min_dist = 0.3,          # UMAP: how tightly points may be packed (0 = tight clumps, 1 = spread out)
  umap_threads = 4,             # UMAP: threads
  umap_seed   = 1,              # UMAP: random seed (the embedding is reproducible)
  umap_facet  = TRUE,           # panel C: small multiples (all alpine vs lowland, then one panel per continent with its ranges coloured) instead of one crowded panel
  umap_ellipse_level = 0.5,     # panel C: ellipses contain this share of the samples of a range ...
  umap_ellipse_type = "t",      # ... "t" = robust (multivariate t, not pulled out by outlying studies), "norm" = Gaussian
  umap_test_rep = 20,           # UMAP test: replicates (one random sample per cluster in each)
  umap_test_pairs = 5000,       # UMAP test: random pairs of samples of different clusters per group and replicate
  umap_boot   = 1000,           # UMAP spread: bootstrap resamples of the clusters
  gmba_path   = file.path(data_folder, "GMBA"),   # panels C, D: GMBA Mountain Inventory v2, a folder with the unzipped shapefile (or the .shp itself): earthenv.org/mountains
                                                  # (without it the regions are the continents)
  range_level = NA,             # panels C, D: which GMBA polygons are the ranges: NA = the 292 of the column Select_300 (made for global comparisons); a number 1-10 = all of that Hier_Lvl
  net_min_clusters = 2,         # panel D: a range needs at least this many alpine clusters to be a node of the network (the coloured ranges are always nodes)
  net_stretch = 0.9,            # panel D: the dissimilarities between ranges are reduced by this share of the smallest one before the MDS (they lie in a narrow band, so the layout would be a blob); 0 = none
  net_k       = 3,              # panel D: every range is linked to its this many most similar ranges (plus the minimum spanning tree)
  range_max_km = 50,            # panels C, D: an alpine sample outside every polygon is given to the nearest one if it is closer than this (km); otherwise it is left out
  region_warn_clusters = 5,     # S8: ranges with fewer alpine clusters are marked with a dagger (interpret with care)
  region_axis = "log1m",        # S8: vertical scale of the violins: "log1m" = stretched near 1 (-log10(1 - d + region_eps)), or "linear"
  region_eps  = 0.01,           # S8: the smaller, the more the scale is stretched near 1 (a dissimilarity of exactly 1 is shown at -log10(eps))
  n_bins      = 20,             # panels A, B: the open circles are the means of this many bins of distance (equal numbers of pairs per bin, pooled over all replicates)
  rep_dots    = TRUE,           # panels A, B: also a small dot per replicate and bin (the n_rep replicates)
  points_max  = 40000,          # pairs drawn as points (random subset; the models use all pairs)
  point_alpha = 0.05,
  cols        = c(Lowland = "#7C9885", Alpine = "#2B6CB0"),
  n_top_ranges = 14,            # panels C, D: the ranges with the most alpine samples get a colour (and links on the map); the others are "Other ranges"
  range_each_continent = TRUE,  # ... and so does the largest range of every continent (if it has at least range_min_samples alpine samples), so that no continent is missing
  range_min_samples = 30,
  cols_continent_shades = list(Europe = c("#8B0000", "#F4796B"), Asia = c("#08306B", "#6BAED6"), "North America" = c("#00441B", "#74C476"), "South America" = c("#7F2704", "#FDAE6B"),
                               Oceania = c("#3F007D", "#B3A2D6"), Africa = c("#6B4F00", "#E6C229"), Antarctica = c("#1F3B4D", "#8FB5C9"), Other = c("#00695C", "#80CBC4")),   # panels C, D: dark and light end of the shades of the ranges of a continent (the largest range gets the darkest)
  col_other_ranges = "#5F5F5F", col_lowland = "#C9C9C9"
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
get_legend <- function(p, pos = "top") {
  g <- ggplotGrob(p + theme(legend.position = pos))
  for (i in grep("guide-box", g$layout$name))
    if (!inherits(g$grobs[[i]], "zeroGrob")) return(g$grobs[[i]])
  NULL
}

# Assemble rows of panels, a shared legend ON TOP, and a caption below into one PNG. row_widths: optional list with the relative widths of the panels of
# every row (NULL = equal widths; a row with one panel uses the full width)
assemble_fig <- function(file, rows, leg, caption, width, height, row_heights = rep(1, length(rows)), row_widths = NULL) {
  nr <- length(rows)
  cap_lines <- strwrap(caption, 165)
  png(file, width = width, height = height, units = "in", res = 300)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nr + 2, 1, heights = unit.c(unit(0.45, "in"), unit(row_heights, "null"), unit(0.17 * length(cap_lines) + 0.2, "in")))))
  if (!is.null(leg)) { pushViewport(viewport(layout.pos.row = 1)); grid.draw(leg); popViewport() }
  for (i in seq_len(nr)) {
    k <- length(rows[[i]]); w <- if (!is.null(row_widths) && length(row_widths) >= i && !is.null(row_widths[[i]])) row_widths[[i]] else rep(1, k)
    pushViewport(viewport(layout.pos.row = i + 1)); pushViewport(viewport(layout = grid.layout(1, k, widths = unit(w, "null"))))
    for (j in seq_len(k)) print(rows[[i]][[j]], vp = viewport(layout.pos.col = j))
    popViewport(2)
  }
  pushViewport(viewport(layout.pos.row = nr + 2))
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

# UMAP samples: ALL alpine samples (umap_n = Inf) and as many lowland samples, spread evenly over latitude / longitude cells (at most one per thin_m
# metres); a smaller umap_n cuts the alpine samples down evenly as well
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

set.seed(9000)
# the samples of one group for the UMAP: an even cut over latitude / longitude cells (one per cell in turn); the lowland group is first cut
# to 3 x n (thinning is slow on many samples) and thinned to at most one per thin_m metres
umap_pick <- function(pool, n, thin = TRUE) {
  if (!nrow(pool)) return(pool)
  if (thin) {
    if (nrow(pool) > 3 * n) pool <- pool[even_order(pool$lat, pool$lon)[seq_len(3 * n)], ]
    pool <- pool[thin_sites(pool$lat, pool$lon, CFG$thin_m), ]
  }
  pool[even_order(pool$lat, pool$lon)[seq_len(min(n, nrow(pool)))], ]
}
n_umap   <- min(CFG$umap_n, sum(meta$group == "Alpine"))
alp_umap <- umap_pick(meta[meta$group == "Alpine", ], n_umap, thin = FALSE)
low_umap <- umap_pick(meta[meta$group == "Lowland", ], n_umap)
# the mountain range of every alpine sample; the n_top_ranges ranges with the most samples get a colour (UMAP, map), the others are "Other ranges"
alp_umap$region <- region_of(alp_umap$lon, alp_umap$lat)
range_n <- sort(table(alp_umap$region[alp_umap$region != "Other"]), decreasing = TRUE)
top_ranges <- names(range_n)[seq_len(min(CFG$n_top_ranges, length(range_n)))]
# the continent of every range (the continent of the centre of its samples)
cen <- alp_umap[alp_umap$region != "Other", ] %>% group_by(region) %>% summarise(lon = mean(lon), lat = mean(lat), .groups = "drop") %>% as.data.frame()
cen$continent <- continent_of(cen$lat, cen$lon); cen$n <- as.integer(range_n[cen$region])
# ... plus the largest range of every continent that is not in yet, if it has enough samples
if (isTRUE(CFG$range_each_continent) && length(range_n)) {
  have <- cen$continent[cen$region %in% top_ranges]
  for (cn in setdiff(unique(cen$continent), c("Other", have))) {
    cand <- cen[cen$continent == cn & cen$n >= CFG$range_min_samples, ]
    if (nrow(cand)) { add <- cand$region[which.max(cand$n)]; top_ranges <- c(top_ranges, add); message(sprintf("  added as the largest range of %s: %s (%d samples)", cn, add, cand$n[which.max(cand$n)])) }
  }
}
top_ranges <- top_ranges[order(-as.integer(range_n[top_ranges]))]
rest <- setdiff(names(range_n), top_ranges)
if (length(rest)) message(sprintf("  next ranges, not coloured (raise n_top_ranges to include them): %s", paste(sprintf("%s (%d)", head(rest, 12), as.integer(range_n[head(rest, 12)])), collapse = ", ")))
alp_umap$range_class <- ifelse(alp_umap$region %in% top_ranges, alp_umap$region, "Other ranges")
range_k <- tapply(alp_umap$cl, alp_umap$region, function(x) length(unique(x)))                # alpine clusters per range (the size of the nodes of the network)
range_continent <- setNames(cen$continent, cen$region)[top_ranges]                          # one colour family per continent; the ranges of a continent are shades of it, the largest the darkest
range_cols <- setNames(rep(NA_character_, length(top_ranges)), top_ranges)
for (cn in unique(range_continent)) {
  rg <- top_ranges[range_continent == cn]; ends <- CFG$cols_continent_shades[[if (cn %in% names(CFG$cols_continent_shades)) cn else "Other"]]
  range_cols[rg] <- if (length(rg) == 1) colorRampPalette(ends)(3)[2] else colorRampPalette(ends)(length(rg))
}
message(sprintf("UMAP samples: %d alpine samples in %d %s, %d without a range; the %d with the most samples: %s; and %d lowland samples (evenly spaced over latitude / longitude, at most one per %g m)",
                nrow(alp_umap), length(range_n), region_kind, sum(alp_umap$region == "Other"), length(top_ranges),
                paste(sprintf("%s (%d)", top_ranges, range_n[top_ranges]), collapse = ", "), nrow(low_umap), CFG$thin_m))

# Range comparison: in every replicate ONE random alpine sample per cluster (ALL clusters), at most one per thin_m metres, with its range
region_sets <- lapply(seq_len(CFG$n_rep), function(r) {
  set.seed(7000 + r)
  a <- meta[meta$group == "Alpine", ]; a <- a[sample.int(nrow(a)), ]; a <- a[!duplicated(a$cl), ]; a <- a[thin_sites(a$lat, a$lon, CFG$thin_m), ]
  a$region <- region_of(a$lon, a$lat); a
})
reg_n <- table(region_sets[[1]]$region)                                                    # alpine clusters per range
reg_keep <- top_ranges[top_ranges %in% names(reg_n)]                                       # the ranges compared: the top ranges
reg_warn <- reg_keep[reg_n[reg_keep] < CFG$region_warn_clusters]
reg_all <- unique(c(reg_keep, setdiff(names(reg_n)[reg_n >= CFG$net_min_clusters], "Other")))                # the nodes of the network: the coloured ranges and the others with enough clusters
message(sprintf("Range comparison: %d alpine clusters; the top ranges and their clusters: %s", nrow(region_sets[[1]]), if (length(reg_keep)) paste(sprintf("%s (%d)", reg_keep, reg_n[reg_keep]), collapse = ", ") else "none"))
all_ids <- unique(c(unlist(lapply(analyses, function(a) unlist(lapply(a$sel, `[[`, "id")))), alp_umap$id, low_umap$id, unlist(lapply(region_sets, `[[`, "id"))))
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
    adj[[metric]] <- data.frame(metric = metric, group = as.character(P$group), rep_id = P$rep_id, x = P$log10_geo, y = P[[metric]] - fitted(m) + pred_at(m, as.character(P$group), P$log10_geo))
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

## ---- 4b. How many replicates? ----------------------------------------------------------------------------------------------------
# The replicates redraw which sample represents a cluster and which pairs are used; they add no independent information. The number is enough when the
# estimate stops moving: the alpine - lowland level (at the reference distance) and rise per 10-fold distance, Bray-Curtis, from the first k replicates and
# from the odd and the even replicates separately. The two halves should differ by much less than the standard error of the full set.
conv <- NULL
if (CFG$df_dist == 1 && all(c("Alpine", "Lowland") %in% main$groups) && length(unique(main$P$rep_id)) >= 4) {
  P_all <- main$P; reps_all <- sort(unique(P_all$rep_id))
  c0m <- if (is.finite(main$ref_main)) log10(main$ref_main) else median(P_all$log10_geo)
  conv_one <- function(rs) {
    d <- P_all[P_all$rep_id %in% rs, ]
    m <- lm(as.formula(sprintf("bray ~ group * I(log10_geo - %s) + d_depth + d_unm", format(c0m, digits = 10))), data = d)
    V <- dyadic_vcov(model.matrix(m), resid(m), d$ca, d$cb); b <- coef(m); sl <- grep("^groupAlpine:", names(b), value = TRUE)
    data.frame(level = unname(b["groupAlpine"]), level_se = sqrt(V["groupAlpine", "groupAlpine"]), rise = unname(b[sl]), rise_se = sqrt(V[sl, sl]))
  }
  ks <- sort(unique(pmin(c(5, 10, 20, 30, 40, length(reps_all)), length(reps_all)))); ks <- ks[ks >= 2]
  idx_odd <- seq_along(reps_all) %% 2 == 1
  conv <- bind_rows(lapply(ks, function(k) cbind(set = paste("first", k, "replicates"), replicates = k, conv_one(reps_all[seq_len(k)]))),
                    cbind(set = "odd replicates", replicates = sum(idx_odd), conv_one(reps_all[idx_odd])), cbind(set = "even replicates", replicates = sum(!idx_odd), conv_one(reps_all[!idx_odd])))
  write.csv(conv, out_t("replicate_convergence.csv"), row.names = FALSE)
  fr <- conv[grepl("^first", conv$set), ]; full <- fr[which.max(fr$replicates), ]; od <- conv[conv$set == "odd replicates", ]; ev <- conv[conv$set == "even replicates", ]
  d_lv <- abs(od$level - ev$level) / full$level_se; d_ri <- abs(od$rise - ev$rise) / full$rise_se
  message("\nReplicates: alpine - lowland (Bray-Curtis) from the first k replicates, level at the reference distance: ", paste(sprintf("k=%d %s", fr$replicates, fmt(fr$level)), collapse = ", "),
          "; rise per 10-fold distance: ", paste(sprintf("k=%d %s", fr$replicates, fmt(fr$rise)), collapse = ", "))
  message(sprintf("  odd vs even replicates differ by %.0f%% (level) and %.0f%% (rise) of the standard error of the full set (the noise of the full set is about half of that): %s",
                  100 * d_lv, 100 * d_ri, if (max(d_lv, d_ri) < 0.2) "enough replicates" else "more replicates would still move the estimate"))
}

## ---- 5. UMAP ---------------------------------------------------------------------------------------------------------------------------
# The relative abundances of the umap_top_otus most abundant OTUs (as in the full community), Hellinger-transformed (square root), centred and
# reduced to umap_pcs principal components (a partial SVD, RSpectra), then embedded with UMAP (uwot; Euclidean distance on the components; the first
# two components are the starting layout). No distance matrix between all samples is computed, so it takes a few minutes.
for (pkg in c("uwot", "RSpectra")) if (!requireNamespace(pkg, quietly = TRUE)) stop("the UMAP needs the R package ", pkg, " (conda: r-", tolower(pkg), "; add it to envs/r_analyses/environment.yml)")
umap_file <- file.path(out_folder, "umap_result.rds")
umap_ids <- c(alp_umap$id, low_umap$id)
umap_key <- list(n = length(umap_ids), chars = sum(nchar(umap_ids)), top = CFG$umap_top_otus, pcs = CFG$umap_pcs, nn = CFG$umap_neighbors, md = CFG$umap_min_dist,
                 seed = CFG$umap_seed, rare = CFG$rare_depth, test = TEST_MODE, v = 2)
nm <- NULL; pcs_keep <- NULL
if (REUSE_FITS && file.exists(umap_file)) { umap_cached <- readRDS(umap_file); if (identical(umap_cached$key, umap_key)) { nm <- umap_cached$nm; pcs_keep <- umap_cached$pcs; message("Reusing the saved UMAP: ", umap_file) } }
if (is.null(nm)) {
  S <- bind_rows(alp_umap, low_umap) %>% filter(id %in% colnames(X_all))
  X <- X_all[, S$id, drop = FALSE]; tot <- Matrix::colSums(X); ok <- tot >= CFG$rare_depth
  S <- S[ok, ]; X <- X[, ok, drop = FALSE]; tot <- tot[ok]
  Prel <- X %*% Matrix::Diagonal(x = 1 / tot)
  rs <- Matrix::rowSums(Prel); top <- order(rs, decreasing = TRUE)[seq_len(min(CFG$umap_top_otus, nrow(Prel)))]
  Hraw <- sqrt(as.matrix(Matrix::t(Prel[top, , drop = FALSE])))                           # samples x the most abundant OTUs, Hellinger-transformed
  kpc <- min(CFG$umap_pcs, ncol(Hraw) - 1, nrow(Hraw) - 1)
  message(sprintf("UMAP: %d samples (%d alpine, %d lowland) x %d OTUs (%.1f%% of all reads), %d principal components",
                  nrow(Hraw), sum(S$group == "Alpine"), sum(S$group == "Lowland"), ncol(Hraw), 100 * sum(rs[top]) / sum(rs), kpc))
  Hc <- sweep(Hraw, 2, colMeans(Hraw)); rm(Hraw)
  sv <- RSpectra::svds(Hc, k = kpc); pcs <- sv$u %*% diag(sv$d, nrow = kpc); rm(Hc, sv)
  emb <- uwot::umap(pcs, n_neighbors = min(CFG$umap_neighbors, nrow(pcs) - 1), min_dist = CFG$umap_min_dist, metric = "euclidean",
                    init = pcs[, 1:2], init_sdev = 1e-4, n_threads = CFG$umap_threads, seed = CFG$umap_seed, verbose = FALSE)
  nm <- data.frame(id = S$id, group = as.character(S$group), lon = S$lon, lat = S$lat, UMAP1 = emb[, 1], UMAP2 = emb[, 2])
  pcs_keep <- pcs; rownames(pcs_keep) <- S$id
  saveRDS(list(key = umap_key, nm = nm, pcs = pcs_keep), umap_file)
}
nm$group <- factor(nm$group, levels = CFG$groups)
nm$cl <- meta$cl[match(nm$id, meta$id)]                                                    # the cluster of every sample
rc <- setNames(alp_umap$range_class, alp_umap$id); rg <- setNames(alp_umap$region, alp_umap$id)          # the range of every alpine sample (the embedding does not depend on it)
nm$range <- ifelse(nm$group == "Alpine", unname(rg[nm$id]), NA_character_)
nm$range_class <- ifelse(nm$group == "Alpine", unname(rc[nm$id]), "Lowland")
write.csv(nm, out_t("umap_samples.csv"), row.names = FALSE)
write.csv(data.frame(category = names(table(nm$range_class)), samples = as.integer(table(nm$range_class))), out_t("umap_summary.csv"), row.names = FALSE)

## ---- 5b. Region grid: the dissimilarities of ALL pairs of alpine clusters within and between regions -----------------------------------
# Per replicate: the Bray-Curtis dissimilarity of every pair of the alpine samples (one per cluster) of the regions kept, from sparse matrices
region_pairs <- function(S0) {
  S <- S0 %>% filter(id %in% colnames(X_all), region %in% reg_all)
  X <- X_all[, S$id, drop = FALSE]; tot <- Matrix::colSums(X); ok <- tot >= CFG$rare_depth
  S <- S[ok, ]; X <- X[, ok, drop = FALSE]; tot <- tot[ok]
  Prel <- X %*% Matrix::Diagonal(x = 1 / tot)
  ij <- which(upper.tri(matrix(0, nrow(S), nrow(S))), arr.ind = TRUE)
  bray <- as.numeric(Matrix::colSums(abs(Prel[, ij[, 1], drop = FALSE] - Prel[, ij[, 2], drop = FALSE]))) / 2
  data.frame(ra = S$region[ij[, 1]], rb = S$region[ij[, 2]], ca = S$cl[ij[, 1]], cb = S$cl[ij[, 2]], bray = bray,
             geo_km = geosphere::distHaversine(cbind(S$lon[ij[, 1]], S$lat[ij[, 1]]), cbind(S$lon[ij[, 2]], S$lat[ij[, 2]])) / 1000,
             d_depth = abs(S$log10_reads[ij[, 1]] - S$log10_reads[ij[, 2]]), d_unm = abs(S$unmapped[ij[, 1]] - S$unmapped[ij[, 2]]), stringsAsFactors = FALSE)
}
rp <- NULL
if (length(reg_all) >= 2) {
  rp_file <- file.path(out_folder, "region_pairs.rds")
  rp_key <- list(ids = sum(nchar(unlist(lapply(region_sets, `[[`, "id")))), n = CFG$n_rep, keep = reg_all, src = region_kind, rare = CFG$rare_depth, test = TEST_MODE, v = 3)
  if (REUSE_FITS && file.exists(rp_file)) { rp_c <- readRDS(rp_file); if (identical(rp_c$key, rp_key)) { rp <- rp_c$rp; message("Reusing the saved region pairs: ", rp_file) } }
  if (is.null(rp)) {
    message("Region grid: dissimilarities of all pairs of alpine clusters, ", CFG$n_rep, " replicates")
    rp <- bind_rows(lapply(seq_along(region_sets), function(r) cbind(region_pairs(region_sets[[r]]), rep_id = r)))
    saveRDS(list(key = rp_key, rp = rp), rp_file)
  }
} else message("Range comparison: fewer than 2 ranges: no network (panel D), no S7, S8")
cell_stats <- NULL; reg_order <- NULL
if (!is.null(rp)) {
  # every pair once, ordered so that region a comes before region b alphabetically (the grid is symmetric)
  sw <- rp$ra > rp$rb; tmp <- rp$ra[sw]; rp$ra[sw] <- rp$rb[sw]; rp$rb[sw] <- tmp
  # the dissimilarity adjusted to the median differences in sequencing depth and unmapped proportion (and, in the second version, to the median geographic
  # distance of the pairs); the coefficients are estimated WITHIN pairs of ranges (every cell is demeaned), a common slope for all cells
  # ONE value per pair of clusters: the mean over the replicates (the replicates only redraw which sample represents a cluster), so that every pair of clusters, and
  # with it every cluster, counts equally; the adjustment, the medians, the network and the violins all work on these pairs of clusters
  rp$cpkey <- paste(pmin(rp$ca, rp$cb), pmax(rp$ca, rp$cb), sep = "|")
  cp <- rp %>% group_by(ra, rb, cpkey) %>% summarise(ca = first(ca), cb = first(cb), bray = mean(bray), d_depth = mean(d_depth), d_unm = mean(d_unm), geo_km = mean(geo_km), n_rep = n(), .groups = "drop") %>% as.data.frame()
  cp$cell <- paste(cp$ra, cp$rb, sep = " | "); cp$log10_geo <- log10(pmax(cp$geo_km, 0.1))
  adjust_pairs <- function(with_space) {
    dm <- function(x) x - ave(x, cp$cell)
    X <- cbind(d_depth = dm(cp$d_depth), d_unm = dm(cp$d_unm)); if (with_space) X <- cbind(X, log10_geo = dm(cp$log10_geo))
    b <- qr.solve(X, dm(cp$bray))
    y <- cp$bray - b["d_depth"] * (cp$d_depth - median(cp$d_depth)) - b["d_unm"] * (cp$d_unm - median(cp$d_unm))
    if (with_space) y <- y - b["log10_geo"] * (cp$log10_geo - median(cp$log10_geo))
    clip01(y)
  }
  cp$bray_adj <- adjust_pairs(FALSE); cp$bray_adj_space <- adjust_pairs(TRUE)
  cell_stats <- cp %>% group_by(ra, rb) %>% summarise(n_cluster_pairs = n(), n_pairs_pooled = sum(n_rep), median = median(bray), mean = mean(bray), q25 = quantile(bray, 0.25), q75 = quantile(bray, 0.75),
                                                      median_adj = median(bray_adj), mean_adj = mean(bray_adj), median_adj_space = median(bray_adj_space), mean_adj_space = mean(bray_adj_space),
                                                      median_km = median(geo_km), .groups = "drop") %>% as.data.frame()
  cell_stats$clusters_a <- as.integer(range_k[cell_stats$ra]); cell_stats$clusters_b <- as.integer(range_k[cell_stats$rb])
  cell_stats$few_clusters <- cell_stats$ra %in% reg_warn | cell_stats$rb %in% reg_warn                       # a region with fewer than region_warn_clusters clusters
  write.csv(cell_stats, out_t("range_pair_dissimilarity.csv"), row.names = FALSE)
  # rows and columns in the order of a hierarchical clustering of the median dissimilarities
  # (the coloured ranges only: they are the violin grid of S8)
  cell_stats_top <- cell_stats[cell_stats$ra %in% reg_keep & cell_stats$rb %in% reg_keep, ]
  Mmed <- matrix(NA_real_, length(reg_keep), length(reg_keep), dimnames = list(reg_keep, reg_keep))
  for (i in seq_len(nrow(cell_stats_top))) { Mmed[cell_stats_top$ra[i], cell_stats_top$rb[i]] <- cell_stats_top$median[i]; Mmed[cell_stats_top$rb[i], cell_stats_top$ra[i]] <- cell_stats_top$median[i] }
  Mfill <- Mmed; Mfill[is.na(Mfill)] <- max(Mmed, na.rm = TRUE)
  reg_order <- reg_keep[if (length(reg_keep) > 2) hclust(as.dist(Mfill), method = "average")$order else seq_along(reg_keep)]
  message(sprintf("Range comparison: %s pairs of %d ranges; median dissimilarity from %.3f to %.3f (adjusted for depth / unmapped: %.3f to %.3f; also for distance: %.3f to %.3f)", format(nrow(rp), big.mark = ","), length(reg_all),
                  min(cell_stats$median), max(cell_stats$median), min(cell_stats$median_adj), max(cell_stats$median_adj), min(cell_stats$median_adj_space), max(cell_stats$median_adj_space)))
}

## ---- 5c. UMAP: are alpine samples more similar to each other than lowland samples, and how do the spreads compare? -----------------------------------
# BEGIN UMAPSTATS
# ---- Are alpine samples more similar to each other than lowland samples, and how do the spreads compare? ---------------------------------------------------
# Random pairs of samples of one group from DIFFERENT clusters (in every replicate ONE random sample per cluster, as in the main analysis), with their distances in
# the given spaces and their geographic distance. d: data.frame with group, cl, lon, lat (rows aligned with the rows of the matrices in `spaces`).
umap_pairs <- function(d, spaces, n_rep, n_pairs, seed) {
  set.seed(seed)
  ix_by_group <- split(seq_len(nrow(d)), as.character(d$group))
  bind_rows(lapply(seq_len(n_rep), function(r) bind_rows(lapply(names(ix_by_group), function(g) {
    ix <- ix_by_group[[g]]; ix <- ix[sample.int(length(ix))]; ix <- ix[!duplicated(d$cl[ix])]                 # one random sample per cluster
    if (length(ix) < 3) return(NULL)
    np <- min(n_pairs, choose(length(ix), 2)); a <- sample(ix, np, TRUE); b <- sample(ix, np, TRUE); k <- a != b; a <- a[k]; b <- b[k]
    out <- data.frame(rep = r, group = g, ca = d$cl[a], cb = d$cl[b], stringsAsFactors = FALSE,
                      geo_km = geosphere::distHaversine(cbind(d$lon[a], d$lat[a]), cbind(d$lon[b], d$lat[b])) / 1000)
    for (sp in names(spaces)) out[[sp]] <- sqrt(rowSums((spaces[[sp]][a, , drop = FALSE] - spaces[[sp]][b, , drop = FALSE])^2))
    out
  }))))
}
# The difference alpine - lowland in the mean distance between samples of different clusters (< 0: alpine samples are more similar), crude and at the same
# geographic distance; cluster-robust standard errors for pair data (the clusters are the independent units)
umap_sim_test <- function(P, resp) {
  P$group <- factor(P$group, levels = c("Lowland", "Alpine")); P$log10_geo <- log10(pmax(P$geo_km, 0.1))
  mn <- tapply(P[[resp]], P$group, mean)
  bind_rows(lapply(c(crude = FALSE, adjusted = TRUE), function(adj) {
    m <- lm(as.formula(paste(resp, "~ group", if (adj) "+ log10_geo" else "")), data = P)
    V <- dyadic_vcov(model.matrix(m), resid(m), P$ca, P$cb); dfree <- max(attr(V, "G") - 1, 1)
    est <- unname(coef(m)["groupAlpine"]); se <- sqrt(V["groupAlpine", "groupAlpine"])
    data.frame(mean_lowland = unname(mn["Lowland"]), mean_alpine = unname(mn["Alpine"]), estimate = est, se = se, lo = est - qt(0.975, dfree) * se, hi = est + qt(0.975, dfree) * se,
               p = 2 * pt(-abs(est / se), dfree), clusters = attr(V, "G"))
  }), .id = "model")
}
# Per cluster: number of samples, centre, sum of squares and mean distance of its samples to the centre
cluster_summary <- function(X, cl) {
  X <- as.matrix(X); cl <- as.character(cl); n <- tapply(rep(1, length(cl)), cl, sum)
  C <- rowsum(X, cl); n <- as.numeric(n[rownames(C)]); C <- C / n
  D2 <- rowSums((X - C[cl, , drop = FALSE])^2)
  list(n = n, C = C, ssw = rowsum(D2, cl)[rownames(C), 1], md_w = rowsum(sqrt(D2), cl)[rownames(C), 1] / n)
}
# Spread of one group from its cluster summary (idx = the clusters used: a bootstrap sample or all). The clusters count equally, except in the variance share.
spread_metrics <- function(cs, idx = seq_along(cs$n)) {
  C <- cs$C[idx, , drop = FALSE]; n <- cs$n[idx]; ssw <- cs$ssw[idx]
  d_cent <- mean(sqrt(rowSums(sweep(C, 2, colMeans(C))^2)))                                     # distance of the cluster centres to their centre
  var_cent <- sum(apply(C, 2, var))                                                             # total variance of the cluster centres
  area <- if (ncol(C) == 2 && nrow(C) > 2) pi * qchisq(0.95, 2) * sqrt(max(det(cov(C)), 0)) else NA_real_    # 95% ellipse of the cluster centres (2D only)
  ok <- n >= 2; within <- if (any(ok)) mean(cs$md_w[idx][ok]) else NA_real_                     # distance of the samples to their cluster centre (clusters with 2+ samples)
  gw <- colSums(C * n) / sum(n); ssb <- sum(n * rowSums(sweep(C, 2, gw)^2))                     # between-cluster share of the total variance of the samples
  c(centroid_dispersion = d_cent, centroid_variance = var_cent, ellipse_area = area, within_cluster_distance = within, between_cluster_share = ssb / (ssb + sum(ssw)))
}
# Alpine vs lowland: ratio of the spread measures (difference for the variance share), 95% interval and p from a bootstrap that resamples whole clusters
spread_boot <- function(csA, csL, B) {
  fA <- spread_metrics(csA); fL <- spread_metrics(csL); is_ratio <- names(fA) != "between_cluster_share"
  comb <- function(a, l) ifelse(is_ratio, a / l, a - l)
  bs <- t(replicate(B, comb(spread_metrics(csA, sample.int(length(csA$n), replace = TRUE)), spread_metrics(csL, sample.int(length(csL$n), replace = TRUE)))))
  null <- ifelse(is_ratio, 1, 0)
  data.frame(metric = names(fA), kind = ifelse(is_ratio, "ratio alpine / lowland", "difference alpine - lowland"), alpine = unname(fA), lowland = unname(fL), estimate = unname(comb(fA, fL)),
             lo = apply(bs, 2, function(b) if (any(is.finite(b))) unname(quantile(b[is.finite(b)], 0.025)) else NA_real_),
             hi = apply(bs, 2, function(b) if (any(is.finite(b))) unname(quantile(b[is.finite(b)], 0.975)) else NA_real_),
             p = vapply(seq_along(null), function(j) { b <- bs[, j]; b <- b[is.finite(b)]; if (!length(b)) NA_real_ else min(1, 2 * (1 + min(sum(b <= null[j]), sum(b >= null[j]))) / (length(b) + 1)) }, 1), row.names = NULL)   # a bootstrap p cannot be below ~2 / resamples
}
# END UMAPSTATS

umap_space_lab <- c(umap2d = "UMAP (2D)", pcs = sprintf("Hellinger PCs (%dD)", ncol(pcs_keep)))
spaces <- list(umap2d = as.matrix(nm[, c("UMAP1", "UMAP2")]), pcs = pcs_keep[nm$id, , drop = FALSE])
d_grp <- data.frame(group = as.character(nm$group), cl = nm$cl, lon = nm$lon, lat = nm$lat, stringsAsFactors = FALSE)
UP <- umap_pairs(d_grp, spaces, CFG$umap_test_rep, CFG$umap_test_pairs, 600)
sim_tab <- bind_rows(lapply(names(spaces), function(sp) cbind(space = sp, umap_sim_test(UP, sp))))
write.csv(sim_tab, out_t("umap_similarity_test.csv"), row.names = FALSE)
cs_list <- lapply(spaces, function(X) list(Alpine = cluster_summary(X[nm$group == "Alpine", , drop = FALSE], nm$cl[nm$group == "Alpine"]),
                                           Lowland = cluster_summary(X[nm$group == "Lowland", , drop = FALSE], nm$cl[nm$group == "Lowland"])))
set.seed(601)
spread_tab <- bind_rows(lapply(names(spaces), function(sp) cbind(space = sp, spread_boot(cs_list[[sp]]$Alpine, cs_list[[sp]]$Lowland, CFG$umap_boot))))
write.csv(spread_tab, out_t("umap_spread.csv"), row.names = FALSE)
message(sprintf("\nUMAP test: %d alpine clusters and %d lowland clusters (the independent units); %s pairs of samples of different clusters", length(cs_list$umap2d$Alpine$n), length(cs_list$umap2d$Lowland$n), format(nrow(UP), big.mark = ",")))
for (sp in names(spaces)) {
  d <- sim_tab[sim_tab$space == sp, ]
  message(sprintf("  %s: mean distance between samples of different clusters: alpine %.3g, lowland %.3g; alpine - lowland %s [%s, %s], p %s; at the same geographic distance %s [%s, %s], p %s",
                  umap_space_lab[[sp]], d$mean_alpine[1], d$mean_lowland[1], fmt(d$estimate[1]), fmt(d$lo[1]), fmt(d$hi[1]), pfmt(d$p[1]), fmt(d$estimate[2]), fmt(d$lo[2]), fmt(d$hi[2]), pfmt(d$p[2])))
  d <- spread_tab[spread_tab$space == sp & !is.na(spread_tab$estimate), ]
  for (i in seq_len(nrow(d))) message(sprintf("     %-26s %s %.3g [%.3g, %.3g], p %s  (alpine %.3g, lowland %.3g)", d$metric[i], if (grepl("^ratio", d$kind[i])) "ratio" else "difference", d$estimate[i], d$lo[i], d$hi[i], pfmt(d$p[i]), d$alpine[i], d$lowland[i]))
}
# the note on panel C (the 2D UMAP) and the sentence of the caption
u1 <- sim_tab[sim_tab$space == "umap2d", ]; u2 <- spread_tab[spread_tab$space == "umap2d", ]; ud <- u2[u2$metric == "centroid_dispersion", ]; uw <- u2[u2$metric == "within_cluster_distance", ]
umap_note <- paste0("Samples of different clusters, mean distance: alpine ", signif(u1$mean_alpine[1], 3), ", lowland ", signif(u1$mean_lowland[1], 3), "\n",
                    "alpine - lowland ", fmt(u1$estimate[1]), " [", fmt(u1$lo[1]), ", ", fmt(u1$hi[1]), "], p ", pfmt(u1$p[1]), ";  same geographic distance: ", fmt(u1$estimate[2]), " [", fmt(u1$lo[2]), ", ", fmt(u1$hi[2]), "], p ", pfmt(u1$p[2]), "\n",
                    "Spread of the cluster centres, alpine / lowland: ", sprintf("%.2f [%.2f, %.2f]", ud$estimate, ud$lo, ud$hi), ";  of the samples within a cluster: ", sprintf("%.2f [%.2f, %.2f]", uw$estimate, uw$lo, uw$hi))
umap_cap <- sprintf("Text in C: the mean distance in the UMAP between two samples of different clusters (one random sample per cluster, %d replicates, %s pairs), the difference alpine - lowland (negative: alpine samples closer together) with its 95%% interval, estimated with standard errors that are cluster-robust for pair data and also at the same geographic distance, and the ratio alpine / lowland of the spread of the cluster centres (distance to their centre) and of the samples within their clusters, with 95%% intervals from a bootstrap of the clusters. The UMAP axes are not metric: the same tests in the space the UMAP was made from are in Supplementary Figure 9 and tables/umap_similarity_test.csv, umap_spread.csv.",
                    CFG$umap_test_rep, format(nrow(UP), big.mark = ","))

## ---- 6. Figures -----------------------------------------------------------------------------------------------------------------
x_scale <- scale_x_continuous(breaks = log10(c(1, 10, 100, 1000, 10000)), labels = c("1", "10", "100", "1,000", "10,000"))
col_scales <- list(scale_colour_manual(values = CFG$cols, name = NULL, breaks = CFG$groups, labels = c(Lowland = "Lowland", Alpine = "Alpine")),
                   scale_fill_manual(values = CFG$cols, guide = "none"),
                   guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1, linewidth = 1.2, linetype = "solid"))))
set.seed(5)

# the differences between groups of one analysis and metric as two lines of text for a panel (A = alpine, L = lowland):
# the estimates and the stars of the interaction model's terms (level at the reference distance; rise per 10-fold distance)
annot_text <- function(rs, metric) {
  d <- rs$tests[rs$tests$metric == metric, ]
  if (is.null(d) || !nrow(d)) return("")
  line <- function(dd) paste(sprintf("%s-%s %s %s", short[dd$g1], short[dd$g2], fmt(dd$estimate), stars(dd$p)), collapse = "   ")
  paste(c(paste0("Level at ", fmt_km(rs$ref_main), " km:  ", line(d[d$quantity == "level", ])), paste0("Rise per 10-fold distance:  ", line(d[d$quantity == "slope", ]))), collapse = "\n")
}

decay_panel <- function(rs, metric, tag, title = NULL, show_x = TRUE) {
  cv <- rs$curves_sum[rs$curves_sum$metric == metric, ]; cv$group <- factor(cv$group, levels = CFG$groups)
  ad <- rs$adj[rs$adj$metric == metric, ] %>% group_by(group) %>% mutate(bin = ntile(x, CFG$n_bins)) %>% ungroup()               # equal numbers of pairs per bin, in every group
  bn <- ad %>% group_by(group, bin) %>% summarise(x = mean(x), y = mean(y), .groups = "drop") %>% as.data.frame(); bn$group <- factor(bn$group, levels = CFG$groups)
  br <- if (isTRUE(CFG$rep_dots)) ad %>% group_by(group, bin, rep_id) %>% summarise(x = mean(x), y = mean(y), .groups = "drop") %>% as.data.frame() else NULL     # the same bins in every replicate
  if (!is.null(br)) br$group <- factor(br$group, levels = CFG$groups)
  pts <- rs$P[sample.int(nrow(rs$P), min(nrow(rs$P), CFG$points_max)), ]; pts$group <- factor(pts$group, levels = CFG$groups)
  ylo <- max(0, floor((min(bn$y, cv$lo, if (!is.null(br)) quantile(br$y, 0.01), na.rm = TRUE) - 0.06) * 20) / 20)               # zoom: pairs below the window are not drawn (they are in the models)
  ggplot() +
    geom_point(data = pts, aes(log10_geo, .data[[metric]], colour = group), size = 0.3, alpha = CFG$point_alpha, shape = 16) +
    geom_ribbon(data = cv, aes(x, ymin = lo, ymax = hi, fill = group), alpha = 0.2) +
    geom_line(data = cv, aes(x, fit, colour = group), linewidth = 1.1) +
    (if (!is.null(br)) geom_point(data = br, aes(x, y, colour = group), shape = 16, size = 0.8, alpha = 0.4)) +
    geom_point(data = bn, aes(x, y, colour = group), shape = 21, fill = "white", size = 2.4, stroke = 1.1) +
    annotate("text", x = -Inf, y = -Inf, label = annot_text(rs, metric), hjust = -0.03, vjust = -0.5, size = 2.9, lineheight = 1.1) +
    x_scale + scale_y_continuous(breaks = function(l) { b <- scales::extended_breaks()(l); b[b <= 1 + 1e-9] }) + coord_cartesian(ylim = c(ylo, 1.005)) + col_scales +
    labs(tag = tag, title = title, x = if (show_x) "Geographic distance between samples (km)" else NULL, y = paste(metrics[[metric]], "dissimilarity")) + theme_fig +
    theme(plot.title = element_text(size = 10, face = "bold"))
}
pA <- decay_panel(main, "bray", "A"); pB <- decay_panel(main, "jaccard", "B")

# C: the UMAP as small multiples (umap_facet): the first panel all alpine samples (blue, with an ellipse) against the lowland samples (grey, with an ellipse), then one panel per
# continent in which the ranges of that continent are coloured (one colour family per continent) and everything else is a pale grey background; the same coordinates in all
# panels. Ellipses: the central half of the samples of a range (robust multivariate t), only for ranges with at least 10 samples. The legend, grouped by continent, and the
# note on the test are a panel of their own, to the right.
pt_size  <- function(n) if (n < 1500) 1.5 else 0.9
pt_alpha <- function(n) min(0.8, max(0.3, 4000 / n))                                       # visible with few points, not a blob with many
cat_levels <- c(top_ranges, if (any(nm$range_class == "Other ranges")) "Other ranges", "Lowland")
n_cat <- table(factor(nm$range_class, levels = cat_levels))
cols_cat <- c(range_cols, "Other ranges" = CFG$col_other_ranges, Lowland = CFG$col_lowland)[cat_levels]
k_cat <- tapply(nm$cl, factor(nm$range_class, levels = cat_levels), function(x) length(unique(x)))
lab_cat <- setNames(sprintf("%s (%s samples, %d clusters)", cat_levels, format(as.integer(n_cat), big.mark = ",", trim = TRUE), as.integer(k_cat)), cat_levels)
nm$range_class <- factor(nm$range_class, levels = cat_levels)
cont_n <- tapply(as.integer(n_cat[top_ranges]), range_continent, sum); cont_order <- names(sort(cont_n, decreasing = TRUE))
rc_chr <- as.character(nm$range_class)
cont_s <- ifelse(nm$group == "Alpine", ifelse(!is.na(range_continent[rc_chr]), range_continent[rc_chr], continent_of(nm$lat, nm$lon)), "Lowland")   # the continent of every sample
extra <- setdiff(names(which(table(cont_s[nm$group == "Alpine"]) >= 10)), cont_order)                    # continents that have only samples of the uncoloured ranges
overview <- if (isTRUE(CFG$umap_facet)) "Alpine vs lowland" else "All samples"
facet_levels <- if (isTRUE(CFG$umap_facet)) c(overview, cont_order, extra) else overview
alp_col <- unname(CFG$cols["Alpine"])
fac <- bind_rows(lapply(facet_levels, function(f) {
  fg <- if (f == overview) nm$group == "Alpine" else nm$group == "Alpine" & cont_s == f
  d <- data.frame(UMAP1 = nm$UMAP1, UMAP2 = nm$UMAP2, rc = rc_chr, group = as.character(nm$group), fg = fg, facet = f, stringsAsFactors = FALSE)
  d$col <- NA_character_; d$grp <- NA_character_
  by_range <- !(f == overview && isTRUE(CFG$umap_facet))                                         # the overview panel has one colour for all alpine samples
  d$col[fg] <- if (by_range) unname(cols_cat[d$rc[fg]]) else alp_col
  d$grp[fg] <- if (by_range) d$rc[fg] else "Alpine"
  d
}))
fac$facet <- factor(fac$facet, levels = facet_levels)
ell_ranges <- top_ranges[top_ranges %in% names(which(table(nm$range_class) >= 10))]               # an ellipse needs enough samples
ell_d <- fac[fac$fg & (fac$grp %in% ell_ranges | (fac$grp == "Alpine" & fac$facet == overview)), ]
if (isTRUE(CFG$umap_facet)) { low_e <- fac[fac$facet == overview & fac$group == "Lowland", ]; low_e$grp <- "Lowland"; low_e$col <- "#6B6B6B"; ell_d <- bind_rows(ell_d, low_e) }
fg_d <- fac[fac$fg, ]; fg_d <- fg_d[order(match(fg_d$rc, c("Other ranges", top_ranges))), ]                    # the large ranges below, the small ones on top
pC <- ggplot() +
  geom_point(data = fac[!fac$fg, ], aes(UMAP1, UMAP2), colour = "#D2D2D2", size = pt_size(nrow(nm)) * 0.8, alpha = 0.35, shape = 16) +
  geom_point(data = fg_d, aes(UMAP1, UMAP2, colour = col), size = pt_size(nrow(nm)), alpha = pt_alpha(nrow(nm)), shape = 16) +
  stat_ellipse(data = ell_d, aes(UMAP1, UMAP2, group = grp, colour = col), level = CFG$umap_ellipse_level, type = CFG$umap_ellipse_type, linewidth = 0.9, show.legend = FALSE) +
  scale_colour_identity() +
  facet_wrap(~ facet, ncol = if (length(facet_levels) > 4) 3 else 2) +
  labs(tag = "C", x = "UMAP 1", y = "UMAP 2") + theme_fig + theme(legend.position = "none", panel.spacing = unit(0.6, "lines"))
# the legend of C: grouped by continent (a ggplot legend has no group headings), with the note on the test, drawn as a panel of its own
rows_l <- list(data.frame(label = "Alpine samples by mountain range", kind = "title", col = NA_character_, stringsAsFactors = FALSE))
for (cn in cont_order) {
  rows_l[[length(rows_l) + 1]] <- data.frame(label = cn, kind = "head", col = NA_character_, stringsAsFactors = FALSE)
  rg <- top_ranges[range_continent == cn]
  rows_l[[length(rows_l) + 1]] <- data.frame(label = lab_cat[rg], kind = "item", col = unname(range_cols[rg]), stringsAsFactors = FALSE)
}
rest_l <- c(if ("Other ranges" %in% cat_levels) "Other ranges", "Lowland")
rows_l[[length(rows_l) + 1]] <- data.frame(label = "Not coloured", kind = "head", col = NA_character_, stringsAsFactors = FALSE)
rows_l[[length(rows_l) + 1]] <- data.frame(label = lab_cat[rest_l], kind = "item", col = unname(cols_cat[rest_l]), stringsAsFactors = FALSE)
note_l <- unlist(lapply(strsplit(umap_note, "\n")[[1]], function(l) strwrap(l, width = 54)))
rows_l[[length(rows_l) + 1]] <- data.frame(label = "Test in the UMAP (clusters as units)", kind = "head", col = NA_character_, stringsAsFactors = FALSE)
rows_l[[length(rows_l) + 1]] <- data.frame(label = note_l, kind = "note", col = NA_character_, stringsAsFactors = FALSE)
LG <- bind_rows(rows_l); LG$y <- rev(seq_len(nrow(LG))); txt_size <- min(4.1, 4.1 * 30 / nrow(LG))
pC_leg <- ggplot(LG) +
  geom_text(data = LG[LG$kind == "title", ], aes(x = 0, y = y, label = label), hjust = 0, fontface = "bold", size = txt_size * 1.05) +
  geom_text(data = LG[LG$kind == "head", ], aes(x = 0, y = y, label = label), hjust = 0, fontface = "bold", size = txt_size) +
  geom_point(data = LG[LG$kind == "item", ], aes(x = 0.05, y = y, colour = col), size = txt_size * 1.15) +
  geom_text(data = LG[LG$kind == "item", ], aes(x = 0.13, y = y, label = label), hjust = 0, size = txt_size) +
  geom_text(data = LG[LG$kind == "note", ], aes(x = 0, y = y, label = label), hjust = 0, size = txt_size * 0.8) +
  scale_colour_identity() + scale_x_continuous(limits = c(0, 1.3), expand = c(0, 0)) + scale_y_continuous(limits = c(0.4, nrow(LG) + 0.6), expand = c(0, 0)) + theme_void()

# S8: the grid of violins (lower triangle, diagonal = within a region): one vertical scale for all cells, fill = median dissimilarity (blue = similar,
# red = dissimilar), rows and columns ordered by the clustering above. Drawn as polygons on a common coordinate system: cell (row i, column j) has its centre
# at x = j, y = n - i + 1. The dissimilarities pile up near 1, so the vertical scale is stretched there (region_axis = "log1m": -log10(1 - d + eps), the same
# in every cell); the dissimilarity values are the ticks at the left of every row, with a light guide line in every cell.
pG <- NULL
if (!is.null(rp) && length(reg_keep) >= 2) {
  nR <- length(reg_order); lo <- floor(quantile(cp$bray, 0.005) * 20) / 20; lo <- min(lo, floor(min(cell_stats_top$q25) * 20) / 20 - 0.05)
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
  polys <- list(); frames <- list(); meds <- list(); labs_n <- list(); guides_l <- list(); dots <- list()
  for (i in seq_len(nR)) for (j in seq_len(i)) {
    a <- reg_order[j]; b <- reg_order[i]; ra <- min(a, b); rb <- max(a, b)                      # (a < b alphabetically, as stored)
    v <- cp$bray[cp$ra == ra & cp$rb == rb]; cx <- j; cy <- nR - i + 1                         # one value per pair of clusters
    frames[[length(frames) + 1]] <- data.frame(xmin = cx - 0.47, xmax = cx + 0.47, ymin = cy - 0.47, ymax = cy + 0.47)
    guides_l[[length(guides_l) + 1]] <- data.frame(x = cx - 0.47, xend = cx + 0.47, y = ymap(tr(tick_vals), cy))
    if (!length(v)) next
    mv <- median(v)
    meds[[length(meds) + 1]] <- data.frame(x = cx - 0.38, xend = cx + 0.38, y = ymap(tr(mv), cy), yend = ymap(tr(mv), cy))
    labs_n[[length(labs_n) + 1]] <- data.frame(x = cx - 0.44, y = cy - 0.43, label = format(length(v), big.mark = ","))
    if (length(v) < 5) {                                                                        # too few pairs of clusters for a violin: the values themselves
      dots[[length(dots) + 1]] <- data.frame(x = cx + seq(-0.2, 0.2, length.out = max(length(v), 2))[seq_along(v)], y = ymap(tr(v), cy), median = mv); next
    }
    dn <- density(tr(v), from = u_lo, to = u_hi, n = 100, adjust = 1.3)                         # the density on the stretched scale
    hw <- 0.45 * dn$y / max(dn$y)                                                               # equal maximum width: the shape, not the number of pairs
    polys[[length(polys) + 1]] <- data.frame(cell = paste(i, j), x = c(cx - hw, rev(cx + hw)), y = c(ymap(dn$x, cy), rev(ymap(dn$x, cy))), median = mv)
  }
  polys <- bind_rows(polys); frames <- bind_rows(frames); meds <- bind_rows(meds); labs_n <- bind_rows(labs_n); guides_l <- bind_rows(guides_l); dots <- bind_rows(dots)
  ticks <- bind_rows(lapply(seq_len(nR), function(i) data.frame(y = ymap(tr(tick_vals), nR - i + 1), label = tick_lab)))     # the values at the left of every row
  reg_lab <- sprintf("%s (%d)%s", reg_order, as.integer(reg_n[reg_order]), ifelse(reg_order %in% reg_warn, " \u2020", ""))
  note <- paste0("One value per pair of clusters\nBar: median\nSmall number: pairs of clusters\nBrackets: alpine clusters of the region",
                 if (any(reg_order %in% reg_warn)) sprintf("\n\u2020 fewer than %d clusters: interpret with care", CFG$region_warn_clusters) else "")      # the rest is in the caption
  pG <- ggplot() +
    geom_rect(data = frames, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "white", colour = "grey80", linewidth = 0.3) +
    geom_segment(data = guides_l, aes(x = x, xend = xend, y = y, yend = y), colour = "grey88", linewidth = 0.25, linetype = "dotted") +
    geom_polygon(data = polys, aes(x, y, group = cell, fill = median), colour = "grey25", linewidth = 0.25) +
    (if (nrow(dots)) geom_point(data = dots, aes(x, y, fill = median), shape = 21, colour = "grey15", size = 2.2, stroke = 0.5)) +
    geom_segment(data = meds, aes(x = x, xend = xend, y = y, yend = yend), colour = "black", linewidth = 0.5) +
    geom_segment(data = ticks, aes(x = 0.49, xend = 0.53, y = y, yend = y), colour = "grey30", linewidth = 0.3) +
    geom_text(data = ticks, aes(x = 0.47, y = y, label = label), hjust = 1, size = 2.5, colour = "grey30") +
    geom_text(data = labs_n, aes(x, y, label = label), hjust = 0, size = 2.2, colour = "grey35") +
    annotate("text", x = nR + 0.45, y = nR + 0.45, hjust = 1, vjust = 1, size = 2.6, lineheight = 1.1, colour = "grey25", label = note) +
    scale_fill_gradientn(colours = c("#2166AC", "#F7F7F7", "#B2182B"), name = "Median Bray-Curtis dissimilarity", limits = range(cell_stats_top$median),
                         values = scales::rescale(c(min(cell_stats_top$median), mean(range(cell_stats_top$median)), max(cell_stats_top$median)), from = range(cell_stats_top$median))) +
    scale_x_continuous(breaks = seq_len(nR), labels = reg_lab, limits = c(0, nR + 0.5), expand = c(0, 0)) +
    scale_y_continuous(breaks = nR - seq_len(nR) + 1, labels = reg_lab, limits = c(0.5, nR + 0.5), expand = c(0, 0)) +
    coord_fixed() + labs(tag = "", x = NULL, y = NULL) + theme_fig +
    guides(fill = guide_colourbar(title.position = "top", barwidth = unit(3, "in"), barheight = unit(0.13, "in"))) +
    theme(legend.position = "top", legend.title = element_text(size = 9), legend.text = element_text(size = 8), panel.grid = element_blank(), panel.border = element_blank(),
          axis.text.x = element_text(angle = 40, hjust = 1, size = 9), axis.text.y = element_text(size = 9), axis.ticks = element_blank())
}

# D: the network of the mountain ranges (see the header). Nodes placed by classical MDS of the median dissimilarities between the ranges; edges = for every range its
# net_k most similar ranges + the minimum spanning tree; colour of an edge = the median dissimilarity.
mst_edges <- function(Dm) {                                                                # minimum spanning tree (Prim) of a full dissimilarity matrix: the pairs of indices
  n <- nrow(Dm); inT <- 1L; ed <- list(); best <- Dm[1, ]; from <- rep(1L, n)
  while (length(inT) < n) {
    cand <- setdiff(seq_len(n), inT); j <- cand[which.min(best[cand])]
    ed[[length(ed) + 1]] <- c(from[j], j); inT <- c(inT, j)
    upd <- Dm[j, ] < best; upd[inT] <- FALSE; best[upd] <- Dm[j, upd]; from[upd] <- j
  }
  do.call(rbind, ed)
}
net_nodes <- data.frame(region = reg_all, stringsAsFactors = FALSE)
net_nodes$n <- as.integer(range_k[net_nodes$region]); miss <- is.na(net_nodes$n); net_nodes$n[miss] <- as.integer(reg_n[net_nodes$region[miss]]); net_nodes$n[is.na(net_nodes$n)] <- 1L
net_nodes$class <- ifelse(net_nodes$region %in% top_ranges, net_nodes$region, "Other ranges")
net_panel_in <- c(7.4, 5.8)                                                                # the panel of the network in inches (width, height): the scale of the node sizes
# Moves nodes apart where they overlap or are crowded, and pulls them gently back to where the MDS put them (so the layout keeps its structure): the radius of a node is
# in mm (the dot, plus room for the name of the named ones)
declump <- function(xy, n, named, panel_in = net_panel_in, iters = 400, pad = 1.5, pull = 0.04) {
  K <- nrow(xy); if (K < 2) return(xy)
  R <- 5.5 * sqrt(n / max(n)) + ifelse(named, 3, 0)
  cx <- mean(range(xy[, 1])); cy <- mean(range(xy[, 2]))
  sx <- 0.82 * panel_in[1] * 25.4 / max(diff(range(xy[, 1])), 1e-9); sy <- 0.82 * panel_in[2] * 25.4 / max(diff(range(xy[, 2])), 1e-9)
  P0 <- cbind((xy[, 1] - cx) * sx, (xy[, 2] - cy) * sy); P <- P0
  for (it in seq_len(iters)) {
    for (i in seq_len(K - 1)) for (j in (i + 1):K) {
      d <- P[j, ] - P[i, ]; dist <- sqrt(sum(d^2)); if (dist < 1e-6) { d <- c(cos(i * 2.4), sin(i * 2.4)); dist <- 1 }
      need <- R[i] + R[j] + pad; soft <- 1.8 * need
      if (dist < soft) { push <- if (dist < need) (need - dist) / 2 + 0.15 * (soft - need) else 0.15 * (soft - dist); P[i, ] <- P[i, ] - push * d / dist; P[j, ] <- P[j, ] + push * d / dist }
    }
    P <- P + pull * (P0 - P)
  }
  cbind(P[, 1] / sx + cx, P[, 2] / sy + cy)
}
net_layout <- function(value_col, ref = NULL) {
  cs <- cell_stats[cell_stats$ra != cell_stats$rb & cell_stats$ra %in% reg_all & cell_stats$rb %in% reg_all, ]
  Dm <- matrix(NA_real_, length(reg_all), length(reg_all), dimnames = list(reg_all, reg_all))
  for (i in seq_len(nrow(cs))) { Dm[cs$ra[i], cs$rb[i]] <- cs[[value_col]][i]; Dm[cs$rb[i], cs$ra[i]] <- cs[[value_col]][i] }
  Df <- Dm; Df[is.na(Df)] <- max(Dm, na.rm = TRUE); diag(Df) <- 0                                  # a pair of ranges without a pair of clusters: as far as the farthest
  Dst <- pmax(Df - CFG$net_stretch * min(Df[upper.tri(Df)]), 0); diag(Dst) <- 0                     # the differences between the dissimilarities, not their common level, place the nodes
  mds <- cmdscale(as.dist(Dst), k = 2, eig = TRUE); xy <- mds$points; rownames(xy) <- reg_all
  if (!is.null(ref)) xy <- vegan::procrustes(ref, xy, scale = TRUE, symmetric = FALSE)$Yrot            # rotated, reflected and scaled onto the reference layout
  xy_mds <- xy
  xy <- declump(xy, net_nodes$n, reg_all %in% top_ranges)
  K <- length(reg_all); ed <- list()
  for (i in seq_len(K)) for (j in setdiff(order(Df[i, ]), i)[seq_len(min(CFG$net_k, K - 1))]) ed[[length(ed) + 1]] <- c(min(i, j), max(i, j))
  tr <- mst_edges(Df); for (r in seq_len(nrow(tr))) ed[[length(ed) + 1]] <- c(min(tr[r, ]), max(tr[r, ]))
  ed <- unique(do.call(rbind, ed))
  edges <- data.frame(x = xy[ed[, 1], 1], y = xy[ed[, 1], 2], xend = xy[ed[, 2], 1], yend = xy[ed[, 2], 2], value = Df[cbind(ed[, 1], ed[, 2])])
  list(nodes = cbind(net_nodes, x = xy[, 1], y = xy[, 2]), edges = edges, xy = xy, xy_mds = xy_mds, share = mds$GOF[1], n_edges = nrow(edges))
}
# Names without overlaps: for every node (large ones first) the first free position among eight directions and growing distances (a free position does not overlap a name
# already placed, any node or the edge of the panel); estimated text size in mm; ux, uy = plot units per mm of the panel. Returns the anchor and the alignment of every name.
place_labels <- function(pts, names_to_place, bx, by, ux, uy) {
  X <- pts$x / ux; Y <- pts$y / uy; R <- 5.5 * sqrt(pts$n / max(pts$n))                      # nodes in mm (max_size 11 = diameter)
  box <- c(bx[1] / ux, bx[2] / ux, by[1] / uy, by[2] / uy)
  dirs <- list(c(0, 1, 0.5, 0), c(0, -1, 0.5, 1), c(1, 0, 0, 0.5), c(-1, 0, 1, 0.5), c(0.7, 0.7, 0, 0), c(-0.7, 0.7, 1, 0), c(0.7, -0.7, 0, 1), c(-0.7, -0.7, 1, 1))
  placed <- matrix(numeric(0), 0, 4); out <- list()
  for (k in order(-pts$n)) if (pts$region[k] %in% names_to_place) {
    w <- 2.05 * nchar(pts$region[k]); h <- 3.9; best <- NULL
    for (d in c(0.8, 2.5, 4.5, 7, 10, 14, 19)) for (dr in dirs) {
      if (!is.null(best)) break
      ax <- X[k] + dr[1] * (R[k] + d); ay <- Y[k] + dr[2] * (R[k] + d)
      x0 <- ax - dr[3] * w; x1 <- x0 + w; y0 <- ay - dr[4] * h; y1 <- y0 + h
      ok <- x0 >= box[1] && x1 <= box[2] && y0 >= box[3] && y1 <= box[4]
      if (ok && nrow(placed)) ok <- !any(x0 < placed[, 2] & x1 > placed[, 1] & y0 < placed[, 4] & y1 > placed[, 3])
      if (ok) { ddx <- pmax(x0 - X, 0, X - x1); ddy <- pmax(y0 - Y, 0, Y - y1); ok <- !any(ddx^2 + ddy^2 < R^2) }
      if (ok) best <- list(ax = ax, ay = ay, hj = dr[3], vj = dr[4], b = c(x0, x1, y0, y1))
    }
    if (is.null(best)) { ax <- X[k]; ay <- Y[k] + R[k] + 0.8; best <- list(ax = ax, ay = ay, hj = 0.5, vj = 0, b = c(ax - w / 2, ax + w / 2, ay, ay + h)) }   # nothing free: above the node
    placed <- rbind(placed, best$b)
    out[[length(out) + 1]] <- data.frame(region = pts$region[k], ax = best$ax * ux, ay = best$ay * uy, hj = best$hj, vj = best$vj)
  }
  bind_rows(out)
}
net_plot <- function(L, tag, legend_title, panel_in = net_panel_in) {
  nodes <- L$nodes[order(-L$nodes$n), ]; e <- L$edges
  pad <- 0.16; rx <- range(nodes$x); ry <- range(nodes$y); bx <- rx + c(-1, 1) * pad * diff(rx); by <- ry + c(-1, 1) * pad * diff(ry)
  lab <- place_labels(nodes, top_ranges, bx, by, diff(bx) / (panel_in[1] * 25.4), diff(by) / (panel_in[2] * 25.4))
  ggplot() +
    geom_segment(data = e, aes(x, y, xend = xend, yend = yend), colour = "grey35", linewidth = 1.7, alpha = 0.3, lineend = "round") +          # an outline: pale colours stay visible
    geom_segment(data = e, aes(x, y, xend = xend, yend = yend, colour = value), linewidth = 1.1, alpha = 0.95, lineend = "round") +
    geom_point(data = nodes, aes(x, y, size = n, fill = class), shape = 21, colour = "grey15", stroke = 0.5) +
    geom_text(data = lab, aes(ax, ay, label = region, hjust = hj, vjust = vj), size = 3.3, fontface = "bold") +
    scale_fill_manual(values = c(range_cols, "Other ranges" = CFG$col_other_ranges), guide = "none") +
    scale_size_area(max_size = 11, name = "Alpine clusters\nin the range", labels = scales::comma) +
    scale_colour_gradientn(colours = c("#2C7BB6", "#FFFFBF", "#D7191C"), name = legend_title, limits = range(e$value)) +
    guides(colour = guide_colourbar(title.position = "top", barwidth = unit(0.2, "in"), barheight = unit(2.4, "in")), size = guide_legend(title.position = "top", ncol = 1)) +
    scale_x_continuous(limits = bx, expand = c(0, 0)) + scale_y_continuous(limits = by, expand = c(0, 0)) +
    labs(tag = tag, x = NULL, y = NULL) + theme_fig +
    theme(legend.position = "none", panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.border = element_rect(colour = "grey70"))
}
leg_cell <- function(p) {                                                                  # the legend of a panel as a panel of its own (to its right)
  g <- get_legend(p + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10)), "right")
  ggplot() + annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + theme_void()
}
L_main <- NULL; pD <- NULL; pD_leg <- NULL
if (!is.null(cell_stats) && length(reg_all) >= 3) {
  L_main <- net_layout("median_adj")
  message(sprintf("Range network: %d ranges (%d named), %d edges (the %d most similar ranges of each + the minimum spanning tree); the two MDS axes show %.0f%% of the variation of the median dissimilarities",
                  length(reg_all), length(top_ranges), L_main$n_edges, CFG$net_k, 100 * L_main$share))
  pD <- net_plot(L_main, "D", "Median Bray-Curtis\ndissimilarity of the pairs\nof clusters between\nranges (adjusted for\ndepth / unmapped,\nnot for distance)")
  pD_leg <- leg_cell(pD)
} else message("Range network: fewer than 3 ranges: panel D is left out")

stat_txt <- paste("Statistics: one linear regression per metric on the pairs of all replicates with a group x distance interaction (distance centred at the reference distance); the stars are the tests of its terms (level at the reference distance and rise per 10-fold distance of A = alpine, L = lowland; A-L = alpine minus lowland), with standard errors that are cluster-robust for pair data (pairs sharing a cluster are dependent; Aronow et al. 2015; the clusters are the independent units); *** p < 0.001, ** p < 0.01, * p < 0.05, ns; not adjusted for multiple tests. All coefficients: tables/model_coefficients.csv.")
caption <- paste(
  "Distance decay of community dissimilarity in alpine and lowland soils (samples within 20 km of an alpine area are not used), the ordination of the soils and the dissimilarity between mountain ranges.",
  sprintf("A, B: Bray-Curtis and Jaccard dissimilarity between pairs of samples of the same group against geographic distance (one straight line per group%s). Lines = regression, ribbons = 95%% interval, open circles = observed dissimilarity adjusted to the median differences in sequencing depth and unmapped proportion (mean of %d bins of distance with equal numbers of pairs, pooled over the replicates)%s.",
          if (CFG$df_dist > 1) sprintf("; natural spline with %d df", CFG$df_dist) else "", CFG$n_bins, if (isTRUE(CFG$rep_dots)) sprintf("; small dots = the same bin means in each of the %d replicates", CFG$n_rep) else ""),
  sprintf("Each of %d replicates: %d samples per group (one random sample per BioProject / location cluster, at least %g m apart, cut down evenly over %g-degree latitude / longitude cells) and up to %d pairs per group (random within classes of geographic distance); the pairs of all replicates are analysed together.", CFG$n_rep, n_main, CFG$thin_m, CFG$even_cell_deg, CFG$n_pairs),
  stat_txt,
  sprintf("C: UMAP of the community composition (Hellinger-transformed relative abundances of the %d most abundant OTUs, %d principal components) of %s alpine samples, coloured by %s (one colour family per continent: Europe reds, Asia blues, North America greens, South America oranges, Oceania purples, Africa golds; the largest range of a continent the darkest; the number of samples in brackets; the other ranges together as Other ranges), and %s lowland samples spread evenly over latitude / longitude (grey); small multiples: the first panel all alpine samples against the lowland, then one panel per continent with its ranges coloured and all other samples pale grey; ellipses = the central half of the samples (robust) of a range. UMAP keeps neighbourhoods, not distances: the axes, the size of the clusters and the gaps between them are not quantitative.",
          CFG$umap_top_otus, CFG$umap_pcs, format(sum(nm$group == "Alpine"), big.mark = ","), region_kind, format(sum(nm$group == "Lowland"), big.mark = ",")),
  umap_cap,
  if (!is.null(pD)) sprintf("D: network of the mountain ranges (%s; %d ranges, the %d coloured ones named; colours as in C, size of a node = number of alpine clusters). The nodes are placed by classical multidimensional scaling of the median Bray-Curtis dissimilarity of all pairs of alpine clusters of two ranges (%d replicates of one random alpine sample per cluster; every pair of clusters has one value, the mean over the replicates, so that every cluster counts equally), adjusted to the median differences in sequencing depth and unmapped proportion but NOT for geographic distance: ranges with similar communities lie close together, distance included (the two axes show %.0f%% of the variation; the dissimilarities are first reduced by a common floor so that their differences, not their common level, place the nodes, and nodes that overlap are moved slightly apart). Edges: every range linked to its %d most similar ranges and the minimum spanning tree, coloured by the median dissimilarity. Supplementary Figure 7 also shows the dissimilarity at the same geographic distance.", region_kind, length(reg_all), length(top_ranges), CFG$n_rep, 100 * L_main$share, CFG$net_k))
rows3 <- list(list(pA, pB), list(pC, pC_leg)); if (!is.null(pD)) rows3[[3]] <- list(pD, pD_leg)
assemble_fig(file.path(out_folder, "fig3_distance_decay.png"), rows3, get_legend(pA), caption, width = 14, height = if (is.null(pD)) 17 else 23.5, row_heights = c(1, 1.7, 1.3)[seq_along(rows3)],
             row_widths = list(NULL, c(2, 1), c(1.5, 1)))
saveRDS(list(A = pA, B = pB, C = pC, D = pD), file.path(out_folder, "fig3_panels.rds"))

# Supplementary 7: the network of the ranges without and with the adjustment for geographic distance (the layout of B rotated onto that of A)
if (!is.null(L_main)) {
  L_sp <- net_layout("median_adj_space", ref = L_main$xy_mds)
  pS7a <- net_plot(L_main, "A", "Median dissimilarity of\nthe pairs of clusters\nbetween ranges (adjusted\nfor depth / unmapped,\nnot for distance)")
  pS7b <- net_plot(L_sp, "B", "Median dissimilarity of\nthe pairs of clusters\nbetween ranges (adjusted\nfor depth / unmapped\nAND for distance)")
  cap7 <- paste("Supplementary Figure 7. The network of the mountain ranges (nodes: ranges, colours as in Figure 3C, size = number of alpine clusters; the coloured ranges are named) without (A, as in Figure 3D) and with (B) the adjustment for geographic distance.",
                "The nodes are placed by classical multidimensional scaling of the median Bray-Curtis dissimilarity of all pairs of alpine clusters of two ranges (replicates of one random alpine sample per cluster; every pair of clusters has one value, the mean over the replicates), edges link every range to its most similar ranges (and the minimum spanning tree), coloured by the median dissimilarity. In both panels the dissimilarity is adjusted to the median differences in sequencing depth and unmapped proportion; in B also to the median geographic distance of the pairs (a common slope of the dissimilarity on log10 distance, estimated within pairs of ranges), so B shows how different the ranges are at the same distance; its layout is rotated and scaled onto that of A to make the two comparable. Compare the two: ranges that are far apart in A and close together in B are far apart in space rather than different.",
                sprintf("The two axes show %.0f%% (A) and %.0f%% (B) of the variation. The numbers (raw and adjusted medians, median distance, number of pairs of clusters) are in tables/range_pair_dissimilarity.csv.", 100 * L_main$share, 100 * L_sp$share))
  assemble_fig(file.path(out_folder, "figS7_range_dissimilarity_space.png"), list(list(pS7a, leg_cell(pS7a)), list(pS7b, leg_cell(pS7b))), NULL, cap7, width = 14, height = 14.5, row_heights = c(1, 1), row_widths = list(c(1.5, 1), c(1.5, 1)))
}
# Supplementary 8: the distributions behind the medians, a grid of violins (the top ranges)
if (!is.null(pG)) {
  cap8 <- paste("Supplementary Figure 8. The distributions of the Bray-Curtis dissimilarities (unadjusted) of all pairs of alpine clusters within and between the mountain ranges with the most alpine samples (", length(reg_keep), " ranges; one value per pair of clusters = the mean over ", CFG$n_rep, " replicates of one random alpine sample per cluster, so that every cluster counts equally): one violin per cell (a cell with fewer than 5 pairs of clusters shows the values themselves; diagonal: within a range), the same vertical scale in all cells",
                if (identical(CFG$region_axis, "log1m")) ", stretched near 1" else "", " (values at the left of every row), coloured by the median (red = more dissimilar); rows and columns ordered by hierarchical clustering of the medians. The median distance of the pairs of every cell is in tables/range_pair_dissimilarity.csv.",
                if (any(reg_order %in% reg_warn)) paste0(" A dagger marks ranges with fewer than ", CFG$region_warn_clusters, " alpine clusters (", paste(reg_order[reg_order %in% reg_warn], collapse = ", "), "): their distributions rest on few pairs of clusters, and a range with 2 clusters on a single within-range pair, whose violin only shows the variation between the samples drawn for those clusters; interpret with care.") else "", sep = "")
  assemble_fig(file.path(out_folder, "figS8_range_violins.png"), list(list(pG)), NULL, cap8, width = 10, height = 11)
}

# Supplementary 9: the UMAP test and the spread, in the UMAP and in the space it was made from
{
  up <- UP[sample.int(nrow(UP), min(nrow(UP), 40000)), ]
  pairs_long <- bind_rows(lapply(names(spaces), function(sp) data.frame(space = factor(umap_space_lab[[sp]], levels = umap_space_lab), group = factor(up$group, levels = CFG$groups), d = up[[sp]])))
  cent_long <- bind_rows(lapply(names(spaces), function(sp) bind_rows(lapply(c("Lowland", "Alpine"), function(g) { C <- cs_list[[sp]][[g]]$C
    data.frame(space = factor(umap_space_lab[[sp]], levels = umap_space_lab), group = factor(g, levels = CFG$groups), d = sqrt(rowSums(sweep(C, 2, colMeans(C))^2))) }))))
  vbox <- list(geom_violin(scale = "width", colour = "grey30", linewidth = 0.3, alpha = 0.85), geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.3),
               scale_fill_manual(values = CFG$cols, guide = "none"), scale_x_discrete(labels = c(Lowland = "Lowland", Alpine = "Alpine")), facet_wrap(~ space, nrow = 1, scales = "free_y"))
  pS9a <- ggplot(pairs_long, aes(group, d, fill = group)) + vbox + stat_summary(fun = mean, geom = "point", shape = 23, size = 2.4, fill = "white", colour = "black") +
    labs(tag = "A", x = NULL, y = "Distance between two samples of different clusters") + theme_fig
  pS9b <- ggplot(cent_long, aes(group, d, fill = group)) + vbox + stat_summary(fun = mean, geom = "point", shape = 23, size = 2.4, fill = "white", colour = "black") +
    labs(tag = "B", x = NULL, y = "Distance of a cluster centre to the centre of its group") + theme_fig
  met_lab <- c(centroid_dispersion = "Spread of the cluster centres\n(distance to the group centre)", centroid_variance = "Total variance of the cluster centres",
               ellipse_area = "Area of the 95% ellipse of the\ncluster centres (2D only)", within_cluster_distance = "Distance of the samples to their\ncluster centre (clusters with 2+ samples)")
  fr <- spread_tab[spread_tab$metric %in% names(met_lab) & is.finite(spread_tab$estimate), ]; fr$metric <- factor(met_lab[fr$metric], levels = rev(met_lab)); fr$space <- factor(umap_space_lab[fr$space], levels = umap_space_lab)
  pS9c <- ggplot(fr, aes(estimate, metric, colour = space)) + geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_linerange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.55), linewidth = 0.7) + geom_point(position = position_dodge(width = 0.55), size = 2.2) +
    geom_text(aes(x = hi, label = stars(p)), position = position_dodge(width = 0.55), hjust = -0.25, vjust = 0.4, size = 2.6, show.legend = FALSE) +
    scale_x_log10(expand = expansion(mult = c(0.05, 0.12))) + scale_colour_manual(values = c("#0072B2", "#D55E00"), name = NULL) +
    labs(tag = "C", x = "Ratio alpine / lowland (< 1: alpine less spread)", y = NULL) + theme_fig + theme(legend.position = "top")
  sh <- spread_tab[spread_tab$metric == "between_cluster_share", ]
  sh_long <- bind_rows(data.frame(space = factor(umap_space_lab[sh$space], levels = umap_space_lab), group = factor("Alpine", levels = CFG$groups), share = sh$alpine), data.frame(space = factor(umap_space_lab[sh$space], levels = umap_space_lab), group = factor("Lowland", levels = CFG$groups), share = sh$lowland))
  pS9d <- ggplot(sh_long, aes(group, share, fill = group)) + geom_col(width = 0.6) + scale_fill_manual(values = CFG$cols, guide = "none") + scale_x_discrete(labels = c(Lowland = "Lowland", Alpine = "Alpine")) +
    facet_wrap(~ space, nrow = 1) + scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    labs(tag = "D", x = NULL, y = "Share of the variance between clusters") + theme_fig
  cap9 <- paste("Supplementary Figure 9. Are alpine samples more similar to each other than lowland samples, and how do the spreads compare? Left facets / points in orange: the 2D UMAP of Figure 3C; right facets / points in blue: the space the UMAP was made from (the principal components of the Hellinger-transformed relative abundances; Euclidean distance there = Hellinger distance between communities, not distorted like a UMAP).",
                "A: the distance between two samples of the same group from different clusters (one random sample per cluster in each of replicates; violin, box = quartiles, diamond = mean); the test of the difference of the means, with cluster-robust standard errors for pair data, is in tables/umap_similarity_test.csv and in Figure 3C. B: the distance of every cluster centre (the mean position of the samples of a cluster) to the centre of its group: the alpine and lowland clusters are the independent units, one point each.",
                "C: ratios alpine / lowland of the spread measures with 95% intervals and stars from a bootstrap that resamples whole clusters (*** p < 0.001, ** p < 0.01, * p < 0.05, ns). D: the share of the total variance of the samples that lies between clusters rather than within them: where it is high, a group is a collection of clumps (studies, places) and the pooled samples give little information beyond the clusters.")
  assemble_fig(file.path(out_folder, "figS9_umap_spread.png"), list(list(pS9a, pS9b), list(pS9c, pS9d)), NULL, cap9, width = 13, height = 11.5, row_heights = c(1, 1))
}

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
s6cap <- paste("Supplementary Figure 6. The overall dissimilarity of alpine and lowland soils, crude and adjusted.",
               sprintf("A: dissimilarity of the pairs of samples of each group (all replicates and distances; violin, box = quartiles, diamond = mean; %s pairs). B: the geographic distance of the same pairs (shaded: the range common to all groups): the groups do not cover the same distances, so the crude comparison in A mixes differences in the communities with differences in where the pairs are.", format(nrow(P0), big.mark = ",")),
               "C, D: differences between the groups (A = alpine, L = lowland) estimated by four linear regressions that add, one after the other, the differences in sequencing depth and unmapped proportion, the geographic distance (log10), and a distance slope per group (the main model of Figure 3, evaluated at the mean distance of the pairs in the common range); points = estimates, bars = 95% intervals, standard errors cluster-robust for pair data (pairs sharing a cluster are dependent; the clusters are the independent units); *** p < 0.001, ** p < 0.01, * p < 0.05, ns; not adjusted for multiple tests. A difference that shrinks towards 0 when distance is added reflects where the pairs are, not the communities.")
assemble_fig(file.path(out_folder, "figS6_overall_dissimilarity.png"), list(list(pS1, pS2), list(pS3, pS4)), get_legend(pS3 + theme(legend.position = "top")), s6cap, width = 13, height = 11.5)

message("Done. Outputs in ", out_folder)
writeLines(capture.output(sessionInfo()), file.path(out_folder, "sessionInfo.txt"))
file.create(file.path(out_folder, ".done"))   # flag file for Snakemake