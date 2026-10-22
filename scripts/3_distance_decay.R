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
#   * One LINEAR model per replicate and metric on the dissimilarity itself (0-1; a beta regression does not converge, many
#     dissimilarities are exactly 1): a distance curve per group (natural spline of log10 km, df_dist), an intercept per group and
#     linear nuisance terms for the differences in sequencing depth and unmapped proportion.
#   * What is compared, in every replicate, and then summarised over the replicates (median and range 2.5-97.5%):
#       - the dissimilarity at the SAME distances (ref_km), group against group  -> "more dissimilar across distance"
#       - the TURNOVER RATE: the rise of the dissimilarity per 10-fold increase of distance, over the distance range that all three
#         groups have in common (common_q), reported separately for its NEAR and its FAR half (split at the geometric mean).
#         A dissimilarity cannot exceed 1: a group whose dissimilarity is high saturates at short distances and then looks FLAT at
#         long distances, so only the near half tells whether dissimilarity rises faster       -> "rises faster with distance"
#     Contrasts: alpine - lowland, peri-alpine - lowland, alpine - peri-alpine; with the share of replicates in which each is > 0.
#
# CAVEATS
#   * The range between replicates shows how much the result depends on WHICH samples are drawn (alpine sample of each cluster, even
#     cut, pairs). It is NOT a confidence interval of the population (that would need the sampling of clusters themselves) and
#     pairs reuse samples: judge robustness by the size of the effects and their consistency between replicates.
#   * The groups differ in geography (alpine samples sit in mountain ranges), so "at the same distance" is not "in the same
#     environment": a higher alpine dissimilarity can also reflect environmental turnover between ranges. Distances are only compared
#     where all three groups have pairs (distance_range_per_group.csv).
#   * Pairs reuse samples (a sample is in several pairs and in several replicates).
#
# FIGURE (in <out_folder>)   fig3_distance_decay.png : A, B = Bray-Curtis / Jaccard dissimilarity against distance per group (lines =
#   median over the replicates, ribbons = range between replicates, open circles = observed dissimilarity adjusted to the median depth /
#   unmapped differences, mean in 8 bins of distance); C = dissimilarity at the reference distances; D = turnover rate per group
#   (near and far half of the common distance range).
# TABLES (in <out_folder>/tables)   curves.csv, levels_at_reference_distances.csv, level_contrasts.csv, turnover.csv,
#   turnover_contrasts.csv, per_replicate_summary.csv, model_checks.csv, samples_per_replicate.csv, pair_counts.csv,
#   distance_range_per_group.csv, share_of_pairs_at_dissimilarity_1.csv
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
  df_dist     = 3,              # degrees of freedom of the natural spline of log10 distance (3 = a flexible but smooth decay)
  ref_km      = c(100, 1000),   # reference distances (km) at which the groups are compared (only where all groups have pairs)
  common_q    = c(0.05, 0.95),  # the distance range common to all groups = from the highest 5th to the lowest 95th percentile
  points_max  = 40000,          # pairs drawn as points (random subset; the models use all pairs)
  point_alpha = 0.05,
  cols        = c(Lowland = "#7C9885", Perialpine = "#CC79A7", Alpine = "#2B6CB0")
)
if (TEST_MODE) {
  out_folder <- paste0(out_folder, "_TEST")
  CFG$n_per_group <- 60; CFG$n_rep <- 1; CFG$n_pairs <- 300
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

# Assemble rows of 2 panels + a shared legend + a caption into one PNG
assemble_fig <- function(file, rows, leg, caption, width, height, row_heights = rep(1, length(rows))) {
  nr <- length(rows)
  cap_lines <- strwrap(caption, 160)
  png(file, width = width, height = height, units = "in", res = 300)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nr + 2, 2,
                        heights = unit.c(unit(row_heights, "null"), unit(0.5, "in"),
                                         unit(0.17 * length(cap_lines) + 0.2, "in")))))
  for (i in seq_len(nr)) for (j in seq_along(rows[[i]])) {
    vp <- viewport(layout.pos.row = i, layout.pos.col = if (length(rows[[i]]) == 1) 1:2 else j)
    print(rows[[i]][[j]], vp = vp)
  }
  if (!is.null(leg)) {
    pushViewport(viewport(layout.pos.row = nr + 1, layout.pos.col = 1:2)); grid.draw(leg); popViewport()
  }
  pushViewport(viewport(layout.pos.row = nr + 2, layout.pos.col = 1:2))
  grid.text(paste(cap_lines, collapse = "\n"), x = 0.01, just = "left", gp = gpar(fontsize = 9))
  popViewport(); dev.off()
}

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
cat("\nSamples and clusters per group:\n")
print(meta %>% group_by(group) %>% summarise(samples = n(), clusters = n_distinct(cl), .groups = "drop") %>% as.data.frame())
if (any(table(meta$group) == 0)) stop("a group has no samples: ", paste(names(table(meta$group))[table(meta$group) == 0], collapse = ", "))

## ---- 2. Replicate subsamples: the same number of samples in every group ----------------------------------------------
# One random sample per cluster (a different one in every replicate), at most one per thin_m metres, then every group is cut down to
# the same n, evenly over latitude / longitude cells: n = the smallest group (the group that limits the comparison is reported).
# The pool sizes reported below are therefore capped at 3 x n_per_group for the largest groups.
sel <- lapply(seq_len(CFG$n_rep), function(r) {
  set.seed(100 + r)
  pools <- lapply(CFG$groups, function(g) {
    x <- meta[meta$group == g, ]; x <- x[sample.int(nrow(x)), ]; x <- x[!duplicated(x$cl), ]
    # a large group (lowland: thousands of clusters) is first cut evenly to 3 x n_per_group: thinning is slow on many samples and only
    # n_per_group of them are needed (the cut is even, so it does not favour any region)
    if (nrow(x) > 3 * CFG$n_per_group) x <- x[even_order(x$lat, x$lon)[seq_len(3 * CFG$n_per_group)], ]
    x[thin_sites(x$lat, x$lon, CFG$thin_m), ]
  })
  names(pools) <- CFG$groups
  n <- min(CFG$n_per_group, vapply(pools, nrow, 1))
  out <- bind_rows(lapply(pools, function(x) x[even_order(x$lat, x$lon)[seq_len(n)], ]))
  attr(out, "pool_sizes") <- vapply(pools, nrow, 1)
  out
})
cells_of <- function(x) length(unique(paste(floor(x$lat / CFG$even_cell_deg), floor(x$lon / CFG$even_cell_deg))))
samples_per_rep <- bind_rows(lapply(seq_along(sel), function(r) {
  x <- sel[[r]]; ps <- attr(x, "pool_sizes")
  bind_rows(lapply(CFG$groups, function(g) data.frame(replicate = r, group = g, pool_one_per_cluster = unname(ps[g]), samples = sum(x$group == g),
                                                      clusters = length(unique(x$cl[x$group == g])), cells = cells_of(x[x$group == g, ]))))
}))
write.csv(samples_per_rep, out_t("samples_per_replicate.csv"), row.names = FALSE)
ps1 <- attr(sel[[1]], "pool_sizes"); n1 <- sum(sel[[1]]$group == CFG$groups[1])
message(sprintf("Every replicate: %d samples per group (one per cluster, at least %g m apart, evenly cut over %g-degree cells). Pool of replicate 1 (one per cluster): %s. The limiting group: %s%s",
                n1, CFG$thin_m, CFG$even_cell_deg, paste(sprintf("%s %d", names(ps1), ps1), collapse = ", "),
                paste(names(ps1)[ps1 == min(ps1)], collapse = "/"), if (n1 == CFG$n_per_group) " (none: n_per_group is the limit)" else ""))
message(sprintf("  cells covered per group (median over replicates): %s",
                paste(sprintf("%s %d", CFG$groups, round(tapply(samples_per_rep$cells, factor(samples_per_rep$group, levels = CFG$groups), median))), collapse = ", ")))
stopifnot(all(tapply(samples_per_rep$samples, samples_per_rep$replicate, function(x) length(unique(x))) == 1))      # the groups really are the same size
stopifnot(all(vapply(sel, function(x) !any(duplicated(x$cl[x$group == "Alpine"])), TRUE)))                           # one alpine sample per cluster
if (n1 < 100) warning("fewer than 100 samples per group: the distance decays will be poorly determined")

all_ids <- unique(unlist(lapply(sel, `[[`, "id")))
message(length(all_ids), " distinct samples across ", CFG$n_rep, " replicate(s)")

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
make_pairs <- function(r) {
  set.seed(200 + r)
  S <- sel[[r]] %>% filter(id %in% colnames(X_all))
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
  pairs <- bind_rows(lapply(CFG$groups, function(g) {
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

  pairs$rep_id <- r
  pairs$geo_km <- geosphere::distHaversine(cbind(S$lon[pairs$i], S$lat[pairs$i]), cbind(S$lon[pairs$j], S$lat[pairs$j])) / 1000
  pairs$d_depth <- abs(S$log10_reads[pairs$i] - S$log10_reads[pairs$j])   # difference in sequencing depth
  pairs$d_unm   <- abs(S$unmapped[pairs$i] - S$unmapped[pairs$j])          # difference in unmapped proportion
  pairs[, c("rep_id", "group", "bray", "jaccard", "geo_km", "d_depth", "d_unm")]
}

# All replicates. Saved (pairs_all.rds) and reused on a rerun if the settings and the drawn samples are the same, so a problem in
# the models does not cost the distance computation again.
pairs_key <- list(dist_breaks = CFG$dist_breaks_km, n_rep = CFG$n_rep, n_per_group = CFG$n_per_group, n_pairs = CFG$n_pairs,
                  rare_depth = CFG$rare_depth, n_ids = length(all_ids), ids_chars = sum(nchar(all_ids)), thin_m = CFG$thin_m,
                  even_cell_deg = CFG$even_cell_deg, test = TEST_MODE, design = "decay-v1")
pairs_file <- file.path(out_folder, "pairs_all.rds")
P_all <- NULL
if (REUSE_FITS && file.exists(pairs_file)) {
  cached <- readRDS(pairs_file)
  if (identical(cached$key, pairs_key)) { P_all <- cached$pairs; message("Reusing saved pairs: ", pairs_file) }
}
if (is.null(P_all)) {
  pr <- list()
  for (r in seq_len(CFG$n_rep)) {
    message("Replicate ", r, ": computing dissimilarities")
    x <- make_pairs(r)
    if (is.null(x)) { message("  no pairs drawn: skipping replicate"); next }
    pr[[length(pr) + 1]] <- x
  }
  P_all <- bind_rows(pr); rm(pr)
  saveRDS(list(key = pairs_key, pairs = P_all), pairs_file)
}

P <- as.data.frame(P_all)
P$group <- factor(P$group, levels = CFG$groups)
P$log10_geo <- log10(pmax(P$geo_km, 0.1))
reps <- sort(unique(P$rep_id))
message(sprintf("%d pairs from %d replicates in %d groups", nrow(P), length(reps), nlevels(droplevels(P$group))))
print(table(group = P$group))
write.csv(P %>% count(rep_id, group), out_t("pair_counts.csv"), row.names = FALSE)
# share of pairs with dissimilarity exactly 1 (no shared OTU): the pile at the boundary
at_one <- bind_rows(lapply(names(metrics), function(mt)
  P %>% group_by(group) %>% summarise(share_at_1 = mean(.data[[mt]] >= 1 - 1e-12), n_pairs = n(), .groups = "drop") %>% mutate(metric = mt)))
write.csv(at_one, out_t("share_of_pairs_at_dissimilarity_1.csv"), row.names = FALSE)
for (i in which(at_one$metric == "bray")) message(sprintf("  %s: %.1f%% of the Bray-Curtis pairs are exactly 1", at_one$group[i], 100 * at_one$share_at_1[i]))

## ---- 4. Distance range, reference distances, turnover window -----------------------------------------------------------
# Each group only has pairs over some distances (alpine samples sit in mountain ranges, evenly cut lowland samples are far apart):
# groups are compared only where ALL of them have pairs
rng <- P %>% group_by(group) %>% summarise(min_km = min(geo_km), q05_km = quantile(geo_km, CFG$common_q[1]), median_km = median(geo_km),
                                           q95_km = quantile(geo_km, CFG$common_q[2]), max_km = max(geo_km), n_pairs = n(), .groups = "drop") %>% as.data.frame()
write.csv(rng, out_t("distance_range_per_group.csv"), row.names = FALSE)
common <- c(max(rng$q05_km), min(rng$q95_km))
message(sprintf("Distance range common to all groups (%g-%g percentile): %s to %s km", 100 * CFG$common_q[1], 100 * CFG$common_q[2], fmt_km(common[1]), fmt_km(common[2])))
refs <- CFG$ref_km[CFG$ref_km >= common[1] & CFG$ref_km <= common[2]]
if (!length(refs) && common[2] > common[1]) {
  refs <- signif(10^(log10(common[1]) + c(1, 2) / 3 * (log10(common[2]) - log10(common[1]))), 2)
  message("  none of ref_km lies in that range: the reference distances are ", paste(refs, collapse = " and "), " km")
}
use_turnover <- is.finite(common[1]) && common[2] / common[1] >= 10                  # at least one decade
if (!use_turnover) message("  the common range is narrower than one decade: no turnover rate")
mid <- sqrt(prod(common))                                                             # geometric mean: the near and the far half
windows <- list(near = c(common[1], mid), far = c(mid, common[2]))
win_lab <- c(near = sprintf("near: %s-%s km", fmt_km(common[1]), fmt_km(mid)), far = sprintf("far: %s-%s km", fmt_km(mid), fmt_km(common[2])))
med_depth <- median(P$d_depth); med_unm <- median(P$d_unm)

## ---- 5. One linear model per replicate and metric ------------------------------------------------------------------------
# dissimilarity ~ group * ns(log10 distance) + depth + unmapped differences (see the header). Predictions at median depth / unmapped
# differences, the same for all groups.
fit_decay <- function(metric, d) lm(as.formula(sprintf("%s ~ group * ns(log10_geo, df = %d) + d_depth + d_unm", metric, CFG$df_dist)), data = d)
pred_at <- function(m, groups, x) predict(m, data.frame(group = factor(groups, levels = levels(m$model$group)), log10_geo = x, d_depth = med_depth, d_unm = med_unm))
grid_g <- bind_rows(lapply(CFG$groups, function(g) { x <- P$log10_geo[P$group == g]
  data.frame(group = g, x = seq(quantile(x, 0.01), quantile(x, 0.99), length.out = 80)) }))
curves <- list(); levs <- list(); turn <- list(); checks <- list(); adj <- list()
for (metric in names(metrics)) {
  message("  fitting ", metric, " (", length(reps), " replicates)")
  for (r in reps) {
    Pr <- P[P$rep_id == r, ]; Pr$group <- droplevels(Pr$group)
    if (nlevels(Pr$group) < length(CFG$groups) || any(table(Pr$group) < 50)) { message("    replicate ", r, ": a group has too few pairs, skipped"); next }
    m <- tryCatch(fit_decay(metric, Pr), error = function(e) { message("    fit failed: ", conditionMessage(e)); NULL })
    if (is.null(m)) next
    key <- paste(metric, r)
    curves[[key]] <- data.frame(replicate = r, metric = metric, group = grid_g$group, x = grid_g$x, fit = clip01(pred_at(m, grid_g$group, grid_g$x)))
    if (length(refs)) {
      nd <- expand.grid(group = CFG$groups, ref_km = refs, stringsAsFactors = FALSE)
      levs[[key]] <- data.frame(replicate = r, metric = metric, nd, level = clip01(pred_at(m, nd$group, log10(nd$ref_km))))
    }
    if (use_turnover) turn[[key]] <- bind_rows(lapply(names(windows), function(w) {
      lo <- pred_at(m, CFG$groups, log10(windows[[w]][1])); hi <- pred_at(m, CFG$groups, log10(windows[[w]][2]))
      data.frame(replicate = r, metric = metric, group = CFG$groups, window = w, turnover = (hi - lo) / (log10(windows[[w]][2]) - log10(windows[[w]][1])))
    }))
    # observed dissimilarity of every pair adjusted to the median depth / unmapped differences (residual + the curve at the pair's distance):
    # binned, these are the open circles that the curve has to run through
    adj[[key]] <- data.frame(metric = metric, group = as.character(Pr$group), x = Pr$log10_geo, y = Pr[[metric]] - fitted(m) + pred_at(m, as.character(Pr$group), Pr$log10_geo))
    fv <- tapply(fitted(m), Pr$group, mean); ob <- tapply(Pr[[metric]], Pr$group, mean)
    checks[[key]] <- data.frame(replicate = r, metric = metric, group = names(ob), observed_mean = as.numeric(ob), fitted_mean = as.numeric(fv[names(ob)]), n_pairs = as.numeric(table(Pr$group)[names(ob)]))
  }
}
if (!length(curves)) stop("no model could be fitted: see the messages above")
checks <- bind_rows(checks); checks$difference <- checks$fitted_mean - checks$observed_mean
write.csv(checks, out_t("model_checks.csv"), row.names = FALSE)
worst <- max(abs(checks$difference), na.rm = TRUE)
message(sprintf("Model check: the fitted mean dissimilarity differs from the observed one by at most %.3f (largest over %d fits x groups)%s", worst, nrow(checks),
                if (worst > 0.02) "  -> LARGER THAN 0.02: look at tables/model_checks.csv, a fit may be wrong" else ""))
curves <- bind_rows(curves); levs <- bind_rows(levs); turn <- bind_rows(turn); adj_all <- bind_rows(adj)

## ---- 6. Summaries over the replicates and the contrasts between groups --------------------------------------------------
curves_sum <- curves %>% group_by(metric, group, x) %>% summarise(lo = qlo(fit), hi = qhi(fit), fit = median(fit), .groups = "drop") %>% as.data.frame()
write.csv(curves_sum %>% mutate(distance_km = 10^x), out_t("curves.csv"), row.names = FALSE)
cmp <- list(c("Alpine", "Lowland"), c("Perialpine", "Lowland"), c("Alpine", "Perialpine"))      # first minus second
summ <- function(d, val, by) d %>% group_by(across(all_of(by))) %>%
  summarise(lo = qlo(.data[[val]]), hi = qhi(.data[[val]]), share_of_replicates_above_0 = mean(.data[[val]] > 0), n_replicates = n(), median = median(.data[[val]]), .groups = "drop") %>% as.data.frame()
contrast <- function(d, val, by) bind_rows(lapply(cmp, function(cp) {
  a <- d[d$group == cp[1], ]; b <- d[d$group == cp[2], ]
  x <- merge(a, b, by = c("replicate", by), suffixes = c("_1", "_2")); x$diff <- x[[paste0(val, "_1")]] - x[[paste0(val, "_2")]]
  x$contrast <- sprintf("%s minus %s", cp[1], cp[2]); x
}))
if (nrow(levs)) {
  lv_sum <- levs %>% summ("level", c("metric", "group", "ref_km")); lv_sum$level <- lv_sum$median; lv_sum$median <- NULL
  write.csv(lv_sum, out_t("levels_at_reference_distances.csv"), row.names = FALSE)
  lv_con <- contrast(levs, "level", c("metric", "ref_km")) %>% summ("diff", c("metric", "contrast", "ref_km")); lv_con$difference <- lv_con$median; lv_con$median <- NULL
  write.csv(lv_con, out_t("level_contrasts.csv"), row.names = FALSE)
}
if (nrow(turn)) {
  tv_sum <- turn %>% summ("turnover", c("metric", "group", "window")); tv_sum$turnover <- tv_sum$median; tv_sum$median <- NULL
  tv_sum$window_km <- win_lab[tv_sum$window]
  write.csv(tv_sum, out_t("turnover.csv"), row.names = FALSE)
  tv_con <- contrast(turn, "turnover", c("metric", "window")) %>% summ("diff", c("metric", "window", "contrast")); tv_con$difference <- tv_con$median; tv_con$median <- NULL
  tv_con$window_km <- win_lab[tv_con$window]
  write.csv(tv_con, out_t("turnover_contrasts.csv"), row.names = FALSE)
}
write.csv(bind_rows(levs, if (nrow(turn)) turn %>% mutate(level = NA_real_, ref_km = NA_real_)), out_t("per_replicate_summary.csv"), row.names = FALSE)

fmt <- function(x) sprintf("%+.3f", x)
message("\nDo mountains look like islands of diversity? Differences between groups (median over replicates [range 2.5-97.5%], share of replicates > 0):")
for (mt in names(metrics)) {
  message("  ", c(bray = "Bray-Curtis", jaccard = "Jaccard (rarefied)")[[mt]], ":")
  if (nrow(levs)) for (i in which(lv_con$metric == mt)) message(sprintf("    dissimilarity at %g km | %s: %s [%s, %s], above 0 in %.0f%% of the replicates", lv_con$ref_km[i], lv_con$contrast[i],
                                                                        fmt(lv_con$difference[i]), fmt(lv_con$lo[i]), fmt(lv_con$hi[i]), 100 * lv_con$share_of_replicates_above_0[i]))
  if (nrow(turn)) for (i in which(tv_con$metric == mt)) message(sprintf("    turnover per 10-fold distance, %s | %s: %s [%s, %s], above 0 in %.0f%% of the replicates", tv_con$window_km[i], tv_con$contrast[i],
                                                                        fmt(tv_con$difference[i]), fmt(tv_con$lo[i]), fmt(tv_con$hi[i]), 100 * tv_con$share_of_replicates_above_0[i]))
}
message("(> 0 for alpine minus lowland = alpine samples MORE dissimilar / dissimilarity rising faster with distance = the island pattern; read the FAR window with care: dissimilarity saturates at 1)")

## ---- 7. Figure -------------------------------------------------------------------------------------------------------------------
x_scale <- scale_x_continuous(breaks = log10(c(1, 10, 100, 1000, 10000)), labels = c("1", "10", "100", "1,000", "10,000"))
col_scales <- list(scale_colour_manual(values = CFG$cols, name = NULL, breaks = CFG$groups, labels = c(Lowland = "Lowland", Perialpine = "Peri-alpine", Alpine = "Alpine")),
                   scale_fill_manual(values = CFG$cols, guide = "none"),
                   guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1, linewidth = 1.2, linetype = "solid"))))
set.seed(5)
pts <- P[sample.int(nrow(P), min(nrow(P), CFG$points_max)), ]
pts$group <- factor(pts$group, levels = CFG$groups)

decay_panel <- function(metric, tag) {
  cv <- curves_sum[curves_sum$metric == metric, ]; cv$group <- factor(cv$group, levels = CFG$groups)
  bn <- adj_all[adj_all$metric == metric, ] %>% group_by(group) %>% mutate(bin = ntile(x, 8)) %>% group_by(group, bin) %>%
    summarise(x = mean(x), y = mean(y), .groups = "drop") %>% as.data.frame(); bn$group <- factor(bn$group, levels = CFG$groups)
  ylo <- max(0, floor((min(bn$y, cv$lo, na.rm = TRUE) - 0.04) * 20) / 20)               # zoom: pairs below the window are not drawn (they are in the models)
  ggplot() +
    geom_point(data = pts, aes(log10_geo, .data[[metric]], colour = group), size = 0.3, alpha = CFG$point_alpha, shape = 16) +
    geom_ribbon(data = cv, aes(x, ymin = lo, ymax = hi, fill = group), alpha = 0.2) +
    geom_line(data = cv, aes(x, fit, colour = group), linewidth = 1.1) +
    geom_point(data = bn, aes(x, y, colour = group), shape = 21, fill = "white", size = 2.6, stroke = 1.1) +
    x_scale + scale_y_continuous(breaks = function(l) { b <- scales::extended_breaks()(l); b[b <= 1 + 1e-9] }) + coord_cartesian(ylim = c(ylo, 1.005)) + col_scales +
    labs(tag = tag, x = "Geographic distance between samples (km)", y = paste(metrics[[metric]], "dissimilarity")) + theme_fig
}
pA <- decay_panel("bray", "A"); pB <- decay_panel("jaccard", "B")

group_x <- scale_x_discrete(limits = CFG$groups, labels = c(Lowland = "Lowland", Perialpine = "Peri-alpine", Alpine = "Alpine"))
pC <- pD <- NULL
if (nrow(levs)) {
  d <- lv_sum; d$group <- factor(d$group, levels = CFG$groups); d$metric <- factor(metrics[d$metric], levels = metrics)
  d$ref <- factor(sprintf("at %s km", formatC(d$ref_km, format = "d", big.mark = ",")), levels = sprintf("at %s km", formatC(sort(unique(d$ref_km)), format = "d", big.mark = ",")))
  pC <- ggplot(d, aes(group, level, colour = group)) + geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.5, linewidth = 0.8) +
    facet_grid(metric ~ ref, scales = "free_y", switch = "y") + group_x + col_scales +
    labs(tag = "C", x = NULL, y = "Dissimilarity at the same distance") + theme_fig + theme(strip.placement = "outside", legend.position = "none")
}
if (nrow(turn)) {
  d <- tv_sum; d$group <- factor(d$group, levels = CFG$groups); d$metric <- factor(metrics[d$metric], levels = metrics)
  d$window_km <- factor(d$window_km, levels = win_lab)
  pD <- ggplot(d, aes(group, turnover, colour = group)) + geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.5, linewidth = 0.8) +
    facet_grid(metric ~ window_km, scales = "free_y", switch = "y") + group_x + col_scales +
    labs(tag = "D", x = NULL, y = "Dissimilarity rise per 10-fold distance") + theme_fig + theme(strip.placement = "outside", legend.position = "none")
}
caption <- paste(
  "Distance decay of community dissimilarity in alpine, peri-alpine (lowland within 20 km of an alpine area) and lowland soils: do mountains behave like islands of diversity (alpine samples more dissimilar at the same distance, or dissimilarity rising faster with distance)?",
  sprintf("A, B: Bray-Curtis and Jaccard dissimilarity between pairs of samples of the same group against geographic distance. Lines = median over %d replicates, ribbons = range between replicates (2.5-97.5%%), open circles = observed dissimilarity adjusted to the median differences in sequencing depth and unmapped proportion (mean in 8 bins of distance).", length(reps)),
  "C: dissimilarity predicted at the same distance; D: turnover = rise of the dissimilarity per 10-fold increase of distance in the near and the far half of the distance range common to all groups (a dissimilarity cannot exceed 1: groups with a high dissimilarity flatten at long distances); points = median, bars = range between replicates.",
  sprintf("Each replicate: %d samples per group (one random sample per BioProject / location cluster, at least %g m apart, cut down evenly over %g-degree latitude / longitude cells) and up to %d pairs per group (random within classes of geographic distance); a linear model per replicate with a natural-spline distance curve per group.",
          n1, CFG$thin_m, CFG$even_cell_deg, CFG$n_pairs),
  "The ranges show how much the result depends on which samples are drawn; they are not confidence intervals of the population (pairs reuse samples). Contrasts: tables/level_contrasts.csv, tables/turnover_contrasts.csv.")
bottom <- if (!is.null(pC) && !is.null(pD)) list(pC, pD) else if (!is.null(pC)) list(pC) else NULL
rows <- c(list(list(pA, pB)), if (!is.null(bottom)) list(bottom))
assemble_fig(file.path(out_folder, "fig3_distance_decay.png"), rows, get_legend(pA), caption, width = 14, height = if (length(rows) == 2) 12.5 else 7.5,
             row_heights = if (length(rows) == 2) c(1, 1.1) else 1)
saveRDS(list(A = pA, B = pB, C = pC, D = pD), file.path(out_folder, "fig3_panels.rds"))

message("Done. Outputs in ", out_folder)
writeLines(capture.output(sessionInfo()), file.path(out_folder, "sessionInfo.txt"))
file.create(file.path(out_folder, ".done"))   # flag file for Snakemake