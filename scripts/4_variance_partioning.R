# =============================================================================
# Script 4 - Variation partitioning: is the alpine soil microbiome better predicted by the ENVIRONMENT or by GEOGRAPHY?
#
# What is partitioned (alpine samples only; column perialpine_alpine == "Alpine")
#   * community COMPOSITION: the Hellinger-transformed relative abundances (square root) of the top_otus most abundant OTUs;
#   * ALPHA diversity: log10 richness, Shannon and Pielou evenness, each as the residual of a regression on log10 sequencing depth.
#   The variation of every response is split into: environment only, geography only, shared by both, and unexplained (adjusted R2, Ezekiel:
#   it corrects for the number of predictors, and its fractions are those of vegan::varpart; checked against it in the log). A second
#   partition splits the environment into climate, soil and terrain (15 fractions, summarised).
#
# Predictor sets (the columns of all_sampled_data.csv; edit env_sets to change them)
#   Climate : CHELSA bio01-bio19, fcf, fgd, scd and the monthly variables at the sampling date (clt, cmi, hurs, pet, pr, rsds, sfcWind, spei12, spi12,
#             tas, tasmax, tasmin, vpd)
#   Soil    : SoilGrids <property>_30cm_avg (bdod, cec, cfvo, clay, sand, silt, nitrogen, phh2o, soc, ocd, ocs)
#   Terrain : DEM elevation_m, slope_deg, aspect_northness, aspect_eastness, tpi_m, roughness_m
#   Space   : spatial eigenvectors (PCNM, vegan::pcnm) of the great-circle distances between the cluster centres: the k_space broadest scales
#             (among clusters); a local polynomial of the coordinates within a cluster (within clusters).
#   Variables with more than max_na missing values are dropped, the others are imputed with the median; right-skewed variables (skewness
#   above skew_log, no negative values) are log1p-transformed; everything is z-scored. Every environmental set is then reduced to its first
#   principal components (up to pcs_max, enough for pca_var of the variance): NO variable selection, so that the set with more candidate
#   variables does not win just by fitting noise, and the response is never used to choose predictors. The sensitivity grid (panel D) shows
#   whether the verdict depends on the number of components / eigenvectors.
#
# CLUSTERS ARE THE UNITS (cluster_id: samples of one BioProject or one place). The environment varies mostly BETWEEN studies, so two levels:
#   * AMONG clusters (the regional / global signal): one point per cluster (at least among_min_n samples): the mean of the responses and of the
#     predictors of its samples, the clusters counting equally; df = clusters - 1.
#   * WITHIN clusters (the local signal, e.g. an elevation gradient inside a study): clusters with at least within_min_n samples; responses and
#     predictors centred on their cluster mean; df = samples - clusters. The composition is represented by its within_pcs principal components
#     (the total variation is that of the full data), which makes the bootstrap affordable.
#   "Better predicted by the environment than by geography" is read from the pure fractions (environment only - geography only, which equals the difference of the marginal adjusted R2 of the two sets) and from the
#   fractions, with 95% intervals from a bootstrap that resamples whole clusters (boot resamples; the predictor sets are
#   those of the full data, so the intervals do not include the uncertainty of the principal components; the bootstrap distribution is centred on the estimate
#   because repeated clusters inflate the R2 of the bootstrap samples).
#
# LIMITS: environment and geography are correlated (that is the point of the shared fraction); the environmental data are 1-250 m grids at the
#   coordinates (and CHELSA at the sampling date, projected after 2010); study / method differences are confounded with place among clusters;
#   a response that is poorly explained by both is unexplained, not "random".
#
# Input : <data_folder>/sample_data_filtered.tsv, results/all_sampled_data.csv (climate + soil + terrain, written by download_earth_data) and the
#         filtered BIOM
# Usage : Rscript 4_variance_partitioning.R [sample_table] [all_sampled_data.csv] [filtered_biom] [out_folder]
# Output: <out_folder>/fig4_variance_partitioning.png and tables/varpart_fractions.csv, varpart_differences.csv, varpart_sensitivity.csv,
#         predictor_sets.csv, sample_counts.csv
# Needs : R packages vegan, geosphere, rbiom, Matrix, dplyr, ggplot2, RSpectra (conda: r-vegan r-geosphere r-rbiom r-matrix r-tidyverse r-rspectra)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(Matrix)
  library(rbiom)
  library(geosphere)
  library(vegan)
  library(grid)
})

TEST_MODE <- FALSE   # TRUE: a tiny run (few OTUs, few bootstraps) to check that it works

args <- commandArgs(trailingOnly = TRUE)
data_folder     <- "/home/renku/work/epfl-altshuler/GlobalAlpineSoilMicrobiome/data"
results_default <- file.path(dirname(data_folder), "results")
sample_table  <- if (length(args) >= 1) as.character(args[1]) else file.path(data_folder, "sample_data_filtered.tsv")
env_table     <- if (length(args) >= 2) as.character(args[2]) else file.path(results_default, "all_sampled_data.csv")
filtered_biom <- if (length(args) >= 3) as.character(args[3]) else
  file.path(data_folder, "microbeatlas", "otuTable_soilCurated_OTU97_nextMAPrelease_filtered.biom")
out_folder    <- if (length(args) >= 4) as.character(args[4]) else file.path(results_default, "4_variance_partitioning")

CFG <- list(
  sample_file = sample_table, env_file = env_table, biom_file = filtered_biom,
  id_table    = "Sample_acc",   # sample ID in the sample table ...
  id_env      = "site",         # ... and in all_sampled_data.csv
  biom_id     = "Accessions",   # column of the sample table that matches the BIOM sample IDs
  group_col   = "perialpine_alpine", alpine_value = "Alpine", cluster_col = "cluster_id",
  rare_depth  = 1000,           # samples with fewer reads are dropped
  top_otus    = 5000,           # the most abundant OTUs of the composition
  env_sets    = list(
    Climate = c("^bio[0-9]{2}$", "^(fcf|fgd|scd)$", "^(clt|cmi|hurs|pet|pr|rsds|sfcWind|spei12|spi12|tas|tasmax|tasmin|vpd)$"),
    Soil    = c("_30cm_avg$"),
    Terrain = c("^elevation_m$", "^slope_deg$", "^aspect_northness$", "^aspect_eastness$", "^tpi_m$", "^roughness_m$")),
  max_na      = 0.2,            # a variable with more missing values among the alpine samples is dropped
  skew_log    = 2,              # right-skewed variables (skewness above this, no negative values) are log1p-transformed
  pca_var     = 0.9,            # principal components of a set: enough for this share of its variance ...
  pcs_max     = 6,              # ... at most this many per set
  k_space     = 12,             # the broadest PCNM eigenvectors (spatial predictors among clusters)
  among_min_n = 3,              # a cluster needs at least this many samples to be a point of the among-cluster analysis
  within_min_n = 10,            # ... and this many for the within-cluster analysis
  within_pcs  = 200,            # principal components of the within-cluster composition
  boot        = 500,            # bootstrap resamples of the clusters
  sens_k_space = c(4, 8, 12, 20), sens_pcs_max = c(2, 4, 6),   # the grid of the sensitivity panel
  cols = c(Environment = "#009E73", Space = "#E69F00", Shared = "#9E9E9E", Unexplained = "#F0F0F0",
           Climate = "#D55E00", Soil = "#8C6D31", Terrain = "#56B4E9", SharedEnv = "#B5D99C")
)
if (TEST_MODE) { out_folder <- paste0(out_folder, "_TEST"); CFG$top_otus <- 300; CFG$boot <- 40; CFG$within_pcs <- 30; CFG$among_min_n <- 2; CFG$within_min_n <- 4
                 CFG$sens_k_space <- c(4, 8); CFG$sens_pcs_max <- c(2, 4) }
dir.create(file.path(out_folder, "tables"), recursive = TRUE, showWarnings = FALSE)
out_t <- function(...) file.path(out_folder, "tables", ...)
set.seed(1)

# BEGIN ENGINE
# ---- The partition engine ---------------------------------------------------------------------------------------------------------------------
# Variation partitioning by regression: the adjusted R2 (Ezekiel) of every union of the sets of predictors, and from them the fractions that belong to
# one set only, or to several (the commonality coefficients; for two sets the fractions [a], [b], [c] of vegan::varpart). Y: n x p (centred), the
# predictors centred, ss_tot = the total sum of squares of Y (needed separately when Y has been truncated to its main principal components), df_tot = the
# degrees of freedom of the total variation (n - 1; n - the number of clusters when the data are centred within clusters).
ss_expl <- function(X, Y) {                                                                   # explained sum of squares of Y by the columns of X, and the rank of X
  if (is.null(X) || ncol(X) == 0) return(c(ss = 0, rank = 0))
  q <- qr(X, tol = 1e-8); r <- q$rank; if (r == 0) return(c(ss = 0, rank = 0))
  c(ss = sum(qr.qty(q, Y)[seq_len(r), , drop = FALSE]^2), rank = r)
}
rsq_adj <- function(X, Y, ss_tot, df_tot) {
  e <- ss_expl(X, Y); if (df_tot - e[["rank"]] <= 0) return(NA_real_)
  1 - (1 - e[["ss"]] / ss_tot) * df_tot / (df_tot - e[["rank"]])
}
commonality <- function(sets, Y, ss_tot, df_tot) {                                            # sets: a named list of matrices
  k <- length(sets); nm <- names(sets); full <- 2^k - 1; bit <- 2^(seq_len(k) - 1); cache <- new.env()
  R <- function(mask) {                                                                       # the adjusted R2 of the union of the sets in `mask`
    if (mask == 0) return(0); key <- as.character(mask)
    if (!is.null(cache[[key]])) return(cache[[key]])
    val <- rsq_adj(do.call(cbind, sets[bitwAnd(mask, bit) > 0]), Y, ss_tot, df_tot); cache[[key]] <- val; val
  }
  fr <- bind_rows(lapply(seq_len(full), function(T) {
    ti <- which(bitwAnd(T, bit) > 0); tot <- 0
    for (U in 0:T) if (bitwAnd(U, T) == U) tot <- tot + (-1)^(length(ti) - sum(bitwAnd(U, bit) > 0)) * R(bitwXor(full, U))
    data.frame(fraction = paste(nm[ti], collapse = " + "), n_sets = length(ti), value = -tot, stringsAsFactors = FALSE)           # (the Moebius inversion of R2(all) - R2(all without U))
  }))
  list(fractions = fr, total = R(full), R = R, names = nm, bit = bit)
}
# Two sets: environment E and space S -> pure E, shared, pure S, and the marginal (total) adjusted R2 of each
partition2 <- function(E, S, Y, ss_tot, df_tot) {
  cm <- commonality(list(E = E, S = S), Y, ss_tot, df_tot); f <- setNames(cm$fractions$value, cm$fractions$fraction)
  c(pure_E = f[["E"]], shared = f[["E + S"]], pure_S = f[["S"]], total_E = cm$R(1), total_S = cm$R(2), explained = cm$total)
}
# END ENGINE

fmt <- function(x) sprintf("%+.3f", x)
theme_fig <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), plot.tag = element_text(face = "bold", size = 16), strip.background = element_rect(fill = "grey95"), strip.text = element_text(face = "bold"))
to_sparse <- function(m) {
  if (inherits(m, "simple_triplet_matrix")) Matrix::sparseMatrix(i = m$i, j = m$j, x = m$v, dims = c(m$nrow, m$ncol), dimnames = list(m$dimnames[[1]], m$dimnames[[2]]))
  else as(m, "dgCMatrix")
}
rowmean <- function(M, g) { M <- as.matrix(M); mu <- rowsum(M, g); mu / as.numeric(table(factor(g, levels = rownames(mu)))) }       # the mean of the rows of every group (sorted by group)
ctr_by <- function(M, g) { M <- as.matrix(M); M - rowmean(M, g)[match(g, rownames(rowmean(M, g))), , drop = FALSE] }                   # centred on the group mean

## ---- 1. Data ------------------------------------------------------------------------------------------------------------------------------
raw <- read.delim(CFG$sample_file, quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE, fill = FALSE)
need <- c(CFG$id_table, CFG$biom_id, CFG$group_col, CFG$cluster_col, "Latitude", "Longitude", "alpha_total_n_reads", "alpha_n_otus", "alpha_shannon", "alpha_pielou")
miss <- setdiff(need, names(raw)); if (length(miss)) stop("Missing columns in ", CFG$sample_file, ": ", paste(miss, collapse = ", "))
tab <- raw[raw[[CFG$group_col]] == CFG$alpine_value, ]
tab <- tab[is.finite(tab$Latitude) & is.finite(tab$Longitude) & is.finite(tab$alpha_total_n_reads) & tab$alpha_total_n_reads >= CFG$rare_depth &
           is.finite(tab$alpha_n_otus) & tab$alpha_n_otus > 0 & is.finite(tab$alpha_shannon) & is.finite(tab$alpha_pielou), ]
env <- read.csv(CFG$env_file, check.names = FALSE, stringsAsFactors = FALSE); env <- env[!duplicated(env[[CFG$id_env]]), ]
mm <- match(tab[[CFG$id_table]], env[[CFG$id_env]])
message(sprintf("%d alpine samples; %d of them have environmental data (%s)", nrow(tab), sum(!is.na(mm)), basename(CFG$env_file)))
tab <- tab[!is.na(mm), ]; env <- env[mm[!is.na(mm)], , drop = FALSE]
if (nrow(tab) < 50) stop("fewer than 50 alpine samples with environmental data: check CFG$id_table / CFG$id_env")

message("Reading the BIOM: ", CFG$biom_file)
biom <- rbiom::read_biom(CFG$biom_file); ids <- intersect(tab[[CFG$biom_id]], biom$samples); biom <- biom[ids]
X <- to_sparse(biom$counts); X <- X[rownames(X) != "Unmapped", , drop = FALSE]; rm(biom); invisible(gc())
keep <- tab[[CFG$biom_id]] %in% colnames(X); tab <- tab[keep, ]; env <- env[keep, , drop = FALSE]
X <- X[, tab[[CFG$biom_id]], drop = FALSE]; tot <- Matrix::colSums(X); ok <- tot >= CFG$rare_depth
tab <- tab[ok, ]; env <- env[ok, , drop = FALSE]; X <- X[, ok, drop = FALSE]; tot <- tot[ok]
Prel <- X %*% Matrix::Diagonal(x = 1 / tot); rs_otu <- Matrix::rowSums(Prel); top <- order(rs_otu, decreasing = TRUE)[seq_len(min(CFG$top_otus, nrow(Prel)))]
H <- sqrt(as.matrix(Matrix::t(Prel[top, , drop = FALSE]))); rm(X, Prel)                                         # samples x the most abundant OTUs, Hellinger-transformed
message(sprintf("%d alpine samples in the analysis; composition: %d OTUs (%.1f%% of all reads)", nrow(H), ncol(H), 100 * sum(rs_otu[top]) / sum(rs_otu)))

cl <- as.character(tab[[CFG$cluster_col]]); no_cl <- is.na(cl) | !nzchar(cl); cl[no_cl] <- paste0("none", seq_len(sum(no_cl)))
lon <- tab$Longitude; lat <- tab$Latitude
reads <- log10(tab$alpha_total_n_reads); dres <- function(y) as.numeric(resid(lm(y ~ poly(reads, 2))))         # the alpha metrics freed of the sequencing depth
A <- cbind(`Richness (log10)` = dres(log10(tab$alpha_n_otus)), Shannon = dres(tab$alpha_shannon), `Pielou evenness` = dres(tab$alpha_pielou))

# the environmental sets: columns, cleaning, z-scores
skew <- function(x) mean(((x - mean(x)) / sd(x))^3)
clean_env <- function(M) {
  M <- as.matrix(M); M <- M[, colMeans(is.na(M)) <= CFG$max_na, drop = FALSE]
  for (j in seq_len(ncol(M))) { x <- M[, j]; x[is.na(x)] <- median(x, na.rm = TRUE); if (sd(x) > 0 && min(x) >= 0 && skew(x) > CFG$skew_log) x <- log1p(x); M[, j] <- x }
  M <- M[, apply(M, 2, sd) > 1e-8, drop = FALSE]; if (ncol(M)) scale(M) else M
}
Z <- list()
for (nm in names(CFG$env_sets)) {
  cols <- unique(unlist(lapply(CFG$env_sets[[nm]], function(p) grep(p, names(env), value = TRUE)))); cols <- cols[vapply(env[cols], is.numeric, TRUE)]
  Z[[nm]] <- clean_env(env[, cols, drop = FALSE]); rownames(Z[[nm]]) <- NULL
  message(sprintf("  %s: %d of %d variables kept (missing values <= %.0f%%, not constant)", nm, ncol(Z[[nm]]), length(cols), 100 * CFG$max_na))
}
Z <- Z[vapply(Z, function(M) ncol(M) > 0, TRUE)]; if (!length(Z)) stop("no environmental variables found: check CFG$env_sets and the columns of ", CFG$env_file)

# a set reduced to its principal components (no use of the response)
env_pcs <- function(M, cap) {
  M <- M[, apply(M, 2, sd) > 1e-8, drop = FALSE]; if (ncol(M) == 0 || nrow(M) < 4) return(NULL)
  pc <- prcomp(M, scale. = TRUE); cum <- cumsum(pc$sdev^2) / sum(pc$sdev^2)
  k <- min(cap, max(2, which(cum >= CFG$pca_var)[1], na.rm = TRUE), ncol(pc$x)); P <- scale(pc$x[, seq_len(k), drop = FALSE], scale = FALSE)
  attr(P, "var") <- cum[k]; attr(P, "n_vars") <- ncol(M); P
}

## ---- 2. The two levels -----------------------------------------------------------------------------------------------------------------------
cl_n <- table(cl)
# AMONG clusters: one point per cluster
cl_a <- names(cl_n)[cl_n >= CFG$among_min_n]; sel_a <- cl %in% cl_a
Ha <- rowmean(H[sel_a, ], cl[sel_a]); Aa <- rowmean(A[sel_a, , drop = FALSE], cl[sel_a]); Za <- lapply(Z, function(M) rowmean(M[sel_a, , drop = FALSE], cl[sel_a]))
ll_a <- rowmean(cbind(lon, lat)[sel_a, ], cl[sel_a]); n_a <- nrow(Ha)
pcnm_all <- tryCatch(vegan::pcnm(as.dist(geosphere::distm(ll_a, fun = geosphere::distHaversine) / 1000))$vectors, error = function(e) NULL)
if (is.null(pcnm_all) || ncol(pcnm_all) < 2) stop("the spatial eigenvectors (PCNM) of the cluster centres could not be computed")
among_sets <- function(k_space = CFG$k_space, cap = CFG$pcs_max) {
  E <- lapply(Za, env_pcs, cap = cap); E <- E[!vapply(E, is.null, TRUE)]
  list(E = E, S = scale(pcnm_all[, seq_len(min(k_space, ncol(pcnm_all))), drop = FALSE], scale = FALSE))
}
as_ <- among_sets()
L_among <- list(level = "Among clusters", units = n_a, rows = as.list(seq_len(n_a)), E_sets = as_$E, S = as_$S,
                prep = function(Y, rs, idx, n_units) { Yb <- scale(Y[idx, , drop = FALSE], scale = FALSE)
                  list(Y = Yb, ss_tot = sum(Yb^2), df_tot = length(idx) - 1, E_sets = lapply(as_$E, function(M) scale(M[idx, , drop = FALSE], scale = FALSE)), S = scale(as_$S[idx, , drop = FALSE], scale = FALSE)) })
resp_a <- c(list(Composition = list(Y = Ha, rs = NULL)), lapply(setNames(seq_len(ncol(Aa)), colnames(Aa)), function(j) list(Y = Aa[, j, drop = FALSE], rs = NULL)))

# WITHIN clusters: everything centred on the cluster mean
cl_w <- names(cl_n)[cl_n >= CFG$within_min_n]; sel_w <- cl %in% cl_w; L_within <- NULL; resp_w <- NULL
if (length(cl_w) >= 5) {
  clw <- cl[sel_w]; n_w <- sum(sel_w)
  ll_c <- rowmean(cbind(lon, lat)[sel_w, ], clw)[match(clw, rownames(rowmean(cbind(lon, lat)[sel_w, ], clw))), , drop = FALSE]
  lx <- (lon[sel_w] - ll_c[, 1]) * cos(ll_c[, 2] * pi / 180) * 111.32; ly <- (lat[sel_w] - ll_c[, 2]) * 110.57                  # km from the centre of the cluster
  Sw <- ctr_by(cbind(lx, ly, lx^2, lx * ly, ly^2), clw); Sw <- Sw[, apply(Sw, 2, sd) > 1e-8, drop = FALSE]; Sw <- scale(Sw, center = FALSE)
  Ew <- lapply(Z, function(M) ctr_by(M[sel_w, , drop = FALSE], clw)); Hw <- ctr_by(H[sel_w, ], clw); Aw <- ctr_by(A[sel_w, , drop = FALSE], clw)
  sv <- RSpectra::svds(Hw, k = min(CFG$within_pcs, min(dim(Hw)) - 1)); Hred <- sv$u %*% diag(sv$d, nrow = length(sv$d))
  within_sets <- function(cap = CFG$pcs_max) { E <- lapply(Ew, env_pcs, cap = cap); E[!vapply(E, is.null, TRUE)] }
  ws_ <- within_sets()
  L_within <- list(level = "Within clusters", units = length(cl_w), rows = split(seq_len(n_w), factor(clw, levels = cl_w)), E_sets = ws_, S = Sw,
                   prep = function(Y, rs, idx, n_units) list(Y = Y[idx, , drop = FALSE], ss_tot = sum(rs[idx]), df_tot = length(idx) - n_units,
                                                             E_sets = lapply(ws_, function(M) M[idx, , drop = FALSE]), S = Sw[idx, , drop = FALSE]))
  resp_w <- c(list(Composition = list(Y = Hred, rs = rowSums(Hw^2))), lapply(setNames(seq_len(ncol(Aw)), colnames(Aw)), function(j) list(Y = Aw[, j, drop = FALSE], rs = Aw[, j]^2)))
} else message("Within clusters: fewer than 5 clusters with at least ", CFG$within_min_n, " samples: that level is left out")
message(sprintf("Among clusters: %d clusters (>= %d samples; %d samples); within clusters: %s", n_a, CFG$among_min_n, sum(sel_a),
                if (is.null(L_within)) "not analysed" else sprintf("%d clusters (>= %d samples; %d samples)", length(cl_w), CFG$within_min_n, sum(sel_w))))
set_info <- bind_rows(lapply(list(L_among, L_within), function(L) if (is.null(L)) NULL else bind_rows(
  lapply(names(L$E_sets), function(nm) data.frame(level = L$level, set = nm, variables = attr(L$E_sets[[nm]], "n_vars"), components = ncol(L$E_sets[[nm]]), variance_of_the_variables = attr(L$E_sets[[nm]], "var"))),
  data.frame(level = L$level, set = "Space", variables = NA, components = ncol(L$S), variance_of_the_variables = NA))))
write.csv(set_info, out_t("predictor_sets.csv"), row.names = FALSE)
write.csv(data.frame(level = c("Among clusters", "Within clusters"), clusters = c(n_a, length(cl_w)), samples = c(sum(sel_a), sum(sel_w))), out_t("sample_counts.csv"), row.names = FALSE)

## ---- 3. The partitions --------------------------------------------------------------------------------------------------------------------------
# A bootstrap that resamples whole clusters repeats some of them: the repeated rows make every regression look better than it is, and more so for the set with more predictors,
# so the bootstrap distribution of a difference of adjusted R2 is shifted. The interval keeps the SPREAD of the bootstrap but is centred on the estimate (the bias of the bootstrap mean removed).
bc_ci <- function(est, b) est + (unname(quantile(b, c(0.025, 0.975), na.rm = TRUE)) - mean(b, na.rm = TRUE))
analyse_level <- function(L, resp) {
  fr <- list(); df <- list()
  for (rn in names(resp)) {
    r <- resp[[rn]]; pr <- L$prep(r$Y, r$rs, unlist(L$rows), length(L$rows))
    E <- do.call(cbind, pr$E_sets); p2 <- partition2(E, pr$S, pr$Y, pr$ss_tot, pr$df_tot)
    cm4 <- commonality(c(pr$E_sets, list(Space = pr$S)), pr$Y, pr$ss_tot, pr$df_tot)
    bs <- vapply(seq_len(CFG$boot), function(i) { u <- sample(length(L$rows), replace = TRUE); idx <- unlist(L$rows[u]); b <- L$prep(r$Y, r$rs, idx, length(u))
      p <- partition2(do.call(cbind, b$E_sets), b$S, b$Y, b$ss_tot, b$df_tot); p[["pure_E"]] - p[["pure_S"]] }, 1)      # (= the difference of the marginal adjusted R2: the shared part cancels)
    fr[[length(fr) + 1]] <- bind_rows(
      data.frame(level = L$level, response = rn, partition = "Environment vs space", fraction = c("Environment only", "Shared", "Space only", "Environment (marginal)", "Space (marginal)", "All sets together"),
                 value = unname(p2[c("pure_E", "shared", "pure_S", "total_E", "total_S", "explained")]), stringsAsFactors = FALSE),
      data.frame(level = L$level, response = rn, partition = "Climate / soil / terrain / space", fraction = cm4$fractions$fraction, value = cm4$fractions$value, stringsAsFactors = FALSE))
    est <- p2[["pure_E"]] - p2[["pure_S"]]; ci <- bc_ci(est, bs)
    df[[length(df) + 1]] <- data.frame(level = L$level, response = rn, measure = "Environment - space", estimate = est, lo = ci[1], hi = ci[2], stringsAsFactors = FALSE)
  }
  list(fractions = bind_rows(fr), diffs = bind_rows(df))
}
message("Partitions and bootstraps (", CFG$boot, " resamples of the clusters per response and level) ...")
res_a <- analyse_level(L_among, resp_a); res_w <- if (!is.null(L_within)) analyse_level(L_within, resp_w) else NULL
fractions <- bind_rows(res_a$fractions, res_w$fractions); diffs <- bind_rows(res_a$diffs, res_w$diffs)
write.csv(fractions, out_t("varpart_fractions.csv"), row.names = FALSE); write.csv(diffs, out_t("varpart_differences.csv"), row.names = FALSE)

# the check against vegan::varpart (among clusters, composition, the two sets)
vp <- vegan::varpart(Ha, do.call(cbind, as_$E), as_$S)$part$indfract$Adj.R.squared
mine <- fractions[fractions$level == "Among clusters" & fractions$response == "Composition" & fractions$partition == "Environment vs space", ]
message(sprintf("Check against vegan::varpart (among clusters, composition): max |difference| of the three fractions = %.2e", max(abs(c(mine$value[mine$fraction == "Environment only"] - vp[1], mine$value[mine$fraction == "Space only"] - vp[2], mine$value[mine$fraction == "Shared"] - vp[3])))))

# sensitivity (among clusters, composition): the verdict for other numbers of components / eigenvectors
sens <- bind_rows(lapply(CFG$sens_k_space, function(ks) bind_rows(lapply(CFG$sens_pcs_max, function(cp) {
  st <- among_sets(ks, cp); Yc <- scale(Ha, scale = FALSE); p <- partition2(do.call(cbind, st$E), st$S, Yc, sum(Yc^2), nrow(Yc) - 1)
  data.frame(k_space = ks, pcs_max = cp, env_components = ncol(do.call(cbind, st$E)), pure_E = p[["pure_E"]], pure_S = p[["pure_S"]], shared = p[["shared"]], difference = p[["pure_E"]] - p[["pure_S"]], stringsAsFactors = FALSE) }))))
write.csv(sens, out_t("varpart_sensitivity.csv"), row.names = FALSE)

# the log
for (lv in unique(diffs$level)) for (rn in unique(diffs$response)) {
  f <- fractions[fractions$level == lv & fractions$response == rn & fractions$partition == "Environment vs space", ]; d <- diffs[diffs$level == lv & diffs$response == rn, ]
  v <- setNames(f$value, f$fraction)
  message(sprintf("  %s | %-16s environment only %5.1f%%, shared %5.1f%%, space only %5.1f%%, unexplained %5.1f%% | environment - space %s [%s, %s]", lv, rn, 100 * v[["Environment only"]], 100 * v[["Shared"]], 100 * v[["Space only"]],
                  100 * (1 - max(0, v[["All sets together"]])), fmt(d$estimate[1]), fmt(d$lo[1]), fmt(d$hi[1])))
}

## ---- 4. Figure --------------------------------------------------------------------------------------------------------------------------------------
resp_lab <- unique(fractions$response); lv_lab <- unique(fractions$level)
fr2 <- fractions[fractions$partition == "Environment vs space" & fractions$fraction %in% c("Environment only", "Shared", "Space only", "All sets together"), ]
barA <- bind_rows(lapply(split(fr2, list(fr2$level, fr2$response), drop = TRUE), function(d) { v <- setNames(d$value, d$fraction); pos <- pmax(v[c("Environment only", "Shared", "Space only")], 0)
  data.frame(level = d$level[1], response = d$response[1], part = c("Environment only", "Shared", "Space only", "Unexplained"), value = c(pos, max(0, 1 - sum(pos))), stringsAsFactors = FALSE) }))
barA$part <- factor(barA$part, levels = c("Unexplained", "Space only", "Shared", "Environment only")); barA$response <- factor(barA$response, levels = resp_lab); barA$level <- factor(barA$level, levels = lv_lab)
colA <- c("Environment only" = CFG$cols[["Environment"]], "Shared" = CFG$cols[["Shared"]], "Space only" = CFG$cols[["Space"]], "Unexplained" = CFG$cols[["Unexplained"]])
pA <- ggplot(barA, aes(response, value, fill = part)) + geom_col(width = 0.7, colour = "grey40", linewidth = 0.25) +
  geom_text(data = barA[barA$part != "Unexplained" & barA$value >= 0.02, ], aes(label = sprintf("%.0f%%", 100 * value)), position = position_stack(vjust = 0.5), size = 3, colour = "grey10") +
  scale_fill_manual(values = colA, breaks = names(colA), name = NULL) + facet_wrap(~ level) + scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
  labs(tag = "A", x = NULL, y = "Share of the variation (adjusted R2)") + theme_fig + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top")
f4 <- fractions[fractions$partition == "Climate / soil / terrain / space", ]
cat4 <- function(fr) { sets <- strsplit(fr, " \\+ "); vapply(sets, function(s) if (length(s) == 1) paste(s, "only") else if ("Space" %in% s) "Shared: environment + space" else "Shared: among environment sets", "") }
f4$cat <- cat4(f4$fraction)
barB <- f4 %>% group_by(level, response, cat) %>% summarise(value = sum(pmax(value, 0)), .groups = "drop") %>% as.data.frame()
barB <- bind_rows(barB, barB %>% group_by(level, response) %>% summarise(cat = "Unexplained", value = max(0, 1 - sum(value)), .groups = "drop"))
catlev <- c("Unexplained", "Space only", "Shared: environment + space", "Shared: among environment sets", "Terrain only", "Soil only", "Climate only")
barB <- barB[barB$cat %in% catlev, ]; barB$cat <- factor(barB$cat, levels = catlev); barB$response <- factor(barB$response, levels = resp_lab); barB$level <- factor(barB$level, levels = lv_lab)
colB <- c("Climate only" = CFG$cols[["Climate"]], "Soil only" = CFG$cols[["Soil"]], "Terrain only" = CFG$cols[["Terrain"]], "Shared: among environment sets" = CFG$cols[["SharedEnv"]],
          "Shared: environment + space" = CFG$cols[["Shared"]], "Space only" = CFG$cols[["Space"]], "Unexplained" = CFG$cols[["Unexplained"]])
pB <- ggplot(barB, aes(response, value, fill = cat)) + geom_col(width = 0.7, colour = "grey40", linewidth = 0.25) +
  geom_text(data = barB[barB$cat != "Unexplained" & barB$value >= 0.02, ], aes(label = sprintf("%.0f%%", 100 * value)), position = position_stack(vjust = 0.5), size = 2.8, colour = "grey10") +
  scale_fill_manual(values = colB, breaks = rev(catlev), name = NULL) + facet_wrap(~ level) + scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
  guides(fill = guide_legend(nrow = 3)) + labs(tag = "B", x = NULL, y = "Share of the variation (adjusted R2)") + theme_fig + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top", legend.text = element_text(size = 8))
dd <- diffs; dd$response <- factor(dd$response, levels = rev(resp_lab)); dd$level <- factor(dd$level, levels = lv_lab)
pC <- ggplot(dd, aes(estimate, response)) + geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_linerange(aes(xmin = lo, xmax = hi), colour = "#6A51A3", linewidth = 0.8) + geom_point(colour = "#6A51A3", size = 2.4) +
  geom_text(aes(x = hi, label = ifelse(abs(estimate) < 0.005, "0.00", sprintf("%+.2f", estimate))), hjust = -0.2, vjust = -0.6, size = 2.8, colour = "grey25") +
  facet_wrap(~ level, scales = "free_x") + scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(tag = "C", x = "Environment - space (adjusted R2; > 0: the environment is the better predictor)", y = NULL) + theme_fig
pD <- ggplot(sens, aes(factor(k_space), factor(pcs_max), fill = difference)) + geom_tile(colour = "white") + geom_text(aes(label = ifelse(abs(difference) < 0.005, "0.00", sprintf("%+.2f", difference))), size = 3.2) +
  scale_fill_gradient2(low = "#E69F00", mid = "white", high = "#009E73", midpoint = 0, name = "Environment only\n- space only") +
  labs(tag = "D", x = "Spatial eigenvectors (PCNM)", y = "Principal components per environmental set (cap)") + theme_fig + theme(panel.grid = element_blank())
cap <- paste("Figure 4. Is the alpine soil microbiome better predicted by the environment or by geography? Variation partitioning (adjusted R2) of the community composition (Hellinger-transformed relative abundances) and of the alpha diversity (richness, Shannon, Pielou evenness, freed of sequencing depth) of alpine soils, with the environment (CHELSA climate, SoilGrids soil, DEM terrain: principal components of each set) and geography (spatial eigenvectors among clusters, a local polynomial within clusters) as predictors.",
             "The clusters (one study / place) are the units: AMONG clusters, one point per cluster (the mean of its samples), the regional signal; WITHIN clusters, everything centred on the cluster mean, the local signal. A: environment only, shared and space only; B: the environment split into climate, soil and terrain (shared fractions summed; fractions below 0 are shown as 0). C: the difference between the environment and space (adjusted R2 of the environment minus that of space; it is the same for the pure fractions and for the marginal ones; > 0: the environment is the better predictor) with 95% intervals from a bootstrap that resamples whole clusters (the predictor sets are those of the full data; the spread of the bootstrap, centred on the estimate).",
             "D: the same difference for the composition among clusters with other numbers of components and spatial eigenvectors: the verdict should not depend on them. Environment and space are correlated (the shared fraction); the unexplained part is what neither explains, not necessarily noise.")
cap_lines <- strwrap(cap, 160)
png(file.path(out_folder, "fig4_variance_partitioning.png"), width = 13, height = 11.5 + 0.17 * length(cap_lines) * 0.4, units = "in", res = 300)
grid.newpage(); pushViewport(viewport(layout = grid.layout(3, 2, heights = unit.c(unit(1, "null"), unit(1, "null"), unit(0.17 * length(cap_lines) + 0.2, "in")))))
print(pA, vp = viewport(layout.pos.row = 1, layout.pos.col = 1)); print(pB, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(pC, vp = viewport(layout.pos.row = 2, layout.pos.col = 1)); print(pD, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1:2)); grid.text(paste(cap_lines, collapse = "\n"), x = 0.01, just = "left", gp = gpar(fontsize = 9)); popViewport(2); dev.off()
saveRDS(list(A = pA, B = pB, C = pC, D = pD), file.path(out_folder, "fig4_panels.rds"))
message("Done. Outputs in ", out_folder)
writeLines(capture.output(sessionInfo()), file.path(out_folder, "sessionInfo.txt"))
file.create(file.path(out_folder, ".done"))   # flag file for Snakemake