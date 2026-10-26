# =============================================================================
# Script 4 - Variation partitioning: is the soil microbiome better predicted by the ENVIRONMENT or by GEOGRAPHY, and is that different in
# ALPINE and in LOWLAND soils?
#
# What is partitioned
#   * community COMPOSITION: the Hellinger-transformed relative abundances (square root) of the top_otus most abundant OTUs (one OTU set for both
#     groups: the OTUs that are most abundant in the alpine and in the lowland samples, the two groups weighing equally);
#   * ALPHA diversity: log10 richness, Shannon and Pielou evenness, each as the residual of a regression on log10 sequencing depth.
#   The variation of every response is split into: environment only, geography only, shared by both, and unexplained (adjusted R2, Ezekiel:
#   it corrects for the number of predictors, and its fractions are those of vegan::varpart; checked against it in the log). A second
#   partition splits the environment into climate, soil and terrain (15 fractions, summarised).
#
# GROUPS: the alpine samples, and the lowland samples of the lowland REPLICATES of define_lowland_sets.py (columns rep_01, rep_02, ... of
#   all_sampled_data.csv: lowland samples spread evenly over 5-degree cells; soil is only downloaded for them). Every lowland replicate is
#   analysed on its own (lowland_reps of them) and summarised by the median (and the range) over replicates; the alpine analysis is repeated for
#   alpine_draws random draws of the samples that form the cluster means (below), for the same reason. The question of the comparison:
#   "environment - space" (adjusted R2 of the environment minus that of geography) in alpine minus the same in lowland soils.
#
# Predictor sets (the columns of all_sampled_data.csv; edit env_sets to change them)
#   Climate : CHELSA bio01-bio19, fcf, fgd, scd and the monthly variables at the sampling date (clt, cmi, hurs, pet, pr, rsds, sfcWind, spei12, spi12,
#             tas, tasmax, tasmin, vpd)
#   Soil    : SoilGrids <property>_30cm_avg (bdod, cec, cfvo, clay, sand, silt, nitrogen, phh2o, soc, ocd, ocs)
#   Terrain : DEM elevation_m, slope_deg, aspect_northness, aspect_eastness, tpi_m, roughness_m
#   Space   : spatial eigenvectors (PCNM, vegan::pcnm) of the great-circle distances between the cluster centres: the k_space broadest scales
#             (among clusters); a local polynomial of the coordinates within a cluster (within clusters).
#   Variables with more than max_na missing values are dropped, the others are imputed with the median; right-skewed variables (skewness
#   above skew_log, no negative values) are log1p-transformed; everything is z-scored (in each group on its own). Every environmental set is then
#   reduced to its first principal components (up to pcs_max, enough for pca_var of the variance): NO variable selection, so that the set with more
#   candidate variables does not win just by fitting noise, and the response is never used to choose predictors. The sensitivity grid (panel D)
#   shows whether the verdict depends on the number of components / eigenvectors.
#
# CLUSTERS ARE THE UNITS (cluster_id: samples of one BioProject or one place). The environment varies mostly BETWEEN studies, so two levels:
#   * AMONG clusters (the regional / global signal): one point per cluster with at least among_n_per_cluster samples: the mean of the responses and
#     of the predictors of exactly that many RANDOM samples of the cluster (so that every unit, in both groups, has the same sampling noise: the
#     lowland clusters are small, and the mean of 1-2 samples is noisier than the mean of 30); df = clusters - 1.
#   * WITHIN clusters (the local signal, e.g. an elevation gradient inside a study): clusters with at least within_min_n samples; responses and
#     predictors centred on their cluster mean; df = samples - clusters. The composition is represented by its within_pcs principal components
#     (the total variation is that of the full data), which makes the bootstrap affordable. This level may not exist for the lowland (few
#     clusters with many samples): it is then left out for that group.
#   "Better predicted by the environment than by geography" is read from the pure fractions (environment only - geography only, which equals the
#   difference of the marginal adjusted R2 of the two sets), with 95% intervals from a bootstrap that resamples whole clusters (and, for a
#   group with several replicates, picks a replicate at random every time; boot resamples). The predictor sets are those of the full data, so the
#   intervals do not include the uncertainty of the principal components; the bootstrap distribution is centred on the estimate because repeated
#   clusters inflate the R2 of the bootstrap samples. The difference between the groups uses independent bootstraps of the two groups.
#
# DIAGNOSTICS (log + tables/diagnostics_noise_ceiling_technical.csv; they change nothing in the analysis): the NOISE CEILING of the among-cluster level (the share of
#   the variation that is reproducible between two independent sets of among_n_per_cluster samples of the same cluster: no predictor can explain more), and the
#   share that goes with technical covariates (log10 reads and the share of unmapped reads, cluster means) alone and beyond environment + space. The unmapped READS are
#   not in the composition (the 'Unmapped' row of the BIOM is removed: relative abundances are shares of the mapped reads).
#
# OTU FILTER: an OTU is kept if it is present (at least one read) in at least min_samples (default 5) of the analysed samples (alpine + lowland replicates pooled), and in
#   at least min_prevalence of them (default 0); top_otus then caps their number (the most abundant ones; Inf = all).
#
# RESPONSE SENSITIVITY (fig4_response_sensitivity.png, tables/response_sensitivity.csv): the noise ceiling and the explained variance for other sets of OTUs (the n most abundant;
#   the OTUs present in at least x of the samples) and other transformations. Rare OTUs are mostly sampling noise: if the explained share rises when they are removed, set top_otus
#   and / or min_samples / min_prevalence accordingly.
#
# LIMITS: environment and geography are correlated (that is the point of the shared fraction); the environmental data are 1-250 m grids at the
#   coordinates (and CHELSA at the sampling date, projected after 2010); study / method differences are confounded with place among clusters;
#   the lowland replicates are evenly spread in space by design, the alpine samples are not, so the groups differ in spatial extent and in the
#   number of clusters (both are in tables/sample_counts.csv; the adjusted R2 corrects for the number of units, the intervals show the precision);
#   a response that is poorly explained by both is unexplained, not "random".
#
# Input : <data_folder>/sample_data_filtered.tsv, results/all_sampled_data.csv (climate + soil + terrain + the lowland replicates, written by
#         download_earth_data) and the filtered BIOM
# Usage : Rscript 4_variance_partitioning.R [sample_table] [all_sampled_data.csv] [filtered_biom] [out_folder]
# Output: <out_folder>/fig4_variance_partitioning.png and tables/varpart_fractions.csv, varpart_differences.csv (per group and alpine - lowland),
#         varpart_sensitivity.csv, predictor_sets.csv, sample_counts.csv, diagnostics_noise_ceiling_technical.csv, response_sensitivity.csv (+ fig4_response_sensitivity.png)
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
  group_col   = "perialpine_alpine", alpine_value = "Alpine", lowland_value = "Lowland", cluster_col = "cluster_id",
  lowland_reps = 5,             # lowland replicates analysed (rep_01, rep_02, ...; 0 = no lowland comparison)
  alpine_draws = 5,             # random draws of the samples that form the alpine cluster means
  rare_depth  = 1000,           # samples with fewer reads are dropped
  top_otus    = 5000,           # at most this many OTUs (the most abundant of those that pass the filters below; Inf = all of them)
  env_sets    = list(
    Climate = c("^bio[0-9]{2}$", "^(fcf|fgd|scd)$", "^(clt|cmi|hurs|pet|pr|rsds|sfcWind|spei12|spi12|tas|tasmax|tasmin|vpd)$"),
    Soil    = c("_30cm_avg$"),
    Terrain = c("^elevation_m$", "^slope_deg$", "^aspect_northness$", "^aspect_eastness$", "^tpi_m$", "^roughness_m$")),
  max_na      = 0.2,            # a variable with more missing values among the samples of a group is dropped
  skew_log    = 2,              # right-skewed variables (skewness above this, no negative values) are log1p-transformed
  pca_var     = 0.9,            # principal components of a set: enough for this share of its variance ...
  pcs_max     = 6,              # ... at most this many per set
  k_space     = 12,             # the broadest PCNM eigenvectors (spatial predictors among clusters)
  among_n_per_cluster = 3,      # among clusters: a cluster needs this many samples, and its unit is the mean of exactly this many random ones
  within_min_n = 10,            # a cluster needs at least this many samples for the within-cluster analysis
  within_pcs  = 200,            # principal components of the within-cluster composition
  boot        = 500,            # bootstrap resamples of the clusters
  min_samples = 5,              # the OTUs of the composition must be present (at least one read) in at least this many of the analysed samples (alpine + lowland pooled); 0 = no filter. Rare OTUs are mostly noise: see response_sensitivity
  min_prevalence = 0,           # ... and in at least this share of them (0 = no filter)
  response_sensitivity = TRUE,  # how the noise ceiling and the explained variance depend on the OTUs used (number, prevalence) and on the transformation (fig4_response_sensitivity.png)
  sens_otus   = c(50, 200, 1000, 5000, Inf), sens_prev = c(0.02, 0.05, 0.10, 0.25),   # ... the grid: the most abundant OTUs / the OTUs present in at least this share of the samples of the group
  sens_k_space = c(4, 8, 12, 20), sens_pcs_max = c(2, 4, 6),   # the grid of the sensitivity panel
  cols = c(Environment = "#009E73", Space = "#E69F00", Shared = "#9E9E9E", Unexplained = "#F0F0F0",
           Climate = "#D55E00", Soil = "#8C6D31", Terrain = "#56B4E9", SharedEnv = "#B5D99C"),
  cols_group = c(Alpine = "#2B6CB0", Lowland = "#7C9885", `Alpine - lowland` = "#222222")
)
if (TEST_MODE) { out_folder <- paste0(out_folder, "_TEST"); CFG$top_otus <- 300; CFG$boot <- 40; CFG$within_pcs <- 30; CFG$among_n_per_cluster <- 2; CFG$within_min_n <- 4
                 CFG$sens_k_space <- c(4, 8); CFG$sens_pcs_max <- c(2, 4); CFG$lowland_reps <- 2; CFG$alpine_draws <- 2; CFG$sens_otus <- c(50, Inf); CFG$sens_prev <- c(0.05, 0.25) }
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
# A bootstrap that resamples whole clusters repeats some of them: the repeated rows make every regression look better than it is, and more so for the set with more predictors,
# so the bootstrap distribution of a difference of adjusted R2 is shifted. The interval keeps the SPREAD of the bootstrap but is centred on the estimate (the bias of the bootstrap mean removed).
bc_ci <- function(est, b) est + (unname(quantile(b, c(0.025, 0.975), na.rm = TRUE)) - mean(b, na.rm = TRUE))

## ---- 1. Data ------------------------------------------------------------------------------------------------------------------------------
raw <- read.delim(CFG$sample_file, quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE, fill = FALSE)
need <- c(CFG$id_table, CFG$biom_id, CFG$group_col, CFG$cluster_col, "Latitude", "Longitude", "alpha_total_n_reads", "alpha_n_otus", "alpha_shannon", "alpha_pielou")
miss <- setdiff(need, names(raw)); if (length(miss)) stop("Missing columns in ", CFG$sample_file, ": ", paste(miss, collapse = ", "))
tab <- raw[raw[[CFG$group_col]] %in% c(CFG$alpine_value, CFG$lowland_value), ]
tab <- tab[is.finite(tab$Latitude) & is.finite(tab$Longitude) & is.finite(tab$alpha_total_n_reads) & tab$alpha_total_n_reads >= CFG$rare_depth &
           is.finite(tab$alpha_n_otus) & tab$alpha_n_otus > 0 & is.finite(tab$alpha_shannon) & is.finite(tab$alpha_pielou), ]
env <- read.csv(CFG$env_file, check.names = FALSE, stringsAsFactors = FALSE); env <- env[!duplicated(env[[CFG$id_env]]), ]
mm <- match(tab[[CFG$id_table]], env[[CFG$id_env]]); tab <- tab[!is.na(mm), ]; env <- env[mm[!is.na(mm)], , drop = FALSE]
is_alp <- tab[[CFG$group_col]] == CFG$alpine_value
message(sprintf("%d alpine and %d lowland samples have environmental data (%s)", sum(is_alp), sum(!is_alp), basename(CFG$env_file)))
if (sum(is_alp) < 50) stop("fewer than 50 alpine samples with environmental data: check CFG$id_table / CFG$id_env")
yes <- function(x) tolower(trimws(as.character(x))) %in% c("yes", "true", "1", "y")
rep_cols <- head(sort(grep("^rep_[0-9]+$", names(env), value = TRUE)), CFG$lowland_reps)
grp_idx <- list(Alpine = which(is_alp))
for (k in seq_len(CFG$alpine_draws - 1)) grp_idx[[paste0("Alpine|draw", k + 1)]] <- grp_idx$Alpine
for (rc in rep_cols) { ix <- which(!is_alp & yes(env[[rc]])); if (length(ix) >= 50) grp_idx[[paste0("Lowland|", rc)]] <- ix else message("  ", rc, ": fewer than 50 lowland samples, skipped") }
names(grp_idx)[1] <- "Alpine|draw1"
has_low <- any(grepl("^Lowland", names(grp_idx)))
if (!has_low) message("No lowland replicate with environmental data (columns rep_XX of ", basename(CFG$env_file), "): the analysis is alpine only")

message("Reading the BIOM: ", CFG$biom_file)
use <- sort(unique(unlist(grp_idx)))
biom <- rbiom::read_biom(CFG$biom_file); ids <- intersect(tab[[CFG$biom_id]][use], biom$samples); biom <- biom[ids]
X <- to_sparse(biom$counts); X <- X[rownames(X) != "Unmapped", , drop = FALSE]; rm(biom); invisible(gc())
col_of <- match(tab[[CFG$biom_id]], colnames(X)); tot <- Matrix::colSums(X)
good <- !is.na(col_of) & tot[pmax(col_of, 1)] >= CFG$rare_depth
grp_idx <- lapply(grp_idx, function(ix) ix[good[ix]]); grp_idx <- grp_idx[vapply(grp_idx, length, 1L) >= 50]
if (!any(grepl("^Alpine", names(grp_idx)))) stop("fewer than 50 alpine samples with reads in the BIOM")
Prel <- X %*% Matrix::Diagonal(x = 1 / tot); rm(X)
cols_of <- function(pat) unique(col_of[unlist(grp_idx[grepl(pat, names(grp_idx))])])
m_alp <- Matrix::rowMeans(Prel[, cols_of("^Alpine"), drop = FALSE])
m_low <- if (any(grepl("^Lowland", names(grp_idx)))) Matrix::rowMeans(Prel[, cols_of("^Lowland"), drop = FALSE]) else m_alp
rank_all <- order((m_alp + m_low) / 2, decreasing = TRUE)                                              # the OTUs by mean relative abundance (alpine and lowland weighing equally)
cols_used <- unique(col_of[unlist(grp_idx)])                                                              # the columns of the analysed samples
n_pres_all <- Matrix::rowSums(Prel[, cols_used, drop = FALSE] > 0)                                       # in how many of the analysed samples every OTU is present
keep_otu <- n_pres_all >= CFG$min_samples & n_pres_all / length(cols_used) >= CFG$min_prevalence
rank_all_f <- rank_all[keep_otu[rank_all]]
message(sprintf("OTU filter: %d of %d OTUs are present in at least %d of the %d analysed samples%s", length(rank_all_f), length(rank_all), CFG$min_samples, length(cols_used),
                if (CFG$min_prevalence > 0) sprintf(" and in at least %.0f%% of them", 100 * CFG$min_prevalence) else ""))
top <- rank_all_f[seq_len(min(CFG$top_otus, length(rank_all_f)))]                                       # ONE set of OTUs for both groups
message(sprintf("Composition: %d OTUs (%.1f%% of the alpine reads, %.1f%% of the lowland reads)", length(top), 100 * sum(m_alp[top]) / sum(m_alp), 100 * sum(m_low[top]) / sum(m_low)))
hell <- function(ix) sqrt(as.matrix(Matrix::t(Prel[top, col_of[ix], drop = FALSE])))                    # samples x OTUs, Hellinger-transformed

skew <- function(x) mean(((x - mean(x)) / sd(x))^3)
clean_env <- function(M) {
  M <- as.matrix(M); M <- M[, colMeans(is.na(M)) <= CFG$max_na, drop = FALSE]
  for (j in seq_len(ncol(M))) { x <- M[, j]; x[is.na(x)] <- median(x, na.rm = TRUE); if (sd(x) > 0 && min(x) >= 0 && skew(x) > CFG$skew_log) x <- log1p(x); M[, j] <- x }
  M <- M[, apply(M, 2, sd) > 1e-8, drop = FALSE]; if (ncol(M)) scale(M) else M
}
# a set reduced to its principal components (no use of the response)
env_pcs <- function(M, cap) {
  M <- M[, apply(M, 2, sd) > 1e-8, drop = FALSE]; if (ncol(M) == 0 || nrow(M) < 4) return(NULL)
  pc <- prcomp(M, scale. = TRUE); cum <- cumsum(pc$sdev^2) / sum(pc$sdev^2)
  k <- min(cap, max(2, which(cum >= CFG$pca_var)[1], na.rm = TRUE), ncol(pc$x)); P <- scale(pc$x[, seq_len(k), drop = FALSE], scale = FALSE)
  attr(P, "var") <- cum[k]; attr(P, "n_vars") <- ncol(M); P
}

## ---- 2. The two levels of one group ------------------------------------------------------------------------------------------------------------
build_group <- function(ix, with_within = TRUE) {
  t_g <- tab[ix, , drop = FALSE]; e_g <- env[ix, , drop = FALSE]; Hh <- hell(ix)
  cl <- as.character(t_g[[CFG$cluster_col]]); no_cl <- is.na(cl) | !nzchar(cl); cl[no_cl] <- paste0("none", seq_len(sum(no_cl)))
  lon <- t_g$Longitude; lat <- t_g$Latitude
  reads <- log10(t_g$alpha_total_n_reads); dres <- function(y) as.numeric(resid(lm(y ~ poly(reads, 2))))         # the alpha metrics freed of the sequencing depth
  A <- cbind(`Richness (log10)` = dres(log10(t_g$alpha_n_otus)), Shannon = dres(t_g$alpha_shannon), `Pielou evenness` = dres(t_g$alpha_pielou))
  Z <- list()
  for (nm in names(CFG$env_sets)) {
    cols <- unique(unlist(lapply(CFG$env_sets[[nm]], function(p) grep(p, names(e_g), value = TRUE)))); cols <- cols[vapply(e_g[cols], is.numeric, TRUE)]
    Z[[nm]] <- clean_env(e_g[, cols, drop = FALSE]); rownames(Z[[nm]]) <- NULL
  }
  Z <- Z[vapply(Z, function(M) ncol(M) > 0, TRUE)]; if (!length(Z)) stop("no environmental variables found: check CFG$env_sets and the columns of ", CFG$env_file)
  cl_n <- table(cl); m <- CFG$among_n_per_cluster
  # AMONG clusters: one unit per cluster = the mean of m random samples of it
  cl_a <- names(cl_n)[cl_n >= m]; ok_a <- which(cl %in% cl_a)
  pick <- unlist(lapply(split(ok_a, cl[ok_a]), function(i) sample(i, m))); cl_p <- cl[pick]
  Ha <- rowmean(Hh[pick, ], cl_p); Aa <- rowmean(A[pick, , drop = FALSE], cl_p); Za <- lapply(Z, function(M) rowmean(M[pick, , drop = FALSE], cl_p))
  ll_a <- rowmean(cbind(lon, lat)[pick, ], cl_p); n_a <- nrow(Ha)
  pcnm_all <- tryCatch(vegan::pcnm(as.dist(geosphere::distm(ll_a, fun = geosphere::distHaversine) / 1000))$vectors, error = function(e) NULL)
  if (n_a < 15 || is.null(pcnm_all) || ncol(pcnm_all) < 2) return(NULL)
  among_sets <- function(k_space = CFG$k_space, cap = CFG$pcs_max) {
    E <- lapply(Za, env_pcs, cap = cap); E <- E[!vapply(E, is.null, TRUE)]
    list(E = E, S = scale(pcnm_all[, seq_len(min(k_space, ncol(pcnm_all))), drop = FALSE], scale = FALSE))
  }
  as_ <- among_sets()
  # DIAGNOSTICS (nothing in the analysis depends on them): (1) the NOISE CEILING = the share of the among-cluster variation that is reproducible at all: two independent
  # sets of m samples of the same cluster (clusters with >= 2m samples) give two units; 1 - (noise variance) / (total variance across clusters); no predictor can explain
  # more than that; (2) the TECHNICAL covariates (log10 reads, share of unmapped reads; cluster means of the same samples): their adjusted R2 alone and beyond environment + space
  seed_keep <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)                 # the diagnostics draw random samples: the random stream of the analysis is restored after them
  diag <- tryCatch({
    resp_all <- c(list(Composition = Hh), lapply(setNames(seq_len(ncol(A)), colnames(A)), function(j) A[, j, drop = FALSE]))
    cl_2 <- names(cl_n)[cl_n >= 2 * m]
    ceil <- vapply(names(resp_all), function(rn) {
      if (length(cl_2) < 15) return(NA_real_); Yr <- resp_all[[rn]]; pcol <- ncol(Yr)
      two <- lapply(cl_2, function(cc) { sm <- sample(which(cl == cc)); list(sm[seq_len(m)], sm[m + seq_len(m)]) })
      U1 <- matrix(vapply(two, function(z) colMeans(Yr[z[[1]], , drop = FALSE]), numeric(pcol)), ncol = pcol, byrow = TRUE)
      U2 <- matrix(vapply(two, function(z) colMeans(Yr[z[[2]], , drop = FALSE]), numeric(pcol)), ncol = pcol, byrow = TRUE)
      Tv <- (sum(apply(U1, 2, var)) + sum(apply(U2, 2, var))) / 2; Nv <- mean(rowSums((U1 - U2)^2)) / 2
      max(0, min(1, 1 - Nv / Tv)) }, 1)
    tech <- NULL
    if (all(c("unmapped_raw_counts", "total_raw_counts") %in% names(t_g))) {
      tv <- cbind(log10_reads = log10(t_g$alpha_total_n_reads), unmapped = t_g$unmapped_raw_counts / t_g$total_raw_counts)
      for (j in 1:2) tv[!is.finite(tv[, j]), j] <- median(tv[is.finite(tv[, j]), j])
      Ta <- scale(rowmean(tv[pick, , drop = FALSE], cl_p), scale = FALSE); Eall <- do.call(cbind, as_$E)
      ra <- c(list(Composition = Ha), lapply(setNames(seq_len(ncol(Aa)), colnames(Aa)), function(j) Aa[, j, drop = FALSE]))
      tech <- bind_rows(lapply(names(ra), function(rn) { Y <- scale(ra[[rn]], scale = FALSE); cm <- commonality(list(E = Eall, S = as_$S, T = Ta), Y, sum(Y^2), nrow(Y) - 1)
        data.frame(response = rn, env_space = cm$R(3), tech_alone = cm$R(4), tech_beyond_env_space = cm$fractions$value[cm$fractions$fraction == "T"]) }))
    }
    list(ceiling = data.frame(response = names(ceil), noise_ceiling = unname(ceil), clusters_used = length(cl_2)), tech = tech)
  }, error = function(e) NULL)
  if (!is.null(seed_keep)) assign(".Random.seed", seed_keep, envir = globalenv())
  L_a <- list(level = "Among clusters", rows = as.list(seq_len(n_a)), E_sets = as_$E, S = as_$S,
              prep = function(Y, rs, idx, n_units) { Yb <- scale(Y[idx, , drop = FALSE], scale = FALSE)
                list(Y = Yb, ss_tot = sum(Yb^2), df_tot = length(idx) - 1, E_sets = lapply(as_$E, function(M) scale(M[idx, , drop = FALSE], scale = FALSE)), S = scale(as_$S[idx, , drop = FALSE], scale = FALSE)) })
  resp_a <- c(list(Composition = list(Y = Ha, rs = NULL)), lapply(setNames(seq_len(ncol(Aa)), colnames(Aa)), function(j) list(Y = Aa[, j, drop = FALSE], rs = NULL)))
  # WITHIN clusters: everything centred on the cluster mean
  cl_w <- names(cl_n)[cl_n >= CFG$within_min_n]; sel_w <- cl %in% cl_w; wi <- NULL
  if (with_within && length(cl_w) >= 5) {
    clw <- cl[sel_w]; n_w <- sum(sel_w); llw <- cbind(lon, lat)[sel_w, , drop = FALSE]; ll_c <- rowmean(llw, clw)[match(clw, rownames(rowmean(llw, clw))), , drop = FALSE]
    lx <- (llw[, 1] - ll_c[, 1]) * cos(ll_c[, 2] * pi / 180) * 111.32; ly <- (llw[, 2] - ll_c[, 2]) * 110.57                  # km from the centre of the cluster
    Sw <- ctr_by(cbind(lx, ly, lx^2, lx * ly, ly^2), clw); Sw <- Sw[, apply(Sw, 2, sd) > 1e-8, drop = FALSE]; Sw <- scale(Sw, center = FALSE)
    Ew <- lapply(Z, function(M) ctr_by(M[sel_w, , drop = FALSE], clw)); Hw <- ctr_by(Hh[sel_w, ], clw); Aw <- ctr_by(A[sel_w, , drop = FALSE], clw)
    sv <- RSpectra::svds(Hw, k = min(CFG$within_pcs, min(dim(Hw)) - 1)); Hred <- sv$u %*% diag(sv$d, nrow = length(sv$d))
    ws_ <- lapply(Ew, env_pcs, cap = CFG$pcs_max); ws_ <- ws_[!vapply(ws_, is.null, TRUE)]
    if (length(ws_) && ncol(Sw)) {
      L_w <- list(level = "Within clusters", rows = split(seq_len(n_w), factor(clw, levels = cl_w)), E_sets = ws_, S = Sw,
                  prep = function(Y, rs, idx, n_units) list(Y = Y[idx, , drop = FALSE], ss_tot = sum(rs[idx]), df_tot = length(idx) - n_units,
                                                            E_sets = lapply(ws_, function(M) M[idx, , drop = FALSE]), S = Sw[idx, , drop = FALSE]))
      resp_w <- c(list(Composition = list(Y = Hred, rs = rowSums(Hw^2))), lapply(setNames(seq_len(ncol(Aw)), colnames(Aw)), function(j) list(Y = Aw[, j, drop = FALSE], rs = Aw[, j]^2)))
      wi <- list(L = L_w, resp = resp_w, clusters = length(cl_w), samples = n_w)
    }
  }
  info <- bind_rows(lapply(names(as_$E), function(nm) data.frame(level = "Among clusters", set = nm, variables = attr(as_$E[[nm]], "n_vars"), components = ncol(as_$E[[nm]]), variance_of_the_variables = attr(as_$E[[nm]], "var"))),
                    data.frame(level = "Among clusters", set = "Space", variables = NA, components = ncol(as_$S), variance_of_the_variables = NA))
  rm(list = intersect(c("Hh", "Hw", "Aw", "Ew"), ls()))             # the closures below keep this whole environment alive: free the big sample x OTU matrices first
  list(among = list(L = L_a, resp = resp_a), within = wi, Ha = Ha, among_sets = among_sets, info = info, diag = diag, units = list(ix = ix, cl = cl, pick = pick, cl_p = cl_p),
       counts = data.frame(samples = length(ix), clusters_among = n_a, mean_cluster_size = mean(as.numeric(cl_n[cl_a])), clusters_within = if (is.null(wi)) 0 else wi$clusters, samples_within = if (is.null(wi)) 0 else wi$samples))
}
message("Building the groups (", length(grp_idx), " data sets: alpine draws and lowland replicates) ...")
G <- list(); for (nm in names(grp_idx)) { g <- build_group(grp_idx[[nm]], with_within = !grepl("draw[2-9]", nm)); if (is.null(g)) message("  ", nm, ": fewer than 15 clusters with ", CFG$among_n_per_cluster, " samples: skipped") else G[[nm]] <- g }
grp_of <- sub("\\|.*", "", names(G)); if (!"Alpine" %in% grp_of) stop("the alpine group has fewer than 15 clusters with ", CFG$among_n_per_cluster, " samples")
counts <- bind_rows(lapply(names(G), function(nm) cbind(group = sub("\\|.*", "", nm), data_set = nm, G[[nm]]$counts)))
write.csv(counts, out_t("sample_counts.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(names(G), function(nm) cbind(group = sub("\\|.*", "", nm), data_set = nm, G[[nm]]$info))), out_t("predictor_sets.csv"), row.names = FALSE)
for (gr in unique(grp_of)) { cc <- counts[counts$group == gr, ]
  message(sprintf("  %s: %d data set(s); among clusters: %.0f clusters (median) with the mean of %d samples; within clusters: %s", gr, nrow(cc), median(cc$clusters_among), CFG$among_n_per_cluster,
                  if (all(cc$clusters_within == 0)) "no cluster with enough samples" else sprintf("%.0f clusters, %.0f samples (median of the data sets that have them)", median(cc$clusters_within[cc$clusters_within > 0]), median(cc$samples_within[cc$clusters_within > 0])))) }

## ---- 3. The partitions --------------------------------------------------------------------------------------------------------------------------
message("Partitions and bootstraps (", CFG$boot, " resamples of the clusters per group, level and response) ...")
fit_p2 <- function(Lr, rn, resample = FALSE) {                                                  # the two-set partition of one response, for the full data or a bootstrap sample
  L <- Lr$L; r <- Lr$resp[[rn]]
  if (resample) { u <- sample(length(L$rows), replace = TRUE); idx <- unlist(L$rows[u]); nu <- length(u) } else { idx <- unlist(L$rows); nu <- length(L$rows) }
  pr <- L$prep(r$Y, r$rs, idx, nu); list(pr = pr, p2 = partition2(do.call(cbind, pr$E_sets), pr$S, pr$Y, pr$ss_tot, pr$df_tot))
}
analyse_group <- function(Gs, level_key, group_name) {
  Ls <- lapply(Gs, function(g) g[[level_key]]); Ls <- Ls[!vapply(Ls, is.null, TRUE)]
  if (!length(Ls)) return(NULL)
  fr <- list(); ds <- list(); bsl <- list()
  for (rn in names(Ls[[1]]$resp)) {
    per <- lapply(Ls, function(Lr) { f <- fit_p2(Lr, rn); cm4 <- commonality(c(f$pr$E_sets, list(Space = f$pr$S)), f$pr$Y, f$pr$ss_tot, f$pr$df_tot); list(p2 = f$p2, f4 = cm4$fractions) })
    P2 <- do.call(rbind, lapply(per, `[[`, "p2")); nrep <- nrow(P2)
    f2 <- data.frame(fraction = c("Environment only", "Shared", "Space only", "Environment (marginal)", "Space (marginal)", "All sets together"),
                     value = apply(P2[, c("pure_E", "shared", "pure_S", "total_E", "total_S", "explained"), drop = FALSE], 2, median),
                     vmin = apply(P2[, c("pure_E", "shared", "pure_S", "total_E", "total_S", "explained"), drop = FALSE], 2, min),
                     vmax = apply(P2[, c("pure_E", "shared", "pure_S", "total_E", "total_S", "explained"), drop = FALSE], 2, max), partition = "Environment vs space", stringsAsFactors = FALSE)
    f4 <- bind_rows(lapply(per, `[[`, "f4")) %>% group_by(fraction) %>% summarise(value = median(value), vmin = min(value), vmax = max(value), .groups = "drop") %>% as.data.frame(); f4$partition <- "Climate / soil / terrain / space"
    fr[[length(fr) + 1]] <- cbind(group = group_name, level = Ls[[1]]$L$level, response = rn, n_data_sets = nrep, rbind(f2, f4))
    bs <- vapply(seq_len(CFG$boot), function(i) { f <- fit_p2(Ls[[sample.int(length(Ls), 1)]], rn, resample = TRUE)$p2; f[["pure_E"]] - f[["pure_S"]] }, 1)   # (= the difference of the marginal adjusted R2: the shared part cancels)
    est <- median(P2[, "pure_E"] - P2[, "pure_S"]); ci <- bc_ci(est, bs)
    ds[[length(ds) + 1]] <- data.frame(group = group_name, level = Ls[[1]]$L$level, response = rn, estimate = est, lo = ci[1], hi = ci[2], stringsAsFactors = FALSE)
    bsl[[paste(Ls[[1]]$L$level, rn, sep = "|")]] <- list(bs = bs, est = est)
  }
  list(fractions = bind_rows(fr), diffs = bind_rows(ds), boots = bsl)
}
res <- list()
for (gr in unique(grp_of)) for (lk in c("among", "within")) { r <- analyse_group(G[grp_of == gr], lk, gr); if (!is.null(r)) res[[paste(gr, lk)]] <- r }
fractions <- bind_rows(lapply(res, `[[`, "fractions")); diffs <- bind_rows(lapply(res, `[[`, "diffs"))
# alpine - lowland: independent bootstraps of the two groups
if (has_low && any(grp_of == "Lowland")) {
  dd2 <- list()
  for (lk in c("among", "within")) { ra <- res[[paste("Alpine", lk)]]; rl <- res[[paste("Lowland", lk)]]
    if (is.null(ra) || is.null(rl)) next
    for (key in intersect(names(ra$boots), names(rl$boots))) { ea <- ra$boots[[key]]$est - rl$boots[[key]]$est; ci <- bc_ci(ea, ra$boots[[key]]$bs - rl$boots[[key]]$bs)
      dd2[[length(dd2) + 1]] <- data.frame(group = "Alpine - lowland", level = sub("\\|.*", "", key), response = sub(".*\\|", "", key), estimate = ea, lo = ci[1], hi = ci[2], stringsAsFactors = FALSE) } }
  diffs <- bind_rows(diffs, dd2)
}
write.csv(fractions, out_t("varpart_fractions.csv"), row.names = FALSE); write.csv(diffs, out_t("varpart_differences.csv"), row.names = FALSE)

# the check against vegan::varpart (alpine, among clusters, composition, first draw, the two sets)
g1 <- G[[which(grp_of == "Alpine")[1]]]; f1 <- fit_p2(list(L = g1$among$L, resp = g1$among$resp), "Composition")
vp <- vegan::varpart(g1$Ha, do.call(cbind, g1$among$L$E_sets), g1$among$L$S)$part$indfract$Adj.R.squared
message(sprintf("Check against vegan::varpart (alpine, among clusters, composition, first draw): max |difference| of the three fractions = %.2e",
                max(abs(c(f1$p2[["pure_E"]] - vp[1], f1$p2[["pure_S"]] - vp[2], f1$p2[["shared"]] - vp[3])))))

# sensitivity (among clusters, composition): the verdict for other numbers of components / eigenvectors, median over the data sets of a group
sens <- bind_rows(lapply(unique(grp_of), function(gr) bind_rows(lapply(CFG$sens_k_space, function(ks) bind_rows(lapply(CFG$sens_pcs_max, function(cp) {
  v <- vapply(G[grp_of == gr], function(g) { st <- g$among_sets(ks, cp); Yc <- scale(g$Ha, scale = FALSE); p <- partition2(do.call(cbind, st$E), st$S, Yc, sum(Yc^2), nrow(Yc) - 1); p[["pure_E"]] - p[["pure_S"]] }, 1)
  data.frame(group = gr, k_space = ks, pcs_max = cp, difference = median(v), stringsAsFactors = FALSE) }))))))
write.csv(sens, out_t("varpart_sensitivity.csv"), row.names = FALSE)

# the log
for (gr in c(unique(grp_of))) for (lv in c("Among clusters", "Within clusters")) for (rn in unique(diffs$response)) {
  f <- fractions[fractions$group == gr & fractions$level == lv & fractions$response == rn & fractions$partition == "Environment vs space", ]; d <- diffs[diffs$group == gr & diffs$level == lv & diffs$response == rn, ]
  if (!nrow(f)) next; v <- setNames(f$value, f$fraction)
  message(sprintf("  %-7s | %s | %-16s environment only %5.1f%%, shared %5.1f%%, space only %5.1f%%, unexplained %5.1f%% | environment - space %s [%s, %s]", gr, lv, rn, 100 * v[["Environment only"]], 100 * v[["Shared"]], 100 * v[["Space only"]],
                  100 * (1 - max(0, v[["All sets together"]])), fmt(d$estimate[1]), fmt(d$lo[1]), fmt(d$hi[1])))
}
if (any(diffs$group == "Alpine - lowland")) for (i in which(diffs$group == "Alpine - lowland")) message(sprintf("  Alpine - lowland | %s | %-16s (environment - space) in alpine minus in lowland: %s [%s, %s]", diffs$level[i], diffs$response[i], fmt(diffs$estimate[i]), fmt(diffs$lo[i]), fmt(diffs$hi[i])))

# the diagnostics: the noise ceiling and the technical covariates (medians over the data sets of a group)
diag_tab <- bind_rows(lapply(names(G), function(nm) { d <- G[[nm]]$diag; if (is.null(d)) return(NULL); x <- d$ceiling; if (!is.null(d$tech)) x <- merge(x, d$tech, by = "response", all.x = TRUE, sort = FALSE)
  cbind(group = sub("\\|.*", "", nm), data_set = nm, x) }))
if (nrow(diag_tab)) {
  dsum <- diag_tab %>% group_by(group, response) %>% summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE)), .groups = "drop") %>% as.data.frame()
  dsum$explained_share_of_ceiling <- ifelse(is.finite(dsum$noise_ceiling) & dsum$noise_ceiling > 0, dsum$env_space / dsum$noise_ceiling, NA)
  write.csv(dsum, out_t("diagnostics_noise_ceiling_technical.csv"), row.names = FALSE)
  pc <- function(x) ifelse(is.finite(x), sprintf("%.1f%%", 100 * x), "NA")
  message("Diagnostics, among clusters (the share of the variation that is reproducible at all, and what the technical covariates take):")
  for (i in seq_len(nrow(dsum))) { d <- dsum[i, ]; message(sprintf("  %-7s | %-16s noise ceiling %s (from %.0f clusters with >= %d samples); environment + space explain %s (%s of the ceiling); read depth + unmapped share: alone %s, beyond environment + space %s",
    d$group, d$response, pc(d$noise_ceiling), d$clusters_used, 2 * CFG$among_n_per_cluster, pc(d$env_space), ifelse(is.finite(d$explained_share_of_ceiling), sprintf("%.0f%%", 100 * d$explained_share_of_ceiling), "NA"),
    if ("tech_alone" %in% names(d)) pc(d$tech_alone) else "not computed", if ("tech_beyond_env_space" %in% names(d)) pc(d$tech_beyond_env_space) else "not computed")) }
}

# RESPONSE SENSITIVITY (among clusters, composition; the first alpine draw and the first lowland replicate, with the same units and predictors as the analysis):
# for several sets of OTUs (the n most abundant; the OTUs present in at least x of the samples of the group) and three transformations (square root = Hellinger; relative
# abundance, which lets the abundant OTUs dominate; log(1 + 1000 p)), the NOISE CEILING (the share of the among-cluster variation that two independent sets of
# among_n_per_cluster samples of the same cluster have in common) and the variance that environment + space explain.
resp_sens <- NULL
if (isTRUE(CFG$response_sensitivity)) resp_sens <- tryCatch({
  m <- CFG$among_n_per_cluster; seed_keep <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  tfs <- list(`Hellinger (sqrt)` = function(M) { M <- as(M, "dgCMatrix"); M@x <- sqrt(M@x); M }, `Relative abundance` = function(M) as(M, "dgCMatrix"),
              `log(1 + 1000 p)` = function(M) { M <- as(M, "dgCMatrix"); M@x <- log1p(1000 * M@x); M })
  avg_mat <- function(cl_vec, rows_n) { u <- sort(unique(cl_vec)); Matrix::sparseMatrix(i = seq_len(rows_n), j = match(cl_vec, u), x = 1 / m, dims = c(rows_n, length(u))) }
  one_group <- function(gname, g) {
    u <- g$units; cols_g <- col_of[u$ix]; W <- avg_mat(u$cl_p, length(u$pick)); Eall <- do.call(cbind, g$among$L$E_sets); S <- g$among$L$S
    cl_tab <- table(u$cl); cl_2 <- names(cl_tab)[cl_tab >= 2 * m]
    two <- lapply(cl_2, function(cc) { sm <- sample(which(u$cl == cc)); list(sm[seq_len(m)], sm[m + seq_len(m)]) })
    r1 <- unlist(lapply(two, `[[`, 1)); r2 <- unlist(lapply(two, `[[`, 2)); W12 <- avg_mat(rep(cl_2, each = m), length(r1))
    prev_g <- Matrix::rowMeans(Prel[, cols_g, drop = FALSE] > 0)
    variant <- function(sel, tf_name, label, type) {
      tfn <- tfs[[tf_name]]; U <- as.matrix(Matrix::t(tfn(Prel[sel, cols_g[u$pick], drop = FALSE]) %*% W)); Y <- scale(U, scale = FALSE)
      p <- partition2(Eall, S, Y, sum(Y^2), nrow(Y) - 1); ceil <- NA_real_
      if (length(cl_2) >= 15) { U1 <- as.matrix(Matrix::t(tfn(Prel[sel, cols_g[r1], drop = FALSE]) %*% W12)); U2 <- as.matrix(Matrix::t(tfn(Prel[sel, cols_g[r2], drop = FALSE]) %*% W12))
        Tv <- (sum(apply(U1, 2, var)) + sum(apply(U2, 2, var))) / 2; Nv <- mean(rowSums((U1 - U2)^2)) / 2; ceil <- max(0, min(1, 1 - Nv / Tv)) }
      data.frame(group = gname, transform = tf_name, selection = label, type = type, n_otus = length(sel), noise_ceiling = ceil, explained = p[["explained"]],
                 pure_E = p[["pure_E"]], shared = p[["shared"]], pure_S = p[["pure_S"]], stringsAsFactors = FALSE)
    }
    bind_rows(lapply(names(tfs), function(tn) bind_rows(
      lapply(CFG$sens_otus, function(K) { sel <- rank_all[seq_len(min(K, length(rank_all)))]; variant(sel, tn, if (is.finite(K)) paste("top", K) else "all OTUs", "most abundant") }),
      lapply(CFG$sens_prev, function(x) { sel <- rank_all[prev_g[rank_all] >= x]; if (length(sel) < 10) NULL else variant(sel, tn, sprintf("present in >= %.0f%%", 100 * x), "prevalence") }),
      if (CFG$min_samples > 0) list({ sel <- rank_all[n_pres_all[rank_all] >= CFG$min_samples]; if (length(sel) < 10) NULL else variant(sel, tn, sprintf("present in >= %d samples (the filter)", CFG$min_samples), "min_samples") }))))
  }
  gs <- list(Alpine = G[["Alpine|draw1"]], Lowland = G[[which(grp_of == "Lowland")[1]]])
  out <- bind_rows(lapply(names(gs), function(gn) if (is.null(gs[[gn]])) NULL else one_group(gn, gs[[gn]])))
  if (!is.null(seed_keep)) assign(".Random.seed", seed_keep, envir = globalenv()); out
}, error = function(e) { message("Response sensitivity not computed: ", conditionMessage(e)); NULL })
if (!is.null(resp_sens) && nrow(resp_sens)) {
  write.csv(resp_sens, out_t("response_sensitivity.csv"), row.names = FALSE)
  pc1 <- function(x) ifelse(is.finite(x), sprintf("%.1f%%", 100 * x), "NA")
  message(sprintf("Response sensitivity (among clusters, composition; the analysis uses %d OTUs, Hellinger transformation; the filter variant is the one with >= %d samples):", length(top), CFG$min_samples))
  for (gn in unique(resp_sens$group)) for (tn in unique(resp_sens$transform)) { d <- resp_sens[resp_sens$group == gn & resp_sens$transform == tn, ]; b <- d[which.max(d$explained), ]; a <- d[d$selection == "all OTUs", ]
    message(sprintf("  %-7s | %-18s best: %-18s (%d OTUs) explains %s of a noise ceiling of %s%s", gn, tn, b$selection, b$n_otus, pc1(b$explained), pc1(b$noise_ceiling),
                    if (nrow(a)) sprintf("; all %d OTUs: %s of %s", a$n_otus[1], pc1(a$explained[1]), pc1(a$noise_ceiling[1])) else "")) }
  rs <- resp_sens %>% tidyr::pivot_longer(c(noise_ceiling, explained, pure_E, pure_S), names_to = "quantity", values_to = "value") %>% as.data.frame()
  rs$quantity <- factor(c(noise_ceiling = "Noise ceiling", explained = "Environment + space", pure_E = "Environment only", pure_S = "Space only")[rs$quantity], levels = c("Noise ceiling", "Environment + space", "Environment only", "Space only"))
  rs$transform <- factor(rs$transform, levels = names(tfs)); rs$value <- pmax(rs$value, 0)
  pR <- ggplot(rs, aes(n_otus, value, colour = quantity, linetype = quantity)) + geom_line(data = rs[rs$type != "min_samples", ], aes(group = interaction(quantity, type)), linewidth = 0.5, alpha = 0.6) + geom_point(aes(shape = type), size = 2) +
    scale_x_log10(labels = scales::comma) + scale_y_continuous(labels = scales::percent) +
    scale_colour_manual(values = c("Noise ceiling" = "grey40", "Environment + space" = "black", "Environment only" = CFG$cols[["Environment"]], "Space only" = CFG$cols[["Space"]]), name = NULL) +
    scale_linetype_manual(values = c("Noise ceiling" = "dashed", "Environment + space" = "solid", "Environment only" = "solid", "Space only" = "solid"), name = NULL) + scale_shape_manual(values = c("most abundant" = 16, prevalence = 17, min_samples = 15), name = "OTUs chosen by") +
    facet_grid(transform ~ group) + labs(x = "Number of OTUs in the composition", y = "Share of the among-cluster variation (adjusted R2 for the predictors)") + theme_fig + theme(legend.position = "top")
  cap_r <- "Response sensitivity. Among-cluster composition, with the same units and predictors as Figure 4, for sets of OTUs chosen by abundance (the n most abundant) or by prevalence (present in at least x of the samples of the group; the square = present in at least min_samples samples, the filter of the analysis) and three transformations. Dashed: the noise ceiling, the share of the variation that two independent sets of samples of the same cluster have in common (no predictor can explain more); solid: the variance explained by the environment and space together, by the environment only and by space only (adjusted R2, below 0 shown as 0). If the explained share rises when rare OTUs are removed, they are noise; the analysis uses the OTUs of top_otus and min_prevalence."
  cl_r <- strwrap(cap_r, 150); png(file.path(out_folder, "fig4_response_sensitivity.png"), width = 11, height = 9.5 + 0.17 * length(cl_r), units = "in", res = 300)
  grid.newpage(); pushViewport(viewport(layout = grid.layout(2, 1, heights = unit.c(unit(1, "null"), unit(0.17 * length(cl_r) + 0.2, "in"))))); print(pR, vp = viewport(layout.pos.row = 1))
  pushViewport(viewport(layout.pos.row = 2)); grid.text(paste(cl_r, collapse = "\n"), x = 0.01, just = "left", gp = gpar(fontsize = 9)); popViewport(2); dev.off()
}

## ---- 4. Figure --------------------------------------------------------------------------------------------------------------------------------------
grp_lev <- intersect(c("Alpine", "Lowland"), unique(fractions$group)); resp_lab <- unique(fractions$response); lv_lab <- intersect(c("Among clusters", "Within clusters"), unique(fractions$level))
prep_facets <- function(d) { d$group <- factor(d$group, levels = grp_lev); d$response <- factor(d$response, levels = resp_lab); d$level <- factor(d$level, levels = lv_lab); d }
fr2 <- fractions[fractions$partition == "Environment vs space" & fractions$fraction %in% c("Environment only", "Shared", "Space only", "All sets together"), ]
barA <- bind_rows(lapply(split(fr2, list(fr2$group, fr2$level, fr2$response), drop = TRUE), function(d) { v <- setNames(d$value, d$fraction); pos <- pmax(v[c("Environment only", "Shared", "Space only")], 0)
  data.frame(group = d$group[1], level = d$level[1], response = d$response[1], part = c("Environment only", "Shared", "Space only", "Unexplained"), value = c(pos, max(0, 1 - sum(pos))), stringsAsFactors = FALSE) }))
barA <- prep_facets(barA); barA$part <- factor(barA$part, levels = c("Unexplained", "Space only", "Shared", "Environment only"))
colA <- c("Environment only" = CFG$cols[["Environment"]], "Shared" = CFG$cols[["Shared"]], "Space only" = CFG$cols[["Space"]], "Unexplained" = CFG$cols[["Unexplained"]])
pA <- ggplot(barA, aes(group, value, fill = part)) + geom_col(width = 0.75, colour = "grey40", linewidth = 0.25) +
  geom_text(data = barA[barA$part != "Unexplained" & barA$value >= 0.03, ], aes(label = sprintf("%.0f%%", 100 * value)), position = position_stack(vjust = 0.5), size = 2.8, colour = "grey10") +
  scale_fill_manual(values = colA, breaks = names(colA), name = NULL) + facet_grid(level ~ response) + scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
  labs(tag = "A", x = NULL, y = "Share of the variation (adjusted R2)") + theme_fig + theme(legend.position = "top")
f4 <- fractions[fractions$partition == "Climate / soil / terrain / space", ]
cat4 <- function(fr) { sets <- strsplit(fr, " \\+ "); vapply(sets, function(s) if (length(s) == 1) paste(s, "only") else if ("Space" %in% s) "Shared: environment + space" else "Shared: among environment sets", "") }
f4$cat <- cat4(f4$fraction)
barB <- f4 %>% group_by(group, level, response, cat) %>% summarise(value = sum(pmax(value, 0)), .groups = "drop") %>% as.data.frame()
barB <- bind_rows(barB, barB %>% group_by(group, level, response) %>% summarise(cat = "Unexplained", value = max(0, 1 - sum(value)), .groups = "drop"))
catlev <- c("Unexplained", "Space only", "Shared: environment + space", "Shared: among environment sets", "Terrain only", "Soil only", "Climate only")
barB <- prep_facets(barB[barB$cat %in% catlev, ]); barB$cat <- factor(barB$cat, levels = catlev)
colB <- c("Climate only" = CFG$cols[["Climate"]], "Soil only" = CFG$cols[["Soil"]], "Terrain only" = CFG$cols[["Terrain"]], "Shared: among environment sets" = CFG$cols[["SharedEnv"]],
          "Shared: environment + space" = CFG$cols[["Shared"]], "Space only" = CFG$cols[["Space"]], "Unexplained" = CFG$cols[["Unexplained"]])
pB <- ggplot(barB, aes(group, value, fill = cat)) + geom_col(width = 0.75, colour = "grey40", linewidth = 0.25) +
  geom_text(data = barB[barB$cat != "Unexplained" & barB$value >= 0.03, ], aes(label = sprintf("%.0f%%", 100 * value)), position = position_stack(vjust = 0.5), size = 2.6, colour = "grey10") +
  scale_fill_manual(values = colB, breaks = rev(catlev), name = NULL) + facet_grid(level ~ response) + scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
  guides(fill = guide_legend(nrow = 2)) + labs(tag = "B", x = NULL, y = "Share of the variation (adjusted R2)") + theme_fig + theme(legend.position = "top", legend.text = element_text(size = 8))
dd <- diffs; dd$level <- factor(dd$level, levels = lv_lab); dd$response <- factor(dd$response, levels = rev(resp_lab)); dd$group <- factor(dd$group, levels = c("Lowland", "Alpine", "Alpine - lowland"))
pC <- ggplot(dd, aes(estimate, response, colour = group, shape = group)) + geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_linerange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.7), linewidth = 0.7) + geom_point(position = position_dodge(width = 0.7), size = 2.3) +
  scale_colour_manual(values = CFG$cols_group, breaks = c("Alpine", "Lowland", "Alpine - lowland"), name = NULL) + scale_shape_manual(values = c(Alpine = 16, Lowland = 16, `Alpine - lowland` = 18), breaks = c("Alpine", "Lowland", "Alpine - lowland"), name = NULL) +
  facet_wrap(~ level, scales = "free_x") + labs(tag = "C", x = "Environment - space (adjusted R2; > 0: the environment is the better predictor)", y = NULL) + theme_fig + theme(legend.position = "top")
sens$group <- factor(sens$group, levels = grp_lev)
pD <- ggplot(sens, aes(factor(k_space), factor(pcs_max), fill = difference)) + geom_tile(colour = "white") + geom_text(aes(label = ifelse(abs(difference) < 0.005, "0.00", sprintf("%+.2f", difference))), size = 3) +
  scale_fill_gradient2(low = "#E69F00", mid = "white", high = "#009E73", midpoint = 0, name = "Environment\n- space") + facet_wrap(~ group) +
  labs(tag = "D", x = "Spatial eigenvectors (PCNM)", y = "Principal components per environmental set (cap)") + theme_fig + theme(panel.grid = element_blank())
cap <- paste("Figure 4. Is the soil microbiome better predicted by the environment or by geography, in alpine and in lowland soils? Variation partitioning (adjusted R2) of the community composition (Hellinger-transformed relative abundances) and of the alpha diversity (richness, Shannon, Pielou evenness, freed of sequencing depth), with the environment (CHELSA climate, SoilGrids soil, DEM terrain: principal components of each set) and geography (spatial eigenvectors among clusters, a local polynomial within clusters) as predictors.",
             "Alpine: all alpine samples; lowland: the lowland replicates (samples evenly spread over 5-degree cells), each analysed on its own and summarised by the median (the alpine analysis is repeated over random draws of the samples that form the cluster means). The clusters (one study / place) are the units: AMONG clusters, one point per cluster (the mean of the same number of random samples in every cluster and group), the regional signal; WITHIN clusters, everything centred on the cluster mean, the local signal (only where clusters have enough samples).",
             "A: environment only, shared and space only; B: the environment split into climate, soil and terrain (shared fractions summed; fractions below 0 are shown as 0). C: the difference between the environment and space (adjusted R2 of the environment minus that of space; > 0: the environment is the better predictor) in each group, and the difference between the groups (alpine minus lowland; > 0: the environment is relatively better in alpine soils), with 95% intervals from bootstraps that resample whole clusters (and pick a replicate at random; the predictor sets are those of the full data; the spread of the bootstrap, centred on the estimate).",
             "D: the same difference for the composition among clusters with other numbers of components and spatial eigenvectors: the verdict should not depend on them. Environment and space are correlated (the shared fraction); the unexplained part is what neither explains, not necessarily noise. The groups differ in spatial extent and in the number of clusters (tables/sample_counts.csv).")
cap_lines <- strwrap(cap, 170); cap_h <- 0.17 * length(cap_lines) + 0.2
png(file.path(out_folder, "fig4_variance_partitioning.png"), width = 14, height = 17 + cap_h, units = "in", res = 300)
grid.newpage(); pushViewport(viewport(layout = grid.layout(4, 2, heights = unit.c(unit(1, "null"), unit(1, "null"), unit(0.9, "null"), unit(cap_h, "in")))))
print(pA, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2)); print(pB, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))
print(pC, vp = viewport(layout.pos.row = 3, layout.pos.col = 1)); print(pD, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 1:2)); grid.text(paste(cap_lines, collapse = "\n"), x = 0.01, just = "left", gp = gpar(fontsize = 9)); popViewport(2); dev.off()
saveRDS(list(A = pA, B = pB, C = pC, D = pD), file.path(out_folder, "fig4_panels.rds"))
message("Done. Outputs in ", out_folder)
writeLines(capture.output(sessionInfo()), file.path(out_folder, "sessionInfo.txt"))
file.create(file.path(out_folder, ".done"))   # flag file for Snakemake