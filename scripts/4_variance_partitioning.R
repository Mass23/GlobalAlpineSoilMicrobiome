# =============================================================================
# Script 4 - What explains the composition of the soil microbiome across clusters (studies / places)? Alpine, lowland, both
#
# ONE UNIT = ONE CLUSTER (cluster_id: samples of one BioProject or one place). In every random DRAW each cluster contributes ONE random sample:
#   Alpine : one sample of every alpine cluster;
#   Lowland: as many lowland clusters (from the lowland replicates rep_XX, spread evenly over 5-degree cells), one random sample each;
#   Both   : the two together.
# The composition is that of the GENERA (the OTUs aggregated to tax_rank with the taxonomy file; OTUs of chloroplasts, mitochondria and eukaryotes removed; genera present in too few
# samples / clusters removed): relative abundances among the kept genera.
#
# PART 1 - VARIATION PARTITIONING (panel A): the Hellinger composition against the sets Climate (CHELSA), Soil (SoilGrids), Terrain (DEM), Technical (log10 reads, share of unmapped reads)
#   and Space (a spherical spline of the coordinates, mgcv bs = "sos", k_space functions). Every environmental set is reduced to its first principal components (up to pcs_max; no use of
#   the response). Adjusted R2 (Ezekiel; the fractions are those of vegan::varpart: checked against it in the log), fractions averaged over the draws. Shown: the share of the EXPLAINED
#   variation that belongs to every set alone, and to two or more sets together ("shared"); the total adjusted R2 is printed above every bar.
# PART 2 - NMDS + ENVFIT (panel B): in every draw an NMDS (Bray-Curtis, 2 axes) of the units and vegan::envfit of EVERY variable (the r2 of the variable with the ordination: how well
#   it explains the position of the samples) with permutation p-values; summarised over the draws (median r2, 95% range, the share of draws with p < 0.05) for Alpine, Lowland and Both.
#   The variables: all the cleaned CHELSA / SoilGrids / terrain variables, the technical ones, latitude and longitude. Panel B shows the top envfit_top_per_set of every set.
#
# SI FIGURE (figS_envfit_shape.png): for every variable of panel B the position of the samples along the direction that envfit fitted, against the variable, with a GAM smooth per group: envfit
#   assumes a straight line, the figure shows where it is curved; tables/envfit_shape.csv has the R2 of the line and of the GAM.
#
# Variables: more than max_na missing values: dropped; the others imputed with the median; right-skewed (skewness > skew_log, no negative values): log1p; everything z-scored.
# The full version of the script (adjusted-R2 + cross-validated ridge, noise ceiling, OTU overlap, rank test, alpha metrics, within-cluster level) is 4_variance_partitioning_full.R.
#
# Input : sample_data_filtered.tsv, results/all_sampled_data.csv (climate + soil + terrain; rep_XX = the lowland replicates), the filtered BIOM, the OTU taxonomy file
# Usage : Rscript 4_variance_partitioning.R [sample_table] [all_sampled_data.csv] [filtered_biom] [out_folder]
# Output: <out_folder>/fig4_variance_partitioning_envfit.png, figS_envfit_shape.png (SI: the shape of the envfit relationships) and tables/vp_fractions_by_draw.csv, vp_summary.csv,
#         envfit_by_draw.csv, envfit_summary.csv, envfit_shape.csv, nmds_stress.csv
# Needs : R packages vegan, mgcv, rbiom, Matrix, dplyr, tidyr, ggplot2
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(Matrix)
  library(rbiom)
  library(vegan)
  library(mgcv)
  library(grid)
})

TEST_MODE <- FALSE   # TRUE: a tiny run (few draws, few permutations)

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
  id_table = "Sample_acc", id_env = "site", biom_id = "Accessions",
  group_col = "perialpine_alpine", alpine_value = "Alpine", lowland_value = "Lowland", cluster_col = "cluster_id",
  rare_depth = 1000,            # samples with fewer (mapped) reads are dropped
  lowland_reps = 20,            # the lowland replicates (rep_01 ... ) whose samples form the lowland pool
  n_draws = 10,                 # random draws (one sample per cluster); the VP is averaged and the NMDS / envfit summarised over them
  env_sets = list(
    Climate = c("^bio[0-9]{2}$", "^(fcf|fgd|scd)$", "^(clt|cmi|hurs|pet|pr|rsds|sfcWind|spei12|spi12|tas|tasmax|tasmin|vpd)$"),
    Soil    = c("_30cm_avg$"),
    Terrain = c("^elevation_m$", "^slope_deg$", "^aspect_northness$", "^aspect_eastness$", "^tpi_m$", "^roughness_m$")),
  use_technical = TRUE,         # log10 reads and the share of unmapped reads as a set of their own (they vary by study and method)
  max_na = 0.2, skew_log = 2, pca_var = 0.9, pcs_max = 6, k_space = 15,
  taxonomy_file = file.path(data_folder, "microbeatlas", "taxonomy_soilCurated_OTU97_nextMAPrelease.tsv"),   # NULL or missing: the OTUs are used
  tax_rank = "family",           # domain, phylum, class, order, family, genus; NULL: OTUs
  tax_fill = TRUE,              # an OTU whose lineage stops above tax_rank goes to an 'unclassified' taxon of its lowest known rank (FALSE: it is left out)
  tax_exclude = "chloroplast|mitochond|Eukaryota|Eukarya",   # OTUs whose lineage matches are removed (plant / organelle DNA); NULL: kept
  tax_id_col = NULL, tax_lineage_col = NULL, tax_skip = 0,
  min_samples = 5, min_clusters = 5,   # a taxon must be present in at least this many samples and clusters (of the analysed ones)
  nmds_trymax = 200, envfit_perm = 199, envfit_top_per_set = 6,
  cols = c(Climate = "#D55E00", Soil = "#8C6D31", Terrain = "#56B4E9", Technical = "#7B6FB0", Space = "#E69F00", Shared = "#9E9E9E",
           Geography = "#E69F00", Alpine = "#2C5F9E", Lowland = "#7A9A7E", Both = "#222222")
)
if (TEST_MODE) { out_folder <- paste0(out_folder, "_TEST"); CFG$n_draws <- 2; CFG$nmds_trymax <- 5; CFG$envfit_perm <- 49; CFG$lowland_reps <- 2 }
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
space_spline <- function(ll, k) {                                                                       # a spherical spline basis of the coordinates (longitude, latitude)
  ll <- ll + matrix(runif(length(ll), -1e-4, 1e-4), nrow(ll)); k <- max(4, min(k, nrow(ll) - 2))
  X <- tryCatch(mgcv::smoothCon(eval(substitute(mgcv::s(lat, lon, bs = "sos", k = KK), list(KK = k))), data = data.frame(lat = pmin(pmax(ll[, 2], -89.9), 89.9), lon = ll[, 1]), absorb.cons = TRUE)[[1]]$X, error = function(e) NULL)
  if (is.null(X)) return(NULL); X <- X[, apply(X, 2, sd) > 1e-8, drop = FALSE]; scale(X, scale = FALSE)
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
  list(Agg = Matrix::sparseMatrix(i = match(lab[assigned], taxa), j = which(assigned), x = 1, dims = c(length(taxa), length(ids)), dimnames = list(taxa, ids)), assigned = assigned, filled = filled, n_taxa = length(taxa))
}


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
Prel <- X %*% Matrix::Diagonal(x = 1 / tot); rm(X)
if (!is.null(tax) && !is.null(CFG$tax_rank)) {                                                           # the OTUs aggregated to a taxonomic rank
  ag <- tax_agg(tax, CFG$tax_rank, rownames(Prel)); mass <- Matrix::rowSums(Prel)
  message(sprintf("Aggregated to %s: %d taxa from %d OTUs; %.1f%% of the reads have a %s, %.1f%% %s; %.1f%% have no lineage and are left out", CFG$tax_rank, ag$n_taxa, sum(ag$assigned), 100 * sum(mass[ag$assigned & !ag$filled]) / sum(mass), CFG$tax_rank,
                  100 * sum(mass[ag$filled]) / sum(mass), if (isTRUE(CFG$tax_fill)) "are in 'unclassified' taxa of a higher rank" else "are left out", 100 * sum(mass[!ag$assigned]) / sum(mass)))
  Prel <- ag$Agg %*% Prel; rm(ag)
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
set_of <- c(unlist(lapply(names(Zs), function(nm) setNames(rep(nm, ncol(Zs[[nm]])), colnames(Zs[[nm]])))), if (!is.null(tech)) setNames(rep("Technical", 2), colnames(tech)), Latitude = "Geography", Longitude = "Geography")

## ---- 2. The random draws: one sample per cluster ------------------------------------------------------------------------------------------
alp_cl <- split(alp_ix, cl_all[alp_ix]); n_u <- length(alp_cl); pick1 <- function(i) i[sample.int(length(i), 1)]
draw_alpine <- function() unname(vapply(alp_cl, pick1, 1L))
if (has_low) { low_cl <- split(low_ix, cl_all[low_ix]); cen <- data.frame(lat = vapply(low_cl, function(i) mean(tab$Latitude[i]), 1), lon = vapply(low_cl, function(i) mean(tab$Longitude[i]), 1)); cell <- paste(floor(cen$lat / 5), floor(cen$lon / 5))
  draw_lowland <- function() { o <- sample(nrow(cen)); rk <- ave(seq_along(o), cell[o], FUN = seq_along); ord <- o[order(rk, match(cell[o], sample(unique(cell[o]))))][seq_len(min(n_u, nrow(cen)))]; unname(vapply(low_cl[ord], pick1, 1L)) } }
message(sprintf("Draws: %d; every draw: %d alpine clusters%s", CFG$n_draws, n_u, if (has_low) sprintf(" + %d lowland clusters (of %d), one random sample each", min(n_u, length(low_cl)), length(low_cl)) else ""))
unit_data <- function(idx) {
  r <- match(idx, an_ix); P <- Prel[taxa, col_of[idx], drop = FALSE]; s <- Matrix::colSums(P); ok <- s > 0; idx <- idx[ok]; r <- r[ok]
  P <- as.matrix(Matrix::t(P[, ok, drop = FALSE] %*% Matrix::Diagonal(x = 1 / s[ok])))                  # units x taxa: relative abundances among the kept taxa
  list(idx = idx, P = P, Z = lapply(Zs, function(M) M[r, , drop = FALSE]), tech = if (!is.null(tech)) tech[r, , drop = FALSE], ll = cbind(tab$Longitude[idx], tab$Latitude[idx]))
}

## ---- 3. Variation partitioning and NMDS + envfit of every draw ----------------------------------------------------------------------------------
vp_run <- function(u) {
  sets <- lapply(u$Z, env_pcs, cap = CFG$pcs_max); sets <- sets[!vapply(sets, is.null, TRUE)]
  if (!is.null(u$tech)) sets$Technical <- scale(u$tech, scale = FALSE)
  S <- space_spline(u$ll, CFG$k_space); if (!is.null(S)) sets$Space <- S
  Y <- scale(sqrt(u$P), scale = FALSE); cm <- commonality(sets, Y, sum(Y^2), nrow(Y) - 1)
  list(fractions = cm$fractions, total = cm$total, sets = sets, Y = Y)
}
nmds_run <- function(u) {
  ord <- tryCatch(suppressWarnings(vegan::metaMDS(vegan::vegdist(u$P, "bray"), k = 2, trymax = CFG$nmds_trymax, trace = 0, autotransform = FALSE)), error = function(e) NULL); if (is.null(ord)) return(NULL)
  V <- data.frame(do.call(cbind, u$Z), check.names = FALSE); if (!is.null(u$tech)) V <- cbind(V, as.data.frame(u$tech)); V$Latitude <- u$ll[, 2]; V$Longitude <- u$ll[, 1]
  V <- V[, vapply(V, function(x) sd(x) > 1e-8, TRUE), drop = FALSE]; ef <- vegan::envfit(ord, V, permutations = CFG$envfit_perm)
  list(env = data.frame(variable = names(ef$vectors$r), r2 = as.numeric(ef$vectors$r), p = as.numeric(ef$vectors$pvals), stress = ord$stress, stringsAsFactors = FALSE),
       scores = unname(vegan::scores(ord, display = "sites")[, 1:2, drop = FALSE]), V = V)                    # the scores and the variables of the units are kept for the SI figure
}
vp_all <- list(); ef_all <- list(); stress_all <- list(); shape_store <- list(); t0 <- Sys.time()
for (r in seq_len(CFG$n_draws)) {
  a <- draw_alpine(); ul <- list(Alpine = a); if (has_low) { l <- draw_lowland(); ul$Lowland <- l; ul$Both <- c(a, l) }
  for (g in names(ul)) {
    u <- unit_data(ul[[g]]); v <- vp_run(u)
    if (r == 1 && g == "Alpine") {                                                                       # the check against vegan::varpart (all the environmental sets against space)
      Eall <- do.call(cbind, v$sets[setdiff(names(v$sets), "Space")]); mine <- partition2(Eall, v$sets$Space, v$Y, sum(v$Y^2), nrow(v$Y) - 1); vpo <- vegan::varpart(v$Y, Eall, v$sets$Space)$part$indfract$Adj.R.squared
      message(sprintf("Check against vegan::varpart (first alpine draw: all environmental sets vs space): max |difference| of the fractions = %.1e", max(abs(c(mine[["pure_E"]] - vpo[1], mine[["pure_S"]] - vpo[2], mine[["shared"]] - vpo[3]))))) }
    vp_all[[length(vp_all) + 1]] <- cbind(group = g, draw = r, n_units = nrow(u$P), v$fractions, total = v$total)
    e <- nmds_run(u); if (!is.null(e)) { ef_all[[length(ef_all) + 1]] <- cbind(group = g, draw = r, e$env); stress_all[[length(stress_all) + 1]] <- data.frame(group = g, draw = r, stress = e$env$stress[1])
      shape_store[[length(shape_store) + 1]] <- list(group = g, draw = r, scores = e$scores, V = e$V) }
  }
  if (r == 1 || r %% 5 == 0 || r == CFG$n_draws) message(sprintf("  draw %d of %d done (%.1f min)", r, CFG$n_draws, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
vp <- bind_rows(vp_all); ef <- bind_rows(ef_all); stress <- bind_rows(stress_all)
write.csv(vp, out_t("vp_fractions_by_draw.csv"), row.names = FALSE); write.csv(ef, out_t("envfit_by_draw.csv"), row.names = FALSE); write.csv(stress, out_t("nmds_stress.csv"), row.names = FALSE)

## ---- 4. Summaries ---------------------------------------------------------------------------------------------------------------------------
grp_lev <- intersect(c("Alpine", "Lowland", "Both"), unique(vp$group))
vp_sum <- vp %>% group_by(group, fraction, n_sets) %>% summarise(value = mean(value), .groups = "drop") %>% as.data.frame()
tot_sum <- vp %>% group_by(group, draw) %>% summarise(total = first(total), .groups = "drop") %>% group_by(group) %>% summarise(R2 = mean(total), R2_lo = quantile(total, 0.025), R2_hi = quantile(total, 0.975), .groups = "drop") %>% as.data.frame()
vp_sum$category <- ifelse(vp_sum$n_sets == 1, paste(vp_sum$fraction, "only"), "Shared (two or more sets)")
cat_sum <- vp_sum %>% group_by(group, category) %>% summarise(value = sum(value), .groups = "drop") %>% as.data.frame()
cat_sum <- cat_sum %>% group_by(group) %>% mutate(share = pmax(value, 0) / sum(pmax(value, 0))) %>% ungroup() %>% as.data.frame()
write.csv(left_join(cat_sum, tot_sum, by = "group"), out_t("vp_summary.csv"), row.names = FALSE)
ef_sum <- ef %>% group_by(group, variable) %>% summarise(r2m = median(r2), lo = quantile(r2, 0.025), hi = quantile(r2, 0.975), sig = mean(p < 0.05), .groups = "drop") %>% as.data.frame()
ef_sum$set <- unname(set_of[ef_sum$variable]); write.csv(ef_sum, out_t("envfit_summary.csv"), row.names = FALSE)
for (g in grp_lev) { t_ <- tot_sum[tot_sum$group == g, ]; cs <- cat_sum[cat_sum$group == g, ]; cs <- cs[order(-cs$share), ]; e_ <- ef_sum[ef_sum$group == g, ]; e_ <- head(e_[order(-e_$r2m), ], 5)
  message(sprintf("  %-7s | total adjusted R2 %.1f%% [%.1f, %.1f] | share of the explained variation: %s | NMDS stress %.3f | top envfit r2: %s", g, 100 * t_$R2, 100 * t_$R2_lo, 100 * t_$R2_hi,
                  paste(sprintf("%s %.0f%%", sub(" only", "", cs$category[1:min(4, nrow(cs))]), 100 * cs$share[1:min(4, nrow(cs))]), collapse = ", "), mean(stress$stress[stress$group == g]), paste(sprintf("%s %.2f", e_$variable, e_$r2m), collapse = ", "))) }

## ---- 5. Figure ----------------------------------------------------------------------------------------------------------------------------------
catlev <- c("Shared (two or more sets)", "Space only", "Technical only", "Terrain only", "Soil only", "Climate only")
colA <- c("Climate only" = CFG$cols[["Climate"]], "Soil only" = CFG$cols[["Soil"]], "Terrain only" = CFG$cols[["Terrain"]], "Technical only" = CFG$cols[["Technical"]], "Space only" = CFG$cols[["Space"]], "Shared (two or more sets)" = CFG$cols[["Shared"]])
cs <- cat_sum[cat_sum$category %in% catlev, ]; cs$category <- factor(cs$category, levels = catlev); cs$group <- factor(cs$group, levels = grp_lev)
tt <- tot_sum; tt$group <- factor(tt$group, levels = grp_lev)
pA <- ggplot(cs, aes(group, share, fill = category)) + geom_col(width = 0.7, colour = "grey30", linewidth = 0.25) +
  geom_text(data = cs[cs$share >= 0.04, ], aes(label = sprintf("%.0f%%", 100 * share)), position = position_stack(vjust = 0.5), size = 3.2, colour = "grey10") +
  geom_text(data = tt, aes(group, 1.0, label = sprintf("R2 = %.0f%%", 100 * R2)), inherit.aes = FALSE, vjust = -0.5, size = 3.4, fontface = "bold") +
  scale_fill_manual(values = colA, breaks = rev(catlev), name = NULL) + scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.08))) +
  labs(tag = "A", x = NULL, y = "Share of the explained variation (adjusted R2)") + theme_fig + theme(legend.position = "bottom") + guides(fill = guide_legend(ncol = 2))
sel <- ef_sum %>% group_by(set, variable) %>% summarise(m = max(r2m), .groups = "drop") %>% group_by(set) %>% slice_max(m, n = CFG$envfit_top_per_set, with_ties = FALSE) %>% ungroup() %>% as.data.frame()
eb <- ef_sum[ef_sum$variable %in% sel$variable, ]; eb$group <- factor(eb$group, levels = grp_lev)
ord_v <- sel[order(sel$set, sel$m), ]; eb$variable <- factor(eb$variable, levels = ord_v$variable)
setlev <- intersect(c("Climate", "Soil", "Terrain", "Technical", "Geography"), unique(eb$set)); eb$set <- factor(eb$set, levels = setlev)
eb$signif <- factor(ifelse(eb$sig >= 0.5, "p < 0.05 in at least half of the draws", "otherwise"), levels = c("p < 0.05 in at least half of the draws", "otherwise"))
pB <- ggplot(eb, aes(r2m, variable, colour = group)) +
  geom_linerange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.7), linewidth = 0.5, alpha = 0.7) + geom_point(aes(shape = signif), position = position_dodge(width = 0.7), size = 2.2, stroke = 0.8) +
  scale_colour_manual(values = CFG$cols[grp_lev], name = NULL) + scale_shape_manual(values = c(16, 1), name = NULL) + facet_grid(set ~ ., scales = "free_y", space = "free_y") +
  labs(tag = "B", x = "envfit r2 with the NMDS (median and 95% range over the draws)", y = NULL) + theme_fig + theme(legend.position = "bottom", strip.text.y = element_text(angle = 0))
cap <- paste("Figure 4. What explains the composition of the soil microbiome across clusters? Units: clusters (studies / places), one random sample each, in", CFG$n_draws, "random draws; Alpine: the alpine clusters; Lowland: as many lowland clusters (from the lowland replicates, spread evenly over 5-degree cells); Both: the two together. The composition is that of the", unit_lab, "(chloroplast, mitochondria and eukaryote OTUs removed; present in at least", CFG$min_samples, "samples and", CFG$min_clusters, "clusters).",
             "A: variation partitioning (adjusted R2, mean over the draws) of the Hellinger composition among the sets Climate (CHELSA), Soil (SoilGrids), Terrain (DEM), Technical (read depth and unmapped share) and Space (spherical spline); every environmental set is reduced to its first principal components. The bars show the share of the explained variation that belongs to every set alone and to two or more sets together (negative fractions as 0); the total adjusted R2 is above every bar.",
             "B: NMDS (Bray-Curtis, 2 axes) of every draw and envfit of every variable: the r2 of the variable with the ordination (median and 95% range over the draws); the top", CFG$envfit_top_per_set, "variables of every set are shown, filled when the permutation p-value is below 0.05 in at least half of the draws. Environment and space are correlated: the shared fraction and the envfit r2 of correlated variables overlap.")
cl_ <- strwrap(cap, 170); n_rows <- nrow(sel)
png(file.path(out_folder, "fig4_variance_partitioning_envfit.png"), width = 14, height = max(8, 0.27 * n_rows + 3.5) + 0.17 * length(cl_), units = "in", res = 300)
grid.newpage(); pushViewport(viewport(layout = grid.layout(2, 2, widths = unit(c(1, 1.5), "null"), heights = unit.c(unit(1, "null"), unit(0.17 * length(cl_) + 0.2, "in")))))
print(pA, vp = viewport(layout.pos.row = 1, layout.pos.col = 1)); print(pB, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1:2)); grid.text(paste(cl_, collapse = "\n"), x = 0.01, just = "left", gp = gpar(fontsize = 9)); popViewport(2); dev.off()
saveRDS(list(A = pA, B = pB), file.path(out_folder, "fig4_panels.rds"))

## ---- 6. SI figure: the shape of the envfit relationships --------------------------------------------------------------------------------------
# For every variable of panel B and every draw: the direction in the NMDS that envfit fitted (the coefficients of the variable on the two axes) and the position of every unit along it
# (the projection of its scores, standardised); if the relationship were linear, that position would be a straight line of the variable. A flexible smooth (GAM) per group shows the shape.
shape_long <- bind_rows(lapply(shape_store, function(z) bind_rows(lapply(intersect(as.character(sel$variable), colnames(z$V)), function(v) {
  x <- z$V[[v]]; b <- coef(lm(x ~ z$scores[, 1] + z$scores[, 2]))[2:3]; if (anyNA(b) || sum(b^2) == 0) return(NULL)
  y <- as.numeric(z$scores %*% (b / sqrt(sum(b^2)))); data.frame(group = z$group, draw = z$draw, variable = v, x = x, y = as.numeric(scale(y))) }))))
if (nrow(shape_long)) {
  shape_tab <- bind_rows(lapply(split(shape_long, list(shape_long$variable, shape_long$group), drop = TRUE), function(d) {
    lin <- summary(lm(y ~ x, d))$r.squared; gm <- tryCatch(mgcv::gam(y ~ s(x, k = 5), data = d), error = function(e) NULL); if (is.null(gm)) return(NULL)
    data.frame(variable = d$variable[1], group = d$group[1], r2_linear = lin, r2_gam = summary(gm)$dev.expl, gain_gam_minus_linear = summary(gm)$dev.expl - lin, edf = sum(gm$edf[-1]), n_points = nrow(d)) }))
  shape_tab$set <- unname(set_of[shape_tab$variable]); write.csv(shape_tab, out_t("envfit_shape.csv"), row.names = FALSE)
  top_gain <- head(shape_tab[order(-shape_tab$gain_gam_minus_linear), ], 5)
  message("Shape of the envfit relationships (GAM minus linear R2 of the position along the envfit axis): the largest departures from linearity: ", paste(sprintf("%s (%s) %+.2f", top_gain$variable, top_gain$group, top_gain$gain_gam_minus_linear), collapse = "; "))
  shape_long$group <- factor(shape_long$group, levels = grp_lev); shape_long$variable <- factor(shape_long$variable, levels = rev(ord_v$variable))
  lab_v <- setNames(sprintf("%s  (%s)", ord_v$variable, ord_v$set), ord_v$variable)
  ann <- shape_tab[shape_tab$group %in% grp_lev, ]; ann$g1 <- substr(ann$group, 1, 1); ann <- ann %>% group_by(variable) %>% summarise(txt = paste0("gain: ", paste(sprintf("%s %+.2f", g1, gain_gam_minus_linear), collapse = ", ")), .groups = "drop")
  ann$variable <- factor(ann$variable, levels = levels(shape_long$variable))
  pS <- ggplot(shape_long, aes(x, y)) +
    geom_point(data = shape_long[shape_long$group != "Both", ], aes(colour = group), alpha = 0.10, size = 0.4) +
    geom_smooth(aes(colour = group), method = "gam", formula = y ~ s(x, k = 5), se = FALSE, linewidth = 0.9) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = txt), inherit.aes = FALSE, hjust = -0.05, vjust = 1.4, size = 2.3, colour = "grey25") +
    facet_wrap(~ variable, scales = "free_x", ncol = 6, labeller = labeller(variable = lab_v)) + coord_cartesian(ylim = c(-2.6, 2.6)) +
    scale_colour_manual(values = CFG$cols[grp_lev], name = NULL) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1.5))) +
    labs(x = "Variable (z-score; latitude / longitude in degrees)", y = "Position of the sample along the envfit axis of the variable (standardised)") + theme_fig +
    theme(legend.position = "top", strip.text = element_text(size = 7.5, face = "plain"), axis.text = element_text(size = 7))
  cap_s <- paste("Supplementary Figure. The shape of the envfit relationships of panel B. For every variable and every draw, the direction in the NMDS that envfit fitted for the variable (the coefficients of the variable on the two NMDS axes) was found and every sample was projected onto it; the projection (standardised) is plotted against the variable, for all draws together (points: Alpine and Lowland; the Both ordination has its own axes and is shown by its line only).",
                 "A line is a flexible smooth (GAM, 5 basis functions) per group; envfit assumes a straight line. 'gain' = the R2 of the GAM minus the R2 of a straight line (A alpine, L lowland, B both; tables/envfit_shape.csv): a large gain means that the relationship is curved (for example an optimum) and that envfit r2 understates it; the draws share clusters, so the points are not independent and no test is given.")
  cl_s <- strwrap(cap_s, 200); nr <- ceiling(length(unique(shape_long$variable)) / 6)
  png(file.path(out_folder, "figS_envfit_shape.png"), width = 16, height = 2.5 * nr + 1.6 + 0.17 * length(cl_s), units = "in", res = 300)
  grid.newpage(); pushViewport(viewport(layout = grid.layout(2, 1, heights = unit.c(unit(1, "null"), unit(0.17 * length(cl_s) + 0.2, "in"))))); print(pS, vp = viewport(layout.pos.row = 1))
  pushViewport(viewport(layout.pos.row = 2)); grid.text(paste(cl_s, collapse = "\n"), x = 0.01, just = "left", gp = gpar(fontsize = 9)); popViewport(2); dev.off()
}
message("Done. Outputs in ", out_folder)
writeLines(capture.output(sessionInfo()), file.path(out_folder, "sessionInfo.txt"))
file.create(file.path(out_folder, ".done"))   # flag file for Snakemake