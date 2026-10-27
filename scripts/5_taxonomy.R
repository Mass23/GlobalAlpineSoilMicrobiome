# =============================================================================
# Script 5 - Which families are alpine? Taxonomy (family level) and differential abundance, alpine vs lowland
#
# Samples : all alpine samples and the lowland samples of the lowland replicates (rep_XX of all_sampled_data.csv), at most max_per_cluster random samples of every cluster (so that big studies do
#           not dominate and the fit stays light); the prevalence is computed on ALL the alpine samples.
# Families: the OTUs aggregated to tax_rank = family with the taxonomy file (OTUs of chloroplasts, mitochondria and eukaryotes removed; an OTU whose lineage stops above family goes to an
#           'unclassified' taxon, which stays in the CLR but is not shown). Families present in at least min_prev of the alpine OR of the lowland samples are tested.
# Abundance: the centred log-ratio (CLR) of the family counts per sample, log2 scale: log2((count + clr_pseudocount) / geometric mean of the families of the sample). A difference of CLR values is
#           a log2 fold change in the abundance relative to the rest of the community (compositional data: it is not an absolute change).
# Model   : one GAM per family (mgcv::bam, fREML, discrete): CLR ~ group + s(log10 reads) + s(share of unmapped reads) + s(latitude, longitude; spherical spline, k_space) + s(cluster; random effect).
#           The coefficient of group (Alpine - Lowland) is the log2 fold change accounting for sequencing depth, unmapped reads, geography (the contrast is made at the same place) and the study.
#           p-values from the t test of the coefficient, Benjamini-Hochberg FDR across the families.
# Figure  : A the taxonomic tree (domain > phylum > class > order > family) with the prevalence of every family in the alpine samples (bar; the grey point: in the lowland samples); B the log2 fold
#           change (95% interval) of every family, aligned with the tree; filled: FDR < fdr; C the link with the ENVFIT results of script 4: the Spearman correlation of the CLR abundance of every
#           family with the c_n_vars variables of the highest envfit r2 (tables/envfit_summary.csv of script 4, group c_group), across the alpine CLUSTERS (cluster means, so a study counts once);
#           dots: FDR < fdr (per variable, across the families); the numbers on top are the envfit r2 of the variables. Without that file panel C is left out.
#
# SUPPLEMENTARY FIGURE: figS5_family_tree_differential_abundance_all_pageNN_of_KK.png: ALL the tested families in taxonomic order, at most si_per_page per page (a phylum is kept on one page when
#           it fits), same palette and axes as the main figure; tables/family_figure_pages.csv says on which page every family is.
#
# FIGURE SIZE: with hundreds of families the figure shows fig_max_families of them, chosen by fig_select ("both": half the most prevalent in the alpine samples, half the largest fold changes, the
#           significant ones first); the tree stays taxonomic; the tables keep all families (column shown_in_figure).
#
# Input : sample_data_filtered.tsv, results/all_sampled_data.csv (only its rep_XX columns), the filtered BIOM, the OTU taxonomy file
# Usage : Rscript 5_family_differential_abundance.R [sample_table] [all_sampled_data.csv] [filtered_biom] [out_folder] [envfit_summary.csv of script 4]
# Output: <out_folder>/fig5_family_tree_differential_abundance.png and tables/family_prevalence.csv, family_differential_abundance.csv, family_envfit_correlation.csv, da_samples.csv
# Needs : R packages mgcv, rbiom, Matrix, dplyr, tidyr, ggplot2 (parallel is part of R)
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

args <- commandArgs(trailingOnly = TRUE)
data_folder     <- "/home/renku/work/epfl-altshuler/GlobalAlpineSoilMicrobiome/data"
results_default <- file.path(dirname(data_folder), "results")
sample_table  <- if (length(args) >= 1) as.character(args[1]) else file.path(data_folder, "sample_data_filtered.tsv")
env_table     <- if (length(args) >= 2) as.character(args[2]) else file.path(results_default, "all_sampled_data.csv")
filtered_biom <- if (length(args) >= 3) as.character(args[3]) else
  file.path(data_folder, "microbeatlas", "otuTable_soilCurated_OTU97_nextMAPrelease_filtered.biom")
out_folder    <- if (length(args) >= 4) as.character(args[4]) else file.path(results_default, "5_family_differential_abundance")
envfit_file   <- if (length(args) >= 5) as.character(args[5]) else file.path(results_default, "4_variance_partitioning", "tables", "envfit_summary.csv")

CFG <- list(
  sample_file = sample_table, env_file = env_table, biom_file = filtered_biom,
  id_table = "Sample_acc", id_env = "site", biom_id = "Accessions",
  group_col = "perialpine_alpine", alpine_value = "Alpine", lowland_value = "Lowland", cluster_col = "cluster_id",
  rare_depth = 1000,            # samples with fewer (mapped) reads are dropped
  lowland_reps = 5,             # the lowland replicates (rep_01 ...) whose samples form the lowland pool
  max_per_cluster = 10,         # at most this many random samples of every cluster in the models
  taxonomy_file = file.path(data_folder, "microbeatlas", "taxonomy_soilCurated_OTU97_nextMAPrelease.tsv"),
  tax_rank = "family", tax_fill = TRUE, tax_exclude = "chloroplast|mitochond|Eukaryota|Eukarya", tax_id_col = NULL, tax_lineage_col = NULL, tax_skip = 0,
  min_prev = 0.10,              # a family is tested if it is present in at least this share of the alpine OR of the lowland samples
  clr_pseudocount = 0.5,        # added to the counts before the log
  k_space = 30,                 # basis functions of the spherical spline of the coordinates
  n_cores = 4,                  # families fitted in parallel
  fdr = 0.05,
  envfit_file = envfit_file,    # tables/envfit_summary.csv of script 4 (panel C); a missing file: no panel C
  c_group = "Alpine", c_n_vars = 10,   # panel C: the variables with the highest envfit r2 in this group of script 4, at most this many
  fig_max_families = 60,        # the figure shows at most this many families (the tables have all of them); NULL / Inf: all
  si_all_families = TRUE,       # supplementary figure(s) with ALL the tested families (in pages of si_per_page), when the main figure shows fewer
  si_per_page = 70,
  fig_select = "both",          # which: "prevalence" (the most prevalent in the alpine samples), "effect" (the largest |fold change|, the significant first) or "both" (half of each)
  min_tips_label = 4,           # phyla with at least this many families in the figure get a name
  cols = c(Alpine = "#2C5F9E", Lowland = "#7A9A7E", ns = "grey55")
)
if (TEST_MODE) { out_folder <- paste0(out_folder, "_TEST"); CFG$lowland_reps <- 2; CFG$n_cores <- 2; CFG$k_space <- 15 }
dir.create(file.path(out_folder, "tables"), recursive = TRUE, showWarnings = FALSE)
out_t <- function(...) file.path(out_folder, "tables", ...)
set.seed(1)

theme_fig <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), plot.tag = element_text(face = "bold", size = 16))
to_sparse <- function(m) {
  if (inherits(m, "simple_triplet_matrix")) Matrix::sparseMatrix(i = m$i, j = m$j, x = m$v, dims = c(m$nrow, m$ncol), dimnames = list(m$dimnames[[1]], m$dimnames[[2]]))
  else as(m, "dgCMatrix")
}
cluster_ids <- function(t_g) { cl <- as.character(t_g[[CFG$cluster_col]]); no <- is.na(cl) | !nzchar(cl); cl[no] <- paste0("none", which(no)); cl }
yes <- function(x) tolower(trimws(as.character(x))) %in% c("yes", "true", "1", "y")

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


## ---- 1. Data ------------------------------------------------------------------------------------------------------------------------------
raw <- read.delim(CFG$sample_file, quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE, fill = FALSE)
need <- c(CFG$id_table, CFG$biom_id, CFG$group_col, CFG$cluster_col, "Latitude", "Longitude", "alpha_total_n_reads", "unmapped_raw_counts", "total_raw_counts")
miss <- setdiff(need, names(raw)); if (length(miss)) stop("Missing columns in ", CFG$sample_file, ": ", paste(miss, collapse = ", "))
tab <- raw[raw[[CFG$group_col]] %in% c(CFG$alpine_value, CFG$lowland_value), ]
tab <- tab[is.finite(tab$Latitude) & is.finite(tab$Longitude) & is.finite(tab$alpha_total_n_reads) & tab$alpha_total_n_reads >= CFG$rare_depth, ]
env <- read.csv(CFG$env_file, check.names = FALSE, stringsAsFactors = FALSE); env <- env[!duplicated(env[[CFG$id_env]]), ]
mm <- match(tab[[CFG$id_table]], env[[CFG$id_env]]); tab <- tab[!is.na(mm), ]; env <- env[mm[!is.na(mm)], , drop = FALSE]
is_alp <- tab[[CFG$group_col]] == CFG$alpine_value
rep_cols <- head(sort(grep("^rep_[0-9]+$", names(env), value = TRUE)), CFG$lowland_reps)
alp_ix <- which(is_alp); low_ix <- sort(unique(unlist(lapply(rep_cols, function(rc) which(!is_alp & yes(env[[rc]]))))))
message(sprintf("%d alpine samples and %d lowland samples (the union of %d replicates)", length(alp_ix), length(low_ix), length(rep_cols)))
if (length(alp_ix) < 50 || length(low_ix) < 50) stop("fewer than 50 alpine or lowland samples: check the replicate columns rep_XX of ", basename(CFG$env_file))

message("Reading the BIOM: ", CFG$biom_file)
biom <- rbiom::read_biom(CFG$biom_file); ids <- intersect(tab[[CFG$biom_id]][c(alp_ix, low_ix)], biom$samples); biom <- biom[ids]
X <- to_sparse(biom$counts); X <- X[rownames(X) != "Unmapped", , drop = FALSE]; rm(biom); invisible(gc())
col_of <- match(tab[[CFG$biom_id]], colnames(X))
if (is.null(CFG$taxonomy_file) || !file.exists(CFG$taxonomy_file)) stop("this script needs the taxonomy file (taxonomy_file): ", CFG$taxonomy_file)
tax <- read_taxonomy(CFG$taxonomy_file, rownames(X)); if (is.null(tax)) stop("the taxonomy file could not be matched to the BIOM (see the lines above)")
if (!is.null(CFG$tax_exclude) && nzchar(CFG$tax_exclude)) {
  ex <- rownames(X) %in% tax$ids[grepl(CFG$tax_exclude, tax$excl, ignore.case = TRUE)]
  for (gp in c("Alpine", "Lowland")) { cols <- unique(col_of[if (gp == "Alpine") alp_ix else low_ix]); cols <- cols[!is.na(cols)]
    v <- Matrix::colSums(X[ex, cols, drop = FALSE]) / pmax(Matrix::colSums(X[, cols, drop = FALSE]), 1)
    message(sprintf("  chloroplast / mitochondria / eukaryote OTUs (%d): %s samples: median %.1f%% of the reads", sum(ex), gp, 100 * median(v))) }
  X <- X[!ex, , drop = FALSE]
}
tot <- Matrix::colSums(X); good <- !is.na(col_of) & tot[pmax(col_of, 1)] >= CFG$rare_depth
alp_ix <- alp_ix[good[alp_ix]]; low_ix <- low_ix[good[low_ix]]; sel_ix <- c(alp_ix, low_ix)
ag <- tax_agg(tax, CFG$tax_rank, rownames(X)); mass <- Matrix::rowSums(X)
message(sprintf("Aggregated to %s: %d taxa from %d OTUs; %.1f%% of the reads have a %s, %.1f%% are in 'unclassified' taxa", CFG$tax_rank, ag$n_taxa, sum(ag$assigned), 100 * sum(mass[ag$assigned & !ag$filled]) / sum(mass), CFG$tax_rank, 100 * sum(mass[ag$filled]) / sum(mass)))
Cfam <- ag$Agg %*% X[, col_of[sel_ix], drop = FALSE]; rm(X)                                              # families x samples: counts
is_a <- seq_along(sel_ix) <= length(alp_ix)
prev_alp <- Matrix::rowMeans(Cfam[, is_a, drop = FALSE] > 0); prev_low <- Matrix::rowMeans(Cfam[, !is_a, drop = FALSE] > 0)
keep <- prev_alp >= CFG$min_prev | prev_low >= CFG$min_prev
message(sprintf("%d of %d %s taxa are present in at least %.0f%% of the alpine or of the lowland samples", sum(keep), nrow(Cfam), CFG$tax_rank, 100 * CFG$min_prev))
if (sum(keep) < 5) stop("fewer than 5 taxa pass min_prev")

# the centred log-ratio (log2) of every sample over the kept taxa
L <- log(as.matrix(Cfam[keep, , drop = FALSE]) + CFG$clr_pseudocount); clr <- sweep(L, 2, colMeans(L)) / log(2); rm(L)
disp <- !grepl("unclassified", rownames(clr), fixed = TRUE)                                              # the 'unclassified' taxa stay in the CLR but are not tested or shown

## ---- 2. The samples of the models and the GAM of every family ----------------------------------------------------------------------------------
cl_all <- cluster_ids(tab); cap_pick <- function(ix) unlist(lapply(split(ix, cl_all[ix]), function(i) if (length(i) > CFG$max_per_cluster) sample(i, CFG$max_per_cluster) else i), use.names = FALSE)
fit_ix <- c(cap_pick(alp_ix), cap_pick(low_ix)); pos <- match(fit_ix, sel_ix)
us <- tab$unmapped_raw_counts[fit_ix] / tab$total_raw_counts[fit_ix]; us[!is.finite(us)] <- median(us[is.finite(us)]); us <- pmin(pmax(us, 0), 1)
dat <- data.frame(group = factor(ifelse(is_alp[fit_ix], "Alpine", "Lowland"), levels = c("Lowland", "Alpine")), log10_reads = log10(tab$alpha_total_n_reads[fit_ix]), unmapped_share = us,
                  lat = tab$Latitude[fit_ix], lon = tab$Longitude[fit_ix], cluster = factor(cl_all[fit_ix]))
da_samples <- data.frame(group = c("Alpine", "Lowland"), samples = as.integer(table(dat$group)[c("Alpine", "Lowland")]), clusters = tapply(as.character(dat$cluster), dat$group, function(z) length(unique(z)))[c("Alpine", "Lowland")])
write.csv(da_samples, out_t("da_samples.csv"), row.names = FALSE)
message(sprintf("Models: %d alpine samples (%d clusters) and %d lowland samples (%d clusters), at most %d per cluster", da_samples$samples[1], da_samples$clusters[1], da_samples$samples[2], da_samples$clusters[2], CFG$max_per_cluster))
k_sp <- max(5, min(CFG$k_space, nrow(unique(dat[, c("lat", "lon")])) - 2))
terms <- c("group", "s(log10_reads, k = 4)", if (sd(dat$unmapped_share) > 1e-6) "s(unmapped_share, k = 4)", sprintf("s(lat, lon, bs = \"sos\", k = %d)", k_sp), if (nlevels(dat$cluster) > 1) "s(cluster, bs = \"re\")")
form <- as.formula(paste("y ~", paste(terms, collapse = " + "))); message("Model: ", deparse(form, width.cutoff = 200))
fit_one <- function(f) {
  d <- dat; d$y <- clr[f, pos]
  m <- tryCatch(mgcv::bam(form, data = d, method = "fREML", discrete = TRUE, nthreads = 1), error = function(e) NULL)
  if (is.null(m)) return(NULL); pt <- summary(m)$p.table; if (!"groupAlpine" %in% rownames(pt)) return(NULL)
  data.frame(family = f, lfc = pt["groupAlpine", 1], se = pt["groupAlpine", 2], p = pt["groupAlpine", 4], stringsAsFactors = FALSE)
}
todo <- rownames(clr)[disp]; message(sprintf("Fitting %d GAMs (%d cores) ...", length(todo), CFG$n_cores)); t0 <- Sys.time()
res <- bind_rows(parallel::mclapply(todo, fit_one, mc.cores = CFG$n_cores))
message(sprintf("  done in %.1f min; %d of %d fits succeeded", as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(res), length(todo)))
res$lo <- res$lfc - 1.96 * res$se; res$hi <- res$lfc + 1.96 * res$se; res$fdr <- p.adjust(res$p, "BH")

## ---- 3. The taxonomy of the families and the tables ---------------------------------------------------------------------------------------------
rk_names <- c("domain", "phylum", "class", "order", "family"); rkc <- paste0("r_", rk_names)                       # (r_ columns: "family" is also the name of the label column)
fam <- data.frame(family = rownames(clr), prevalence_alpine = as.numeric(prev_alp[keep]), prevalence_lowland = as.numeric(prev_low[keep]), stringsAsFactors = FALSE)
R5 <- ag$R[match(fam$family, ag$lab), seq_len(5), drop = FALSE]; R5[is.na(R5)] <- ""
for (k in 1:5) fam[[rkc[k]]] <- R5[, k]
fam$name <- ifelse(nzchar(R5[, 5]), R5[, 5], sub(".*;", "", fam$family))
fam <- inner_join(fam, res, by = "family"); fam <- fam[order(fam$r_domain, fam$r_phylum, fam$r_class, fam$r_order, fam$name), ]
fam$sig <- fam$fdr < CFG$fdr
X_fig <- if (is.null(CFG$fig_max_families) || !is.finite(CFG$fig_max_families)) nrow(fam) else CFG$fig_max_families
if (nrow(fam) > X_fig) {
  by_prev <- order(-fam$prevalence_alpine, -abs(fam$lfc)); by_eff <- order(-fam$sig, -abs(fam$lfc))                    # the significant families first, by |fold change|
  idx <- switch(CFG$fig_select, prevalence = by_prev[seq_len(X_fig)], effect = by_eff[seq_len(X_fig)], both = head(unique(c(rbind(by_prev[seq_len(ceiling(X_fig / 2))], by_eff[seq_len(ceiling(X_fig / 2))]), by_eff, by_prev)), X_fig), stop("fig_select must be prevalence, effect or both"))
  fam$shown_in_figure <- seq_len(nrow(fam)) %in% idx
} else fam$shown_in_figure <- TRUE
message(sprintf("Figure: %d of the %d tested families are shown (%s)", sum(fam$shown_in_figure), nrow(fam), if (all(fam$shown_in_figure)) "all" else switch(CFG$fig_select, prevalence = "the most prevalent in the alpine samples", effect = "the largest fold changes, the significant first", both = "half the most prevalent in the alpine samples, half the largest fold changes")))
write.csv(setNames(fam[, c("family", rkc[1:4], "name", "prevalence_alpine", "prevalence_lowland", "lfc", "se", "lo", "hi", "p", "fdr", "shown_in_figure")], c("taxon", rk_names[1:4], "family", "prevalence_alpine", "prevalence_lowland", "log2_fold_change_alpine_minus_lowland", "se", "ci_lo", "ci_hi", "p", "fdr", "shown_in_figure")), out_t("family_differential_abundance.csv"), row.names = FALSE)
write.csv(data.frame(family = rownames(Cfam), prevalence_alpine = as.numeric(prev_alp), prevalence_lowland = as.numeric(prev_low)), out_t("family_prevalence.csv"), row.names = FALSE)
up <- fam[fam$sig & fam$lfc > 0, ]; dn <- fam[fam$sig & fam$lfc < 0, ]
message(sprintf("%d families tested; FDR < %.2f: %d enriched in the alpine soils, %d enriched in the lowland soils", nrow(fam), CFG$fdr, nrow(up), nrow(dn)))
top5 <- function(d, dec) paste(sprintf("%s %+.2f", head(d$name[order(if (dec) -d$lfc else d$lfc)], 5), head(sort(d$lfc, decreasing = dec), 5)), collapse = ", ")
if (nrow(up)) message("  most alpine: ", top5(up, TRUE)); if (nrow(dn)) message("  most lowland: ", top5(dn, FALSE))

## ---- 3b. The link with the envfit results of script 4 -------------------------------------------------------------------------------------------
cor_long <- NULL; vars_c <- NULL
if (!is.null(CFG$envfit_file) && file.exists(CFG$envfit_file)) {
  efs <- read.csv(CFG$envfit_file, stringsAsFactors = FALSE); efs <- efs[efs$group == CFG$c_group, ]; efs <- head(efs[order(-efs$r2m), ], CFG$c_n_vars)
  get_var <- function(v) if (v %in% names(env)) suppressWarnings(as.numeric(env[[v]])) else if (v == "Latitude") tab$Latitude else if (v == "Longitude") tab$Longitude else if (v == "log10_reads") log10(tab$alpha_total_n_reads) else
    if (v == "unmapped_share") tab$unmapped_raw_counts / tab$total_raw_counts else NULL
  efs <- efs[vapply(efs$variable, function(v) !is.null(get_var(v)), TRUE), ]
  if (nrow(efs) >= 2) {
    cl_a <- cl_all[alp_ix]; ucl <- sort(unique(cl_a)); n_cl <- as.numeric(table(factor(cl_a, levels = ucl)))
    cm <- t(rowsum(t(clr[fam$family, is_a, drop = FALSE]), cl_a)) / rep(n_cl, each = nrow(fam))                      # families x alpine clusters: the mean CLR of the cluster
    cor_long <- bind_rows(lapply(seq_len(nrow(efs)), function(i) {
      v <- efs$variable[i]; vm <- as.numeric(tapply(get_var(v)[alp_ix], factor(cl_a, levels = ucl), function(z) if (any(is.finite(z))) mean(z[is.finite(z)]) else NA_real_)); ok <- is.finite(vm)
      rho <- as.numeric(cor(t(cm[, ok, drop = FALSE]), vm[ok], method = "spearman")); nn <- sum(ok); tt <- rho * sqrt((nn - 2) / pmax(1 - rho^2, 1e-12)); p <- 2 * pt(-abs(tt), nn - 2)
      data.frame(family = fam$family, variable = v, envfit_r2 = efs$r2m[i], rho = rho, p = p, fdr = p.adjust(p, "BH"), n_clusters = nn, stringsAsFactors = FALSE) }))
    vars_c <- efs$variable; write.csv(cor_long, out_t("family_envfit_correlation.csv"), row.names = FALSE)
    message(sprintf("Panel C: correlation of the families with the %d variables of the highest envfit r2 (%s): %s; %d family-variable correlations with FDR < %.2f", length(vars_c), CFG$c_group, paste(sprintf("%s (r2 %.2f)", efs$variable, efs$r2m), collapse = ", "), sum(cor_long$fdr < CFG$fdr, na.rm = TRUE), CFG$fdr))
  }
} else message("No envfit file (", CFG$envfit_file, "): panel C is left out")

## ---- 4. Figures: the tree, the prevalence, the fold change and the link with envfit -----------------------------------------------------------------
# draw_fig() draws the four panels for a set of families (taxonomic order): the main figure (at most fig_max_families) and the supplementary pages (all the families), with the same palette and
# the same axis ranges, so that the pages can be compared.
fam_all <- fam; n_tested <- nrow(fam_all)
phy_all <- unique(ifelse(nzchar(fam_all$r_phylum), fam_all$r_phylum, "(no phylum)")); pal <- setNames(grDevices::hcl.colors(max(3, length(phy_all)), "Dark 3")[seq_along(phy_all)], phy_all)
xlim_B <- range(c(fam_all$lo, fam_all$hi), na.rm = TRUE)                                                   # the same fold-change axis on every page
draw_fig <- function(fam, file, cap_head, cap_sel, rds = NULL) {
  n <- nrow(fam); fam$y <- n:1
  for (k in 1:4) { pk <- fam[[rkc[k]]]; fam[[paste0("p", k)]] <- ifelse(nzchar(pk), pk, paste0("(no ", rk_names[k], ")")) }
  for (k in 1:4) fam[[paste0("k", k)]] <- do.call(paste, c(lapply(1:k, function(j) fam[[paste0("p", j)]]), sep = "|"))
  fam$phy <- fam$p2
  vseg <- list(); hseg <- list(); child <- data.frame(key = fam$family, y = fam$y, parent = fam$k4, phy = fam$phy, stringsAsFactors = FALSE)
  for (d in 4:1) {
    nd <- child %>% group_by(parent) %>% summarise(ymin = min(y), ymax = max(y), yn = mean(y), phy = first(phy), .groups = "drop")
    vseg[[d]] <- data.frame(x = d, ymin = nd$ymin, ymax = nd$ymax, phy = nd$phy); hseg[[d]] <- data.frame(x = d, xend = d + 1, y = child$y, phy = child$phy)
    if (d > 1) { child <- data.frame(key = nd$parent, y = nd$yn, parent = fam[[paste0("k", d - 1)]][match(nd$parent, fam[[paste0("k", d)]])], phy = nd$phy, stringsAsFactors = FALSE) } else root <- nd
  }
  vseg <- bind_rows(vseg); hseg <- bind_rows(hseg); stub <- data.frame(x = 0.4, xend = 1, y = root$yn, phy = root$phy)
  blocks <- fam %>% group_by(phy) %>% summarise(ymin = min(y) - 0.5, ymax = max(y) + 0.5, tips = n(), .groups = "drop") %>% arrange(desc(ymax)); blocks$band <- rep(c(TRUE, FALSE), length.out = nrow(blocks))
  band_layer <- geom_rect(data = blocks[blocks$band, ], aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = "grey94")
  top_pad <- if (!is.null(cor_long)) 8 else 0                                                              # rows of free space at the top of every panel (the headers of panel C)
  ylim_ <- scale_y_continuous(limits = c(0.5, n + top_pad + 0.5), expand = c(0, 0)); no_y <- theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank())
  pm <- theme(plot.margin = margin(6, 4, 6, 4), panel.grid.minor = element_blank())
  pT <- ggplot() + band_layer +
    geom_segment(data = vseg, aes(x = x, xend = x, y = ymin, yend = ymax, colour = phy), linewidth = 0.35) + geom_segment(data = hseg, aes(x = x, xend = xend, y = y, yend = y, colour = phy), linewidth = 0.35) +
    geom_segment(data = stub, aes(x = x, xend = xend, y = y, yend = y), colour = "grey40", linewidth = 0.35) +
    geom_point(data = fam, aes(5, y, colour = phy), size = 1.1) + geom_text(data = fam, aes(5.2, y, label = name), hjust = 0, size = 2.2) +
    geom_rect(data = blocks, aes(xmin = -0.95, xmax = -0.2, ymin = ymin + 0.1, ymax = ymax - 0.1, fill = phy), inherit.aes = FALSE) +
    geom_text(data = blocks[blocks$tips >= CFG$min_tips_label, ], aes(-0.575, (ymin + ymax) / 2, label = phy), angle = 90, size = 2.6, colour = "white", fontface = "bold") +
    scale_colour_manual(values = pal, guide = "none") + scale_fill_manual(values = pal, guide = "none") + scale_x_continuous(limits = c(-0.95, 9), expand = c(0, 0), breaks = 0, labels = "0") + ylim_ +
    labs(x = " ", y = NULL) + theme_fig + pm + no_y + theme(plot.margin = margin(6, 4, 6, 24), axis.text.x = element_text(colour = "white"), axis.ticks.x = element_line(colour = "white"), axis.title.x = element_text(colour = "white"), panel.grid.major.x = element_blank(), panel.border = element_blank())
  pP <- ggplot(fam) + band_layer + geom_col(aes(prevalence_alpine, y), orientation = "y", width = 0.7, fill = CFG$cols[["Alpine"]]) + geom_point(aes(prevalence_lowland, y), colour = "grey35", size = 0.9, shape = 18) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0), labels = scales::percent, breaks = c(0, 0.5, 1)) + ylim_ + labs(x = "Prevalence in alpine samples (bars);\nlowland (diamonds)", y = NULL) + theme_fig + pm + no_y
  fam$dir <- ifelse(!fam$sig, "ns", ifelse(fam$lfc > 0, "Alpine", "Lowland"))
  pB <- ggplot(fam) + band_layer + geom_vline(xintercept = 0, colour = "grey30", linetype = "dashed", linewidth = 0.4) +
    geom_linerange(aes(xmin = lo, xmax = hi, y = y, colour = dir), linewidth = 0.4) + geom_point(aes(lfc, y, colour = dir, shape = sig), size = 1.5, stroke = 0.7) +
    scale_colour_manual(values = c(CFG$cols[c("Alpine", "Lowland")], ns = CFG$cols[["ns"]]), guide = "none") + scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") + scale_x_continuous(limits = xlim_B, expand = expansion(mult = 0.04)) + ylim_ +
    labs(x = "log2 fold change in CLR abundance\n(alpine - lowland; > 0: more abundant in alpine soils)", y = NULL) + theme_fig + pm + no_y + theme(plot.margin = margin(6, 4, 6, 24))
  pC <- NULL
  if (!is.null(cor_long)) {
    cd <- inner_join(cor_long, fam[, c("family", "y")], by = "family"); cd$x <- match(cd$variable, vars_c); hd <- unique(cd[, c("variable", "envfit_r2", "x")])
    pC <- ggplot(cd) + geom_tile(aes(x, y, fill = rho), width = 0.94, height = 1) + geom_point(data = cd[cd$fdr < CFG$fdr, ], aes(x, y), size = 0.55, colour = "grey10") +
      geom_text(data = hd, aes(x, n + 1.1, label = sprintf("%.2f", envfit_r2)), size = 2.4, fontface = "bold") + geom_text(data = hd, aes(x, n + 2.3, label = variable), angle = 90, hjust = 0, size = 2.5) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-0.6, 0.6), oob = scales::squish, name = "rho", guide = guide_colourbar(barheight = 4, barwidth = 0.5)) +
      scale_x_continuous(limits = c(0.4, length(vars_c) + 0.6), expand = c(0, 0), breaks = NULL) + ylim_ + labs(x = "Spearman correlation with the envfit variables\n(numbers on top: envfit r2; dots: FDR < 0.05)", y = NULL) + theme_fig + pm + no_y + theme(plot.margin = margin(6, 4, 6, 24), legend.position = "right", panel.grid.major.x = element_blank())
  }
  hdr <- function(txt) if (top_pad > 0) annotate("text", x = -Inf, y = n + top_pad - 0.3, label = txt, hjust = -0.04, vjust = 1, size = 3, fontface = "bold") else NULL
  pT <- pT + hdr("Taxonomy (domain to family)"); pP <- pP + hdr("Prevalence"); pB <- pB + hdr("Differential abundance"); if (!is.null(pC)) pC <- pC + hdr("Link with envfit")
  cap <- paste(cap_head, "Taxonomy and differential abundance of the", CFG$tax_rank, "level, alpine vs lowland soils. A: taxonomic tree (domain, phylum, class, order, family; branches and phylum blocks coloured by phylum) of", cap_sel, "present in at least", sprintf("%.0f%%", 100 * CFG$min_prev),
               "of the alpine or of the lowland samples, with the prevalence of every family in all alpine samples (bars) and in the lowland samples (diamonds); chloroplast, mitochondria and eukaryote OTUs removed.",
               "B: log2 fold change of the centred log-ratio (CLR) abundance, alpine minus lowland, with its 95% interval, from one GAM per family (mgcv::bam) with the covariates log10 sequencing depth and the share of unmapped reads (smooths), latitude and longitude (spherical spline) and a random effect of the cluster (study / place), on at most", CFG$max_per_cluster,
               "samples per cluster; blue: enriched in the alpine soils, green: in the lowland soils (filled: FDR <", paste0(CFG$fdr, "; Benjamini-Hochberg across the families), grey: not significant. CLR values are relative to the geometric mean of the community: a fold change is relative, not absolute."),
               if (!is.null(pC)) paste("C: link with the envfit results of script 4: Spearman correlation between the CLR abundance of every family and each of the", length(vars_c), "variables with the highest envfit r2 (", CFG$c_group, "; numbers on top: that r2), across alpine clusters (cluster means); red: the family is more abundant where the variable is high; dots: FDR <", paste0(CFG$fdr, " per variable across the families. The envfit r2 comes from the ordination of the whole composition; the correlations show which families carry it.")) else "")
  cl_ <- strwrap(cap, 215); np <- if (is.null(pC)) 3 else 4; w_in <- if (is.null(pC)) 15 else 18; h_in <- max(7, 0.145 * (n + top_pad) + 2.2)
  png(file.path(out_folder, file), width = w_in, height = h_in + 0.17 * length(cl_), units = "in", res = 300)
  grid.newpage(); pushViewport(viewport(layout = grid.layout(2, np, widths = unit(c(5.3, 2.0, 3.6, if (!is.null(pC)) 4.0), "null"), heights = unit.c(unit(1, "null"), unit(0.17 * length(cl_) + 0.2, "in")))))
  plots <- list(pT, pP, pB, pC); plots <- plots[!vapply(plots, is.null, TRUE)]
  gs <- tryCatch({ g <- lapply(plots, ggplotGrob); if (length(unique(vapply(g, function(z) length(z$heights), 1L))) != 1) stop("the plots have different numbers of rows")
    mh <- do.call(grid::unit.pmax, lapply(g, function(z) z$heights)); lapply(g, function(z) { z$heights <- mh; z }) }, error = function(e) { message("  (the panels could not be aligned exactly: ", conditionMessage(e), ")"); NULL })   # the same row heights: the rows of the panels line up
  for (i in seq_along(plots)) { if (!is.null(gs)) { pushViewport(viewport(layout.pos.row = 1, layout.pos.col = i)); grid.draw(gs[[i]]); popViewport() } else print(plots[[i]], vp = viewport(layout.pos.row = 1, layout.pos.col = i)) }
  tag <- function(txt, col) grid.text(txt, x = unit(0.05, "in"), y = unit(1, "npc") - unit(0.15, "in"), just = c("left", "top"), gp = gpar(fontsize = 16, fontface = "bold"), vp = viewport(layout.pos.row = 1, layout.pos.col = col))
  tag("A", 1); tag("B", 3); if (!is.null(pC)) tag("C", 4)
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1:np)); grid.text(paste(cl_, collapse = "\n"), x = 0.01, just = "left", gp = gpar(fontsize = 9)); popViewport(2); dev.off()

  if (!is.null(rds)) saveRDS(list(tree = pT, prevalence = pP, lfc = pB, envfit = pC), rds)
  invisible(NULL)
}
cap_main <- if (sum(fam_all$shown_in_figure) < n_tested) sprintf("%d of the %d families tested (%s; all of them in tables/family_differential_abundance.csv and in the supplementary figure), the tested families being those", sum(fam_all$shown_in_figure), n_tested,
                      switch(CFG$fig_select, prevalence = "the most prevalent in the alpine samples", effect = "the largest fold changes, the significant first", both = "half the most prevalent in the alpine samples and half the largest fold changes, the significant first")) else "the families"
draw_fig(fam_all[fam_all$shown_in_figure, ], "fig5_family_tree_differential_abundance.png", "Figure 5.", cap_main, rds = file.path(out_folder, "fig5_panels.rds"))

# SUPPLEMENTARY: all the tested families, in taxonomic order, a page of at most si_per_page families (a phylum is kept on one page when it fits)
pages <- NULL
if (isTRUE(CFG$si_all_families) && sum(fam_all$shown_in_figure) < n_tested) {
  ph <- ifelse(nzchar(fam_all$r_phylum), fam_all$r_phylum, "(no phylum)"); runs <- split(seq_len(n_tested), factor(ph, levels = unique(ph))); per <- CFG$si_per_page; pages <- list(); cur <- integer(0)
  flush <- function() { if (length(cur)) pages[[length(pages) + 1]] <<- cur; cur <<- integer(0) }
  for (idx in runs) { if (length(idx) > per) { flush(); for (ch in split(idx, ceiling(seq_along(idx) / per))) pages[[length(pages) + 1]] <- ch }
                      else if (length(cur) + length(idx) > per) { flush(); cur <- idx } else cur <- c(cur, idx) }
  flush(); K <- length(pages)
  for (k in seq_len(K)) draw_fig(fam_all[pages[[k]], ], sprintf("figS5_family_tree_differential_abundance_all_page%02d_of_%02d.png", k, K), sprintf("Supplementary Figure 5 (page %d of %d).", k, K),
                                 sprintf("the families %d to %d of the %d tested (all the families, in taxonomic order; the tested families being those", min(pages[[k]]), max(pages[[k]]), n_tested))
  fam_all$si_page <- NA_integer_; for (k in seq_len(K)) fam_all$si_page[pages[[k]]] <- k
  message(sprintf("Supplementary figure: all %d families on %d page(s) of at most %d", n_tested, K, per))
}
write.csv(data.frame(family = fam_all$family, name = fam_all$name, shown_in_main_figure = fam_all$shown_in_figure, supplementary_page = if (!is.null(pages)) fam_all$si_page else NA_integer_), out_t("family_figure_pages.csv"), row.names = FALSE)
message("Done. Outputs in ", out_folder)
writeLines(capture.output(sessionInfo()), file.path(out_folder, "sessionInfo.txt"))
file.create(file.path(out_folder, ".done"))   # flag file for Snakemake