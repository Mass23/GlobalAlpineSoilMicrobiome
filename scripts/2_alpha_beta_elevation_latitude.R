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
# PART B  beta diversity (community dissimilarity) vs DIFFERENCE in elevation, DIFFERENCE
#   in latitude and geographic distance, within lowland and within alpine pairs
#   * Each replicate draws n_per_group random samples from the alpine and from the lowland
#     samples, then n_pairs random pairwise distances for
#     EACH of three comparisons: within lowland, within alpine, between alpine and lowland.
#     The pairs of ALL n_rep replicates are pooled in one model per metric, with a random
#     intercept for the replicate (s(rep_id, bs = "re")): all the data are used, and the
#     differences between the replicate draws are absorbed by the random effect.
#   * Bray-Curtis (relative abundances) and Jaccard (presence/absence after rarefaction to
#     rare_depth reads), computed only for the drawn pairs from sparse matrices
#     (checked against vegan::vegdist in TEST_MODE).
#   * Beta-regression GAM on the pooled pairs: dissimilarity ~ group + s(rep_id, "re") + s(|diff in log10 depth|)
#     + s(|diff in unmapped proportion|)
#     + s(d_elev), s(d_lat), s(log10 geo distance), each with difference smooths
#     (within alpine and alpine-vs-lowland, relative to within lowland).
#   CAVEAT: pairs reuse samples (within a replicate each sample is in ~2 of the n_pairs pairs, and
#   a sample can be drawn in several replicates). The random effect for the replicate does not
#   remove that dependence, so the confidence ribbons and the p-values / stars (beta/models_summary.csv,
#   beta/group_difference_tests.csv) are descriptive and optimistic.
#
# FIGURES (in <out_folder>)
#   fig2_alpha_beta_gradients.png : A, B = alpha diversity vs elevation and latitude;
#                                   C, D = Bray-Curtis vs elevation and latitude difference
#   figS1_beta_jaccard.png        : Supplementary Figure 1 = same as C, D with Jaccard
#   alpha/ and beta/              : full grids (all metrics / variables / specs), tables, diagnostics
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
  abs_lat     = TRUE,           # TRUE: |latitude|; FALSE: signed latitude
  k_space     = 60,             # basis dimension of the sos spatial spline
  k_space_alpine = 30,          # Supplementary Figure 3: basis dimension of the alpine-specific surface
  grid_res    = 2,              # Supplementary Figure 3: grid resolution (degrees)
  min_support = 5,              # Supplementary Figure 3: min. samples within the 3x3 cell neighbourhood to show a cell
  map_ylim    = c(-60, 83),     # Supplementary Figure 3: latitude window of the maps
  n_rug       = 3000,           # rug marks drawn per group (plot only)
  moran_n_max = 1500, moran_nperm = 199,
  # --- part B ---
  n_per_group = 1000,           # random samples per group per replicate
  n_rep       = 20,             # replicate random subsamples (each re-drawn from scratch), pooled with a random effect
  rare_depth  = 1000,           # rarefaction depth for Jaccard (samples below are dropped)
  n_pairs     = 1000,           # random pairwise distances per comparison type, replicate and pair set
  min_pairs   = 100,            # comparison types with fewer pairs are dropped from a model
  k_smooth    = 5,              # basis dimension of every 1-D smooth (alpha and beta models): small = smooth curves
  bs_smooth   = "ts",           # thin plate spline WITH shrinkage (the null space is penalised too)
  select      = FALSE,          # TRUE would add a second shrinkage penalty on top of "ts"
  beta_points = "all",          # beta panels: "all" = every pair of every replicate as a point; "binned" = 15 binned means per group
  beta_point_alpha = 0.08,      # transparency of those points
  # --- main figure ---
  main_alpha_metrics = c("richness", "shannon", "pielou"),   # alpha metrics shown in panels A, B
  main_beta_spec     = "all_pairs",                          # pair set shown in panels C, D
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
  for (i in seq_len(nr)) for (j in seq_along(rows[[i]]))
    print(rows[[i]][[j]], vp = viewport(layout.pos.row = i,
                                        layout.pos.col = if (length(rows[[i]]) == 1) 1:2 else j))
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

sm1 <- function(v, by = NULL) sprintf("s(%s, bs = '%s', k = %d%s)", v, CFG$bs_smooth, CFG$k_smooth,
                                      if (is.null(by)) "" else paste0(", by = ", by))
base_rhs <- paste("Alpine +", sm1("log10_reads"), "+", sm1("elevation"), "+", sm1("elevation", "AlpineO"),
                  "+", sm1("glat"), "+", sm1("glat", "AlpineO"))
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
  excl <- intersect("s(lat,lon)", vapply(m$smooth, function(s) s$label, ""))
  bind_rows(lapply(levels(d$Alpine), function(a) {
    sub <- d[d$Alpine == a, ]
    nd <- data.frame(Alpine = factor(a, levels = levels(d$Alpine)),
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
                                           "+ s(lat, lon, bs = 'sos', by = AlpineO, k = %d)"),
                                     sm1("log10_reads"), CFG$k_space, CFG$k_space_alpine))
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
  nd$fit <- as.numeric(predict(m_sp, nd, type = "response", newdata.guaranteed = TRUE, discrete = FALSE))
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

## ---- B1. Replicate subsamples and counts -------------------------------------------------------
sel <- lapply(seq_len(CFG$n_rep), function(r) {
  set.seed(100 + r)
  meta %>%
    group_by(Alpine) %>% slice_sample(n = CFG$n_per_group) %>% ungroup()
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
  iL <- which(S$Alpine == "Lowland"); iA <- which(S$Alpine == "Alpine")

  # n_pairs random pairs (indices into S). within = TRUE: unordered pairs inside one set
  draw_pairs <- function(i1, i2, within) {
    cand <- data.frame(a = sample.int(length(i1), 4 * CFG$n_pairs, replace = TRUE),
                       b = sample.int(length(i2), 4 * CFG$n_pairs, replace = TRUE))
    if (within) {
      cand <- cand[cand$a != cand$b, ]
      cand <- data.frame(a = pmin(cand$a, cand$b), b = pmax(cand$a, cand$b))
    }
    cand <- head(cand[!duplicated(cand), ], CFG$n_pairs)
    data.frame(i = i1[cand$a], j = i2[cand$b])
  }
  lab_pairs <- function(d, label) data.frame(group = rep(label, nrow(d)), i = d$i, j = d$j)

  pairs <- bind_rows(lab_pairs(draw_pairs(iL, iL, TRUE),  beta_groups[1]),
                     lab_pairs(draw_pairs(iA, iA, TRUE),  beta_groups[2]),
                     lab_pairs(draw_pairs(iA, iL, FALSE), beta_groups[3]))
  pairs$spec <- "all_pairs"
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
  pairs[, c("rep_id", "spec", "group", "bray", "jaccard", "d_elev", "d_lat", "geo_km", "d_depth", "d_unm")]
}

## ---- B3. Models ---------------------------------------------------------------------------------------
beta_groups  <- c("Within lowland", "Within alpine", "Alpine vs lowland")   # reference = first
beta_metrics <- c(bray = "Bray\u2013Curtis", jaccard = "Jaccard (rarefied)")
beta_vars    <- c(d_elev = "Elevation difference (km)", d_lat = "Latitude difference (\u00B0)",
             log10_geo = "Geographic distance (log10 km)")
beta_specs   <- c(all_pairs = "Random pairs")

prep_pairs <- function(P, spec) {
  P <- as.data.frame(P[P$spec == spec, ])
  # drop comparison types with too few pairs
  keep <- intersect(beta_groups, names(which(table(P$group) >= CFG$min_pairs)))
  if (length(setdiff(beta_groups, keep))) message("    dropped (< ", CFG$min_pairs, " pairs): ",
      paste(setdiff(beta_groups, keep), collapse = ", "))
  P <- P[P$group %in% keep, ]
  P$group <- factor(P$group, levels = keep)
  P$groupO <- factor(P$group, levels = keep, ordered = TRUE)        # reference = first level
  P$rep_id <- factor(P$rep_id)                                        # replicate: random effect
  P$log10_geo <- log10(pmax(P$geo_km, 0.1))
  P
}

fit_beta <- function(P, metric) {
  n <- nrow(P); P$y <- (P[[metric]] * (n - 1) + 0.5) / n          # squeeze into (0, 1)
  sm <- unlist(lapply(names(beta_vars), function(v) c(sm1(v), sm1(v, "groupO"))))
  # nuisance covariates shared by the three comparison types: differences in depth and unmapped proportion
  terms <- c("group", sm1("d_depth"), sm1("d_unm"), sm)
  if (nlevels(P$rep_id) > 1) terms <- c(terms, "s(rep_id, bs = 're')")                  # random intercept: replicate
  bam(as.formula(paste("y ~", paste(terms, collapse = " + "))), data = P,
      family = betar(link = "logit"), method = "fREML", discrete = TRUE,
      select = CFG$select, nthreads = CFG$nthreads)
}

# Partial curve of one pair variable per group (others at group medians), from the pooled model: the replicate
# random effect is left out (a typical replicate); the ribbon is the model's 95% confidence interval
predict_var <- function(m, P, v) {
  excl <- intersect("s(rep_id)", vapply(m$smooth, function(s) s$label, ""))
  bind_rows(lapply(levels(P$group), function(g) {
    sub <- P[P$group == g, ]
    x <- seq(quantile(sub[[v]], .01), quantile(sub[[v]], .99), length.out = 80)
    nd <- data.frame(group = factor(g, levels = levels(P$group)),
                     groupO = factor(g, levels = levels(P$group), ordered = TRUE),
                     d_elev = median(sub$d_elev), d_lat = median(sub$d_lat),
                     log10_geo = median(sub$log10_geo),
                     d_depth = median(sub$d_depth), d_unm = median(sub$d_unm),
                     rep_id = factor(levels(P$rep_id)[1], levels = levels(P$rep_id)))
    nd <- nd[rep(1, length(x)), ]; nd[[v]] <- x
    p <- predict(m, nd, type = "link", se.fit = TRUE, exclude = excl, newdata.guaranteed = TRUE, discrete = FALSE)
    li <- family(m)$linkinv
    data.frame(group = g, x = x, fit = li(p$fit), lo = li(p$fit - 1.96 * p$se.fit), hi = li(p$fit + 1.96 * p$se.fit))
  }))
}

# Binned means of the observed dissimilarities (15 quantile bins per group, all replicates): data check
bin_obs <- function(P, metric, v)
  P %>% group_by(group) %>% mutate(bin = ntile(.data[[v]], 15)) %>%
    group_by(group, bin) %>%
    summarise(x = mean(.data[[v]]), y = mean(.data[[metric]]), .groups = "drop")

# all replicates: the random pairs and their dissimilarities. Saved (beta/pairs_all.rds) and reused on a rerun if the
# settings and the drawn samples are the same, so a problem in the models does not cost the distance computation again.
pairs_key <- list(n_rep = CFG$n_rep, n_per_group = CFG$n_per_group, n_pairs = CFG$n_pairs, rare_depth = CFG$rare_depth,
                  n_ids = length(all_ids), ids_chars = sum(nchar(all_ids)), test = TEST_MODE)
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

preds <- list(); obs <- list(); summ <- list(); pair_counts <- list(); beta_raw <- list()
for (spec in names(beta_specs)) {
  P <- prep_pairs(P_all, spec)
  pair_counts[[spec]] <- P %>% count(rep_id, group) %>% mutate(spec = spec)
  beta_raw[[spec]] <- P[, c("group", names(beta_vars), names(beta_metrics))]      # every pair of every replicate (for the plots)
  message(sprintf("%s: %d pairs from %d replicates", spec, nrow(P), nlevels(P$rep_id)))
  for (metric in names(beta_metrics)) {
    message("  fitting ", spec, " / ", metric, " (pooled, random effect: replicate)")
    m <- tryCatch(fit_beta(P, metric), error = function(e) { message("  FAILED: ", conditionMessage(e)); NULL })
    if (is.null(m)) next
    key <- paste(spec, metric)
    preds[[key]] <- bind_rows(lapply(names(beta_vars), function(v)
      cbind(spec = spec, metric = metric, variable = v, predict_var(m, P, v))))
    obs[[key]] <- bind_rows(lapply(names(beta_vars), function(v)
      cbind(spec = spec, metric = metric, variable = v, bin_obs(P, metric, v))))
    sm_m <- summary(m)
    st <- as.data.frame(sm_m$s.table); st$term <- rownames(st)
    one <- function(x) if (length(x) == 1 && is.numeric(x)) as.numeric(x) else NA_real_     # always exactly one number
    # standard deviation of the replicate random effect (NA if this mgcv version returns it in another shape)
    rep_sd <- tryCatch({
      vc <- gam.vcomp(m, rescale = FALSE)
      if (is.list(vc) && !is.matrix(vc)) vc <- vc$vc
      one(if (is.matrix(vc)) vc["s(rep_id)", 1] else vc["s(rep_id)"])
    }, error = function(e) NA_real_)
    summ[[key]] <- data.frame(spec = spec, metric = metric, dev_expl = one(sm_m$dev.expl), n_pairs = nrow(P),
                              replicate_sd = rep_sd, st, check.names = FALSE, stringsAsFactors = FALSE)
    saveRDS(m, out_b(paste0("model_", spec, "_", metric, ".rds")))
  }
}
preds <- bind_rows(preds); obs <- bind_rows(obs)
summ_all <- bind_rows(summ)
write.csv(summ_all, out_b("models_summary.csv"), row.names = FALSE)   # p-values: descriptive only

# Difference smooths vs the within-lowland reference (pooled model; descriptive: pairs reuse samples)
pcol <- grep("^p[-.]value$", names(summ_all), value = TRUE)[1]
beta_tests <- summ_all %>%
  filter(grepl("groupO", term)) %>%
  mutate(variable = ifelse(grepl("d_elev", term), "d_elev", ifelse(grepl("d_lat", term), "d_lat", "log10_geo")),
         contrast = sub(".*groupO", "", term), p = .data[[pcol]]) %>%
  select(spec, metric, variable, contrast, p)
write.csv(beta_tests, out_b("group_difference_tests.csv"), row.names = FALSE)
star_of <- function(p, contrast, which) { v <- p[contrast == which]; if (length(v)) stars(v[1]) else "n/a" }
beta_label <- function(tt) paste0(" vs within lowland:\n within alpine ", star_of(tt$p, tt$contrast, "Within alpine"),
                                   " | alpine-lowland ", star_of(tt$p, tt$contrast, "Alpine vs lowland"))
write.csv(bind_rows(pair_counts), out_b("pair_counts.csv"), row.names = FALSE)
write.csv(preds, out_b("predicted_curves.csv"), row.names = FALSE)

## ---- B4. Beta-diversity figures ----------------------------------------------------------------------------------------
# One panel: one metric against one pair variable (line = pooled model, ribbon = its 95% confidence interval,
# points = every pair of every replicate, or binned means if CFG$beta_points = "binned")
pt_size  <- if (CFG$beta_points == "all") 0.35 else 0.9
pt_alpha <- if (CFG$beta_points == "all") CFG$beta_point_alpha else 0.6
pts_txt  <- if (CFG$beta_points == "all") "points = every pair of every replicate (transparent)" else "points = binned observed means"
beta_panel <- function(mt, vv, sp, tag) {
  sel_rows <- function(d) d[d$spec == sp & d$metric == mt & d$variable == vv, ]
  pm <- sel_rows(preds)
  ob <- if (CFG$beta_points == "all") { raw <- beta_raw[[sp]]; data.frame(group = raw$group, x = raw[[vv]], y = raw[[mt]]) }
        else sel_rows(obs)
  pm$mid <- pm$fit
  for (nm in c("ob", "pm")) {
    d <- get(nm); d$group <- factor(d$group, levels = beta_groups); assign(nm, d)
  }
  tt <- beta_tests[beta_tests$spec == sp & beta_tests$metric == mt & beta_tests$variable == vv, ]
  an <- data.frame(label = beta_label(tt))
  ggplot() +
    geom_ribbon(data = pm, aes(x, ymin = lo, ymax = hi, fill = group), alpha = 0.25) +
    geom_point(data = ob, aes(x, y, colour = group), size = pt_size, alpha = pt_alpha, shape = 16) +
    geom_line(data = pm, aes(x, mid, colour = group), linewidth = 1.2) +
    geom_text(data = an, aes(label = label), x = -Inf, y = Inf, hjust = 0, vjust = 1.3,
              size = 3, lineheight = 0.95, inherit.aes = FALSE) +
    scale_colour_manual(values = CFG$cols_beta, name = NULL) +
    scale_fill_manual(values = CFG$cols_beta, guide = "none") +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1, linewidth = 1.2))) +
    labs(tag = tag, x = beta_vars[[vv]], y = paste(beta_metrics[[mt]], "dissimilarity")) +
    theme_fig
}

# Full grids (all metrics x all pair variables, one file per pair set), kept in beta/
for (sp in names(beta_specs)) {
  pr <- preds %>% filter(spec == sp)
  ob <- if (CFG$beta_points == "all") {
    raw <- beta_raw[[sp]]
    bind_rows(lapply(names(beta_metrics), function(mt) bind_rows(lapply(names(beta_vars), function(v)
      data.frame(metric = mt, variable = v, group = raw$group, x = raw[[v]], y = raw[[mt]])))))
  } else obs %>% filter(spec == sp)
  if (!nrow(pr)) next
  lab <- function(d) {
    d$metric   <- factor(beta_metrics[d$metric], levels = beta_metrics)
    d$variable <- factor(beta_vars[d$variable], levels = beta_vars)
    d$group    <- factor(d$group, levels = beta_groups)
    d
  }
  pr <- lab(pr); ob <- lab(ob)
  pm <- pr; pm$mid <- pm$fit
  ann_g <- beta_tests %>% filter(spec == sp) %>% group_by(metric, variable) %>%
    group_modify(~ data.frame(label = beta_label(.x))) %>% ungroup()
  ann_g$metric   <- factor(beta_metrics[ann_g$metric], levels = beta_metrics)
  ann_g$variable <- factor(beta_vars[ann_g$variable], levels = beta_vars)
  g <- ggplot() +
    geom_ribbon(data = pm, aes(x, ymin = lo, ymax = hi, fill = group), alpha = 0.25) +
    geom_point(data = ob, aes(x, y, colour = group), size = pt_size, alpha = pt_alpha, shape = 16) +
    geom_line(data = pm, aes(x, mid, colour = group), linewidth = 1.2) +
    geom_text(data = ann_g, aes(label = label), x = -Inf, y = Inf, hjust = 0, vjust = 1.3,
              size = 3, lineheight = 0.95, inherit.aes = FALSE) +
    facet_grid(metric ~ variable, scales = "free", switch = "both") +
    scale_colour_manual(values = CFG$cols_beta, name = NULL) +
    scale_fill_manual(values = CFG$cols_beta, guide = "none") +
    labs(x = NULL, y = "Dissimilarity", title = beta_specs[[sp]]) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1, linewidth = 1.2))) +
    theme_bw(base_size = 12) +
    theme(strip.placement = "outside", strip.background = element_blank(),
          strip.text = element_text(face = "bold"), panel.grid.minor = element_blank(),
          legend.position = "bottom")
  ggsave(out_b(paste0("beta_full_grid_", sp, ".png")), g, width = 12, height = 7, dpi = 300)
}

beta_caption <- sprintf(paste(
  "Dissimilarity between pairs of samples (within lowland, within alpine, between alpine and lowland; pair set: %s) against the difference in elevation",
  "and in latitude between the two samples. One beta-regression GAM per metric on the pooled pairs of %d replicate subsamples",
  "(each: %d random samples per group, %d random pairwise distances per comparison type), with a random intercept for the replicate",
  "and adjusted for geographic distance and for differences in sequencing depth and unmapped proportion: lines = fitted curves,",
  "ribbons = 95%% confidence intervals of the pooled model, points = binned observed means. Pairs reuse samples, so the intervals and",
  "the stars (p of the difference smooths vs the within-lowland curve: *** p < 0.001, ** p < 0.01, * p < 0.05, ns) are descriptive",
  "and optimistic."),
  beta_specs[[CFG$main_beta_spec]], CFG$n_rep, CFG$n_per_group, CFG$n_pairs)
beta_caption <- sub("points = binned observed means", pts_txt, beta_caption, fixed = TRUE)

# Supplementary Figure 1: Jaccard (rarefied presence/absence)
pS1A <- beta_panel("jaccard", "d_elev", CFG$main_beta_spec, "A")
pS1B <- beta_panel("jaccard", "d_lat",  CFG$main_beta_spec, "B")
assemble_fig(file.path(out_folder, "figS1_beta_jaccard.png"), list(list(pS1A, pS1B)),
             get_legend(pS1A), paste("Supplementary Figure 1.", beta_caption), 12, 5.6)

}  # end RUN_BETA

## ======================= MAIN FIGURE (needs both parts) =====================================
if (RUN_ALPHA && RUN_BETA) {
  pA <- alpha_panel(grad_lab[["elevation"]], "A")
  pB <- alpha_panel(grad_lab[["glat"]], "B")
  pC <- beta_panel("bray", "d_elev", CFG$main_beta_spec, "C")
  pD <- beta_panel("bray", "d_lat",  CFG$main_beta_spec, "D")
  main_caption <- paste(
    "A, B: alpha diversity along elevation (A) and latitude (B); GAM partial effects adjusted for sequencing depth and unmapped proportion, solid = net of space",
    "(95% CI), dashed = no spatial term; stars = test that the alpine gradient differs in shape from the lowland",
    "one (net model; *** p < 0.001, ** p < 0.01, * p < 0.05, ns). C, D: Bray-Curtis dissimilarity between pairs of samples (within lowland, within alpine, between alpine and lowland) against the",
    sprintf("difference in elevation (C) and latitude (D); one model on the pooled pairs of %d replicate subsamples (each: %d random samples per group,",
            CFG$n_rep, CFG$n_per_group),
    sprintf("%d random pairwise distances per comparison type) with a random intercept for the replicate; lines = fitted curves, ribbons = 95%% confidence intervals, points = binned observed means.",
            CFG$n_pairs),
    "Stars in C, D: p of the difference smooths vs within lowland (descriptive and optimistic: pairs reuse samples).",
    "Jaccard: Supplementary Figure 1.")
  main_caption <- sub("points = binned observed means", pts_txt, main_caption, fixed = TRUE)
  assemble_fig(file.path(out_folder, "fig2_alpha_beta_gradients.png"),
               list(list(pA, pB), list(pC, pD)), list(get_legend(pA), get_legend(pC)), main_caption,
               width = 14, height = 10.5, row_heights = c(1.1, 1))
  saveRDS(list(A = pA, B = pB, C = pC, D = pD), file.path(out_folder, "fig2_panels.rds"))
} else message("Main figure skipped: it needs RUN_ALPHA and RUN_BETA both TRUE.")

message("Done. Outputs in ", out_folder)
writeLines(capture.output(sessionInfo()), file.path(out_folder, "sessionInfo.txt"))
file.create(file.path(out_folder, ".done"))   # flag file for Snakemake