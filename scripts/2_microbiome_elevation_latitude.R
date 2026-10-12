#!/usr/bin/env Rscript
# =============================================================================
# 2_alpha_beta_elevation_latitude.R  --  Figure 2
# Alpha AND beta diversity along ELEVATION and LATITUDE, alpine vs lowland.
#
# PART A  alpha diversity (Observed OTUs, Shannon, Pielou J) vs elevation and latitude
#   Two GAM specifications per metric:
#     net : alpine + s(log10 reads)
#           + s(elevation) + s(elevation, by = alpine-ordered)   [lowland = reference,
#           + s(latitude)  + s(latitude,  by = alpine-ordered)    by-smooths = alpine - lowland]
#           + sos spatial spline            -> gradients net of spatial structure
#     raw : the same WITHOUT the spatial spline -> raw gradients
#   The ordered-factor difference smooths test whether each gradient differs in shape
#   between alpine and lowland (alpha/group_difference_tests.csv, printed on the figure).
#   NB: the sos spline absorbs much of the latitudinal signal (concurvity).
#
# PART B  beta diversity (community dissimilarity) vs DIFFERENCE in elevation, DIFFERENCE
#   in latitude and geographic distance, within lowland and within alpine pairs
#   * Each replicate draws n_per_group random samples from the alpine and from the lowland
#     samples, then n_pairs random pairwise distances for
#     EACH of three comparisons: within lowland, within alpine, between alpine and lowland.
#     n_rep replicates show sampling stability.
#   * Bray-Curtis (relative abundances) and Jaccard (presence/absence after rarefaction to
#     rare_depth reads), computed only for the drawn pairs from sparse matrices
#     (checked against vegan::vegdist in TEST_MODE).
#   * Beta-regression GAM on pairs: dissimilarity ~ group
#     + s(d_elev), s(d_lat), s(log10 geo distance), each with difference smooths
#     (within alpine and alpine-vs-lowland, relative to within lowland).
#   CAVEAT: pairs reuse samples (each sample is in ~2 of the n_pairs pairs), so they are
#   not fully independent: p-values in beta/models_summary.csv are descriptive only; judge
#   robustness by the stability across replicates.
#
# FIGURES (in <out_folder>)
#   fig2_alpha_beta_gradients.png : A, B = alpha diversity vs elevation and latitude;
#                                   C, D = Bray-Curtis vs elevation and latitude difference
#   figS1_beta_jaccard.png        : Supplementary Figure 1 = same as C, D with Jaccard
#   alpha/ and beta/              : full grids (all metrics / variables / specs), tables, diagnostics
#
# Input : <data_folder>/sample_data_filtered.csv, elevation from env_file and, for part B,
#         the filtered BIOM written by microbeatlas_filter_data.R
# Output: <out_folder>/alpha/ and <out_folder>/beta/
# Usage : Rscript 2_alpha_beta_elevation_latitude.R [data_folder] [out_folder]
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

args <- commandArgs(trailingOnly = TRUE)
data_folder <- if (length(args) >= 1) as.character(args[1]) else
  "/home/renku/work/epfl-altshuler/GlobalAlpineSoilMicrobiome/data"
out_folder  <- if (length(args) >= 2) as.character(args[2]) else
  file.path(dirname(data_folder), "results", "2_alpha_beta_elevation_latitude")

CFG <- list(
  # --- input ---
  alpha_file  = file.path(data_folder, "sample_data_filtered.csv"),
  env_file    = file.path(dirname(data_folder), "results", "all_sampled_data.csv"),  # elevation (DEM)
  alpha_id    = "Sample_acc",   # sample ID column in sample_data_filtered.csv   (check!)
  env_id      = "site",         # matching ID column in env_file                  (check!)
  env_elev    = "elevation",    # elevation column (metres) in env_file           (check!)
  biom_file   = file.path(data_folder, "microbeatlas",
                          "samples-otus.97.mapped.metag.minfilter.refilt_filtered.biom"),
  biom_id     = "Accessions",   # column of alpha_file that matches the BIOM sample IDs
  # --- part A ---
  abs_lat     = TRUE,           # TRUE: |latitude|; FALSE: signed latitude
  k_space     = 60,             # basis dimension of the sos spatial spline
  n_rug       = 3000,           # rug marks drawn per group (plot only)
  moran_n_max = 1500, moran_nperm = 199,
  # --- part B ---
  n_per_group = 1000,           # random samples per group per replicate
  n_rep       = 20,             # replicate random subsamples (each re-drawn from scratch); the 95% band needs many
  rare_depth  = 1000,           # rarefaction depth for Jaccard (samples below are dropped)
  n_pairs     = 1000,           # random pairwise distances per comparison type, replicate and pair set
  min_pairs   = 100,            # comparison types with fewer pairs are dropped from a model
  k_beta      = 6,              # basis dimension of the pair smooths
  # --- main figure ---
  main_alpha_metrics = c("richness", "shannon", "pielou"),   # alpha metrics shown in panels A, B
  main_beta_spec     = "all_pairs",                          # pair set shown in panels C, D
  # --- common ---
  nthreads    = 4,
  cols        = c(Lowland = "#3A8F4B", Alpine = "#2B6CB0"),
  cols_beta   = c("Within lowland" = "#3A8F4B", "Within alpine" = "#2B6CB0",
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
  png(file, width = width, height = height, units = "in", res = 300)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nr + 2, 2,
                                             heights = unit(c(row_heights, 0.09, 0.14), "null"))))
  for (i in seq_len(nr)) for (j in 1:2)
    print(rows[[i]][[j]], vp = viewport(layout.pos.row = i, layout.pos.col = j))
  if (!is.null(leg)) {
    legs <- if (inherits(leg, "grob")) list(leg) else leg
    for (k in seq_along(legs)) {
      pushViewport(viewport(layout.pos.row = nr + 1, layout.pos.col = if (length(legs) == 1) 1:2 else k))
      grid.draw(legs[[k]]); popViewport()
    }
  }
  pushViewport(viewport(layout.pos.row = nr + 2, layout.pos.col = 1:2))
  grid.text(paste(strwrap(caption, 160), collapse = "\n"), x = 0.01, just = "left", gp = gpar(fontsize = 9))
  popViewport(); dev.off()
}

## ---- 1. Data (shared by both parts) -----------------------------------------------------------------------
# Reads sample_data_filtered.csv even when it was written WITHOUT quotes (quote = F).
#  * commas inside {...} are protected: everything within braces is ONE field;
#  * rows that still have extra fields (unbraced commas in free text such as
#    "Washington state, USA") get the extras merged back into the geo_col column;
#  * rows that still do not fit (non-numeric coordinates, Alpine not Alpine/Lowland, too few
#    fields) are dropped and counted.
# A properly quoted CSV (first data line starts with a quote) is read with read.csv as usual.
read_filtered_table <- function(path, geo_col = "GPT_Geotext") {
  first <- readLines(path, n = 2)
  if (length(first) >= 2 && substr(first[2], 1, 1) == "\"")
    return(read.csv(path, stringsAsFactors = FALSE))
  
  lines <- readLines(path)
  hdr <- strsplit(lines[1], ",", fixed = TRUE)[[1]]; n <- length(hdr); lines <- lines[-1]
  
  m <- gregexpr("\\{[^{}]*\\}", lines)                       # protect commas inside {...}
  regmatches(lines, m) <- lapply(regmatches(lines, m), function(z) gsub(",", "\001", z, fixed = TRUE))
  parts <- strsplit(lines, ",", fixed = TRUE); rm(lines, m)
  
  g <- match(geo_col, hdr)
  extra_rows <- which(lengths(parts) > n)
  if (length(extra_rows)) {
    if (is.na(g)) stop("Rows with extra fields, but column '", geo_col, "' not found in the header")
    parts[extra_rows] <- lapply(parts[extra_rows], function(p) {
      e <- length(p) - n
      rest <- if (g + e < length(p)) p[(g + e + 1):length(p)] else character(0)
      c(p[seq_len(g - 1)], paste(p[g:(g + e)], collapse = ","), rest)
    })
  }
  ok_len <- lengths(parts) == n
  mat <- matrix(unlist(parts[ok_len], use.names = FALSE), ncol = n, byrow = TRUE); rm(parts)
  colnames(mat) <- hdr
  
  valid <- rep(TRUE, nrow(mat))
  if (all(c("Latitude", "Longitude", "Alpine") %in% hdr))
    valid <- is.finite(suppressWarnings(as.numeric(mat[, "Latitude"]))) &
    is.finite(suppressWarnings(as.numeric(mat[, "Longitude"]))) &
    mat[, "Alpine"] %in% c("Alpine", "Lowland")
  message(sprintf("%s: %d lines | %d with extra fields merged into %s | %d with too few fields | %d invalid rows dropped",
                  basename(path), length(ok_len), length(extra_rows), geo_col, sum(!ok_len), sum(!valid)))
  mat <- mat[valid, , drop = FALSE]
  
  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  df[] <- lapply(df, function(x) type.convert(gsub("\001", ",", x, fixed = TRUE), as.is = TRUE))
  df
}

raw <- read_filtered_table(CFG$alpha_file)
required <- c("Longitude", "Latitude", "Alpine", "alpha_total_n_reads",
              "alpha_n_otus", "alpha_shannon", "alpha_pielou", CFG$alpha_id,
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
      shannon     = alpha_shannon, pielou = alpha_pielou
    ) %>%
    filter(is.finite(lon), is.finite(lat), abs(lon) <= 180, abs(lat) <= 90,
           is.finite(elevation), is.finite(reads), reads > 0,
           is.finite(richness), is.finite(shannon), shannon > 0, is.finite(pielou),
           !is.na(Alpine))
  
  # Pielou J must be strictly in (0, 1) for beta regression (squeeze only if needed)
  if (any(d$pielou <= 0 | d$pielou >= 1)) {
    n <- nrow(d); d$pielou <- (d$pielou * (n - 1) + 0.5) / n
    message("Pielou J contained 0/1 values: applied the Smithson-Verkuilen squeeze.")
  }
  d$AlpineO <- factor(d$Alpine, levels = c("Lowland", "Alpine"), ordered = TRUE)   # for difference smooths
  print(table(alpine = d$Alpine))
  write.csv(select(d, any_of(c(CFG$alpha_id, "Alpine", "lon", "lat", "glat", "elevation", "reads",
                               "richness", "shannon", "pielou"))),
            out_a("analysis_table.csv"), row.names = FALSE)
  
  ## ---- A1. Models -----------------------------------------------------------------------
  families <- list(richness = nb(), shannon = Gamma(link = "log"), pielou = betar(link = "logit"))
  resp_lab <- c(richness = "Observed OTUs", shannon = "Shannon index", pielou = "Pielou J")
  
  base_rhs <- paste("Alpine + s(log10_reads, k = 6)",
                    "+ s(elevation, k = 6) + s(elevation, by = AlpineO, k = 6)",
                    "+ s(glat, k = 6) + s(glat, by = AlpineO, k = 6)")
  specs <- c(net = sprintf("%s + s(lat, lon, bs = 'sos', k = %d)", base_rhs, CFG$k_space),
             raw = base_rhs)
  spec_lab <- c(net = "Net of space (sos)", raw = "No spatial term")
  
  fits <- list()
  for (resp in names(families)) for (sp in names(specs)) {
    nm <- paste(resp, sp, sep = "__"); message("Fitting ", nm)
    fits[[nm]] <- bam(as.formula(paste(resp, "~", specs[[sp]])), data = d, family = families[[resp]],
                      method = "fREML", discrete = TRUE, select = TRUE, nthreads = CFG$nthreads)
  }
  saveRDS(fits, out_a("fits.rds"))
  
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
    writeLines(c(paste(parts[1], "~", specs[[parts[2]]]), "", capture.output(print(sm)),
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
  # group median, depth at the overall median, spatial term excluded.
  predict_gradient <- function(m, var) {
    excl <- intersect("s(lat,lon)", vapply(m$smooth, function(s) s$label, ""))
    bind_rows(lapply(levels(d$Alpine), function(a) {
      sub <- d[d$Alpine == a, ]
      nd <- data.frame(Alpine = factor(a, levels = levels(d$Alpine)),
                       AlpineO = factor(a, levels = levels(d$Alpine), ordered = TRUE),
                       elevation = median(sub$elevation), glat = median(sub$glat),
                       log10_reads = median(d$log10_reads), lat = median(sub$lat), lon = median(sub$lon))
      nd <- nd[rep(1, 100), ]
      nd[[var]] <- seq(quantile(sub[[var]], .01), quantile(sub[[var]], .99), length.out = 100)
      p <- predict(m, nd, type = "link", se.fit = TRUE, exclude = excl, newdata.guaranteed = TRUE)
      li <- family(m)$linkinv
      data.frame(Alpine = nd$Alpine, x = nd[[var]], fit = li(p$fit),
                 lo = li(p$fit - 1.96 * p$se.fit), hi = li(p$fit + 1.96 * p$se.fit))
    }))
  }
  
  grad_lab <- c(elevation = "Elevation (km)",
                glat = if (CFG$abs_lat) "Absolute latitude (\u00B0)" else "Latitude (\u00B0)")
  pred <- bind_rows(lapply(names(fits), function(nm) {
    parts <- strsplit(nm, "__")[[1]]
    bind_rows(lapply(names(grad_lab), function(v)
      cbind(response = parts[1], spec = parts[2], gradient = grad_lab[[v]],
            predict_gradient(fits[[nm]], v))))
  }))
  pred$response <- factor(resp_lab[pred$response], levels = resp_lab)
  pred$gradient <- factor(pred$gradient, levels = grad_lab)
  pred$spec     <- factor(spec_lab[pred$spec], levels = spec_lab)
  write.csv(pred, out_a("predicted_gradients.csv"), row.names = FALSE)
  
  # Alpine-vs-lowland difference smooths: does the gradient shape differ between groups?
  diff_tests <- bind_rows(smooths) %>%
    filter(grepl("AlpineO", term)) %>%
    mutate(gradient = factor(ifelse(grepl("elevation", term), grad_lab[["elevation"]], grad_lab[["glat"]]),
                             levels = grad_lab),
           response = factor(resp_lab[response], levels = resp_lab))
  write.csv(diff_tests, out_a("group_difference_tests.csv"), row.names = FALSE)
  ann <- diff_tests %>% filter(spec == "net") %>%
    mutate(label = sprintf("Alpine vs lowland shape: p %s",
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
  ggsave(out_a("alpha_full_grid.png"), g, width = 10, height = 10, dpi = 300)
  saveRDS(g, out_a("alpha_full_grid_plot.rds"))
  
  # Panel builder for the main figure: one gradient, the chosen alpha metrics as facets
  main_alpha_labels <- resp_lab[CFG$main_alpha_metrics]
  ann_short <- ann %>% mutate(label = sub("Alpine vs lowland shape: ", "", label))
  alpha_panel <- function(gr_label, tag) {
    pr <- pred[pred$gradient == gr_label & pred$response %in% main_alpha_labels, ]
    pr$response <- droplevels(pr$response)
    an <- ann_short[ann_short$gradient == gr_label & ann_short$response %in% main_alpha_labels, ]
    rg <- rug[rug$gradient == gr_label, ]
    ggplot(pr, aes(x, fit, colour = Alpine)) +
      geom_ribbon(data = pr[pr$spec == spec_lab[["net"]], ],
                  aes(ymin = lo, ymax = hi, fill = Alpine), alpha = 0.25, colour = NA) +
      geom_line(aes(linetype = spec), linewidth = 1) +
      geom_rug(data = rg, aes(x = x, colour = Alpine), inherit.aes = FALSE,
               sides = "b", alpha = 0.05, length = unit(0.03, "npc")) +
      geom_text(data = an, aes(label = label), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
                size = 3, inherit.aes = FALSE) +
      facet_wrap(~response, nrow = 1, scales = "free_y") +
      scale_colour_manual(values = CFG$cols, name = NULL) +
      scale_fill_manual(values = CFG$cols, guide = "none") +
      scale_linetype_manual(values = c("solid", "dashed"), name = NULL) +
      labs(tag = tag, x = gr_label, y = NULL) + theme_fig
  }
  
}  # end RUN_ALPHA

## =============================== PART B: BETA DIVERSITY ================================
if (RUN_BETA) {
  
  meta <- raw %>%
    transmute(
      id        = .data[[CFG$biom_id]],
      Alpine    = factor(Alpine, levels = c("Lowland", "Alpine")),
      lon       = Longitude, lat = Latitude,
      elevation = .elev_m / 1000                                      # km
    ) %>%
    filter(is.finite(lon), is.finite(lat), abs(lon) <= 180, abs(lat) <= 90,
           is.finite(elevation), !is.na(Alpine), !is.na(id))
  print(table(group = meta$Alpine))
  
  ## ---- B1. Replicate subsamples and counts -------------------------------------------------------
  sel <- lapply(seq_len(CFG$n_rep), function(r) {
    set.seed(100 + r)
    meta %>%
      group_by(Alpine) %>% slice_sample(n = CFG$n_per_group) %>% ungroup()
  })
  if (CFG$n_rep < 10) message("Note: a 95% interval across replicate curves is poorly defined with n_rep < 10")
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
    
    pairs$rep <- r
    pairs$d_elev <- abs(S$elevation[pairs$i] - S$elevation[pairs$j])
    pairs$d_lat  <- abs(S$lat[pairs$i] - S$lat[pairs$j])
    pairs$geo_km <- geosphere::distHaversine(cbind(S$lon[pairs$i], S$lat[pairs$i]),
                                             cbind(S$lon[pairs$j], S$lat[pairs$j])) / 1000
    pairs[, c("rep", "spec", "group", "bray", "jaccard", "d_elev", "d_lat", "geo_km")]
  }
  
  ## ---- B3. Models ---------------------------------------------------------------------------------------
  beta_groups  <- c("Within lowland", "Within alpine", "Alpine vs lowland")   # reference = first
  beta_metrics <- c(bray = "Bray\u2013Curtis", jaccard = "Jaccard (rarefied)")
  beta_vars    <- c(d_elev = "Elevation difference (km)", d_lat = "Latitude difference (\u00B0)",
                    log10_geo = "Geographic distance (log10 km)")
  beta_specs   <- c(all_pairs = "Random pairs")
  grid_cache <- new.env()
  
  prep_pairs <- function(P, spec) {
    P <- as.data.frame(P[P$spec == spec, ])
    # drop comparison types with too few pairs
    keep <- intersect(beta_groups, names(which(table(P$group) >= CFG$min_pairs)))
    if (length(setdiff(beta_groups, keep))) message("    dropped (< ", CFG$min_pairs, " pairs): ",
                                                    paste(setdiff(beta_groups, keep), collapse = ", "))
    P <- P[P$group %in% keep, ]
    P$group <- factor(P$group, levels = keep)
    P$groupO <- factor(P$group, levels = keep, ordered = TRUE)        # reference = first level
    P$log10_geo <- log10(pmax(P$geo_km, 0.1))
    P
  }
  
  fit_beta <- function(P, metric, spec) {
    n <- nrow(P); P$y <- (P[[metric]] * (n - 1) + 0.5) / n          # squeeze into (0, 1)
    sm <- unlist(lapply(names(beta_vars), function(v)
      c(sprintf("s(%s, k = %d)", v, CFG$k_beta), sprintf("s(%s, by = groupO, k = %d)", v, CFG$k_beta))))
    terms <- c("group", sm)
    bam(as.formula(paste("y ~", paste(terms, collapse = " + "))), data = P,
        family = betar(link = "logit"), method = "fREML", discrete = TRUE,
        select = TRUE, nthreads = CFG$nthreads)
  }
  
  # Partial curve of one pair variable per group (others at group medians), on a grid
  # fixed once per spec/group/variable so replicates can be averaged
  predict_var <- function(m, P, v, spec) {
    bind_rows(lapply(levels(P$group), function(g) {
      sub <- P[P$group == g, ]; key <- paste(spec, g, v)
      if (is.null(grid_cache[[key]]))
        grid_cache[[key]] <- seq(quantile(sub[[v]], .01), quantile(sub[[v]], .99), length.out = 80)
      x <- grid_cache[[key]]
      nd <- data.frame(group = factor(g, levels = levels(P$group)),
                       groupO = factor(g, levels = levels(P$group), ordered = TRUE),
                       d_elev = median(sub$d_elev), d_lat = median(sub$d_lat),
                       log10_geo = median(sub$log10_geo))
      nd <- nd[rep(1, length(x)), ]; nd[[v]] <- x
      data.frame(group = g, x = x, fit = as.numeric(predict(m, nd, type = "response")))
    }))
  }
  
  # Binned means of the observed dissimilarities (15 quantile bins per group): data check
  bin_obs <- function(P, metric, v)
    P %>% group_by(group) %>% mutate(bin = ntile(.data[[v]], 15)) %>%
    group_by(group, bin) %>%
    summarise(x = mean(.data[[v]]), y = mean(.data[[metric]]), .groups = "drop")
  
  preds <- list(); obs <- list(); summ <- list(); pair_counts <- list()
  for (r in seq_len(CFG$n_rep)) {
    message("Replicate ", r, ": computing dissimilarities")
    P_all <- make_pairs(r)
    if (is.null(P_all)) { message("  no pairs drawn: skipping replicate"); next }
    for (spec in names(beta_specs)) {
      P <- prep_pairs(P_all, spec)
      pair_counts[[paste(r, spec)]] <- P %>% count(group) %>% mutate(rep = r, spec = spec)
      for (metric in names(beta_metrics)) {
        message("  fitting ", spec, " / ", metric)
        m <- tryCatch(fit_beta(P, metric, spec), error = function(e) { message("  FAILED: ", conditionMessage(e)); NULL })
        if (is.null(m)) next
        key <- paste(r, spec, metric)
        preds[[key]] <- bind_rows(lapply(names(beta_vars), function(v)
          cbind(rep = r, spec = spec, metric = metric, variable = v, predict_var(m, P, v, spec))))
        obs[[key]] <- bind_rows(lapply(names(beta_vars), function(v)
          cbind(rep = r, spec = spec, metric = metric, variable = v, bin_obs(P, metric, v))))
        st <- as.data.frame(summary(m)$s.table); st$term <- rownames(st)
        summ[[key]] <- cbind(rep = r, spec = spec, metric = metric, dev_expl = summary(m)$dev.expl,
                             n_pairs = nrow(P), st)
      }
    }
  }
  preds <- bind_rows(preds); obs <- bind_rows(obs)
  write.csv(bind_rows(summ), out_b("models_summary.csv"), row.names = FALSE)   # p-values: descriptive only
  write.csv(bind_rows(pair_counts), out_b("pair_counts.csv"), row.names = FALSE)
  write.csv(preds, out_b("predicted_curves.csv"), row.names = FALSE)
  
  ## ---- B4. Beta-diversity figures ----------------------------------------------------------------------------------------
  # One panel: one metric against one pair variable (line = median, ribbon = 95% interval across replicates,
  # points = binned observed means)
  beta_panel <- function(mt, vv, sp, tag) {
    sel_rows <- function(d) d[d$spec == sp & d$metric == mt & d$variable == vv, ]
    pr <- sel_rows(preds); ob <- sel_rows(obs)
    # median curve and 95% interval (2.5-97.5th percentile) across the replicate curves
    pm <- pr %>% group_by(group, x) %>%
      summarise(mid = median(fit), lo = quantile(fit, 0.025), hi = quantile(fit, 0.975), .groups = "drop")
    # binned observed means, averaged over replicates (one point per group and bin)
    ob <- ob %>% group_by(group, bin) %>% summarise(x = mean(x), y = mean(y), .groups = "drop")
    for (nm in c("ob", "pm")) {
      d <- get(nm); d$group <- factor(d$group, levels = beta_groups); assign(nm, d)
    }
    ggplot() +
      geom_ribbon(data = pm, aes(x, ymin = lo, ymax = hi, fill = group), alpha = 0.25) +
      geom_point(data = ob, aes(x, y, colour = group), size = 0.9, alpha = 0.6, shape = 16) +
      geom_line(data = pm, aes(x, mid, colour = group), linewidth = 1.2) +
      scale_colour_manual(values = CFG$cols_beta, name = NULL) +
      scale_fill_manual(values = CFG$cols_beta, guide = "none") +
      guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1, linewidth = 1.2))) +
      labs(tag = tag, x = beta_vars[[vv]], y = paste(beta_metrics[[mt]], "dissimilarity")) +
      theme_fig
  }
  
  # Full grids (all metrics x all pair variables, one file per pair set), kept in beta/
  for (sp in names(beta_specs)) {
    pr <- preds %>% filter(spec == sp); ob <- obs %>% filter(spec == sp)
    if (!nrow(pr)) next
    lab <- function(d) {
      d$metric   <- factor(beta_metrics[d$metric], levels = beta_metrics)
      d$variable <- factor(beta_vars[d$variable], levels = beta_vars)
      d$group    <- factor(d$group, levels = beta_groups)
      d
    }
    pr <- lab(pr); ob <- lab(ob)
    pm <- pr %>% group_by(metric, variable, group, x) %>%
      summarise(mid = median(fit), lo = quantile(fit, 0.025), hi = quantile(fit, 0.975), .groups = "drop")
    ob <- ob %>% group_by(metric, variable, group, bin) %>%
      summarise(x = mean(x), y = mean(y), .groups = "drop")
    g <- ggplot() +
      geom_ribbon(data = pm, aes(x, ymin = lo, ymax = hi, fill = group), alpha = 0.25) +
      geom_point(data = ob, aes(x, y, colour = group), size = 0.9, alpha = 0.6, shape = 16) +
      geom_line(data = pm, aes(x, mid, colour = group), linewidth = 1.2) +
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
    "and in latitude between the two samples. Beta-regression GAM partial curves: lines = median and ribbons =",
    "95% interval (2.5-97.5th percentile) across %d replicate subsamples (each: %d random samples per group,",
    "%d random pairwise distances per comparison type); points = binned observed means. Pairs reuse samples,",
    "so the ribbon reflects subsampling variability, not a model-based confidence interval."),
    beta_specs[[CFG$main_beta_spec]], CFG$n_rep, CFG$n_per_group, CFG$n_pairs)
  
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
    "A, B: alpha diversity along elevation (A) and latitude (B); GAM partial effects, solid = net of space",
    "(95% CI), dashed = no spatial term; p = test that the alpine gradient differs in shape from the lowland",
    "one (net model). C, D: Bray-Curtis dissimilarity between pairs of samples (within lowland, within alpine, between alpine and lowland) against the",
    sprintf("difference in elevation (C) and latitude (D); lines = median, ribbons = 95%% interval across %d replicate subsamples",
            CFG$n_rep),
    sprintf("(each: %d random samples per group, %d random pairwise distances per comparison type), points = binned observed means.",
            CFG$n_per_group, CFG$n_pairs),
    "Jaccard: Supplementary Figure 1.")
  assemble_fig(file.path(out_folder, "fig2_alpha_beta_gradients.png"),
               list(list(pA, pB), list(pC, pD)), list(get_legend(pA), get_legend(pC)), main_caption,
               width = 14, height = 10.5, row_heights = c(1.1, 1))
  saveRDS(list(A = pA, B = pB, C = pC, D = pD), file.path(out_folder, "fig2_panels.rds"))
} else message("Main figure skipped: it needs RUN_ALPHA and RUN_BETA both TRUE.")

message("Done. Outputs in ", out_folder)
writeLines(capture.output(sessionInfo()), file.path(out_folder, "sessionInfo.txt"))