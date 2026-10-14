#!/usr/bin/env Rscript
# =============================================================================
# 1_alpha_overview.R  --  Figure 1
#   A (top, full width) : world map of all samples, lowland (green) / alpine (blue; drawn on top)
#   B (middle left)     : observed OTUs vs sequencing depth, GAM fit per group, adjusted for
#                         space (sos spline); tests of alpine vs lowland (level and curve shape)
#   C (middle right)    : alpine vs lowland boxplots of Shannon, total reads (log10) and the
#                         proportion of unmapped reads, with Wilcoxon rank-sum tests
#   D (bottom left)     : richness and Shannon vs the proportion of unmapped reads (GAM, per group)
#   E (bottom right)    : alpine-vs-lowland difference in richness and Shannon with and without
#                         adjusting for the unmapped proportion: does it explain the discrepancy?
#
# Input : <data_folder>/sample_data_filtered.tsv  (from microbeatlas_filter_data.R)
# Output: <out_folder>/fig1_alpha_overview.png, model reports, test tables, panel objects
# Usage : Rscript 1_alpha_overview.R [sample_table] [out_folder]
#         Snakemake passes the exact paths (all file names are defined in the Snakefile only)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mgcv)
  library(maps)
  library(grid)
})

## ---- 0. Folders and configuration --------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
# Snakemake passes: args[1] = sample table (TSV), args[2] = output folder. The defaults are for running by hand.
data_folder  <- "/home/renku/work/epfl-altshuler/GlobalAlpineSoilMicrobiome/data"
sample_table <- if (length(args) >= 1) as.character(args[1]) else file.path(data_folder, "sample_data_filtered.tsv")
out_folder   <- if (length(args) >= 2) as.character(args[2]) else file.path(dirname(data_folder), "results", "1_alpha_overview")
dir.create(out_folder, recursive = TRUE, showWarnings = FALSE)

CFG <- list(
  alpha_file  = sample_table,
  map_ylim    = c(-60, 83),     # latitude window of the map (Antarctica is cut by default)
  pt_size     = c(Lowland = 0.5,  Alpine = 0.8),    # map point sizes
  pt_alpha    = c(Lowland = 0.35, Alpine = 0.8),    # map point transparencies
  depth_xmax  = 5e5,            # panel B: x axis is cut here (reads)
  log_y       = TRUE,           # panel B: log10 y axis for observed OTUs
  k_space     = 60,             # basis dimension of the sos spatial spline
  k_smooth    = 5,              # basis dimension of every 1-D smooth: small = smooth curves
  bs_smooth   = "ts",           # thin plate spline WITH shrinkage (the null space is penalised too)
  select      = FALSE,          # TRUE would add a second shrinkage penalty on top of "ts"
  n_plot_lowland = 20000,       # lowland points drawn behind the GAM curves (plot only)
  fig_width   = 13,             # inches
  row_heights = c(6, 4.2, 3.8), # inches: map | B, C | D, E
  nthreads    = 4,
  cols        = c(Lowland = "#7C9885", Alpine = "#2B6CB0")   # greyish green vs blue
)
set.seed(1)

## ---- 1. Data -------------------------------------------------------------------
# The sample table is a TSV without quoting (its text fields contain commas): read it with quote = "" so that
# a stray double quote cannot swallow the following lines; fill = FALSE makes a misaligned row an error.
read_sample_table <- function(path)
  read.delim(path, quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE, fill = FALSE)

raw <- read_sample_table(CFG$alpha_file)
required <- c("Longitude", "Latitude", "Alpine", "alpha_total_n_reads", "alpha_n_otus",
              "alpha_shannon", "total_raw_counts", "unmapped_raw_counts")
missing_cols <- setdiff(required, names(raw))
if (length(missing_cols)) stop("Missing columns in ", CFG$alpha_file, ": ",
                               paste(missing_cols, collapse = ", "))

d <- raw %>%
  mutate(
    Alpine        = factor(Alpine, levels = c("Lowland", "Alpine")),
    lon           = Longitude,
    lat           = Latitude,
    reads         = alpha_total_n_reads,       # mapped reads = depth behind the alpha metrics
    log10_reads   = log10(reads),
    richness      = alpha_n_otus,
    shannon       = alpha_shannon,
    unmapped_prop = unmapped_raw_counts / total_raw_counts
  ) %>%
  filter(is.finite(lon), is.finite(lat), abs(lon) <= 180, abs(lat) <= 90,
         is.finite(reads), reads > 0, is.finite(richness), !is.na(Alpine))
d$AlpineO <- factor(d$Alpine, levels = c("Lowland", "Alpine"), ordered = TRUE)   # for difference smooths

n_grp <- table(d$Alpine)
message(sprintf("%d samples: %d lowland, %d alpine", nrow(d), n_grp["Lowland"], n_grp["Alpine"]))
grp_labels <- setNames(sprintf("%s (n = %s)", names(n_grp), format(as.integer(n_grp), big.mark = ",", trim = TRUE)),
                       names(n_grp))
fmt_p <- function(p) ifelse(p < 1e-4, "< 0.0001", sprintf("= %.4f", p))
# Significance stars: *** p < 0.001, ** p < 0.01, * p < 0.05, ns otherwise
lpad <- function(x) paste0(" ", gsub("\n", "\n ", x, fixed = TRUE))   # small left padding for annotations
stars <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

theme_fig <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "none",
        plot.tag = element_text(face = "bold", size = 16),
        strip.background = element_rect(fill = "grey95"))

## ---- 2. Panel A: world map ---------------------------------------------------------
world <- map_data("world")
n_outside <- sum(d$lat < CFG$map_ylim[1] | d$lat > CFG$map_ylim[2])
if (n_outside > 0) message(n_outside, " samples fall outside the map latitude window (CFG$map_ylim)")

pA <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group),
               fill = "grey93", colour = "grey75", linewidth = 0.15) +
  # lowland first, alpine on top so the (few) alpine points are not hidden
  geom_point(data = filter(d, Alpine == "Lowland"), aes(lon, lat, colour = Alpine),
             size = CFG$pt_size[["Lowland"]], alpha = CFG$pt_alpha[["Lowland"]], shape = 16) +
  geom_point(data = filter(d, Alpine == "Alpine"), aes(lon, lat, colour = Alpine),
             size = CFG$pt_size[["Alpine"]], alpha = CFG$pt_alpha[["Alpine"]], shape = 16) +
  scale_colour_manual(values = CFG$cols, labels = grp_labels, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  # not coord_quickmap: the panel fills the full width, so A, B and C share the same edges
  coord_cartesian(xlim = c(-180, 180), ylim = CFG$map_ylim, expand = FALSE) +
  labs(tag = "A", x = "Longitude", y = "Latitude") +
  theme_fig + theme(legend.position = "bottom")

## ---- 3. Models (panels B, D, E) ---------------------------------------------------------
# Modelling data: finite unmapped proportion, positive richness and Shannon
dm <- d %>% filter(is.finite(unmapped_prop), richness > 0, shannon > 0)
message(nrow(dm), " samples in the models")

fam <- list(richness = nb(), shannon = Gamma(link = "log"))
resp_lab <- c(richness = "Observed OTUs", shannon = "Shannon index")

# Lowland = reference level: the by-smooths (ordered factor) are alpine MINUS lowland, so their
# p-values test whether the alpine curve differs in shape from the lowland one
sm1 <- function(v, by = NULL) sprintf("s(%s, bs = '%s', k = %d%s)", v, CFG$bs_smooth, CFG$k_smooth,
                                      if (is.null(by)) "" else paste0(", by = ", by))
base_rhs <- paste("Alpine +", sm1("log10_reads"), "+", sm1("log10_reads", "AlpineO"),
                  sprintf("+ s(lat, lon, bs = 'sos', k = %d)", CFG$k_space))
rhs_u <- c(m0  = base_rhs,                                                    # space + depth
           m1  = paste(base_rhs, "+", sm1("unmapped_prop")),                  # + unmapped (common)
           m1g = paste(base_rhs, "+", sm1("unmapped_prop"),                   # + unmapped (per group)
                       "+", sm1("unmapped_prop", "AlpineO")))

fit_gam <- function(resp, rhs)
  bam(as.formula(paste(resp, "~", rhs)), data = dm, family = fam[[resp]],
      method = "fREML", discrete = TRUE, select = CFG$select, nthreads = CFG$nthreads)

fits_u <- list()
for (resp in names(fam)) for (mm in names(rhs_u)) {
  key <- paste(resp, mm, sep = "__"); message("Fitting ", key)
  fits_u[[key]] <- fit_gam(resp, rhs_u[[mm]])
}
saveRDS(fits_u, file.path(out_folder, "fig1_gam_fits.rds"))
writeLines(unlist(lapply(names(fits_u), function(k)
  c(paste0("==== ", k, ": ", strsplit(k, "__")[[1]][1], " ~ ", rhs_u[[strsplit(k, "__")[[1]][2]]]),
    capture.output(print(summary(fits_u[[k]]))), "", capture.output(print(k.check(fits_u[[k]]))), ""))),
  file.path(out_folder, "fig1_gam_summaries.txt"))

# tests of alpine vs lowland: level (parametric term) and curve shape (difference smooth)
gam_tests <- function(m, term = "log10_reads") {
  sm <- summary(m); pt <- sm$p.table; st <- sm$s.table
  i <- which(grepl("AlpineO", rownames(st)) & grepl(term, rownames(st), fixed = TRUE))
  c(level_p = pt["AlpineAlpine", ncol(pt)], shape_p = if (length(i)) st[i[1], ncol(st)] else NA_real_)
}

## ---- 4. Panel B: richness vs sequencing depth ------------------------------------------
m0 <- fits_u[["richness__m0"]]
tB <- gam_tests(m0)
annB <- data.frame(label = sprintf("Alpine vs lowland\nLevel: %s (p %s)\nDepth-curve shape: %s (p %s)",
                                   stars(tB[["level_p"]]), fmt_p(tB[["level_p"]]),
                                   stars(tB[["shape_p"]]), fmt_p(tB[["shape_p"]])))

grid_depth <- bind_rows(lapply(levels(dm$Alpine), function(a) {
  x <- dm$log10_reads[dm$Alpine == a]
  data.frame(Alpine  = factor(a, levels = levels(dm$Alpine)),
             AlpineO = factor(a, levels = levels(dm$Alpine), ordered = TRUE),
             log10_reads = seq(quantile(x, .01), quantile(x, .99), length.out = 100))
}))
grid_depth$lat <- median(dm$lat); grid_depth$lon <- median(dm$lon)    # spatial term is excluded
pr <- predict(m0, grid_depth, type = "link", se.fit = TRUE, exclude = "s(lat,lon)",
              newdata.guaranteed = TRUE, discrete = FALSE)
li <- family(m0)$linkinv
grid_depth$reads <- 10^grid_depth$log10_reads
grid_depth$fit <- li(pr$fit)
grid_depth$lo  <- li(pr$fit - 1.96 * pr$se.fit)
grid_depth$hi  <- li(pr$fit + 1.96 * pr$se.fit)

# background points (all alpine, a random subset of lowland: for drawing only)
low_idx <- which(dm$Alpine == "Lowland")
bg <- bind_rows(dm[dm$Alpine == "Alpine", ],
                dm[sample(low_idx, min(CFG$n_plot_lowland, length(low_idx))), ])

pB <- ggplot() +
  geom_point(data = bg, aes(reads, richness, colour = Alpine), size = 0.4, alpha = 0.15, shape = 16) +
  geom_ribbon(data = grid_depth, aes(reads, ymin = lo, ymax = hi, fill = Alpine), alpha = 0.3) +
  geom_line(data = grid_depth, aes(reads, fit, colour = Alpine), linewidth = 1.1) +
  geom_text(data = annB, aes(label = lpad(label)), x = -Inf, y = Inf, hjust = 0, vjust = 1.3,
            size = 3.2, lineheight = 0.95, inherit.aes = FALSE) +
  scale_colour_manual(values = CFG$cols, name = NULL) +
  scale_fill_manual(values = CFG$cols, guide = "none") +
  scale_x_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  (if (CFG$log_y) scale_y_log10(breaks = c(10, 30, 100, 300, 1000, 3000)) else NULL) +
  coord_cartesian(xlim = c(900, CFG$depth_xmax), ylim = quantile(dm$richness, c(0.001, 0.999))) +
  labs(tag = "B", x = "Sequencing depth (mapped reads)", y = "Observed OTUs",
       subtitle = "GAM fit (95% CI), adjusted for space") +
  theme_fig

## ---- 5. Panel C: boxplots and Wilcoxon tests ------------------------------------------------
metric_lab <- c(shannon = "Shannon index", reads = "Total reads (log10)",
                unmapped = "Unmapped reads (proportion)")
long <- d %>%
  transmute(Alpine,
            shannon  = alpha_shannon,
            reads    = log10(total_raw_counts),
            unmapped = unmapped_prop) %>%
  pivot_longer(-Alpine, names_to = "metric", values_to = "value") %>%
  filter(is.finite(value)) %>%
  mutate(metric = factor(metric_lab[metric], levels = metric_lab))

# Wilcoxon rank-sum test; Cliff's delta = 2W/(n1 n2) - 1 (> 0: alpine tends to be higher).
# With this many samples p-values are tiny: read the effect size.
wil <- long %>% group_by(metric) %>% group_modify(~ {
  x <- .x$value[.x$Alpine == "Alpine"]; y <- .x$value[.x$Alpine == "Lowland"]
  w <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
  data.frame(n_alpine = length(x), n_lowland = length(y), W = unname(w$statistic), p = w$p.value,
             cliffs_delta = 2 * unname(w$statistic) / (length(x) * length(y)) - 1,
             median_alpine = median(x), median_lowland = median(y))
}) %>% ungroup()
write.csv(wil, file.path(out_folder, "fig1_panelC_wilcoxon.csv"), row.names = FALSE)
print(wil)
# significance bracket above each pair of boxes: stars + effect size
brk <- long %>% group_by(metric) %>%
  summarise(ymax = max(value), r = diff(range(value)), .groups = "drop") %>%
  left_join(wil, by = "metric") %>%
  mutate(y_line = ymax + 0.03 * r, y_delta = ymax + 0.13 * r, star = stars(p),
         delta_lab = sprintf("Cliff's \u03B4 = %.2f", cliffs_delta))

pC <- ggplot(long, aes(Alpine, value, fill = Alpine)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.25, width = 0.5, linewidth = 0.4) +
  geom_segment(data = brk, aes(x = 1, xend = 2, y = y_line, yend = y_line), linewidth = 0.4, inherit.aes = FALSE) +
  geom_text(data = filter(brk, star != "ns"), aes(x = 1.5, y = y_line, label = star),
            vjust = -0.15, size = 5.5, inherit.aes = FALSE) +
  geom_text(data = filter(brk, star == "ns"), aes(x = 1.5, y = y_line, label = star),
            vjust = -0.45, size = 3.5, inherit.aes = FALSE) +
  geom_text(data = brk, aes(x = 1.5, y = y_delta, label = delta_lab), vjust = 0, size = 3, inherit.aes = FALSE) +
  facet_wrap(~metric, nrow = 1, scales = "free_y", labeller = label_wrap_gen(14)) +
  scale_fill_manual(values = CFG$cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.25))) +
  labs(tag = "C", x = NULL, y = NULL,
       subtitle = "Wilcoxon rank-sum test; \u03B4 = Cliff's delta (> 0: alpine higher)") +
  theme_fig

## ---- 6. Panel D: richness and Shannon vs the unmapped proportion ----------------------------------
# Per-group curves (m1g models), depth held at the overall median, space excluded
pred_unm <- function(resp) {
  m <- fits_u[[paste(resp, "m1g", sep = "__")]]
  bind_rows(lapply(levels(dm$Alpine), function(a) {
    x <- dm$unmapped_prop[dm$Alpine == a]
    nd <- data.frame(Alpine  = factor(a, levels = levels(dm$Alpine)),
                     AlpineO = factor(a, levels = levels(dm$Alpine), ordered = TRUE),
                     unmapped_prop = seq(quantile(x, .01), quantile(x, .99), length.out = 100),
                     log10_reads = median(dm$log10_reads), lat = median(dm$lat), lon = median(dm$lon))
    p <- predict(m, nd, type = "link", se.fit = TRUE, exclude = "s(lat,lon)", newdata.guaranteed = TRUE, discrete = FALSE)
    li <- family(m)$linkinv
    data.frame(response = resp_lab[[resp]], Alpine = nd$Alpine, x = nd$unmapped_prop,
               fit = li(p$fit), lo = li(p$fit - 1.96 * p$se.fit), hi = li(p$fit + 1.96 * p$se.fit))
  }))
}
pu <- bind_rows(lapply(names(fam), pred_unm))
pu$response <- factor(pu$response, levels = resp_lab)

# does the curve against the unmapped proportion differ between alpine and lowland? (difference smooth)
unm_tests <- bind_rows(lapply(names(fam), function(resp) {
  m <- fits_u[[paste(resp, "m1g", sep = "__")]]; st <- summary(m)$s.table
  i <- which(grepl("AlpineO", rownames(st)) & grepl("unmapped_prop", rownames(st), fixed = TRUE))
  data.frame(response = resp, unmapped_shape_diff_p = if (length(i)) st[i[1], ncol(st)] else NA_real_)
}))
annD <- data.frame(response = factor(resp_lab[unm_tests$response], levels = resp_lab),
                   label = sprintf("Curves differ between groups: %s", stars(unm_tests$unmapped_shape_diff_p)))

pD <- ggplot(pu, aes(x, fit, colour = Alpine)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = Alpine), alpha = 0.3, colour = NA) +
  geom_line(linewidth = 1.1) +
  geom_text(data = annD, aes(label = lpad(label)), x = -Inf, y = Inf, hjust = 0, vjust = 1.4,
            size = 3.2, inherit.aes = FALSE) +
  facet_wrap(~response, nrow = 1, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  scale_colour_manual(values = CFG$cols, name = NULL) +
  scale_fill_manual(values = CFG$cols, guide = "none") +
  labs(tag = "D", x = "Unmapped reads (proportion)", y = "Predicted value",
       subtitle = "GAM fits (95% CI), adjusted for space and depth") +
  theme_fig

## ---- 7. Panel E: does the unmapped proportion explain the alpine-lowland difference? ---------------------
# (a) Coefficient of the Alpine term (log scale), kept in a table for reference. With a COMMON curve
#     for the unmapped proportion (m1) it is misleading when the groups' curves differ (panel D).
coef_tab <- bind_rows(lapply(names(fam), function(resp) bind_rows(lapply(c("m0", "m1"), function(mm) {
  m <- fits_u[[paste(resp, mm, sep = "__")]]; sm <- summary(m); pt <- sm$p.table
  b <- pt["AlpineAlpine", "Estimate"]; se <- pt["AlpineAlpine", "Std. Error"]
  data.frame(response = resp, model = mm, beta = b, se = se,
             pct = 100 * (exp(b) - 1), pct_lo = 100 * (exp(b - 1.96 * se) - 1),
             pct_hi = 100 * (exp(b + 1.96 * se) - 1), p = pt["AlpineAlpine", ncol(pt)],
             depth_shape_p = gam_tests(m)[["shape_p"]], AIC = AIC(m), dev_expl = sm$dev.expl)
}))))
write.csv(coef_tab, file.path(out_folder, "fig1_panelE_alpine_coefficients.csv"), row.names = FALSE)
write.csv(unm_tests, file.path(out_folder, "fig1_panelD_unmapped_group_difference.csv"), row.names = FALSE)

# (b) Plotted: effect STANDARDISED over the alpine samples. For every alpine sample, the model
#     prediction with its own (alpine) status vs the same sample (same depth, unmapped proportion
#     and location) treated as lowland, averaged on the log scale. This stays valid when the two
#     groups have different curves, which the single Alpine coefficient does not.
alpine_att <- function(m) {
  nd1 <- dm[dm$Alpine == "Alpine", ]
  nd0 <- nd1
  nd0$Alpine  <- factor("Lowland", levels = levels(dm$Alpine))
  nd0$AlpineO <- factor("Lowland", levels = levels(dm$Alpine), ordered = TRUE)
  X1 <- predict(m, nd1, type = "lpmatrix", discrete = FALSE)
  X0 <- predict(m, nd0, type = "lpmatrix", discrete = FALSE)
  a <- colMeans(X1 - X0)                                     # the spatial terms cancel out
  est <- sum(a * coef(m)); se <- sqrt(drop(t(a) %*% vcov(m) %*% a))
  data.frame(beta = est, se = se, pct = 100 * (exp(est) - 1),
             pct_lo = 100 * (exp(est - 1.96 * se) - 1), pct_hi = 100 * (exp(est + 1.96 * se) - 1),
             p = 2 * pnorm(-abs(est / se)))
}
att_tab <- bind_rows(lapply(names(fam), function(resp) bind_rows(lapply(c("m0", "m1g"), function(mm)
  cbind(response = resp, model = mm, alpine_att(fits_u[[paste(resp, mm, sep = "__")]]))))))
write.csv(att_tab, file.path(out_folder, "fig1_panelE_alpine_effect.csv"), row.names = FALSE)
print(att_tab)

att_tab <- att_tab %>%
  mutate(star = stars(p),
         response = factor(resp_lab[response], levels = resp_lab),
         model = factor(c(m0 = "Space + depth", m1g = "+ unmapped proportion")[model],
                        levels = c("+ unmapped proportion", "Space + depth")))

pE <- ggplot(att_tab, aes(pct, model)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = pct_lo, xmax = pct_hi), size = 0.6) +
  geom_text(aes(x = pct_hi, label = star), hjust = -0.4, vjust = 0.3, size = 5) +
  facet_wrap(~response, ncol = 1, scales = "free_x") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  labs(tag = "E", x = "Alpine vs same samples as lowland (%, 95% CI)", y = NULL,
       subtitle = "Stars: differs from 0") +
  theme_fig + theme(strip.text = element_text(size = 9))

## ---- 8. Assemble and save --------------------------------------------------------------------
png(file.path(out_folder, "fig1_alpha_overview.png"), width = CFG$fig_width,
    height = sum(CFG$row_heights), units = "in", res = 300)
grid.newpage()

# align the left edge of the panel area of A, B and D (their y-axis labels differ in width)
left_mm <- function(p) {
  g <- ggplotGrob(p)
  i <- min(g$layout$l[grepl("^panel", g$layout$name)])
  convertWidth(sum(g$widths[seq_len(i - 1)]), "mm", valueOnly = TRUE)
}
L <- c(A = left_mm(pA), B = left_mm(pB), D = left_mm(pD))
pad <- function(p, delta_mm)
  p + theme(axis.title.y = element_text(margin = margin(r = 2.75 + delta_mm * 72.27 / 25.4)))
pA <- pad(pA, max(L) - L[["A"]]); pB <- pad(pB, max(L) - L[["B"]]); pD <- pad(pD, max(L) - L[["D"]])

pushViewport(viewport(layout = grid.layout(3, 2, heights = unit(CFG$row_heights, "null"),
                                           widths = unit(c(1.15, 1), "null"))))
print(pA, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
print(pB, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(pC, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(pD, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
print(pE, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
dev.off()

saveRDS(list(A = pA, B = pB, C = pC, D = pD, E = pE), file.path(out_folder, "fig1_panels.rds"))

summary_tab <- long %>% group_by(metric, Alpine) %>%
  summarise(n = n(), median = median(value), q25 = quantile(value, .25),
            q75 = quantile(value, .75), .groups = "drop")
write.csv(summary_tab, file.path(out_folder, "fig1_panelC_summary.csv"), row.names = FALSE)

message("Done. Outputs in ", out_folder)
file.create(file.path(out_folder, ".done"))   # flag file for Snakemake