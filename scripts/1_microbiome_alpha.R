#!/usr/bin/env Rscript
# =============================================================================
# 1_alpha_overview.R  --  Figure 1
#   A (top, wide) : world map of all samples, lowland (green) / alpine (blue)
#   B (bottom L)  : observed OTUs vs sequencing depth, GAM fit per group
#                   (alpine vs lowland), accounting for space (sos spline)
#   C (bottom R)  : alpine vs lowland boxplots of Shannon, total reads (log10)
#                   and the proportion of unmapped reads
#
# Input : <data_folder>/sample_data_filtered.csv  (from microbeatlas_filter_data.R)
# Output: <out_folder>/fig1_alpha_overview.png, GAM summary, summary table,
#         ggplot objects (fig1_panels.rds) for later tweaking
# Usage : Rscript 1_alpha_overview.R [data_folder] [out_folder]
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
data_folder <- if (length(args) >= 1) as.character(args[1]) else
  "/home/renku/work/epfl-altshuler/GlobalAlpineSoilMicrobiome/data"
out_folder  <- if (length(args) >= 2) as.character(args[2]) else
  file.path(dirname(data_folder), "results", "1_alpha_overview")
dir.create(out_folder, recursive = TRUE, showWarnings = FALSE)

CFG <- list(
  alpha_file  = file.path(data_folder, "sample_data_filtered.csv"),
  map_ylim    = c(-60, 83),     # latitude window of the map (Antarctica is cut by default)
  k_space     = 60,             # basis dimension of the sos spatial spline
  n_plot_lowland = 20000,       # lowland points drawn behind the GAM curves (plot only)
  nthreads    = 4,
  cols        = c(Lowland = "#3A8F4B", Alpine = "#2B6CB0")
)
set.seed(1)

## ---- 1. Data -------------------------------------------------------------------
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
    unmapped_prop = unmapped_raw_counts / total_raw_counts
  ) %>%
  filter(is.finite(lon), is.finite(lat), abs(lon) <= 180, abs(lat) <= 90,
         is.finite(reads), reads > 0, is.finite(richness), !is.na(Alpine))

n_grp <- table(d$Alpine)
message(sprintf("%d samples: %d lowland, %d alpine", nrow(d), n_grp["Lowland"], n_grp["Alpine"]))
grp_labels <- setNames(sprintf("%s (n = %s)", names(n_grp), format(as.integer(n_grp), big.mark = ",")),
                       names(n_grp))

theme_fig <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.tag = element_text(face = "bold", size = 16))

## ---- 2. Panel A: world map ---------------------------------------------------------
world <- map_data("world")
n_outside <- sum(d$lat < CFG$map_ylim[1] | d$lat > CFG$map_ylim[2])
if (n_outside > 0) message(n_outside, " samples fall outside the map latitude window (CFG$map_ylim)")

pA <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group),
               fill = "grey93", colour = "grey75", linewidth = 0.15) +
  # lowland first, alpine on top so the (few) alpine points are not hidden
  geom_point(data = filter(d, Alpine == "Lowland"), aes(lon, lat, colour = Alpine),
             size = 0.25, alpha = 0.35, shape = 16) +
  geom_point(data = filter(d, Alpine == "Alpine"), aes(lon, lat, colour = Alpine),
             size = 0.4, alpha = 0.8, shape = 16) +
  scale_colour_manual(values = CFG$cols, labels = grp_labels, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  coord_quickmap(xlim = c(-180, 180), ylim = CFG$map_ylim, expand = FALSE) +
  labs(tag = "A", x = "Longitude", y = "Latitude") +
  theme_fig + theme(legend.position = "bottom")

## ---- 3. Panel B: GAM of richness vs sequencing depth ------------------------------------
dm <- d

rhs <- sprintf("Alpine + s(log10_reads, by = Alpine, k = 8) + s(lat, lon, bs = 'sos', k = %d)", CFG$k_space)
m <- bam(as.formula(paste("richness ~", rhs)), data = dm, family = nb(),
         method = "fREML", discrete = TRUE, select = TRUE, nthreads = CFG$nthreads)
writeLines(c(paste("richness ~", rhs), "", capture.output(print(summary(m))),
             "", capture.output(print(k.check(m)))),
           file.path(out_folder, "fig1_panelB_gam_summary.txt"))

# Predictions along sequencing depth (1-99% range per group); space held at the
# median location (spatial term excluded)
grid_depth <- bind_rows(lapply(levels(dm$Alpine), function(a) {
  x <- dm$log10_reads[dm$Alpine == a]
  data.frame(Alpine = factor(a, levels = levels(dm$Alpine)),
             log10_reads = seq(quantile(x, .01), quantile(x, .99), length.out = 100))
}))
grid_depth$lat <- median(dm$lat); grid_depth$lon <- median(dm$lon)

pr <- predict(m, grid_depth, type = "link", se.fit = TRUE, exclude = "s(lat,lon)",
              newdata.guaranteed = TRUE)
li <- family(m)$linkinv
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
  scale_colour_manual(values = CFG$cols, labels = grp_labels, name = NULL) +
  scale_fill_manual(values = CFG$cols, guide = "none") +
  scale_x_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  coord_cartesian(ylim = c(0, quantile(dm$richness, .995))) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(tag = "B", x = "Sequencing depth (mapped reads)", y = "Observed OTUs",
       subtitle = "GAM fit (95% CI), adjusted for space") +
  theme_fig + theme(legend.position = "bottom")

## ---- 4. Panel C: boxplots ---------------------------------------------------------------
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

pC <- ggplot(long, aes(Alpine, value, fill = Alpine)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.25, width = 0.6, linewidth = 0.4) +
  facet_wrap(~metric, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = CFG$cols, guide = "none") +
  labs(tag = "C", x = NULL, y = NULL) +
  theme_fig + theme(strip.background = element_rect(fill = "grey95"))

## ---- 5. Assemble and save -------------------------------------------------------------------
png(file.path(out_folder, "fig1_alpha_overview.png"), width = 13, height = 10, units = "in", res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2,
                                           heights = unit(c(1.1, 1), "null"),
                                           widths  = unit(c(1, 1.35), "null"))))
print(pA, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
print(pB, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(pC, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()

saveRDS(list(A = pA, B = pB, C = pC), file.path(out_folder, "fig1_panels.rds"))

summary_tab <- long %>% group_by(metric, Alpine) %>%
  summarise(n = n(), median = median(value), q25 = quantile(value, .25),
            q75 = quantile(value, .75), .groups = "drop")
write.csv(summary_tab, file.path(out_folder, "fig1_panelC_summary.csv"), row.names = FALSE)
print(summary_tab)

message("Done. Outputs in ", out_folder)