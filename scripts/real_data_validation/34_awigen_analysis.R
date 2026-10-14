library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(patchwork)

repo_root <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("--file=", "", args[grep("^--file=", args)])
  dirname(dirname(dirname(normalizePath(file_arg))))
})

manifest_path <- file.path(repo_root, "results/awigen_validation_samples.tsv")
results_dir   <- file.path(repo_root, "benchmark_results/awigen_benchmark_results")
outdir_fig    <- file.path(repo_root, "figures")
outdir_stats  <- file.path(repo_root, "stats")
dir.create(outdir_fig,   showWarnings = FALSE, recursive = TRUE)
dir.create(outdir_stats, showWarnings = FALSE, recursive = TRUE)

# ── Load per-individual results ──────────────────────────────────────────────
# Each result file is a single unheadered row: metaG_id,total_reads,default_hits,african_hits
result_files <- list.files(results_dir, pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(result_files) > 0)

results <- map_dfr(result_files, function(p) {
  read_csv(p, col_names = c("metaG_id", "total_reads", "default_hits", "african_hits"),
            col_types = "cddd", show_col_types = FALSE)
})

manifest <- read_tsv(manifest_path, show_col_types = FALSE)

dat <- results %>%
  inner_join(manifest, by = "metaG_id") %>%
  mutate(
    default_rate = default_hits / total_reads,
    african_rate = african_hits / total_reads,
    gain_rate    = african_rate - default_rate,
    extra_reads  = african_hits - default_hits,
    relative_gain_pct = ifelse(default_hits > 0, extra_reads / default_hits * 100, NA_real_)
  )

cat(sprintf("Loaded %d individuals (expected 600).\n", nrow(dat)))
if (nrow(dat) != 600) cat("WARNING: individual count does not match expected 600.\n")
if (any(dat$default_rate == 0 & dat$african_rate == 0)) {
  cat("WARNING: some individuals have zero hits for both conditions - check before trusting results.\n")
}

write_csv(dat, file.path(outdir_stats, "awigen_per_individual.csv"))

# ── Paired Wilcoxon signed-rank test, per-individual hit rate ───────────────
# Rate (hits/total_reads), not raw counts - normalizes for the large spread in
# per-individual read depth (12M-68M+), same reasoning as every other
# per-individual test in this project. conf.int=TRUE gives a Hodges-Lehmann
# estimate and 95% CI on the median paired difference alongside the test.
run_wilcoxon <- function(df, label) {
  # exact = TRUE, not the normal approximation: wilcox.test's approximate
  # (uniroot-based) CI/estimate search can converge to the wrong root under
  # heavy skew in the paired differences - confirmed on this data, where the
  # exact=FALSE South Africa estimate (0.0153 pp) fell outside its own
  # reported CI (0.0104-0.0113 pp), an impossible result. The exact method
  # matches a direct manual computation of the Hodges-Lehmann estimator
  # (median of all pairwise Walsh averages) in every group here, and is
  # computationally exact anyway - n<=600 with zero tied differences in every
  # group, well within reach of the exact algorithm (no approximation needed
  # or wanted).
  wt <- wilcox.test(df$african_rate, df$default_rate, paired = TRUE,
                     conf.int = TRUE, conf.level = 0.95, exact = TRUE)
  tibble(
    group = label,
    n = nrow(df),
    median_default_pct = median(df$default_rate) * 100,
    median_african_pct  = median(df$african_rate) * 100,
    median_gain_pct_points = median(df$gain_rate) * 100,
    hodges_lehmann_estimate_pct_points = unname(wt$estimate) * 100,
    ci_lower_pct_points = wt$conf.int[1] * 100,
    ci_upper_pct_points = wt$conf.int[2] * 100,
    median_relative_gain_pct = median(df$relative_gain_pct),
    mean_relative_gain_pct = mean(df$relative_gain_pct),
    median_extra_reads = median(df$extra_reads),
    mean_extra_reads = mean(df$extra_reads),
    V_statistic = unname(wt$statistic),
    V_max_possible = nrow(df) * (nrow(df) + 1) / 2,
    p_value = wt$p.value
  )
}

overall_stats <- run_wilcoxon(dat, "Overall")

country_stats <- dat %>%
  group_by(country) %>%
  group_map(~ run_wilcoxon(.x, .y$country)) %>%
  bind_rows()

stats_table <- bind_rows(overall_stats, country_stats)

cat("\n=== Per-individual paired Wilcoxon signed-rank test (African vs Default hit rate) ===\n")
print(stats_table %>% mutate(across(where(is.numeric), ~round(.x, 6))), n = Inf, width = Inf)

write_csv(stats_table, file.path(outdir_stats, "awigen_validation_wilcoxon.csv"))
cat("\nWrote awigen_validation_wilcoxon.csv\n")

# ── Figure ────────────────────────────────────────────────────────────────────
# No panel/figure titles or subtitles - only A/B/C corner tags (via
# plot_annotation(tag_levels)). Full description (what each panel shows, the
# significance test, summary numbers) lives in the figure legend and results
# text instead, not baked into the image.
# Panel A: for each individual, the fold-change African_rate / Default_rate -
# "African DB caught Nx as many human reads as Default did, for this person."
# One histogram of all 600 individuals, coloured by country, log10 x-axis
# (fold-change ranges 1.01x-11.5x - multiplicative, not additive, so a log
# axis is the natural scale). A dashed reference line at 1x marks "no
# difference"; since every individual's fold-change is >1, the entire
# distribution sits to the right of that line - a single visual answering
# "did African DB ever do worse" (no) and "by what multiple did it typically
# do better" (median 1.25x, mean 1.75x, long right tail to >11x) at once.
# Replaces an earlier per-individual scatter plot against a 1:1 line, which
# accurately showed the same "always at or above" result but required
# reading two log-scaled percentage axes at once - a single fold-change axis
# is a more direct answer to "how much better," which is what this panel is
# for (Panels B/C already cover relative-% and absolute-count gain).
# Panel B: per-individual relative gain (%) - extra reads caught by African
# DB as a fraction of Default's own catch. Only the median is annotated
# (matches the rank-based Wilcoxon test/Hodges-Lehmann estimator used
# throughout, and is far less skew-sensitive than a mean; repeating n and p
# per box was redundant with the one-line figure-level significance
# statement).
# Panel C: per-individual absolute extra reads caught (count, not %) - the
# tangible "how many more actual reads per person" number, alongside the
# relative view in Panel B.

country_cols <- c("Burkina Faso" = "#66C2A5", "Ghana" = "#FC8D62",
                   "Kenya" = "#8DA0CB", "South Africa" = "#E78AC3")

base_theme <- theme_classic(base_size = 10.5) +
  theme(
    strip.text         = element_text(face = "bold", size = 9.5),
    strip.background   = element_rect(fill = "grey94", colour = NA),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
    plot.title         = element_text(face = "bold", size = 11.5, hjust = 0),
    plot.subtitle      = element_text(size = 8.5, colour = "grey35")
  )

# Median, not mean: matches the Wilcoxon test itself (rank-based, and its
# Hodges-Lehmann estimate is a median of Walsh averages) and is far less
# sensitive to the right-skew already visible in panels B/C than a mean
# would be.
panel_labels <- country_stats %>%
  rename(country = group) %>%
  mutate(
    gain_label  = sprintf("median +%.0f%%", median_relative_gain_pct),
    reads_label = sprintf("median %s reads", scales::comma(round(median_extra_reads)))
  )

dat <- dat %>% mutate(fold_change = african_rate / default_rate)

panel_a <- ggplot(dat, aes(x = fold_change, fill = country)) +
  geom_histogram(bins = 40, colour = "white", linewidth = 0.15) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
  annotate("text", x = 1, y = Inf, label = "no difference", angle = 90,
           vjust = -0.6, hjust = 1.1, size = 2.9, colour = "grey30") +
  # Doubling breaks (1x,2x,4x,8x) land at equal intervals on a log axis
  # (equal ratio = equal log-distance, regardless of log base) - irregular
  # breaks like 1/2/3/5/8/12 have uneven step ratios (2x, 1.5x, 1.67x, 1.6x,
  # 1.5x) and so look inconsistently spaced even though the axis is correct.
  scale_x_log10(breaks = c(1, 2, 4, 8), labels = function(x) paste0(x, "x")) +
  scale_fill_manual(values = country_cols, name = NULL) +
  labs(x = "Fold-change (African DB rate ÷ Default rate)", y = "Number of individuals") +
  base_theme +
  theme(legend.position = "bottom")

# log10 y-axis: relative gain spans ~3 orders of magnitude (1.2%-1052%,
# all strictly positive - no zero/negative values to handle) and is
# right-skewed by a handful of large-gain individuals, so a linear axis
# crushes the bulk of the distribution into an unreadable band near zero.
panel_b <- ggplot(dat, aes(x = country, y = relative_gain_pct, fill = country)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.4, colour = "grey30") +
  geom_text(data = panel_labels, aes(x = country, y = 2200, label = gain_label),
            inherit.aes = FALSE, size = 2.8, colour = "grey20") +
  scale_y_log10(labels = scales::label_percent(scale = 1, accuracy = 1), limits = c(NA, 3000)) +
  scale_fill_manual(values = country_cols, guide = "none") +
  labs(x = NULL, y = "Relative gain (%)") +
  base_theme

# extra_reads (African hits - Default hits) is strictly positive for all 600
# individuals and spans >2 orders of magnitude (1.2K-537K) - log scale again.
panel_c <- ggplot(dat, aes(x = country, y = extra_reads, fill = country)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.4, colour = "grey30") +
  geom_text(data = panel_labels, aes(x = country, y = 1.1e6, label = reads_label),
            inherit.aes = FALSE, size = 2.8, colour = "grey20") +
  scale_y_log10(labels = scales::label_comma(), limits = c(NA, 2e6)) +
  scale_fill_manual(values = country_cols, guide = "none") +
  labs(x = NULL, y = "Extra reads caught (count)") +
  base_theme

# No panel titles/subtitles or figure-level title - full description (what
# each panel shows, the significance test, the summary numbers) goes in the
# figure legend and results text instead. tag_levels = "A" auto-labels the
# three panels A/B/C in the corner, which is all the image itself carries.
combined_fig <- (panel_a / (panel_b | panel_c)) +
  plot_layout(heights = c(1.1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))

ggsave(file.path(outdir_fig, "awigen_validation.pdf"), combined_fig, width = 11, height = 8.5)
ggsave(file.path(outdir_fig, "awigen_validation.png"), combined_fig, width = 11, height = 8.5, dpi = 300)
cat("\nWrote awigen_validation\n")
