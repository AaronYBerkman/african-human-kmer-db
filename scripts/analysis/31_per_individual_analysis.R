library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(patchwork)

repo_root <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("--file=", "", args[grep("^--file=", args)])
  dirname(dirname(dirname(normalizePath(file_arg))))
})

per_indiv_path       <- file.path(repo_root, "benchmark_results/per_individual_results/per_individual_summary.csv")
per_indiv_afr25_path <- file.path(repo_root, "benchmark_results/per_individual_african25_results/per_individual_african25_summary.csv")
outdir   <- file.path(repo_root, "figures")
statsdir <- file.path(repo_root, "stats")
dir.create(outdir,   showWarnings = FALSE, recursive = TRUE)
dir.create(statsdir, showWarnings = FALSE, recursive = TRUE)

raw <- read_csv(per_indiv_path, show_col_types = FALSE)

# Cap sensitivity at 1.0 — aligns_to can report duplicate rows for a read
# that maps to multiple DB entries; the raw line-count then exceeds n_reads.
# The correct treatment is deduplication at the bash level (fixed in script 38
# via sort -u); here we cap the rate so one duplicated individual doesn't
# distort medians or tests, and flag it in the annotated CSV.
raw <- raw %>%
  mutate(
    dup_hit_flag        = african_sensitivity > 1 | combined_sensitivity > 1,
    african_sensitivity  = pmin(african_sensitivity,  1),
    combined_sensitivity = pmin(combined_sensitivity, 1)
  )
if (any(raw$dup_hit_flag))
  cat("NOTE: capped sensitivity at 1.0 for",
      paste(raw$run_accession[raw$dup_hit_flag], collapse = ", "),
      "(duplicate aligns_to rows — rerun script 38 with sort -u to fix properly)\n\n")

# ── Population annotation ────────────────────────────────────────────────────
superpop <- tribble(
  ~population_code, ~superpopulation, ~ancestry_group,
  "YRI", "AFR", "Continental-African",
  "LWK", "AFR", "Continental-African",
  "GWD", "AFR", "Continental-African",
  "MSL", "AFR", "Continental-African",
  "ESN", "AFR", "Continental-African",
  "ASW", "AFR", "African diaspora",
  "ACB", "AFR", "African diaspora",
  "CEU", "EUR", "Non-African",
  "TSI", "EUR", "Non-African",
  "FIN", "EUR", "Non-African",
  "GBR", "EUR", "Non-African",
  "IBS", "EUR", "Non-African",
  "CHB", "EAS", "Non-African",
  "JPT", "EAS", "Non-African",
  "CHS", "EAS", "Non-African",
  "CDX", "EAS", "Non-African",
  "KHV", "EAS", "Non-African",
  "GIH", "SAS", "Non-African",
  "PJL", "SAS", "Non-African",
  "BEB", "SAS", "Non-African",
  "STU", "SAS", "Non-African",
  "ITU", "SAS", "Non-African",
  "MXL", "AMR", "Non-African",
  "PUR", "AMR", "Non-African",
  "CLM", "AMR", "Non-African",
  "PEL", "AMR", "Non-African"
)

dat <- raw %>%
  left_join(superpop, by = "population_code") %>%
  mutate(
    gain_african  = african_sensitivity  - default_sensitivity,
    gain_combined = combined_sensitivity - default_sensitivity,
    ancestry_group = factor(ancestry_group,
      levels = c("Continental-African", "African diaspora", "Non-African"))
  )

write_csv(dat, file.path(statsdir, "per_individual_summary_annotated.csv"))
cat("=== Per-individual summary ===\n")
print(dat %>% select(population_code, run_accession, default_sensitivity,
                     african_sensitivity, combined_sensitivity, gain_african) %>%
        mutate(across(where(is.numeric), ~round(.x, 4))), n = 27)

# ── Fix 1: Wilcoxon signed-rank (paired, n=27) ──────────────────────────────
wt_african  <- wilcox.test(dat$african_sensitivity,  dat$default_sensitivity, paired = TRUE, exact = FALSE)
wt_combined <- wilcox.test(dat$combined_sensitivity, dat$default_sensitivity, paired = TRUE, exact = FALSE)

med_default  <- median(dat$default_sensitivity)
med_african  <- median(dat$african_sensitivity)
med_combined <- median(dat$combined_sensitivity)

cat(sprintf("\n=== Fix 1: Wilcoxon signed-rank (n=%d paired individuals) ===\n", nrow(dat)))
cat(sprintf("Median: Default=%.2f%%  African DB=%.2f%%  Combined=%.2f%%\n",
            med_default*100, med_african*100, med_combined*100))
cat(sprintf("African DB vs Default:  V=%g, p=%.4g\n", wt_african$statistic,  wt_african$p.value))
cat(sprintf("Combined vs Default:    V=%g, p=%.4g\n", wt_combined$statistic, wt_combined$p.value))

forbes27_wilcoxon <- tibble(
  panel          = "Global",
  comparison     = c("African DB vs Default", "Combined vs Default"),
  n_pairs        = nrow(dat),
  median_default = med_default,
  median_method  = c(med_african, med_combined),
  median_gain_pp = c(med_african - med_default, med_combined - med_default) * 100,
  V_max_possible = nrow(dat) * (nrow(dat) + 1) / 2,
  V_statistic    = c(wt_african$statistic, wt_combined$statistic),
  p_value        = c(wt_african$p.value,   wt_combined$p.value)
)
write_csv(forbes27_wilcoxon, file.path(statsdir, "wilcoxon_per_individual_forbes27.csv"))

# ── African-25 per-individual Wilcoxon ───────────────────────────────────────
afr25_file_exists <- file.exists(per_indiv_afr25_path)
if (afr25_file_exists) {
  raw25 <- read_csv(per_indiv_afr25_path, show_col_types = FALSE)
  raw25 <- raw25 %>%
    mutate(
      dup_hit_flag        = african_sensitivity > 1 | combined_sensitivity > 1,
      african_sensitivity  = pmin(african_sensitivity,  1),
      combined_sensitivity = pmin(combined_sensitivity, 1)
    )
  if (any(raw25$dup_hit_flag))
    cat("NOTE (African-25): capped sensitivity at 1.0 for",
        paste(raw25$run_accession[raw25$dup_hit_flag], collapse = ", "), "\n\n")

  wt25_african  <- wilcox.test(raw25$african_sensitivity,  raw25$default_sensitivity, paired = TRUE, exact = FALSE)
  wt25_combined <- wilcox.test(raw25$combined_sensitivity, raw25$default_sensitivity, paired = TRUE, exact = FALSE)

  med25_default  <- median(raw25$default_sensitivity)
  med25_african  <- median(raw25$african_sensitivity)
  med25_combined <- median(raw25$combined_sensitivity)

  cat(sprintf("\n=== African-25 per-individual Wilcoxon (n=%d) ===\n", nrow(raw25)))
  cat(sprintf("Median: Default=%.2f%%  African DB=%.2f%%  Combined=%.2f%%\n",
              med25_default*100, med25_african*100, med25_combined*100))
  cat(sprintf("African DB vs Default:  V=%g (max=%g), p=%.4g\n",
              wt25_african$statistic,  nrow(raw25)*(nrow(raw25)+1)/2, wt25_african$p.value))
  cat(sprintf("Combined vs Default:    V=%g (max=%g), p=%.4g\n",
              wt25_combined$statistic, nrow(raw25)*(nrow(raw25)+1)/2, wt25_combined$p.value))

  afr25_wilcoxon <- tibble(
    panel          = "Continental-African",
    comparison     = c("African DB vs Default", "Combined vs Default"),
    n_pairs        = nrow(raw25),
    median_default = med25_default,
    median_method  = c(med25_african, med25_combined),
    median_gain_pp = c(med25_african - med25_default, med25_combined - med25_default) * 100,
    V_max_possible = nrow(raw25) * (nrow(raw25) + 1) / 2,
    V_statistic    = c(wt25_african$statistic, wt25_combined$statistic),
    p_value        = c(wt25_african$p.value,   wt25_combined$p.value)
  )

  supp_perind <- bind_rows(forbes27_wilcoxon, afr25_wilcoxon)
  cat("\n=== tab:supp-perind — combined supplementary table ===\n")
  print(supp_perind %>% mutate(across(where(is.numeric), ~round(.x, 4))), n = Inf)
  write_csv(supp_perind, file.path(statsdir, "supp_tab_per_individual_wilcoxon.csv"))
  cat("\nWrote supp_tab_per_individual_wilcoxon.csv\n")
} else {
  cat("\nNOTE: African-25 per-individual summary not yet available —\n",
      "run script 40 first, then re-run this script to get the combined table.\n",
      "Forbes-27 results written to wilcoxon_per_individual_forbes27.csv\n\n")
  write_csv(forbes27_wilcoxon, file.path(statsdir, "supp_tab_per_individual_wilcoxon.csv"))
}

# ── Fix 2: within-panel ancestry comparison ──────────────────────────────────
gain_cont_afr <- dat %>% filter(ancestry_group == "Continental-African") %>% pull(gain_african)
gain_diaspora <- dat %>% filter(ancestry_group == "African diaspora")    %>% pull(gain_african)
gain_non_afr  <- dat %>% filter(ancestry_group == "Non-African")         %>% pull(gain_african)

mw_cont_vs_non   <- wilcox.test(gain_cont_afr, gain_non_afr,  exact = FALSE)
mw_diasp_vs_non  <- wilcox.test(gain_diaspora, gain_non_afr,  exact = FALSE)
mw_cont_vs_diasp <- wilcox.test(gain_cont_afr, gain_diaspora, exact = FALSE)

cat("\n=== Fix 2: within-panel ancestry gain comparison ===\n")
dat %>%
  group_by(ancestry_group) %>%
  summarise(n = n(), median_gain_pp = round(median(gain_african)*100, 2),
            iqr_pp = round(IQR(gain_african)*100, 2), .groups = "drop") %>%
  print()
cat(sprintf("Continental-African vs Non-African: W=%g, p=%.4g\n",
            mw_cont_vs_non$statistic, mw_cont_vs_non$p.value))
cat(sprintf("Diaspora vs Non-African:            W=%g, p=%.4g\n",
            mw_diasp_vs_non$statistic, mw_diasp_vs_non$p.value))

ancestry_stats <- dat %>%
  group_by(ancestry_group) %>%
  summarise(n = n(),
            median_gain_pp = median(gain_african)*100,
            iqr_gain_pp    = IQR(gain_african)*100,
            median_default  = median(default_sensitivity),
            median_african  = median(african_sensitivity),
            median_combined = median(combined_sensitivity),
            .groups = "drop")

write_csv(ancestry_stats, file.path(statsdir, "ancestry_gain_summary.csv"))
write_csv(
  tibble(comparison = c("Continental-African vs Non-African",
                        "African diaspora vs Non-African",
                        "Continental-African vs African diaspora"),
         W_statistic = c(mw_cont_vs_non$statistic, mw_diasp_vs_non$statistic, mw_cont_vs_diasp$statistic),
         p_value     = c(mw_cont_vs_non$p.value,   mw_diasp_vs_non$p.value,   mw_cont_vs_diasp$p.value)),
  file.path(statsdir, "ancestry_mannwhitney.csv")
)

# ── Shared theme ─────────────────────────────────────────────────────────────
group_cols <- c("Continental-African" = "#7B3294",
                "African diaspora"    = "#C2A5CF",
                "Non-African"         = "#A6DBA0")

method_cols <- c("Default (Human RefSeq)" = "#E64B35",
                 "African DB only"         = "#4DBBD5",
                 "Combined (Stage 1+2)"    = "#00A087")

base_theme <- theme_classic(base_size = 10.5) +
  theme(
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
    legend.text        = element_text(size = 10),
    plot.title         = element_text(face = "bold", size = 11.5, hjust = 0),
    plot.subtitle      = element_text(size = 8.5, colour = "grey35", hjust = 0,
                                      margin = margin(b = 8))
  )

# ── Figure 3: per-individual sensitivity (Fix 1) ─────────────────────────────
long_sens <- dat %>%
  pivot_longer(cols = c(default_sensitivity, african_sensitivity, combined_sensitivity),
               names_to  = "method_key", values_to = "sensitivity") %>%
  mutate(
    method = recode(method_key,
      default_sensitivity  = "Default (Human RefSeq)",
      african_sensitivity  = "African DB only",
      combined_sensitivity = "Combined (Stage 1+2)"),
    method = factor(method, levels = names(method_cols))
  )

medians3 <- long_sens %>%
  group_by(method) %>%
  summarise(med = median(sensitivity, na.rm = TRUE), .groups = "drop")

y_lo <- floor(min(long_sens$sensitivity, na.rm = TRUE) * 1000) / 1000
y_hi <- 1.002

sub3 <- sprintf(
  "Wilcoxon signed-rank, n=%d paired individuals  |  African DB vs Default: p = %.2g  |  Combined vs Default: p = %.2g",
  nrow(dat), wt_african$p.value, wt_combined$p.value)

set.seed(42)
fig3 <- ggplot(long_sens, aes(x = method, y = sensitivity)) +
  geom_line(aes(group = run_accession), colour = "grey78", linewidth = 0.35) +
  geom_point(aes(colour = ancestry_group),
             position = position_jitter(width = 0.07), size = 2.7) +
  geom_crossbar(data = medians3,
                aes(x = method, y = med, ymin = med, ymax = med),
                inherit.aes = FALSE, colour = "grey15",
                width = 0.5, linewidth = 0.75) +
  scale_colour_manual(values = group_cols, name = "Ancestry group") +
  scale_x_discrete(labels = c("Default\n(Human RefSeq)", "African\nDB only",
                               "Combined\n(Stage 1+2)")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                     limits = c(y_lo, y_hi),
                     expand = expansion(mult = c(0.02, 0.01))) +
  labs(title    = "Per-individual sensitivity by method  (Global panel, n = 27)",
       subtitle = sub3,
       x = NULL, y = "Sensitivity") +
  base_theme +
  theme(legend.position  = "right",
        legend.title     = element_text(size = 9.5, face = "bold"),
        axis.text.x      = element_text(size = 10))

ggsave(file.path(outdir, "per_individual_sensitivity.pdf"), fig3, width = 8, height = 5.5)
ggsave(file.path(outdir, "per_individual_sensitivity.png"), fig3, width = 8, height = 5.5, dpi = 300)
cat("\nWrote fig3_per_individual_sensitivity\n")

# ── Figure 4: per-individual gain by ancestry group (Fix 2) ──────────────────
# Superpopulation label for each point, kept inside the panel using ggrepel
# if available, falling back to geom_text with nudge otherwise.
use_repel <- requireNamespace("ggrepel", quietly = TRUE)

gain_lo <- floor(min(dat$gain_african) * 200) / 200   # round to nearest 0.5 pp
gain_hi <- ceiling(max(dat$gain_african) * 200) / 200

# Group median bars
group_medians <- dat %>%
  group_by(ancestry_group) %>%
  summarise(med = median(gain_african), .groups = "drop")

# Significance annotation: place above the highest point in each group
sig_df <- tibble(
  ancestry_group = factor(c("Continental-African", "Non-African"),
                          levels = levels(dat$ancestry_group)),
  label = c(
    sprintf("n = %d\nmedian = %.1f pp", length(gain_cont_afr), median(gain_cont_afr)*100),
    sprintf("n = %d\nmedian = %.1f pp", length(gain_non_afr),  median(gain_non_afr)*100)
  ),
  y = gain_hi * 0.98
)
diasp_df <- tibble(
  ancestry_group = factor("African diaspora", levels = levels(dat$ancestry_group)),
  label = sprintf("n = %d\nmedian = %.1f pp", length(gain_diaspora), median(gain_diaspora)*100),
  y = gain_hi * 0.98
)
ann_df <- bind_rows(sig_df, diasp_df)

mw_subtitle <- sprintf(
  "Continental-African vs Non-African: p = %.3g (Mann-Whitney U, n = 5 vs 20)  |  not significant at α = 0.05",
  mw_cont_vs_non$p.value)

set.seed(42)
fig4_base <- ggplot(dat, aes(x = ancestry_group, y = gain_african, colour = ancestry_group)) +
  geom_crossbar(data = group_medians,
                aes(x = ancestry_group, y = med, ymin = med, ymax = med),
                inherit.aes = FALSE, colour = "grey15", width = 0.45, linewidth = 0.7) +
  geom_jitter(size = 3.2, width = 0.14) +
  geom_text(data = ann_df, aes(x = ancestry_group, y = y, label = label),
            inherit.aes = FALSE, size = 3, colour = "grey30", vjust = 1, lineheight = 0.9) +
  scale_colour_manual(values = group_cols, guide = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                     limits = c(gain_lo, gain_hi),
                     expand = expansion(mult = c(0.04, 0.08))) +
  labs(title    = "Per-individual sensitivity gain: African DB minus Default  (Global panel)",
       subtitle = mw_subtitle,
       x = NULL,
       y = "Sensitivity gain (African DB − Default)") +
  base_theme +
  theme(legend.position = "none", axis.text.x = element_text(size = 10.5))

# Add population labels (ggrepel preferred to avoid overlap)
if (use_repel) {
  fig4 <- fig4_base +
    ggrepel::geom_text_repel(aes(label = population_code), size = 2.8,
                             colour = "grey25", box.padding = 0.3,
                             min.segment.length = 0.2, seed = 42, max.overlaps = 20)
} else {
  fig4 <- fig4_base +
    geom_text(aes(label = population_code), size = 2.6, colour = "grey25",
              vjust = -0.8, hjust = 0.5, check_overlap = FALSE)
}

ggsave(file.path(outdir, "ancestry_gain.pdf"), fig4, width = 8, height = 5.5)
ggsave(file.path(outdir, "ancestry_gain.png"), fig4, width = 8, height = 5.5, dpi = 300)
cat("Wrote fig4_ancestry_gain\n")

cat("\nStats written to", statsdir, "\n")
