library(dplyr)
library(readr)
library(tidyr)

repo_root <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("--file=", "", args[grep("^--file=", args)])
  dirname(dirname(dirname(normalizePath(file_arg))))
})

wilson_ci <- function(x, n, conf_level = 0.95) {
  z <- qnorm(1 - (1 - conf_level) / 2)
  p_hat <- x / n
  denom <- 1 + z^2 / n
  center <- (p_hat + z^2 / (2 * n)) / denom
  half_width <- (z / denom) * sqrt((p_hat * (1 - p_hat) / n) + (z^2 / (4 * n^2)))
  data.frame(estimate = p_hat, lower = pmax(0, center - half_width), upper = pmin(1, center + half_width))
}

viral_path     <- file.path(repo_root, "benchmark_results/viral_mixed_benchmark_results/viral_mixed_benchmark_summary.csv")
bacterial_path <- file.path(repo_root, "benchmark_results/bacterial_benchmark_results/bacterial_benchmark_summary.csv")
african25_path <- file.path(repo_root, "benchmark_results/african25_benchmark_results/african25_benchmark_summary.csv")
outdir <- file.path(repo_root, "stats")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

viral_raw <- read_csv(viral_path, show_col_types = FALSE)
bacterial_raw <- read_csv(bacterial_path, show_col_types = FALSE)
african25_raw <- read_csv(african25_path, show_col_types = FALSE)

viral_pooled <- viral_raw %>%
  transmute(
    human_pct = human_titre_pct, true_human, true_other,
    default_tp, default_fp, african_tp, african_fp, combined_tp, combined_fp,
    background = "Viral", human_pop = "Global"
  )

bacterial_pooled <- bacterial_raw %>%
  transmute(
    human_pct = human_titre_pct, true_human, true_other = true_bacterial,
    default_tp, default_fp, african_tp, african_fp, combined_tp, combined_fp,
    background = "Microbial", human_pop = "Global"
  )

african25_pooled <- african25_raw %>%
  transmute(
    human_pct = human_titre_pct, true_human, true_other,
    default_tp, default_fp, african_tp, african_fp, combined_tp, combined_fp,
    background = ifelse(arm == "microbial", "Microbial", "Viral"),
    human_pop = "Continental-African"
  )

combined <- bind_rows(viral_pooled, bacterial_pooled, african25_pooled)

# Build one row per (background, titre, method) with sensitivity + specificity + both CIs
build_method_table <- function(df, method_name, tp_col, fp_col) {
  df %>%
    rowwise() %>%
    mutate(
      sens_ci = list(wilson_ci(.data[[tp_col]], true_human)),
      fp_ci   = list(wilson_ci(.data[[fp_col]], true_other))
    ) %>%
    ungroup() %>%
    mutate(method = method_name) %>%
    unnest(sens_ci, names_sep = "_") %>%
    unnest(fp_ci, names_sep = "_") %>%
    transmute(
      background, human_pop, human_pct, method,
      true_human, true_other,
      tp = .data[[tp_col]], fp = .data[[fp_col]],
      sensitivity = sens_ci_estimate, sens_ci_lower = sens_ci_lower, sens_ci_upper = sens_ci_upper,
      specificity = 1 - fp_ci_estimate, spec_ci_lower = 1 - fp_ci_upper, spec_ci_upper = 1 - fp_ci_lower
    )
}

full_table <- bind_rows(
  build_method_table(combined, "Default (human RefSeq)", "default_tp", "default_fp"),
  build_method_table(combined, "African DB only", "african_tp", "african_fp"),
  build_method_table(combined, "Combined (Stage 1+2)", "combined_tp", "combined_fp")
)

# McNemar's test: Default vs African DB (the comparison reported in the main text)
mcnemar_table <- combined %>%
  rowwise() %>%
  mutate(
    b = max(0, african_tp - default_tp),
    c = max(0, default_tp - african_tp),
    mcnemar_stat = ifelse(b + c > 0, (abs(b - c) - 1)^2 / (b + c), NA_real_),
    mcnemar_p = ifelse(b + c > 0, pchisq(mcnemar_stat, df = 1, lower.tail = FALSE), NA_real_)
  ) %>%
  ungroup() %>%
  transmute(background, human_pop, human_pct,
            discordant_pairs_b_african_only = b,
            discordant_pairs_c_default_only = c,
            mcnemar_chi_sq = mcnemar_stat,
            mcnemar_p_value = mcnemar_p)

# Merge: one comprehensive row per background x human_pop x titre x method, with
# McNemar result attached only where relevant (Default and African DB rows both
# show the same McNemar p-value for that cell, since it's a test of those two
# methods specifically)
final_table <- full_table %>%
  left_join(mcnemar_table, by = c("background", "human_pop", "human_pct")) %>%
  arrange(human_pop, background, human_pct, factor(method, levels = c("Default (human RefSeq)", "African DB only", "Combined (Stage 1+2)")))

write_csv(final_table, file.path(outdir, "supplementary_table_S1_full_stats.csv"))

cat("=== Supplementary Table S1: full statistics ===\n")
print(final_table %>%
        mutate(across(where(is.numeric), ~round(.x, 6))) %>%
        select(human_pop, background, human_pct, method, tp, fp, sensitivity, sens_ci_lower, sens_ci_upper,
               specificity, spec_ci_lower, spec_ci_upper, mcnemar_p_value),
      n = 60)

cat("\nWrote full table to", file.path(outdir, "supplementary_table_S1_full_stats.csv"), "\n")

# ============================================================
# Cross-population comparison: is sensitivity/specificity for a given
# method different between Global and Continental-African
# ancestry? The original benchmarking design never had more than one
# human population, so there's no established paired test for this -
# unlike Default-vs-African DB (same reads, same individuals, hence
# McNemar's test on paired counts), the two population panels are built
# from different reads/individuals entirely - an independent-samples
# comparison, not a paired one. Newcombe's method (Newcombe 1998,
# "Method 10") extends the same Wilson-score framework already used for
# each individual proportion above to a CI on the difference of two
# independent proportions, so this reuses wilson_ci() rather than
# introducing an unrelated statistical framework.
# ============================================================
newcombe_diff_ci <- function(x1, n1, x2, n2, conf_level = 0.95) {
  ci1 <- wilson_ci(x1, n1, conf_level)
  ci2 <- wilson_ci(x2, n2, conf_level)
  diff <- ci1$estimate - ci2$estimate
  lower <- diff - sqrt((ci1$estimate - ci1$lower)^2 + (ci2$upper - ci2$estimate)^2)
  upper <- diff + sqrt((ci1$upper - ci1$estimate)^2 + (ci2$estimate - ci2$lower)^2)
  data.frame(diff = diff, lower = lower, upper = upper)
}

population_wide <- full_table %>%
  filter(human_pop %in% c("Global", "Continental-African")) %>%
  mutate(pop_key = ifelse(human_pop == "Global", "mixed", "african")) %>%
  select(background, pop_key, human_pct, method, tp, fp, true_human, true_other) %>%
  pivot_wider(names_from = pop_key, values_from = c(tp, fp, true_human, true_other))

population_comparison <- population_wide %>%
  rowwise() %>%
  mutate(
    sens_diff_ci = list(newcombe_diff_ci(tp_mixed, true_human_mixed, tp_african, true_human_african)),
    spec_diff_ci = list(newcombe_diff_ci(true_other_mixed - fp_mixed, true_other_mixed,
                                          true_other_african - fp_african, true_other_african))
  ) %>%
  ungroup() %>%
  unnest(sens_diff_ci, names_sep = "_") %>%
  unnest(spec_diff_ci, names_sep = "_") %>%
  transmute(
    background, human_pct, method,
    sensitivity_mixed = tp_mixed / true_human_mixed,
    sensitivity_continental_african = tp_african / true_human_african,
    sensitivity_diff = sens_diff_ci_diff,
    sensitivity_diff_ci_lower = sens_diff_ci_lower,
    sensitivity_diff_ci_upper = sens_diff_ci_upper,
    specificity_mixed = 1 - fp_mixed / true_other_mixed,
    specificity_continental_african = 1 - fp_african / true_other_african,
    specificity_diff = spec_diff_ci_diff,
    specificity_diff_ci_lower = spec_diff_ci_lower,
    specificity_diff_ci_upper = spec_diff_ci_upper
  ) %>%
  arrange(background, human_pct, factor(method, levels = c("Default (human RefSeq)", "African DB only", "Combined (Stage 1+2)")))

write_csv(population_comparison, file.path(outdir, "population_comparison_wilson_ci.csv"))

cat("\n=== Cross-population comparison: Global vs Continental-African ===\n")
cat("(diff = Global minus Continental-African; CI excludes 0 => significant at the stated level)\n")
print(population_comparison %>% mutate(across(where(is.numeric), ~round(.x, 6))), n = 60)

cat("\nWrote population comparison to", file.path(outdir, "population_comparison_wilson_ci.csv"), "\n")
