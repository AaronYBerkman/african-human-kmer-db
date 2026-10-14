library(dplyr)
library(readr)

# Supplementary Tables S5 and S7, formatted for LaTeX, from the full statistics
# written by 28_supplementary_table_S1.R.
#   S5 (stats_display.csv): sensitivity and specificity with Wilson 95% CIs
#   S7 (gap_display.csv):   sensitivity gain of the African database over the
#                           default, in percentage points, with a Newcombe 95% CI
#                           (same method as the population comparison in script
#                           28) and the McNemar p-value

repo_root <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("--file=", "", args[grep("^--file=", args)])
  dirname(dirname(dirname(normalizePath(file_arg))))
})
statsdir <- file.path(repo_root, "stats")

s1 <- read_csv(file.path(statsdir, "supplementary_table_S1_full_stats.csv"), show_col_types = FALSE)

pct  <- function(x, digits) formatC(100 * x, format = "f", digits = digits)
pfmt <- function(p) ifelse(p < 0.001, "$<$0.001", formatC(p, format = "f", digits = 3))
methods <- c("Default (human RefSeq)", "African DB only", "Combined (Stage 1+2)")

# ── Table S5 ─────────────────────────────────────────────────────────────────
s5_order <- tibble(background = c("Microbial", "Viral", "Microbial", "Viral"),
                   human_pop  = c("Continental-African", "Continental-African", "Global", "Global"),
                   block = 1:4)
s5 <- s1 %>%
  inner_join(s5_order, by = c("background", "human_pop")) %>%
  mutate(method_rank = match(method, methods)) %>%
  arrange(block, human_pct, method_rank) %>%
  transmute(background, human_pop, human_pct, method,
            sensitivity_ci = sprintf("%s\\%% (%s--%s)", pct(sensitivity, 2), pct(sens_ci_lower, 2), pct(sens_ci_upper, 2)),
            specificity_ci = sprintf("%s\\%% (%s--%s)", pct(specificity, 3), pct(spec_ci_lower, 3), pct(spec_ci_upper, 3)))

# ── Table S7 ─────────────────────────────────────────────────────────────────
# Newcombe hybrid score interval from the two Wilson intervals (Newcombe 1998, method 10)
s7_order <- tibble(background = c("Viral", "Microbial", "Microbial", "Viral"),
                   human_pop  = c("Global", "Global", "Continental-African", "Continental-African"),
                   block = 1:4)
wide <- s1 %>%
  filter(method %in% methods[1:2]) %>%
  mutate(m = ifelse(method == methods[1], "def", "afr")) %>%
  select(background, human_pop, human_pct, m, sensitivity, sens_ci_lower, sens_ci_upper, mcnemar_p_value) %>%
  tidyr::pivot_wider(names_from = m, values_from = c(sensitivity, sens_ci_lower, sens_ci_upper, mcnemar_p_value))
s7 <- wide %>%
  inner_join(s7_order, by = c("background", "human_pop")) %>%
  arrange(block, human_pct) %>%
  mutate(gap   = sensitivity_afr - sensitivity_def,
         lower = gap - sqrt((sensitivity_afr - sens_ci_lower_afr)^2 + (sens_ci_upper_def - sensitivity_def)^2),
         upper = gap + sqrt((sens_ci_upper_afr - sensitivity_afr)^2 + (sensitivity_def - sens_ci_lower_def)^2)) %>%
  transmute(background, human_pop, human_pct,
            gap_ci = sprintf("%s (%s--%s)", pct(gap, 2), pct(lower, 2), pct(upper, 2)),
            mcnemar_p = pfmt(mcnemar_p_value_afr))

write_csv(s5, file.path(statsdir, "stats_display.csv"))
write_csv(s7, file.path(statsdir, "gap_display.csv"))
cat("Wrote stats_display.csv (", nrow(s5), "rows) and gap_display.csv (", nrow(s7), "rows)\n")
