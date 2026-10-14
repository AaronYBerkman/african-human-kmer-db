library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

# Repo root, resolved from this script's own path (scripts/analysis/<this file>)
# so paths below work regardless of the working directory Rscript is invoked from.
repo_root <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("--file=", "", args[grep("^--file=", args)])
  dirname(dirname(dirname(normalizePath(file_arg))))
})

viral_path     <- file.path(repo_root, "benchmark_results/viral_mixed_benchmark_results/viral_mixed_benchmark_summary.csv")
microbial_path <- file.path(repo_root, "benchmark_results/bacterial_benchmark_results/bacterial_benchmark_summary.csv")
african25_path <- file.path(repo_root, "benchmark_results/african25_benchmark_results/african25_benchmark_summary.csv")

outdir <- file.path(repo_root, "figures")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

viral_raw    <- read_csv(viral_path,    show_col_types = FALSE)
microbial_raw <- read_csv(microbial_path, show_col_types = FALSE)
african25_raw <- read_csv(african25_path, show_col_types = FALSE)

cell1_viral <- viral_raw %>%
  transmute(human_pct = human_titre_pct, true_human, true_other,
            default_tp, african_tp, combined_tp,
            default_fp, african_fp, combined_fp,
            background = "Viral", human_pop = "Global")

cell2_microbial <- microbial_raw %>%
  transmute(human_pct = human_titre_pct, true_human, true_other = true_bacterial,
            default_tp, african_tp, combined_tp,
            default_fp, african_fp, combined_fp,
            background = "Microbial", human_pop = "Global")

cell34_african <- african25_raw %>%
  transmute(human_pct = human_titre_pct, true_human, true_other,
            default_tp, african_tp, combined_tp,
            default_fp, african_fp, combined_fp,
            background = ifelse(arm == "microbial", "Microbial", "Viral"),
            human_pop = "Continental-African")

raw <- bind_rows(cell1_viral, cell2_microbial, cell34_african) %>%
  mutate(
    background = factor(background, levels = c("Viral", "Microbial")),
    human_pop  = factor(human_pop,  levels = c("Global", "Continental-African"))
  )

method_cols   <- c("Default (Human RefSeq)" = "#E64B35",
                   "African DB only"         = "#4DBBD5",
                   "Combined (Stage 1+2)"    = "#00A087")
method_shapes <- c("Default (Human RefSeq)" = 16, "African DB only" = 17, "Combined (Stage 1+2)" = 15)
method_ltype  <- c("Default (Human RefSeq)" = "solid", "African DB only" = "dashed", "Combined (Stage 1+2)" = "dotdash")

base_theme <- theme_classic(base_size = 10.5) +
  theme(
    strip.text         = element_text(face = "bold", size = 9.5),
    strip.background   = element_rect(fill = "grey94", colour = NA),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
    panel.spacing.y    = unit(1.5, "lines"),
    plot.margin        = margin(4, 30, 2, 4),
    legend.position    = "bottom",
    legend.text        = element_text(size = 10),
    plot.title         = element_text(face = "bold", size = 11.5, hjust = 0)
  )

x_scale <- scale_x_continuous(
  breaks = seq(10, 90, by = 20),
  labels = paste0(seq(10, 90, by = 20), "%"),
  expand = expansion(mult = c(0.05, 0.15))
)

mcnemar_p_paired <- function(x1, x2) {
  b <- pmax(0, x1 - x2); c <- pmax(0, x2 - x1)
  stat <- ifelse(b + c > 0, (abs(b - c) - 1)^2 / (b + c), NA_real_)
  ifelse(b + c > 0, pchisq(stat, df = 1, lower.tail = FALSE), NA_real_)
}

# ── Figure 1: Absolute sensitivity ──────────────────────────────────────────

sens_long <- raw %>%
  pivot_longer(cols = c(default_tp, african_tp, combined_tp),
               names_to = "method_key", values_to = "tp") %>%
  mutate(
    method = case_when(
      method_key == "default_tp"  ~ "Default (Human RefSeq)",
      method_key == "african_tp"  ~ "African DB only",
      method_key == "combined_tp" ~ "Combined (Stage 1+2)"
    ),
    sensitivity = tp / true_human,
    method = factor(method, levels = names(method_cols))
  )

mean_sens <- sens_long %>%
  group_by(background, human_pop, method) %>%
  summarise(mean_sens = mean(sensitivity), .groups = "drop") %>%
  mutate(human_pct = 90, label = scales::percent(mean_sens, accuracy = 0.1))

sig_sens <- raw %>%
  mutate(p_african  = mcnemar_p_paired(african_tp,  default_tp),
         p_combined = mcnemar_p_paired(combined_tp, default_tp),
         max_p = pmax(p_african, p_combined, na.rm = TRUE)) %>%
  group_by(human_pop) %>%
  summarise(max_p = max(max_p), .groups = "drop") %>%
  mutate(
    human_pop = factor(human_pop, levels = levels(raw$human_pop)),
    human_pct = ifelse(human_pop == "Global", 12, 88),
    label     = ifelse(max_p < 0.001, "p < 0.001", paste0("max p = ", signif(max_p, 2))),
    hjust_val = ifelse(human_pop == "Global", 0, 1)
  )

fig1 <- ggplot(sens_long, aes(x = human_pct, y = sensitivity,
                               colour = method, linetype = method, shape = method, group = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4) +
  geom_text(data = mean_sens, aes(x = human_pct, y = mean_sens, label = label),
            hjust = -0.2, size = 3, fontface = "bold", show.legend = FALSE) +
  geom_text(data = sig_sens,
            mapping = aes(x = human_pct, y = -Inf, label = label, hjust = hjust_val),
            inherit.aes = FALSE, size = 3.2, fontface = "italic",
            colour = "grey30", vjust = -0.5) +
  facet_grid(background ~ human_pop, scales = "free_y") +
  scale_colour_manual(values = method_cols,   name = NULL) +
  scale_linetype_manual(values = method_ltype, name = NULL) +
  scale_shape_manual(values = method_shapes,   name = NULL) +
  x_scale +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(title = "Absolute sensitivity by method, background type, and human population",
       x = "Human read content", y = "Sensitivity") +
  base_theme

ggsave(file.path(outdir, "sensitivity.pdf"), fig1, width = 9.5, height = 6.8)
ggsave(file.path(outdir, "sensitivity.png"), fig1, width = 9.5, height = 6.8, dpi = 300)
cat("Wrote sensitivity\n")

# ── Figure 2: Absolute specificity ──────────────────────────────────────────

spec_long <- raw %>%
  pivot_longer(cols = c(default_fp, african_fp, combined_fp),
               names_to = "method_key", values_to = "fp") %>%
  mutate(
    method = case_when(
      method_key == "default_fp"  ~ "Default (Human RefSeq)",
      method_key == "african_fp"  ~ "African DB only",
      method_key == "combined_fp" ~ "Combined (Stage 1+2)"
    ),
    specificity = 1 - fp / true_other,
    method = factor(method, levels = names(method_cols))
  )

fig2 <- ggplot(spec_long, aes(x = human_pct, y = specificity,
                               colour = method, linetype = method, shape = method, group = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4) +
  facet_grid(background ~ human_pop) +
  scale_colour_manual(values = method_cols,   name = NULL) +
  scale_linetype_manual(values = method_ltype, name = NULL) +
  scale_shape_manual(values = method_shapes,   name = NULL) +
  x_scale +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.01),
                     limits = c(0.999, 1.0001)) +
  labs(title = "Absolute specificity by method, background type, and human population",
       x = "Human read content", y = "Specificity") +
  base_theme

ggsave(file.path(outdir, "specificity.pdf"), fig2, width = 9.5, height = 6.8)
ggsave(file.path(outdir, "specificity.png"), fig2, width = 9.5, height = 6.8, dpi = 300)
cat("Wrote specificity\n")
