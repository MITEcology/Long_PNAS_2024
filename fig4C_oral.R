# R code for analyze empirical gut microbiota data
# reproduce figure 4C or figure Sxx

rm(list = ls())
source("code/toolbox.R")

# --- import time series data ---
data_donor <- read_delim("data/empirical_time_series/oral_subject_a.csv", delim = ",")

# --- coarse grain based on taxonomy ---
data_family <- data_donor %>%
  group_by(D0,D1,D2,D3,D4) %>%
  summarize(across(c(-1, -2), sum, .names = "s_{.col}"), .groups = "keep") %>%
  unite(eff_taxo, D0,D1,D2,D3,D4, sep = ":")
data_genus <- data_donor %>%
  group_by(D0,D1,D2,D3,D4,D5) %>%
  summarize(across(c(-1), sum, .names = "s_{.col}"), .groups = "keep") %>%
  unite(eff_taxo, D0,D1,D2,D3,D4,D5, sep = ":")
data_species <- data_donor %>%
  group_by(D0,D1,D2,D3,D4,D5,D6) %>%
  summarize(across(everything(), sum, .names = "s_{.col}"), .groups = "keep") %>%
  unite(eff_taxo, D0,D1,D2,D3,D4,D5,D6, sep = ";")


# --- generate and analyze empirical data ---
data_plot <- data_family
data_plot <- data_plot %>%
  mutate(across(c(-1), ~ . / sum(.)))

sizes <- c(3:12) # choose the size range
steps <- 200
n_size <- max(sizes) - min(sizes) + 1
set.seed(180)
emp_oral_tbl <- analyze_empirical_timeseries(sizes, steps, data_plot, can_thre = 1e-5, pa_thre = 1e-6, tele_alpha = 1e-1, rep_num = 1)


# --- plot: Figure 4C or Figure Sxx ---
set.seed(150)
p_oral_size <- ggplot() +
  geom_boxplot(data = filter(emp_oral_tbl, type == "randomized"), aes(size, ent_adap, color = type, group = cut_interval(size, 18)), size = 1) +
  geom_jitter(data = filter(emp_oral_tbl, type == "randomized"), aes(size, ent_adap, color = type), size = 2, alpha = 0.45) +
  geom_boxplot(data = filter(emp_oral_tbl, type == "empirical"), aes(size, ent_adap, color = type, group = cut_interval(size, 18)), size = 1, alpha = 0.4) +
  geom_jitter(data = filter(emp_oral_tbl, type == "empirical"), aes(size, ent_adap, color = type), size = 2, alpha = 0.45) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = expression("Size of Taxon Pool" ~ (S)), y = expression("Switching Capacity" ~ (H))) +
  # labs(x = NULL, y = NULL) +
  theme_my(
    axis.text.x = element_text(size = 25, vjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 35),
    axis.title.x = element_text(size = 50, vjust = 0),
    axis.title.y = element_text(size = 50),
    legend.position = "none"
  ) +
  scale_color_manual(values = c("dodgerblue", "darkorange"), labels = c("Empirical", "Randomized")) +
  scale_x_continuous(breaks = seq(3, 12, by = 1))
print(p_oral_size)

ggsave("figs/fig4C.pdf", p_oral_size, width = 12.5, height = 10, dpi = 100, device = pdf)