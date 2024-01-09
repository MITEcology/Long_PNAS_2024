# R code for analyze empirical ocean microbiota data
# reproduce figure 4D or figure Sxx

rm(list = ls())
source("code/toolbox.R")

# --- import time series data ---
data_donor <- read_delim("data/empirical_time_series/ocean_site_37.csv", delim = ",")

# --- coarse grain based on taxonomy ---
data_family <- data_donor %>%
  group_by(D0,D1,D2,D3,D4) %>%
  summarize(across(c(-1, -2), sum, .names = "s_{.col}"), .groups = "keep") %>%
  unite(eff_taxo, D0,D1,D2,D3,D4, sep = ";")

data_species <- data_donor %>%
  group_by(D0,D1,D2,D3,D4,D5,D6) %>%
  summarize(across(everything(), sum, .names = "s_{.col}"), .groups = "keep") %>% # make it standard
  unite(eff_taxo, D0,D1,D2,D3,D4,D5,D6, sep = ";")

# --- generate and analyze empirical data ---
data_plot <- data_family
data_plot <- data_plot %>%
  mutate(across(c(-1), ~ . / sum(.)))


sizes <- c(3:12) # choose the size range
steps <- 200
n_size <- max(sizes) - min(sizes) + 1
set.seed(180)
emp_ocean_tbl <- analyze_empirical_timeseries(sizes, steps, data_plot, can_thre = 1e-6, pa_thre = 1e-6, tele_alpha = 1e-1, rep_num = 1)

# --- plot: Figure 4D or Figure Sxx ---
set.seed(150)
p_ocean_size <- ggplot() +
  geom_boxplot(data = filter(emp_ocean_tbl, type == "randomized"), aes(size, ent_adap, color = type, group = cut_interval(size, 18)), size = 1) +
  geom_jitter(data = filter(emp_ocean_tbl, type == "randomized"), aes(size, ent_adap, color = type), size = 2, alpha = 0.45) +
  geom_boxplot(data = filter(emp_ocean_tbl, type == "empirical"), aes(size, ent_adap, color = type, group = cut_interval(size, 18)), size = 1, alpha = 0.4) +
  geom_jitter(data = filter(emp_ocean_tbl, type == "empirical"), aes(size, ent_adap, color = type), size = 2, alpha = 0.45) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = expression("Size of Taxon Pool" ~ (S)), y = expression("Switching Capacity" ~ (H))) +
  theme_my(
    axis.text.x = element_text(size = 25, vjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 35),
    axis.title.x = element_text(size = 50, vjust = 0),
    axis.title.y = element_text(size = 50),
    legend.position = "none",
    legend.text = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 30, face = "bold")
  ) +
  guides(color = guide_legend(title = "Type of time series")) +
  scale_color_manual(values = c("dodgerblue", "darkorange"), labels = c("Observed (Structured)", "Randomized (Unstructured)")) +
  scale_x_continuous(breaks = seq(3, 20, by = 1))
print(p_ocean_size)

ggsave("figs/fig4D.pdf", p_ocean_size, width = 12.5, height = 10, dpi = 300, device = pdf)