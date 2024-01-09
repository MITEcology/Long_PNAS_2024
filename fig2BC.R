rm(list = ls())
source("code/toolbox.R")

# load("data/theoretical_data/uniform_omega_tbl.RData")
load("data/theoretical_data/lnorm1_omega_tbl.RData")
# --- generate plots for figure 2B and 2C ---
size_range <- unique(gm_total$size)
size_cut <- max(size_range) - min(size_range) + 1

# --- Figure 2B: optimal structural parameter vs. size of taxon pool ---

beta_opt_size_tbl <- enframe(beta_opt_size, name = "size", value = "beta_opt") %>%
  mutate(size = as.integer(size), beta_opt = as.numeric(beta_opt))
p_beta_opt <- ggplot(beta_opt_size_tbl) +
  geom_abline(slope = 0.5, intercept = -0.5, linetype = "dashed") +
  geom_point(aes(size, beta_opt), color = "grey30", size = 6) +
  scale_x_continuous(breaks = c(3:12)) +
  scale_y_continuous(breaks = seq(0, 8, by = 1), limits = c(0, 8)) +
  labs(x = expression("Size of Taxon Pool" ~ (S)), y = expression("Optimal Struct. Param." ~ (beta["opt"]))) +
  annotate("text", x = 3.8, y = 3, label = "beta[opt]==frac(S,2)-frac(1,2)", size = 13.5, color = "black", parse = TRUE) +
  coord_cartesian(xlim = c(3, 12), ylim = c(0, 8)) +
  theme_my(
    axis.text.x = element_text(size = 25, vjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 50, vjust = 0),
    axis.title.y = element_text(size = 40)
  )
print(p_beta_opt)
ggsave("figs/sfig-lnorm1-b.pdf", width = 12.5, height = 10, units = "in", dpi = 300, device = pdf)


# --- Figure 2C: switching capacity vs. size of taxon pool ---
gm_total <- gm_total %>%
  mutate(type = if_else(beta == beta_opt_size[as.character(size)], "Optimized", type)) #mark the beta_optimal boxplots
p_gm_size <- ggplot() +
  geom_boxplot(
    data = gm_total %>% filter(beta == "0.00"), aes(size, ent_adap, color = type, group = cut_interval(size, size_cut)), size = 0.8
  ) +
  geom_jitter(
    data = gm_total %>% filter(beta == "0.00"), aes(size, ent_adap, color = type), size = 1, alpha = 0.35) +
  geom_boxplot(
    data = gm_total %>% filter(beta == "2.00"), aes(size, ent_adap, color = type, group = cut_interval(size, size_cut)), size = 0.8
  ) +
  geom_jitter(
    data = gm_total %>% filter(beta == "2.00"), aes(size, ent_adap, color = type), size = 1, alpha = 0.35) +
  geom_boxplot(
    data = gm_total %>% filter(type == "Optimized"), aes(size, ent_adap, color = type, group = cut_interval(size, size_cut)), size = 0.8
  ) +
  geom_jitter(
    data = gm_total %>% filter(type == "Optimized"), aes(size, ent_adap, color = type), size = 1, alpha = 0.35) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = expression("Size of Taxon Pool" ~ (S)), y = expression("Switching Capacity" ~ (H))) +
  scale_color_manual(
    values = c("grey50", "dodgerblue", "darkorange"),
    labels = c(expression(beta ~ "=" ~ beta["opt"] ~ (S)), expression(beta ~ "=" ~ 2), expression(beta ~ "=" ~ 0))
  ) +
  scale_x_continuous(breaks = seq(min(size_range), max(size_range), by = 1)) +
  guides(color = guide_legend(title = "Structural\nParameters")) +
  theme_my(
    axis.text.x = element_text(size = 25, vjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 50, vjust = 0),
    axis.title.y = element_text(size = 50),
    # legend.position = c(0.78, 0.69),
    legend.position = "none",
    legend.text = element_text(size = 35, hjust = 0),
    legend.key.height = unit(2, "lines"),
    legend.title = element_text(size = 35, face = "bold")
  )
print(p_gm_size)
ggsave("figs/sfig-lnorm1-c.pdf", width = 12.5, height = 10, units = "in", dpi = 300, device = pdf)
