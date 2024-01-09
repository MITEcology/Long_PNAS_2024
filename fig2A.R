rm(list = ls())
source("code/toolbox.R")

# --- initialize ---
num <- 5 # set the size of the taxon pool
NUM <- 2^num # the number of compositions
steps <- 200 # set the number of data points

# --- generate feasible partitions ---
set.seed(200)
# romegass <- unif_omegass(num, steps)
romegass <- lnorm_omegass(num, steps, logmean = -0.8369882, logsd = 0.53636)
rcentrds <- random_sphere(num, steps)

# --- compute transition matrices and their switching capacities ---
gm_tbl <- tibble(
  beta = character(), type = character(), ent_adap = double(), size = integer()
) # a tibble to store H(β) relationships
for (step in c(seq_len(steps))) {
  r_center <- rcentrds[[step]]
  r_omegas <- romegass[[step]]
  # For unstructured transitions, directly compute ent_adap via unst_ent_adap()
  for (beta in c(0)) {
    gm_tbl <- rbind(gm_tbl,
      tibble(beta = sprintf("%.2f", beta), type = "Unstructured", ent_adap = unst_ent_adap(r_omegas), size = num)
    )
  }
  # For structured transitions, use the markov_analysis function to compute ent_adap
  for (beta in c(1:32) * 0.25) {
    grav_tmat <- gravity_tmat(r_omegas, r_center, beta)
    mkv <- markov_analysis(grav_tmat, r_omegas, TRUE)
    gm_tbl <- rbind(gm_tbl,
      tibble(beta = sprintf("%.2f", beta), type = "Structured", ent_adap = mkv$ent_adap, size = num)
    )
  }
  # print the computation progress, can be commented out
  print(step)
}
#rbind all gm_tbl with sizes range from c(3:12) will be gm_total used in Figure 2B,C

# --- generate plots for Figure 2A ---
# visualization of changing Structural parameter (β) with fixed pool size (S)
set.seed(150)
p_gm_beta <- ggplot() +
  geom_boxplot(
    data = gm_tbl %>% filter(type == "Unstructured"), aes(beta, ent_adap, color = type), size = 1) +
  geom_jitter(
    data = gm_tbl %>% filter(type == "Unstructured"), aes(beta, ent_adap, color = type), size = 1, alpha = 0.35) +
  geom_boxplot(
    data = gm_tbl %>% filter(type == "Structured"), aes(beta, ent_adap, color = type), size = 1) +
  geom_jitter(
    data = gm_tbl %>% filter(type == "Structured"), aes(beta, ent_adap, color = type), size = 1, alpha = 0.35) +
  geom_vline(xintercept = "2.50", linetype = "dashed", linewidth = 1, color = "grey30") + # beta_opt is specified manually
  annotate("text", x = 13, y = 0.06, label = "beta[opt]", size = 16.5, color = "grey30", parse = TRUE) +
  annotate("text", x = 26, y = 1, label = "Size of Taxon Pool", size = 13.5, color = "black", fontface = "bold") +
  annotate("text", x = 26, y = 0.93, label = "S==5", size = 13.5, color = "black", parse = TRUE) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = expression("Structural Parameter" ~ (beta)), y = expression("Switching Capacity" ~ (H))) +
  scale_color_manual(values = c("dodgerblue", "darkorange"), labels = c("Structured", "Unstructured")) +
  scale_x_discrete(breaks = sprintf("%.2f", seq(0, 8, by = 0.5))) +
  guides(NULL) +
  theme_my(
    axis.text.x = element_text(size = 25, angle = 45, hjust = 0.9),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 50, vjust = 0),
    axis.title.y = element_text(size = 50),
    legend.position = "none",
    legend.text = element_text(size = 25, face = "bold"),
    legend.title = element_text(size = 30),
  )
print(p_gm_beta)
ggsave("figs/sfig-lnorm2-a.pdf", width = 12.5, height = 10, units = "in", dpi = 300, device = pdf)
