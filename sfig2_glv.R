rm(list = ls())
source("code/toolbox.R")

# --- generate partitions from generalized Lotka-Volterra model ---
# set the ensemble parameters for random interaction matrces
num <- 4
NUM <- 2 ^ num
stren <- 1
conne <- 1
steps <- 1000
record1 <- get_compo(num, -1)
record0 <- get_compo(num, 0)

set.seed(200)
interaction_mats <- ensem_inmat_rand("norm", num, steps, stren, conne)
load("data/theoretical_data/glv_omegass_seed200_uniform.RData") # can be computed from `omegass <- ensem_geometry(interaction_mats)`

# --- compute transition matrices and their switching capacities ---
gm_glv <- tibble(
  beta = character(), type = character(), ent_adap = double()
)
for (beta in c(0)) {
  for (step in c(1:200)) {
    omegas <- omegass[[step]]
    gm_glv <- rbind(gm_glv,
      tibble(beta = sprintf("%.2f", beta), type = "Unstructured", ent_adap = unst_ent_adap(omegas))
    )
  }
}
for (beta in c(1:32) * 0.25) {
  grav_tmats <- ensem_tmat_LV(interaction_mats[1:200], omegass[1:200], beta)
  results <- ensem_trans_analysis(grav_tmats, omegass[1:200], TRUE)
  gm_glv <- rbind(gm_glv,
    tibble(beta = sprintf("%.2f", beta), type = "Structured", ent_adap = results[, 2])
  )
}

# --- generate plots for Figure S1B ---
set.seed(150)
p_gm_glv <- ggplot() +
  geom_boxplot(
    data = gm_glv %>% filter(type == "Unstructured"), aes(beta, ent_adap, color = type), size = 1) +
  geom_jitter(
    data = gm_glv %>% filter(type == "Unstructured"), aes(beta, ent_adap, color = type), size = 1, alpha = 0.35) +
  geom_boxplot(
    data = gm_glv %>% filter(type == "Structured"), aes(beta, ent_adap, color = type), size = 1) +
  geom_jitter(
    data = gm_glv %>% filter(type == "Structured"), aes(beta, ent_adap, color = type), size = 1, alpha = 0.35) +
  # annotate("text", x = 13, y = 0.06, label = "beta[opt]", size = 16.5, color = "grey30", parse = TRUE) +
  annotate("text", x = 26, y = 1, label = "Size of Taxon Pool", size = 13.5, color = "black", fontface = "bold") +
  annotate("text", x = 26, y = 0.93, label = "S==4", size = 13.5, color = "black", parse = TRUE) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = expression("Structural Parameter" ~ (beta)), y = expression("Switching Capacity" ~ (H))) +
  scale_color_manual(values = c("dodgerblue", "darkorange"), labels = c("Structured", "Unstructured")) +
  scale_x_discrete(breaks = sprintf("%.2f", seq(0, 8, by = 0.5))) +
  guides(color = guide_legend(title = "Type of Transitions")) +
  theme_my(
  axis.text.x = element_text(size = 25, angle = 45, hjust = 0.9),
  axis.text.y = element_text(size = 30),
  axis.title.x = element_text(size = 50, vjust = 0),
  axis.title.y = element_text(size = 50),
  legend.position = "none",
  legend.text = element_text(size = 25, face = "bold"),
  legend.title = element_text(size = 30)
  )
print(p_gm_glv)
ggsave("figs/figS2-glv.pdf", width = 12.5, height = 10, units = "in", dpi = 300, device = pdf)
