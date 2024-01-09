# R code for priciple coordinate analysis of ocean microbiota data
# reproduce figure S1D

rm(list = ls())
source("code/toolbox.R")

# --- import time series data and binarize ---
data_donor <- read_delim("data/empirical_time_series/ocean_site_37.csv", delim = ",")
data_family <- data_donor %>%
  group_by(D0,D1,D2,D3,D4) %>%
  summarize(across(c(-1, -2), sum, .names = "s_{.col}"), .groups = "keep") %>%
  unite(eff_taxo, D0,D1,D2,D3,D4, sep = ";")
data_focus <- bin_pa(data_family)
data_focus <- bin_pa_p(data_family)

# --- PCoA and Clustering ---
distance_matrix <- dist(data_focus[, -1])
pcoa <- cmdscale(distance_matrix)

k <- 4 # number of clusters
kmeans_result <- kmeans(pcoa, centers = k)
cluster_assignments <- kmeans_result$cluster
clustered_data <- data.frame(pcoa, Cluster = as.factor(cluster_assignments))

# --- plot the data with clusters ---
p_pcoa <- ggplot(clustered_data, aes(x = X1, y = X2, color = Cluster)) +
  geom_point(size = 8, alpha = 0.5) +
  labs(x = "PCo1", y = "PCo2", title = NULL) +
  theme_my(
    axis.text.x = element_text(size = 35),
    axis.text.y = element_text(size = 35),
    axis.title.x = element_text(size = 45, vjust = 0),
    axis.title.y = element_text(size = 45),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 30)
  )
print(p_pcoa)