library(igraph)
library(tidyverse)

source("scripts/validation/00_load_objects.R")

seed_kras <- result_kras$seeds
seed_nras <- result_nras$seeds
seed_hras <- result_hras$seeds

backbone_nodes <- Reduce(intersect, list(seed_kras, seed_nras, seed_hras))

clust <- result_kras[["clusters"]]
cluster_nodes <- names(clust[clust == 2])
cluster_nodes_noiso <- cluster_nodes[!cluster_nodes %in% c("P01116-2","P01116-1")] 

# cluster density
g_cluster <- induced_subgraph(g_kras, cluster_nodes)

cluster_density <- edge_density(g_cluster)
network_density <- edge_density(g_kras)

density_enrichment <- cluster_density / network_density

cat("Cluster density:", cluster_density, "\n")
cat("Network density:", network_density, "\n")
cat("Fold enrichment:", density_enrichment, "\n")

# cluster density - without KRAS isoforms
g_cluster_noiso <- induced_subgraph(g_kras, cluster_nodes_noiso)

cluster_density_noiso <- edge_density(g_cluster_noiso)

density_enrichment_noiso <- cluster_density_noiso / network_density

cat("Cluster density without KRAS isoforms:", cluster_density_noiso, "\n")
cat("Network density:", network_density, "\n")
cat("Fold enrichment without KRAS isoforms::", density_enrichment_noiso, "\n")



#betweenness centrality
bt <- betweenness(g_kras, normalized = TRUE)

cluster_bt <- bt[names(bt) %in% cluster_nodes]
backbone_bt <- bt[names(bt) %in% backbone_nodes]

wilcox_bt <- wilcox.test(cluster_bt, backbone_bt)

cat("Mean cluster betweenness:", mean(cluster_bt), "\n")
cat("Mean backbone betweenness:", mean(backbone_bt), "\n")
print(wilcox_bt)

#betweenness centrality - without KRAS isoforms
cluster_bt_noiso <- bt[names(bt) %in% cluster_nodes_noiso]

wilcox_bt_noiso <- wilcox.test(cluster_bt_noiso, backbone_bt)

cat("Mean cluster betweenness without KRAS isoforms:", mean(cluster_bt), "\n")
cat("Mean backbone betweenness:", mean(backbone_bt), "\n")
print(wilcox_bt_noiso)



df_bt <- data.frame(
  betweenness = c(cluster_bt, cluster_bt_noiso, backbone_bt),
  group = c(
    rep("KRAS module", length(cluster_bt)),
    rep("KRAS module (no isoforms)", length(cluster_bt_noiso)),
    rep("Conserved backbone", length(backbone_bt))
  )
)

df_bt$bet_log <- log10(df_bt$betweenness + 1e-6)

ggplot(df_bt, aes(x = group, y = bet_log , fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = NA)+
  #geom_boxplot(width = 0.1) +
  geom_jitter(width = 0.12, size = 1.5, alpha = 0.4)+
  scale_fill_manual(values = c(
    "KRAS module" = "#5B2A86",   # viola
    "KRAS module (no isoforms)" = "#9C9EDE", 
    "Conserved backbone" = "#2CA58D"        # verde
  )) +
  labs(
    x = "",
    y = "Betweenness centrality (log10 scale)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  ) +
  coord_cartesian(ylim = c(-7.5, 0.5))+
  ggpubr::stat_compare_means(
    comparisons = list(
      c("Conserved backbone", "KRAS module"),
      c("Conserved backbone", "KRAS module (no isoforms)")
    ),
    label = "p.format",
    label.y = c(-0.6, -0.2),
    step.increase = 0
  )

results <- tibble(
  cluster_density = cluster_density,
  cluster_density_noiso,
  network_density = network_density,
  density_enrichment = density_enrichment,
  density_enrichment_noiso, 
  mean_cluster_bt = mean(cluster_bt),
  mean_backbone_bt = mean(backbone_bt),
  betweenness_p = wilcox_bt$p.value,
  mean_cluster_bt_noiso = mean(cluster_bt_noiso),
  betweenness_p_noiso = wilcox_bt_noiso$p.value
)

write.csv(results, "output/validation_results/module_validation_kras.csv", row.names = FALSE)
