library(igraph)
library(tidyverse)

source("scripts/validation/00_load_objects.R")

seed_kras <- result_kras$seeds
seed_nras <- result_nras$seeds
seed_hras <- result_hras$seeds

backbone_nodes <- Reduce(intersect, list(seed_kras, seed_nras, seed_hras))

clust <- result_kras[["clusters"]]
cluster_nodes <- names(clust[clust == 2])

bridges_nodes <- explore_kras$bridge_direct$node

library(tidygraph)
library(ggraph)
set.seed(123)

g_tbl <- as_tbl_graph(g_kras)

g_tbl <- g_tbl %>%
  mutate(type = case_when(
    name %in% "P01116"~ "KRAS",
    name %in% c("P01116-2","P01116-1") ~ "KRAS isoforms",
    name %in% cluster_nodes ~ "KRAS module",
    name %in% bridges_nodes ~ "Bridge nodes",
    TRUE ~ "Backbone"
  ))


ggraph(g_tbl, layout = "fr") +
  
  geom_edge_link(
    color = "grey85",
    width = 0.2,
    alpha = 0.4
  ) +
  
  ggforce::geom_mark_hull(
    data = function(x) x[x$type == "KRAS module",],
    aes(x = x, y = y),
    fill = "#5B2A86",
    alpha = 0.1,
    color = "#5B2A86",
    expand = unit(3, "mm")
  ) +
  
  geom_node_point(
    aes(color = type, size = type),
    alpha = 0.7
  ) +
  
  scale_color_manual(values = c(
    "Backbone" = "grey75",
    "KRAS" = "black",
    "KRAS module" = "#5B2A86",
    "KRAS isoforms" = "#E64B35",
    "Bridge nodes" = "#2CA58D"
  )) +
  
  scale_size_manual(values = c(
    "Backbone" = 1.2,
    "KRAS module" = 2.5,
    "KRAS isoforms" = 3.5,
    "Bridge nodes" = 1.8,
    "KRAS" = 5
  )) +
  
  geom_node_text(
    aes(label = ifelse(type %in% c("KRAS isoforms", "KRAS"), name, "")),
    repel = TRUE,
    size = 3.5
  ) +
  
  theme_void()


# patway barplot
library(dplyr)

module_ptw <- read.csv("output/propagated_networks/enrichment_results/cluster2_ptw.csv")
top_module <- module_ptw %>% arrange(Entities.FDR) %>% pull(Pathway.identifier) %>% head(20)

module_ptw_top <- module_ptw %>% filter(Pathway.identifier %in% top_module)


embl_palette <- c("#5B2A86", "#B2ABD2", "#E6F5D0", "#1B7837")

module_plot <- ggplot(
  module_ptw_top %>%
    mutate(Pathway.name = forcats::fct_reorder(Pathway.name, -log(Entities.pValue))),
  aes(x = -log10(Entities.pValue),
      y = Pathway.name,
      fill = -log10(Entities.pValue))
) +
  geom_col(width = 0.8) +
  scale_fill_gradientn(colors = embl_palette) +
  labs(
    title = "Pathway enrichment of KRAS module nodes",
    x = expression(-log[10](p-value)),
    y = NULL,
    fill = expression(-log[10](p-value))
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 11),
    legend.position = "right"
  )

