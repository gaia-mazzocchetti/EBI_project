
library(igraph)
library(ggraph)
# Plot del network
source("scripts/validation/00_load_objects.R")

seed_kras <- result_kras$seeds
seed_nras <- result_nras$seeds
seed_hras <- result_hras$seeds

core_proteins <- c("P55011","P05023","Q8NF37","P98172","Q00587","P52799","Q9P2W9")

# KRAS 
V(g_kras)$type <- "other"

V(g_kras)$type[V(g_kras)$name %in% seed_kras] <- "seed"
V(g_kras)$type[V(g_kras)$name %in% core_proteins] <- "core"
V(g_kras)$type[V(g_kras)$name == "P01116"] <- "KRAS"
V(g_kras)$type[V(g_kras)$name == "P01116-1"] <- "KRAS-1"
V(g_kras)$type[V(g_kras)$name == "P01116-2"] <- "KRAS-2"


set.seed(123)

a <- ggraph(g_kras, layout = "stress") +
  
  geom_edge_link(
    alpha = 0.25,
    color = "grey80"
  ) +
  
  geom_node_point(
    aes(color = type),
    size = 4,
    alpha= 0.6
  ) +
  
  geom_node_text(
    aes(label = ifelse(type %in% c("core","KRAS", "KRAS-1", "KRAS-2"), name, "")),
    repel = TRUE,
    size = 3
  ) +
  
  scale_color_manual(values = c(
    KRAS = "#862d86",
    core = "#9999ff",
    seed = "#40bf40",
    other = "grey85"
  )) +
  
  theme_void() +
  
  labs(
    title = "KRAS network with seed distribution",
    color = "Node type"
  )


seed_nodes <- intersect(c(ssed_kras, "P01116"), V(g_kras)$name)
seed_ids <- which(V(g_kras)$name %in% seed_nodes)

neighbors_seed <- unlist(neighborhood(g_kras, order = 1, nodes = seed_ids))
nodes_keep <- unique(c(seed_ids, neighbors_seed))

g_seed <- induced_subgraph(g_kras, nodes_keep)

deg <- degree(g_seed)

g_seed <- delete_vertices(g_seed, V(g_seed)[deg == 1 & !(name %in% seed_nodes)])

V(g_seed)$label <- ifelse(
  V(g_seed)$type %in% c("seed","core","KRAS", "KRAS-1", "KRAS-2"),
  V(g_seed)$name,
  ""
)
set.seed(123)

a <- ggraph(g_seed, layout = "stress") +
  
  geom_edge_link(
    color="grey80",
    alpha=0.4
  ) +
  
  geom_node_point(
    aes(color=type),
    size=4
  ) +
  
  geom_node_text(
    aes(label = ifelse(type %in% c("core","KRAS", "KRAS-1", "KRAS-2"), name, "")),
    repel=TRUE,
    size=3
  ) +
  
  scale_color_manual(values=c(
    KRAS="black",
    core="#E69F00",
    seed="#D55E00",
    other="grey85"
  )) +
  
  theme_void() +
  labs(
    title = "KRAS network with seed distribution",
    color = "Node type"
  )


# NRAS 
V(g_nras)$type <- "other"

V(g_nras)$type[V(g_nras)$name %in% seed_nras] <- "seed"
V(g_nras)$type[V(g_nras)$name %in% core_proteins] <- "core"
V(g_nras)$type[V(g_nras)$name == "P01111"] <- "NRAS"

set.seed(123)

b <- ggraph(g_nras, layout = "stress") +
  
  geom_edge_link(
    alpha = 0.25,
    color = "grey80"
  ) +
  
  geom_node_point(
    aes(color = type),
    size = 4,
    alpha= 0.6
  ) +
  
  geom_node_text(
    aes(label = ifelse(type %in% c("core","NRAS"), name, "")),
    repel = TRUE,
    size = 3
  ) +
  
  scale_color_manual(values = c(
    NRAS = "#862d86",
    core = "#9999ff",
    seed = "#40bf40",
    other = "grey85"
  )) +
  
  theme_void() +
  
  labs(
    title = "NRAS network with seed distribution",
    color = "Node type"
  )


# HRAS 
V(g_hras)$type <- "other"

V(g_hras)$type[V(g_hras)$name %in% seed_hras] <- "seed"
V(g_hras)$type[V(g_hras)$name %in% core_proteins] <- "core"
V(g_hras)$type[V(g_hras)$name == "P01112"] <- "HRAS"

set.seed(123)

c <- ggraph(g_hras, layout = "stress") +
  
  geom_edge_link(
    alpha = 0.25,
    color = "grey80"
  ) +
  
  geom_node_point(
    aes(color = type),
    size = 4,
    alpha= 0.6
  ) +
  
  geom_node_text(
    aes(label = ifelse(type %in% c("core","HRAS"), name, "")),
    repel = TRUE,
    size = 3
  ) +
  
  scale_color_manual(values = c(
    HRAS = "#862d86",
    core = "#9999ff",
    seed = "#40bf40",
    other = "grey85"
  )) +
  
  theme_void() +
  
  labs(
    title = "HRAS network with seed distribution",
    color = "Node type"
  )


 a + b + c

 
 # Venn diagram seeds
 seed_list <- list(
   KRAS = seed_kras,
   NRAS = seed_nras,
   HRAS = seed_hras
 )

 library(ggVennDiagram)
 library(ggplot2)
 library(patchwork)
 
 p_venn <- ggVennDiagram(seed_list) +
   scale_fill_gradient(high = "#008B45", low = "#7CCD7C") +
   theme_void()+
   theme(legend.position = "none")
 
 library(patchwork)
 (a | b | c) /
   p_venn
 
 
 # patway barplot
 library(dplyr)
 core_ptw <- read.csv("output/propagated_networks/enrichment_results/core_enrichment_ptw.csv")
 top_core <- core_ptw %>% arrange(Entities.FDR) %>% pull(Pathway.identifier) %>% head(20)
 
 core_ptw_top <- core_ptw %>% filter(Pathway.identifier %in% top_core)
 
 
 embl_palette <- c("#5E3C99", "#B2ABD2", "#E6F5D0", "#1B7837")
 
 core_plot <- ggplot(
   core_ptw_top %>%
     mutate(Pathway.name = forcats::fct_reorder(Pathway.name, -log(Entities.pValue))),
   aes(x = -log10(Entities.pValue),
       y = Pathway.name,
       fill = -log10(Entities.pValue))
 ) +
   geom_col(width = 0.8) +
   scale_fill_gradientn(colors = embl_palette) +
   labs(
     title = "Pathway enrichment of core nodes",
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
 