################################################################################
# Unified Personalized PageRank Workflow for KRAS, NRAS, HRAS
################################################################################

library(igraph)
library(dplyr)
library(ggraph)
library(ggplot2)

set.seed(1234)

# -------------------------------------------------------------------------------
# 1) Load seeds
# -------------------------------------------------------------------------------

seeds <- read.csv("output/seed_weights_ann.csv")

# Assign UniProt IDs as names
seed_kras <- seeds$seed_kras; names(seed_kras) <- seeds$uniprotswissprot
seed_nras <- seeds$seed_nras; names(seed_nras) <- seeds$uniprotswissprot
seed_hras <- seeds$seed_hras; names(seed_hras) <- seeds$uniprotswissprot

# -------------------------------------------------------------------------------
# 2) Load networks
# -------------------------------------------------------------------------------

load_network <- function(file){
  nt <- read.csv(file)
  edges <- nt %>%
    dplyr::select(interactor_A, interactor_B) %>%
    filter(interactor_A != interactor_B)
  
  graph_from_data_frame(edges, directed = FALSE)
}

g_kras <- load_network("input/networks/kras_nt.csv")
g_nras <- load_network("input/networks/nras_nt.csv")
g_hras <- load_network("input/networks/hras_nt.csv")

# -------------------------------------------------------------------------------
# 3) PageRank + clustering function
# -------------------------------------------------------------------------------

run_ras_ppr <- function(g, seed_vector, topN = 300,
                        quant_thr = 0.85,
                        min_cluster_size = 10,
                        prefix = NULL){
  
  # ---------------------------------------------------------------------------
  # Robust seed selection
  # ---------------------------------------------------------------------------
  
  seed_vector_in_graph <- seed_vector[names(seed_vector) %in% V(g)$name]
  seed_vector_in_graph <- sort(seed_vector_in_graph, decreasing = TRUE)
  seed_vector_in_graph <- seed_vector_in_graph[1:min(topN, length(seed_vector_in_graph))]
  
  seed_nodes <- names(seed_vector_in_graph)
  
  cat(prefix, "- Seeds in graph:", length(seed_nodes), "\n")
  
  if(length(seed_nodes) == 0){
    stop(paste("No seed nodes found in graph for", prefix))
  }
  
  # ---------------------------------------------------------------------------
  # Normalize seed weights (robust for single seed)
  # ---------------------------------------------------------------------------
  
  if(length(seed_vector_in_graph) > 1){
    seed_vector_in_graph <- seed_vector_in_graph - min(seed_vector_in_graph)
    seed_vector_in_graph <- seed_vector_in_graph / max(seed_vector_in_graph)
  } else {
    seed_vector_in_graph <- rep(1, length(seed_vector_in_graph))
    names(seed_vector_in_graph) <- seed_nodes
  }
  
  # ---------------------------------------------------------------------------
  # Personalized vector
  # ---------------------------------------------------------------------------
  pers <- rep(0, vcount(g))
  names(pers) <- V(g)$name
  pers[seed_nodes] <- seed_vector_in_graph
  
  if(sum(pers) == 0){
    stop("Personalized vector is empty.")
  }
  
  pers <- pers / sum(pers)
  
  # ---------------------------------------------------------------------------
  # Personalized PageRank
  # ---------------------------------------------------------------------------
  
  ppr <- page_rank(g, personalized = pers)$vector
  
  # ---------------------------------------------------------------------------
  # Top subgraph selection
  # ---------------------------------------------------------------------------
  
  thr <- quantile(ppr, quant_thr, na.rm = TRUE)
  top_nodes <- names(ppr)[ppr >= thr]
  g_top <- induced_subgraph(g, vids = top_nodes)
  
  # ---------------------------------------------------------------------------
  # Clustering
  # ---------------------------------------------------------------------------
  
  wc <- cluster_walktrap(g_top)
  cl <- membership(wc)
  
  cl_nodes <- names(cl)
  seed_in_top <- intersect(seed_nodes, cl_nodes)
  cluster_list <- split(cl_nodes, cl)
  ppr_top <- ppr[cl_nodes]
  
  # ---------------------------------------------------------------------------
  # KS test per cluster
  # ---------------------------------------------------------------------------
  
  ks_p <- sapply(cluster_list, function(nodes_i){
    if(length(nodes_i) < min_cluster_size) return(NA_real_)
    ks.test(ppr_top[nodes_i], ppr_top, alternative="greater")$p.value
  })
  
  ks_fdr <- p.adjust(ks_p, method="BH")
  
  summary_tbl <- data.frame(
    cluster    = names(cluster_list),
    size       = sapply(cluster_list, length),
    n_seed     = sapply(cluster_list, function(x) sum(x %in% seed_in_top)),
    mean_ppr   = sapply(cluster_list, function(x) mean(ppr_top[x], na.rm=TRUE)),
    median_ppr = sapply(cluster_list, function(x) median(ppr_top[x], na.rm=TRUE)),
    ks_p       = ks_p[names(cluster_list)],
    ks_fdr     = ks_fdr[names(cluster_list)]
  )
  
  summary_tbl <- summary_tbl[summary_tbl$n_seed >= 2, ]
  
  return(list(
    ppr = ppr,
    g_top = g_top,
    clusters = cl,
    summary = summary_tbl,
    seeds = seed_nodes
  ))
}

# -------------------------------------------------------------------------------
# 4) Run workflow
# -------------------------------------------------------------------------------

result_kras <- run_ras_ppr(g_kras, seed_kras, topN=100, prefix="KRAS")
result_nras <- run_ras_ppr(g_nras, seed_nras, topN=100, prefix="NRAS")
result_hras <- run_ras_ppr(g_hras, seed_hras, topN=100, prefix="HRAS")

# -------------------------------------------------------------------------------
# 5) Network exploration function
# -------------------------------------------------------------------------------

explore_network <- function(g, result_obj, prefix, topN = 100,
                            hub_quantile = 0.95){
  
  seed_nodes <- result_obj$seeds
  
  # Full PPR ranking
  ppr_df <- data.frame(
    node = names(result_obj$ppr),
    ppr  = as.numeric(result_obj$ppr)
  ) %>% arrange(desc(ppr))
  
  # Direct bridge nodes
  neighbors_list <- neighborhood(g, order = 1, nodes = seed_nodes)
  neighbors_names <- lapply(neighbors_list, function(x) V(g)$name[x])
  bridge_direct <- setdiff(unique(unlist(neighbors_names)), seed_nodes)
  
  bridge_dir_df <- data.frame(
    node = bridge_direct,
    ppr  = result_obj$ppr[bridge_direct]
  ) %>% arrange(desc(ppr))
  
  # Indirect bridge nodes (robust)
  if(length(seed_nodes) >= 2){
    
    sp_matrix <- shortest_paths(g, from = seed_nodes,
                                to = V(g), output = "vpath")
    
    node_counts <- table(unlist(lapply(sp_matrix$vpath, names)))
    bridge_indirect <- setdiff(names(node_counts[node_counts >= 2]),
                               seed_nodes)
    
  } else {
    bridge_indirect <- character(0)
  }
  
  bridge_ind_df <- data.frame(
    node = bridge_indirect,
    ppr  = result_obj$ppr[bridge_indirect]
  ) %>% arrange(desc(ppr))
  
  # High-ranked nodes
  ppr_threshold <- quantile(result_obj$ppr, hub_quantile, na.rm = TRUE)
  
  high_ranked_df <- ppr_df %>%
    filter(ppr >= ppr_threshold)
  
  secondary_hubs_df <- high_ranked_df %>%
    filter(!(node %in% seed_nodes))
  
  return(list(
    ppr = ppr_df,
    bridge_direct = bridge_dir_df,
    bridge_indirect = bridge_ind_df,
    high_ranked = high_ranked_df,
    secondary_hubs = secondary_hubs_df
  ))
}

# -------------------------------------------------------------------------------
# 6) Run exploration
# -------------------------------------------------------------------------------

explore_kras <- explore_network(g_kras, result_kras, "kras")
explore_nras <- explore_network(g_nras, result_nras, "nras")
explore_hras <- explore_network(g_hras, result_hras, "hras")


clust <- result_kras[["clusters"]]
names(clust[clust == 2])
writeLines(names(clust[clust == 2]), "output/results/kras_cluster2_proteins.txt")
writeLines(explore_kras$bridge_direct$node, "output/results/kras_direct_bridge.txt")

write.csv(explore_hras[["bridge_direct"]], "output/results/hras_bridge_direct.csv")
