################################################################################
# Comparative Network Propagation: Top-100 vs Top-300 Seeds
################################################################################

library(igraph)
library(dplyr)
library(biomaRt)
set.seed(1234)
# -------------------------------------------------------------------------------
# 1) Load seed weights and map HUGO symbols to UniProt IDs
# -------------------------------------------------------------------------------

seeds <- read.csv("output/seed_weights_ann.csv")

# -------------------------------------------------------------------------------
# 2) Load the KRAS interaction network from IntAct
# -------------------------------------------------------------------------------

kras_nt <- read.csv("input//networks/kras_nt.csv")

# Edges for the igraph object
edges <- kras_nt %>% dplyr::select(interactor_A, interactor_B)

# Remove self‑edges
edges <- edges[edges$interactor_A != edges$interactor_B, ]

# Create undirected graph from edge list
g_kras <- graph_from_data_frame(d = edges, directed = FALSE)

# -------------------------------------------------------------------------------
# 3) Prepare personalized PageRank seed vector
# -------------------------------------------------------------------------------

# Extract the seed weights specific to KRAS
seed_vector <- seeds$seed_kras
names(seed_vector) <- seeds$uniprotswissprot  # assign UniProt IDs

# ----------------------------
# Function: PPR + clustering with stability controls
# ----------------------------
run_stable_ppr <- function(g, seed_vector, min_top=100, max_top=200, step=20, 
                           quant_thr=0.85, min_cluster_size=10){
  
  results <- list()
  
  for(topN in seq(min_top, max_top, by=step)){
    
    # Take top-N seeds
    seed_vector_top <- sort(seed_vector, decreasing=TRUE)[1:min(topN, length(seed_vector))]
    common_genes <- intersect(names(seed_vector_top), V(g)$name)
    seed_vector_top <- seed_vector_top[common_genes]
    
    if(length(seed_vector_top)==0){
      message(sprintf("Top-%d seeds not present in graph, skipping...", topN))
      next
    }
    
    # Normalize to 0-1
    seed_vector_top <- seed_vector_top - min(seed_vector_top)
    seed_vector_top <- seed_vector_top / max(seed_vector_top)
    
    # Build personalized vector
    pers <- rep(0, vcount(g))
    names(pers) <- V(g)$name
    pers[common_genes] <- seed_vector_top
    pers <- pers / sum(pers)
    
    # PageRank
    ppr <- page_rank(g, personalized = pers)$vector
    
    # Top nodes subgraph
    thr <- quantile(ppr, quant_thr, na.rm=TRUE)
    top_nodes <- names(ppr)[ppr >= thr]
    g_top <- induced_subgraph(g, vids=top_nodes)
    
    # Clustering
    wc <- cluster_walktrap(g_top)
    cl <- membership(wc)
    cl_nodes <- names(cl)
    
    # Count seed per cluster
    seed_in_top <- intersect(names(seed_vector_top), cl_nodes)
    n_seed_per_cluster <- tapply(cl_nodes %in% seed_in_top, cl, sum)
    
    # Keep only clusters with at least min_cluster_size nodes and 2 seeds
    keep_clusters <- names(n_seed_per_cluster[n_seed_per_cluster >= 2])
    
    # Summary
    cluster_list <- split(cl_nodes, cl)
    summary_tbl <- data.frame(
      cluster     = names(cluster_list),
      size        = sapply(cluster_list, length),
      n_seed      = sapply(cluster_list, function(x) sum(x %in% seed_in_top)),
      mean_ppr    = sapply(cluster_list, function(x) mean(ppr[x], na.rm=TRUE)),
      median_ppr  = sapply(cluster_list, function(x) median(ppr[x], na.rm=TRUE))
    )
    summary_tbl <- summary_tbl[summary_tbl$cluster %in% keep_clusters, ]
    
    results[[paste0("topN_", topN)]] <- list(
      ppr=ppr,
      g_top=g_top,
      clusters=cl,
      summary=summary_tbl
    )
    
    message(sprintf("Top-%d seeds: %d clusters retained", topN, nrow(summary_tbl)))
  }
  
  return(results)
}



# Run the stable propagation workflow
ppr_results <- run_stable_ppr(
  g = g_kras,
  seed_vector = seed_vector,
  min_top = 100,
  max_top = 350,
  step = 25,
  quant_thr = 0.85,
  min_cluster_size = 10
)


# Example: compare clusters between topN=100 and topN=125
summary_250 <- ppr_results$topN_250$summary
summary_300 <- ppr_results$topN_300$summary

