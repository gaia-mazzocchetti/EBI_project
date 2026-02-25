################################################################################
# Personalized PageRank Network Propagation (KRAS)
################################################################################

library(igraph)
library(dplyr)
library(biomaRt)

# -------------------------------------------------------------------------------
# 1) Load seed weights and map HUGO symbols to UniProt IDs
# -------------------------------------------------------------------------------

seeds <- read.csv("output/seed_weights.csv")

# Extract Hugo gene symbols
hugo_genes <- seeds$X

# Connect to Ensembl for mapping
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Query UniProt Swiss‑Prot IDs for Hugo symbols
mapping <- getBM(
  attributes = c("hgnc_symbol", "uniprotswissprot"),
  filters = "hgnc_symbol",
  values = hugo_genes,
  mart = ensembl
)
mapping <- mapping[mapping$uniprotswissprot != "",]
# Join mapping back to seeds
seeds <- seeds %>% 
  left_join(mapping, by = c("X" = "hgnc_symbol"))

# Keep only seeds with valid UniProt IDs (needed for network alignment)
seeds2 <- seeds[!is.na(seeds$uniprotswissprot),]

write_csv(seeds2, "output/seed_weights_ann.csv")
# -------------------------------------------------------------------------------
# 2) Load the KRAS interaction network from IntAct
# -------------------------------------------------------------------------------

kras_nt <- read.csv("output/networks/kras_nt.csv")

# Edges for the igraph object
edges <- kras_nt %>% select(interactor_A, interactor_B)

# Remove self‑edges
edges <- edges[edges$interactor_A != edges$interactor_B, ]

# Create undirected graph from edge list
g_kras <- graph_from_data_frame(d = edges, directed = FALSE)

# -------------------------------------------------------------------------------
# 3) Prepare personalized PageRank seed vector
# -------------------------------------------------------------------------------

# Extract the seed weights specific to KRAS
seed_vector <- seeds2$seed_kras
names(seed_vector) <- seeds2$uniprotswissprot  # assign UniProt IDs

# Only keep seeds present in the KRAS graph
common_genes <- intersect(names(seed_vector), V(g_kras)$name)
seed_vector <- seed_vector[common_genes]

# If too few seeds remain, you may choose top N seeds to sharpen the signal
topN <- 300
if(length(seed_vector) > topN){
  seed_vector <- sort(seed_vector, decreasing = TRUE)[1:topN]
  common_genes  <- names(seed_vector)
}

# Rescale seed weights to a positive range [0,1]
seed_vector <- seed_vector - min(seed_vector)
seed_vector <- seed_vector / max(seed_vector)

# Build the personalization vector with zeros for non‑seed nodes
pers <- rep(0, vcount(g_kras))
names(pers) <- V(g_kras)$name

# Assign seed weights only to common genes
pers[common_genes] <- seed_vector

# Normalize to sum equal 1 (required for igraph personalized PageRank)
pers <- pers / sum(pers)

# -------------------------------------------------------------------------------
# 4) Compute Personalized PageRank
# -------------------------------------------------------------------------------

ppr <- page_rank(g_kras, personalized = pers)$vector

# -------------------------------------------------------------------------------
# 5) Extract high‑ranking subgraph (e.g., top 15%)
# -------------------------------------------------------------------------------

thr <- quantile(ppr, 0.85, na.rm = TRUE)
top_nodes <- names(ppr)[ppr >= thr]
g_top <- induced_subgraph(g_kras, vids = top_nodes)

# -------------------------------------------------------------------------------
# 6) Cluster the high‑score subgraph
# -------------------------------------------------------------------------------

wc <- cluster_walktrap(g_top)
cl  <- membership(wc)

# Prepare a data frame of PageRank scores and cluster membership
out <- data.frame(node  = names(ppr),
                  ppr   = as.numeric(ppr),
                  cluster = NA_integer_)

out$cluster[match(names(cl), out$node)] <- cl

# Only keep clusters with at least k nodes (e.g., 10)
cl_size        <- table(cl)
keep_clusters  <- names(cl_size[cl_size >= 10])
out2           <- out[out$cluster %in% keep_clusters, ]

# -------------------------------------------------------------------------------
# 7) Assess cluster significance (optional)
# -------------------------------------------------------------------------------

# Extract only the PPR values for nodes in the top subgraph
ppr_top <- ppr[top_nodes]

ks_p <- sapply(split(top_nodes, cl), function(nodes_i) {
  if (length(nodes_i) < 10) return(NA_real_)
  ks.test(ppr_top[nodes_i], ppr_top, alternative = "greater")$p.value
})

ks_fdr <- p.adjust(ks_p, method = "BH")

# -------------------------------------------------------------------------------
# 8) Build summary table of clusters
# -------------------------------------------------------------------------------

cluster_list <- split(top_nodes, cl)

summary_tbl <- data.frame(
  cluster     = names(cluster_list),
  size        = sapply(cluster_list, length),
  n_seed      = sapply(cluster_list, function(x) sum(x %in% common_genes)),
  mean_ppr    = sapply(cluster_list, function(x) mean(ppr_top[x], na.rm = TRUE)),
  median_ppr  = sapply(cluster_list, function(x) median(ppr_top[x], na.rm = TRUE)),
  ks_p        = ks_p[names(cluster_list)],
  ks_fdr      = ks_fdr[names(cluster_list)]
)

# Order by adjusted p‑value and effect size
summary_tbl <- summary_tbl[order(summary_tbl$ks_fdr, -summary_tbl$mean_ppr), ]

# Print summary
print(summary_tbl)


