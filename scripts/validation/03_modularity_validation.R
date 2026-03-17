################################################################################
# 02_modularity_validation.R
################################################################################

source("scripts/validation/00_load_objects.R")

compute_modularity <- function(g){
  
  cl <- cluster_louvain(g)
  
  data.frame(
    nodes = vcount(g),
    edges = ecount(g),
    modularity = modularity(cl),
    n_clusters = length(unique(membership(cl)))
  )
}

mod_kras <- compute_modularity(g_kras)
mod_nras <- compute_modularity(g_nras)
mod_hras <- compute_modularity(g_hras)

modularity_table <- rbind(
  data.frame(paralog="KRAS", mod_kras),
  data.frame(paralog="NRAS", mod_nras),
  data.frame(paralog="HRAS", mod_hras)
)

write.table(
  modularity_table,
  "output/validation_results/network_modularity.tsv",
  sep="\t",
  row.names=FALSE,
  quote=FALSE
)
