################################################################################
# 03_bridge_validation.R
################################################################################

source("scripts/validation/00_load_objects.R")

seed_kras <- result_kras$seeds
seed_nras <- result_nras$seeds
seed_hras <- result_hras$seeds

shared_all <- Reduce(intersect, list(seed_kras, seed_nras, seed_hras))

shared_kn <- setdiff(intersect(seed_kras, seed_nras), shared_all)
shared_kh <- setdiff(intersect(seed_kras, seed_hras), shared_all)
shared_nh <- setdiff(intersect(seed_nras, seed_hras), shared_all)

bridge_enrichment <- function(g, bridge, core, outfile, save = T){
  
  universe <- V(g)$name
  bridge <- intersect(bridge, universe)
  core <- intersect(core, universe)
  
  bridge_core <- intersect(bridge, core)
  
  non_core <- setdiff(universe, core)
  
  p <- phyper(
    length(bridge_core) - 1,
    m = length(core),
    n = length(non_core),
    k = length(bridge),
    lower.tail = FALSE
  )
  
  res <- data.frame(
    universe_size = length(universe),
    n_bridge = length(bridge),
    n_core = length(core),
    bridge_in_core = length(bridge_core),
    p_value = p
  )
  
  if(save) {    
    write.table(
      res,
      outfile,
      sep="\t",
      row.names=FALSE,
      quote=FALSE
    )
  }
  
return(res)
  
}

core_kras <- union(shared_all, union(shared_kn, shared_kh))

brdge_kras <- bridge_enrichment(
  g_kras,
  explore_kras$bridge_direct$node,
  shared_all,
  "output/validation_results/kras_bridge_enrichment.tsv", save =F
)

bridge_nras <- bridge_enrichment(
  g_nras,
  explore_nras$bridge_direct$node,
  shared_all,
  "output/validation_results/nras_bridge_enrichment.tsv",  save =F
)

bridge_hras <- bridge_enrichment(
  g_hras,
  explore_hras$bridge_direct$node,
  shared_all,
  "output/validation_results/hras_bridge_enrichment.tsv",  save =F
)



bridge_connection_test <- function(g, bridge, core, seeds){
  
  peripheral <- setdiff(seeds, core)
  
  core_links <- sapply(bridge, function(n){
    
    neigh <- neighbors(g, n)
    neigh_names <- V(g)$name[neigh]
    
    sum(neigh_names %in% core)
    
  })
  
  peripheral_links <- sapply(bridge, function(n){
    
    neigh <- neighbors(g, n)
    neigh_names <- V(g)$name[neigh]
    
    sum(neigh_names %in% peripheral)
    
  })
  
  df <- data.frame(
    bridge = bridge,
    core_links = core_links,
    peripheral_links = peripheral_links
  )
  
  
  df$core_links_norm <- df$core_links / length(core)
  df$peripheral_links_norm <- df$peripheral_links / length(peripheral)
  
  wl <- wilcox.test(df$core_links_norm, df$peripheral_links_norm)
  
  res <- list(test =wl, links = df)
  return(res)
  
}

kras_neigh <- bridge_connection_test(
  g_kras,
  explore_kras$bridge_direct$node,
  shared_all,
  seed_kras
)

core_kras <- union(shared_all, union(shared_kn, shared_kh))
kras_neigh_more <- bridge_connection_test(
  g_kras,
  explore_kras$bridge_direct$node,
  core_kras,
  seed_kras
)

sink("output/validation_results/mean_links_bridge_nodes.txt")
mean(kras_neigh$links$core_links)
mean(kras_neigh$links$peripheral_links)
kras_neigh$test

mean(kras_neigh_more$links$core_links)
mean(kras_neigh_more$links$peripheral_links)
kras_neigh_more$test
sink()
