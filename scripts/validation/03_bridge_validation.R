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
core_nras <- union(shared_all, union(shared_kn, shared_nh))
core_hras <- union(shared_all, union(shared_kh, shared_nh))


brdge_kras <- bridge_enrichment(
  g_kras,
  explore_kras$bridge_direct$node,
  core_kras,
  "output/validation_results/kras_bridge_enrichment.tsv"
)

bridge_nras <- bridge_enrichment(
  g_nras,
  explore_nras$bridge_direct$node,
  core_nras,
  "output/validation_results/nras_bridge_enrichment.tsv"
)

bridge_hras <- bridge_enrichment(
  g_hras,
  explore_hras$bridge_direct$node,
  core_hras,
  "output/validation_results/hras_bridge_enrichment.tsv"
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
  
  wilcox.test(df$core_links, df$peripheral_links)
  
}

bridge_connection_test(
  g_hras,
  explore_hras$bridge_direct$node,
  core_hras,
  seed_hras
)

bridge_connection_test(
  g_kras,
  explore_kras$bridge_direct$node,
  core_kras,
  seed_kras
)

bridge_connection_test(
  g_nras,
  explore_nras$bridge_direct$node,
  core_nras,
  seed_nras
)


ecount(g_kras)
ecount(g_nras)
ecount(g_hras)

vcount(g_kras)
vcount(g_nras)
vcount(g_hras)


seed_density <- function(g, seeds){
  
  sub <- induced_subgraph(g, seeds)
  
  edge_density(sub)
  
}

seed_density(g_kras, seed_kras)
seed_density(g_nras, seed_nras)
seed_density(g_hras, seed_hras)
