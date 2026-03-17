################################################################################
# 00_load_objects.R
# Load networks and workflow outputs
################################################################################

library(igraph)

load_network <- function(file){
  nt <- readr::read_csv(file)
  edges <- nt %>%
    dplyr::select(interactor_A, interactor_B) %>%
    dplyr::filter(interactor_A != interactor_B)
  
  graph_from_data_frame(edges, directed = FALSE)
}

# Networks
g_kras <- load_network("input/networks/kras_nt.csv")
g_nras <- load_network("input/networks/nras_nt.csv")
g_hras <- load_network("input/networks/hras_nt.csv")

# PPR results
result_kras <- readRDS("output/propagated_networks/result_kras.rds")
result_nras <- readRDS("output/propagated_networks/result_nras.rds")
result_hras <- readRDS("output/propagated_networks/result_hras.rds")

# Exploration results
explore_kras <- readRDS("output/propagated_networks/explore_kras.rds")
explore_nras <- readRDS("output/propagated_networks/explore_nras.rds")
explore_hras <- readRDS("output/propagated_networks/explore_hras.rds")