library(igraph)
library(dplyr)

# -------------------------------------------------------------------------------
# 1) Load networks
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
# 2) Compute global degree
# -------------------------------------------------------------------------------

deg_kras <- degree(g_kras)
deg_hras <- degree(g_hras)
deg_nras <- degree(g_nras)

degree_df_kras <- data.frame(
  protein = names(deg_kras),
  degree = as.numeric(deg_kras)
)

degree_df_hras <- data.frame(
  protein = names(deg_hras),
  degree = as.numeric(deg_hras)
)

degree_df_nras <- data.frame(
  protein = names(deg_nras),
  degree = as.numeric(deg_nras)
)

# -------------------------------------------------------------------------------
# 3) Category definition
# -------------------------------------------------------------------------------
seed_hras_100 <- result_hras$seeds
seed_kras_100 <- result_kras$seeds
seed_nras_100 <- result_nras$seeds
# Intersezioni
shared_all <- Reduce(intersect, list(seed_kras_100,
                                     seed_nras_100,
                                     seed_hras_100))

shared_kh <- intersect(seed_kras_100, seed_hras_100)
shared_kn <- intersect(seed_kras_100, seed_nras_100)
shared_nh <- intersect(seed_nras_100, seed_hras_100)

# Rimuoviamo quelli shared_all dalle pairwise
shared_kh <- setdiff(shared_kh, shared_all)
shared_kn <- setdiff(shared_kn, shared_all)
shared_nh <- setdiff(shared_nh, shared_all)

# Unique
unique_k <- setdiff(seed_kras_100,
                    union(union(shared_kh, shared_kn), shared_all))

unique_n <- setdiff(seed_nras_100,
                    union(union(shared_kn, shared_nh), shared_all))

unique_h <- setdiff(seed_hras_100,
                    union(union(shared_kh, shared_nh), shared_all))


seed_category <- data.frame(
  protein = c(shared_all,
              shared_kh,
              shared_kn,
              shared_nh,
              unique_k,
              unique_n,
              unique_h),
  category = c(rep("shared_all", length(shared_all)),
               rep("shared_kh", length(shared_kh)),
               rep("shared_kn", length(shared_kn)),
               rep("shared_nh", length(shared_nh)),
               rep("unique_k", length(unique_k)),
               rep("unique_n", length(unique_n)),
               rep("unique_h", length(unique_h)))
)

# -------------------------------------------------------------------------------
# 4) Statistical tests
# -------------------------------------------------------------------------------
seed_degree_kras <- seed_category %>%
  left_join(degree_df_kras, by = "protein")

seed_degree_kras <- seed_degree_kras %>%
  mutate(group = ifelse(grepl("shared", category),
                        "shared",
                        "unique"))

seed_degree_hras <- seed_category %>%
  left_join(degree_df_hras, by = "protein")

seed_degree_hras <- seed_degree_hras %>%
  mutate(group = ifelse(grepl("shared", category),
                        "shared",
                        "unique"))

seed_degree_nras <- seed_category %>%
  left_join(degree_df_nras, by = "protein")

seed_degree_nras <- seed_degree_nras %>%
  mutate(group = ifelse(grepl("shared", category),
                        "shared",
                        "unique"))

wilcox.test(degree ~ group, data = seed_degree_kras)
wilcox.test(degree ~ group, data = seed_degree_nras)
wilcox.test(degree ~ group, data = seed_degree_hras)

kruskal.test(degree ~ category, data = seed_degree_kras)
pairwise.wilcox.test(seed_degree_kras$degree,
                     seed_degree_kras$category,
                     p.adjust.method = "BH")

kruskal.test(degree ~ category, data = seed_degree_nras)
pairwise.wilcox.test(seed_degree_nras$degree,
                     seed_degree_nras$category,
                     p.adjust.method = "BH")

kruskal.test(degree ~ category, data = seed_degree_hras)
pairwise.wilcox.test(seed_degree_hras$degree,
                     seed_degree_hras$category,
                     p.adjust.method = "BH")


