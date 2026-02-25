

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

core_all <- shared_all
core_pairwise_k <- union(shared_kh, shared_kn)

core_kras <- union(core_all, core_pairwise_k)

peripheral_kras <- setdiff(seed_kras_100, core_kras)

library(data.table)
bridge_kras <- explore_kras[["bridge_direct"]]$node
bridge_nras <- explore_nras[["bridge_direct"]]$node
bridge_hras <- explore_hras[["bridge_direct"]]$node

table(bridge_kras %in% core_kras)

bridge_core_k <- intersect(bridge_kras, core_kras)
bridge_peripheral_k <- intersect(bridge_kras, peripheral_kras)

length(bridge_core_k)
length(bridge_peripheral_k)


phyper(k-1,
       m = length(peripheral_kras),
       n = length(core_kras),
       k = length(bridge_kras),
       lower.tail = FALSE)


core_pairwise_h <- union(shared_kh, shared_nh)
core_hras <- union(core_all, core_pairwise_h)
bridge_core_h <- intersect(bridge_hras, core_hras)
peripheral_hras <- setdiff(seed_hras_100, core_hras)
bridge_peripheral_h <- intersect(bridge_hras, peripheral_hras)

length(bridge_core_h)
length(bridge_peripheral_h)

core_pairwise_n <- union(shared_kn, shared_nh)
core_nras <- union(core_all, core_pairwise_n)
bridge_core_n <- intersect(bridge_nras, core_nras)
peripheral_nras <- setdiff(seed_nras_100, core_nras)
bridge_peripheral_n <- intersect(bridge_nras, peripheral_nras)

length(bridge_core_n)
length(bridge_peripheral_n)

