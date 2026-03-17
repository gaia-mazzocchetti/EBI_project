################################################################################
# 01_seed_validation.R
################################################################################

source("scripts/validation/00_load_objects.R")
library(dplyr)

# ------------------------------------------------------------------
# Seed sets
# ------------------------------------------------------------------

seed_kras <- result_kras$seeds
seed_nras <- result_nras$seeds
seed_hras <- result_hras$seeds

shared_all <- Reduce(intersect, list(seed_kras, seed_nras, seed_hras))

shared_kn <- setdiff(intersect(seed_kras, seed_nras), shared_all)
shared_kh <- setdiff(intersect(seed_kras, seed_hras), shared_all)
shared_nh <- setdiff(intersect(seed_nras, seed_hras), shared_all)

unique_k <- setdiff(seed_kras, c(shared_all, shared_kn, shared_kh))
unique_n <- setdiff(seed_nras, c(shared_all, shared_kn, shared_nh))
unique_h <- setdiff(seed_hras, c(shared_all, shared_kh, shared_nh))

# ------------------------------------------------------------------
# Permutation overlap test
# ------------------------------------------------------------------

perm_overlap_test <- function(seed1, seed2, g1, g2, n_perm = 5000){
  
  universe <- intersect(V(g1)$name, V(g2)$name)
  
  seed1 <- intersect(seed1, universe)
  seed2 <- intersect(seed2, universe)
  
  obs <- length(intersect(seed1, seed2))
  
  rand <- replicate(n_perm,{
    r1 <- sample(universe, length(seed1))
    r2 <- sample(universe, length(seed2))
    length(intersect(r1,r2))
  })
  
  p <- (sum(rand >= obs)+1)/(n_perm+1)
  
  data.frame(
    observed_overlap = obs,
    mean_random = mean(rand),
    sd_random = sd(rand),
    p_value = p
  )
}

overlap_kn <- perm_overlap_test(seed_kras, seed_nras, g_kras, g_nras)
overlap_kh <- perm_overlap_test(seed_kras, seed_hras, g_kras, g_hras)
overlap_nh <- perm_overlap_test(seed_nras, seed_hras, g_nras, g_hras)

overlap_table <- rbind(
  data.frame(comparison="KRAS-NRAS", overlap_kn),
  data.frame(comparison="KRAS-HRAS", overlap_kh),
  data.frame(comparison="NRAS-HRAS", overlap_nh)
)

write.table(
  overlap_table,
    "output/validation_results/seed_overlap_results.tsv",
  sep="\t",
  row.names=FALSE,
  quote=FALSE
)

# ------------------------------------------------------------------
# Centrality analysis
# ------------------------------------------------------------------

compute_centrality <- function(g){
  
  data.frame(
    protein = V(g)$name,
    degree = degree(g),
    betweenness = betweenness(g, normalized = TRUE)
  )
}

analyze_seed_centrality <- function(g, seeds, outfile, save = T){
  
  cent <- compute_centrality(g)
  
  cent$group <- ifelse(
    cent$protein %in% seeds, "seed", "non_seed"
  )
  
  #cent <- cent[cent$protein %in% seeds, ]
  
  if(save) {
   write.table(
    cent,
    outfile,
    sep="\t",
    row.names=FALSE,
    quote=FALSE
  )
  }
  
  return(cent)
}

#KRAS
cent_kras <- analyze_seed_centrality(
  g_kras,
  seed_kras,
  "output/validation_results/kras_seed_centrality.tsv"
)

cent_kras$group2 <- ifelse(cent_kras$protein %in% unique_k, "paralog-specific",
                           ifelse(cent_kras$protein %in% shared_all, "shared", "other"))

cent_kras2 <- cent_kras %>% filter(group2 != "other") 
kras_test <- wilcox.test(betweenness ~ group2, data= cent_kras2 )
kras_btw <- cent_kras2 %>% group_by(group2) %>% summarize(mean_btw = mean(betweenness)) %>% mutate(prot = "KRAS",
                                                                                                   wilcox_p= kras_test$p.value)

#NRAS
cent_nras <- analyze_seed_centrality(
  g_nras,
  seed_nras,
  "output/validation_results/nras_seed_centrality.tsv"
)

cent_nras$group2 <- ifelse(cent_nras$protein %in% unique_n, "paralog-specific",
                           ifelse(cent_nras$protein %in% shared_all, "shared", "other"))

cent_nras2 <- cent_nras %>% filter(group2 != "other") 
nras_test <- wilcox.test(betweenness ~ group2, data= cent_nras2 )
nras_btw <- cent_nras2 %>% group_by(group2) %>% summarize(mean_btw = mean(betweenness)) %>% mutate(prot = "NRAS",
                                                                                                   wilcox_p= nras_test$p.value)


#HRAS
cent_hras <- analyze_seed_centrality(
  g_hras,
  seed_hras,
  "output/validation_results/hras_seed_centrality.tsv"
)


cent_hras$group2 <- ifelse(cent_hras$protein %in% unique_h, "paralog-specific",
                           ifelse(cent_hras$protein %in% shared_all, "shared", "other"))

cent_hras2 <- cent_hras %>% filter(group2 != "other") 
hras_test <- wilcox.test(betweenness ~ group2, data= cent_hras2 )
hras_btw <- cent_hras2 %>% group_by(group2) %>% summarize(mean_btw = mean(betweenness)) %>% mutate(prot = "HRAS",
                                                                                                   wilcox_p= hras_test$p.value)



btw <- rbind(kras_btw, nras_btw, hras_btw)

write.table(
  btw,
  "output/validation_results/core_shared_seed_betweenness_means.tsv",
  sep="\t",
  row.names=FALSE,
  quote=FALSE
)

sink("output/validation_results/wilcoxon_seed_betweenness.txt")
kras_test <- wilcox.test(betweenness ~ group, data= cent_kras)
wilcox.test(betweenness ~ group, data= cent_nras )
wilcox.test(betweenness ~ group, data= cent_hras )
sink()
