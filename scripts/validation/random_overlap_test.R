


universe_kn <- intersect(V(g_kras)$name,
                         V(g_nras)$name)
universe_kh <- intersect(V(g_kras)$name,
                         V(g_hras)$name)
universe_nh <- intersect(V(g_nras)$name,
                         V(g_hras)$name)
universe_all <- Reduce(intersect,
                       list(V(g_kras)$name,
                            V(g_nras)$name,
                            V(g_hras)$name))

length(V(g_kras))
length(V(g_nras))
length(V(g_hras))

length(universe_kn)
length(universe_kh)
length(universe_nh)
length(universe_all)


perm_overlap_test <- function(seed1, seed2, g1, g2, 
                              n_perm = 5000, 
                              seed_random = 123) {
  
  set.seed(seed_random)
  
  # Universo corretto = intersezione dei nodi
  universe <- intersect(V(g1)$name, V(g2)$name)
  
  # Restringi seed all’universo comune
  seed1_u <- intersect(seed1, universe)
  seed2_u <- intersect(seed2, universe)
  
  # Overlap osservato
  obs_overlap <- length(intersect(seed1_u, seed2_u))
  
  # Permutazioni
  rand_overlap <- numeric(n_perm)
  
  for(i in seq_len(n_perm)){
    rand1 <- sample(universe, length(seed1_u))
    rand2 <- sample(universe, length(seed2_u))
    rand_overlap[i] <- length(intersect(rand1, rand2))
  }
  
  # p-value empirico
  p_emp <- (sum(rand_overlap >= obs_overlap) + 1) / (n_perm + 1)
  
  # Output strutturato
  list(
    observed_overlap = obs_overlap,
    mean_random = mean(rand_overlap),
    sd_random = sd(rand_overlap),
    p_value = p_emp,
    random_distribution = rand_overlap
  )
}

# KRAS-NRAS
res_kn <- perm_overlap_test(seed_kras_100,
                            seed_nras_100,
                            g_kras,
                            g_nras)

# KRAS-HRAS
res_kh <- perm_overlap_test(seed_kras_100,
                            seed_hras_100,
                            g_kras,
                            g_hras)

# NRAS-HRAS
res_nh <- perm_overlap_test(seed_nras_100,
                            seed_hras_100,
                            g_nras,
                            g_hras)

# ALL THREE PARALOGS
perm_overlap_test_3 <- function(seed1, seed2, seed3,
                                g1, g2, g3,
                                n_perm = 5000,
                                seed_random = 123) {
  
  set.seed(seed_random)
  
  universe <- Reduce(intersect,
                     list(V(g1)$name,
                          V(g2)$name,
                          V(g3)$name))
  
  seed1_u <- intersect(seed1, universe)
  seed2_u <- intersect(seed2, universe)
  seed3_u <- intersect(seed3, universe)
  
  obs_overlap <- length(Reduce(intersect,
                               list(seed1_u,
                                    seed2_u,
                                    seed3_u)))
  
  rand_overlap <- numeric(n_perm)
  
  for(i in seq_len(n_perm)){
    rand1 <- sample(universe, length(seed1_u))
    rand2 <- sample(universe, length(seed2_u))
    rand3 <- sample(universe, length(seed3_u))
    
    rand_overlap[i] <- length(Reduce(intersect,
                                     list(rand1, rand2, rand3)))
  }
  
  p_emp <- (sum(rand_overlap >= obs_overlap) + 1) / (n_perm + 1)
  
  list(
    observed_overlap = obs_overlap,
    mean_random = mean(rand_overlap),
    sd_random = sd(rand_overlap),
    p_value = p_emp,
    random_distribution = rand_overlap
  )
}

res_all <- perm_overlap_test_3(seed_kras_100,
                               seed_nras_100,
                               seed_hras_100,
                               g_kras,
                               g_nras,
                               g_hras)
print_result <- function(res){
  cat("Observed overlap:", res$observed_overlap, "\n")
  cat("Mean random:", round(res$mean_random,2), "\n")
  cat("SD random:", round(res$sd_random,2), "\n")
  cat("Empirical p-value:", res$p_value, "\n\n")
}

print_result(res_kn)
print_result(res_kh)
print_result(res_nh)
print_result(res_all)
