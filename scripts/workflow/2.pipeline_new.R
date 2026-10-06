################################################################################
# Unified Personalized PageRank Workflow for KRAS, NRAS, HRAS  (revised)
#
# Changes with respect to the previous version:
#  1. KS test direction fixed: alternative = "less" (= cluster values LARGER),
#     background excludes the cluster itself, seeds excluded from both sides.
#     A one-sided Wilcoxon test is reported alongside as a sanity check.
#  2. Bridge nodes: shortest paths are now computed between PAIRS OF SEEDS
#     (shortest_paths() only accepts a single 'from' vertex). First neighbours
#     of seeds are reported separately as "seed_neighbours".
#  3. Seed weights: normalised by their sum (the old min-max scaling gave the
#     lowest seed a weight of 0). Option to select seeds by |score| so that
#     down-regulated genes are not silently discarded.
#  4. Network cleaning: UniProt isoform / chain IDs mapped to the canonical
#     accession, non-protein interactors removed, multi-edges and self-loops
#     removed with simplify().
#  5. Hub bias: degree-matched seed permutations give, for every node, a
#     z-score and an empirical p-value. Node selection and cluster tests use
#     the z-score instead of the raw PPR.
#  6. Cluster statistics are reported on non-seed nodes as well.
#  7. Clusters are filtered (size, n_seed) BEFORE multiple-testing correction.
#  8. New: Jaccard overlap of clusters between paralogs.
################################################################################

suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
})

set.seed(1234)

# -------------------------------------------------------------------------------
# 0) Parameters
# -------------------------------------------------------------------------------

params <- list(
  seeds_file       = "output/DEG/seed_weights_ann.csv",
  network_file     = "input/all_cancer_interactome.csv",
  out_dir          = "output2/propagated_networks/",
  topN             = 100,    # number of seeds per paralog
  seed_use_abs     = TRUE,   # rank seeds by |score| (keeps down-regulated genes)
  damping          = 0.85,
  n_perm           = 500,    # degree-matched permutations (>= 1000 for the paper)
  n_degree_bins    = 20,
  score_type       = "z",    # "z" (degree-corrected) or "ppr" (raw)
  quant_thr        = 0.85,   # keep top 15% of nodes by score
  min_cluster_size = 10,
  min_cluster_seed = 2
)

dir.create(params$out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------------
# 1) Helpers: UniProt IDs
# -------------------------------------------------------------------------------

uniprot_regex <- "^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"

# P01116-1 -> P01116 ; P05067-PRO_0000000092 -> P05067 ; uniprotkb:P01112 -> P01112
canonical_uniprot <- function(x){
  x <- trimws(as.character(x))
  x <- sub("^uniprotkb:", "", x)
  sub("-(PRO_)?[0-9]+$", "", x)
}

# -------------------------------------------------------------------------------
# 2) Load and clean the network
# -------------------------------------------------------------------------------

load_network <- function(file){
  nt <- read.csv(file, stringsAsFactors = FALSE)
  
  edges <- nt %>%
    dplyr::transmute(A = canonical_uniprot(interactor_A),
                     B = canonical_uniprot(interactor_B)) %>%
    dplyr::filter(grepl(uniprot_regex, A), grepl(uniprot_regex, B), A != B)
  
  n_raw <- nrow(nt)
  g <- graph_from_data_frame(edges, directed = FALSE)
  g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  
  comp <- components(g)
  cat(sprintf("Network: %d rows -> %d nodes, %d unique edges, %d components (largest = %d nodes)\n",
              n_raw, vcount(g), ecount(g), comp$no, max(comp$csize)))
  g
}

# -------------------------------------------------------------------------------
# 3) Seeds
# -------------------------------------------------------------------------------

prepare_seeds <- function(scores, ids, g, topN, use_abs = TRUE, prefix = ""){
  df <- data.frame(id = canonical_uniprot(ids), score = scores,
                   stringsAsFactors = FALSE) %>%
    dplyr::filter(!is.na(id), id != "", !is.na(score)) %>%
    dplyr::mutate(w = if (use_abs) abs(score) else score) %>%
    dplyr::filter(id %in% V(g)$name, w > 0) %>%
    dplyr::group_by(id) %>%
    dplyr::slice_max(w, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(w)) %>%
    dplyr::slice_head(n = topN)
  
  if (nrow(df) == 0) stop("No seed nodes found in graph for ", prefix)
  
  cat(sprintf("%s - seeds in graph: %d (positive score: %d, negative score: %d)\n",
              prefix, nrow(df), sum(df$score > 0), sum(df$score < 0)))
  
  list(weights = setNames(df$w, df$id),
       sign    = setNames(sign(df$score), df$id))
}

# -------------------------------------------------------------------------------
# 4) PPR and degree-matched permutation null
# -------------------------------------------------------------------------------

run_ppr <- function(g, weights, damping){
  w <- tapply(weights, names(weights), sum)        # merges duplicated nodes
  pers <- setNames(numeric(vcount(g)), V(g)$name)
  pers[names(w)] <- w
  pers <- pers / sum(pers)                         # sum normalisation, no seed lost
  page_rank(g, personalized = pers, damping = damping)$vector
}

make_degree_bins <- function(g, n_bins){
  deg <- degree(g)
  bin_of <- setNames(dplyr::ntile(deg, n_bins), V(g)$name)
  list(bin_of = bin_of, members = split(names(bin_of), bin_of))
}

# Each seed is replaced by a random node from the same degree bin and keeps its weight
sample_matched_seeds <- function(weights, bins){
  seed_bins <- bins$bin_of[names(weights)]
  new_names <- character(length(weights))
  for (b in unique(seed_bins)){
    idx  <- which(seed_bins == b)
    pool <- bins$members[[as.character(b)]]
    new_names[idx] <- if (length(pool) == 1) rep(pool, length(idx)) else
      sample(pool, length(idx), replace = length(pool) < length(idx))
  }
  setNames(weights, new_names)
}

ppr_with_null <- function(g, weights, bins, n_perm, damping){
  obs <- run_ppr(g, weights, damping)
  s1 <- s2 <- n_ge <- numeric(length(obs))
  for (i in seq_len(n_perm)){
    p  <- run_ppr(g, sample_matched_seeds(weights, bins), damping)
    s1 <- s1 + p
    s2 <- s2 + p^2
    n_ge <- n_ge + (p >= obs)
  }
  mu <- s1 / n_perm
  sdv <- sqrt(pmax(s2 / n_perm - mu^2, 0))
  z <- ifelse(sdv > 0, (obs - mu) / sdv, NA_real_)
  data.frame(node = names(obs), ppr = obs, null_mean = mu, null_sd = sdv,
             z = z, emp_p = (n_ge + 1) / (n_perm + 1), degree = degree(g)[names(obs)],
             stringsAsFactors = FALSE)
}

# -------------------------------------------------------------------------------
# 5) Cluster statistics
# -------------------------------------------------------------------------------

summarise_clusters <- function(cl, seeds, node_df, min_size, min_seed){
  z   <- setNames(node_df$z,   node_df$node)
  ppr <- setNames(node_df$ppr, node_df$node)
  
  cl_list <- split(names(cl), cl)
  nonseed_top <- setdiff(names(cl), seeds)
  
  tbl <- dplyr::bind_rows(lapply(names(cl_list), function(k){
    x  <- cl_list[[k]]
    ns <- setdiff(x, seeds)                  # non-seed nodes of the cluster
    bg <- setdiff(nonseed_top, ns)           # non-seed nodes of all OTHER clusters
    ks_p <- wx_p <- NA_real_
    if (length(x) >= min_size && length(ns) >= 3 && length(bg) >= 3){
      # "less" = CDF of cluster lies BELOW background = cluster values are LARGER
      ks_p <- suppressWarnings(ks.test(z[ns], z[bg], alternative = "less")$p.value)
      wx_p <- suppressWarnings(wilcox.test(z[ns], z[bg], alternative = "greater")$p.value)
    }
    data.frame(cluster = k,
               size = length(x),
               n_seed = sum(x %in% seeds),
               n_nonseed = length(ns),
               mean_ppr_all = mean(ppr[x]),
               median_ppr_nonseed = if (length(ns)) median(ppr[ns]) else NA_real_,
               median_z_nonseed   = if (length(ns)) median(z[ns], na.rm = TRUE) else NA_real_,
               median_z_background = if (length(bg)) median(z[bg], na.rm = TRUE) else NA_real_,
               ks_p = ks_p, wilcox_p = wx_p,
               stringsAsFactors = FALSE)
  }))
  
  # filter first, then correct only the clusters that will be interpreted
  tbl <- tbl %>% dplyr::filter(size >= min_size, n_seed >= min_seed)
  tbl$ks_fdr     <- p.adjust(tbl$ks_p, method = "BH")
  tbl$wilcox_fdr <- p.adjust(tbl$wilcox_p, method = "BH")
  dplyr::arrange(tbl, wilcox_fdr)
}

# -------------------------------------------------------------------------------
# 6) Main function: PPR + null + top subgraph + clustering
# -------------------------------------------------------------------------------

run_ras_ppr <- function(g, seed_scores, seed_ids, bins, p, prefix){
  
  sd_obj <- prepare_seeds(seed_scores, seed_ids, g, p$topN, p$seed_use_abs, prefix)
  seed_nodes <- names(sd_obj$weights)
  
  node_df <- ppr_with_null(g, sd_obj$weights, bins, p$n_perm, p$damping)
  node_df$is_seed <- node_df$node %in% seed_nodes
  
  score <- if (p$score_type == "z") node_df$z else node_df$ppr
  ok  <- is.finite(score) & node_df$ppr > 0
  thr <- quantile(score[ok], p$quant_thr)
  top_nodes <- node_df$node[ok & score >= thr]
  
  g_top <- induced_subgraph(g, vids = top_nodes)
  cl <- membership(cluster_walktrap(g_top))
  
  cat(sprintf("%s - top subgraph: %d nodes, %d edges, %d seeds inside, %d clusters\n",
              prefix, vcount(g_top), ecount(g_top),
              sum(seed_nodes %in% top_nodes), length(unique(cl))))
  
  summary_tbl <- summarise_clusters(cl, seed_nodes, node_df,
                                    p$min_cluster_size, p$min_cluster_seed)
  
  list(node_scores = dplyr::arrange(node_df, dplyr::desc(z)),
       g_top = g_top,
       clusters = cl,
       summary = summary_tbl,
       seeds = seed_nodes,
       seed_sign = sd_obj$sign)
}

# -------------------------------------------------------------------------------
# 7) Network exploration: seed neighbours, bridges between seeds, top nodes
# -------------------------------------------------------------------------------

find_seed_bridges <- function(g, seeds){
  seeds <- intersect(seeds, V(g)$name)
  n <- length(seeds)
  if (n < 2) return(data.frame(node = character(0), n_seed_pairs = integer(0)))
  
  inner <- character(0)
  for (i in seq_len(n - 1)){
    paths <- suppressWarnings(
      shortest_paths(g, from = seeds[i], to = seeds[(i + 1):n], output = "vpath")$vpath)
    for (pth in paths){
      v <- as_ids(pth)
      if (length(v) > 2) inner <- c(inner, v[-c(1, length(v))])
    }
  }
  inner <- inner[!inner %in% seeds]
  if (!length(inner)) return(data.frame(node = character(0), n_seed_pairs = integer(0)))
  tab <- table(inner)
  data.frame(node = names(tab), n_seed_pairs = as.integer(tab), stringsAsFactors = FALSE)
}

explore_network <- function(g, res, min_pairs = 2, top_quantile = 0.95){
  seeds <- res$seeds
  ns <- res$node_scores
  
  # first neighbours of seeds (previously called "bridge_direct")
  neigh <- unique(unlist(lapply(neighborhood(g, order = 1, nodes = seeds), as_ids)))
  seed_neigh <- ns %>% dplyr::filter(node %in% setdiff(neigh, seeds))
  
  # nodes lying on shortest paths between >= min_pairs seed pairs.
  # NB: shortest paths favour hubs -> always read together with degree and z.
  bridges <- find_seed_bridges(g, seeds) %>%
    dplyr::filter(n_seed_pairs >= min_pairs) %>%
    dplyr::left_join(ns, by = "node") %>%
    dplyr::arrange(dplyr::desc(n_seed_pairs))
  
  z_thr <- quantile(ns$z, top_quantile, na.rm = TRUE)
  top_nonseed <- ns %>% dplyr::filter(!is_seed, !is.na(z), z >= z_thr)
  
  list(seed_neighbours = seed_neigh, bridges = bridges, top_nonseed = top_nonseed)
}

# -------------------------------------------------------------------------------
# 8) Cross-paralog overlap of clusters
# -------------------------------------------------------------------------------

cluster_jaccard <- function(resA, resB){
  getcl <- function(r){
    l <- split(names(r$clusters), r$clusters)
    l[as.character(r$summary$cluster)]
  }
  A <- getcl(resA); B <- getcl(resB)
  if (!length(A) || !length(B)) return(matrix(numeric(0), 0, 0))
  sapply(B, function(b) sapply(A, function(a)
    length(intersect(a, b)) / length(union(a, b))))
}

# -------------------------------------------------------------------------------
# 9) Run workflow
# -------------------------------------------------------------------------------

if (sys.nframe() == 0L || interactive()) {
  
  seeds_tbl   <- read.csv(params$seeds_file, stringsAsFactors = FALSE)
  interactome <- load_network(params$network_file)
  bins        <- make_degree_bins(interactome, params$n_degree_bins)
  
  paralog_cols <- c(KRAS = "seed_kras", NRAS = "seed_nras", HRAS = "seed_hras")
  
  results <- lapply(names(paralog_cols), function(pr)
    run_ras_ppr(interactome, seeds_tbl[[paralog_cols[[pr]]]],
                seeds_tbl$uniprotswissprot, bins, params, pr))
  names(results) <- names(paralog_cols)
  
  explore <- lapply(results, function(r) explore_network(interactome, r))
  
  jaccard <- list(KRAS_NRAS = cluster_jaccard(results$KRAS, results$NRAS),
                  KRAS_HRAS = cluster_jaccard(results$KRAS, results$HRAS),
                  NRAS_HRAS = cluster_jaccard(results$NRAS, results$HRAS))
  
  for (pr in names(results)){
    cat("\n====", pr, "====\n")
    print(results[[pr]]$summary)
    lc <- tolower(pr)
    saveRDS(results[[pr]], file.path(params$out_dir, paste0("result_", lc, ".rds")))
    saveRDS(explore[[pr]], file.path(params$out_dir, paste0("explore_", lc, ".rds")))
    write.csv(results[[pr]]$summary,
              file.path(params$out_dir, paste0("cluster_summary_", lc, ".csv")), row.names = FALSE)
    write.csv(results[[pr]]$node_scores,
              file.path(params$out_dir, paste0("node_scores_", lc, ".csv")), row.names = FALSE)
    cl <- results[[pr]]$clusters
    for (k in results[[pr]]$summary$cluster){
      writeLines(names(cl[cl == as.integer(k)]),
                 file.path(params$out_dir, sprintf("%s_cluster%s_proteins.txt", lc, k)))
    }
  }
  
  cat("\n==== Jaccard overlap between paralog clusters ====\n")
  print(lapply(jaccard, round, 2))
  saveRDS(jaccard, file.path(params$out_dir, "cluster_jaccard.rds"))
}


mod_info <- function(res, k){
  nodes <- names(res$clusters)[res$clusters == as.integer(k)]
  ns <- res$node_scores[res$node_scores$node %in% nodes, ]
  list(seeds   = intersect(nodes, res$seeds),
       n_sig   = sum(ns$emp_p < 0.05 & !ns$is_seed),
       top     = head(ns[!ns$is_seed, c("node", "z", "emp_p", "degree")], 15))
}
mod_info(results$HRAS, 12); mod_info(results$KRAS, 43)   # modulo C
mod_info(results$KRAS, 33); mod_info(results$HRAS, 3)    # moduli D ed E

#----------------- check for 10k permutations -----------------------

if (sys.nframe() == 0L || interactive()) {
  
  params$n_perm <- 10000
  results_10Kperm <- lapply(names(paralog_cols), function(pr)
    run_ras_ppr(interactome, seeds_tbl[[paralog_cols[[pr]]]],
                seeds_tbl$uniprotswissprot, bins, params, pr))
  names(results_10Kperm) <- names(paralog_cols)
  
  explore <- lapply(results_10Kperm, function(r) explore_network(interactome, r))
  
  jaccard <- list(KRAS_NRAS = cluster_jaccard(results_10Kperm$KRAS, results_10Kperm$NRAS),
                  KRAS_HRAS = cluster_jaccard(results_10Kperm$KRAS, results_10Kperm$HRAS),
                  NRAS_HRAS = cluster_jaccard(results_10Kperm$NRAS, results_10Kperm$HRAS))
  
  for (pr in names(results_10Kperm)){
    cat("\n====", pr, "10k perm ====\n")
    print(results_10Kperm[[pr]]$summary)
    lc <- tolower(pr)
    saveRDS(results_10Kperm[[pr]], file.path(params$out_dir, paste0("result_", lc, ".rds")))
    saveRDS(explore[[pr]], file.path(params$out_dir, paste0("explore_", lc, ".rds")))
    write.csv(results_10Kperm[[pr]]$summary,
              file.path(params$out_dir, paste0("cluster_summary_", lc, ".csv")), row.names = FALSE)
    write.csv(results_10Kperm[[pr]]$node_scores,
              file.path(params$out_dir, paste0("node_scores_", lc, ".csv")), row.names = FALSE)
    cl <- results_10Kperm[[pr]]$clusters
    for (k in results_10Kperm[[pr]]$summary$cluster){
      writeLines(names(cl[cl == as.integer(k)]),
                 file.path(params$out_dir, sprintf("%s_cluster%s_proteins.txt", lc, k)))
    }
  }
  
  cat("\n==== Jaccard overlap between paralog clusters, 10k perm ====\n")
  print(lapply(jaccard, round, 2))
  saveRDS(jaccard, file.path(params$out_dir, "cluster_jaccard_10kperm.rds"))
}


#----------------- check for different seed -----------------------
set.seed(5678)

if (sys.nframe() == 0L || interactive()) {
  
  results_newseed <- lapply(names(paralog_cols), function(pr)
    run_ras_ppr(interactome, seeds_tbl[[paralog_cols[[pr]]]],
                seeds_tbl$uniprotswissprot, bins, params, pr))
  names(results_newseed) <- names(paralog_cols)
  
  explore <- lapply(results_newseed, function(r) explore_network(interactome, r))
  
  jaccard <- list(KRAS_NRAS = cluster_jaccard(results_newseed$KRAS, results_newseed$NRAS),
                  KRAS_HRAS = cluster_jaccard(results_newseed$KRAS, results_newseed$HRAS),
                  NRAS_HRAS = cluster_jaccard(results_newseed$NRAS, results_newseed$HRAS))
  
  for (pr in names(results_newseed)){
    cat("\n====", pr, "10k perm ====\n")
    print(results_newseed[[pr]]$summary)
    lc <- tolower(pr)
    saveRDS(results_newseed[[pr]], file.path(params$out_dir, paste0("result_", lc, ".rds")))
    saveRDS(explore[[pr]], file.path(params$out_dir, paste0("explore_", lc, ".rds")))
    write.csv(results_newseed[[pr]]$summary,
              file.path(params$out_dir, paste0("cluster_summary_", lc, ".csv")), row.names = FALSE)
    write.csv(results_newseed[[pr]]$node_scores,
              file.path(params$out_dir, paste0("node_scores_", lc, ".csv")), row.names = FALSE)
    cl <- results_newseed[[pr]]$clusters
    for (k in results_newseed[[pr]]$summary$cluster){
      writeLines(names(cl[cl == as.integer(k)]),
                 file.path(params$out_dir, sprintf("%s_cluster%s_proteins.txt", lc, k)))
    }
  }
  
  cat("\n==== Jaccard overlap between paralog clusters, 10k perm ====\n")
  print(lapply(jaccard, round, 2))
  saveRDS(jaccard, file.path(params$out_dir, "cluster_jaccard_newseed.rds"))
}


#--------------------- Check if the module z-score between the paralogs -----------------------------------
mod_nodes <- function(res, k) names(res$clusters)[res$clusters == as.integer(k)]
eval_fixed_module <- function(res, nodes){
  z  <- setNames(res$node_scores$z, res$node_scores$node)
  m  <- setdiff(nodes, res$seeds)
  bg <- setdiff(res$node_scores$node, c(nodes, res$seeds))
  c(n_nonseed = length(m), median_z = median(z[m], na.rm = TRUE),
    wilcox_p = wilcox.test(z[m], z[bg], alternative = "greater")$p.value)
}
modA <- mod_nodes(results_10Kperm$KRAS, 18); eval_fixed_module(results_newseed$HRAS, modA)

mods <- list(A = mod_nodes(results_10Kperm$KRAS, 18),
             C = mod_nodes(results_10Kperm$HRAS, 11),
             D = mod_nodes(results_10Kperm$KRAS, 32))
runs <- list(perm10k = results_10Kperm, seed5678 = results_newseed)
tab <- do.call(rbind, lapply(names(runs), function(r) do.call(rbind, lapply(c("KRAS","NRAS","HRAS"), function(p)
  do.call(rbind, lapply(names(mods), function(m)
    data.frame(run = r, paralog = p, module = m, t(eval_fixed_module(runs[[r]][[p]], mods[[m]])))))))))
tab

for (p in c("NRAS", "HRAS")) {
  cl <- results_10Kperm[[p]]$clusters
  cat(p, ": nel top 15% ", sum(mods$D %in% names(cl)), " nodi su ", length(mods$D), "\n")
  print(table(cl[intersect(mods$D, names(cl))]))
}


seedsB <- intersect(mods$C, results_10Kperm$HRAS$seeds)   # mods$C nel codice = nuovo modulo B
genesB <- seeds_tbl$hgnc_symbol[match(seedsB, seeds_tbl$uniprotswissprot)]
wt  <- !(colData$hras_mut | colData$kras_mut | colData$nras_mut)
mut <- which(colData$hras_mut)
zmat <- sapply(genesB, function(g) sapply(mut, function(i) {
  ref <- expr[g, wt & colData$CANCER_TYPE == colData$CANCER_TYPE[i]]
  (expr[g, i] - mean(ref)) / sd(ref) }))
rownames(zmat) <- paste(colnames(expr)[mut], colData$CANCER_TYPE[mut]); round(zmat, 2)
