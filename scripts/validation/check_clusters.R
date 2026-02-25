library(igraph)
g_kras   # igraph oggetto propagato KRAS
g_nras   # igraph oggetto propagato NRAS
g_hras   # igraph oggetto propagato HRAS

cluster_kras <- cluster_louvain(g_kras)
cluster_nras <- cluster_louvain(g_nras)
cluster_hras <- cluster_louvain(g_hras)

top100_kras <- seed_kras_100  # vettore nomi nodi seed
top100_nras <- seed_nras_100
top100_hras <- seed_hras_100

cluster2_nodes <- names(result_kras[["clusters"]][result_kras[["clusters"]] == 2])  # vettore nodi del modulo KRAS (Cluster 2)
backbone_nodes <- shared_all # 6 nodi condivisi o lista shared core


# observed modularity
Q_kras <- modularity(cluster_kras)
Q_nras <- modularity(cluster_nras)
Q_hras <- modularity(cluster_hras)

df_mod <- data.frame(
  paralog = c("KRAS","NRAS","HRAS"),
  Q = c(Q_kras, Q_nras, Q_hras)
)


vcount(g_kras)
vcount(g_nras)
vcount(g_hras)

ecount(g_kras)
ecount(g_nras)
ecount(g_hras)

length(intersect(V(g_nras)$name, V(g_hras)$name)) / vcount(g_nras)

print(cluster_louvain(g_kras))
sizes(cluster_louvain(g_hras))
sizes(cluster_louvain(g_nras))


mem <- membership(cluster_kras)
cluster2_nodes_lov <- names(mem[mem == 2])
cluster1_nodes_lov <- names(mem[mem == 1])
cluster3_nodes_lov <- names(mem[mem == 3]) 

length(intersect(cluster2_nodes, cluster2_nodes_lov)) /
  length(cluster2_nodes)

length(intersect(cluster2_nodes, cluster1_nodes_lov)) /
  length(cluster2_nodes)
length(intersect(cluster2_nodes, cluster3_nodes_lov)) /
  length(cluster2_nodes)

sub_walk <- induced_subgraph(g_kras, cluster2_nodes)
edge_density(sub_walk)
edge_density(g_kras)

mean(node_df$bet[node_df$node %in% cluster2_nodes])
mean(node_df$bet[node_df$module == "Backbone"])

intersect(cluster2_nodes, c("P04637","P42345","P78527","O14980","P21796","P45880"))

edge_density(g_nras)
edge_density(g_hras)

# permutation test
nperm <- 1000
Q_perm <- numeric(nperm)

for(i in 1:nperm){
  
  rand_nodes <- sample(V(g_global)$name, length(top100_kras))
  sub_rand <- induced_subgraph(g_global, vids = rand_nodes)
  
  cl_rand <- cluster_louvain(sub_rand)
  Q_perm[i] <- modularity(cl_rand)
}


# centrality
bet_kras <- betweenness(g_kras, normalized = TRUE)
eig_kras <- eigen_centrality(g_kras)$vector
close_kras <- closeness(g_kras, normalized = TRUE)

node_df <- data.frame(
  node = names(bet_kras),
  bet = bet_kras,
  eig = eig_kras,
  close = close_kras,
  module = ifelse(names(bet_kras) %in% cluster2_nodes, "Module2",
                  ifelse(names(bet_kras) %in% backbone_nodes, "Backbone","Other"))
)

wilcox.test(bet ~ module, 
            data = subset(node_df, module %in% c("Module2","Backbone")))

