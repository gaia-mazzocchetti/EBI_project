library(tidyverse)
 
deg_kras <- read.csv("output/deg_kras.csv")
deg_hras <- read.csv("output/deg_hras.csv")
deg_nras <- read.csv("output/deg_nras.csv")

intersect(deg_kras$X, deg_hras$X)
intersect(deg_kras$X, deg_nras$X)
intersect(deg_hras$X, deg_nras$X)

writeLines(unique(deg_kras$X))
writeLines(unique(deg_hras$X))
writeLines(unique(deg_nras$X))



seed_hras_100 <- result_hras$seeds
seed_kras_100 <- result_kras$seeds
seed_nras_100 <- result_nras$seeds

intersect(seed_kras_100, seed_hras_100) |> length()
intersect(seed_kras_100, seed_nras_100) |> length() 
intersect(seed_nras_100, seed_hras_100) |> length()

intersect(seed_kras_100, seed_hras_100) |> writeLines()
intersect(seed_kras_100, seed_nras_100) |> writeLines()
intersect(seed_nras_100, seed_hras_100) |> writeLines()

shared_all <- Reduce(intersect, list(seed_kras_100,
                                     seed_nras_100,
                                     seed_hras_100))
