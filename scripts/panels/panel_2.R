

library(tidyverse)
library(forcats)

kras_ptw <- read.csv("output/DEG/kras_reactome.csv")
top_kras <- kras_ptw %>% arrange(Entities.FDR) %>% pull(Pathway.identifier) %>% head(20)

kras_ptw_top <- kras_ptw %>% filter(Pathway.identifier %in% top_kras)

hras_ptw <- read.csv("output/DEG/hras_reactome.csv")
top_hras <- hras_ptw %>% arrange(Entities.FDR) %>% pull(Pathway.identifier) %>% head(20)

hras_ptw_top <- hras_ptw %>% filter(Pathway.identifier %in% top_hras)

nras_ptw <- read.csv("output/DEG/nras_reactome.csv")
top_nras <- nras_ptw %>% arrange(Entities.FDR) %>% pull(Pathway.identifier) %>% head(20)

nras_ptw_top <- nras_ptw %>% filter(Pathway.identifier %in% top_nras)


embl_palette <- c("#5E3C99", "#B2ABD2", "#E6F5D0", "#1B7837")

kras_plot <- ggplot(
  kras_ptw_top %>%
    mutate(Pathway.name = fct_reorder(Pathway.name, -log(Entities.pValue))),
  aes(x = -log10(Entities.pValue),
      y = Pathway.name,
      fill = -log10(Entities.pValue))
) +
  geom_col(width = 0.8) +
  scale_fill_gradientn(colors = embl_palette) +
  labs(
    title = "Pathway enrichment in KRAS-mutated samples",
    x = expression(-log[10](p-value)),
    y = NULL,
    fill = expression(-log[10](p-value))
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 11),
    legend.position = "right"
  )

hras_plot <- ggplot(
  hras_ptw_top %>%
    mutate(Pathway.name = fct_reorder(Pathway.name, -log(Entities.pValue))),
  aes(x = -log10(Entities.pValue),
      y = Pathway.name,
      fill = -log10(Entities.pValue))
) +
  geom_col(width = 0.8) +
  scale_fill_gradientn(colors = embl_palette) +
  labs(
    title = "Pathway enrichment in HRAS-mutated samples",
    x = expression(-log[10](p-value)),
    y = NULL,
    fill = expression(-log[10](p-value))
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 11),
    legend.position = "right"
  )

nras_plot <- ggplot(
  nras_ptw_top %>%
    mutate(Pathway.name = fct_reorder(Pathway.name, -log(Entities.pValue))),
  aes(x = -log10(Entities.pValue),
      y = Pathway.name,
      fill = -log10(Entities.pValue))
) +
  geom_col(width = 0.8) +
  scale_fill_gradientn(colors = embl_palette) +
  labs(
    title = "Pathway enrichment in NRAS-mutated samples",
    x = expression(-log[10](p-value)),
    y = NULL,
    fill = expression(-log[10](p-value))
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 11),
    legend.position = "right"
  )

ggpubr::ggarrange(kras_plot, hras_plot, nras_plot, labels = c("A", "B", "C"))
ggsave("../../plot/Figure_pathway_enrichment.pdf",
       width = 20, height = 15)



setwd("~/Desktop/paper_analysis/paper_repository")
deg_kras <- read.csv("output/DEG/deg_kras.csv")
deg_hras <- read.csv("output/DEG/deg_hras.csv")
deg_nras <- read.csv("output/DEG/deg_nras.csv")

plot_db <- data.frame(Paralog = factor(c("KRAS", "HRAS", "NRAS"), levels = c("KRAS","NRAS", "HRAS")),
                      degs = c(nrow(deg_kras), nrow(deg_hras), nrow(deg_nras)))
library(ggplot2)

ggplot(plot_db, aes(Paralog, degs, fill= Paralog)) + 
  geom_bar(stat ="identity")+
  scale_fill_manual(values = c("#40bf40","#862d86","#9999ff"))+
  theme_bw()


deg_kras <- deg_kras %>% mutate(dir = ifelse(logFC >= 0.5, "up", ifelse(logFC <= -0.5, "down", "norm"))) %>% arrange(adj.P.Val, abs(logFC))
table(deg_kras$dir)
head(deg_kras, n=10)


deg_nras <- deg_nras %>% mutate(dir = ifelse(logFC >= 0.5, "up", ifelse(logFC <= -0.5, "down", "norm"))) %>% arrange(adj.P.Val, abs(logFC))
table(deg_nras$dir)
head(deg_nras, n=10)


deg_hras <- deg_hras %>% mutate(dir = ifelse(logFC >= 0.5, "up", ifelse(logFC <= -0.5, "down", "norm"))) %>% arrange(adj.P.Val, abs(logFC))
table(deg_hras$dir)
head(deg_hras, n=10)