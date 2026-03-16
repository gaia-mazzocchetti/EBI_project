

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
    title = "Pathways of KRAS-mutated samples",
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
    title = "Pathways of HRAS-mutated samples",
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
    title = "Pathways of NRAS-mutated samples",
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

 ggpubr::ggarrange(kras_plot, hras_plot, nras_plot, nrow = 1)
print(panel_ptw)
ggsave("../../plot/Figure_pathway_enrichment.pdf",
       width = 30, height = 10)



deg_kras <- read.csv("output/DEG/deg_kras.csv")
deg_hras <- read.csv("output/DEG/deg_hras.csv")
deg_nras <- read.csv("output/DEG/deg_nras.csv")

plot_db <- data.frame(Paralog = factor(c("KRAS", "HRAS", "NRAS"), levels = c("KRAS","NRAS", "HRAS")),
                      degs = c(nrow(deg_kras), nrow(deg_hras), nrow(deg_nras)))
library(ggplot2)
library(ggrepel)

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


volcano_plot <- function(df, title, highlight_genes){
  
  df <- df %>%
    mutate(
      significance = case_when(
        adj.P.Val < 0.05 & logFC > 0 ~ "Up",
        adj.P.Val < 0.05 & logFC < 0 ~ "Down",
        TRUE ~ "NS"
      ),
      label = ifelse(X %in% highlight_genes, X, NA)
    )
  
  ggplot(df, aes(x = logFC, y = -log10(adj.P.Val))) +
    
    geom_point(aes(color = significance), size = 2, alpha = 0.8) +
    
    scale_color_manual(values = c(
      "Up" = "#2E8B57",   # verde
      "Down" = "#7A3E9D", # viola
      "NS" = "grey80"
    )) +
    
    geom_text_repel(
      aes(label = label),
      size = 5,
      max.overlaps = 20
    ) +
    
    labs(
      title = title,
      x = "log2 Fold Change",
      y = "-log10 adjusted p-value"
    ) +
    
    theme_classic() +
    geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", colour="grey60") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey60")+
    
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 12)
    )
}

kras_genes <- c(
  "MUC5AC","MUC2","MUC5B","TFF1", "REG4", "DKK4", "KRT7",
  "DUSP4", "DUSP6","LEFTY1","RP11-350J20.12", "SLC26A3"
)

nras_genes <- c(
  "DEFA4","AZU1","LY6G6D","AC104135.4",
  "RP11-386M24.8","RP11-366L20.2", "SLC26A3"
)

hras_genes <- c(
  "DPPA5","IFITM5","KIR3DX1","LIPK",
  "AC005725.1","RP5-991C6.2","LINC00676", "RP11-59D5__B.2",
  "AC006262.10", "LINC00704", "RP11-259O2.1"
)


p_kras <- volcano_plot(deg_kras, "KRAS DEGs", kras_genes)
p_nras <- volcano_plot(deg_nras, "NRAS DEGs", nras_genes)
p_hras <- volcano_plot(deg_hras, "HRAS DEGs", hras_genes)

library(patchwork)

final_plot <- (p_kras | p_nras | p_hras) / (kras_plot |hras_plot| nras_plot) +
  plot_layout(heights = c(1, 1))


final_plot + plot_annotation(tag_levels = "A")

