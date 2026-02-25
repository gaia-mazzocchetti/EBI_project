
mutation_data <- read.delim("data_mutations.txt", skip = 2)

ras_paralog_mut <- mutation_data %>% filter(Hugo_Symbol %in% c("KRAS", "NRAS", "HRAS"),
                                            DNA_VAF >= 0.10)

sample_info <- read.csv("../../paper_repository/input/colData.csv")
sample_info$Tumor_Sample_Barcode <- sample_info$SAMPLE_ID
ras_paralog_mut2 <- left_join(ras_paralog_mut, sample_info[,c(1,4,11)], by = c("Tumor_Sample_Barcode" = "SAMPLE_ID"))
ras_paralog_mut2 <- ras_paralog_mut2 %>% filter(!is.na(CANCER_TYPE))

ras_maf <- maftools::read.maf(maf = ras_paralog_mut2, clinicalData = sample_info)

library(pals)
library(maftools)

tumor_types <- unique(ras_maf@clinical.data$CANCER_TYPE)

cols <- setNames(
  pals::cols25(length(tumor_types)),
  tumor_types
)

pdf("../../plot/oncoplot.pdf", width = 12, height = 10)
maftools::oncoplot(
  maf = ras_maf,
  genes = c("KRAS","NRAS","HRAS"),
  clinicalFeatures = "CANCER_TYPE",
  annotationColor = list(CANCER_TYPE = cols),
  sortByAnnotation = TRUE,
  draw_titv = TRUE,
  drawRowBar = TRUE,
  drawColBar = TRUE,
  showTumorSampleBarcodes = FALSE
)

dev.off()


maftools::lollipopPlot(
  maf = ras_maf,
  gene = "KRAS",
  AACol = "Protein_position",
  showMutationRate = TRUE
)

maftools::lollipopPlot(
  maf = ras_maf,
  gene = "NRAS",
  AACol = "Protein_position"
)

maftools::lollipopPlot(
  maf = ras_maf,
  gene = "HRAS",
  AACol = "Protein_position"
)

pdf("../../plot/lollipop_combined.pdf", width = 10, height = 5)


# ---- KRAS ----
lollipopPlot(
  maf = ras_maf,
  gene = "KRAS",
  AACol = "Protein_position",
)

# ---- NRAS ----
lollipopPlot(
  maf = ras_maf,
  gene = "NRAS",
  AACol = "Protein_position"
)

# ---- HRAS ----
lollipopPlot(
  maf = ras_maf,
  gene = "HRAS",
  AACol = "Protein_position"
)

dev.off()

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


# Plot del network
set.seed(123)  # per layout riproducibile
ggraph(g_kras, layout = "fr") + 
  geom_edge_link(alpha = 0.3, color = "gray") +
  geom_node_point(aes(color = type, size = ppr)) +
  # Aggiungi i nomi dei nodi solo se PPR alto (hub)
  geom_node_text(aes(label = ifelse(is_hub, name, NA)),
                 repel = TRUE,   # evita sovrapposizioni
                 size = 3,
                 color = "black") +
  scale_color_manual(values = c("seed" = "red", "bridge" = "orange", "other" = "lightgray")) +
  scale_size_continuous(range = c(3, 8)) +
  theme_void() +
  labs(title = "KRAS Mutation-enriched Network",
       color = "Node type",
       size = "PageRank") +
  theme(legend.position = "right"

