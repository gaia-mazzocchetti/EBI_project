
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


