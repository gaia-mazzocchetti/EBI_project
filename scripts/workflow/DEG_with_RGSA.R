library(ReactomeGSA)
library(Biobase)
library(dplyr)

colData <- read.csv("input/colData.csv")
counts <- data.table::fread("input/counts_log2fpkm.csv")

counts_collapsed <- counts %>%
  group_by(gene_name) %>%
  summarise(across(where(is.numeric), mean))

counts_collapsed <- as.data.frame(counts_collapsed)
rownames(counts_collapsed) <- counts_collapsed$gene_name
counts_collapsed <- counts_collapsed[,-1]

rownames(colData) <- colData$SAMPLE_ID
phenoData <- AnnotatedDataFrame(colData)

expr <- ExpressionSet(assayData = as.matrix(counts_collapsed), phenoData)


my_request <-ReactomeAnalysisRequest(method = "Camera")

# do not create a visualization for this example
my_request <- set_parameters(request = my_request, create_reactome_visualization = F)

# add the dataset using the loaded object
my_request <- add_dataset(request = my_request, 
                          expression_values = expr, 
                          name = "KRAS analysis", 
                          type = "microarray_norm",
                          comparison_factor = "kras_mut", 
                          comparison_group_1 = "TRUE", 
                          comparison_group_2 = "FALSE",
                          additional_factors = c("CANCER_TYPE", "PURITY"))

res <- perform_reactome_analysis(my_request)
