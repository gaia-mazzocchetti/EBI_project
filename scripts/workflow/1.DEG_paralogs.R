# ===============================================================
# Paralog-specific Differential Expression and Weighted Seed Generation
# for Network Propagation (log2(FPKM+1) RNA-seq)
# ===============================================================

suppressPackageStartupMessages({
  library(limma)   # For linear modeling on expression data
  library(edgeR)   # Provides voom function and DGE utilities
  library(dplyr) # For data manipulation
  library(ggplot2)
})

# ----------------------------
# 1) INPUT: FPKM and metadata
# ----------------------------
# fpkm: matrix (genes x samples), numeric FPKM values
# colData: data.frame with at least:
#   sample_id, tumor_type, batch (optional), purity (optional)
#   KRAS_mut, NRAS_mut, HRAS_mut (0/1 or TRUE/FALSE)


setwd("~/Desktop/paper_analysis/paper_repository")
colData <- read.csv("input/colData.csv")
counts <- data.table::fread("../input/counts_log2fpkm.csv")

rownames(colData) <- colData$SAMPLE_ID

# ----------------------------
# 2) Function for paralog-specific DE
# ----------------------------
run_limma_paralog <- function(paralog_col, counts, colData,
                              tumor_type_col = "CANCER_TYPE",
                              purity_col = "PURITY",
                              min_mut = 5,
                              pca = F) {
  
  stopifnot(paralog_col %in% colnames(colData))
  
  # ----------------------------
  # 1) Prepare numeric counts matrix
  # ----------------------------
  # Separate gene names from numeric matrix
  gene_names <- counts[[1]]             # first column contains gene symbols
  counts_numeric <- counts[, -1]        # remove first column
  counts_numeric <- as.data.frame(counts_numeric)
  
  # Make rownames unique (needed for limma and PCA)
  rownames(counts_numeric) <- make.unique(as.character(gene_names))
  
  # Ensure numeric matrix
  counts_numeric <- as.matrix(counts_numeric)
  storage.mode(counts_numeric) <- "numeric"
  
  # ----------------------------
  # 2) Define MUT vs WT-all
  # ----------------------------
  muts <- colData[[paralog_col]] == TRUE
  wts  <- colData$mparalogs == 0
  selected <- muts | wts
  
  counts_sel <- counts_numeric[, rownames(colData)[selected], drop=FALSE]
  colData_sel <- colData[selected, , drop=FALSE]
  
  # Status factor
  status <- ifelse(colData_sel[[paralog_col]], "MT", "WT")
  colData_sel$Status <- factor(status, levels = c("WT","MT"))
  
  n_mut <- sum(status=="MT")
  n_wt  <- sum(status=="WT")
  
  if(n_mut < min_mut){
    warning(sprintf("[%s] Too few MUT samples (%d). DE may be unstable.", paralog_col, n_mut))
  }
  if(length(unique(status)) < 2){
    stop(sprintf("[%s] Cannot run DE: need at least one MUT and one WT sample.", paralog_col))
  }
  
  # ----------------------------
  # 3) Covariates and design
  # ----------------------------
  colData_sel$TumorType <- factor(colData_sel[[tumor_type_col]])
  colData_sel$Purity    <- as.numeric(colData_sel[[purity_col]])
  
  design <- model.matrix(~ TumorType + Purity + Status, data = colData_sel)
  
  # ----------------------------
  # 4) limma DE
  # ----------------------------
  fit <- lmFit(counts_sel, design)
  fit <- eBayes(fit)
  
  # Find Status coefficient
  coef_name <- grep("^Status", colnames(fit$coefficients), value=TRUE)
  if(length(coef_name) == 0){
    stop(sprintf("[%s] Coefficient for MUT vs WT not found.", paralog_col))
  }
  
  res <- topTable(fit, coef = coef_name[1], number=Inf, sort.by="P")
  
  # ----------------------------
  # 5) Weighted seed
  # ----------------------------
  pval <- if("P.Value" %in% colnames(res)) res$P.Value else res$adj.P.Val
  seed <- with(res, sign(logFC) * -log10(pmax(pval, 1e-300)))
  names(seed) <- rownames(res)
  seed[is.na(seed)] <- 0
  
  # ----------------------------
  # 6) PCA diagnostic (robust)
  # ----------------------------
  if(pca) {
  gene_var <- apply(counts_sel, 1, var)
  counts_sel_filtered <- counts_sel[gene_var > 0, , drop=FALSE]
  
  if(ncol(counts_sel_filtered) >= 2 && nrow(counts_sel_filtered) >= 2 && n_mut >= 2 && n_wt >= 2){
    top_genes <- names(sort(gene_var[gene_var>0], decreasing=TRUE))[1:min(500, nrow(counts_sel_filtered))]
    pca_mat <- t(counts_sel_filtered[top_genes, ])
    pca_res <- prcomp(pca_mat, scale.=TRUE)
    
    pca_df <- data.frame(
      PC1 = pca_res$x[,1],
      PC2 = pca_res$x[,2],
      Status = status,
      TumorType = colData_sel$TumorType
    )
    
    p <- ggplot(pca_df, aes(x=PC1, y=PC2, color=TumorType, shape=Status)) +
      geom_point(size=3, alpha=0.8) +
      theme_bw() +
      labs(title=paste0("PCA: ", paralog_col, " MUT vs WT-all (top 500 variable genes)"),
           x="PC1", y="PC2") +
      theme(legend.position="right", legend.title=element_text(size=10))
    print(p)
    
  } else {
    message(sprintf("[%s] PCA skipped: too few non-zero variance genes or samples (MUT=%d, WT=%d).",
                    paralog_col, n_mut, n_wt))
  }
  }
  # ----------------------------
  # 7) Return results
  # ----------------------------
  return(list(res = res, seed = seed, n_MUT = n_mut, n_WT = n_wt))
}


# ----------------------------
# 3) Run DE for each paralog
# ----------------------------

kras <- run_limma_paralog("kras_mut", counts, colData)
nras <- run_limma_paralog("nras_mut", counts, colData)
hras <- run_limma_paralog("hras_mut", counts, colData)

# ----------------------------
# 4) Extract top DE genes (optional filter)
# ----------------------------
deg_kras <- kras$res %>% filter(!is.na(adj.P.Val) & adj.P.Val < 0.05 & abs(logFC) >= 0.5)
deg_nras <- nras$res %>% filter(!is.na(adj.P.Val) & adj.P.Val < 0.05 & abs(logFC) >= 0.5)
deg_hras <- hras$res %>% filter(!is.na(adj.P.Val) & adj.P.Val < 0.05 & abs(logFC) >= 0.5)

# Weighted seed vectors ready for network propagation
seed_kras <- kras$seed
seed_nras <- nras$seed
seed_hras <- hras$seed

seed_db <- data.frame(seed_kras= kras$seed,
                      seed_nras = nras$seed,
                      seed_hras= hras$seed)
seed_db <- rownames_to_column(seed_db)

write.csv(seed_db,"output/seed_weights.csv", row.names = T, quote = F)

write.csv(deg_kras,"output/deg_kras.csv", row.names = T, quote = F)
write.csv(deg_nras,"output/deg_nras.csv", row.names = T, quote = F)
write.csv(deg_hras,"output/deg_hras.csv", row.names = T, quote = F)

summary(seed_db)

plot(density(seed_db$seed_kras))
plot(density(seed_db$seed_nras))
plot(density(seed_db$seed_hras))

