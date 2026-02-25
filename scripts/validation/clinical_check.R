clinical <- fread("../data/pancan_pcawg_2020/data_clinical_patient.txt", skip = 4)
mut_info <- fread("input/colData.csv")

complete_db <- left_join(mut_info, clinical[,1:7])

clinical$resistance_status <- NA

complete_db$resistance_status[
  complete_db$`FIRST THERAPY_RESPONSE` %in% 
    c("Disease Progression", "Stable Disease")
] <- "Resistant"

complete_db$resistance_status[
  complete_db$`FIRST THERAPY_RESPONSE` %in% 
    c("Complete Response", "Partial Response")
] <- "Sensitive"

table(complete_db$resistance_status, useNA = "ifany")

expr <- fread("input/counts_log2fpkm.csv")

expr_scaled <- t(scale(t(expr[,-1])))
rownames(expr_scaled) <- expr$gene_name

cluster2_uniprot <- fread("output/results/kras_cluster2_proteins.txt") %>% unlist() # tutta la tua lista

library(org.Hs.eg.db)
library(AnnotationDbi)

cluster2_genes_hugo <- mapIds(org.Hs.eg.db,
                              keys = cluster2_uniprot,
                              column = "SYMBOL",
                              keytype = "UNIPROT",
                              multiVals = "first")

cluster2_genes <- intersect(cluster2_genes_hugo, rownames(expr_scaled))
cluster2_score <- colMeans(expr_scaled[cluster2_genes, ], na.rm = TRUE)

complete_db$cluster2_score <- cluster2_score[match(complete_db$SAMPLE_ID,
                                                colnames(expr_scaled))]
summary(complete_db$cluster2_score)

complete_db_mut <- complete_db[complete_db$group == "MT", ]
table(complete_db_mut$kras_mut, complete_db_mut$OS_STATUS)

kras_patients <- complete_db[complete_db$kras_mut == T &
                            !is.na(complete_db$resistance_status), ]
wilcox.test(cluster2_score ~ OS_STATUS,
            data = complete_db_mut)

library(ggplot2)

ggplot(complete_db_mut,
       aes(x = OS_STATUS,
           y = cluster2_score)) +
  geom_boxplot() +
  theme_bw()


complete_db_kmut <- complete_db[complete_db$kras_mut == T, ]

wilcox.test(cluster2_score ~ OS_STATUS,
            data = complete_db_kmut)

library(ggplot2)

ggplot(complete_db_kmut,
       aes(x = OS_STATUS,
           y = cluster2_score)) +
  geom_boxplot() +
  theme_bw()


mdl <- glm(as.factor(OS_STATUS) ~ cluster2_score,
           data = complete_db_mut,
           family = binomial)

summary(mdl)

mdl2 <- glm(as.factor(OS_STATUS) ~ cluster2_score + CANCER_TYPE + SEX + AGE,
    data = complete_db_mut,
    family = binomial)

summary(mdl2)

mdl3 <- glm(as.factor(OS_STATUS) ~ cluster2_score,
           data = complete_db_kmut,
           family = binomial)

summary(mdl3)

mdl4 <- glm(as.factor(OS_STATUS) ~ cluster2_score + CANCER_TYPE + SEX + AGE,
            data = complete_db_kmut,
            family = binomial)

summary(mdl4)

subset_colon <- subset(complete_db_mut,
                       CANCER_TYPE == "Colorectal Cancer")

cln_model <- glm(as.factor(OS_STATUS) ~ cluster2_score,
    data = subset_colon,
    family = binomial)
summary(cln_model)

complete_db_kmut$colon_cancer <- ifelse(complete_db_kmut$CANCER_TYPE == "Colorectal Cancer", 1, 0)

ggplot(complete_db_kmut,
       aes(x = as.factor(colon_cancer),
           y = cluster2_score)) +
  geom_boxplot() +
  theme_bw()


wilcox.test(cluster2_score ~ as.factor(colon_cancer),
            data = complete_db_kmut)
table(complete_db_kmut$colon_cancer)
