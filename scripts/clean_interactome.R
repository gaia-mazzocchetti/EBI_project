library(tidyverse)
cancer_inter <- read.delim("../../Documents/projects/EBI/Cancer.txt")

head(cancer_inter)
table(cancer_inter$Host.organism.s.) %>% sort(decreasing = T) %>% head()

cancer_inter$Host.organism.s.[grepl(pattern = "human", cancer_inter$Host.organism.s.)]

table(cancer_inter$Interaction.type.s.) %>% sort(decreasing = T) %>% head()

cancer_inter_cln <- cancer_inter %>% filter(grepl("human", Host.organism.s.),
                        grepl("association|interaction", Interaction.type.s.)) %>% 
  mutate(MI_score = as.numeric(str_remove(string = Confidence.value.s., pattern = "intact-miscore:"))) %>% 
  filter(MI_score >= 0.6)

table(cancer_inter_cln$Interaction.type.s.) 

cancer_ready <- cancer_inter %>% select(X.ID.s..interactor.A, ID.s..interactor.B) %>% 
  mutate(interactor_A = gsub(pattern ="uniprotkb:" , x = X.ID.s..interactor.A, replacement = ""),
         interactor_B = gsub(pattern = "uniprotkb:", x = ID.s..interactor.B, replacement = "")) %>% 
  select(interactor_A, interactor_B)

write.csv(cancer_ready, "input/all_cancer_interactome.csv", quote = F, row.names = F)
