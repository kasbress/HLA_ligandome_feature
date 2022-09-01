# title: "XGBoost for ligands trial"
# author: "BP Nicolet"
# date: "14/12/2021"

library(tidyverse)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(colorRamps)
library(ggfortify)
library(reshape)
library(circlize)
library(stringr)
library(biomaRt)
library(doMC)
library(caret)
library(randomForest)
library(e1071)
library(xgboost)
library(Matrix)
#install.packages("xgboost")

set.seed(12345)

setwd("/DATA/users/k.bresser/ligandome_project/analysis_new/Ligandome_models")




### importing Full_data ###
read_tsv("Data/Protein_per_Uniprot_entry_library_v2_RBP_GC_length_codon_AA_m6A_m5C_AtoI_m1A_m7G_CD8miRDB_PTM.csv") %>% 
  mutate(across(everything(), replace_na, 0)) %>% 
  inner_join(read_tsv("Output/test_train_tables_new/rf_train_peptides.tsv"),
             by = c("Entry" = "swissprot_id")) -> Full_data
dim(Full_data)


Full_data %>% 
  group_by(ligand) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(data = map2(data, c(40000,10000), sample_n)) %>% 
  unnest(data) -> Full_data


# Scale rna data
file.names <- list.files("./Data/rna_seq_data", pattern = "MEL") 

expression.data <- tibble(
  
  expression.info = map( paste0("./Data/rna_seq_data/", file.names), read_tsv ),
  tumor = map_chr(str_split(file.names, "_"), 1)
  
)

expression.data %>% 
  unnest(expression.info) %>% 
  group_by(tumor) %>% 
  mutate(rna = scale(log10(rna))[,1]) %>% 
  dplyr::select(swissprot_id, rna) %>% 
  nest() %>% 
  dplyr::rename(expr = data) -> expression.data

write_rds(expression.data, "./Data/rna_seq_data/all_lines_final.rds")
  
Full_data %>% 
  mutate(tumor = fct_recode(tumor, MEL1 = "M026", MEL2 = "mel95", MEL3 = "RTE")) %>% 
  dplyr::select(-rna) %>% 
  group_by(tumor) %>% 
  nest() %>% 
  inner_join(expression.data) %>% 
  transmute(tumor = tumor,
            data = map2(data, expr, inner_join, by = c("Entry" = "swissprot_id"))) %>% 
  unnest(data) %>% 
  ungroup() -> Full_data
 
write_tsv(Full_data, "./Data/Train_set_final.tsv")

dim(Full_data)

Full_data$swissprot_id <- NULL
Full_data$tumor <- NULL
Full_data$peptide <- NULL
Full_data$external_gene_name <- NULL
Full_data$ID <- NULL
Full_data$peptide_length <- NULL
Full_data$Entry <- NULL




###-------------------------------------------------------------------###
###-------------------------------------------------------------------###

### Modeling ###

###---------###
## Set up ##
###---------###
gc()
registerDoMC(1)
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 2,
                        search="grid",
                        verboseIter = TRUE,
                        sampling = "down",
                        allowParallel = T, 
                        returnData = F)

xgbGrid <- expand.grid(nrounds = 1000,
                       max_depth = 1,
                       colsample_bytree = 0.5,
                       eta = 0.3,
                       gamma=1,
                       min_child_weight = 0.9,
                       subsample = 1)

set.seed(12345)


#objective = "binary:logistic"


##-----------------------###
# Random model ##
##-----------------------###
print("Random model")

start_time <- Sys.time()


Full_data %>% 
  dplyr::mutate(ligand = sample(ligand)) %>% 
  train(as.factor(ligand)~. ,
        data=.,
        method="xgbTree",
        trControl=control,
        metric="Accuracy",
        tuneGrid= xgbGrid,
        na.action = na.omit,
        nthread = 12,
        verbose = TRUE) -> XGB.model

end_time <- Sys.time()
print(end_time - start_time) # Time difference of

print(XGB.model)

saveRDS(XGB.model,"./XGB_models_final/xgb_ligands|random|24_04_2022.RDS")

rm(XGB.model)
gc()






###-----------------------###
## RNA only model ##
###-----------------------###
print("RNA only model")

start_time <- Sys.time()


Full_data %>% 
  dplyr::select(c(ligand, rna)) %>% 
  train(as.factor(ligand)~. ,
        data=.,
        method="xgbTree",
        trControl=control,
        metric="Accuracy",
        tuneGrid= xgbGrid,
        na.action = na.omit,
        nthread = 12,
        verbose = TRUE) -> XGB.model

end_time <- Sys.time()
print(end_time - start_time) # Time difference of

print(XGB.model)

saveRDS(XGB.model,"./XGB_models_final/xgb_ligands|rna_only|24_04_2022.RDS")

rm(XGB.model)
gc()




###-----------------------###
## ribo only model ##
###-----------------------###
print("ribo only model")

start_time <- Sys.time()


Full_data %>% 
  dplyr::select(c(ligand, ribo)) %>% 
  train(as.factor(ligand)~. ,
        data=.,
        method="xgbTree",
        trControl=control,
        metric="Accuracy",
        tuneGrid= xgbGrid,
        na.action = na.omit,
        nthread = 12,
        verbose = TRUE) -> XGB.model

end_time <- Sys.time()
print(end_time - start_time) # Time difference of

print(XGB.model)

saveRDS(XGB.model,"./XGB_models_final/xgb_ligands|ribo_only|24_04_2022.RDS")

rm(XGB.model)
gc()




###-----------------------###
## rna + ribo only model ##
###-----------------------###
print("rna + ribo only model")

start_time <- Sys.time()


Full_data %>% 
  dplyr::select(c(ligand, ribo, rna)) %>% 
  train(as.factor(ligand)~. ,
        data=.,
        method="xgbTree",
        trControl=control,
        metric="Accuracy",
        tuneGrid= xgbGrid,
        na.action = na.omit,
        nthread = 12,
        verbose = TRUE) -> XGB.model

end_time <- Sys.time()
print(end_time - start_time) # Time difference of

print(XGB.model)

saveRDS(XGB.model,"./XGB_models_final/xgb_ligands|ribo_only|24_04_2022.RDS")

rm(XGB.model)
gc()


###-----------------------###
## Aff + chop only model ##
###-----------------------###
print("Aff + chop model")

start_time <- Sys.time()

Full_data %>% 
  dplyr::select(c(ligand, affMIN, chop)) %>% 
  train(as.factor(ligand)~. ,
        data=.,
        method="xgbTree",
        trControl=control,
        metric="Accuracy",
        tuneGrid= xgbGrid,
        na.action = na.omit,
        nthread = 12,
        verbose = TRUE) -> XGB.model

end_time <- Sys.time()
print(end_time - start_time) # Time difference of

print(XGB.model)

saveRDS(XGB.model,"./XGB_models_final/xgb_ligands|aff_chop_only|24_04_2022.RDS")



rm(XGB.model)
gc()



###-----------------------###
## Aff only model ##
###-----------------------###
print("Aff model")

start_time <- Sys.time()

Full_data %>% 
  dplyr::select(c(ligand, affMIN)) %>% 
  train(as.factor(ligand)~. ,
        data=.,
        method="xgbTree",
        trControl=control,
        metric="Accuracy",
        tuneGrid= xgbGrid,
        na.action = na.omit,
        nthread = 12,
        verbose = TRUE) -> XGB.model

end_time <- Sys.time()
print(end_time - start_time) # Time difference of

print(XGB.model)

saveRDS(XGB.model,"./XGB_models_final/xgb_ligands|aff_only|24_04_2022.RDS")



rm(XGB.model)
gc()





###-----------------------###
## Library only model ##
###-----------------------###
print("Library only model")

start_time <- Sys.time()

Full_data %>% 
  dplyr::select(-c(affMIN, chop, rna, ribo)) %>% 
  train(as.factor(ligand)~. ,
        data=.,
        method="xgbTree",
        trControl=control,
        metric="Accuracy",
        tuneGrid= xgbGrid,
        na.action = na.omit,
        nthread = 12,
        verbose = TRUE) -> XGB.model

end_time <- Sys.time()
print(end_time - start_time) # Time difference of

print(XGB.model)

saveRDS(XGB.model,"./XGB_models_final/xgb_ligands|library_only|23_12_21.RDS")



rm(XGB.model)
gc()


###---------------------------------------------------###
## Modeling with only RNA +  chop + affMIN ## 
###---------------------------------------------------###


print("Modeling with only RNA +  chop + affMIN ")

start_time <- Sys.time()

Full_data %>% 
  dplyr::select(c(ligand, affMIN, chop, rna)) %>% 
  train(as.factor(ligand)~. ,
        data=.,
        method="xgbTree",
        trControl=control,
        metric="Accuracy",
        tuneGrid= xgbGrid,
        na.action = na.omit,
        nthread = 12,
        verbose = TRUE) -> XGB.model

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 31.67485 mins

print(XGB.model)

saveRDS(XGB.model,"./XGB_models_final/xgb_ligands|aff_chop_rna|24_04_2022.RDS")



rm(XGB.model)
gc()



###---------------------------------------------------###
## Modeling with only ribo +  chop + affMIN ## 
###---------------------------------------------------###


print("Modeling with only ribo +  chop + affMIN ")

start_time <- Sys.time()

Full_data %>% 
  dplyr::select(c(ligand, affMIN, chop, ribo)) %>% 
  train(as.factor(ligand)~. ,
        data=.,
        method="xgbTree",
        trControl=control,
        metric="Accuracy",
        tuneGrid= xgbGrid,
        na.action = na.omit,
        nthread = 12,
        verbose = TRUE) -> XGB.model

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 31.67485 mins

print(XGB.model)

saveRDS(XGB.model,"./XGB_models_final/xgb_ligands|aff_chop_ribo|24_04_2022.RDS")



rm(XGB.model)
gc()


###---------------------------------------------------###
## Modeling with only rna + ribo +  chop + affMIN ## 
###---------------------------------------------------###


print("Modeling with only rna + ribo +  chop + affMIN ")

start_time <- Sys.time()

Full_data %>% 
  dplyr::select(c(ligand, affMIN, chop, ribo, rna)) %>% 
  train(as.factor(ligand)~. ,
        data=.,
        method="xgbTree",
        trControl=control,
        metric="Accuracy",
        tuneGrid= xgbGrid,
        na.action = na.omit,
        nthread = 12,
        verbose = TRUE) -> XGB.model

end_time <- Sys.time()
print(end_time - start_time) # Time difference of 31.67485 mins

print(XGB.model)

saveRDS(XGB.model,"./XGB_models_final/xgb_ligands|aff_chop_rna_ribo|24_04_2022.RDS")



rm(XGB.model)
gc()



###---------------------------------------###
## Modeling with library + chop + affMIN ## 
###---------------------------------------###

print("Modeling with library + chop + affMIN")


start_time <- Sys.time()

Full_data %>% 
  dplyr::select(-c(rna, ribo)) %>% 
  train(as.factor(ligand)~. ,
        data=.,
        method="xgbTree",
        trControl=control,
        metric="Accuracy",
        tuneGrid= xgbGrid,
        na.action = na.omit,
        nthread = 12,
        verbose = TRUE) -> XGB.model

end_time <- Sys.time()
print(end_time - start_time) # Time difference of

print(XGB.model)

saveRDS(XGB.model,"./XGB_models_final/xgb_ligands|aff_chop_Library|23_12_21.RDS")


rm(XGB.model)
gc()


###---------------------------------------###
## Modeling with library + chop + affMIN ## 
###---------------------------------------###

print("Modeling with library + chop + affMIN + rna")


start_time <- Sys.time()

Full_data %>% 
  dplyr::select(-c(ribo)) %>% 
  train(as.factor(ligand)~. ,
        data=.,
        method="xgbTree",
        trControl=control,
        metric="Accuracy",
        tuneGrid= xgbGrid,
        na.action = na.omit,
        nthread = 12,
        verbose = TRUE) -> XGB.model

end_time <- Sys.time()
print(end_time - start_time) # Time difference of

print(XGB.model)

saveRDS(XGB.model,"./XGB_models_final/xgb_ligands|aff_chop_Library_rna|23_12_21.RDS")


rm(XGB.model)
gc()


###---------------------------------------###
## Modeling with library + chop + affMIN ## 
###---------------------------------------###

print("Modeling with library + chop + affMIN + ribo")


start_time <- Sys.time()

Full_data %>% 
  dplyr::select(-c(rna)) %>% 
  train(as.factor(ligand)~. ,
        data=.,
        method="xgbTree",
        trControl=control,
        metric="Accuracy",
        tuneGrid= xgbGrid,
        na.action = na.omit,
        nthread = 12,
        verbose = TRUE) -> XGB.model

end_time <- Sys.time()
print(end_time - start_time) # Time difference of

print(XGB.model)

saveRDS(XGB.model,"./XGB_models_final/xgb_ligands|aff_chop_Library_ribo|23_12_21.RDS")


rm(XGB.model)
gc()

###---------------------------------------------###
## Modeling with library + RNA + chop + affMIN ## 
###---------------------------------------------###

print("Modeling with everything")


start_time <- Sys.time()


Full_data %>% 
  train(as.factor(ligand)~. ,
        data=.,
        method="xgbTree",
        trControl=control,
        metric="Accuracy",
        tuneGrid= xgbGrid,
        na.action = na.omit,
        nthread = 12,
        verbose = TRUE) -> XGB.model

end_time <- Sys.time()
print(end_time - start_time) # Time difference of

print(XGB.model)

saveRDS(XGB.model,"./XGB_models_final/xgb_ligands|aff_chop_Library_rna_ribo|24_04_2022.RDS")



rm(XGB.model)
gc()


