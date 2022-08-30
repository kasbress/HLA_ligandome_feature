Test XGB models
================
Kaspar Bresser
24/12/2021

``` r
library(here)
library(caret)
library(xgboost)
library(tidyverse)
```

Import the feature library

``` r
feature.table <- read_tsv(here("Data", "Protein_per_Uniprot_entry_library_v2_RBP_GC_length_codon_AA_m6A_m5C_AtoI_m1A_m7G_CD8miRDB_PTM.csv")) %>% 
  mutate(across(everything(), replace_na, 0))

feature.table <- rename(feature.table, swissprot_id = Entry)
```

Get the file list, extract the names that will be used as column names

``` r
files <- list.files(here("XGB_models_final"), pattern = "xgb")
#files <- c("xgb_ligands|aff_chop_only|11_01_2022.RDS", "xgb_ligands|rna_aff_chop|11_01_2022.RDS")



files %>% 
  str_extract("\\|.*\\|") %>% 
  str_remove_all("\\|") -> names

files <- set_names(files, names)

files
```

    ##                                  aff_chop_Library_ribo 
    ##       "xgb_ligands|aff_chop_Library_ribo|23_12_21.RDS" 
    ##                              aff_chop_Library_rna_ribo 
    ## "xgb_ligands|aff_chop_Library_rna_ribo|24_04_2022.RDS" 
    ##                                   aff_chop_Library_rna 
    ##        "xgb_ligands|aff_chop_Library_rna|23_12_21.RDS" 
    ##                                       aff_chop_Library 
    ##            "xgb_ligands|aff_chop_Library|23_12_21.RDS" 
    ##                                          aff_chop_only 
    ##             "xgb_ligands|aff_chop_only|24_04_2022.RDS" 
    ##                                          aff_chop_ribo 
    ##             "xgb_ligands|aff_chop_ribo|24_04_2022.RDS" 
    ##                                      aff_chop_rna_ribo 
    ##         "xgb_ligands|aff_chop_rna_ribo|24_04_2022.RDS" 
    ##                                           aff_chop_rna 
    ##              "xgb_ligands|aff_chop_rna|24_04_2022.RDS" 
    ##                                               aff_only 
    ##                  "xgb_ligands|aff_only|24_04_2022.RDS" 
    ##                                           library_only 
    ##                "xgb_ligands|library_only|23_12_21.RDS" 
    ##                                                 random 
    ##                    "xgb_ligands|random|24_04_2022.RDS" 
    ##                                              ribo_only 
    ##                 "xgb_ligands|ribo_only|24_04_2022.RDS" 
    ##                                               rna_only 
    ##                  "xgb_ligands|rna_only|24_04_2022.RDS"

Import the XGB models.

``` r
files %>% 
  map(~read_rds(here("XGB_models_final", .))) %>% 
  set_names(names) -> XGB.models

XGB.models
```

    ## $aff_chop_Library_ribo
    ## eXtreme Gradient Boosting 
    ## 
    ## No pre-processing
    ## Resampling: Cross-Validated (10 fold, repeated 2 times) 
    ## Summary of sample sizes: 45000, 45000, 45000, 45000, 45000, 45000, ... 
    ## Addtional sampling using down-sampling
    ## 
    ## Resampling results:
    ## 
    ##   Accuracy  Kappa    
    ##   0.89807   0.7183847
    ## 
    ## Tuning parameter 'nrounds' was held constant at a value of 1000
    ## 
    ## Tuning parameter 'min_child_weight' was held constant at a value of 0.9
    ## 
    ## Tuning parameter 'subsample' was held constant at a value of 1
    ## 
    ## $aff_chop_Library_rna_ribo
    ## eXtreme Gradient Boosting 
    ## 
    ## No pre-processing
    ## Resampling: Cross-Validated (10 fold, repeated 2 times) 
    ## Summary of sample sizes: 45000, 45000, 45000, 45000, 45000, 45000, ... 
    ## Addtional sampling using down-sampling
    ## 
    ## Resampling results:
    ## 
    ##   Accuracy  Kappa    
    ##   0.90184   0.7275291
    ## 
    ## Tuning parameter 'nrounds' was held constant at a value of 1000
    ## 
    ## Tuning parameter 'min_child_weight' was held constant at a value of 0.9
    ## 
    ## Tuning parameter 'subsample' was held constant at a value of 1
    ## 
    ## $aff_chop_Library_rna
    ## eXtreme Gradient Boosting 
    ## 
    ## No pre-processing
    ## Resampling: Cross-Validated (10 fold, repeated 2 times) 
    ## Summary of sample sizes: 45000, 45000, 45000, 45000, 45000, 45000, ... 
    ## Addtional sampling using down-sampling
    ## 
    ## Resampling results:
    ## 
    ##   Accuracy  Kappa    
    ##   0.90131   0.7262852
    ## 
    ## Tuning parameter 'nrounds' was held constant at a value of 1000
    ## 
    ## Tuning parameter 'min_child_weight' was held constant at a value of 0.9
    ## 
    ## Tuning parameter 'subsample' was held constant at a value of 1
    ## 
    ## $aff_chop_Library
    ## eXtreme Gradient Boosting 
    ## 
    ## No pre-processing
    ## Resampling: Cross-Validated (10 fold, repeated 2 times) 
    ## Summary of sample sizes: 45000, 45000, 45000, 45000, 45000, 45000, ... 
    ## Addtional sampling using down-sampling
    ## 
    ## Resampling results:
    ## 
    ##   Accuracy  Kappa    
    ##   0.89114   0.7012391
    ## 
    ## Tuning parameter 'nrounds' was held constant at a value of 1000
    ## 
    ## Tuning parameter 'min_child_weight' was held constant at a value of 0.9
    ## 
    ## Tuning parameter 'subsample' was held constant at a value of 1
    ## 
    ## $aff_chop_only
    ## eXtreme Gradient Boosting 
    ## 
    ## No pre-processing
    ## Resampling: Cross-Validated (10 fold, repeated 2 times) 
    ## Summary of sample sizes: 45000, 45000, 45000, 45000, 45000, 45000, ... 
    ## Addtional sampling using down-sampling
    ## 
    ## Resampling results:
    ## 
    ##   Accuracy  Kappa    
    ##   0.85906   0.6285168
    ## 
    ## Tuning parameter 'nrounds' was held constant at a value of 1000
    ## 
    ## Tuning parameter 'min_child_weight' was held constant at a value of 0.9
    ## 
    ## Tuning parameter 'subsample' was held constant at a value of 1
    ## 
    ## $aff_chop_ribo
    ## eXtreme Gradient Boosting 
    ## 
    ## No pre-processing
    ## Resampling: Cross-Validated (10 fold, repeated 2 times) 
    ## Summary of sample sizes: 45000, 45000, 45000, 45000, 45000, 45000, ... 
    ## Addtional sampling using down-sampling
    ## 
    ## Resampling results:
    ## 
    ##   Accuracy  Kappa    
    ##   0.88788   0.6951733
    ## 
    ## Tuning parameter 'nrounds' was held constant at a value of 1000
    ## 
    ## Tuning parameter 'min_child_weight' was held constant at a value of 0.9
    ## 
    ## Tuning parameter 'subsample' was held constant at a value of 1
    ## 
    ## $aff_chop_rna_ribo
    ## eXtreme Gradient Boosting 
    ## 
    ## No pre-processing
    ## Resampling: Cross-Validated (10 fold, repeated 2 times) 
    ## Summary of sample sizes: 45000, 45000, 45000, 45000, 45000, 45000, ... 
    ## Addtional sampling using down-sampling
    ## 
    ## Resampling results:
    ## 
    ##   Accuracy  Kappa    
    ##   0.8985    0.7200657
    ## 
    ## Tuning parameter 'nrounds' was held constant at a value of 1000
    ## 
    ## Tuning parameter 'min_child_weight' was held constant at a value of 0.9
    ## 
    ## Tuning parameter 'subsample' was held constant at a value of 1
    ## 
    ## $aff_chop_rna
    ## eXtreme Gradient Boosting 
    ## 
    ## No pre-processing
    ## Resampling: Cross-Validated (10 fold, repeated 2 times) 
    ## Summary of sample sizes: 45000, 45000, 45000, 45000, 45000, 45000, ... 
    ## Addtional sampling using down-sampling
    ## 
    ## Resampling results:
    ## 
    ##   Accuracy  Kappa    
    ##   0.89723   0.7168193
    ## 
    ## Tuning parameter 'nrounds' was held constant at a value of 1000
    ## 
    ## Tuning parameter 'min_child_weight' was held constant at a value of 0.9
    ## 
    ## Tuning parameter 'subsample' was held constant at a value of 1
    ## 
    ## $aff_only
    ## eXtreme Gradient Boosting 
    ## 
    ## No pre-processing
    ## Resampling: Cross-Validated (10 fold, repeated 2 times) 
    ## Summary of sample sizes: 45000, 45000, 45000, 45000, 45000, 45000, ... 
    ## Addtional sampling using down-sampling
    ## 
    ## Resampling results:
    ## 
    ##   Accuracy  Kappa    
    ##   0.83528   0.5756165
    ## 
    ## Tuning parameter 'nrounds' was held constant at a value of 1000
    ## 
    ## Tuning parameter 'min_child_weight' was held constant at a value of 0.9
    ## 
    ## Tuning parameter 'subsample' was held constant at a value of 1
    ## 
    ## $library_only
    ## eXtreme Gradient Boosting 
    ## 
    ## No pre-processing
    ## Resampling: Cross-Validated (10 fold, repeated 2 times) 
    ## Summary of sample sizes: 45000, 45000, 45000, 45000, 45000, 45000, ... 
    ## Addtional sampling using down-sampling
    ## 
    ## Resampling results:
    ## 
    ##   Accuracy  Kappa    
    ##   0.70766   0.3168004
    ## 
    ## Tuning parameter 'nrounds' was held constant at a value of 1000
    ## 
    ## Tuning parameter 'min_child_weight' was held constant at a value of 0.9
    ## 
    ## Tuning parameter 'subsample' was held constant at a value of 1
    ## 
    ## $random
    ## eXtreme Gradient Boosting 
    ## 
    ## No pre-processing
    ## Resampling: Cross-Validated (10 fold, repeated 2 times) 
    ## Summary of sample sizes: 45000, 45000, 45000, 45000, 45000, 45000, ... 
    ## Addtional sampling using down-sampling
    ## 
    ## Resampling results:
    ## 
    ##   Accuracy  Kappa        
    ##   0.49949   -0.0003384882
    ## 
    ## Tuning parameter 'nrounds' was held constant at a value of 1000
    ## 
    ## Tuning parameter 'min_child_weight' was held constant at a value of 0.9
    ## 
    ## Tuning parameter 'subsample' was held constant at a value of 1
    ## 
    ## $ribo_only
    ## eXtreme Gradient Boosting 
    ## 
    ## No pre-processing
    ## Resampling: Cross-Validated (10 fold, repeated 2 times) 
    ## Summary of sample sizes: 45000, 45000, 45000, 45000, 45000, 45000, ... 
    ## Addtional sampling using down-sampling
    ## 
    ## Resampling results:
    ## 
    ##   Accuracy  Kappa    
    ##   0.66329   0.2903702
    ## 
    ## Tuning parameter 'nrounds' was held constant at a value of 1000
    ## 
    ## Tuning parameter 'min_child_weight' was held constant at a value of 0.9
    ## 
    ## Tuning parameter 'subsample' was held constant at a value of 1
    ## 
    ## $rna_only
    ## eXtreme Gradient Boosting 
    ## 
    ## No pre-processing
    ## Resampling: Cross-Validated (10 fold, repeated 2 times) 
    ## Summary of sample sizes: 45000, 45000, 45000, 45000, 45000, 45000, ... 
    ## Addtional sampling using down-sampling
    ## 
    ## Resampling results:
    ## 
    ##   Accuracy  Kappa   
    ##   0.64806   0.276133
    ## 
    ## Tuning parameter 'nrounds' was held constant at a value of 1000
    ## 
    ## Tuning parameter 'min_child_weight' was held constant at a value of 0.9
    ## 
    ## Tuning parameter 'subsample' was held constant at a value of 1

get column names.

``` r
here("Output", "test_train_tables_new", "rf_test_peptides.tsv") %>% 
  read_tsv( n_max = 1) %>% 
  names() -> column.names
```

``` r
here("Data", "rna_seq_data", "all_lines_final.rds") %>% 
  read_rds() %>% 
  unnest(expr) -> expr.data
```

Will import the test set in chunks, perform the predictions with each
model and write out the scores in a single pipe to keep from depleting
memory. The test set contains 2540538 entries.

``` r
2540538/33
```

    ## [1] 76986

``` r
2540538-76986
```

    ## [1] 2463552

``` r
seq(0, 2463552, 76986)
```

    ##  [1]       0   76986  153972  230958  307944  384930  461916  538902  615888
    ## [10]  692874  769860  846846  923832 1000818 1077804 1154790 1231776 1308762
    ## [19] 1385748 1462734 1539720 1616706 1693692 1770678 1847664 1924650 2001636
    ## [28] 2078622 2155608 2232594 2309580 2386566 2463552

``` r
start.points <- seq(0, 2463552, 76986)+1
```

Define a function that will import a portion of the test data and join
with the feature library.

Next, loop over the start points. In each loop use the subset of the
test data to perform predictions with all imported models. Column-bind
this to the ligand column of the same partition of the data, and write
out.

``` r
importstuff <- function(start){
  here("Output", "test_train_tables_new", "rf_test_peptides.tsv") %>% 
    read_tsv(skip = start , n_max = 76986, col_names = column.names) %>% 
    mutate(tumor = fct_recode(tumor, M026.X1 = "MEL1", SKMEL95 = "MEL2", NKIRTIL006 = "MEL3")) %>% 
    dplyr::select(-rna) %>% 
    left_join(expr.data, by = c("tumor", "swissprot_id")) %>% 
    inner_join(feature.table)  -> out
  return(out)
}


for (i in start.points){
    XGB.models %>% 
      map2_dfc(list(importstuff(i)) , ~predict(object = .x, newdata = .y, type = "prob")$`TRUE`) %>% 
      bind_cols(read_tsv(here("Output", "test_train_tables_new", "rf_test_peptides.tsv"),
                         skip = i, 
                         n_max = 76986, 
                         col_names = column.names, 
                         col_select = c(tumor, ligand))) %>% 
      write_tsv(here("Output", "XGB_final", paste0("prediction_results_", i, ".tsv")))
  gc()
}
```

