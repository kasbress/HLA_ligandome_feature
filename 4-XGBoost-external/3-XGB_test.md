Test XGB models on Broad data
================
Kaspar Bresser
11/01/2021

``` r
library(here)
library(caret)
library(xgboost)
library(tidyverse)
```

Import the feature library

``` r
feature.table <- read_tsv(here( "Protein_per_Uniprot_entry_library_v2_RBP_GC_length_codon_AA_m6A_m5C_AtoI_m1A_m7G_CD8miRDB_PTM.csv")) %>% 
  mutate(across(everything(), replace_na, 0))

feature.table <- rename(feature.table, ID = Entry)
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

    ##                            aff_chop_Library_rna 
    ## "xgb_ligands|aff_chop_Library_rna|23_12_21.RDS" 
    ##                                aff_chop_Library 
    ##     "xgb_ligands|aff_chop_Library|23_12_21.RDS" 
    ##                                   aff_chop_only 
    ##      "xgb_ligands|aff_chop_only|24_04_2022.RDS" 
    ##                                    aff_chop_rna 
    ##       "xgb_ligands|aff_chop_rna|24_04_2022.RDS" 
    ##                                        aff_only 
    ##           "xgb_ligands|aff_only|24_04_2022.RDS" 
    ##                                Library_chop_aff 
    ##     "xgb_ligands|Library_chop_aff|23_12_21.RDS"

Import the XGB models.

``` r
files %>% 
  map(~read_rds(here("XGB_models_final", .))) %>% 
  set_names(names) -> XGB.models

XGB.models
```

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
    ## $Library_chop_aff
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
    ##   0.89374   0.7073627
    ## 
    ## Tuning parameter 'nrounds' was held constant at a value of 1000
    ## 
    ## Tuning parameter 'min_child_weight' was held constant at a value of 0.9
    ## 
    ## Tuning parameter 'subsample' was held constant at a value of 1

get column names.

``` r
here("Output", "Broad_test_set_forXGB.tsv") %>% 
  read_tsv( n_max = 1) %>% 
  names() -> column.names

column.names
```

    ## [1] "ID"               "rna"              "Peptide"          "processing_score"
    ## [5] "ligand"           "Rank"             "allele"

Will import the test set in chunks, perform the predictions with each
model and write out the scores in a single pipe to keep from depleting
memory. The test set contains 2,802,729 entries. Read in by chunks of
75,000

``` r
2802729/75000
```

    ## [1] 37.36972

``` r
75000*38
```

    ## [1] 2850000

``` r
start.points <- seq(0, 2802729, 75000)+1

start.points
```

    ##  [1]       1   75001  150001  225001  300001  375001  450001  525001  600001
    ## [10]  675001  750001  825001  900001  975001 1050001 1125001 1200001 1275001
    ## [19] 1350001 1425001 1500001 1575001 1650001 1725001 1800001 1875001 1950001
    ## [28] 2025001 2100001 2175001 2250001 2325001 2400001 2475001 2550001 2625001
    ## [37] 2700001 2775001

Define a function that will import a portion of the test data,
rename/retain the columns needed, and join with the feature library.

Next, loop over the start points. In each loop use the subset of the
test data to perform predictions with all imported models. Column-bind
this to the ligand column of the same partition of the data, and write
out.

``` r
importstuff <- function(start){
  here("Output", "Broad_test_set_forXGB.tsv") %>% 
    read_tsv(skip = start , n_max = 75000, col_names = column.names) %>% 
    transmute(rna = rna, chop = processing_score, affMIN = Rank, ID = ID) %>% 
    left_join(feature.table, by = "ID") -> out
  return(out)
}


for(i in start.points){
    XGB.models %>% 
      map2_dfc(list(importstuff(i)) , ~predict(object = .x, newdata = .y, type = "prob")$`TRUE`) %>% 
      bind_cols(read_tsv(here("Output", "Broad_test_set_forXGB.tsv"),
                         skip = i, 
                         n_max = 75000, 
                         col_names = column.names, 
                         col_select = c(allele, ligand))) %>% 
      write_tsv(here("Output", "XGB_final_plus", paste0("prediction_results_", i, ".tsv")))
  gc()
}
```
