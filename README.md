# Multiple Myeloma DREAM Challenge (Through RandomForest, scRNAseq and survival analysis)

### Clinical context
Multiple Myeloma (MM) is a type of bone marrow cancer, called plasma cell myeloma and simply myeloma as well. The cause of this disease is the Carcinogenesis of effector B cell, which mature in bone marrow. MM have shown the existence of genetic heterogeneity in precursor stages, which leads to the need to classify the high-risk stage patient. Therefore this project specified this particular needs for MM.


### Data
Expression data:
```MMRF_CoMMpass_IA9_E74GTF_Salmon_entrezID_TPM_hg19.csv```
(https://www.synapse.org/#!Synapse:syn10573789)

Clinical data and labels:
```sc3_Training_ClinAnnotations.csv```
(https://www.synapse.org/#!Synapse:syn9926878)

###  Description of analysis

- The first part of the chanllenge is to construct a propre machine learning model for our multiple myeloma to prognosis the risk of patients. Through the previous research, the Random forest is recommanded[1], and some specialists recommend as well the lasso regression method for RNAseq expressed data, the class label that we use is the HR_FLAG in clinical informations. For prudence consideration, I also tried another model to see the final performance in this train model procedure.

- Then, for the interest of Biologie, I accomplished the scRNAseq analysis by using other public data and survival analysis.

## The utilization and description of the different script (based on Jupyter notebook and R)

- ```Xinyue-survivalPrediction-MultipleMyeloma.ipynb```
This script works with Jupyter notebook, which is a script of construction of RandomForest and other models, the function of specify step is commanded. The PDF format of this work is availble as well, but the interactive plot could not be visualized in PDF. (PS: Those interactive plots take more time than other commend, I will suggest not running this part of code.)

```Machine_learning_MM.R```
This script is based on R, with RandomForest model to confirm the result of jupyter analysis, please use this command if there is message error about protect(). The utilization of each function is written in script and the predict process is also possible as long as the same structure as input data. (df_preprocessing_and_merged.csv)

```sh
Rscript --max-ppsize=500000 Machine_learning_MM.R
```

- ```MM_scRNAseq_xinyue.R``` 
This script is based on scRNAseq data of MM-->  [GSE161195](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161195).[2] 

Command for decompression of data file :
```sh
tar -xf GSE161195_RAW.tar -C ./GSM_file/
``` 
Since what we interested mainly about the cluster of B cells, so the application of gene lists based on one part accroding to bibiography[3][4], and other part based on the top 20 important items of Randomforest model except the RNA5S items, since they are universal items in nucleous.

Attention!!! this analysis needs more than 100Go memory to run.

- ``` survival_analysis_MM_xinyue.Rmd```
This script contains the survival analysis according to the result of scRNAseq analysis, we select genes which are likely significantly presented at B cells (pig. gene_biblio and gene_ML). So that we have a gene list that we selected for survival analysis. 

## Result 
gene_list=c("SNORD89",  "AURKAIP1", "RNVU1-18","ELANE", "LINC02108","ANKRD11","NCOA3","CYLD","RNY1","RFESD","RNU1-1","RNU1-3")
The result of this gene list is based on randomforest.

gene_biblio=c ("KRAS","NRAS","FAM46C","DIS3","TP53", "BRAF", "TRAF3", "PRDM1", "CYLD", "RB1", "IRF4", "EGR1","HIST1H1E","ACTG1")
This list based on articles.[3][4]

After the observation with scRNAseq, we have final those genes:
gene_list=c("AURKAIP1","ANKRD11", "NCOA3", "CYLD", "PRDM1","IRF4","ACTG1"), then with survival analysis, we might conclude AURKAIP1 and PRDM1 could be the important feautures that improve performance of machine learn. Otherwise, Age, gender, PFS, OS and CYTO might help for machine learning model.

AURKAIP1 has positive regulation of proteolysis. Located in mitochondrion and nucleoplasm. has not a specify related with immunity, but with RNA5S, which might universal for every cell.

PRDM1 has relation with humoral immunity and plasma cells, biologically could be the biomarker of MM, which confirmed again with this analysis


## Conclusion and Discussion

We have around 66% accuracy with Randomforest model by using the expression RNAseq data, but it's not an ideal model for the selection of biomarkers, there is CNN might be a better choice; but due to the limit of time, the application doesn't succeed.

The selection of scRNAseq for MM is a tricky problem since there aren't much data available. After trying 3 datasets, I chose one that is likely to approach to this challenge, and also related to the treatment of MM. During the analysis, the data chosen contained around 200 samples, which make some difficulties for calculating with a personal PC, the effet batch hasn't been applied, so the cluster might not be ideal.

The selection of features is considered to use the RNAseq differential expressed analysis, but there isn't any interests result since the samples are all concerned about MM samples.

## Installation

All the dependency packages needed are written in every script, jupyter notebook needs mostly sklean, pandas, numpy, matplot. ect

## Work environement for R
R version 4.1.1 (2021-08-10)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Big Sur 11.2.3

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

locale:
[1] fr_FR.UTF-8/fr_FR.UTF-8/fr_FR.UTF-8/C/fr_FR.UTF-8/fr_FR.UTF-8

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods  
[9] base     

other attached packages:
 [1] caret_6.0-93          sp_1.5-0              SeuratObject_4.1.1   
 [4] Seurat_4.1.1          rms_6.3-0             SparseM_1.81         
 [7] Hmisc_4.7-1           Formula_1.2-4         lattice_0.20-45      
[10] GGally_2.1.2          devtools_2.4.4        usethis_2.1.6        
[13] survminer_0.4.9       ggpubr_0.4.0          survival_3.4-0       
[16] forcats_0.5.2         stringr_1.4.1         purrr_0.3.4          
[19] readr_2.1.2           tidyr_1.2.1           tibble_3.1.8         
[22] ggplot2_3.3.6         tidyverse_1.3.2       clusterProfiler_4.0.5
[25] org.Hs.eg.db_3.13.0   AnnotationDbi_1.54.1  IRanges_2.26.0       
[28] S4Vectors_0.30.2      Biobase_2.52.0        BiocGenerics_0.38.0  
[31] randomForest_4.7-1.1  dplyr_1.0.10         

loaded via a namespace (and not attached):
  [1] scattermore_0.8        ModelMetrics_1.2.2.2   bit64_4.0.5           
  [4] knitr_1.40             irlba_2.3.5            multcomp_1.4-20       
  [7] data.table_1.14.2      rpart_4.1.16           KEGGREST_1.32.0       
 [10] hardhat_1.2.0          RCurl_1.98-1.8         generics_0.1.3        
 [13] callr_3.7.2            cowplot_1.1.1          TH.data_1.1-1         
 [16] RSQLite_2.2.17         RANN_2.6.1             polspline_1.1.20      
 [19] shadowtext_0.1.2       future_1.28.0          bit_4.0.4             
 [22] tzdb_0.3.0             enrichplot_1.12.3      spatstat.data_2.2-0   
 [25] xml2_1.3.3             lubridate_1.8.0        httpuv_1.6.6          
 [28] assertthat_0.2.1       viridis_0.6.2          gargle_1.2.1          
 [31] gower_1.0.0            xfun_0.33              hms_1.1.2             
 [34] evaluate_0.16          promises_1.2.0.1       fansi_1.0.3           
 [37] dbplyr_2.2.1           readxl_1.4.1           km.ci_0.5-6           
 [40] igraph_1.3.4           DBI_1.1.3              htmlwidgets_1.5.4     
 [43] reshape_0.8.9          spatstat.geom_2.4-0    googledrive_2.0.0     
 [46] ellipsis_0.3.2         backports_1.4.1        deldir_1.0-6          
 [49] vctrs_0.4.1            ROCR_1.0-11            remotes_2.4.2         
 [52] quantreg_5.94          abind_1.4-5            cachem_1.0.6          
 [55] withr_2.5.0            ggforce_0.3.4          progressr_0.11.0      
 [58] checkmate_2.1.0        sctransform_0.3.4      treeio_1.16.2         
 [61] prettyunits_1.1.1      goftest_1.2-3          cluster_2.1.4         
 [64] DOSE_3.18.3            ape_5.6-2              lazyeval_0.2.2        
 [67] crayon_1.5.1           recipes_1.0.1          pkgconfig_2.0.3       
 [70] tweenr_2.0.2           GenomeInfoDb_1.28.4    nlme_3.1-159          
 [73] pkgload_1.3.0          nnet_7.3-17            rlang_1.0.5           
 [76] globals_0.16.1         lifecycle_1.0.2        miniUI_0.1.1.1        
 [79] MatrixModels_0.5-1     sandwich_3.0-2         downloader_0.4        
 [82] modelr_0.1.9           cellranger_1.1.0       polyclip_1.10-0       
 [85] matrixStats_0.62.0     lmtest_0.9-40          Matrix_1.5-1          
 [88] aplot_0.1.7            KMsurv_0.1-5           carData_3.0-5         
 [91] zoo_1.8-10             reprex_2.0.2           base64enc_0.1-3       
 [94] ggridges_0.5.3         processx_3.7.0         googlesheets4_1.0.1   
 [97] pheatmap_1.0.12        png_0.1-7              viridisLite_0.4.1     
[100] bitops_1.0-7           KernSmooth_2.23-20     pROC_1.18.0           
[103] Biostrings_2.60.2      blob_1.2.3             qvalue_2.24.0         
[106] spatstat.random_2.2-0  parallelly_1.32.1      jpeg_0.1-9            
[109] rstatix_0.7.0          gridGraphics_0.5-1     ggsignif_0.6.3        
[112] scales_1.2.1           ica_1.0-3              memoise_2.0.1         
[115] magrittr_2.0.3         plyr_1.8.7             zlibbioc_1.38.0       
[118] compiler_4.1.1         scatterpie_0.1.8       RColorBrewer_1.1-3    
[121] fitdistrplus_1.1-8     cli_3.4.0              urlchecker_1.0.1      
[124] XVector_0.32.0         listenv_0.8.0          pbapply_1.5-0         
[127] patchwork_1.1.2        ps_1.7.1               htmlTable_2.4.1       
[130] mgcv_1.8-40            MASS_7.3-58.1          tidyselect_1.1.2      
[133] stringi_1.7.8          yaml_2.3.5             GOSemSim_2.18.1       
[136] latticeExtra_0.6-30    ggrepel_0.9.1          survMisc_0.5.6        
[139] grid_4.1.1             fastmatch_1.1-3        tools_4.1.1           
[142] future.apply_1.9.1     rstudioapi_0.14        foreach_1.5.2         
[145] foreign_0.8-82         gridExtra_2.3          prodlim_2019.11.13    
[148] Rtsne_0.16             farver_2.1.1           ggraph_2.0.6          
[151] rgeos_0.5-9            digest_0.6.29          BiocManager_1.30.18   
[154] shiny_1.7.2            lava_1.6.10            Rcpp_1.0.9            
[157] car_3.1-0              broom_1.0.1            later_1.3.0           
[160] RcppAnnoy_0.0.19       httr_1.4.4             colorspace_2.0-3      
[163] tensor_1.5             reticulate_1.26        rvest_1.0.3           
[166] fs_1.5.2               splines_4.1.1          uwot_0.1.14           
[169] yulab.utils_0.0.5      spatstat.utils_2.3-1   tidytree_0.4.0        
[172] graphlayouts_0.8.1     ggplotify_0.1.0        plotly_4.10.0         
[175] sessioninfo_1.2.2      xtable_1.8-4           jsonlite_1.8.0        
[178] ggtree_3.0.4           tidygraph_1.2.2        timeDate_4021.104     
[181] ggfun_0.0.7            ipred_0.9-13           R6_2.5.1              
[184] profvis_0.3.7          pillar_1.8.1           htmltools_0.5.3       
[187] mime_0.12              glue_1.6.2             fastmap_1.1.0         
[190] BiocParallel_1.26.2    class_7.3-20           codetools_0.2-18      
[193] fgsea_1.18.0           pkgbuild_1.3.1         mvtnorm_1.1-3         
[196] utf8_1.2.2             spatstat.sparse_2.1-1  leiden_0.4.3          
[199] GO.db_3.13.0           interp_1.1-3           rmarkdown_2.16        
[202] munsell_0.5.0          DO.db_2.9              GenomeInfoDbData_1.2.6
[205] iterators_1.0.14       haven_2.5.1            reshape2_1.4.4        
[208] gtable_0.3.1           spatstat.core_2.4-4  

## Bibliography
[1] Allegra A, Tonacci A, Sciaccotta R, Genovese S, Musolino C, Pioggia G, Gangemi S. Machine Learning and Deep Learning Applications in Multiple Myeloma Diagnosis, Prognosis, and Treatment Selection. Cancers (Basel). 2022 Jan 25;14(3):606. doi: 10.3390/cancers14030606. PMID: 35158874; PMCID: PMC8833500.

[2]Cohen, Y.C., Zada, M., Wang, SY. et al. Identification of resistance pathways and therapeutic targets in relapsed multiple myeloma patients through single-cell sequencing. Nat Med 27, 491–503 (2021). https://doi.org/10.1038/s41591-021-01232-w

[3] Dutta, A.K., Alberge, JB., Sklavenitis-Pistofidis, R. et al. Single-cell profiling of tumour evolution in multiple myeloma — opportunities for precision medicine. Nat Rev Clin Oncol 19, 223–236 (2022). https://doi.org/10.1038/s41571-021-00593-y

[4]Kumar SK, Rajkumar V, Kyle RA, van Duin M, Sonneveld P, Mateos MV, Gay F, Anderson KC. Multiple myeloma. Nat Rev Dis Primers. 2017 Jul 20;3:17046. doi: 10.1038/nrdp.2017.46. PMID: 28726797.
