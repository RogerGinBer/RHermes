---
title: "R Notebook"
output:
  html_document:
    df_print: paged
---




```r
# setwd("../") #For unknitted use of the chunk
myHermes <- RHermesExp()
myHermes <- setExpParam(myHermes,
                        params = ExpParam(ppm = 3, res = 120000,
                                          instr = "Orbitrap", minmz = 50,
                                          maxmz = 1200))
negHermes <- RHermesExp()
negHermes <- setExpParam(negHermes,
                         params = ExpParam(ppm = 3, res = 120000,
                                          instr = "Orbitrap", minmz = 50,
                                          maxmz = 1200, ion = "-"))
```



```r
myHermes <- setDB(myHermes, "custom", filename = "D:/merge_KEGG_ECMDB.csv") 
```

```
## Parsing the custom formula database
```

```r
myHermes <- remAd(myHermes, myHermes@metadata@ExpParam@adlist$adduct[-c(1:5)])
```

```
## This is the new adduct list:
```

```r
negHermes <- setDB(negHermes, "custom", filename = "D:/merge_KEGG_ECMDB.csv")
```

```
## Parsing the custom formula database
```

```r
negHermes <- remAd(negHermes, negHermes@metadata@ExpParam@adlist$adduct[-c(1,7)])
```

```
## This is the new adduct list:
```




```r
dir <- "D:/Ecoli_GaryMichi/20201109_MSH_ID-X_pHILIC_E-coli_MS1_Yanes-Lab/MS1_positive/"
fil <- list.files(dir, pattern = ".*pos.*.mzML", full.names = TRUE)[c(6,9)] #Blank and Sample
myHermes <- FileProc(myHermes, files = fil)
```

```
## Preprocessing steps, calculating possible ionic formulas and isotopic distributions
```

```
## 
```

```
##  done.
```

```
## 
```

```
## 
```

```
## 
```

```r
dir <- "D:/Ecoli_GaryMichi/20201109_MSH_ID-X_pHILIC_E-coli_MS1_Yanes-Lab/MS1_negative/"
fil <- list.files(dir, pattern = ".*neg.*.mzML", full.names = TRUE)[c(6,9)]
negHermes <- FileProc(negHermes, files = fil)
```

```
## 
```

```
##  done.
```

```
## 
```

```
## 
```

```
## 
```


```r
s <- getSOIpar("double")
myHermes <- SOIfinder(myHermes, params = s, fileID = c(1,2,2,2,2),
                      against = c(0,0,1,1,1))
```

```
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
```

```r
myHermes <- SOIcleaner(myHermes, soiid = 4, minint = 10000, FALSE)
```

```
## 
## 
## 
```

```r
myHermes <- SOIcleaner(myHermes, soiid = 5, minint = 10000, TRUE)
```

```
## 
## 
## 
## 
## 
## 
## 
```

```r
myHermes <- ISFproc(myHermes, id = 4)
```

```
## 
```

```r
myHermes <- ISFproc(myHermes, id = 5)
```

```
## 
```

```r
saveRDS(myHermes, file = "D:/HermesResults/Ecoli/PosSOIFilter.rds")


negHermes <- SOIfinder(negHermes, params = s, fileID = c(1,2,2,2,2),
                      against = c(0,0,1,1,1))
```

```
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
## 
```

```r
negHermes <- SOIcleaner(negHermes, soiid = 4, 10000, FALSE)
```

```
## 
## 
## 
```

```r
negHermes <- SOIcleaner(negHermes, soiid = 5, 10000, TRUE)
```

```
## 
## 
## 
## 
## 
## 
## 
```

```r
negHermes <- ISFproc(negHermes, id = 4)
```

```
## 
```

```r
negHermes <- ISFproc(negHermes, id = 5)
```

```
## 
```

```r
saveRDS(negHermes, file = "D:/HermesResults/Ecoli/NegSOIFilter.rds")
```



```r
myHermes <- genIL(myHermes, id = 4, par = ILParam(priorization = "yes", ad = c("M+H", "M+NH4", "M+CH3OH+H")))
```

```
## Now processing the IL:
```

```
## Prioritizing M+H
```

```
## Prioritizing M+NH4
```

```
## Prioritizing M+CH3OH+H
```

```
## Grouping the entries
```

```
## Tidying the entries
```

```
## Done!
```

```r
# exportIL(myHermes, 1, maxOver = 15, sepFiles = TRUE)

negHermes <- genIL(negHermes, id = 4, par = ILParam(priorization = "yes", ad = c("M-H")))
```

```
## Now processing the IL:
```

```
## Prioritizing M-H
```

```
## Grouping the entries
```

```
## Tidying the entries
```

```
## Done!
```

```r
# exportIL(negHermes, 1, maxOver = 15, sepFiles = TRUE)
```


```r
ms2f <- list.files("D:/Ecoli_GaryMichi/tMS2_positive/", pattern = ".*pos.*.mzML", full.names = TRUE)
myHermes <- MS2Proc(myHermes, 1, MS2files = ms2f, useDB = TRUE, referenceDB = "D:/MS2ID_20200824_202808.rds")
```

```
## Error in (function (classes, fdef, mtable) : unable to find an inherited method for function 'MS2Proc' for signature '"RHermesExp", "numeric", "character", "missing"'
```

```r
saveRDS(myHermes, file = "D:/HermesResults/Ecoli/PosResults.rds")


ms2f <- list.files("D:/Ecoli_GaryMichi/tMS2_negative/", pattern = ".*neg.*.mzML", full.names = TRUE)
negHermes <- MS2Proc(negHermes, 1, MS2files = ms2f, useDB = TRUE, referenceDB = "D:/MS2ID_20200824_202808.rds")
```

```
## Error in (function (classes, fdef, mtable) : unable to find an inherited method for function 'MS2Proc' for signature '"RHermesExp", "numeric", "character", "missing"'
```

```r
saveRDS(negHermes, file = "D:/HermesResults/Ecoli/NegResults.rds")
```

<!--  ```{r} -->
<!--  IL <- myHermes@data@MS2Exp[[1]]@IL -->

<!--  refined_entries <- c(1378, 1679, 2006, 2044, 2068, 2040, 428, 965, 1301, 1266, 446, 1021, 879, 1469, 1556, 1655, 487,551, 1629, 501, 582, 1908, 1456, 1309, 730, 1760, 1832, 1599, 1645, 1540, 1792, 227, 487, 551, 1371, 1540, 1608, 1141, 1348, 469, 564, 2046, 1375, 1491, 1435, 1953, 1268, 2087, 1339, 1929, 210, 2109) -->

<!--  any(duplicated(refined_entries)) -->

<!--  curatedIL <- IL@IL[unique(refined_entries), ] -->

<!--  lowintIL <- curatedIL[curatedIL$MaxInt < 100000,] -->

<!-- ``` -->



<!-- ## 3 m/z range analysis -->
<!-- ```{r warning=FALSE} -->
<!-- # setwd("../") #For unknitted use of the chunk -->
<!-- myHermes <- RHermesExp() -->
<!-- myHermes <- setExpParam(myHermes, -->
<!--                         params = ExpParam(ppm = 3, res = 120000, -->
<!--                                           instr = "Orbitrap", minmz = 50, -->
<!--                                           maxmz = 1200)) -->

<!-- negHermes <- RHermesExp() -->
<!-- negHermes <- setExpParam(negHermes, -->
<!--                          params = ExpParam(ppm = 3, res = 120000, -->
<!--                                           instr = "Orbitrap", minmz = 50, -->
<!--                                           maxmz = 1200, ion = "-")) -->

<!-- ``` -->


<!-- ```{r} -->
<!-- myHermes <- setDB(myHermes, "custom", filename = "D:/merge_KEGG_ECMDB.csv") -->
<!-- myHermes <- remAd(myHermes, c("M+IsoProp+H", "M+IsoProp+Na+H", "M+DMSO+H")) -->

<!-- negHermes <- setDB(negHermes, "custom", filename = "D:/merge_KEGG_ECMDB.csv") -->
<!-- negHermes <- remAd(negHermes, c("M+IsoProp+H", "M+IsoProp+Na+H", "M+DMSO+H")) -->

<!-- ``` -->



<!-- ```{r, message = 1, warning = FALSE} -->
<!-- dir <- "D:/Ecoli_GaryMichi/20201109_MSH_ID-X_pHILIC_E-coli_MS1_Yanes-Lab/MS1_positive/" -->
<!-- fil <- list.files(dir, pattern = ".*pos.*.mzML", full.names = TRUE)[5:10] -->
<!-- myHermes <- FileProc(myHermes, files = fil) -->
<!-- saveRDS(myHermes, file ="D:/Ecoli_GaryMichi/Ranges_MS1_pos_PL.rds") -->


<!-- dir <- "D:/Ecoli_GaryMichi/20201109_MSH_ID-X_pHILIC_E-coli_MS1_Yanes-Lab/MS1_negative/" -->
<!-- fil <- list.files(dir, pattern = ".*neg.*.mzML", full.names = TRUE)[5:10] -->
<!-- negHermes <- FileProc(negHermes, files = fil) -->
<!-- saveRDS(negHermes, file ="D:/Ecoli_GaryMichi/Ranges_MS1_neg_PL.rds") -->

<!-- ``` -->

<!-- ```{r, message = FALSE, warning=FALSE} -->
<!-- s <- getSOIpar("double") -->
<!-- myHermes <- SOIfinder(myHermes, params = s, fileID = c(1,2,3,4,5,6,3,5,5), -->
<!--                       against = c(0,0,0,0,0,0,5,3,6)) -->
<!-- saveRDS(myHermes, file = "D:/Ecoli_GaryMichi/Ranges_MS1_pos_SOI_filter.rds") -->


<!-- negHermes <- SOIfinder(negHermes, params = s, fileID = c(1,2,3,4,5,6,3,5,5), -->
<!--                       against = c(0,0,0,0,0,0,5,3,6)) -->
<!-- saveRDS(negHermes, file = "D:/Ecoli_GaryMichi/Ranges_MS1_neg_SOI_filter.rds") -->

<!-- ``` -->

## Low intensity High IT study


```r
myHermes <- MS2Proc(myHermes, 1, MS2files = "D:/Ecoli_GaryMichi/Curated_lowint/Curated_lowint/UNL_Ecoli_low-int_pos.mzML",
                    useDB = T, mincos = 0.1, referenceDB = "D:/MS2ID_B2R_20201113_083214.rds", sstype = "1scan")
```

```
## Starting MS/MS data importation, merging and sorting within the IL entries
```

```
## A total of 1208 (98.37%) entries were not covered in the experiment
```

```
## Starting superspectra generation. This may take a while...
```

```
## Retrieving MS2 spectra from the reference database
```

```
## Calculating Cosine similarities
```

```
## Done!
```

```r
myHermes@data@MS2Exp[[2]] <- myHermes@data@MS2Exp[[1]] 
myHermes <- MS2Proc(myHermes, 2, MS2files = "D:/Ecoli_GaryMichi/Curated_lowint/Curated_lowint/UNL_Ecoli_low-int_HCD50_pos.mzML",
                    useDB = T, mincos = 0.1, referenceDB = "D:/MS2ID_B2R_20201113_083214.rds")
```

```
## Error in (function (classes, fdef, mtable) : unable to find an inherited method for function 'MS2Proc' for signature '"RHermesExp", "numeric", "character", "missing"'
```

```r
saveRDS(myHermes, "D:/HermesResults/Ecoli/Lowint.rds")
```

## Labelled samples study


```r
myHermes <- RHermesExp()
myHermes <- setExpParam(myHermes,
                        params = ExpParam(ppm = 3, res = 120000,
                                          instr = "Orbitrap", minmz = 50,
                                          maxmz = 1200))
negHermes <- RHermesExp()
negHermes <- setExpParam(negHermes,
                         params = ExpParam(ppm = 3, res = 120000,
                                          instr = "Orbitrap", minmz = 50,
                                          maxmz = 1200, ion = "-"))
```



```r
myHermes <- setDB(myHermes, "custom", filename = "D:/merge_KEGG_ECMDB.csv") 
```

```
## Parsing the custom formula database
```

```r
myHermes <- remAd(myHermes, myHermes@metadata@ExpParam@adlist$adduct[-c(1:5)])
```

```
## This is the new adduct list:
```

```r
negHermes <- setDB(negHermes, "custom", filename = "D:/merge_KEGG_ECMDB.csv")
```

```
## Parsing the custom formula database
```

```r
negHermes <- remAd(negHermes, negHermes@metadata@ExpParam@adlist$adduct[-c(1,7)])
```

```
## This is the new adduct list:
```



```r
dir <- "D:/Ecoli_GaryMichi/20201109_MSH_ID-X_pHILIC_E-coli_MS1_Yanes-Lab/MS1_positive/"
fil <- list.files(dir, pattern = ".*pos.*.mzML", full.names = TRUE)[c(6,3)] #Blank and Sample
myHermes <- FileProc(myHermes, files = fil, labelled = TRUE)
```

```
## Preprocessing steps, calculating possible ionic formulas and isotopic distributions
```

```
## 
```

```
##  done.
```

```
## 
```

```
## 
```

```
## 
```

```r
dir <- "D:/Ecoli_GaryMichi/20201109_MSH_ID-X_pHILIC_E-coli_MS1_Yanes-Lab/MS1_negative/"
fil <- list.files(dir, pattern = ".*neg.*.mzML", full.names = TRUE)[c(6,3)]
negHermes <- FileProc(negHermes, files = fil, labelled = TRUE)
```

```
## 
```

```
##  done.
```

```
## 
```

```
## 
```

```
## 
```


```r
unlabPos <- readRDS("D:/HermesResults/Ecoli/PosResults.rds")
myHermes@data@SOI <- c(unlabPos@data@SOI[[2]], unlabPos@data@SOI[[4]], unlabPos@data@SOI[[5]])
myHermes@data@SOI[[1]]@filename <- myHermes@metadata@filenames[2]
myHermes@data@SOI[[2]]@filename <- myHermes@metadata@filenames[2]
myHermes@data@SOI[[3]]@filename <- myHermes@metadata@filenames[2]
saveRDS(myHermes, file ="D:/HermesResults/Ecoli/Pos13CSOI.rds")


unlabNeg <- readRDS("D:/HermesResults/Ecoli/NegResults.rds")
negHermes@data@SOI <- c(unlabNeg@data@SOI[[2]], unlabNeg@data@SOI[[4]], unlabNeg@data@SOI[[5]])
negHermes@data@SOI[[1]]@filename <- negHermes@metadata@filenames[2]
negHermes@data@SOI[[2]]@filename <- negHermes@metadata@filenames[2]
negHermes@data@SOI[[3]]@filename <- negHermes@metadata@filenames[2]
saveRDS(negHermes, file ="D:/HermesResults/Ecoli/Neg13CSOI.rds")
```


```r
sink("D:/HermesResults/Ecoli/sessionInfo.txt")
sessionInfo()
```

```
## R version 3.6.3 (2020-02-29)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 18363)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=Spanish_Spain.1252  LC_CTYPE=Spanish_Spain.1252    LC_MONETARY=Spanish_Spain.1252 LC_NUMERIC=C                  
## [5] LC_TIME=Spanish_Spain.1252    
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] plotly_4.9.2.1      BiocParallel_1.20.1 keras_2.2.5.0       data.table_1.12.8   RCy3_2.6.3          gtools_3.8.1       
##  [7] forcats_0.4.0       stringr_1.4.0       dplyr_0.8.4         purrr_0.3.3         readr_1.3.1         tidyr_1.0.2        
## [13] tibble_3.0.1        tidyverse_1.3.0     viridis_0.5.1       viridisLite_0.3.0   ggplot2_3.3.0       CHNOSZ_1.3.4       
## [19] magrittr_1.5        mzR_2.20.0          Rcpp_1.0.4.6        igraph_1.2.4.2      enviPat_2.4         RHermes_0.99.9007  
## 
## loaded via a namespace (and not attached):
##  [1] colorspace_1.4-1     ellipsis_0.3.0       XVector_0.26.0       base64enc_0.1-3      fs_1.3.1             rstudioapi_0.11     
##  [7] DT_0.13              fansi_0.4.1          lubridate_1.7.4      xml2_1.2.2           codetools_0.2-16     R.methodsS3_1.7.1   
## [13] ncdf4_1.17           knitr_1.28           zeallot_0.1.0        jsonlite_1.6.1       packrat_0.5.0        broom_0.5.6         
## [19] dbplyr_1.4.2         png_0.1-7            R.oo_1.23.0          tfruns_1.4           shinydashboard_0.7.1 graph_1.64.0        
## [25] shiny_1.5.0          compiler_3.6.3       httr_1.4.1           backports_1.1.7      assertthat_0.2.1     fastmap_1.0.1       
## [31] lazyeval_0.2.2       cli_2.0.2            later_1.0.0          visNetwork_2.0.9     htmltools_0.5.0      tools_3.6.3         
## [37] gtable_0.3.0         glue_1.4.0           rappdirs_0.3.1       Biobase_2.46.0       cellranger_1.1.0     slickR_0.4.9        
## [43] vctrs_0.2.4          Biostrings_2.54.0    RJSONIO_1.3-1.4      nlme_3.1-144         iterators_1.0.12     xfun_0.12           
## [49] rvest_0.3.5          mime_0.9             lifecycle_0.2.0      XML_3.99-0.3         zlibbioc_1.32.0      scales_1.1.0        
## [55] hms_0.5.3            doSNOW_1.0.18        promises_1.1.0       ProtGenerics_1.18.0  parallel_3.6.3       yaml_2.2.1          
## [61] reticulate_1.14      gridExtra_2.3        stringi_1.4.5        S4Vectors_0.24.3     tensorflow_2.0.0     foreach_1.4.7       
## [67] BiocGenerics_0.32.0  rlang_0.4.6          pkgconfig_2.0.3      evaluate_0.14        lattice_0.20-38      htmlwidgets_1.5.1   
## [73] tidyselect_1.0.0     R6_2.4.1             IRanges_2.20.2       snow_0.4-3           generics_0.0.2       DBI_1.1.0           
## [79] pillar_1.4.3         haven_2.2.0          whisker_0.4          withr_2.2.0          KEGGREST_1.26.1      modelr_0.1.5        
## [85] crayon_1.3.4         shinyWidgets_0.5.3   rmarkdown_2.1        grid_3.6.3           readxl_1.3.1         reprex_0.3.0        
## [91] digest_0.6.25        xtable_1.8-4         httpuv_1.5.2         R.utils_2.9.2        stats4_3.6.3         munsell_0.5.0
```

```r
paste("Date: ", date())
```

```
## [1] "Date:  Fri Jan 01 13:13:43 2021"
```

```r
sink()
```



