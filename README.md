---
title: "README: scripts related to manuscript Addressing unmeasured confounders in cohort studies: Instrumental Variable method for a time-fixed exposure on an outcome trajectory"

output:
  html_document:
    toc: true
    toc_depth: 4
    toc_title: "Table of Contents"
date: "2023-11-02"
---


## Author 

Corresponding author : Kateline Le Bourdonnec,
kateline.le-bourdonnec@u-bordeaux.fr

## Documentation

The "Biom_journal.Proj" project contains the R codes used to write the paper "Addressing unmeasured confounders in cohort studies: Instrumental Variable method for a time-fixed exposure on an outcome trajectory" to find the same results. 

All the scripts are available on the GitHub link: 

``` r
https://github.com/KatelineLeBourdonnec/CodePaper_InstrumentalVariable
```
## Packages 

load packages :
```{r, message=F, warning=F, results = 'hide'}
options(repos = "https://cran.r-project.org")
install.packages("splines")
install.packages("MASS")
install.packages("ggplot2")
install.packages("fmsb")
install.packages("doParallel")
install.packages("foreach")
install.packages("marqLevAlg")
install.packages("mvtnorm")
install.packages("here")
install.packages("lcmm")
```

## How to use scripts 

### 1.For the application: [Script : application data](script_figure3.R)

As the real data used for our application is confidential and requires investigator agreements to be shared, we have not provided them. However, the data set called "data_long.txt" is a simulated dataset based on real application set.

The "script_figure3.R" file contains the script used to replicate figure3.

### 2. For the simulations: 

#### 2.a. Running the simulations

To run the simulations, you need to load geneData function (in [geneData](geneData.R) file) and bootstrap_param function (in [boostrap_param](boostrap_param.R) for binary exposure or [boostrap_param_lin](boostrap_param_lin.R) for continuous exposure). 

The simulations were run in parallel mode on a computing server to reduce calculation time. Four scripts were used to run the simulations:

- [Function_for_Cont_Log-Sub](CODE_BJ_Function_Simu_Cont_Log-Sub.R) for a continuous exposure and the substitution method; 

- [Function_for_Linear-Sub](CODE_BJ_Function_Simu_Log-Linear.R) for a binary exposure, a linear first stage and the substitution method;

- [Function_for_Log-Sub](CODE_BJ_Function_Simu_Log-Sub.R") for a binary exposure, a logistic first stage and the substitution method;

- [Function_for_Log-Res](CODE_BJ_Function_Simu_Res.R) for a binary exposure, a logistic first stage and the residual method. 

The outputs of these scripts (once run on a computing server), are provided in files:

- [cont_log_sub](cont_log_sub)  for a continuous exposure and the substitution method; 

- [linear_sub](linear_sub)  for a binary exposure, a linear first stage and the substitution method;

- [log_res](log_res) for a binary exposure, a logistic first stage and the residual method;

- [log_sub](log_sub) for a binary exposure, a logistic first stage and the substitution method. 

In these file, each .Rdata returns a list with in this order:

* CoefTrue,CoefTrueInter,sdTrue, sdTrueInter <- coefficients and standard errors of the simple effect of the exposure and the effect of the interaction with time (Inter) when the true regression model is estimated (except for log_res)

* CoefNaifInter0,CoefNaif0,sdNaif0,sdNaifInter0 <- coefficients and standard errors of the simple effect of the exposure and the effect of the interaction with time (Inter) in the naive regression model

* Coef2stageM,Coef2stageInter,sd2stageM,sd2stageInter <- coefficients and standard errors  of the simple effect of the exposure and the effect of the interaction with time (Inter) in the IV regression model

* r2_1,F_1 <- R² statistic and Fisher statistic for the first stage model (F stat only for linear first stage)

* var_corrig, se_corr <- lists of the corrected variance and of the corrected standard error: each one contains two elements: for the simple effect of the exposure and for the effect of the interaction with time (Inter)

* TC, TCint <- Coverage rate for coef2stageM and Coef2stageInter - not to be used. 

* VarIntra,VarInter,VarTot,VarIntraInt,VarInterInt,VarTotInt <- variance intra-bootstrap, inter-bootstrap, and total variance for coef2stageM and Coef2stageInter

* resu[[1]],resu[[2]] <- bootstrap coefficients for coef2stageM and Coef2stageInter - not to be used.

<font color="red"> All these statistical indicators were used for the analyses conducted, but only the naive and two-stages coefficients over-time are employed for the reproducibility of Figure 2 in script [Script : simulation data](script_figure2.R). </font>



If one wants to replicate the exact same results, the exact same seed should be used. The following table provides the seeds to be specified in the scripts that run the simulations (in  "For reproducibility" section): 

```{r, echo=F}
# table of job & seed
BINARY <- data.frame(
  ALPHA = c(2,3,4,2,3,4,2,3,4),
  SIZE =  rep(c("2000","6000","20K"),3),
  LINEAR = c("set.seed(543*rep+254*7958158)","set.seed(543*rep+254*7958259)","set.seed(543*rep+254*7958310)","set.seed(543*rep+254*7958768)","set.seed(543*rep+254*7958790)","set.seed(543*rep+254*7958790)","set.seed(543*rep+254*7958794)","set.seed(543*rep+254*7958792)","set.seed(543*rep+254*7958791)"),
  LOGSUB = c(rep("set.seed(543*rep)",9)),
  LOGRES = c(rep("set.seed(543*rep)",9))
)
BINARY
# for "rep" number of replication ranges 1 to 500

CONTINUOUS <- data.frame(
  ALPHA = c(0.5,1,0.5,1,0.5,1),
  SIZE =  rep(c("2000","6000","20K"),2),
  CONTINUOUS = c("/","set.seed(1984211)","set.seed(543*rep+254)",
                 "set.seed(94321)","set.seed(560)","set.seed(94321)"))
CONTINUOUS
```

The rep argument should be varying from 1 to 500. 

#### 2.b. Summarizing the simulation results : [Script : simulation data](Script_table2.R)

This script combines the simulation outputs in order to create the Tables and Figures displayed in the manuscript.

Different scenarios are represented :
            * Binary / Continuous
            
            * N = 2000, 6000, 20K
            
            * alpha = 2, 3, 4 / 0.5, 1
            
            * 4 methods : Naive, Log/Res, Line/Subs / Log/Subs
            

## Session info  

### 1. Personal Computer 

The code was written in R with the following software versions : 

R version 4.2.1 (2022-06-23 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default

locale:
[1] LC_COLLATE=French_France.utf8  LC_CTYPE=French_France.utf8    LC_MONETARY=French_France.utf8
[4] LC_NUMERIC=C                   LC_TIME=French_France.utf8    

attached base packages:
[1] parallel  splines   stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] lcmm_2.0.2        marqLevAlg_2.0.8  randtoolbox_2.0.4 rngWELL_0.10-9    mvtnorm_1.2-2    
 [6] survival_3.3-1    MASS_7.3-57       dplyr_1.1.3       here_1.0.1        gridExtra_2.3    
[11] ggpubr_0.6.0      ggplot2_3.4.4    

loaded via a namespace (and not attached):
 [1] pillar_1.9.0      compiler_4.2.1    iterators_1.0.14  tools_4.2.1       nlme_3.1-157     
 [6] lattice_0.20-45   lifecycle_1.0.3   tibble_3.2.1      gtable_0.3.3      pkgconfig_2.0.3  
[11] rlang_1.1.1       foreach_1.5.2     Matrix_1.5-4.1    cli_3.6.1         rstudioapi_0.14  
[16] withr_2.5.0       generics_0.1.3    vctrs_0.6.3       rprojroot_2.0.3   grid_4.2.1       
[21] tidyselect_1.2.0  glue_1.6.2        R6_2.5.1          rstatix_0.7.2     fansi_1.0.4      
[26] carData_3.0-5     farver_2.1.1      purrr_1.0.1       tidyr_1.3.0       car_3.1-2        
[31] magrittr_2.0.3    codetools_0.2-18  scales_1.2.1      backports_1.4.1   abind_1.4-5      
[36] colorspace_2.1-0  ggsignif_0.6.4    labeling_0.4.2    utf8_1.2.3        doParallel_1.0.17
[41] munsell_0.5.0     broom_1.0.5       crayon_1.5.2

### 2. Computing server: 

Computer time for this study was provided by the computing facilities of the MCIA (Mésocentre de Calcul Intensif Aquitain). 
