# GPAM
## This repository contains the R codes to replicate the results in the paper "Generalized Point Process Additive Models".


### (1) Make sure to install all required packages (see below), in addition to the included penReg_0.1.0.tar.gz package.

#### To install penReg_0.1.0.tar.gz package, run the following codes in R 

install.packages('penReg_0.1.0.tar.gz',type='source',repos=NULL)

#### To install all the other required packages, run the following codes in R 

require_packages = c("MCMCpack","data.table","survival","splines",
                    "Matrix", "Rcpp","RcppEnsmallen","RcppArmadillo",
                    "methods","parallel","pracma","RSpectra",
                    "lokern","glmnet","mgcv","gamsel","fdapace",
                    "PtProcess", "MASS","xtable","ggplot2")

new_packages = lapply(require_packages, require, character.only = TRUE)

new_packages = require_packages[!unlist(new_packages)]

if(length(new_packages)) {install.packages(new_packages)}


### (2) The following R codes can be used to replicate the results in simulation studies in the paper. 

#### After downloading all R codes into one folder and setting the working directory in R to this folder, each of the following R codes can be run independently to reproduce the corresponding results.

Table1.R: run this code to generate the results in Table 1 of the paper.

Table2.R: run this code to generate the results in Table 2 of the paper.

Table3.R: run this code to generate the results in Table 3 of the paper.

TableS1.R: run this code to generate the results in Table S1 of the paper.

TableS2.R: run this code to generate the results in Table S2 of the paper.

FigureS1andS2.R: run this code to generate Figures S1 and S2 of the paper.

### (3) The following R codes can be used to replicate the results in real data application in the paper. 
#### The ICU EHR dataset used in the paper can be downloaded from the website: https://physionet.org/content/challenge-2012/1.0.0/. Specifically, the following files are required for the analysis: Outcomes-a.txt, Outcomes-b.txt, Outcomes-c.txt, set-a.tar.gz, set-b.tar.gz, and set-c.tar.gz. These files can be downloaded individually from the *Folder Navigation* section at the bottom of the webpage (there is no need to download the large ZIP archive). After downloading the files and extracting all the .tar.gz archives, the following R code can be run to reproduce the corresponding results.

preprocess_mimic.R: this code has to be run **first** to prprocess the dataset before running the codes below. Remember to change the **path_to_data** variable in the code to specify the correct path to the folder containing all the data files, that is Outcomes-a.txt, Outcomes-b.txt, Outcomes-c.txt, set-a/, set-b/, and set-c/.

Figure1.R: after running preprocess_mimic.R, this code can be run to generate Figure 1 of the paper. 

analyze_mimic.R: after running preprocess_mimic.R, this code can be run to generate the AUCs and beta estimates for different methods. 

Figure2.R: after running analyze_mimic.R, this code can be run to generate Figure 2 of the paper. 











