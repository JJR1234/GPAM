# GPAM
### This repository contains the R codes for the paper "Generalized Point Process Additive Models".

The following R codes can be used to obtain the simulation results in the paper. 

run_Model1.R: run the simulation setting "Model 1" for one replication.

run_Model2.R: run the simulation setting "Model 2" for one replication.

run_Model3.R: run the simulation setting "Model 3" for one replication.

The results returned by these R codes include the AUC estimate and Beta estimate for each method.




### Make sure to install all required packages (see below), in addition to the included penReg_0.1.0.tar.gz package.

## To install penReg_0.1.0.tar.gz package, run the following codes in R 

install.packages('penReg_0.1.0.tar.gz',type='source',repos=NULL)

## To install all the other required packages, run the following codes in R 

require_packages = c("MCMCpack","data.table","survival","splines",
                    "Matrix", "Rcpp","RcppEnsmallen","RcppArmadillo",
                    "methods","parallel","pracma","RSpectra",
                    "lokern","glmnet","mgcv","gamsel","fdapace",
                    "PtProcess", "MASS","adestr","ADCT")

new_packages = lapply(require_packages, require, character.only = TRUE)

new_packages = require_packages[!unlist(new_packages)]

if(length(new_packages)) {install.packages(new_packages)}








