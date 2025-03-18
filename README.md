# GPAM
This repository contains the R codes for the paper "Generalized Point Process Additive Models".

The following R codes can be used to obtain the simulation results in the paper. 

run_Model1.R: run the simulation setting "Model 1" for one replication.

run_Model2.R: run the simulation setting "Model 2" for one replication.

run_Model3.R: run the simulation setting "Model 3" for one replication.

The results returned by these R codes include the AUC estimate and Beta estimate for each method.



Make sure to install all required packages (see below), in addition to the included penReg_0.1.0.tar.gz package.

To install penReg_0.1.0.tar.gz package, use install.packages('penReg_0.1.0.tar.gz',type='source',repos=NULL)

require(MCMCpack)

require(data.table)

require(survival)

require(splines)

require(Matrix)

require(Rcpp)

require(RcppEnsmallen)

require(RcppArmadillo)

require(methods)

require(parallel)

require(data.table)

require(pracma)

require(RSpectra)

require(lokern)

require(glmnet)

require(mgcv)

require(gamsel)

require(fdapace)

require(PtProcess)

require(MASS)







