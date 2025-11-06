# package needed #
require_packages = c("MCMCpack","data.table","survival","splines",
                    "Matrix", "Rcpp","RcppEnsmallen","RcppArmadillo",
                    "methods","parallel","pracma","RSpectra",
                    "lokern","glmnet","mgcv","gamsel","fdapace",
                    "PtProcess", "MASS","xtable","ggplot2")

new_packages = lapply(require_packages, require, character.only = TRUE)
new_packages = require_packages[!unlist(new_packages)]
if(length(new_packages)) {install.packages(new_packages)}

