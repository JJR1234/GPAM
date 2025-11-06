## run Model T ##
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
print(args)
ffs = as.integer(args[1])
Tmax = as.integer(args[2])
set.seed(ffs)


require(Rcpp)
source("SS_FUN.R")
require(penReg)
source("ROCFUN.R")

setting = "ModelT"

lam_lo = 5
lam_up = 15

ncores = 1

source("simuT.R")

patient_sel = 1:n
feature_sel = 1:p

#### FPC scores ####

gram_distr = get_gram_distribution_vec(data_all=data_df,
                                       patient_sel=patient_sel,
                                       feature_sel=feature_sel, 
                                       ncores=ncores,
                                       Gaus_1=TRUE, Gaus_2=TRUE)

d_distr = 10

alpha_x_distr_list = lapply(seq_along(feature_sel), function(j){
    gram_distr_j = gram_distr[[j]]$gram
    get_PC_score(gram_distr_j,d=d_distr)
})


## regression analysis 
res_list = list()
test_id = 501:1000
get_grpid <- function(perG=5,p=100){
    G = p/perG
    grpid = cbind(seq(0, p-perG, by =perG), seq(perG-1, p-1, by =perG))
    rbind(c(0,0), grpid+1)
}

eps = 1e-8
maxiter = 100

lam_seq = exp(seq(log(1),log(0.0001),length=50))
ridge_seq = lam_seq/100

ntrain = 300
p = 50

iter_name = paste("ntrain=",ntrain,";p=",p,sep="")
train_id = 1:ntrain
print(iter_name)

## prepare the data
alpha_x_distr = do.call(cbind, alpha_x_distr_list[1:p])

Ytest = Y[test_id]
alpha_x_distr_test = alpha_x_distr[test_id, ]

Ytrain = Y[train_id]
alpha_x_distr_train = alpha_x_distr[train_id,]



## run regression with PC scores
W = rep(1,length(Ytrain))
foldid = get_foldid(Ytrain)
gvec = c(0, rep(1,p))
ridge_gvec = c(0, rep(1,p))


## distribution (GPAM)

Xtrain = cbind(1,alpha_x_distr_train)
Xtest = cbind(1,alpha_x_distr_test)
grpid = get_grpid(d_distr, ncol(alpha_x_distr) )

res = get_beta_auc_binom_wtest( W,  Xtrain, Ytrain,  Xtest,Ytest, 
                                foldid,  grpid, lam_seq, 
                                ridge_seq, ridge_gvec, gvec, 
                                eps, maxiter, ncores)

res_list[[paste(iter_name,";distr",sep="")]] = list(auc=res$auc, beta=res$beta_norm)



## results

sapply(res_list, function(x){x$auc})

filename = paste("resultSimu/", setting, Tmax, "_", ffs, ".rda", sep="" )
save(res_list, file=filename)
quit(save="no")



