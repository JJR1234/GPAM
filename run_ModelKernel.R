## run simu: sensitivity to kernel choice ##
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
print(args)
ffs = as.integer(args[1])
set.seed(ffs)

require(Rcpp)
source("SS_FUN.R")
require(penReg)
source("ROCFUN.R")

setting = "ModelKernel"

lam_lo = 5
lam_up = 15

ncores = 1

source("simu1.R")

patient_sel = 1:n
feature_sel = 1:p


gram_distr1 = get_gram_distribution_vec(data_all=data_df,
                                        patient_sel=patient_sel,
                                        feature_sel=feature_sel, 
                                        ncores=ncores,
                                        Gaus_1=TRUE, Gaus_2=TRUE)

gram_distr2 = get_gram_distribution_vec(data_all=data_df,
                                        patient_sel=patient_sel,
                                        feature_sel=feature_sel, 
                                        ncores=ncores,
                                        Gaus_1=FALSE, Gaus_2=TRUE)

gram_distr3 = get_gram_distribution_vec(data_all=data_df,
                                        patient_sel=patient_sel,
                                        feature_sel=feature_sel, 
                                        ncores=ncores,
                                        Gaus_1=TRUE, Gaus_2=FALSE)

gram_distr4 = get_gram_distribution_vec(data_all=data_df,
                                        patient_sel=patient_sel,
                                        feature_sel=feature_sel, 
                                        ncores=ncores,
                                        Gaus_1=FALSE, Gaus_2=FALSE)

d_distr = 10


alpha_x_distr1_list = lapply(seq_along(feature_sel), function(j){
    gram_distr_j = gram_distr1[[j]]$gram
    get_PC_score(gram_distr_j,d=d_distr)
})

alpha_x_distr2_list = lapply(seq_along(feature_sel), function(j){
    gram_distr_j = gram_distr2[[j]]$gram
    get_PC_score(gram_distr_j,d=d_distr)
})

alpha_x_distr3_list = lapply(seq_along(feature_sel), function(j){
    gram_distr_j = gram_distr3[[j]]$gram
    get_PC_score(gram_distr_j,d=d_distr)
})

alpha_x_distr4_list = lapply(seq_along(feature_sel), function(j){
    gram_distr_j = gram_distr4[[j]]$gram
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
alpha_x_distr1 = do.call(cbind, alpha_x_distr1_list[1:p])
alpha_x_distr2 = do.call(cbind, alpha_x_distr2_list[1:p])
alpha_x_distr3 = do.call(cbind, alpha_x_distr3_list[1:p])
alpha_x_distr4 = do.call(cbind, alpha_x_distr4_list[1:p])

Ytest = Y[test_id]
alpha_x_distr1_test = alpha_x_distr1[test_id, ]
alpha_x_distr2_test = alpha_x_distr2[test_id, ]
alpha_x_distr3_test = alpha_x_distr3[test_id, ]
alpha_x_distr4_test = alpha_x_distr4[test_id, ]

Ytrain = Y[train_id]
alpha_x_distr1_train = alpha_x_distr1[train_id,]
alpha_x_distr2_train = alpha_x_distr2[train_id,]
alpha_x_distr3_train = alpha_x_distr3[train_id,]
alpha_x_distr4_train = alpha_x_distr4[train_id,]


## run regression with PC scores

W = rep(1,length(Ytrain))
foldid = get_foldid(Ytrain)
gvec = c(0, rep(1,p))
ridge_gvec = c(0, rep(1,p))


# distribution 

Xtrain = cbind(1,alpha_x_distr1_train)
Xtest = cbind(1,alpha_x_distr1_test)
grpid = get_grpid(d_distr, ncol(alpha_x_distr1) )

res = get_beta_auc_binom_wtest( W,  Xtrain, Ytrain,  Xtest,Ytest, 
                                foldid,  grpid, lam_seq, 
                                ridge_seq, ridge_gvec, gvec, 
                                eps, maxiter, ncores)

res_list[[paste(iter_name,";distr1",sep="")]] = list(auc=res$auc, beta=res$beta_norm)

# distribution 
Xtrain = cbind(1,alpha_x_distr2_train)
Xtest = cbind(1,alpha_x_distr2_test)
grpid = get_grpid(d_distr, ncol(alpha_x_distr2) )

res = get_beta_auc_binom_wtest( W,  Xtrain, Ytrain,  Xtest,Ytest, 
                                foldid,  grpid, lam_seq, 
                                ridge_seq, ridge_gvec, gvec, 
                                eps, maxiter, ncores)

res_list[[paste(iter_name,";distr2",sep="")]] = list(auc=res$auc, beta=res$beta_norm)

# distribution 
Xtrain = cbind(1,alpha_x_distr3_train)
Xtest = cbind(1,alpha_x_distr3_test)
grpid = get_grpid(d_distr, ncol(alpha_x_distr3) )

res = get_beta_auc_binom_wtest( W,  Xtrain, Ytrain,  Xtest,Ytest, 
                                foldid,  grpid, lam_seq, 
                                ridge_seq, ridge_gvec, gvec, 
                                eps, maxiter, ncores)

res_list[[paste(iter_name,";distr3",sep="")]] = list(auc=res$auc, beta=res$beta_norm)


# distribution 
Xtrain = cbind(1,alpha_x_distr4_train)
Xtest = cbind(1,alpha_x_distr4_test)
grpid = get_grpid(d_distr, ncol(alpha_x_distr4) )

res = get_beta_auc_binom_wtest( W,  Xtrain, Ytrain,  Xtest,Ytest, 
                                foldid,  grpid, lam_seq, 
                                ridge_seq, ridge_gvec, gvec, 
                                eps, maxiter, ncores)

res_list[[paste(iter_name,";distr4",sep="")]] = list(auc=res$auc, beta=res$beta_norm)


## results

sapply(res_list, function(x){x$auc})

filename = paste("resultSimu/", setting, "_", ffs, ".rda", sep="" )
save(res_list, file=filename)
quit(save="no")



