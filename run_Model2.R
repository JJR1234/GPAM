## run Model 2 ##
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
print(args)
ffs = as.integer(args[1])
set.seed(ffs)


require(Rcpp)
source("SS_FUN.R")
require(penReg)
source("ROCFUN.R")

setting = "Model2"

lam_lo = 10
lam_up = 10

ncores = 1

source("simu1.R")

patient_sel = 1:n
feature_sel = 1:p

##### count data ######
count_value = data_df[,.(count=.N), by=c("subject_num","feature_id")]
count_value[, i:= match(subject_num, patient_sel) ]
count_value[, j:= match(feature_id, feature_sel) ]

mat_dx = sparseMatrix(i=count_value[,i], j=count_value[,j], x=count_value[,count],
                      dims = c(length(patient_sel), length(feature_sel)))
## mat_dx@x = log(mat_dx@x+1)
colnames(mat_dx) = paste("code",feature_sel,sep="")
mat_dx = as.matrix(mat_dx)

#### FPC scores ####

kde_score = get_KDE_PCscore(data_all=data_df, patient_sel=patient_sel,
                            feature_sel=feature_sel, ncores=ncores)

gram_distr = get_gram_distribution_vec(data_all=data_df,
                                       patient_sel=patient_sel,
                                       feature_sel=feature_sel, 
                                       ncores=ncores,
                                       Gaus_1=TRUE, Gaus_2=TRUE)

d_fpca = 5
d_distr = 10

alpha_x_distr_list = lapply(seq_along(feature_sel), function(j){
    gram_distr_j = gram_distr[[j]]$gram
    get_PC_score(gram_distr_j,d=d_distr)
})


alpha_fpca_list = lapply(seq_along(feature_sel), function(j){
    if(ncol(kde_score[[j]]$score) >= d_fpca){
        scale(cbind(kde_score[[j]]$score[,1:d_fpca], mat_dx[,j]) )
    }else{
        scale(cbind(kde_score[[j]]$score, mat_dx[,j]) )
    }
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

for(ntrain in c(100,200,300)){
    for(p in c(50,100,150)){
        
        iter_name = paste("ntrain=",ntrain,";p=",p,sep="")
        train_id = 1:ntrain
        print(iter_name)
        
        ## prepare the data
        grpid_fpca = sapply(alpha_fpca_list[1:p], function(x){
            ncol(x)
        })
        grpid_fpca = c(1, grpid_fpca)
        grpid_fpca = cbind(c(0, cumsum(grpid_fpca)[-length(grpid_fpca)]), 
                           cumsum(grpid_fpca)-1)
        
        
        alpha_x_distr = do.call(cbind, alpha_x_distr_list[1:p])
        alpha_fpca = do.call(cbind, alpha_fpca_list[1:p])
        
        Ytest = Y[test_id]
        alpha_x_distr_test = alpha_x_distr[test_id, ]
        mat_dx_test = mat_dx[test_id,1:p]
        alpha_fpca_test = alpha_fpca[test_id,]
        
        Ytrain = Y[train_id]
        alpha_x_distr_train = alpha_x_distr[train_id,]
        mat_dx_train = mat_dx[train_id,1:p]
        alpha_fpca_train = alpha_fpca[train_id,]
        
        
        ### run glmnet with count
        
        res = cv.glmnet(x=mat_dx_train, y=Ytrain, family='binomial',  alpha=0.9,
                        lambda.min.ratio=0.001,  nfolds=5)
        beta_glmnet = coef(res)[,1]
        
        xb =  cbind(1,mat_dx_test) %*% beta_glmnet
        #loglik_glmnet = get_binom_loglik(rep(1,length(Ytest))/length(Ytest), xb, Ytest )
        if(length(unique(xb))<2){
            auc_glmnet = 0.5 
        }else{
            auc_glmnet = ROC.Est.FUN(Ytest, xb, fpr0=0.05)[1]
        }
        
        res_list[[paste(iter_name,";glmnet",sep="")]] = list(auc=auc_glmnet, beta=beta_glmnet)
        
        ## run GAM with count
        degrees = get_df_gamsel(mat_dx_train)
        dfs =  degrees - 2
        dfs[dfs<1] = 1
        gamsel.cv = cv.gamsel(mat_dx_train,Ytrain,degrees=degrees,dfs=dfs,family = "binomial")
        #gamsel.cv = cv.gamsel(mat_dx,Y,bases=bases,family = "binomial")
        lam_min = gamsel.cv$lambda[gamsel.cv$index.1se]
        #gamsel.out = gamsel(mat_dx,Y,bases=bases,lambda =lam_min,family = "binomial" )
        gamsel.out = gamsel(mat_dx_train,Ytrain,degrees=degrees,dfs=dfs,lambda =lam_min,family = "binomial" )
        
        xb = predict(gamsel.out, mat_dx_test)
        beta_gam = get_beta_gamsel(gamsel.out)
        
        if(length(unique(xb))<2){
            auc_gam = 0.5 
        }else{
            auc_gam = ROC.Est.FUN(Ytest, xb, fpr0=0.05)[1]
        }
        
        res_list[[paste(iter_name,";gam",sep="")]] = list(auc=auc_gam, beta=beta_gam)
        
        
        ## run regression with PC scores
        
        W = rep(1,length(Ytrain))
        foldid = get_foldid(Ytrain)
        gvec = c(0, rep(1,p))
        ridge_gvec = c(0, rep(1,p))
        
        
        ## kde_fpca + count
        Xtrain = cbind(1,alpha_fpca_train)
        Xtest = cbind(1, alpha_fpca_test)
        #grpid = get_grpid(d_fpca+1, ncol(alpha_fpca) )
        grpid = grpid_fpca
        
        res = get_beta_auc_binom_wtest( W,  Xtrain, Ytrain,  Xtest,Ytest, 
                                        foldid,  grpid, lam_seq, 
                                        ridge_seq, ridge_gvec, gvec, 
                                        eps, maxiter, ncores)
        
        res_list[[paste(iter_name,";fpca",sep="")]] = list(auc=res$auc, beta=res$beta_norm)
        
        ## distribution (GPAM)
        
        Xtrain = cbind(1,alpha_x_distr_train)
        Xtest = cbind(1,alpha_x_distr_test)
        grpid = get_grpid(d_distr, ncol(alpha_x_distr) )
        
        res = get_beta_auc_binom_wtest( W,  Xtrain, Ytrain,  Xtest,Ytest, 
                                        foldid,  grpid, lam_seq, 
                                        ridge_seq, ridge_gvec, gvec, 
                                        eps, maxiter, ncores)
        
        res_list[[paste(iter_name,";distr",sep="")]] = list(auc=res$auc, beta=res$beta_norm)
        
    }
}



## results

sapply(res_list, function(x){x$auc})

filename = paste("resultSimu/", setting, "_", ffs, ".rda", sep="" )
save(res_list, file=filename)
quit(save="no")



