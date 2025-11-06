## mimic data ##
rm(list=ls())
require(data.table)
require(Matrix)
require(glmnet)
require(gamsel)
require(Rcpp)
require(penReg)
source("SS_FUN.R")
source("ROCFUN.R")

get_grpid <- function(perG=5,p=100){
    G = p/perG
    grpid = cbind(seq(0, p-perG, by =perG), seq(perG-1, p-1, by=perG))
    rbind(c(0,0), grpid+1)
}

eps = 1e-8
maxiter = 100

load("data_mimic.rda")
patient_sel = data_base$subject_num
feature_sel = seq_along(codes)

d_distr = 10
d_fpca = 5

## get PC scores ##
cat("get PC scores.... \n ")
alpha_x_distr1_list = list()
alpha_x_distr2_list = list()
alpha_x_distr3_list = list()
alpha_x_distr4_list = list()
kde_score_list = list()

for(ffs in feature_sel){
    cat("get PC scores for ", ffs, "th point process ..\n")
    gram_distr = get_gram_distribution_vec(data_all=data_long,
                                            patient_sel=patient_sel,
                                            feature_sel=ffs, 
                                            ncores=ncores,
                                            Gaus_1=TRUE, Gaus_2=TRUE)
    alpha_x_distr1 = get_PC_score(gram_distr[[1]]$gram,d=d_distr)
    
    gram_distr = get_gram_distribution_vec(data_all=data_long,
                                            patient_sel=patient_sel,
                                            feature_sel=ffs, 
                                            ncores=ncores,
                                            Gaus_1=FALSE, Gaus_2=TRUE)
    
    alpha_x_distr2 = get_PC_score(gram_distr[[1]]$gram,d=d_distr)
    
    gram_distr = get_gram_distribution_vec(data_all=data_long,
                                            patient_sel=patient_sel,
                                            feature_sel=ffs, 
                                            ncores=ncores,
                                            Gaus_1=TRUE, Gaus_2=FALSE)
    alpha_x_distr3 = get_PC_score(gram_distr[[1]]$gram,d=d_distr)
    
    gram_distr = get_gram_distribution_vec(data_all=data_long,
                                            patient_sel=patient_sel,
                                            feature_sel=ffs, 
                                            ncores=ncores,
                                            Gaus_1=FALSE, Gaus_2=FALSE)
    alpha_x_distr4 = get_PC_score(gram_distr[[1]]$gram,d=d_distr)
    
    kde_score = get_KDE_PCscore(data_all=data_long, patient_sel=patient_sel,
                                feature_sel=ffs, ncores=ncores)
    
    alpha_x_distr1_list[[ffs]] = alpha_x_distr1[,1:d_distr]
    alpha_x_distr2_list[[ffs]] = alpha_x_distr2[,1:d_distr]
    alpha_x_distr3_list[[ffs]] = alpha_x_distr3[,1:d_distr]
    alpha_x_distr4_list[[ffs]] = alpha_x_distr4[,1:d_distr]
    kde_score_list = c(kde_score_list, kde_score)
    
}

alpha_x_distr1 = alpha_x_distr1_list
alpha_x_distr2 = alpha_x_distr2_list
alpha_x_distr3 = alpha_x_distr3_list
alpha_x_distr4 = alpha_x_distr4_list
kde_score = kde_score_list


#### scores ####
alpha_x_distr1 = do.call(cbind, alpha_x_distr1)
alpha_x_distr2 = do.call(cbind, alpha_x_distr2)
alpha_x_distr3 = do.call(cbind, alpha_x_distr3)
alpha_x_distr4 = do.call(cbind, alpha_x_distr4)

##### count data ######
count_value = data_long[,.(count=.N), by=c("subject_num","feature_id")]
count_value[, i:= match(subject_num, patient_sel) ]
count_value[, j:= match(feature_id, feature_sel) ]

mat_dx = sparseMatrix(i=count_value[,i], j=count_value[,j], x=count_value[,count],
                      dims = c(length(patient_sel), length(feature_sel)))
colnames(mat_dx) = codes
mat_dx = as.matrix(mat_dx)

## FPCA data ##

alpha_fpca = lapply(seq_along(feature_sel), function(j){
    if(ncol(kde_score[[j]]$score) >= d_fpca){
        scale(cbind(kde_score[[j]]$score[,1:d_fpca], mat_dx[,j]) )
    }else{
        scale(cbind(kde_score[[j]]$score, mat_dx[,j]) )
    }
    
})

grpid_fpca = sapply(alpha_fpca, function(x){
    ncol(x)
})
grpid_fpca = c(1,rep(1,3), grpid_fpca)
grpid_fpca = cbind(c(0, cumsum(grpid_fpca)[-length(grpid_fpca)]), 
                   cumsum(grpid_fpca)-1)

alpha_fpca = do.call(cbind, alpha_fpca)

## run regression ##
cat("run regressions .... \n ")

Y = data_base$`In-hospital_death`
set.seed(1)
foldid = get_foldid(Y)

## run cross-validation

auc_all_list = list()
beta_list = list()
res_list = list()

## glmnet on count
X = scale(cbind(data_base[,c("Age","Gender","BMI")], mat_dx))

##
x = as.matrix(X) 
y = Y

res = cv.glmnet(x=x, y=y, family='binomial',  alpha=0.99,
                lambda.min.ratio=0.001,  nfolds=5)
beta_glmnet = coef(res)[,1]

## cv error
score.all = NULL
y.all = NULL
for(foldi in 1:5){
    ## prepare x and y ##
    cat("prepare x and y at fold =",foldi,'\n')
    x = as.matrix(X[foldid!=foldi,]) 
    y = Y[foldid!=foldi]
    
    res = cv.glmnet(x=x, y=y, family='binomial',  alpha=0.99,
                    lambda.min.ratio=0.001,  nfolds=5)
    beta = coef(res)[,1]
    
    ## compute scores ## 
    x = as.matrix(X[foldid==foldi,]) 
    y = Y[foldid==foldi]
    
    score = as.numeric(cbind(1,x)%*%beta)
    score.all = c(score.all, score)
    y.all = c(y.all, y)
}

auc = ROC.Est.FUN(y.all, score.all, fpr0=0.05)[c(1,4:7)]
names(auc) = c("auc","spec","sens","ppv","npv")


auc_all_list[["GLMnet"]] = auc
beta_list[["GLMnet"]] = beta_glmnet

## GAM on count
X = scale(cbind(data_base[,c("Age","Gender","BMI")], mat_dx))

degrees = get_df_gamsel(X)
dfs =  degrees - 2
dfs[dfs<1] = 1

##
x = as.matrix(X) 
y = Y

gamsel.cv = cv.gamsel(x,y,degrees=degrees,dfs=dfs,family = "binomial")
lam_min = gamsel.cv$lambda[gamsel.cv$index.1se]
gamsel.out = gamsel(x,y,degrees=degrees,dfs=dfs,lambda =lam_min,family = "binomial" )

beta_gam = get_beta_gamsel(gamsel.out)
beta_norm_gam = beta_gam

## cv error
score.all = NULL
y.all = NULL
for(foldi in 1:5){
    ## prepare x and y ##
    cat("prepare x and y at fold =",foldi,'\n')
    x = as.matrix(X[foldid!=foldi,]) 
    y = Y[foldid!=foldi]
    
    gamsel.cv = cv.gamsel(x,y,degrees=degrees,dfs=dfs,family = "binomial")
    lam_min = gamsel.cv$lambda[gamsel.cv$index.1se]
    gamsel.out = gamsel(x,y,degrees=degrees,dfs=dfs,lambda=lam_min,family = "binomial")
    
    ## compute scores ## 
    x = as.matrix(X[foldid==foldi,]) 
    y = Y[foldid==foldi]
    
    score = predict(gamsel.out, x)[,1]
    score.all = c(score.all, score)
    y.all = c(y.all, y)
}

auc = ROC.Est.FUN(y.all, score.all, fpr0=0.05)[c(1,4:7)]
names(auc) = c("auc","spec","sens","ppv","npv")

auc_all_list[["GAM"]] = auc
beta_list[["GAM"]] = beta_norm_gam


## kde_fpca + count
ZZ = scale(data_base[,c("Age","Gender","BMI")])
W = rep(1, nrow(alpha_fpca))

gvec = c(0,rep(1,ncol(ZZ)) ,rep(1,length(codes)))
ridge_gvec = c(0, rep(1,ncol(ZZ)), rep(1,length(codes)))

lam_seq = exp(seq(log(1),log(0.0001),length=100))
ridge_seq = lam_seq/100

X = as.matrix(cbind(1,ZZ,alpha_fpca))
grpid = grpid_fpca

## 
res = get_beta_auc_binom_CV( W,  X, Y, foldid,  grpid, lam_seq,
                             ridge_seq, ridge_gvec, gvec, 
                             eps, maxiter, ncores)
auc_all_list[["GFLM"]] = res$auc_all
beta_list[["GFLM"]] = res$beta_norm
res_list[["GFLM"]] = res


### GAM on distribution 

## GPAMgg
X = cbind(1,ZZ, alpha_x_distr1)
grpid = get_grpid(d_distr, ncol(alpha_x_distr1))
grpid = rbind(matrix(0:ncol(ZZ),nrow=ncol(ZZ)+1,ncol = 2), grpid[-1,]+ncol(ZZ))


res = get_beta_auc_binom_CV( W,  X, Y, foldid,  grpid, lam_seq,
                             ridge_seq, ridge_gvec, gvec, 
                             eps, maxiter, ncores)
auc_all_list[["GPAMgg"]] = res$auc_all
beta_list[["GPAMgg"]] = res$beta_norm
res_list[["GPAMgg"]] = res

## GPAMlg
X = cbind(1,ZZ, alpha_x_distr2)
grpid = get_grpid(d_distr, ncol(alpha_x_distr2))
grpid = rbind(matrix(0:ncol(ZZ),nrow=ncol(ZZ)+1,ncol = 2), grpid[-1,]+ncol(ZZ))

res = get_beta_auc_binom_CV( W,  X, Y, foldid,  grpid, lam_seq,
                             ridge_seq, ridge_gvec, gvec, 
                             eps, maxiter, ncores)
auc_all_list[["GPAMlg"]] = res$auc_all
beta_list[["GPAMlg"]] = res$beta_norm
res_list[["GPAMlg"]] = res

## GPAMgl
X = cbind(1,ZZ, alpha_x_distr3)
grpid = get_grpid(d_distr, ncol(alpha_x_distr3))
grpid = rbind(matrix(0:ncol(ZZ),nrow=ncol(ZZ)+1,ncol = 2), grpid[-1,]+ncol(ZZ))

res = get_beta_auc_binom_CV( W,  X, Y, foldid,  grpid, lam_seq,
                             ridge_seq, ridge_gvec, gvec, 
                             eps, maxiter, ncores)
auc_all_list[["GPAMgl"]] = res$auc_all
beta_list[["GPAMgl"]] = res$beta_norm
res_list[["GPAMgl"]] = res

## GPAMll
X = cbind(1,ZZ, alpha_x_distr4)
grpid = get_grpid(d_distr, ncol(alpha_x_distr4))
grpid = rbind(matrix(0:ncol(ZZ),nrow=ncol(ZZ)+1,ncol = 2), grpid[-1,]+ncol(ZZ))

res = get_beta_auc_binom_CV( W,  X, Y, foldid,  grpid, lam_seq,
                             ridge_seq, ridge_gvec, gvec, 
                             eps, maxiter, ncores)
auc_all_list[["GPAMll"]] = res$auc_all
beta_list[["GPAMll"]] = res$beta_norm
res_list[["GPAMll"]] = res


#### save results

folder_path <- "resultMIMIC/"

# Check if the folder exists
if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
    cat("Folder created:", folder_path, "\n")
} else {
    cat("Folder already exists:", folder_path, "\n")
}

print(auc_all_list)

AUC = as.data.frame(do.call(rbind, auc_all_list))
AUC = cbind(Method = rownames(AUC), AUC)

write.table(AUC, file="resultMIMIC/AUC.csv",
            row.names = F, sep=",",quote=F)

beta = data.table(codes,beta=beta_list[["GPAMgg"]][-(1:4)])
setorder(beta, -beta)

write.table(beta, file="resultMIMIC/beta.csv",
            row.names = F, sep=",",quote=F)

save(auc_all_list, beta_list, res_list, 
     codes, alpha_x_distr1, 
     file="resultMIMIC/mimic_res.rda")

# quit(save="no")






