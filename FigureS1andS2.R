## run simu: to check estimated f(x) ##

rm(list=ls())
ffs = 1
set.seed(ffs)

require(Rcpp)
source("SS_FUN.R")
require(penReg)
source("ROCFUN.R")

setting = "Model_fx"

lam_lo = 5
lam_up = 15

ncores = 1

source("simu1.R")

patient_sel = 1:n
feature_sel = 1:p


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
    grpid = cbind(seq(0, p-perG, by=perG), seq(perG-1, p-1, by =perG))
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


# distribution 

Xtrain = cbind(1,alpha_x_distr_train)
Xtest = cbind(1,alpha_x_distr_test)
grpid = get_grpid(d_distr, ncol(alpha_x_distr) )

res = get_beta_auc_binom_wtest( W,  Xtrain, Ytrain,  Xtest,Ytest, 
                                foldid,  grpid, lam_seq, 
                                ridge_seq, ridge_gvec, gvec, 
                                eps, maxiter, ncores)

beta = res$beta

p0 = 4 # informative markers

fx_est_pen = sapply(2:(p0+1), function(ii){
    col1 = grpid[ii,1]+1
    col2 = grpid[ii,2]+1
    Xtest[,col1:col2] %*% beta[col1:col2]
})


res = get_beta_auc_binom_wtest( W,  Xtrain[,1:(1+d_distr*p0)], Ytrain,  
                                Xtest[,1:(1+d_distr*p0)],Ytest, 
                                foldid,  grpid[1:(p0+1),], c(0,0), 
                                c(0,0), ridge_gvec, gvec, 
                                eps, maxiter, ncores)

beta = res$beta
fx_est_unpen = sapply(2:(p0+1), function(ii){
    col1 = grpid[ii,1]+1
    col2 = grpid[ii,2]+1
    Xtest[,col1:col2] %*% beta[col1:col2]
})


fx_true = fx_true[test_id,]


#### make plots 
require(ggplot2)

##Figure S1
plot_df = data.frame(fx_est = as.numeric(fx_est_pen),
                     fx_true = as.numeric(scale(fx_true,scale=F) ),
                     type = rep(paste("f[",1:4,"]","(", "X[",1:4,"]"  ,")",sep=""),
                                each=nrow(fx_est_pen))
)

ggplot(data=plot_df, aes(x=fx_true, y=fx_est)) +
    geom_point()+
    geom_abline() +
    xlab("True value") + ylab("Estimated value") +
    theme(plot.title = element_text(hjust = 0.5,size=16),
          # legend.position.inside = c(0.8, 0.8),
          legend.title = element_text(size=15),
          legend.text = element_text(size=15),
          legend.key.size = unit(1, 'cm'),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))  + 
    facet_wrap(.~type,nrow=2,labeller = label_parsed)

filename = "FigureS1.eps"
ggsave(file=filename, width = 6, height = 6)



## Figure S2
plot_df = data.frame(fx_est = as.numeric(fx_est_unpen),
                     fx_true = as.numeric(scale(fx_true,scale=F) ),
                     type = rep(paste("f[",1:4,"]","(", "X[",1:4,"]"  ,")",sep=""),
                                each=nrow(fx_est_pen))
)

ggplot(data=plot_df, aes(x=fx_true, y=fx_est)) +
    geom_point()+
    geom_abline() +
    xlab("True value") + ylab("Estimated value") +
    theme(plot.title = element_text(hjust = 0.5,size=16),
          # legend.position.inside = c(0.8, 0.8),
          legend.title = element_text(size=15),
          legend.text = element_text(size=15),
          legend.key.size = unit(1, 'cm'),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))  + 
    facet_wrap(.~type,nrow=2,labeller = label_parsed)

filename = "FigureS2.eps"
ggsave(file=filename, width = 6, height = 6)

