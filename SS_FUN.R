## Functions for additive kernel regression ##

require(Matrix)
require(Rcpp)
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

## calc gamma
get_gamma <- function(tseq=NULL){
    n = length(tseq)
    if(n>2& n<1000){
        res = 1/(sum(dist(tseq))*2/n/(n-1))^2
    }else if(n >= 1000){
        tseq = sample(tseq, 1000)
        n = length(tseq)
        res = 1/(sum(dist(tseq))*2/n/(n-1))^2
    }else{
        res = NaN
    }
    res
}

## canl gram for Time
sqfun <- function(x=NULL,y=NULL,gamma=NULL){
    exp(-gamma*(x-y)^2)
}

KTfun <- function(x=NULL, y=NULL, gamma=1){
    mean(outer(x,y,sqfun,gamma=gamma))
}

KT_gram <- function(xlist = NULL){
    
    gamma_list = sapply(xlist, function(xx){
        get_gamma(xx)
    })
    gamma_med = median(gamma_list, na.rm = TRUE)
    
    res = matrix(NA, length(xlist), length(xlist))
    for(i in seq_along(xlist)){
        for(j in i:length(xlist)){
            res[i,j] = KTfun(xlist[[i]], xlist[[j]], gamma_med)
        }
    }
    
    res[lower.tri(res)] = t(res)[lower.tri(res)]
    DXX = outer(diag(res), diag(res), "+")
    res = DXX - 2*res
    
    gamma_est = 0.5 * sum(sqrt(res)) / choose(nrow(res),2)
    gamma_est = 1/gamma_est^2
    exp(-gamma_est* res)
}


#### gram 
# X , X_test: NxD matrices
cal_gram <- function(X_test=NULL, X=NULL, KerType= 1, 
                     polyD=NULL, gamma=NULL){
    
    Lx=NULL
    gamma_est = NULL
    XX = X_test%*%t(X)
    
    if(KerType==1){
        # KerType = 1 -- Gaussian kernel
        # k(x,y)=exp(-gamma(x-y)^2)
        X_test_diag = diag(X_test%*%t(X_test))
        X_diag = diag(X%*%t(X))
        
        DXX = outer(X_test_diag, X_diag, "+")
        XX = DXX - 2*XX
        if(!is.null(gamma)){
            Lx = exp(-gamma* XX)
        }else{
            gamma_est = 0.5 * sum(sqrt(XX)) / choose(nrow(XX),2)
            gamma_est = 1/gamma_est^2
            Lx = exp(-gamma_est* XX)
        }
        
    }else if(KerType==2){
        # KerType = 2 -- polynomial kernel
        # k(x,y)=(1+xy)^d
        Lx = (1+XX)^polyD
    }else{
        ## other -- linear kernel
        Lx = XX
    }
    
    list(Lx=Lx, gamma_est=gamma_est)
}


### centering gram matrix

QXQ_FUN <- function(XX=NULL){
    cmean = apply(XX,2,mean)
    XX = t(XX - cmean)
    cmean = apply(XX,1,mean)
    XX- cmean
}


## get gram matrix

get_gram_count_vec <- function(data_all=NULL,patient_sel=NULL,
                                  feature_sel=NULL, ncores=3){
    
    data_all = data_all[is.element(feature_id, feature_sel),]
    
    res = lapply(feature_sel, function(i){
        
        data_i = data_all[feature_id==i,]
        count_value = data_i[,.(count=.N), by=c("subject_num")]
        
        patient_now = unique(count_value$subject_num)
        
        scores = rep(0, length(patient_sel) )
        
        scores[match(patient_now,patient_sel)] = count_value$count
        rm(count_value)
        
        ## no scaling?
        scores = matrix(scores,ncol=1)
        
        ## 
        gram = cal_gram(scores, scores)
        
        list(gram=gram)
    })
    res
}


get_gram_shape_vec <- function(data_all=NULL,patient_sel=NULL,
                               feature_sel=NULL, ncores=3){
    
    data_all = data_all[is.element(feature_id, feature_sel),]
    
    res = lapply(feature_sel, function(i){
        
        data_i = data_all[feature_id==i,]
        count_value = data_i[,.(count=.N), by=c("subject_num")]
        count_med = median(count_value$count)
        
        x_list = lapply(patient_sel, function(id){
            ind = which(data_i$subject_num == id)
            tt = NULL
            if(length(ind) > 0){
                tt =  data_i$time[ind]
            }else{
                tt = runif(50)
                # tt = runif(count_med)
            }
            tt
        })
        ##
        min_len = min(c(500,length(x_list)))
        gamma = get_gamma_list(x_list[1:min_len])
        gram =  get_KT_gram(x_list, gamma)
        #gram = KT_gram(x_list)
        
        list(gram=gram)
    })
    res
}



get_gram_distribution_vec <- function(data_all=NULL,patient_sel=NULL,
                                      feature_sel=NULL, ncores=3){
    
    data_all = data_all[is.element(feature_id, feature_sel),]
    
    res = lapply(feature_sel, function(i){
        
        data_i = data_all[feature_id==i,]
        count_value = data_i[,.(count=.N), by=c("subject_num")]
        count_med = median(count_value$count)
        
        x_list = lapply(patient_sel, function(id){
            ind = which(data_i$subject_num == id)
            tt = NA
            if(length(ind) > 0){
                tt =  data_i$time[ind]
            }
            tt
        })
        ##
        na_list = sapply(x_list, function(x){any(is.na(x))})
        min_len = min(c(500,sum(!na_list)))
        gamma = get_gamma_list(x_list[!na_list][1:min_len])
        gram =  get_KT_gram_sum(x_list, gamma)
        #gram = KT_gram(x_list)
        
        list(gram=gram)
    })
    res
}


get_KDE_PCscore <- function(data_all=NULL,patient_sel=NULL,
                            feature_sel=NULL, ncores=3){
    
    data_all = data_all[is.element(feature_id, feature_sel),]
    
    res = lapply(feature_sel, function(i){
        
        data_i = data_all[feature_id==i,]
        x_seq = seq(0,1,length=101)
        y_list = lapply(patient_sel, function(id){
            ind = which(data_i$subject_num == id)
            yy = NULL
            if(length(ind) > 0){
                tt = data_i$time[ind]
                if(length(tt)==1){
                    den = density(tt,from=0,to=1,n=101,bw=0.1)
                }else{
                    den = density(tt,from=0,to=1,n=101)
                }
                yy = den$y
            }else{
                yy = rep(1,101)
            }
            #cbind(id,x_seq, yy)
            yy
        })
        
        x_list = rep(list(x_seq), length(patient_sel))
        
        fitted_fpca = FPCA(Ly = y_list, Lt=x_list,optns=list(maxK=10))
        
        list(scores = fitted_fpca$xiEst,cumFVE = fitted_fpca$cumFVE)
    })
    
    res
}


### get pc scores

get_PC_score <- function(Lx=NULL,d=5, d_max=20, var_exp_threshold=0.9){
    # #EnK = apply(gram$Lx, 1, mean)
    # Fik = Lx - apply(Lx, 1, mean) ## C_n
    # EnC = apply(Fik, 2 , mean)
    # C = t(t(Fik) -  EnC) ## C_n - EnC
    # Gx = QXQ_FUN(Lx)/nrow(Lx)
    Gx = QXQ_FUN(Lx)
    pca_svd = eigs_sym(Gx,d_max)
    eig_val =  pca_svd$values
    eig_val[eig_val < 1e-16] = 0
    var_exp = cumsum(eig_val)/sum(diag(Gx))
    if(d<1){
        if(max(var_exp)<var_exp_threshold){
            d = d_max
            warning("d_max does not reach 95% Var\n")
        }else{
            d = which(var_exp>var_exp_threshold)[1] 
        } 
    }
    
    eig_vec = pca_svd$vectors[,1:d]
    t(t(eig_vec)*sqrt(eig_val)[1:d])
    
    ## eta_X = t(t(eig_vec)*sqrt(eig_val)[1:d])
    ## alpha_X = C %*% eta_X
    ## alpha_X
}


### CV validation for binomial

# get_binom_loglik <- function(W=NULL, X=NULL, Y=NULL, beta=NULL){
#     xb = X%*%beta
#     sum(W*(-log1p(exp(xb)) + Y*xb))
# }


get_foldid <- function(Y=NULL, nfold=5){
    yi = Y
    fd = rep(0,length(yi))
    fd[yi==0] = sample(nfold, sum(yi==0), replace=TRUE)
    fd[yi==1] = sample(nfold, sum(yi==1), replace=TRUE)
    fd
}


# negative loglik for binomial
get_binom_loglik <- function(W=NULL, xb=NULL, Y=NULL){
    - sum(W*(-log1p(exp(xb)) + Y*xb))
}

CV_binom_nestGrplasso <- function( W=NULL,  X=NULL, Y=NULL, foldid=NULL, 
                                   beta_init=NULL, grpid=NULL, lam_seq=NULL, 
                                   ridge_seq=NULL, ridge_gvec=NULL, gvec=NULL, 
                                   eps=1e-8, maxiter=100, ncores=3){
    
    # loglik_cv = mclapply(1:max(foldid), function(fid){
    #     Xtrain = X[foldid!=fid,]
    #     Ytrain = Y[foldid!=fid]
    #     Wtrain = W[foldid!=fid] 
    #     Wtrain = Wtrain/sum(Wtrain)
    #     
    #     Xtest = X[foldid==fid,]
    #     Ytest = Y[foldid==fid]
    #     Wtest = W[foldid==fid] 
    #     Wtest = Wtest/sum(Wtest)
    #     
    #     beta_res = BinomWLM_nestGrplasso_seq( Wtrain,  Xtrain,  Ytrain,
    #                                           beta_init, grpid, lam_seq, ridge_seq, 
    #                                           ridge_gvec, gvec, eps , maxiter)
    #     xb_test = Xtest %*% beta_res
    #     loglik_test = apply(xb_test,2,function(xb){
    #         get_binom_loglik(Wtest, xb, Ytest)
    #     })
    #     loglik_test
    # },
    # mc.cores = ncores)
    
    loglik_cv = lapply(1:max(foldid), function(fid){
        Xtrain = X[foldid!=fid,]
        Ytrain = Y[foldid!=fid]
        Wtrain = W[foldid!=fid] 
        Wtrain = Wtrain/sum(Wtrain)
        
        Xtest = X[foldid==fid,]
        Ytest = Y[foldid==fid]
        Wtest = W[foldid==fid] 
        Wtest = Wtest/sum(Wtest)
        
        beta_res = BinomWLM_nestGrplasso_seq( Wtrain,  Xtrain,  Ytrain,
                                              beta_init, grpid, lam_seq, ridge_seq, 
                                              ridge_gvec, gvec, eps , maxiter)
        xb_test = Xtest %*% beta_res
        loglik_test = apply(xb_test,2,function(xb){
            get_binom_loglik(Wtest, xb, Ytest)
        })
        loglik_test
    })
    loglik_cv = do.call(rbind, loglik_cv)
    
    error_mean = apply(loglik_cv, 2,  mean)
    error_sd = apply(loglik_cv, 2,  sd)/sqrt(nrow(loglik_cv))
    error_1se = min(error_mean) + error_sd[which.min(error_mean)]
    
    
    ind_1se = which(error_mean <= error_1se)[1]
    ind_min = which.min(error_mean)
    
    list(loglik_cv = loglik_cv, ind_1se=ind_1se, ind_min=ind_min)
}



CV_binom_Grplasso <- function( W=NULL,  X=NULL, Y=NULL, foldid=NULL, 
                                   beta_init=NULL, grpid=NULL, lam_seq=NULL, 
                                   ridge_seq=NULL, ridge_gvec=NULL, gvec=NULL, 
                                   eps=1e-8, maxiter=100, ncores=3){
    
    # loglik_cv = mclapply(1:max(foldid), function(fid){
    #     Xtrain = X[foldid!=fid,]
    #     Ytrain = Y[foldid!=fid]
    #     Wtrain = W[foldid!=fid] 
    #     Wtrain = Wtrain/sum(Wtrain)
    #     
    #     Xtest = X[foldid==fid,]
    #     Ytest = Y[foldid==fid]
    #     Wtest = W[foldid==fid] 
    #     Wtest = Wtest/sum(Wtest)
    #     
    #     beta_res = BinomWLM_nestGrplasso_seq( Wtrain,  Xtrain,  Ytrain,
    #                                           beta_init, grpid, lam_seq, ridge_seq, 
    #                                           ridge_gvec, gvec, eps , maxiter)
    #     xb_test = Xtest %*% beta_res
    #     loglik_test = apply(xb_test,2,function(xb){
    #         get_binom_loglik(Wtest, xb, Ytest)
    #     })
    #     loglik_test
    # },
    # mc.cores = ncores)
    
    loglik_cv = lapply(1:max(foldid), function(fid){
        Xtrain = X[foldid!=fid,]
        Ytrain = Y[foldid!=fid]
        Wtrain = W[foldid!=fid] 
        Wtrain = Wtrain/sum(Wtrain)
        
        Xtest = X[foldid==fid,]
        Ytest = Y[foldid==fid]
        Wtest = W[foldid==fid] 
        Wtest = Wtest/sum(Wtest)
        
        beta_res = BinomWLM_Grplasso_seq( Wtrain,  Xtrain,  Ytrain,
                                              beta_init, grpid, lam_seq, ridge_seq, 
                                              ridge_gvec, gvec, eps , maxiter)
        xb_test = Xtest %*% beta_res
        loglik_test = apply(xb_test,2,function(xb){
            get_binom_loglik(Wtest, xb, Ytest)
        })
        loglik_test
    })
    loglik_cv = do.call(rbind, loglik_cv)
    
    error_mean = apply(loglik_cv, 2,  mean)
    error_sd = apply(loglik_cv, 2,  sd)/sqrt(nrow(loglik_cv))
    error_1se = min(error_mean) + error_sd[which.min(error_mean)]
    
    
    ind_1se = which(error_mean <= error_1se)[1]
    ind_min = which.min(error_mean)
    
    list(loglik_cv = loglik_cv, ind_1se=ind_1se, ind_min=ind_min)
}


get_beta_auc_binom_wtest <- function( W=NULL,  X=NULL, Y=NULL, 
                                      Xtest=NULL,Ytest=NULL, 
                                      foldid=NULL,  grpid=NULL, lam_seq=NULL, 
                                      ridge_seq=NULL, ridge_gvec=NULL, gvec=NULL, 
                                      eps=1e-8, maxiter=100, ncores=1){
    beta_init = rep(0, ncol(X))
    
    beta_res = BinomWLM_Grplasso_seq( W/sum(W),  X,  Y, beta_init, grpid, lam_seq, ridge_seq, 
                                      ridge_gvec, gvec,eps , maxiter)
    
    stat_res = CV_binom_Grplasso(W,X, Y, foldid,  beta_init, grpid, lam_seq, ridge_seq, 
                                 ridge_gvec, gvec, eps, maxiter, ncores)
    
    beta = beta_res[,stat_res$ind_1se,drop=FALSE]
    beta_norm  = calc_normBeta(beta, grpid)
    
    xb =  Xtest %*% beta
    # loglik_count = get_binom_loglik(rep(1,length(Ytest))/length(Ytest), xb, Ytest )
    if(length(unique(xb))<2){
        auc = 0.5 
    }else{
        auc = ROC.Est.FUN(Ytest, xb, fpr0=0.05)[1]
    }
    
    list(auc=auc, beta=beta, beta_norm=beta_norm, beta_all=beta_res)
    
}


get_beta_auc_binom_CV <- function( W=NULL,  X=NULL, Y=NULL,
                                   foldid=NULL,  grpid=NULL, lam_seq=NULL,
                                   ridge_seq=NULL, ridge_gvec=NULL, gvec=NULL, 
                                   eps=1e-8, maxiter=100, ncores=1){
    
    beta_init = rep(0, ncol(X))
    beta_res = BinomWLM_Grplasso_seq( W/sum(W),  X,  Y, beta_init, grpid, lam_seq, ridge_seq, 
                                      ridge_gvec, gvec,eps , maxiter)
    
    stat_res = CV_binom_Grplasso(W,X, Y, foldid,  beta_init, grpid, lam_seq, ridge_seq, 
                                 ridge_gvec, gvec, eps, maxiter, ncores)
    
    beta = beta_res[,stat_res$ind_1se,drop=FALSE]
    beta_norm  = calc_normBeta(beta, grpid)
    
    ## cv
    score.all = NULL
    y.all = NULL
    for(foldi in 1:max(foldid)){
        ## prepare x and y ##
        cat("prepare x and y at fold =",foldi,'\n')
        x = as.matrix(X[foldid!=foldi,]) 
        y = Y[foldid!=foldi]
        w = W[foldid!=foldi]
        
        beta_res_cv = BinomWLM_Grplasso_seq( w/sum(w),  x,  y, beta_init, grpid, 
                                          lam_seq, ridge_seq, 
                                          ridge_gvec, gvec,eps , maxiter)
        
        foldid_sub = get_foldid(y)
        stat_res = CV_binom_Grplasso(w, x, y, foldid_sub,  beta_init, grpid, 
                                     lam_seq, ridge_seq, 
                                     ridge_gvec, gvec, eps, maxiter, ncores)
        
        beta_cv = beta_res_cv[,stat_res$ind_1se,drop=FALSE]
        
        ## compute scores ## 
        x = as.matrix(X[foldid==foldi,]) 
        y = Y[foldid==foldi]
        score = as.numeric(x%*%beta_cv)
        score.all = c(score.all, score)
        y.all = c(y.all, y)
    }
    
    auc = ROC.Est.FUN(y.all, score.all, fpr0=0.05)[c(1,4:7)]
    names(auc) = c("auc","spec","sens","ppv","npv")
    
    list(auc=auc[1], auc_all=auc, beta=beta, beta_norm=beta_norm, beta_all=beta_res)
    
}

### help function for gamsel

get_df_gamsel <- function(X=NULL){
    uni = apply(X,2,function(x){
        length(unique(x))
    })
    dfs = rep(5, ncol(X))
    dfs[uni <= 5] = 1
    dfs
}

get_beta_gamsel <- function(gamselfit = NULL){
    alphas = abs(gamselfit$alphas)
    df_idx = cumsum(gamselfit$degrees)
    df_idx = cbind(c(1,df_idx[-length(df_idx)]+1), df_idx)
    betas = sapply(1:nrow(df_idx), function(ii){
        sqrt(mean(gamselfit$betas[df_idx[ii,1]:df_idx[ii,2]]^2))
    })
    c(gamselfit$intercept, pmax(alphas, betas))
}




