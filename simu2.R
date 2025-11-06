## simu 2 ###
# require(mable)
# rm(list=ls())
require(MCMCpack)
require(data.table)
require(parallel)
require(survival)
require(splines)
source("simuFun_simple.R")

phi_k <- function(t=NULL, params=NULL,Ti=1){
    
    prob = params[1:ncomp]
    shape1 = params[1:ncomp + (ncomp)]
    shape2 = params[1:ncomp + 2*(ncomp)]
    count = params[length(params)]
    
    res = sapply(seq_along(prob), function(ii){
        prob[ii] * dbeta(t,shape1[ii],shape2[ii])
    })
    res = matrix(res, ncol=length(prob))
    res = count * apply(res, 1, sum)
    res
}

myfun <- function(x=NULL, type=NULL){
    res = matrix(1, nrow=length(x),ncol=p0)
    coordinate = cbind(1:nrow(res),type)
    res[coordinate]
}


n = 1000
p0 = 4
p = 150

beta_density = 1*rep(c(1,-1),each=p0/2)

ncomp = 3
count_mat = matrix(runif(n*p,lam_lo,lam_up), ncol = p)

mixprob_list = lapply(1:n, function(i){
    rdirichlet(p,rep(5,ncomp))
})
beta_shape1_list = lapply(1:n, function(i){
    matrix(runif(p*ncomp, 1, 20),ncol=ncomp)
})
beta_shape2_list = lapply(1:n, function(i){
    matrix(runif(p*ncomp, 1, 20),ncol=ncomp)
})



### generate point process data
data.list = lapply(1:n, function(i){
    
    obstime = 1
    TT = c(0, obstime)
    
    eventlist = lapply(1:p, function(pp){
        
        prob_ip = mixprob_list[[i]][pp,] 
        b_shape1 = beta_shape1_list[[i]][pp,]
        b_shape2 = beta_shape2_list[[i]][pp,]
        
        param_i = c(prob_ip, b_shape1,b_shape2,count_mat[i,pp])
        
        model = mpp(data=NULL, gif = fexp_gif2, mark = list(NULL, NULL),
                    params= param_i, TT = TT, 
                    gmap=expression(params), mmap=NULL)
        
        res = NULL
        max.rate = optimize(bhz2, interval=c(0,obstime),
                            params=param_i, maximum=TRUE)$objective
        max.rate = 1.1 *max.rate
        
        while(any(class(res)=="try-error") | any(is.null(res)) ){
            res = try({simulate(model,  max.rate = max.rate)},silent = TRUE)
            max.rate = max.rate*2
            #cat("increased max rate!\n")
            #cat(xb , "\n")
        }
        #print(max.rate)
        df = NULL
        if(!is.null(res$data)){
            df=data.frame(feature_id=rep(pp,length(res$data$time)),
                          time=res$data$time)
        }
        df
    })
    
    df = do.call(rbind, eventlist)
    #list(eventlist=eventlist,X=Xi,obstime=obstime)
    if(!is.null(df)){
        df = data.table(ID=i, df)
    }else{
        df = data.table(ID=i, feature_id =NA, time=NA)
    }
    
    df
})

data_df = do.call(rbind, data.list)
data_df = data_df[!is.na(time), ]
colnames(data_df)[1] = 'subject_num'


## generate outcome data

data_df_related = data_df[is.element(feature_id, 1:p0),]
data_df_related[,value:= myfun(time, feature_id)]

xb_density = data_df_related[,sum(value),by=c("subject_num","feature_id")]

xb_density = sparseMatrix(i=xb_density[,subject_num], 
                          j=xb_density[,feature_id], x=xb_density[,V1],
                          dims = c(n,p0))

xb_density = as.matrix(xb_density)


mu = (lam_lo+lam_up)/2

xb_density[,1] = 0.2*(xb_density[,1]-mu)
xb_density[,2] = 2*as.numeric(xb_density[,2]>mu)
xb_density[,3] = 0.05*(xb_density[,3]-mu)^2
xb_density[,4] = 2*sin(xb_density[,4]*pi/mu)

# apply(xb_density,2,sd)
xb_density_tmp = xb_density

xb_density = xb_density %*% beta_density
intercept = 0
xb = exp(intercept+xb_density) /( 1+ exp(intercept+xb_density))
Y = sapply(1:n, function(ii){
    sample(0:1, 1, replace = TRUE, c(1-xb[ii], xb[ii]) )
})

fx_true = t(t(xb_density_tmp) * beta_density)



