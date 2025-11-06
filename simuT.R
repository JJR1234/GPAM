## simu for T's effect ###
# require(mable)
# rm(list=ls())
require(MCMCpack)
require(data.table)
require(parallel)
require(survival)
require(splines)

myfun <- function(x=NULL, type=NULL){
    res = matrix(1, nrow=length(x), ncol=p0)
    coordinate = cbind(1:nrow(res), type)
    #apply(res * t(type_arr),1,sum)
    res[coordinate]
}

n = 1000
p0 = 4
p = 50
beta_density = 1*rep(c(1,-1),each=p0/2)

count_mat = matrix(runif(n*p,lam_lo,lam_up), ncol = p)*Tmax

### generate point process data
data.list = lapply(1:n, function(i){
    
    eventlist = lapply(1:p, function(pp){
        
        mu = count_mat[i,pp]
        count = rpois(1, mu)
        if(count!=0){
            times = sort(runif(count, 0, Tmax))
        }
        df = NULL
        if(count!=0){
            df = data.frame(feature_id=rep(pp,length(times)),
                          time=times)
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
data_df[,time:=time/Tmax]

## generate outcome data

data_df_related = data_df[is.element(feature_id, 1:p0),]
data_df_related[,value:= myfun(time, feature_id)]

xb_density = data_df_related[,sum(value),by=c("subject_num","feature_id")]

xb_density = sparseMatrix(i=xb_density[,subject_num], 
                          j=xb_density[,feature_id], x=xb_density[,V1],
                          dims = c(n,p0))

xb_density = as.matrix(xb_density)
xb_density = scale(xb_density)


xb_density[,1] = 1*(xb_density[,1])
xb_density[,2] = 2*as.numeric(xb_density[,2]>0)
xb_density[,3] = 0.8*(xb_density[,3]-0)^2
xb_density[,4] = 1.5*sin(xb_density[,4]*pi)

apply(xb_density,2,sd)

xb_density_tmp = xb_density

xb_density = xb_density %*% beta_density
intercept = 0
xb = exp(intercept+xb_density) /( 1+ exp(intercept+xb_density))
Y = sapply(1:n, function(ii){
    sample(0:1, 1, replace = TRUE, c(1-xb[ii], xb[ii]) )
})

fx_true = t(t(xb_density_tmp) * beta_density)



