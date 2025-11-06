## get results in Table 3 ##

rm(list = ls())

# Define the folder path
folder_path <- "resultSimu/"

# Check if the folder exists
if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
    cat("Folder created:", folder_path, "\n")
} else {
    cat("Folder already exists:", folder_path, "\n")
}

# run Model 3 for 100 replications
# this could take time
# parallelizing it is recommended
for (iter in 1:100) {
    # Build the command
    cmd = paste("Rscript run_Model3.R", iter)
    cat("Running:", cmd, "\n")
    # Run the command
    system(cmd)
}

# summarize results into Table 3

require(xtable)

get_TPFP <- function(stat=NULL,p0=4){
    ind = which(stat!=0)
    TP = length(intersect(ind, 1:p0))
    FP = length(setdiff(ind, 1:p0))
    c(TP, FP)
}


settings = c("Model3")
p0 = 4
#p = 100

res_all_list = list()

for(setting in settings){
    print(setting)
    res_ffs = list()
    
    for(ffs in 1:100){
        errorind = 0
        filename = paste(folder_path, setting, "_", ffs, ".rda", sep="" )
        if(file.exists(filename)){
            #load(filename)
            tryCatch( { load(filename) },
                      error = function(e){
                          errorind = 1
                          cat(filename, "no exists\n")
                      })
            if(errorind==1){next}
            # tmp = matrix(res_list$gam$beta,nrow=3,byrow = FALSE)
            # res_list$gam$beta = c(0,apply(tmp,2,function(x){mean(x^2)}))
            beta_norm = sapply(res_list, function(x){x$beta[-1]})
            auc = sapply(res_list, function(x){x$auc})
            TPFP_res = sapply(beta_norm,function(stat){
                get_TPFP(stat,p0)
            })
            p = strsplit(colnames(TPFP_res),";")
            p = sapply(p, function(x){x[2]})
            p = as.numeric(gsub("\\D", "", p))
            TPFP_res = TPFP_res/rbind(p0,p-p0)
            
            res_ffs[[ffs]] = rbind(TPFP_res, auc)
            rownames(res_ffs[[ffs]]) = c("TPR","FPR","AUC")
        }else{
            cat(filename, "no exists\n")
        }
    }
    res_all_list[[setting]] = res_ffs
}

res_df_mean = lapply(res_all_list, function(x){
    Reduce("+",x)/length(x)
})

res_df_sd = lapply(res_all_list, function(x){
    Ex = Reduce("+",x)/length(x)
    Ex2 = Reduce("+",lapply(x, function(xx){xx^2}))/length(x)
    sqrt(Ex2 - Ex^2)
})

ndigit = 3

res_df = lapply(seq_along(res_df_mean), function(i){
    mu = round(as.numeric(res_df_mean[[i]]),ndigit)
    sd = round(as.numeric(res_df_sd[[i]]),ndigit)
    xx = paste(format(mu,nsmall=ndigit), " (",
               format(sd, nsmall=ndigit),
               ")",sep="")
    xxx = matrix(xx, nrow=nrow(res_df_mean[[i]]))
    rownames(xxx) = rownames(res_df_mean[[i]])
    colnames(xxx) = colnames(res_df_mean[[i]])
    xxx
})


Methods = c("GLMnet","GAM","GFLM","GPAM")

res_xtable = lapply(res_df, function(x){
    xx = t(x)
    rnames = strsplit(rownames(xx),";")
    n = sapply(rnames, function(xxx){
        as.numeric(gsub("\\D", "", xxx[1]))
    })
    p = sapply(rnames, function(xxx){
        as.numeric(gsub("\\D", "", xxx[2]))
    })
    
    xxx = data.frame(Method=Methods,n=n,p=p,xx)
    
    rownames(xxx)=NULL
    xxx
})

res_xtable = xtable(do.call(rbind, res_xtable),
                    display=c("s","s","d","d","s","s","s"))

print(res_xtable,include.rownames=FALSE)

print(res_xtable, file=paste(folder_path,"Table3.txt",sep=""),include.rownames=FALSE)

