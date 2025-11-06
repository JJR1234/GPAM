rm(list=ls())
library(data.table)

#######################################
### modify path_to_data accordingly ###
#######################################
path_to_data = "../../MIMIC data/"
cov_base = c('RecordID','Age','Gender','Height','ICUType','Weight')
codes = NULL
data_base_list = list()
data_long_list = list()

for(name in c("a","b","c")){
    print(name)
    # List all text files in the the folder
    files = list.files(path = paste(path_to_data,"set-",name,"/",sep=""), 
                       pattern = "\\.txt$", full.names = TRUE)
    
    # Read and combine all files 
    combined_dt = rbindlist(
        lapply(files, function(f) {
            id = as.integer(sub("\\.txt$", "", basename(f)))
            dt = fread(f)
            dt[, ID := id]
            dt
        }),
        use.names = TRUE, fill = TRUE, idcol = NULL
    )
    
    # reorder columns (ID first)
    setcolorder(combined_dt, c("ID", setdiff(names(combined_dt), "ID")))
    
    # normalize [0,48] hours to [0,1]
    combined_dt[, time := {
        parts <- tstrsplit(Time, ":", type.convert = TRUE)
        (parts[[1]] * 60 + parts[[2]]) / (48 * 60)
    }]
    
    data_long = combined_dt[!is.element(Parameter, cov_base),]
    data_long = data_long[Parameter!="",]
    if(is.null(codes)){
        codes = unique(data_long$Parameter)
    }
    data_long[,feature_id := match(Parameter,codes) ]
    
    data_base = combined_dt[is.element(Parameter, cov_base),]
    data_base = data_base[time==0,]
    # multiple weight measurements at baseline for a few subjects
    # use mean value
    data_base = dcast(data_base, ID ~ Parameter, value.var = "Value",
                      fun.aggregate = mean)
    
    # replace the missing values ("-1") by sample mean
    cols = c("Gender", "Height", "Weight")
    for (j in cols) {
        mean_val = mean(data_base[get(j) != -1, get(j)], na.rm = TRUE)
        data_base[get(j) == -1, (j) := mean_val]
    }
    
    data_base[,BMI := Weight/(Height/100)^2]
    
    # merge with outcome data
    outcome = fread(paste(path_to_data,'Outcomes-',name,".txt",sep=""))
    
    data_base = merge(data_base, outcome, by="RecordID", all.x = TRUE)
    
    colnames(data_base)[1] = "subject_num"
    colnames(data_long)[1] = "subject_num"
    
    data_base_list[[name]] = data_base
    data_long_list[[name]] = data_long
}

data_base = do.call(rbind, data_base_list)
data_long = do.call(rbind, data_long_list)

save(data_long, data_base, codes, file='data_mimic.rda')


