#### Figure 2#####

rm(list=ls())
require(data.table)
load("resultMIMIC/mimic_res.rda")

name_it = "MechVent" # MechVent,NIDiasABP,AST
idx_it = which(codes==name_it)
d_distri = 10
xx = alpha_x_distr1[,1:d_distri + (idx_it-1)*d_distri]
bb = res_list[["GPAMgg"]]$beta[-(1:4),]
xb = xx %*% bb[1:d_distri + (idx_it-1)*d_distri]

load("data_mimic.rda")
data_long[,Time:=time*48]
# boxplot(xb~data_base$`In-hospital_death`)

count_all = data_long[,.(count=.N), by=c("subject_num","Parameter")]
count_all = count_all[Parameter==name_it,]

matid = match(data_base$subject_num, count_all$subject_num)
count_all = count_all[matid,]
count_all$subject_num = data_base$subject_num

count_all[,count:=ifelse(is.na(count),0,count)]
count_all$Y = data_base$`In-hospital_death`
count_all$xb = xb
# boxplot(count_all$count~count_all$Y)

count_all = count_all[count>5,]
count_all0 = count_all[count_all$Y==0,]
count_all1 = count_all[count_all$Y==1,]
setorder(count_all0, xb)
setorder(count_all1, -xb)

ID1 = count_all1$subject_num[1:50]
ID0 = count_all0$subject_num[1:50]

tmp = data_long[Parameter==name_it,]
dens1 = lapply(ID1, function(id){
    Ntmp = length(tmp$Time[tmp$subject_num==id])
    dens = density(tmp$Time[tmp$subject_num==id],from=0,to=48)
    dens$y = Ntmp*dens$y
    dens
})

dens0 = lapply(ID0, function(id){
    Ntmp = length(tmp$Time[tmp$subject_num==id])
    dens = density(tmp$Time[tmp$subject_num==id],from=0,to=48)
    dens$y = Ntmp*dens$y
    dens
})


## ggplot
require(ggplot2)
dens1_df = lapply(seq_along(dens1), function(ii){
    df = dens1[[ii]]
    data.frame(x=df$x, y=df$y,ID=ii,Outcome="Death")
})
dens1_df = do.call(rbind, dens1_df)

dens0_df = lapply(seq_along(dens0), function(ii){
    df = dens0[[ii]]
    data.frame(x=df$x, y=df$y,ID=ii+100,Outcome="Survive")
})
dens0_df = do.call(rbind, dens0_df)

dens_df = rbind(dens0_df, dens1_df)

ggplot(data=dens_df, aes(x=x, y=y, group=ID)) +
    geom_line(aes(linetype=Outcome,color=Outcome))+xlim(c(0,48)) +
    xlab("Time") + ylab("Intensity") +
    theme(plot.title = element_text(hjust = 0.5,size=16),
          legend.position.inside = c(0.8, 0.8),
          legend.title = element_text(size=15),
          legend.text = element_text(size=15),
          legend.key.size = unit(1, 'cm'),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))  

filename = paste("resultMIMIC/","Figure2.eps",sep='')
ggsave(file=filename, width = 8, height = 4)

