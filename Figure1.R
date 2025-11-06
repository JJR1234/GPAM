#### Figure 1 ####

# id: 139157, 149613 #
rm(list=ls())
require(data.table)
load("data_mimic.rda")
set.seed(3000)
names_plot = c("NIMAP","Na","MechVent")
codes = c("Non-invasive Mean\n Arterial Blood Pressure\n Measurement",
          "Blood Sodium\n Measurement",
          "Mechanical\n Ventilation")
data_long[,Time:= time*48]
count_all = data_long[,.(count=.N), by=c("subject_num","Parameter")]
count_all = count_all[is.element(Parameter,names_plot),]
count_all = count_all[,.(count=.N), by=c("subject_num")]
count_all = count_all[count==length(codes),]
id = sample(count_all$subject_num,2)

filename ="resultMIMIC/Figure1.eps"


setEPS()
postscript(file=filename, width = 12, height = 4)
layout(matrix(1:length(id),ncol=2))
par(mar=c(5,10,3,3))

for(i in seq_along(id)){
    tmp = data_long[ subject_num==id[i], ]
    obs = 48
    plot(0,0,xlim=c(0,obs),ylim=c(0.5,length(codes)+1),cex=0.000001,
         bty='L', ylab = "", xlab ="Time", yaxt="n", xaxt="n",
         main=paste("Patient", i))
    at_lab = seq(0,48,by=6)
    axis(1, at=at_lab,labels=at_lab,  las=1)
    axis(2, at=1:length(codes),labels=codes,  las=2)
    
    segments(0,1,obs,1)
    for(xx in tmp$Time[tmp$Parameter==names_plot[[1]]]){
        segments(xx,1,xx,1.5)
    }
    for(j in 2:length(codes)){
        segments(0,j,obs,j)
        for(xx in tmp$Time[tmp$Parameter==names_plot[[j]] ]){
            segments(xx,0+j,xx,0.5+j)
        }
    }
}
dev.off()


