sum.I <- function(yy,FUN,Yi,Vi=NULL)
{
    if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}
    pos <- rank(c(yy,Yi),ties.method='f')[1:length(yy)]-rank(yy,ties.method='f')
    if (substring(FUN,2,2)=="=") pos <- length(Yi)-pos
    if (!is.null(Vi)) {
        if(substring(FUN,2,2)=="=") tmpind <- order(-Yi) else  tmpind <- order(Yi)
        ##Vi <- cumsum2(as.matrix(Vi)[tmpind,])
        Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)
        return(rbind(0,Vi)[pos+1,])
    } else return(pos)
}

####sum.I.old######
sum.I.old <- function(yy,FUN,Yi,Vi=NULL)   ## sum I(yy FUN Yi)Vi
{  
    if(FUN=="<"|FUN=="<=") { yy <- -yy; Yi <- -Yi}
    if(substring(FUN,2,2)=="=") yy <- yy + 1e-8 else yy <- yy - 1e-8
    pos <- rank(c(yy,Yi))[1:length(yy)] - rank(yy)
    if(is.null(Vi)){return(pos)}else{
        Vi <- cumsum2(as.matrix(Vi)[order(Yi),,drop=F])
        out <- matrix(0, nrow=length(yy), ncol=dim(as.matrix(Vi))[2])
        out[pos!=0,] <- Vi[pos,]
        if(is.null(dim(Vi))) out <- c(out)
        return(out) ## n.y x p
    }
} 




####S.FUN######
S.FUN <- function(yy,Yi,Di,yes.smooth=F,yes.half=T)
{
    if(yes.smooth){
        Y1i = Yi[Di==1]; n1 = sum(Di); bw = bw.nrd(Y1i)/n1^0.6
        out=c(t(rep(1/n1,n1))%*%pnorm((Y1i-VTM(yy,n1))/bw))
    }else{
        out = sum.I(yy,"<=",Yi,Vi=Di)/sum(Di)
        if(yes.half){out=(sum.I(yy,"<",Yi,Vi=Di)/sum(Di)+out)/2}
    }
    out
}

#####Sinv.FUN#####
Sinv.FUN <- function(uu,Yi,Di,yes.smooth=F)
{
    yy0<-unique(sort(Yi,decreasing=T)); ss0 <- S.FUN(yy0,Yi,Di,yes.smooth=yes.smooth) 
    return(approx(ss0[!duplicated(ss0)],yy0[!duplicated(ss0)],uu,method="linear",f=0,rule=2)$y)
}

#####ROC.Est.FUN#####
ROC.Est.FUN <- function(Di,yyi,yy0=NULL,fpr0=NULL,wgti=NULL,yes.smooth=F,yes.half=T)
{
    out.yy <- out.pp <- out.AUC <- out.TPR <- out.FPR <- out.PPV <- out.NPV <- NULL
    if(is.null(wgti)){wgti=rep(1,length(Di))}; yyi = as.matrix(yyi); pp=ncol(as.matrix(yyi));  
    mu0 = sum(wgti*(1-Di))/sum(wgti); mu1 = 1-mu0  
    for(k in 1:pp)
    {
        yy = yy0; 
        if(!is.null(fpr0)){
            tpr.all = S.FUN(yyi[,k],Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth,yes.half=yes.half); 
            fpr.all = S.FUN(yyi[,k],Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth,yes.half=yes.half);
            TPR = approx(c(0,fpr.all,1),c(0,tpr.all,1),fpr0,method="linear",rule=2)$y; 
            if(!is.null(yy0)){TPR = c(S.FUN(yy0,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth,yes.half=yes.half), TPR)} 
            yy = c(yy,Sinv.FUN(fpr0,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth))           
            FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth,yes.half=yes.half)
        }else{
            TPR = S.FUN(yy,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth,yes.half=yes.half); 
            FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth,yes.half=yes.half)
        }
        out.yy = cbind(out.yy, yy)
        out.pp = cbind(out.pp, S.FUN(yy,Yi=yyi[,k],wgti,yes.smooth=yes.smooth,yes.half=yes.half))
        out.TPR = cbind(out.TPR,  TPR);  out.FPR  <- cbind(out.FPR,  FPR)
        PPV <- 1/(1+FPR*mu0/(TPR*mu1)); NPV <- 1/(1+(1-TPR)*mu1/((1-FPR)*mu0))
        out.PPV <- cbind(out.PPV, PPV); out.NPV <- cbind(out.NPV, NPV)
        #AUC <- sum((sum.I(yyi[,k],"<=",Yi=yyi[,k],Vi=Di*wgti)+sum.I(yyi[,k],"<",Yi=yyi[,k],Vi=Di*wgti))*(1-Di)*wgti/2
        #             )/(sum((1-Di)*wgti)*sum(Di*wgti))
        AUC = sum(S.FUN(yyi[,k],Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth)*(1-Di)*wgti)/sum((1-Di)*wgti)
        out.AUC <- c(out.AUC, AUC)
    }
    out = c(out.AUC,out.yy,out.pp,out.FPR,out.TPR,out.PPV,out.NPV)
    out
}


