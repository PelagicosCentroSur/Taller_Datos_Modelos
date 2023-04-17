ProbEstatusproy_sept<-function(admb,nesc,name,pF,nyears){
  
  for(i in 1:length(nesc)){
    std<- read.table(paste(admb,nesc[i],".std",sep=""),header=T,sep="",na="NA",fill=T) 
    rpr_y1     <-subset(std,name=="RPRequ3")$value[nyears-1];
    rpr.std_y1 <-subset(std,name=="RPRequ3")$std[nyears-1]
    rpr_y2     <-subset(std,name=="RPRequ3")$value[nyears];
    rpr.std_y2 <-subset(std,name=="RPRequ3")$std[nyears]
    rpr_y3     <-subset(std,name=="RPR_p0")$value[1];
    rpr.std_y3 <-subset(std,name=="RPR_p0")$std[1]
    #***************************************************************
    # biomasa desovante vs BDrms
    #***************************************************************
    xa1  <- rnorm(1000, mean = rpr_y1, sd = rpr.std_y1)
    xa  <-seq(min(xa1),max(xa1),0.005)
    ya  <-dnorm(xa, mean = rpr_y1, sd =rpr.std_y1)
    ica <-qnorm(c(0.05,0.95,0.5),rpr_y1,rpr.std_y1)
    xxa <-c(xa[xa>=ica[1]&xa<=ica[2]],rev(xa[xa>=ica[1]&xa<=ica[2]]))
    yya <-c(ya[xa>=ica[1]&xa<=ica[2]],rep(0,length(ya[xa>=ica[1]&xa<=ica[2]])))
    
    xb1  <- rnorm(1000, mean = rpr_y2, sd = rpr.std_y2)
    xb  <-seq(min(xb1),max(xb1),0.005)
    yb  <-dnorm(xb, mean = rpr_y2, sd =rpr.std_y2)
    icb <-qnorm(c(0.05,0.95,0.5),rpr_y2,rpr.std_y2)
    xxb <-c(xb[xb>=icb[1]&xb<=icb[2]],rev(xb[xb>=icb[1]&xb<=icb[2]]))
    yyb <-c(yb[xb>=icb[1]&xb<=icb[2]],rep(0,length(yb[xb>=icb[1]&xb<=icb[2]])))
    
    xc1  <- rnorm(1000, mean = rpr_y3, sd = rpr.std_y3)
    xc  <-seq(min(xc1),max(xc1),0.005)
    yc  <-dnorm(xc, mean = rpr_y3, sd =rpr.std_y3)
    icc <-qnorm(c(0.05,0.95,0.5),rpr_y3,rpr.std_y3)
    xxc <-c(xc[xc>=icc[1]&xc<=icc[2]],rev(xc[xc>=icc[1]&xc<=icc[2]]))
    yyc <-c(yc[xc>=icc[1]&xc<=icc[2]],rep(0,length(yc[xc>=icc[1]&xc<=icc[2]])))
    #***************************************************************
    #par(mar=c(3,3,2,1)+0.5)
    ifelse((nchar(name) > 1),upspace <- 1.5,upspace <- 0)
    par(mai=c(0.4,0.4,0.05,0.05),oma=c(0.0,0,upspace,0.0))
    par(cex=0.7, mgp=c(1.35,0.35,0))
    plot(xa,ya,type="n",ylab="Densidad de probabilidad",xaxs="i",cex.lab=0.8,
    xlab="BD/BDrms",las=0,yaxs= "i",ylim=c(0,5),xlim=c(0,2.5),cex.axis=0.8)
    if (upspace > 0) mtext(paste(name,pF[i],sep=" "),side=3,outer=T,cex=.7,font=1)
    
    polygon(xxa,yya,col=gray(0.7,0.5),border="gray85")
    polygon(xxb,yyb,col=gray(0.7,0.5),border="gray80")
    polygon(xxc,yyc,col=gray(0.7,0.5),border="gray95")
    lines(xa,ya,lwd=1,col="red",lty=1)
    lines(xb,yb,lwd=1,col="orange",lty=1)
    lines(xc,yc,lwd=1,col="green",lty=1)
    legend(0.1,5,c(paste("BD2018/19_IC95% = [",round(ica[1],2),"-",round(ica[2],2),"]",sep=" "),
                     paste("BD2019/20_IC95% = [",round(icb[1],2),"-",round(icb[2],2),"]",sep=" "),
                     paste("BD2020/21_IC95% = [",round(icc[1],2),"-",round(icc[2],2),"]",sep=" ")),
           lty=c(1,1,1),col=c("red","orange","green"),bty="n",lwd=1,cex=0.7)
    box()
  }
}

ProbEstatusproy_marzo<-function(admb,nesc,name,pF,nyears){
  
  for(i in 1:length(nesc)){
    std<- read.table(paste(admb,nesc[i],".std",sep=""),header=T,sep="",na="NA",fill=T) 
   # rpr_y1     <-subset(std,name=="RPRequ3")$value[nyears-1];
    #rpr.std_y1 <-subset(std,name=="RPRequ3")$std[nyears-1]
    rpr_y2     <-subset(std,name=="RPRequ3")$value[nyears];
    rpr.std_y2 <-subset(std,name=="RPRequ3")$std[nyears]
    rpr_y3     <-subset(std,name=="RPR_p0")$value[1];
    rpr.std_y3 <-subset(std,name=="RPR_p0")$std[1]
    #***************************************************************
    # biomasa desovante vs BDrms
    #***************************************************************
    #xa1  <- rnorm(1000, mean = rpr_y1, sd = rpr.std_y1)
    #xa  <-seq(min(xa1),max(xa1),0.005)
    #ya  <-dnorm(xa, mean = rpr_y1, sd =rpr.std_y1)
    #ica <-qnorm(c(0.05,0.95,0.5),rpr_y1,rpr.std_y1)
    #xxa <-c(xa[xa>=ica[1]&xa<=ica[2]],rev(xa[xa>=ica[1]&xa<=ica[2]]))
    #yya <-c(ya[xa>=ica[1]&xa<=ica[2]],rep(0,length(ya[xa>=ica[1]&xa<=ica[2]])))
    
    xb1  <- rnorm(1000, mean = rpr_y2, sd = rpr.std_y2)
    xb  <-seq(min(xb1),max(xb1),0.005)
    yb  <-dnorm(xb, mean = rpr_y2, sd =rpr.std_y2)
    icb <-qnorm(c(0.05,0.95,0.5),rpr_y2,rpr.std_y2)
    xxb <-c(xb[xb>=icb[1]&xb<=icb[2]],rev(xb[xb>=icb[1]&xb<=icb[2]]))
    yyb <-c(yb[xb>=icb[1]&xb<=icb[2]],rep(0,length(yb[xb>=icb[1]&xb<=icb[2]])))
    
    xc1  <- rnorm(1000, mean = rpr_y3, sd = rpr.std_y3)
    xc  <-seq(min(xc1),max(xc1),0.005)
    yc  <-dnorm(xc, mean = rpr_y3, sd =rpr.std_y3)
    icc <-qnorm(c(0.05,0.95,0.5),rpr_y3,rpr.std_y3)
    xxc <-c(xc[xc>=icc[1]&xc<=icc[2]],rev(xc[xc>=icc[1]&xc<=icc[2]]))
    yyc <-c(yc[xc>=icc[1]&xc<=icc[2]],rep(0,length(yc[xc>=icc[1]&xc<=icc[2]])))
    #***************************************************************
    #par(mar=c(3,3,2,1)+0.5)
    ifelse((nchar(name) > 1),upspace <- 1.5,upspace <- 0)
    par(mai=c(0.45,0.45,0.05,0.05),oma=c(0.0,0,upspace,0.0))
    par(cex=0.7, mgp=c(1.35,0.35,0))
    plot(xb,yb,type="n",ylab="Densidad de probabilidad",xaxs="i",cex.lab=1,
         xlab="BD/BDrms",las=0,yaxs= "i",ylim=c(0,4.5),xlim=c(0,2.5),cex.axis=1)
    if (upspace > 0) mtext(paste(name,pF[i],sep=" "),side=3,outer=T,cex=.7,font=2)
    
    #polygon(xxa,yya,col=gray(0.2,0.5),border="gray85")
    polygon(xxb,yyb,col=gray(0.8,0.7),border="gray80")
    polygon(xxc,yyc,col=gray(0.5,0.7),border="gray95")
   # lines(xa,ya,lwd=1,col=4,lty=1)
    lines(xb,yb,lwd=1,col=4,lty=1)
    lines(xc,yc,lwd=1,col=5,lty=1)
    legend(0.1,4.5,c(paste("BD2019/20_IC95% = [",round(icb[1],2),"-",round(icb[2],2),"]",sep=" "),
                     paste("BD2020/21_IC95% = [",round(icc[1],2),"-",round(icc[2],2),"]",sep=" ")),
           lty=c(1,1),col=c(4,5),bty="n",lwd=1,cex=0.7)
    box()
  }
}
