Estatusproy_sept<-function(admb,mfyr,escR,year1,year2,year3){
  
  rep<- reptoRlist(paste(admb,".rep",sep=""))
  nyears<-length(rep$years)
  
  nmfyr<-length(mfyr)
  
  rpr_y1     <-matrix(ncol=nmfyr,nrow=nyears);
  rpr.std_y1 <-matrix(ncol=nmfyr,nrow=nyears)
  rpr_y2     <-matrix(ncol=nmfyr,nrow=nyears);
  rpr.std_y2 <-matrix(ncol=nmfyr,nrow=nyears)
  rpr_y3     <-matrix(ncol=nmfyr,nrow=nyears);
  rpr.std_y3 <-matrix(ncol=nmfyr,nrow=nyears)
  
  for(i in 1:nmfyr){
    std<-read.table(paste(admb,mfyr[i],".std",sep=""),header=T,sep="",na="NA",fill=T) 
    rpr_y1[,i]     <-subset(std,name=="RPRequ3")$value[nyears-1];
    rpr.std_y1[,i] <-subset(std,name=="RPRequ3")$std[nyears-1]
    rpr_y2[,i]     <-subset(std,name=="RPRequ3")$value[nyears];
    rpr.std_y2[,i] <-subset(std,name=="RPRequ3")$std[nyears]
    rpr_y3[,i]     <-subset(std,name=="RPR_p0")$value[1];
    rpr.std_y3[,i] <-subset(std,name=="RPR_p0")$std[1]
  }
  
  #c?lculo de probabilidades
  p0.9_y1<-rep(0,nmfyr)
  p0.5_y1<-rep(0,nmfyr)
  p0.9_y2<-rep(0,nmfyr)
  p0.5_y2<-rep(0,nmfyr)
  p0.9_y3<-rep(0,nmfyr)
  p0.5_y3<-rep(0,nmfyr)
  
  for(i in 1:nmfyr){
    p0.9_y1[i]<-pnorm(0.9,rpr_y1[,i],rpr.std_y1[,i],lower.tail = TRUE,log.p = F)
    p0.5_y1[i]<-pnorm(0.5,rpr_y1[,i],rpr.std_y1[,i],lower.tail = TRUE,log.p = F)}
  
  for(i in 1:nmfyr){
    p0.9_y2[i]<-pnorm(0.9,rpr_y2[,i],rpr.std_y2[,i],lower.tail = TRUE,log.p = F)
    p0.5_y2[i]<-pnorm(0.5,rpr_y2[,i],rpr.std_y2[,i],lower.tail = TRUE,log.p = F)}
  
  for(i in 1:nmfyr){
    p0.9_y3[i]<-pnorm(0.9,rpr_y3[,i],rpr.std_y3[,i],lower.tail = TRUE,log.p = F)
    p0.5_y3[i]<-pnorm(0.5,rpr_y3[,i],rpr.std_y3[,i],lower.tail = TRUE,log.p = F)}
  
  probEstatus<-round(rbind(p0.9_y1,p0.5_y1,p0.9_y2,p0.5_y2,p0.9_y3,p0.5_y3),2)
  
  rownames(probEstatus)<-c(paste("p(sobre-explotación)_",year1,sep=""),paste("p(colapso)_",year1,sep=""),
                           paste("p(sobre-explotación)_",year2,sep=""),paste("p(colapso)_",year2,sep=""),
                           paste("p(sobre-explotación)_",year3,sep=""),paste("p(colapso)_",year3,sep=""))
  
  colnames(probEstatus)<-c(paste(escR,"[F~RMS~*1]",sep=""),paste("[F~RMS~*0.9]",sep=""),paste("[F~RMS~*0.7]",sep=""))
  
  kable((probEstatus),align='c')
  
  #rprom<-t(probEstatus[,c(1:3,7:9,13:15)]);rprom
  
  #rquiebres<-t(probEstatus[,c(4:6,10:12,16:18)]);rquiebr_es
  
}

Estatusproy_marzo<-function(admb,mfyr,escR,year2,year3){
  
  rep<- reptoRlist(paste(admb,".rep",sep=""))
  nyears<-length(rep$years)
  
  nmfyr<-length(mfyr)
  
  rpr_y1     <-matrix(ncol=nmfyr,nrow=nyears);
  rpr.std_y1 <-matrix(ncol=nmfyr,nrow=nyears)
  rpr_y2     <-matrix(ncol=nmfyr,nrow=nyears);
  rpr.std_y2 <-matrix(ncol=nmfyr,nrow=nyears)
  rpr_y3     <-matrix(ncol=nmfyr,nrow=nyears);
  rpr.std_y3 <-matrix(ncol=nmfyr,nrow=nyears)
  
  for(i in 1:nmfyr){
    std<-read.table(paste(admb,mfyr[i],".std",sep=""),header=T,sep="",na="NA",fill=T) 
    rpr_y1[,i]     <-subset(std,name=="RPRequ3")$value[nyears-1];
    rpr.std_y1[,i] <-subset(std,name=="RPRequ3")$std[nyears-1]
    rpr_y2[,i]     <-subset(std,name=="RPRequ3")$value[nyears];
    rpr.std_y2[,i] <-subset(std,name=="RPRequ3")$std[nyears]
    rpr_y3[,i]     <-subset(std,name=="RPR_p0")$value[1];
    rpr.std_y3[,i] <-subset(std,name=="RPR_p0")$std[1]
  }
  
  #c?lculo de probabilidades
  p0.9_y1<-rep(0,nmfyr)
  p0.5_y1<-rep(0,nmfyr)
  p0.9_y2<-rep(0,nmfyr)
  p0.5_y2<-rep(0,nmfyr)
  p0.9_y3<-rep(0,nmfyr)
  p0.5_y3<-rep(0,nmfyr)
  
  for(i in 1:nmfyr){
    p0.9_y1[i]<-pnorm(0.9,rpr_y1[,i],rpr.std_y1[,i],lower.tail = TRUE,log.p = F)-pnorm(0.5,rpr_y1[,i],rpr.std_y1[,i],lower.tail = TRUE,log.p = F)
    p0.5_y1[i]<-pnorm(0.5,rpr_y1[,i],rpr.std_y1[,i],lower.tail = TRUE,log.p = F)}
  
  for(i in 1:nmfyr){
    p0.9_y2[i]<-pnorm(0.9,rpr_y2[,i],rpr.std_y2[,i],lower.tail = TRUE,log.p = F)-pnorm(0.5,rpr_y2[,i],rpr.std_y2[,i],lower.tail = TRUE,log.p = F)
    p0.5_y2[i]<-pnorm(0.5,rpr_y2[,i],rpr.std_y2[,i],lower.tail = TRUE,log.p = F)}
  
  for(i in 1:nmfyr){
    p0.9_y3[i]<-pnorm(0.9,rpr_y3[,i],rpr.std_y3[,i],lower.tail = TRUE,log.p = F)-pnorm(0.5,rpr_y3[,i],rpr.std_y3[,i],lower.tail = TRUE,log.p = F)
    p0.5_y3[i]<-pnorm(0.5,rpr_y3[,i],rpr.std_y3[,i],lower.tail = TRUE,log.p = F)}
  
  probEstatus<-round(rbind(p0.9_y2,p0.5_y2,p0.9_y3,p0.5_y3),2)
  
  rownames(probEstatus)<-c(paste("p(BD<0,9BD~RMS~)_",year2,sep=""),paste("p(BD<0,5BD~RMS~)_",year2,sep=""),
                           paste("p(BD<0,9BD~RMS~)_",year3,sep=""),paste("p(BD<0,5BD~RMS~)_",year3,sep=""))
  
  colnames(probEstatus)<-c(paste(escR,"[F~RMS~*1]",sep=""),paste("[F~RMS~*0.9]",sep=""),paste("[F~RMS~*0.7]",sep=""))
  
  kable((probEstatus),align='c')
  
  #rprom<-t(probEstatus[,c(1:3,7:9,13:15)]);rprom
  
  #rquiebres<-t(probEstatus[,c(4:6,10:12,16:18)]);rquiebr_es
  
}

