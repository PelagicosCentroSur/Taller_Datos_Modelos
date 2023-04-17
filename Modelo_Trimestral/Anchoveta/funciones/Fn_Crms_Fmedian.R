Crms_Fmedian<-function(dir,admb,system){
  
  dir<-dir.1
  admb<-"MAE323"
  system<-"mac"
  
  setwd(dir)
  #Asesoría de septiembre
  #para mac
  if(system=="mac"){
    system(paste("~/admb-12.2/admb",admb,sep=" "))
    system(paste("./",admb,sep=""))
  }
  if(system=="windows"){
    system(paste("/ADMB/admb",admb,sep=" "))
    system(admb)
  }
  #------------------------------------------------
  dat.file   = paste(admb,".dat",sep="")
  data.0        <- lisread(paste(dir.1,dat.file, sep='/'));
  names(data.0) <-  str_trim(names(data.0), side="right")
  data.1        <- data.0
  rep      <- reptoRlist(paste(admb,".rep",sep=""))                                               
  std      <- read.table(paste(admb,".std",sep=""),header=T,sep="",na="NA",fill=T) 
  
  years  <- data.1$Ind[,1]                                                              
  nyears <- data.1$nanos      
  
  
  Crms<-subset(std,name=="YTP_r0")$value
  FRMS         <- subset(std,name=="log_Fref")$value[1]  
  Fstatuquo   <- exp(subset(std,name=="log_Ft")$value[nyears]) 
  Fmedian<-exp(median(subset(std,name=="log_Ft")$value,na.rm = T))   
  
  data_update<-round(data.frame(valores=rbind(Crms,FRMS,Fstatuquo,Fmedian)),2)
  
  print(data_update)
  
}
  #------------------------------