
Verosimilitud<-function(admb,dir.0,dir.1,Carpeta,logRo,priorRo,l_log_Ro,l_opt_Ro,l_Fase_desRt,system){
  
  dir<-paste(dir.0,Carpeta,sep="")
  
  dat_admb<-paste(admb,".dat",sep="")
  tpl_admb<-paste(admb,".tpl",sep="")
  rep_admb<-paste(admb,".rep",sep="")
  std_admb<-paste(admb,".std",sep="")
  
  
  unlink(dir,recursive=T) #borra la carpeta "CBA2016"
  dir.create(file.path(dir.0,Carpeta))#crea la carpeta "CBA2016"" nuevamente
  setwd(dir.1);file.copy(c(dat_admb,tpl_admb),dir) #copia los archivos de la carpeta MAE0316
  setwd(dir)
  
  # if(system=="mac"){
  #    system(paste("~/admb-12.2/admb",admb,sep=" "))
  #   system(paste("./",admb,sep=""))}
  
  #  if(system=="windows"){
  #   system(paste("admb",admb,sep=" "))
  #  system(admb)}
  
  data.0 <- readLines(dat_admb,encoding="UTF-8")
  data.1 <- data.0
  
  log_Ro    <- logRo
  prior_Ro  <- priorRo
  casos     <- length(prior_Ro);casos
  
  #===========================================================================================================
  # TERCER PASO: CREA LOS ARCHIVOS .DAT DE CADA CASO
  #===========================================================================================================
  for(i in 1:casos){
    data.1[l_log_Ro]  <- prior_Ro[i]
    data.1[l_opt_Ro]  <- -1
    data.1[l_Fase_desRt]  <- 1
    
    cat(data.1,file=(can<-file(paste(admb,i,".dat", sep=''),"wb",encoding="UTF-8")),sep="\n")
    close(can)
    
    setwd(dir.1)
    file.copy(c(paste(admb,".tpl",sep="")),dir)
    
    setwd(dir)
    file.rename(paste(admb,".tpl",sep=""),paste(admb,i,".tpl",sep="")) 
    
    if(system=="mac"){
      system(paste("~/admb-12.2/admb ",admb,i,sep=""))
      system(paste("./",admb,i,sep="")) }
    
    if(system=="windows"){
      system(paste("/ADMB/admb ",admb,i,sep=""))
      system(paste(admb,i,sep="")) }
    
    
    
    file.remove(paste(admb,i,".htp", sep=""),
                paste(admb,i,".cpp", sep=""),
                paste(admb,i,".obj", sep=""),
                paste(admb,i,".p01", sep=""),
                paste(admb,i,".b01", sep=""),
                paste(admb,i,".r01", sep=""),
                paste(admb,i,".p02", sep=""),
                paste(admb,i,".b02", sep=""),
                paste(admb,i,".r02", sep=""),
                paste(admb,i,".p03", sep=""),
                paste(admb,i,".b03", sep=""),
                paste(admb,i,".r03", sep=""),
                paste(admb,i,".p04", sep=""),
                paste(admb,i,".b04", sep=""),
                paste(admb,i,".r04", sep=""),
                paste(admb,i,".p05", sep=""),
                paste(admb,i,".b05", sep=""),
                paste(admb,i,".r05", sep=""),
                paste(admb,i,".p06", sep=""),
                paste(admb,i,".b06", sep=""),
                paste(admb,i,".r06", sep=""),
                paste(admb,i,".p07", sep=""),
                paste(admb,i,".b07", sep=""),
                paste(admb,i,".r07", sep=""),
                paste(admb,i,".p08", sep=""),
                paste(admb,i,".b08", sep=""),
                paste(admb,i,".r08", sep=""),
                paste(admb,i,".p09", sep=""),
                paste(admb,i,".b09", sep=""),
                paste(admb,i,".r09", sep=""),
                paste(admb,i,".p10", sep=""),
                paste(admb,i,".b10", sep=""),
                paste(admb,i,".r10", sep=""),
                paste(admb,i,".par", sep=""),
                paste(admb,i,".bar", sep=""),
                paste(admb,i,".eva", sep=""),
                paste(admb,i,".cor", sep=""),
                paste(admb,i,".std", sep=""),
                paste(admb,i,".log", sep=""),
                paste(admb,i,".exe", sep=""))
    
}
}
