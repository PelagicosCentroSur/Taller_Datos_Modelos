
CBA<-function(dir.0,dir.1,Carpeta,admb,l_opt_proy,l_opRec,l_mf,opt_proy,system){
  
  dir.5<-paste(dir.0,Carpeta,sep="")
  
  dat_admb<-paste(admb,".dat",sep="")
  exe_admb<-paste(admb,".exe",sep="")
  tpl_admb<-paste(admb,".tpl",sep="")

  unlink(dir.5,recursive=T) #borra la carpeta "CBA_sept"
  dir.create(file.path(dir.0,Carpeta))#crea la carpeta "CBA_sept"" nuevamente
  
  setwd(dir.1);file.copy(c(dat_admb,tpl_admb),dir.5) #copia los archivos de la carpeta calendario
  
  setwd(dir.5);
  
  if(system=="mac"){
    system(paste("~/admb-12.2/admb",admb,sep=" "))
    system(paste("./",admb,sep=""))
  }

  if(system=="windows"){
    system(paste("/ADMB/admb",admb,sep=" "))
    system(admb)
  }
  
  
  std_admb<-paste(admb,".std",sep="")
  rep_admb<-paste(admb,".rep",sep="")
  
  opRec<-seq(1,3,1) #//ESCENARIO DE RECLUTAMIENTO PROMEDIO
  mf<-c(1,0.9,1.1)
  
  data_file<-(paste(admb,".dat",sep="")) #// archivos de datos originales
  S<-readLines(data_file,encoding="UTF-8")
  Se<-S
  
  for( i in 1:3){
    for(j in 1:3){
      Se[l_opt_proy]<-opt_proy #opt_proy 
      Se[l_opRec]<-opRec[i] #opRec
      Se[l_mf]<-mf[j]
 
      cat(Se,file=(can<-file(paste(admb,".dat",sep=""),"wb",encoding="UTF-8")),sep="\n")
      close(can)
      
      if(system=="mac"){
        system(paste("./",admb,sep=""))
      }
      
      if(system=="windows"){
        system(admb)
      }
      
      
      Std <- readLines(std_admb,encoding="UTF-8") 
      cat(Std,file=(can<-file(paste(admb,j,i,".std",sep=""),"wb",encoding="UTF-8")),sep="\n");close(can)
      Rep <- readLines(rep_admb,encoding="UTF-8") 
      cat(Rep,file=(can<-file(paste(admb,j,i,".rep",sep=""),"wb",encoding="UTF-8")),sep="\n");close(can)
      
    }}
  cat(S,file=(can<-file(paste(admb,".dat",sep=""),"wb",encoding="UTF-8")),sep="\n")
  close(can)
}


CreaDataProybase<-function(dir.0,carpetaCBA,admb,yearProy,dir.Rdata,Hito){
  
  reps1a     <- reptoRlist(here(carpetaCBA,paste(admb,"11.rep",sep="")))  
  reps2a     <- reptoRlist(here(carpetaCBA,paste(admb,"12.rep",sep="")))
  reps3a     <- reptoRlist(here(carpetaCBA,paste(admb,"13.rep",sep=""))) 
  
  
  stds1     <- read.table(here(carpetaCBA,paste(admb,"11.std", sep='')),
                          header=T,sep="",na="NA",fill=T) 
  stds2     <- read.table(here(carpetaCBA,paste(admb,"12.std", sep='')),
                          header=T,sep="",na="NA",fill=T) 
  stds3     <- read.table(here(carpetaCBA,paste(admb,"13.std", sep='')),
                          header=T,sep="",na="NA",fill=T) 
  
  # Biomasa desovante
  bds1     <- subset(stds1,name=="BD_p0")$value ; 
  bds1std  <- subset(stds1,name=="BD_p0")$std 
  bds2     <- subset(stds2,name=="BD_p0")$value ; 
  bds2std  <- subset(stds2,name=="BD_p0")$std 
  bds3     <- subset(stds3,name=="BD_p0")$value ; 
  bds3std  <- subset(stds3,name=="BD_p0")$std 
  
  # razon BD/BDrms
  RpRps1     <- subset(stds1,name=="RPR_p0")$value ; 
  RpRps1std  <- subset(stds1,name=="RPR_p0")$std 
  RpRps2     <- subset(stds2,name=="RPR_p0")$value ; 
  RpRps2std  <- subset(stds2,name=="RPR_p0")$std 
  RpRps3     <- subset(stds3,name=="RPR_p0")$value ; 
  RpRps3std  <- subset(stds3,name=="RPR_p0")$std 
  
  # Capturas
  cs1     <- subset(stds1,name=="YTP_p0")$value ; 
  cs1std  <- subset(stds1,name=="YTP_p0")$std 
  cs2     <- subset(stds2,name=="YTP_p0")$value ; 
  cs2std  <- subset(stds2,name=="YTP_p0")$std 
  cs3     <- subset(stds3,name=="YTP_p0")$value ; 
  cs3std  <- subset(stds3,name=="YTP_p0")$std 
  
  #CAPTURA ESTIMADA ULTIMO AÑO EVALUACIÓN
  
  #Captura año biológico actual
  cs0     <- subset(stds1,name=="YTP_r0")$value ; 
  cs0std  <- subset(stds1,name=="YTP_r0")$std 
  
  #Captura primer semestre (Captura año biológico actual - desembarque2dosem)
  cs0D     <- subset(stds1,name=="YTP_r01ersem")$value ; 
  cs0Dstd  <- subset(stds1,name=="YTP_r01ersem")$std 
  
  #Captura primer semestre (año biológico actual - desembarque2dosem - remanente)
  cs0R     <- subset(stds1,name=="YTP_r01ersemR")$value ; 
  cs0Rstd  <- subset(stds1,name=="YTP_r01ersemR")$std 
  
  #######################################################################
  # DESCUENTO DEL DESCARTE
  #######################################################################
  
  #Captura año biológico actual - descarte
  cs0d     <- subset(stds1,name=="YTP_r0d")$value ; 
  cs0dstd  <- subset(stds1,name=="YTP_r0d")$std 
  
  #Captura primer semestre (Captura año biológico actual - desembarque2dosem) - descarte
  cs01erS     <- subset(stds1,name=="YTP_r0d1ersem")$value ; 
  cs01erSstd  <- subset(stds1,name=="YTP_r0d1ersem")$std 
  
  #Captura primer semestre (Captura año biológico actual - desembarque2dosem - remanente) - descarte
  cs01erSR     <- subset(stds1,name=="YTP_r0d1ersemR")$value ; 
  cs01erSRstd  <- subset(stds1,name=="YTP_r0d1ersemR")$std 
  
  # Reclutas proyectatos
  Rtp              <- reps1a$Reclutas
  SSBp             <- reps1a$SSB
  BRMSp            <- reps1a$SSBpbr[3]
  FRMSp            <- reps1a$Fs[2]
  desembarquepredp <- reps1a$desembarquepred
  
  Rs1    <- reps1a$Np[1]
  Rs2    <- reps2a$Np[1]
  Rs3    <- reps3a$Np[1]
  
  years <- reps1a$years
  nyears<- length(years)
  
  yearsp <- c(years ,yearProy)
  nyears1<- length(yearsp)
  
  yearbiol_proy  <-c(paste(years[nyears],"/",str_sub(years[nyears]+1,3,4),sep=""),
                     paste(years[nyears]+1,"/",str_sub(years[nyears]+2,3,4),sep=""))
  
  # Data frame
  dataRproy<-data.frame(S1=c(Rtp,rep(Rs1,2)),
                        S2=c(Rtp,rep(Rs2,2)),
                        S3=c(Rtp,rep(Rs3,2)),
                        base=c(Rtp,rep(NA,2))) %>% 
    mutate(years=yearsp,indicador='Rt') %>% 
    melt(id.vars=c('years','indicador'))
  dataBDproy<-data.frame(S1=c(SSBp,bds1)/1000,
                         S2=c(SSBp,bds2)/1000,
                         S3=c(SSBp,bds3)/1000,
                         base=c(SSBp,rep(NA,2))/1000) %>% 
    mutate(years=yearsp,indicador='BD')%>% 
    melt(id.vars=c('years','indicador'))
  dataBD_rmsproy<-data.frame(S1=round(c(SSBp,bds1)/BRMSp,2),
                             S2=round(c(SSBp,bds2)/BRMSp,2),
                             S3=round(c(SSBp,bds3)/BRMSp,2),
                             base=round(c(SSBp/BRMSp,rep(NA,2)),2)) %>% 
    mutate(years=yearsp,indicador='BD_BDrms')%>% 
    melt(id.vars=c('years','indicador'))
  dataCt_proy<-data.frame(S1=c(desembarquepredp,cs1)/1000,
                          S2=c(desembarquepredp,cs2)/1000,
                          S3=c(desembarquepredp,cs3)/1000,
                          base=c(desembarquepredp,rep(NA,2))/1000) %>% 
    mutate(years=yearsp,indicador='Ct')%>% 
    melt(id.vars=c('years','indicador'))
  
  DataProy<-rbind(dataRproy,dataBDproy,dataBD_rmsproy,dataCt_proy)
  
  
  # RECLUTAMIENTO ESTIMADO ULTIMO AÑO EVALUACIÓN (AÑO ACTUAL) ----
  RTs0s      <- subset(stds1,name=="Reclutas")$value[nyears] ;
  RTs0s_std  <- subset(stds1,name=="Reclutas")$std[nyears]
  
  #BIOMASA DESOVANTE ESTIMADA ULTIMO AÑO EVALUACIÓN ----
  bds0s      <- subset(stds1,name=="SSB")$value[nyears] ;
  bds0s_std  <- subset(stds1,name=="SSB")$std[nyears]
  
  # aporte del grupo de edad 0 (reclutamiento) año actual ----
  C1eryearR1act  <-round(reps1a$YTP_r0W_actual[1]/sum(reps1a$YTP_r0W_actual),2)
  C1eryearR1act2 <-round(reps1a$YTP_r0W_actual[2]/sum(reps1a$YTP_r0W_actual),2)
  # aporte del grupo de edad 0 (reclutamiento) 1er año proyectado ----
  C1eryearR1     <-round(reps1a$YTP_p0W_proyectada[1,1]/sum(reps1a$YTP_p0W_proyectada[1,]),2)
  C1eryearR2     <-round(reps2a$YTP_p0W_proyectada[1,1]/sum(reps2a$YTP_p0W_proyectada[1,]),2)
  C1eryearR3     <-round(reps3a$YTP_p0W_proyectada[1,1]/sum(reps3a$YTP_p0W_proyectada[1,]),2)
  # aporte del grupo de edad 1 (reclutamiento) 1er año proyectado ----
  C1eryearR1a    <-round(reps1a$YTP_p0W_proyectada[1,2]/sum(reps1a$YTP_p0W_proyectada[1,]),2)
  C1eryearR2a    <-round(reps2a$YTP_p0W_proyectada[1,2]/sum(reps2a$YTP_p0W_proyectada[1,]),2)
  C1eryearR3a    <-round(reps3a$YTP_p0W_proyectada[1,2]/sum(reps3a$YTP_p0W_proyectada[1,]),2)
  # aporte del grupo de edad 0 (reclutamiento) 2do año proyectado ----
  C1eryearR12    <-round(reps1a$YTP_p0W_proyectada[2,1]/sum(reps1a$YTP_p0W_proyectada[2,]),2)
  C1eryearR22    <-round(reps2a$YTP_p0W_proyectada[2,1]/sum(reps2a$YTP_p0W_proyectada[2,]),2)
  C1eryearR32    <-round(reps3a$YTP_p0W_proyectada[2,1]/sum(reps3a$YTP_p0W_proyectada[2,]),2)
  # aporte del grupo de edad 1  2do año proyectado ----
  C1eryearR12a   <-round(reps1a$YTP_p0W_proyectada[2,2]/sum(reps1a$YTP_p0W_proyectada[2,]),2)
  C1eryearR22a   <-round(reps2a$YTP_p0W_proyectada[2,2]/sum(reps2a$YTP_p0W_proyectada[2,]),2)
  C1eryearR32a   <-round(reps3a$YTP_p0W_proyectada[2,2]/sum(reps3a$YTP_p0W_proyectada[2,]),2)
  
  # ESTATUS PROYECTADO ----
  
  # PRIMER AÑO PROYECTADO ----
  ### *Probabilidad de estar bajo BRMS* ----
  pa1 <-pnorm(0.9,RpRps1[1],RpRps1std[1],lower.tail = TRUE,log.p = F)
  pa2 <-pnorm(0.9,RpRps2[1],RpRps2std[1],lower.tail = TRUE,log.p = F)
  pa3 <-pnorm(0.9,RpRps3[1],RpRps3std[1],lower.tail = TRUE,log.p = F)
  ### *Probabilidad de estar en zona de sobreexplotacion* ----
  pc1 <-pnorm(0.9,RpRps1[1],RpRps1std[1],lower.tail = TRUE,log.p = F)-pnorm(0.5,
        RpRps1[1],RpRps1std[1],lower.tail = TRUE,log.p = F)
  pc2 <-pnorm(0.9,RpRps2[1],RpRps2std[1],lower.tail = TRUE,log.p = F)-pnorm(0.5,
        RpRps2[1],RpRps2std[1],lower.tail = TRUE,log.p = F)
  pc3 <-pnorm(0.9,RpRps3[1],RpRps3std[1],lower.tail = TRUE,log.p = F)-pnorm(0.5,
        RpRps3[1],RpRps3std[1],lower.tail = TRUE,log.p = F)
  ### *Probabilidad de estar en zona de colapso*----
  pd1 <-pnorm(0.5,RpRps1[1],RpRps1std[1],lower.tail = TRUE,log.p = F)
  pd2 <-pnorm(0.5,RpRps2[1],RpRps2std[1],lower.tail = TRUE,log.p = F)
  pd3 <-pnorm(0.5,RpRps3[1],RpRps3std[1],lower.tail = TRUE,log.p = F)
  
  # SEGUNDO AÑO PROYECTADO ----
  
  ### *Probabilidad de estar bajo BRMS* ----
  pa12 <-pnorm(0.9,RpRps1[2],RpRps1std[2],lower.tail = TRUE,log.p = F)
  pa22 <-pnorm(0.9,RpRps2[2],RpRps2std[2],lower.tail = TRUE,log.p = F)
  pa32 <-pnorm(0.9,RpRps3[2],RpRps3std[2],lower.tail = TRUE,log.p = F)
  ### *Probabilidad de estar en zona de sobreexplotacion* ----
  pc12 <-pnorm(0.9,RpRps1[2],RpRps1std[2],lower.tail = TRUE,log.p = F)-pnorm(0.5,
         RpRps1[2],RpRps1std[2],lower.tail = TRUE,log.p = F)
  pc22 <-pnorm(0.9,RpRps2[2],RpRps2std[2],lower.tail = TRUE,log.p = F)-pnorm(0.5,
         RpRps2[2],RpRps2std[2],lower.tail = TRUE,log.p = F)
  pc32 <-pnorm(0.9,RpRps3[2],RpRps3std[2],lower.tail = TRUE,log.p = F)-pnorm(0.5,
         RpRps3[2],RpRps3std[2],lower.tail = TRUE,log.p = F)
  ### *Probabilidad de estar en zona de colapso* ----
  pd12 <-pnorm(0.5,RpRps1[2],RpRps1std[2],lower.tail = TRUE,log.p = F)
  pd22 <-pnorm(0.5,RpRps2[2],RpRps2std[2],lower.tail = TRUE,log.p = F)
  pd32 <-pnorm(0.5,RpRps3[2],RpRps3std[2],lower.tail = TRUE,log.p = F)
  
  
  # CBA INICIAL ----
  n<-3
  q       <- seq(0.1,0.5,0.1)  # niveles de riesgo (cuantiles)
  nq      <- length(q)
  CBA_sept     <- matrix(ncol=nq,nrow=n)
  CBAp_sept    <- rep(0,n)
  CBApstd_sept <- rep(0,n)
  
  buffer   <- matrix(ncol=nq,nrow=n)
  descarte <- matrix(ncol=nq,nrow=n)
  
  for(i in 1:n){
    std             <-read.table(here(carpetaCBA,paste(admb,"1",i,".std",sep="")),header=T,sep="",na="NA",fill=T)
    CBAp_sept[i]    <-subset(std,name=="CBA_c0")$value[1]
    CBApstd_sept[i] <-subset(std,name=="CBA_c0")$std[1]
    for(j in 1:nq){
      CBA_sept[i,j]   <-qnorm(q[j],CBAp_sept[i],CBApstd_sept[i])
    }
  }
  
  for(i in 1:n){for(j in 1:nq){
    buffer[i,j]<-round(1-CBA_sept[i,j]/CBA_sept[i,5],2)}}
  
  # CBA INICIAL MENOS DESCARTE ----
  n<-3
  q       <- seq(0.1,0.5,0.1)  # niveles de riesgo (cuantiles)
  nq      <- length(q)
  CBAd_sept     <- matrix(ncol=nq,nrow=n)
  CBApd_sept    <- rep(0,n)
  CBApdstd_sept <- rep(0,n)
  
  
  for(i in 1:n){
    std              <-read.table(here(carpetaCBA,paste(admb,"1",i,".std",sep="")),header=T,sep="",na="NA",fill=T)
    CBApd_sept[i]    <-subset(std,name=="CBA_c0d")$value[1]
    CBApdstd_sept[i] <-subset(std,name=="CBA_c0d")$std[1]
    for(j in 1:nq){CBAd_sept[i,j]<-qnorm(q[j],CBApd_sept[i],CBApdstd_sept[i])}}
  
  
  #===============================================================================
  # escenario 1
  #===============================================================================
  xca1 <-rnorm(1000, mean = CBAp_sept[1],sd=CBApstd_sept[1])
  xca  <-seq(min(xca1),max(xca1),0.5)
  yca  <-dnorm(xca, mean = CBAp_sept[1],sd=CBApstd_sept[1])
  icca <-qnorm(c(0.1,0.5,0.5),CBAp_sept[1],CBApstd_sept[1])
  xxca <-c(xca[xca>=icca[1]&xca<=icca[2]],
           rev(xca[xca>=icca[1]&xca<=icca[2]]))
  yyca <-c(yca[xca>=icca[1]&xca<=icca[2]],
           rep(0,length(yca[xca>=icca[1]&xca<=icca[2]])))
  #===============================================================================
  # escenario 2
  #===============================================================================
  xcb1 <-rnorm(1000, mean = CBAp_sept[2],sd= CBApstd_sept[2])
  xcb  <-seq(min(xcb1),max(xcb1),0.5)
  ycb  <-dnorm(xcb, mean = CBAp_sept[2],sd= CBApstd_sept[2])
  iccb <-qnorm(c(0.1,0.5,0.5),CBAp_sept[2],CBApstd_sept[2])
  xxcb <-c(xcb[xcb>=iccb[1]&xcb<=iccb[2]],
           rev(xcb[xcb>=iccb[1]&xcb<=iccb[2]]))
  yycb <-c(ycb[xcb>=iccb[1]&xcb<=iccb[2]],
           rep(0,length(ycb[xcb>=iccb[1]&xcb<=iccb[2]])))
  #===============================================================================
  # escenario 3
  #===============================================================================
  xcc1 <-rnorm(1000, mean = CBAp_sept[3],sd=CBApstd_sept[3])
  xcc  <-seq(min(xcc1),max(xcc1),0.5)
  ycc  <-dnorm(xcc, mean = CBAp_sept[3],sd=CBApstd_sept[3])
  iccc <-qnorm(c(0.1,0.5,0.5),CBAp_sept[3],CBApstd_sept[3])
  xxcc <-c(xcc[xcc>=iccc[1]&xcc<=iccc[2]],
           rev(xcc[xcc>=iccc[1]&xcc<=iccc[2]]))
  yycc <-c(ycc[xcc>=iccc[1]&xcc<=iccc[2]],
           rep(0,length(ycc[xcc>=iccc[1]&xcc<=iccc[2]])))
  
  save(list=ls(all=T),
       file=here("Rdata",paste('Datos_Proy_',Hito,'.RData',sep="")))
  
  
  
  
  
}


CreaDataProybase_HITO2Y3<-function(dir.0,carpetaCBA,admb,yearProy,dir.Rdata,Hito){
  
  
  data1a        <- lisread(here(carpetaCBA,paste(admb,".dat",sep=""))) 
  names(data1a) <- str_trim(names(data1a), side="right")
  
  
  reps1a     <- reptoRlist(here(carpetaCBA,paste(admb,"11.rep",sep="")))  
  reps2a     <- reptoRlist(here(carpetaCBA,paste(admb,"12.rep",sep="")))
  reps3a     <- reptoRlist(here(carpetaCBA,paste(admb,"13.rep",sep=""))) 
  
  
  stds1     <- read.table(here(carpetaCBA,paste(admb,"11.std", sep='')),
                          header=T,sep="",na="NA",fill=T) 
  stds2     <- read.table(here(carpetaCBA,paste(admb,"12.std", sep='')),
                          header=T,sep="",na="NA",fill=T) 
  stds3     <- read.table(here(carpetaCBA,paste(admb,"13.std", sep='')),
                          header=T,sep="",na="NA",fill=T) 
  
  # Biomasa desovante
  bds1     <- subset(stds1,name=="BD_p0")$value[1] ; 
  bds1std  <- subset(stds1,name=="BD_p0")$std[1] 
  bds2     <- subset(stds2,name=="BD_p0")$value[1] ; 
  bds2std  <- subset(stds2,name=="BD_p0")$std[1] 
  bds3     <- subset(stds3,name=="BD_p0")$value[1] ; 
  bds3std  <- subset(stds3,name=="BD_p0")$std[1] 
  
  # razon BD/BDrms
  RpRps1     <- subset(stds1,name=="RPR_p0")$value[1] ; 
  RpRps1std  <- subset(stds1,name=="RPR_p0")$std[1] 
  RpRps2     <- subset(stds2,name=="RPR_p0")$value[1] ; 
  RpRps2std  <- subset(stds2,name=="RPR_p0")$std[1] 
  RpRps3     <- subset(stds3,name=="RPR_p0")$value[1] ; 
  RpRps3std  <- subset(stds3,name=="RPR_p0")$std[1] 
  
  # Capturas
  cs1     <- subset(stds1,name=="YTP_p0")$value[1] ; 
  cs1std  <- subset(stds1,name=="YTP_p0")$std[1] 
  cs2     <- subset(stds2,name=="YTP_p0")$value[1] ; 
  cs2std  <- subset(stds2,name=="YTP_p0")$std[1] 
  cs3     <- subset(stds3,name=="YTP_p0")$value[1] ; 
  cs3std  <- subset(stds3,name=="YTP_p0")$std[1] 
  
  #CAPTURA ESTIMADA ULTIMO AÑO EVALUACIÓN
  desem2doSem <- data1a$Desemb2doSem
  remanente   <- data1a$Remanente
  #Captura año biológico actual
  cs0     <- subset(stds1,name=="YTP_r0")$value ; 
  cs0std  <- subset(stds1,name=="YTP_r0")$std 
  
  #Captura primer semestre (Captura año biológico actual - desembarque2dosem)
  cs0D     <- subset(stds1,name=="YTP_r01ersem")$value ; 
  cs0Dstd  <- subset(stds1,name=="YTP_r01ersem")$std 
  
  #Captura primer semestre (año biológico actual - desembarque2dosem - remanente)
  cs0R     <- subset(stds1,name=="YTP_r01ersemR")$value ; 
  cs0Rstd  <- subset(stds1,name=="YTP_r01ersemR")$std 
  
  #######################################################################
  # DESCUENTO DEL DESCARTE
  #######################################################################
  
  #Captura año biológico actual - descarte
  cs0d     <- subset(stds1,name=="YTP_r0d")$value ; 
  cs0dstd  <- subset(stds1,name=="YTP_r0d")$std 
  
  #Captura primer semestre (Captura año biológico actual - desembarque2dosem) - descarte
  cs01erS     <- subset(stds1,name=="YTP_r0d1ersem")$value ; 
  cs01erSstd  <- subset(stds1,name=="YTP_r0d1ersem")$std 
  
  #Captura primer semestre (Captura año biológico actual - desembarque2dosem - remanente) - descarte
  cs01erSR     <- subset(stds1,name=="YTP_r0d1ersemR")$value ; 
  cs01erSRstd  <- subset(stds1,name=="YTP_r0d1ersemR")$std 
  
  # Reclutas proyectatos
  Rtp              <- reps1a$Reclutas
  SSBp             <- reps1a$SSB
  BRMSp            <- reps1a$SSBpbr[3]
  FRMSp            <- reps1a$Fs[2]
  desembarquepredp <- reps1a$desembarquepred
  
  Rs1    <- reps1a$Np[1]
  Rs2    <- reps2a$Np[1]
  Rs3    <- reps3a$Np[1]
  
  years <- reps1a$years
  nyears<- length(years)
  
  yearsp <- c(years ,yearProy)
  nyears1<- length(yearsp)
  
  yearbiol_proy  <-c(paste(years[nyears],"/",str_sub(years[nyears]+1,3,4),sep=""))
  
  # Data frame
  dataRproy<-data.frame(S1=c(Rtp,rep(Rs1,1)),
                        S2=c(Rtp,rep(Rs2,1)),
                        S3=c(Rtp,rep(Rs3,1)),
                        base=c(Rtp,rep(NA,1))) %>% 
    mutate(years=yearsp,indicador='Rt') %>% 
    melt(id.vars=c('years','indicador'))
  dataBDproy<-data.frame(S1=c(SSBp,bds1)/1000,
                         S2=c(SSBp,bds2)/1000,
                         S3=c(SSBp,bds3)/1000,
                         base=c(SSBp,rep(NA,1))/1000) %>% 
    mutate(years=yearsp,indicador='BD')%>% 
    melt(id.vars=c('years','indicador'))
  dataBD_rmsproy<-data.frame(S1=round(c(SSBp/BRMSp,RpRps1),2),
                             S2=round(c(SSBp/BRMSp,RpRps2),2),
                             S3=round(c(SSBp/BRMSp,RpRps3),2),
                             base=round(c(SSBp/BRMSp,rep(NA,1)),2)) %>% 
    mutate(years=yearsp,indicador='BD_BDrms')%>% 
    melt(id.vars=c('years','indicador'))
  dataCt_proy<-data.frame(S1=c(desembarquepredp,cs1)/1000,
                          S2=c(desembarquepredp,cs2)/1000,
                          S3=c(desembarquepredp,cs3)/1000,
                          base=c(desembarquepredp,rep(NA,1))/1000) %>% 
    mutate(years=yearsp,indicador='Ct')%>% 
    melt(id.vars=c('years','indicador'))
  
  DataProy<-rbind(dataRproy,dataBDproy,dataBD_rmsproy,dataCt_proy)
  
  
  # RECLUTAMIENTO ESTIMADO ULTIMO AÑO EVALUACIÓN (AÑO ACTUAL) ----
  RTs0s      <- subset(stds1,name=="Reclutas")$value[nyears] ;
  RTs0s_std  <- subset(stds1,name=="Reclutas")$std[nyears]
  
  #BIOMASA DESOVANTE ESTIMADA ULTIMO AÑO EVALUACIÓN ----
  bds0s      <- subset(stds1,name=="SSB")$value[nyears] ;
  bds0s_std  <- subset(stds1,name=="SSB")$std[nyears]
  
  # aporte del grupo de edad 0 (reclutamiento) año actual ----
  C1eryearR1act  <-round(reps1a$YTP_r0W_actual[1]/sum(reps1a$YTP_r0W_actual),2)
  C1eryearR1act2 <-round(reps1a$YTP_r0W_actual[2]/sum(reps1a$YTP_r0W_actual),2)
  # aporte del grupo de edad 0 (reclutamiento) 1er año proyectado ----
  C1eryearR1     <-round(reps1a$YTP_p0W_proyectada[1,1]/sum(reps1a$YTP_p0W_proyectada[1,]),2)
  C1eryearR2     <-round(reps2a$YTP_p0W_proyectada[1,1]/sum(reps2a$YTP_p0W_proyectada[1,]),2)
  C1eryearR3     <-round(reps3a$YTP_p0W_proyectada[1,1]/sum(reps3a$YTP_p0W_proyectada[1,]),2)
  # aporte del grupo de edad 1 (reclutamiento) 1er año proyectado ----
  C1eryearR1a    <-round(reps1a$YTP_p0W_proyectada[1,2]/sum(reps1a$YTP_p0W_proyectada[1,]),2)
  C1eryearR2a    <-round(reps2a$YTP_p0W_proyectada[1,2]/sum(reps2a$YTP_p0W_proyectada[1,]),2)
  C1eryearR3a    <-round(reps3a$YTP_p0W_proyectada[1,2]/sum(reps3a$YTP_p0W_proyectada[1,]),2)
 
  # ESTATUS PROYECTADO ----
  
  # PRIMER AÑO PROYECTADO ----
  ### *Probabilidad de estar bajo BRMS* ----
  pa1 <-pnorm(0.9,RpRps1[1],RpRps1std[1],lower.tail = TRUE,log.p = F)
  pa2 <-pnorm(0.9,RpRps2[1],RpRps2std[1],lower.tail = TRUE,log.p = F)
  pa3 <-pnorm(0.9,RpRps3[1],RpRps3std[1],lower.tail = TRUE,log.p = F)
  ### *Probabilidad de estar en zona de sobreexplotacion* ----
  pc1 <-pnorm(0.9,RpRps1[1],RpRps1std[1],lower.tail = TRUE,log.p = F)-pnorm(0.5,
        RpRps1[1],RpRps1std[1],lower.tail = TRUE,log.p = F)
  pc2 <-pnorm(0.9,RpRps2[1],RpRps2std[1],lower.tail = TRUE,log.p = F)-pnorm(0.5,
        RpRps2[1],RpRps2std[1],lower.tail = TRUE,log.p = F)
  pc3 <-pnorm(0.9,RpRps3[1],RpRps3std[1],lower.tail = TRUE,log.p = F)-pnorm(0.5,
        RpRps3[1],RpRps3std[1],lower.tail = TRUE,log.p = F)
  ### *Probabilidad de estar en zona de colapso*----
  pd1 <-pnorm(0.5,RpRps1[1],RpRps1std[1],lower.tail = TRUE,log.p = F)
  pd2 <-pnorm(0.5,RpRps2[1],RpRps2std[1],lower.tail = TRUE,log.p = F)
  pd3 <-pnorm(0.5,RpRps3[1],RpRps3std[1],lower.tail = TRUE,log.p = F)
  
 
  
  # CBA HITO 2  ----
  n<-3
  q       <- seq(0.1,0.5,0.1)  # niveles de riesgo (cuantiles)
  nq      <- length(q)
  CBA_marzo     <- matrix(ncol=nq,nrow=n)
  CBAp_marzo    <- rep(0,n)
  CBApstd_marzo <- rep(0,n)
  
  buffer   <- matrix(ncol=nq,nrow=n)
  descarte <- matrix(ncol=nq,nrow=n)
  
  for(i in 1:n){
    std              <-read.table(here(carpetaCBA,paste(admb,"1",i,".std",sep="")),header=T,sep="",na="NA",fill=T)
    CBAp_marzo[i]    <-subset(std,name=="CBA_c0")$value[1]
    CBApstd_marzo[i] <-subset(std,name=="CBA_c0")$std[1]
    for(j in 1:nq){
      CBA_marzo[i,j] <-qnorm(q[j],CBAp_marzo[i],CBApstd_marzo[i])
    }
  }
  
  for(i in 1:n){for(j in 1:nq){
    buffer[i,j]<-round(1-CBA_marzo[i,j]/CBA_marzo[i,5],2)}}
  
  # CBA MENOS DESCARTE ----
  n<-3
  q       <- seq(0.1,0.5,0.1)  # niveles de riesgo (cuantiles)
  nq      <- length(q)
  CBAd_marzo     <- matrix(ncol=nq,nrow=n)
  CBApd_marzo    <- rep(0,n)
  CBApdstd_marzo <- rep(0,n)
  
  
  for(i in 1:n){
    std               <-read.table(here(carpetaCBA,paste(admb,"1",i,".std",sep="")),header=T,sep="",na="NA",fill=T)
    CBApd_marzo[i]    <-subset(std,name=="CBA_c0d")$value[1]
    CBApdstd_marzo[i] <-subset(std,name=="CBA_c0d")$std[1]
    for(j in 1:nq){
      CBAd_marzo[i,j] <-qnorm(q[j],CBApd_marzo[i],CBApdstd_marzo[i])}}
  
  
  #===============================================================================
  # escenario 1
  #===============================================================================
  xca1 <-rnorm(1000, mean = CBAp_marzo[1],sd=CBApstd_marzo[1])
  xca  <-seq(min(xca1),max(xca1),0.5)
  yca  <-dnorm(xca, mean = CBAp_marzo[1],sd=CBApstd_marzo[1])
  icca <-qnorm(c(0.1,0.5,0.5),CBAp_marzo[1],CBApstd_marzo[1])
  xxca <-c(xca[xca>=icca[1]&xca<=icca[2]],
           rev(xca[xca>=icca[1]&xca<=icca[2]]))
  yyca <-c(yca[xca>=icca[1]&xca<=icca[2]],
           rep(0,length(yca[xca>=icca[1]&xca<=icca[2]])))
  #===============================================================================
  # escenario 2
  #===============================================================================
  xcb1 <-rnorm(1000, mean = CBAp_marzo[2],sd= CBApstd_marzo[2])
  xcb  <-seq(min(xcb1),max(xcb1),0.5)
  ycb  <-dnorm(xcb, mean = CBAp_marzo[2],sd= CBApstd_marzo[2])
  iccb <-qnorm(c(0.1,0.5,0.5),CBAp_marzo[2],CBApstd_marzo[2])
  xxcb <-c(xcb[xcb>=iccb[1]&xcb<=iccb[2]],
           rev(xcb[xcb>=iccb[1]&xcb<=iccb[2]]))
  yycb <-c(ycb[xcb>=iccb[1]&xcb<=iccb[2]],
           rep(0,length(ycb[xcb>=iccb[1]&xcb<=iccb[2]])))
  #===============================================================================
  # escenario 3
  #===============================================================================
  xcc1 <-rnorm(1000, mean = CBAp_marzo[3],sd=CBApstd_marzo[3])
  xcc  <-seq(min(xcc1),max(xcc1),0.5)
  ycc  <-dnorm(xcc, mean = CBAp_marzo[3],sd=CBApstd_marzo[3])
  iccc <-qnorm(c(0.1,0.5,0.5),CBAp_marzo[3],CBApstd_marzo[3])
  xxcc <-c(xcc[xcc>=iccc[1]&xcc<=iccc[2]],
           rev(xcc[xcc>=iccc[1]&xcc<=iccc[2]]))
  yycc <-c(ycc[xcc>=iccc[1]&xcc<=iccc[2]],
           rep(0,length(ycc[xcc>=iccc[1]&xcc<=iccc[2]])))

  ## ----------------------------------------------------------
  # CBA HITO 2  MENOS DESEMBARQUE Y SALDO CUOTA AÑO PREVIO----
  n<-3
  q       <- seq(0.1,0.5,0.1)  # niveles de riesgo (cuantiles)
  nq      <- length(q)
  CBA_marzoD2SR     <- matrix(ncol=nq,nrow=n)
  CBAp_marzoD2SR    <- rep(0,n)
  CBApstd_marzoD2SR <- rep(0,n)
  
  bufferD2SR   <- matrix(ncol=nq,nrow=n)
  descarteD2SR <- matrix(ncol=nq,nrow=n)
  
  for(i in 1:n){
    std              <-read.table(here(carpetaCBA,paste(admb,"1",i,".std",sep="")),header=T,sep="",na="NA",fill=T)
    CBAp_marzoD2SR[i]    <-subset(std,name=="CBA_c0R")$value[1]
    CBApstd_marzoD2SR[i] <-subset(std,name=="CBA_c0R")$std[1]
    for(j in 1:nq){
      CBA_marzoD2SR[i,j] <-qnorm(q[j],CBAp_marzoD2SR[i],CBApstd_marzoD2SR[i])
    }
  }
  
  for(i in 1:n){for(j in 1:nq){
    bufferD2SR[i,j]<-round(1-CBA_marzoD2SR[i,j]/CBA_marzoD2SR[i,5],2)}}
  
  # CBA MENOS DESCARTE ----
  n<-3
  q       <- seq(0.1,0.5,0.1)  # niveles de riesgo (cuantiles)
  nq      <- length(q)
  CBAd_marzoD2SdR     <- matrix(ncol=nq,nrow=n)
  CBApd_marzoD2SdR    <- rep(0,n)
  CBApdstd_marzoD2SdR <- rep(0,n)
  
  
  for(i in 1:n){
    std               <-read.table(here(carpetaCBA,paste(admb,"1",i,".std",sep="")),header=T,sep="",na="NA",fill=T)
    CBApd_marzoD2SdR[i]    <-subset(std,name=="CBA_c0dR")$value[1]
    CBApdstd_marzoD2SdR[i] <-subset(std,name=="CBA_c0dR")$std[1]
    for(j in 1:nq){
      CBAd_marzoD2SdR[i,j] <-qnorm(q[j],CBApd_marzoD2SdR[i],CBApdstd_marzoD2SdR[i])}}
  
  save(list=ls(all=T),
       file=here("Rdata",paste('Datos_Proy_',Hito,'.RData',sep="")))
  
  
  
  
  
}
