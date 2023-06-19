# Tamaño de muestra ----

NM_Ian<-function(dir,admb,system){
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
  
  #---------------------------------------
  # ============================================================================== #
  # I. INDICES DE ABUNDANCIA                                                       #
  # ============================================================================== #
  years  <- data.1$Ind[,1]                                                              
  nyears <- data.1$nanos                                                                
  age    <- seq(0,4,1)                    
  nage   <- data.1$nedades                                                            
  Amax   <- data.1$nedades        
  Age    <- seq(0,4,1) 
  
  
  x  <-c(years,rev(years))
  x1 <-c(years[1],years[nyears]+1,nyears+1/2) #xaxp
  x2 <-c(years[1]-1,years[nyears]+1) #xlim
  
  #============================================================#
  # II. COMPOSICI?N EDAD DE LAS CAPTURAS                       #
  #============================================================#
  #age <-dat$"#Edades"                                         
  age  <-seq(0,4,1)                                            
  nage<-length(age)                                            
  #Proporci?n observada                                        
  pobsF<-rep$pf_obs                                              
  pobsR<-rep$pobs_RECLAS                                       
  pobsP<-rep$pobs_PELACES                                      
  #Proporci?n predicha                                         
  ppredF<-rep$pf_pred                                          
  ppredR<-rep$ppred_RECLAS                                     
  ppredP<-rep$ppred_PELACES                                    
  resfl <-matrix(ncol=nage,nrow=nyears)                         
  for(i in 1:nyears){                                          
    for(j in 1:nage){                                          
      resfl[,j]<-pobsF[,j]-ppredF[,j]}}                        
  #Proporciones                                                
  pF   <- c(pobsF,ppredF); pF[pF==0]  <-NA                     
  pR   <- c(pobsR,ppredR); pR[pR==0]  <-NA                     
  pP   <- c(pobsP,ppredP); pP[pP==0]  <-NA                     
  #arreglos                                                    
  edad <- rep(gl((length(age)),length(years),label=age),2)     
  años <- rep(years,length(age)*2)                             
  ind  <- c(rep("capt_obs",length(years)*length(age)),         
            rep("capt_est",length(years)*length(age)))        
  pro  <- data.frame(años,edad,ind,pF,pR,pP)    
  # ==========================================================================
  
  #=================================================================#
  # M?TODO de Francis
  #=================================================================#
  Nf1 <-60
  Nr1 <-34
  Np1 <-6
  #-------------------------------------------#
  #FLOTA
  fanos<-years
  fobs <-pobsF
  fpre <-ppredF 
  #RECLAS
  ranos<-years
  robs <-pobsR[rowSums(pobsR)>0,]
  rpre <-ppredR[rowSums(pobsR)>0,]
  #PELACES
  panos<-years
  pobs <-pobsP[rowSums(pobsP)>0,]
  ppre <-ppredP[rowSums(pobsP)>0,] 
  
  Of  <- rep(0,length(fanos))
  Ef  <- rep(0,length(fanos))
  vf  <- rep(0,length(fanos))
  vNf <- rep(0,length(fanos))
  
  Or  <- rep(0,length(robs[,1]))
  Er  <- rep(0,length(robs[,1]))
  vr  <- rep(0,length(robs[,1]))
  vNr <- rep(0,length(robs[,1]))
  
  Op  <- rep(0,length(pobs[,1]))
  Ep  <- rep(0,length(pobs[,1]))
  vp  <- rep(0,length(pobs[,1]))
  vNp <- rep(0,length(pobs[,1]))
  #-------------------------------------------#
  for(i in 1:length(fanos)){
    Of[i]  <- sum(fobs[i,]*age)
    Ef[i]  <- sum(fpre[i,]*age)
    vf[i]  <- sum(fpre[i,]*age^2)-Ef[i]^2
    vNf[i] <- vf[i]/Nf1}
  
  for(i in 1:length(robs[,1])){
    Or[i]  <- sum(robs[i,]*age)
    Er[i]  <- sum(rpre[i,]*age)
    vr[i]  <- sum(rpre[i,]*age^2)-Er[i]^2
    vNr[i] <- vr[i]/Nr1}
  
  for(i in 1:length(pobs[,1])){
    Op[i]  <- sum(pobs[i,]*age)
    Ep[i]  <- sum(ppre[i,]*age)
    vp[i]  <- sum(ppre[i,]*age^2)-Ep[i]^2
    vNp[i] <- vp[i]/Np1}
  #--------------------------------------------#
  wf  <- 1/var((Of-Ef)/sqrt(vNf)) 
  wr  <- 1/var((Or-Er)/sqrt(vNr)) 
  wp  <- 1/var((Op-Ep)/sqrt(vNp))
  
  Nf2 <- Nf1*wf                   # NM FLOTA
  Nr2 <- Nr1*wr                   # NM RECLAS
  Np2 <- Np1*wp                   # NM PELACES
  #-----------------------------------------------------------------#
  NM_Fran <- data.frame(nmF=c(Nf1,Nf2),nmR=c(Nr1,Nr2),nmP=c(Np1,Np2))
  #-----------------------------------------------------------------#
  
  #=================================================================#
  # M?todo de Ianelli 2002
  #=================================================================#
  Ofl <-ppredF[rowSums(pobsF)>0,]*(1-ppredF[rowSums(pobsF)>0,])
  Efl <-(pobsF[rowSums(pobsF)>0,]-ppredF[rowSums(pobsF)>0,])^2
  wfl <-rep(0,length(Ofl[,1]))
  for(i in 1:length(Ofl[,1])){
    wfl[i] <-sum(Ofl[i,])/sum(Efl[i,])}
  
  nmf_ari <-mean(wfl)                      # MEDIA ARITMETICA
  nmf_geo <-exp(sum(log(wfl))/length(wfl)) # MEDIA GEOM?TRICA
  nmf_arm <-1/mean(1/wfl)                  # MEDIA ARM?NICA
  
  #------------------------------------------------------------
  Ore <-ppredR[rowSums(pobsR)>0,]*(1-ppredR[rowSums(pobsR)>0,])
  Ere <-(pobsR[rowSums(pobsR)>0,]-ppredR[rowSums(pobsR)>0,])^2
  wre <-rep(0,length(Ore[,1]))
  for(i in 1:length(Ore[,1])){	
    wre[i] <-sum(Ore[i,])/sum(Ere[i,])}
  nmr_ari <-mean(wre)                      # MEDIA ARITMETICA
  nmr_geo <-exp(sum(log(wre))/length(wre)) # MEDIA GEOM?TRICA
  nmr_arm <-1/mean(1/wre)                  # MEDIA ARM?NICA
  #------------------------------------------------------------
  Ope <-ppredP[rowSums(pobsP)>0,]*(1-ppredP[rowSums(pobsP)>0,])
  Epe <-(pobsP[rowSums(pobsP)>0,]-ppredP[rowSums(pobsP)>0,])^2
  wpe <-rep(0,length(Ope[,1]))
  for(i in 1:length(Ope[,1])){
    wpe[i] <-sum(Ope[i,])/sum(Epe[i,])}
  
  nmp_ari <-mean(wpe)                      # MEDIA ARITMETICA
  nmp_geo <-exp(sum(log(wpe))/length(wpe)) # MEDIA GEOM?TRICA
  nmp_arm <-1/mean(1/wpe)                  # MEDIA ARM?NICA
  #------------------------------------------------------------
  #------------------------------------------------------------
  NM_Ian <- data.frame(media=c("Arismética","Geométrica","Armonica"),
                       nmF=c(nmf_ari,nmf_geo,nmf_arm),
                       nmR=c(nmr_ari,nmr_geo,nmr_arm),
                       nmP=c(nmp_ari,nmp_geo,nmp_arm)) 
  
  print(NM_Ian)
}

