#><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
# Función DataHito ----
#><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
Rdata_base<-function(codigos,admb_dat,admb_rep,admb_std,Especie,Region,Hito,MesHito,dir.Rdata){
  #librerias requeridas----
  paquetes <- c("stringr", "tidyverse", "kableExtra","ggplot2",
                "ggthemes","patchwork","dplyr","reshape","here")
  lapply(paquetes, require, character.only = TRUE)
  #><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
  especie     <- Especie
  region      <- Region
  hito        <- Hito
  admb_base   <- str_sub(admb_dat,1,-5)
  meshito     <- MesHito
  data        <- lisread(here(codigos,admb_dat)) 
  names(data) <- str_trim(names(data), side="right")
  dat         <- data
  rep         <- reptoRlist(here(codigos,admb_rep))
  std         <- read.table(here(codigos,admb_std),header=T,sep="",na="NA",fill=T) 
  #datos para Rdata ----
  years  <- rep$years
  nyears   <- length(years)
  yearsbiol      <-paste(years-1,"/",str_sub(years, 3, 4),sep="")
  yearbiol_last  <-paste(years[nyears-1],"/",str_sub(years[nyears],3,4),sep="")
  yearbiol_proy  <-c(paste(years[nyears],"/",str_sub(years[nyears]+1,3,4),sep=""),
                     paste(years[nyears]+1,"/",str_sub(years[nyears]+2,3,4),sep=""))
  
  age     <- seq(0,4,1)                                            
  nage    <- length(age)  
  #pesos medios ----
  WmedF    <- dat$Wmed                                             
  WiniF    <- dat$Wini  
  
  #><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
  #indices ----
  #><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
  names_ind<- c('Crucero_verano', 
                'Crucero_otoño',
                'Crucero_huevos', 
                'Desembarques') 
  #observados 
  reclasobs     <-rep$reclasobs %>% na_if(0)
  pelacesobs    <-rep$pelacesobs %>% na_if(0)
  mphobs        <-rep$mphobs %>% na_if(0)
  desembarqueobs<-rep$desembarqueobs %>% na_if(0)
          #data.frame observados 
          indobs  <- data.frame(
                     reclasobs,
                     pelacesobs,
                     mphobs,
                     desembarqueobs) %>% 
                     #na_if(0) %>% 
                     magrittr::set_colnames(names_ind) %>%
                     mutate(Asesoria=hito,type='observado',yrs= years) %>% 
                     melt(id.var=c('yrs','type', 'Asesoria'))  
  #predichos
  reclaspred     <-rep$reclaspred 
  pelacespred    <-rep$pelacespred 
  mphpred        <-rep$mphpred
  desembarquepred<-rep$desembarquepred
          #data.frame predichos
          indpred  <- data.frame(
                      reclaspred,
                      pelacespred,
                      mphpred,
                      desembarquepred) %>% 
                      #na_if(0) %>% 
                      magrittr::set_colnames(names_ind) %>%
                      mutate(Asesoria=hito,type='predicho',yrs= years) %>% 
                      melt(id.var=c('yrs','type', 'Asesoria')) 
  #><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
  #Cvs indices ----
  #><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
  cvBcV   <-0.30
  cvBcO   <-0.30
  cvdes   <-0.01
  
  #><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
  #composiciones edad ----
  #><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
  #observados
  pfobs<-rep$pf_obs
  pRobs<-rep$pobs_RECLAS
  pPobs<-rep$pobs_PELACES
  #predichos
  pfpred<-rep$pf_pred
  pRpred<-rep$ppred_RECLAS
  pPpred<-rep$ppred_PELACES
  
  #><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
  #Variables poblacionales ----
  #><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
  
  Rt      <- subset(std,name=="Reclutas")$value 
  Rtstd   <- subset(std,name=="Reclutas")$std
  BT      <- subset(std,name=="BT")$value   
  BTstd   <- subset(std,name=="BT")$std
  BD      <- subset(std,name=="SSB")$value   
  BDstd   <- subset(std,name=="SSB")$std
  Ft      <- subset(std,name=="log_Ft")$value   
  Ftstd   <- subset(std,name=="log_Ft")$std
    #data.frame Variables poblacionales
    rt<- data.frame(x=years,
                    y=Rt,
                    lower = (Rt-1.96*Rtstd),
                    upper = (Rt+1.96*Rtstd))%>% 
                    mutate(indicador='Rt')
    bt<- data.frame(x=years,
                    y=BT,
                    lower = (BT-1.96*BTstd),
                    upper = (BT+1.96*BTstd))%>% 
                    mutate(indicador='BT')
    bd<- data.frame(x=years,
                    y=BD,
                    lower = (BD-1.96*BDstd),
                    upper = (BD+1.96*BDstd))%>% 
                    mutate(indicador='BD')
    ft<- data.frame(x=years,
                    y=exp(Ft),
                    lower = exp(Ft-1.96*Ftstd),
                    upper = exp(Ft+1.96*Ftstd))%>% 
                    mutate(indicador='Ft')
    Var<-rbind(rt,bt,bd,ft)%>%
        mutate(Hito=hito)%>% 
        melt(id.vars=c('x','Hito','indicador','lower','upper'))
    
 #><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> 
  #Selectividades ----
 #><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
   
  sel_Flota<-rep$S_f[1,]
  sel_CruV <-rep$Scru_reclas[1,]
  sel_CruO <-rep$Scru_pelaces[1,]
  
  if(especie=="Sardina común"){
  sel_Flota<-rep$Sel_flota[1,]
  sel_CruV <-rep$Sel_reclas[1,]
  sel_CruO <-rep$Sel_pelaces[1,]
  }

  
  #><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
  # Puntos biológicos de referencia ----
  #><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
  
  Amax        <- dat$nedades
  Fmort 	    <- seq(0,3.5,0.02)
  nf          <- length(Fmort)
  R0 		      <- 1
  #datos de entrada
  Dat<-list()
  Dat$M		    <- dat$par[5]
  Dat$Tspw	  <- dat$Dt[3]
  Dat$Mad	    <- dat$madurezsexual
  Dat$Wmed	  <- colMeans(dat$Wmed)
  Dat$Wini	  <- colMeans(dat$Wini)
  Dat$Sel     <- sel_Flota
  Rmed        <-mean(Rt,na.rm = T) 
  Bmed        <- mean(BD,na.rm = T)         
  Fmedian     <- exp(median(Ft,na.rm = T))
  Fstatuquo   <- rep$Ftot[nyears]
  Bo           <- rep$SSBpbr[1]        # Paso 4: Obtenci?n de Bo
  BRMS         <- rep$SSBpbr[3]        # Paso 5: Obtenci?n de Brms = 60%SPRo = 55%Bo
  FRMS         <- rep$Fs[2]
  BLIM         <- Bo*0.275             # Paso 6: Obtenci?n de Blim = 20%Bo 
  FLIM         <- rep$Fs[3]            # Paso 6: Obtenci?n de Flim = 30%SPRo
  SpB          <- BD                   # BD serie hist?rica de evaluaci?n de stock 
  SpBSE        <- BDstd                # desviaci?n estandar BD
  ln_Fyr       <- Ft                   # logaritmo de Ft
  ln_FSE       <- Ftstd                # logaritmo de la desviaci?n standar de Ft
  
  #><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
  # Indicadores de Estatus ----
  #><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
  
  RPR     <- subset(std,name=="RPRequ3")$value
  RPRstd  <- subset(std,name=="RPRequ3")$std
  FRPR    <- subset(std,name=="Frpr")$value
  FRPRstd <- subset(std,name=="Frpr")$std
  Y_BT    <- desembarquepred/BT
  C_N     <- rowSums(rep$pred_Ctot)/rowSums(rep$N)
   #data.frame indicadores de estatus
    Rpr<- data.frame(x=years,
                     y=RPR,
                     lower = (RPR-1.96*RPRstd),
                     upper = (RPR+1.96*RPRstd))%>% 
                     mutate(indicador='BD_BDRMS')
    Frpr<- data.frame(x=years,
                      y=FRPR,
                      lower = (FRPR-1.96*FRPRstd),
                      upper = (FRPR+1.96*FRPRstd))%>% 
                      mutate(indicador='F_FRMS')
    y_bt<-data.frame(x=years,
                     y=Y_BT,
                     lower = rep(NA,nyears),
                     upper = rep(NA,nyears))%>% 
                     mutate(indicador='Y_BT')
    c_n<-data.frame(x=years,
                    y=C_N,
                    lower = rep(NA,nyears),
                    upper = rep(NA,nyears))%>% 
                    mutate(indicador='C_N')
    
    IndStatus<-rbind(Rpr,Frpr,y_bt,c_n)%>%
               mutate(Hito=hito)%>% 
               melt(id.vars=c('x','Hito','indicador','lower','upper'))
  
  #><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
  # Probabilidades de estatus ----
  #><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
  
  #Para densidad de probabilidad último año----
  rprlast     <-subset(std,name=="RPRequ3")$value[nyears]
  rprlaststd  <-subset(std,name=="RPRequ3")$std[nyears]
  Frprlast    <-subset(std,name=="Frpr")$value[nyears]
  Frprlaststd <-subset(std,name=="Frpr")$std[nyears]
  # biomasa desovante vs BDrms----
  xb1 <-rnorm(1000, mean = rprlast, sd = rprlaststd)
  xb  <-seq(min(xb1),max(xb1),0.005)
  yb  <-dnorm(xb, mean = rprlast, sd =rprlaststd)
  icb <-qnorm(c(0.05,0.95,0.5),rprlast,rprlaststd)
  # mortalidad por pesca vs Frms----
  xf1 <- rnorm(1000, mean = Frprlast, sd = Frprlaststd)
  xf  <-seq(min(xf1),max(xf1),0.005)
  yf  <-dnorm(xf, mean = Frprlast, sd =Frprlaststd)
  icf <-qnorm(c(0.05,0.95,0.5),Frprlast,Frprlaststd)
  #distribución probabilidad----
  xxb<- c(xb[xb>=icb[1]&xb<=icb[2]],rev(xb[xb>=icb[1]&xb<=icb[2]]))
  yyb<- c(yb[xb>=icb[1]&xb<=icb[2]],rep(0,length(yb[xb>=icb[1]&xb<=icb[2]])))
  xxf<- c(xf[xf>=icf[1]&xf<=icf[2]],rev(xf[xf>=icf[1]&xf<=icf[2]]))
  yyf<- c(yf[xf>=icf[1]&xf<=icf[2]],rep(0,length(yf[xf>=icf[1]&xf<=icf[2]])))

      #dataframe densidad de probabilidad indicadores de estatus
      densb_b  <- data.frame(x=xxb, y=yyb)%>% 
                  mutate(t='a',
                         r=seq(1,length(xxb),1),
                         indicador='BD_BDRMS')
      
      densb_f  <- data.frame(x=xxf, y=yyf)%>% 
                  mutate(t='a',
                         r=seq(1,length(xxf),1),
                         indicador='F_FRMS')
      densb<-rbind(densb_b,densb_f)%>%
        mutate(Hito=hito)
  
      densb_xy <- data.frame(x=xb,y=yb)%>% 
                 mutate(t='a',
                 indicador='BD_BDRMS')
      densf_xy <- data.frame(x=xf,y=yf)%>% 
                 mutate(t='a',
                 indicador='F_FRMS')
      
      densbf<-rbind(densb_xy,densf_xy)%>%
        mutate(Hito=hito)
      
  # *Probabilidad de estar bajo BRMS* #Asesoría  #P(BD<BDrms)---- 
  pa <-pnorm(0.9,rprlast,rprlaststd,lower.tail = TRUE,log.p = F)
  # *Probabilidad de estar bajo FRMS* #Asesoría  #P(F>Frms)----
  pb <-1-pnorm(1.1,Frprlast,Frprlaststd,lower.tail = TRUE,log.p = F)
  # *Probabilidad de estar en zona de sobreexplotacion*  #Asesoría  #P(BD<BDrms)---- 
  pc <-pnorm(0.9,rprlast,rprlaststd,lower.tail = TRUE,log.p = F)-pnorm(0.5,
                 rprlast,rprlaststd,lower.tail = TRUE,log.p = F)
  # *Probabilidad de estar en zona de colapso* #Asesoría  #P(BD<BDrms) ----
  pd <-pnorm(0.5,rprlast,rprlaststd,lower.tail = TRUE,log.p = F)
  # *Probailidad de sobrepesca* #Asesoría  #P(F>Frms)----
  pe <-1-pnorm(1.1,Frprlast,Frprlaststd,lower.tail = TRUE,log.p = F)
  #guarda Rdata ----
  save(list=ls(all=T),
       file=paste(dir.Rdata,'/Datos_',Hito,'.RData',sep=""))
  
}



# # 4. COMPARACIÓN CON ASESORÍAS PREVIAS ><> ><> ><> ><> ----
Rdata_retrohist<-function(archivo.Rdata,dir,dir.Rdata){
  load(archivo.Rdata)
  
  rep_mar21  <- reptoRlist(here(dir,"MAE0321.rep"))
  rep_jul21  <- reptoRlist(here(dir,"MAE0721.rep"))
  rep_sept21 <- reptoRlist(here(dir,"MAE0921.rep"))
  rep_mar22  <- reptoRlist(here(dir,"MAE322.rep"))
  rep_jul22  <- reptoRlist(here(dir,"MAE722.rep"))
  rep_sept22 <- reptoRlist(here(dir,"MAE922.rep"))
  rep_mar23  <- reptoRlist(here(dir,"MAE323.rep"))
  
  Rtcomp <- data.frame(x=years,
                       Rt_mar21=c(rep_mar21$Reclutas,NA,NA),
                       Rt_jul21=c(rep_jul21$Reclutas,NA,NA),
                       Rt_sept21=c(rep_sept21$Reclutas,NA,NA),
                       Rt_mar22=c(rep_mar22$Reclutas,NA),
                       Rt_jul22=c(rep_jul22$Reclutas,NA),
                       Rt_sept22=c(rep_sept22$Reclutas,NA),
                       Rt_mar23=rep_mar23$Reclutas,
                       Rt_jul23=rep$Reclutas)
  
  SSBtcomp <- data.frame(x=years,
                         SSBt_mar21=c(rep_mar21$SSB,NA,NA),
                         SSBt_jul21=c(rep_jul21$SSB,NA,NA),
                         SSBt_sept21=c(rep_sept21$SSB,NA,NA),
                         SSBt_mar22=c(rep_mar22$SSB,NA),
                         SSBt_jul22=c(rep_jul22$SSB,NA),
                         SSBt_sept22=c(rep_sept22$SSB,NA),
                         SSBt_mar23=rep_mar23$SSB,
                         SSBt_jul23=rep$SSB)
  
  Ftcomp <- data.frame(x=years,
                       Ft_mar21=c(rep_mar21$Ftot,NA,NA),
                       Ft_jul21=c(rep_jul21$Ftot,NA,NA),
                       Ft_sept21=c(rep_sept21$Ftot,NA,NA),
                       Ft_mar22=c(rep_mar22$Ftot,NA),
                       Ft_jul22=c(rep_jul22$Ftot,NA),
                       Ft_sept22=c(rep_sept22$Ftot,NA),
                       Ft_mar23=rep_mar23$Ftot,
                       Ft_jul23=rep$Ftot)
  
  year_retros <- c('Hito_3 (2023)','Hito_2 (2023)','Hito_1 (2022)','Hito_3 (2022)','Hito_2 (2022)','Hito_1 (2021)',
                   "Hito_3 (2021)","Hito_2 (2021)")
  nretros <-8
  
  save(list=ls(all=T),
       file=paste("Rdata",'/RetrosHist_',Hito,'.RData',sep=""))
}

# # 6. ANÁLISIS RETROSPECTIVO ><> ><> ><> ><> ----
Rdata_retroanalitico<-function(archivo.Rdata,dirRDretroA,dir.Rdata){


  load(archivo.Rdata);  Hito<-hito
  
listretros2<-list.files(dirRDretroA,pattern=".rep")  
retros2  <- seq(2,length(listretros2))
nretros2 <- length(retros2)
#year_retros2<-as.factor(years[(nyears-(nretros2-1)):nyears])

retroR      <- matrix(0,nrow=nyears,ncol=nretros2+1)
retroBD     <- matrix(0,nrow=nyears,ncol=nretros2+1)
retroF      <- matrix(0,nrow=nyears,ncol=nretros2+1)

#setwd(dir)
for(i in 1:nretros2){
  x<-retros2[i]
  rep2<- reptoRlist(paste(dirRDretroA,listretros2[x],sep="/"))
  retroR[,i+1] <- c(rep2$Reclutas,rep(NA,i-1))
  retroBD[,i+1] <- c(rep2$SSB,rep(NA,i-1))
  retroF[,i+1]  <- c(rep2$Ftot,rep(NA,i-1)) 
  }

# retrospectivo relativo (cálculo)
mohn.r       <- rep(NA, nretros2)
rel.diff.r   <- matrix(NA, nrow=nyears, ncol=(nretros2))
mohn.ssb     <- rep(NA, nretros2)
rel.diff.ssb <- matrix(NA, nrow=nyears, ncol=(nretros2))
mohn.f       <- rep(NA, nretros2)
rel.diff.f   <- matrix(NA, nrow=nyears, ncol=(nretros2))

for(j in 1:nretros2){
  rel.diff.r[,j]   <- (retroR[,(j+1)]-retroR[,2])/retroR[,2]
  mohn.r[j]        <- rel.diff.r[(nyears-j),j]
  rel.diff.ssb[,j] <- (retroBD[,(j+1)]-retroBD[,2])/retroBD[,2]
  mohn.ssb[j]      <- rel.diff.ssb[(nyears-j),j]
  rel.diff.f[,j]   <- (retroF[,(j+1)]-retroF[,2])/retroF[,2]
  mohn.f[j]        <- rel.diff.f[(nyears-j),j]}

ave.mohn.r    <- mean(mohn.r)
ave.mohn.ssb  <- mean(mohn.ssb)
ave.mohn.f    <- mean(mohn.f)

# Arreglo datos

#Para retrospectivo tradicional
Rt_retro<- data.frame(x=years,
                      y1=retroR[,2],
                      y2=retroR[,3],
                      y3=retroR[,4],
                      y4=retroR[,5],
                      y5=retroR[,6])
BD_retro<- data.frame(x=years,
                      y1=retroBD[,2],
                      y2=retroBD[,3],
                      y3=retroBD[,4],
                      y4=retroBD[,5],
                      y5=retroBD[,6])
Ft_retro<- data.frame(x=years,
                      y1=retroF[,2],
                      y2=retroF[,3],
                      y3=retroF[,4],
                      y4=retroF[,5],
                      y5=retroF[,6])

#Para restrospectivo relativo
Rt_retroRel<- data.frame(x=years,
                         y1=rel.diff.r[,1],
                         y2=rel.diff.r[,2],
                         y3=rel.diff.r[,3],
                         y4=rel.diff.r[,4],
                         y5=rel.diff.r[,5])
BD_retroRel<- data.frame(x=years,
                         y1=rel.diff.ssb[,1],
                         y2=rel.diff.ssb[,2],
                         y3=rel.diff.ssb[,3],
                         y4=rel.diff.ssb[,4],
                         y5=rel.diff.ssb[,5])
Ft_retroRel<- data.frame(x=years,
                         y1=rel.diff.f[,1],
                         y2=rel.diff.f[,2],
                         y3=rel.diff.f[,3],
                         y4=rel.diff.f[,4],
                         y5=rel.diff.f[,5])

year_retros <- as.factor(years[(nyears-4):nyears])
nretros <- length(year_retros)

save(list=ls(all=T),
     file=paste(dir.Rdata,'/RetroAnalitico_',Hito,'.RData',sep=""))

}


# # 5. PERFIL VEROSIMILITUD ><> ><> ><> ><> ----
Rdata_perfil<-function(dir,admb,dir.Rdata,Hito){
setwd(dir)

casos <-23
logRo    <- rep(0,casos)
likeval  <- matrix(ncol=15,nrow=casos)
slikeval <- matrix(ncol=16,nrow=casos)

#
for(i in 1:casos){
  report      <- reptoRlist(paste(dir,admb,i,".rep",sep=""))
  logRo[i]    <- report$log_Ro
  likeval[i,] <- report$likeval}

like    <- data.frame(round(likeval,3),Total=apply(likeval,1,sum))
minLik  <- apply(like,2,min)

# busca el mínimo
for(i in 1:16){
  slikeval[,i]<-like[,i]-minLik[i]
}    # Estandarización

names<-c("Ro","Bio_Reclas","Bio_Pelaces","Desembarques","Bio_Mph","C.Edad_Flota",
         "C.Edad_Recl","C.Edad_Pel","prepPelTall","DesvRt","qreclas","qpela","PenFt",
         "PenFspr","NA","NA","Total")
# Tabla verosimilitud
TLk1 <- data.frame(exp(logRo),like);colnames(TLk1)<-names
# Tabla estandarizada
TLk2<- data.frame(exp(logRo),slikeval);colnames(TLk2)<-names

Ro_reclas  <- TLk2$Ro[TLk2$Bio_Reclas==0]
Ro_pelaces <- TLk2$Ro[TLk2$Bio_Pelaces==0]
Ro_desemb  <- min(TLk2$Ro[TLk2$Desembarques==0])
Ro_MPDH    <- ifelse(sum(TLk2$Bio_Mph)>0,0,TLk2$Ro[TLk2$Bio_Mph==0])
Ro_propF   <- TLk2$Ro[TLk2$C.Edad_Flota==0]
Ro_propRecl<- TLk2$Ro[TLk2$C.Edad_Recl==0]
Ro_propPel <- TLk2$Ro[TLk2$C.Edad_Pel==0]
Ro         <- TLk2$Ro[TLk2$Total==0]

names_res<-c("Bio_Reclas","Bio_Pelaces","Desembarques","Bio_MPDH","C.Edad_Flota","C.Edad_Reclas","C.Edad_Pelaces")
res<-c((Ro_reclas-Ro),
       (Ro_pelaces-Ro),
       (Ro_desemb-Ro),
       (Ro_MPDH-Ro),
       (Ro_propF-Ro),
       (Ro_propRecl-Ro),
       (Ro_propPel-Ro))

residuos<-res/Ro
#datares<-data.frame(names_res,residuos)
save(list=ls(all=T),
     file=paste(dir.Rdata,'/PerfilVero_',Hito,'.RData',sep=""))

}
# 

Rdata_retroanalitico2<-function(Hito,dirRDretroA,dir.Rdata){
  

  
  listretros2<-list.files(dirRDretroA,pattern=".rep")  
  retros2  <- seq(2,length(listretros2))
  nretros2 <- length(retros2)
  #year_retros2<-as.factor(years[(nyears-(nretros2-1)):nyears])
  
  retroR      <- matrix(0,nrow=nyears,ncol=nretros2+1)
  retroBD     <- matrix(0,nrow=nyears,ncol=nretros2+1)
  retroF      <- matrix(0,nrow=nyears,ncol=nretros2+1)
  
  #setwd(dir)
  for(i in 1:nretros2){
    x<-retros2[i]
    rep2<- reptoRlist(paste(dirRDretroA,listretros2[x],sep="/"))
    retroR[,i+1] <- c(rep2$Reclutas,rep(NA,i-1))
    retroBD[,i+1] <- c(rep2$SSB,rep(NA,i-1))
    retroF[,i+1]  <- c(rep2$Ftot,rep(NA,i-1)) 
  }
  
  # retrospectivo relativo (cálculo)
  mohn.r       <- rep(NA, nretros2)
  rel.diff.r   <- matrix(NA, nrow=nyears, ncol=(nretros2))
  mohn.ssb     <- rep(NA, nretros2)
  rel.diff.ssb <- matrix(NA, nrow=nyears, ncol=(nretros2))
  mohn.f       <- rep(NA, nretros2)
  rel.diff.f   <- matrix(NA, nrow=nyears, ncol=(nretros2))
  
  for(j in 1:nretros2){
    rel.diff.r[,j]   <- (retroR[,(j+1)]-retroR[,2])/retroR[,2]
    mohn.r[j]        <- rel.diff.r[(nyears-j),j]
    rel.diff.ssb[,j] <- (retroBD[,(j+1)]-retroBD[,2])/retroBD[,2]
    mohn.ssb[j]      <- rel.diff.ssb[(nyears-j),j]
    rel.diff.f[,j]   <- (retroF[,(j+1)]-retroF[,2])/retroF[,2]
    mohn.f[j]        <- rel.diff.f[(nyears-j),j]}
  
  ave.mohn.r    <- mean(mohn.r)
  ave.mohn.ssb  <- mean(mohn.ssb)
  ave.mohn.f    <- mean(mohn.f)
  
  # Arreglo datos
  
  #Para retrospectivo tradicional
  Rt_retro<- data.frame(x=years,
                        y1=retroR[,2],
                        y2=retroR[,3],
                        y3=retroR[,4],
                        y4=retroR[,5],
                        y5=retroR[,6])
  BD_retro<- data.frame(x=years,
                        y1=retroBD[,2],
                        y2=retroBD[,3],
                        y3=retroBD[,4],
                        y4=retroBD[,5],
                        y5=retroBD[,6])
  Ft_retro<- data.frame(x=years,
                        y1=retroF[,2],
                        y2=retroF[,3],
                        y3=retroF[,4],
                        y4=retroF[,5],
                        y5=retroF[,6])
  
  #Para restrospectivo relativo
  Rt_retroRel<- data.frame(x=years,
                           y1=rel.diff.r[,1],
                           y2=rel.diff.r[,2],
                           y3=rel.diff.r[,3],
                           y4=rel.diff.r[,4],
                           y5=rel.diff.r[,5])
  BD_retroRel<- data.frame(x=years,
                           y1=rel.diff.ssb[,1],
                           y2=rel.diff.ssb[,2],
                           y3=rel.diff.ssb[,3],
                           y4=rel.diff.ssb[,4],
                           y5=rel.diff.ssb[,5])
  Ft_retroRel<- data.frame(x=years,
                           y1=rel.diff.f[,1],
                           y2=rel.diff.f[,2],
                           y3=rel.diff.f[,3],
                           y4=rel.diff.f[,4],
                           y5=rel.diff.f[,5])
  
  year_retros <- as.factor(years[(nyears-4):nyears])
  nretros <- length(year_retros)
  
  save(list=ls(all=T),
       file=paste(dir.Rdata,'/RetroA_',Hito,'.RData',sep=""))
  
}
