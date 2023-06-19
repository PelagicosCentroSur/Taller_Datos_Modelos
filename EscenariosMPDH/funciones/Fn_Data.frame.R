
# ARREGLO DATOS PARA SALIDAS DE TABLAS Y FIGURAS ----

# función datos de composición de edad ----
df_fig8<-function(archivo.Rdata,flota,type){
  load(archivo.Rdata)

  if(flota=='Flota'&type=='observado'){propEdad<-pfobs}
  if(flota=='Flota'&type=='predicho'){propEdad<-pfpred}
  
  if(flota=='Crucero_verano'&type=='observado'){propEdad<-pRobs}
  if(flota=='Crucero_verano'&type=='predicho'){propEdad<-pRpred}
  
  if(flota=='Crucero_otoño'&type=='observado'){propEdad<-pPobs}
  if(flota=='Crucero_otoño'&type=='predicho'){propEdad<-pPpred}
  
  
propE  <- as.data.frame(propEdad) %>% 
         magrittr::set_colnames(age)%>% 
         mutate(yrs=years,
         Asesoria=hito,
         flota=flota,
         type=type) %>% 
        melt(id.vars=c('yrs','Asesoria','flota','type'))
propE
}

# # 8. PUNTOS BIOLÓGICOS DE REFERENCIA ><> ><> ><> ><> ----
# Funcion genera datos para tabla Pasos PBRproxy ----
df_tb4<-function(archivo.Rdata,pSPR_Fmh,pB_Fmh){
  load(archivo.Rdata)
  
  Hito<-c(round(Bmed/10^3,0),
           formatC(round(Fmedian,2), decimal.mark = ","),
           formatC(pSPR_Fmh*100, decimal.mark = ","),
           formatC(60, decimal.mark = ","),
           formatC(pB_Fmh*100, decimal.mark = ","),
           formatC(55, decimal.mark = ","),
           round(Bo/10^3,0),
           round(BRMS/10^3,0),
           round(BLIM/10^3,0))

 return(Hito)
}



# # 10. PROYECCION ><> ><> ><> ><>  ----
# 
# # RESULTADOS DE PROYECCIÓN ASESORÍA DE SEPTIEMBRE (HITO 1) ----
# 
# 
# 
# 
# # Arreglos para figuras proyeccion Fig18 ----
# 
# 
#   datafig18<-function(dir.0,carpetaCBA,admb,yearProy,SSB,BRMS,desembarquepred){
# 
#    # carpetaCBA<-carpetaCBA_sept
#   #  admb<-admb_sept
#     
#     source(paste(dir.fun,"Fn_CBA.R",sep=""))
#     CreaDataProybase(dir.0,carpetaCBA,admb)
#     
#     
#     dira<-paste(dir.0,carpetaCBA,sep="")
#     setwd(dira)
#      load('dataProybase.RData')
#      
# Rtp<-reps1a$Reclutas
# Rs1<-reps1a$Np[1]
# Rs2<-reps2a$Np[1]
# Rs3<-reps3a$Np[1]
# yearsp<-c(reps1a$years,yearProy)
# 
# dataRproy<-data.frame(S1=c(Rtp,rep(Rs1,2)),
#                       S2=c(Rtp,rep(Rs2,2)),
#                       S3=c(Rtp,rep(Rs3,2)),
#                       base=c(Rtp,rep(NA,2))) %>% 
#            mutate(years=yearsp,indicador='Rt') %>% 
#            melt(id.vars=c('years','indicador'))
# 
# 
# dataBDproy<-data.frame(S1=c(SSB,bds1)/1000,
#                        S2=c(SSB,bds2)/1000,
#                        S3=c(SSB,bds3)/1000,
#                        base=c(SSB,rep(NA,2))/1000) %>% 
#             mutate(years=yearsp,indicador='BD')%>% 
#             melt(id.vars=c('years','indicador'))
# 
# 
# dataBD_rmsproy<-data.frame(S1=round(c(SSB,bds1)/BRMS,1),
#                            S2=round(c(SSB,bds2)/BRMS,1),
#                            S3=round(c(SSB,bds3)/BRMS,1),
#                            base=round(c(SSB/BRMS,rep(NA,2)),1)) %>% 
#                 mutate(years=yearsp,indicador='BD_BDrms')%>% 
#                 melt(id.vars=c('years','indicador'))
#  
# dataCt_proy<-data.frame(S1=c(desembarquepred,cs1)/1000,
#                         S2=c(desembarquepred,cs2)/1000,
#                         S3=c(desembarquepred,cs3)/1000,
#                         base=c(desembarquepred,rep(NA,2))/1000) %>% 
#              mutate(years=yearsp,indicador='Ct')%>% 
#              melt(id.vars=c('years','indicador'))
# 
# DataProy<-rbind(dataRproy,dataBDproy,dataBD_rmsproy,dataCt_proy)
# 
# 
# }
# # 
# # 
# # # Arreglos para tablas ----
# # 
# # # RECLUTAMIENTO ESTIMADO ULTIMO AÑO EVALUACIÓN (AÑO ACTUAL) ----
# # RTs0s     <- subset(stds1,name=="Reclutas")$value[nyears1] ; 
# # RTs0s_std  <- subset(stds1,name=="Reclutas")$std[nyears1] 
# # 
# # #BIOMASA DESOVANTE ESTIMADA ULTIMO AÑO EVALUACIÓN ----
# # bds0s     <- subset(stds1,name=="SSB")$value[nyears1] ; 
# # bds0s_std  <- subset(stds1,name=="SSB")$std[nyears1] 
# # 
# # # aporte del grupo de edad 0 (reclutamiento) año actual ----
# # C1eryearR1act<-round(reps1a$YTP_r0W_actual[1]/sum(reps1a$YTP_r0W_actual),2)
# # C1eryearR1act2<-round(reps1a$YTP_r0W_actual[2]/sum(reps1a$YTP_r0W_actual),2)
# # 
# # # aporte del grupo de edad 0 (reclutamiento) 1er año proyectado ----
# # C1eryearR1<-round(reps1a$YTP_p0W_proyectada[1,1]/sum(reps1a$YTP_p0W_proyectada[1,]),2)
# # C1eryearR2<-round(reps2a$YTP_p0W_proyectada[1,1]/sum(reps2a$YTP_p0W_proyectada[1,]),2)
# # C1eryearR3<-round(reps3a$YTP_p0W_proyectada[1,1]/sum(reps3a$YTP_p0W_proyectada[1,]),2)
# # 
# # # aporte del grupo de edad 1 (reclutamiento) 1er año proyectado ----
# # C1eryearR1a<-round(reps1a$YTP_p0W_proyectada[1,2]/sum(reps1a$YTP_p0W_proyectada[1,]),2)
# # C1eryearR2a<-round(reps2a$YTP_p0W_proyectada[1,2]/sum(reps2a$YTP_p0W_proyectada[1,]),2)
# # C1eryearR3a<-round(reps3a$YTP_p0W_proyectada[1,2]/sum(reps3a$YTP_p0W_proyectada[1,]),2)
# # 
# # # aporte del grupo de edad 0 (reclutamiento) 2do año proyectado ----
# # C1eryearR12<-round(reps1a$YTP_p0W_proyectada[2,1]/sum(reps1a$YTP_p0W_proyectada[2,]),2)
# # C1eryearR22<-round(reps2a$YTP_p0W_proyectada[2,1]/sum(reps2a$YTP_p0W_proyectada[2,]),2)
# # C1eryearR32<-round(reps3a$YTP_p0W_proyectada[2,1]/sum(reps3a$YTP_p0W_proyectada[2,]),2)
# # 
# # # aporte del grupo de edad 1  2do año proyectado ----
# # C1eryearR12a<-round(reps1a$YTP_p0W_proyectada[2,2]/sum(reps1a$YTP_p0W_proyectada[2,]),2)
# # C1eryearR22a<-round(reps2a$YTP_p0W_proyectada[2,2]/sum(reps2a$YTP_p0W_proyectada[2,]),2)
# # C1eryearR32a<-round(reps3a$YTP_p0W_proyectada[2,2]/sum(reps3a$YTP_p0W_proyectada[2,]),2)
# # 
# # # ESTATUS PROYECTADO ----
# # 
# # # PRIMER AÑO PROYECTADO ----
# # ### *Probabilidad de estar bajo BRMS* ----
# # pa1<-pnorm(0.9,RpRps1[1],RpRps1std[1],lower.tail = TRUE,log.p = F)
# # pa2<-pnorm(0.9,RpRps2[1],RpRps2std[1],lower.tail = TRUE,log.p = F)
# # pa3<-pnorm(0.9,RpRps3[1],RpRps3std[1],lower.tail = TRUE,log.p = F)
# # ### *Probabilidad de estar en zona de sobreexplotacion* ----
# # pc1<-pnorm(0.9,RpRps1[1],RpRps1std[1],lower.tail = TRUE,log.p = F)-pnorm(0.5,
# #                RpRps1[1],RpRps1std[1],lower.tail = TRUE,log.p = F)
# # pc2<-pnorm(0.9,RpRps2[1],RpRps2std[1],lower.tail = TRUE,log.p = F)-pnorm(0.5,
# #                RpRps2[1],RpRps2std[1],lower.tail = TRUE,log.p = F)
# # pc3<-pnorm(0.9,RpRps3[1],RpRps3std[1],lower.tail = TRUE,log.p = F)-pnorm(0.5,
# #                RpRps3[1],RpRps3std[1],lower.tail = TRUE,log.p = F)
# # ### *Probabilidad de estar en zona de colapso*----
# # pd1<-pnorm(0.5,RpRps1[1],RpRps1std[1],lower.tail = TRUE,log.p = F)
# # pd2<-pnorm(0.5,RpRps2[1],RpRps2std[1],lower.tail = TRUE,log.p = F)
# # pd3<-pnorm(0.5,RpRps3[1],RpRps3std[1],lower.tail = TRUE,log.p = F)
# # 
# # # SEGUNDO AÑO PROYECTADO ----
# # 
# # ### *Probabilidad de estar bajo BRMS* ----
# # pa12<-pnorm(0.9,RpRps1[2],RpRps1std[2],lower.tail = TRUE,log.p = F)
# # pa22<-pnorm(0.9,RpRps2[2],RpRps2std[2],lower.tail = TRUE,log.p = F)
# # pa32<-pnorm(0.9,RpRps3[2],RpRps3std[2],lower.tail = TRUE,log.p = F)
# # ### *Probabilidad de estar en zona de sobreexplotacion* ----
# # pc12<-pnorm(0.9,RpRps1[2],RpRps1std[2],lower.tail = TRUE,log.p = F)-pnorm(0.5,
# #                 RpRps1[2],RpRps1std[2],lower.tail = TRUE,log.p = F)
# # pc22<-pnorm(0.9,RpRps2[2],RpRps2std[2],lower.tail = TRUE,log.p = F)-pnorm(0.5,
# #                 RpRps2[2],RpRps2std[2],lower.tail = TRUE,log.p = F)
# # pc32<-pnorm(0.9,RpRps3[2],RpRps3std[2],lower.tail = TRUE,log.p = F)-pnorm(0.5,
# #                 RpRps3[2],RpRps3std[2],lower.tail = TRUE,log.p = F)
# # ### *Probabilidad de estar en zona de colapso* ----
# # pd12<-pnorm(0.5,RpRps1[2],RpRps1std[2],lower.tail = TRUE,log.p = F)
# # pd22<-pnorm(0.5,RpRps2[2],RpRps2std[2],lower.tail = TRUE,log.p = F)
# # pd32<-pnorm(0.5,RpRps3[2],RpRps3std[2],lower.tail = TRUE,log.p = F)
# # 
# # # CBA INICIAL ----
# # n<-3
# # q       <- seq(0.1,0.5,0.1)  # niveles de riesgo (cuantiles)                                
# # nq      <- length(q)                                                                                   
# # CBA_sept     <- matrix(ncol=nq,nrow=n)
# # CBAp_sept    <- rep(0,n)
# # CBApstd_sept <- rep(0,n)
# # 
# # buffer   <- matrix(ncol=nq,nrow=n)
# # descarte <- matrix(ncol=nq,nrow=n)
# # 
# # for(i in 1:n){
# #   std     <- read.table(paste(admb_sept,"1",i,".std",sep=""),header=T,sep="",na="NA",fill=T) 
# #   CBAp_sept[i]    <-subset(std,name=="CBA_c0")$value[1]
# #   CBApstd_sept[i] <-subset(std,name=="CBA_c0")$std[1]
# #   for(j in 1:nq){CBA_sept[i,j]<-qnorm(q[j],CBAp_sept[i],CBApstd_sept[i])}}
# # 
# # for(i in 1:n){for(j in 1:nq){	
# #   buffer[i,j]<-round(1-CBA_sept[i,j]/CBA_sept[i,5],2)}}
# # 
# # # CBA INICIAL MENOS DESCARTE ----
# # n<-3
# # q       <- seq(0.1,0.5,0.1)  # niveles de riesgo (cuantiles)                                
# # nq      <- length(q)                                                                                   
# # CBAd_sept     <- matrix(ncol=nq,nrow=n)
# # CBApd_sept    <- rep(0,n)
# # CBApdstd_sept <- rep(0,n)
# # 
# # 
# # for(i in 1:n){
# #   std     <- read.table(paste(admb_sept,"1",i,".std",sep=""),header=T,sep="",na="NA",fill=T) 
# #   CBApd_sept[i]    <-subset(std,name=="CBA_c0d")$value[1]
# #   CBApdstd_sept[i] <-subset(std,name=="CBA_c0d")$std[1]
# #   for(j in 1:nq){CBAd_sept[i,j]<-qnorm(q[j],CBApd_sept[i],CBApdstd_sept[i])}}
# # 
# # 
# # 
