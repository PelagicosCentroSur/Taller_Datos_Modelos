
# Tabla índices de abundancia ----
tb1<-function(label_tb,archivo.Rdata){
  
  load(archivo.Rdata)
  
  dataInd<-data.frame(years,reclasobs,pelacesobs,mphobs)
  
  caption1<-"Estimaciones de biomasas  de los cruceros de Verano (RECLAS), 
  Otoño (PELACES) y crucero de huevos (MPDH) utilizadas en la evaluación de stock. "
  
  kbl(dataInd, booktabs = T,format = "latex",position="h!",escape = F,align="c",
      format.args=list(big.mark = '.'),
      col.names = linebreak(c("Año\ncalendario ",
                              "Biomasa crucero\nde verano\n(toneladas)",
                              "Biomasa crucero\nde otoño\n(toneladas)",
                              "Biomasa desovante\nMPDH\n(toneladas)"),align="c"),
      caption = paste(label_tb,caption1,especie,region,sep=" ")) %>% 
    kable_styling(latex_options = c("striped", "condensed"),
                  full_width = FALSE,font_size=10)
}

# tabla capturas ----
tb2<-function(label_tb,bioyear,desembarque,porcDesc_actualizado,CapturaDescartada,CapturaTotal,especie){
  
  dataDes_y_descar<-data.frame("Año"=bioyear,
                               "Desembarques"=round(desembarque,0),
                               "Pdescarte"=porcDesc_actualizado,
                               "Capturadesc"=round(CapturaDescartada,0),
                               "Capturatotal"=round(CapturaTotal,0))
  
  caption1<-"Desembarques en toneladas, porcentaje de descarte supuesto,
captura descartada (toneladas) y captura total (toneladas) estimadas en año biológico. "
  
  
  kbl(dataDes_y_descar, booktabs = T,format = "latex",position="h!",align="c",
      format.args=list(big.mark = '.'),escape=F,
      col.names = linebreak(c("Año\nbiológico ",
                              "Desembarques\n(toneladas)",
                              "Porcentaje\nDescarte",
                              "Captura\ndescartada\n(toneladas)",
                              "Captura\ntotal\n(toneladas)"),align="c"),   
      caption = paste(label_tb,caption1,especie,region,sep=" ")) %>% 
    kable_styling(latex_options = c("striped", "condensed"),
                  full_width = FALSE,font_size=10)
}

# tabla variables poblacionales ----
tb3<-function(label_tb,Rdata_H1,Rdata_H2){
  
  load(Rdata_H1); VarPob_H1<-Var
  load(Rdata_H2); VarPob_H2<-Var
  
  VarPobl1<- cbind(yearsbiol,
                   format(round(c(VarPob_H1[VarPob_H1$indicador=='BD',]$value,NA),0),nsmall=0, big.mark="."),
                   format(round(VarPob_H2[VarPob_H2$indicador=='BD',]$value,0),nsmall=0, big.mark="."),
                   format(round(c(VarPob_H1[VarPob_H1$indicador=='BT',]$value,NA),0),nsmall=0, big.mark="."),
                   format(round(VarPob_H2[VarPob_H2$indicador=='BT',]$value,0),nsmall=0, big.mark="."),
                   format(round(c(VarPob_H1[VarPob_H1$indicador=='Rt',]$value,NA),0),nsmall=0, big.mark="."),
                   format(round(VarPob_H2[VarPob_H2$indicador=='Rt',]$value,0),nsmall=0, big.mark="."),
                   formatC(round(c(VarPob_H1[VarPob_H1$indicador=='Ft',]$value,NA),3),decimal.mark = ",",digits = 3),
                   formatC(round(VarPob_H2[VarPob_H2$indicador=='Ft',]$value,3),decimal.mark = ",",digits = 3))
  
  caption1<-"\\label{Tab17}Estimaciones medias de las biomasa desovante (t), 
  biomasa total (t), reclutamientos (\\#) y  mortalidad por pesca (año-1)  
  estimadas en las asesorías de septiembre 2022 (Hito 1) y marzo 2023 (Hito 2). "
  # TABLA
  VarPobl1 %>%
  kbl(booktabs = T,format = "latex",position="h!",escape = F,align="c",
      col.names = linebreak(c("biológico",
                              "Hito 1","Hito 2","Hito 1","Hito 2","Hito 1","Hito 2","Hito 1","Hito 2"),align="c"),
  caption = paste(label_tb,caption1,especie,region,sep=" ")) %>%
    kable_styling(latex_options = c("striped", "condensed"),
                  full_width = FALSE,font_size=8)%>%
    add_header_above(c('Año'=1,
                        "Biomasa\ndesovante"=2,
                        "Biomasa total"=2,
                        "Reclutamientos"=2,
                        "Mortalidad\npor pesca"=2))
}

# tabla puntos biológicos de referencia ----
tb4<-function(label_tb,pbrHito1,pbrHito2,especie){
  
  etapas<-c("Paso 1", "Paso 1", "Paso 2","Paso 2","Paso 3", "Paso 3", "Paso 4", "Paso 5", "Paso 6")
  
  pbrs<-c("$BD_{promedio}$",
          "$F_{mh}$",
          "$\\%BDPR_{F_{mh}}$",
          "$\\%BDPR_{F_{RMS}}$",
          "$\\%BD_{F_{mh}}$",
          "$\\%BD_{F_{RMS}}$",
          "$BD_{0}$",
          "$BD_{55\\%}$",
          "$BD_{27,5\\%}$")
  
  calculo<-c("Promedio de la serie histórica",
             "Mediana de la serie histórica",
             rep("Cálculo de la curva de biomasa por recluta (BDPR)",2),
             "$\\%BDPR_{F_{mh}}-5\\%$",
             "$\\%BDPR_{F_{RMS}}-5\\%$",
             "$BD_{0}=BD_{promedio}/\\%BD_{F_{mh}}$",
             "$BD_{RMS}=BD_{0}*\\%BD_{F_{RMS}}$",
             "$BD_{LIM}=BD_{0}*\\%BD_{F_{LIM}}$")
 
  tablapbrs<-data.frame(etapas,calculo,pbrs,pbrHito1,pbrHito2)
  caption1 <- "\\label{Tab18}Procedimiento de cálculo de los puntos biológicos de
  referencia de biomasa (miles de t) estimados en las asesorías de septiembre 
  2022 (Hito 1) y marzo 2023 (Hito 2), calculados siguiendo los pasos descritos 
  en la metodología de este informe. "
  
  #-----------
  # TABLA 
  #-----------
  kbl(tablapbrs, booktabs = T,format = "latex",position="h!",align="cllc",escape = FALSE,
      col.names = linebreak(c('Etapas','Cálculo','Aproximación','Hito 1','Hito 2'),align="c"),
      caption = paste(label_tb,caption1,especie,region,sep=" "))%>%
    collapse_rows(columns = 1:2,latex_hline = "major", valign = "middle")%>% 
    kable_styling(latex_options = c("striped", "condensed"),
                  full_width = FALSE,font_size=8)%>% 
    add_header_above(c(" " = 3,"Asesorías" = 2))
}

tb5<-function(label_tb,pbrs,especie){
#attach(pbrs)  
 
Tabla4.1<-rbind("Año biológico"=pbrs$yearStatus,
                "$F_{RMS}$ $(año^{-1})$"=formatC(round(pbrs$FRMS,2), decimal.mark = ","),
                "$BD_{RMS}$ (mil t)"=round(pbrs$BRMS/10^3,0),
                "$BD_{LIM}$ (mil t)"=round(pbrs$BLIM/10^3,0),
                "$p(BD_{last}<BD_{RMS})^1$"=formatC(round(pbrs$pa,3), decimal.mark = ","),
                "$p(F_{last}>F_{RMS})^2$"=formatC(round(pbrs$pb,3), decimal.mark = ","),
                "$p(sobre-explotación)^3$"=formatC(round(pbrs$pc,3), decimal.mark = ","),
                "$p(agotado/colapsado)^4$"=formatC(round(pbrs$pd,3), decimal.mark = ","),
                "$p(sobrepesca)^5$"=formatC(round(pbrs$pe,3), decimal.mark = ","))

colnames(Tabla4.1)<-pbrs$Hitos

footnote<-c("Probabilidad que $BD$ del año más reciente sea menor a $BD_{RMS}$ según el diagrama de fase",
            "Probabilidad que $F$ del año más reciente sea mayor a $F_{RMS}$ según el diagrama de fase",
            "Probabilidad de estar en sobreexplotación = $p(0,5<BD_{last}/BD_{RMS}<0,9)$",
            "Probabilidad de estar en colapso =$p(BD_{last}/BD_{RMS}<0,5)$",
            "Probabilidad de estar en sobrepesca = $p(F_{last}/F_{RMS}>1,1)$")

# TABLA 

caption1<-"\\label{Tab21}Puntos Biológicos de referencia (PBRs) y 
probabilidades de estar bajo $BD_{RMS}$ y sobre $F_{RMS}$ y en sobreexplotación, 
colapsado o sobrepesca. Estimaciones realizadas en asesorías de septiembre (Hito 1) 
y marzo (Hito 2). "

kbl(Tabla4.1, booktabs = T,format = "latex",position="h!",align="c",escape = FALSE,
    caption = paste(label_tb,caption1,especie,region,sep=" ")) %>% 
  kable_styling(latex_options = c("striped"),
                full_width = FALSE,font_size=9) %>%
  footnote(number = footnote,threeparttable = T,escape = FALSE,fixed_small_size = TRUE) %>%
  column_spec(1, width = "25em") %>% 
  add_header_above(c(" " = 1,"Asesorías" = 2))

}

# Tabla indicadores de estatus

tb6<-function(label_tb,Rdata_H1,Rdata_H2){

  load(Rdata_H1); IndStatus_H1<-IndStatus
  load(Rdata_H2); IndStatus_H2<-IndStatus
  
  Estatus2<- data.frame(
    'Años'=yearsbiol,
    "F_FRMS_h1"=formatC(round(c(IndStatus_H1[IndStatus_H1$indicador=='F_FRMS',]$value,NA),3), decimal.mark = ","),
    "F_FRMS_h2"=formatC(round(c(IndStatus_H2[IndStatus_H2$indicador=='F_FRMS',]$value),3), decimal.mark = ","),
    "BD_BDRMS_h1"=formatC(round(c(IndStatus_H1[IndStatus_H1$indicador=='BD_BDRMS',]$value,NA),3), decimal.mark = ","),
    "BD_BDRMS_h2"=formatC(round(c(IndStatus_H2[IndStatus_H2$indicador=='BD_BDRMS',]$value),3), decimal.mark = ","),
    "Y_BT_h1"=formatC(round(c(IndStatus_H1[IndStatus_H1$indicador=='Y_BT',]$value,NA),3), decimal.mark = ","),
    "Y_BT_h2"=formatC(round(c(IndStatus_H2[IndStatus_H2$indicador=='Y_BT',]$value),3), decimal.mark = ","),
    "C_N_h1"=formatC(round(c(IndStatus_H1[IndStatus_H1$indicador=='C_N',]$value,NA),3), decimal.mark = ","),
    "C_N_h2"=formatC(round(c(IndStatus_H2[IndStatus_H2$indicador=='C_N',]$value),3), decimal.mark = ","))
  
  

caption1 <- "\\label{Tab19}Índices de reducción de F respecto de $F_{RMS}$ 
($F/F_{RMS}$), BD respecto de $BD_{RMS}$ ($BD/BD_{RMS}$), tasas de explotación
anual referidos a la biomasa ($Y/BT$) y a la abundancia estimada ($C\\#/N\\#$), 
estimadas en las asesorías de septiembre (Hito 1) y marzo (Hito 2). "

kbl(Estatus2, booktabs = T,format = "latex",position="h!",align="c",escape = F,
    col.names = linebreak(c('biológico','Hito 1','Hito 2','Hito 1','Hito 2','Hito 1','Hito 2','Hito 1','Hito 2'),align="c"),
    caption = paste(label_tb,caption1,especie,region,sep=" ")) %>%
  kable_styling(latex_options = c("striped"),
                full_width = FALSE,font_size=8) %>%
  add_header_above(c("Año" = 1,
                     "$F/F_{RMS}$"=2,
                     "$BD/BD_{RMS}$"=2,
                     "$Y/BT$"=2,
                     "$C\\\\#/N\\\\#$"=2), escape = F)

}

# Tablas proyección hito 1 ----

tb7<-function(label_tb,Rdata_proy_Hito,escRecl,especie,region){
  
  load(Rdata_proy_Hito)
  
  probEstatus<-rbind("$R_{proy}$"=formatC(c(reps1a$Np[1],
                                            reps2a$Np[1],
                                            reps3a$Np[1],
                                            reps1a$Np[1],
                                            reps2a$Np[1],
                                            reps3a$Np[1])/1000, decimal.mark = ","),
                     "$F_{proy}=F_{RMS}$"=formatC(c(rep(FRMSp,6)), decimal.mark = ","),
                     "$BD_{RMS}$"=formatC(c(rep(BRMSp,6))/1000, decimal.mark = ","),
                     '$BD_{proy}$'=formatC(round(c(bds1[1],
                                                   bds2[1],
                                                   bds3[1],
                                                   bds1[2],
                                                   bds2[2],
                                                   bds3[2])/1000,0), decimal.mark = ","),
                     '$BD_{proy}/BD_{RMS}$'=formatC(round(c(RpRps1[1],
                                                            RpRps2[1],
                                                            RpRps3[1],
                                                            RpRps1[2],
                                                            RpRps2[2],
                                                            RpRps3[2]),2), decimal.mark = ","),
                     "$p(sobreexplotación)^1$"=formatC(round(c(pc1,
                                                               pc2,
                                                               pc3,
                                                               pc12,
                                                               pc22,
                                                               pc32),2), decimal.mark = ","),
                     "$p(agotado/colapsado)^2$"=formatC(round(c(pd1,
                                                                pd2,
                                                                pd3,
                                                                pd12,
                                                                pd22,
                                                                pd32),2), decimal.mark = ","))
  
  colnames(probEstatus)<-rep(escRecl,2)
  
  footnote2<-c("Probabilidad de estar en sobreexplotación = $p(0,5<BD_{last}/BD_{RMS}<0,9)$",
               "Probabilidad de estar en colapso =$p(BD_{last}/BD_{RMS}<0,5)$")
  
  caption1 <- "\\label{Tab25s}Proyección de la biomasa desovante hacia el año biológico 2022/23 y 
      2023/24 considerando tres escenarios de reclutamientos y mortalidad por pesca igual 
      al $F_{RMS}$. Probabilidad de estar bajo $BD_{RMS}$ y sobre $F_{RMS}$ y la probabilidad 
      de sobre-explotación y colapso."
  #-----------
  # TABLA 
  #-----------
  kbl(probEstatus, booktabs = T,format = "latex",position="h!",align="c",escape = FALSE,
      caption = paste(label_tb,caption1,especie,region,sep=" ")) %>% 
    kable_styling(latex_options = c("striped", "condensed"),
                  full_width = FALSE,font_size=10)%>% 
    add_header_above(c(" " = 1,
                       "Escenarios de reclutamiento" = 3,
                       "Escenarios de reclutamiento" = 3),line = F)%>% 
    add_header_above(c(" " = 1,
                       "Estatus 2022/23" = 3,"Estatus 2023/24" = 3),line = F)%>%
    footnote(number = footnote2,threeparttable = T,escape = FALSE,fixed_small_size = TRUE)
}



tb8<-function(label_tb,Rdata_proy_Hito,escRecl,especie,region){
  
  load(Rdata_proy_Hito)
  
  Capturaproy<-rbind('Año biológico proyectada 2022/2023'=round(c(cs1[1],cs2[1],cs3[1])/1000,0),
                     'Primer Semestre 2023'=round(c(cs1[1]*0.7,cs2[1]*0.7,cs3[1]*0.7)/1000,0),
                     'Año biológico proyectada 2023/2024'=round(c(cs1[2],cs2[2],cs3[2])/1000,0),
                     'Segundo Semestre 2023'=round(c(cs1[2]*0.3,cs2[2]*0.3,cs3[2]*0.3)/1000,0),
                     'Año calendario 2023'=round(c((cs1[1]*0.7)+(cs1[2]*0.3),
                                                   (cs2[1]*0.7)+(cs2[2]*0.3),
                                                   (cs3[1]*0.7)+(cs3[2]*0.3))/1000,0))
  
  colnames(Capturaproy)<-escRecl
  
  caption1 <-"\\label{Tab26s}Proyección de la captura estimada para los años biológicos 
      2022/23 y 2023/24 considerando una mortalidad por pesca igual al $F_{RMS}$, 
      con sus respectivos escenarios de reclutamientos. El cálculo de la Captura 
      para el año calendario 2023 se obtiene como el promedio ponderado según la 
      estacionalidad semestral de la pesquería que a la fecha se asume 70\\% para el primer 
      semestre y 30\\% para el segundo semestre."
  
  #-----------
  # TABLA 
  #-----------
  kbl(Capturaproy, booktabs = T,format = "latex",position="h!",align="c",escape = FALSE,
      caption = paste(label_tb,caption1,especie,region,sep=" ")) %>% 
    kable_styling(latex_options = c("striped", "condensed"),
                  full_width = FALSE,font_size=10)%>% 
    add_header_above(c("Captura" = 1,
                       "Escenarios de reclutamiento" = 3),line = F)
}


tb9<-function(label_tb,Rdata_proy_Hito,escRecl,especie,region){
  
  load(Rdata_proy_Hito)
  
  percentiles<-c("10\\%","20\\%","30\\%","40\\%","50\\%")
  
  tCBA<-t(data.frame(CBAp_sept,CBApstd_sept,CBA_sept))
  colnames(tCBA)  <-escRecl
  row.names(tCBA) <-c('mean','sd',percentiles)
  
  caption1 <-"\\label{Tab27s}CBA inicial año calendario 2023 calculada bajo 
      $F_{RMS}$ con sus respectivos percentiles de captura entre 10\\% y 50\\% y 
      tres escenarios de reclutamientos."
  #-----------
  # TABLA 
  #-----------
  kbl(round(tCBA,0), booktabs = T,format = "latex",position="h!",align="c",escape = FALSE,
      format.args=list(big.mark = '.'),
      caption = paste(label_tb,caption1,especie,region,sep=" ")) %>% 
    kable_styling(latex_options = c("striped", "condensed"),
                  full_width = FALSE,font_size=10)%>% 
    add_header_above(c(" " = 1,
                       "Escenarios de reclutamiento" = 3),line = F)%>% 
    add_header_above(c(" " = 1,
                       "CBA 2023" = 3))
  
}


tb10_anch<-function(label_tb,Rdata_proy_Hito,escRecl,especie,region){
  
  load(Rdata_proy_Hito)
  percentiles<-c("10\\%","20\\%","30\\%","40\\%","50\\%")
  
  tCBAd<-t(data.frame(CBApd_sept,CBApdstd_sept,CBAd_sept))
  colnames(tCBAd)<-escRecl
  row.names(tCBAd)<-c('mean','sd',percentiles)
  
  caption1 <-"\\label{Tab27s_b}CBA inicial año calendario 2023 considerando el
              descuento del 2\\% de descarte."
  
  #-----------
  # TABLA 
  #-----------
  kbl(round(tCBAd,0), booktabs = T,format = "latex",position="h!",escape = FALSE,align="c",
      format.args=list(big.mark = '.'),
      caption = paste(label_tb,caption1,especie,region,sep=" ")) %>% 
    kable_styling(latex_options = c("striped", "condensed"),
                  full_width = FALSE,font_size=10)%>% 
    add_header_above(c(" " = 1,
                       "Escenarios de reclutamiento" = 3),line = F)%>% 
    add_header_above(c(" " = 1,
                       "CBA 2023 - 2\\% descarte" = 3))
  
}


tb11<-function(label_tb,Rdata_proy_Hito,escRecl,especie,region){
  
  load(Rdata_proy_Hito)
  percentiles<-c("10\\%","20\\%","30\\%","40\\%","50\\%")
  
  buffer2<-data.frame(t(buffer))
  colnames(buffer2)<-escRecl
  row.names(buffer2)<-percentiles
  
  caption1 <-"\\label{Tab28s}Nivel de resguardo de cada percentil de Captura al RMS."
  
  #----------
  # TABLA 
  #-----------
  kbl(buffer2, booktabs = T,format = "latex",position="h!",align="c",escape = FALSE,
      format.args=list(decimal.mark = ','),
      caption = paste(label_tb,caption1,especie,region,sep=" ")) %>% 
    kable_styling(latex_options = c("striped", "condensed"),
                  full_width = FALSE,font_size=10)%>% 
    add_header_above(c(" " = 1,
                       "Escenarios de reclutamiento" = 3),line = F)%>% 
    add_header_above(c(" " = 1,
                       "Resguardo" = 3))
  
  
}




# Tablas proyección hito 2 ----

tb7_H2<-function(label_tb,Rdata,Rdata_proy_Hito,escRecl,especie,region){
  
  load(Rdata);pc;pd
  
  load(Rdata_proy_Hito)
  
  probEstatus<-rbind("$R_{actual/proy}$"=formatC(c(Rtp[nyears],
                                            reps1a$Np[1],
                                            reps2a$Np[1],
                                            reps3a$Np[1])/1000, decimal.mark = ","),
                     "$F_{supuesto/proy}=F_{RMS}$"=formatC(round(c(rep(FRMSp,4)),2), decimal.mark = ","),
                     "$BD_{RMS}$"=formatC(c(rep(BRMSp,4))/1000, decimal.mark = ","),
                     '$BD_{actual/proy}$'=formatC(round(c(SSBp[nyears],
                                                   bds1[1],
                                                   bds2[1],
                                                   bds3[1])/1000,0), decimal.mark = ","),
                     '$BD_{actual/proy}/BD_{RMS}$'=formatC(round(c(SSBp[nyears]/BRMSp,
                                                            RpRps1[1],
                                                            RpRps2[1],
                                                            RpRps3[1]),2), decimal.mark = ","),
                     "$p(sobreexplotación)^1$"=formatC(round(c(pc,
                                                               pc1,
                                                               pc2,
                                                               pc3),2), decimal.mark = ","),
                     "$p(agotado/colapsado)^2$"=formatC(round(c(pd,
                                                                pd1,
                                                                pd2,
                                                                pd3),2), decimal.mark = ","))
  
  colnames(probEstatus)<-c("$R_{2023}$",escRecl)
  
  footnote2<-c("Probabilidad de estar en sobreexplotación = $p(0,5<BD_{last}/BD_{RMS}<0,9)$",
               "Probabilidad de estar en colapso =$p(BD_{last}/BD_{RMS}<0,5)$")
  
  caption1 <- "\\label{Tab25m} Estatus año 2022/23 actualizado en Hito 2 (marzo 2023) y 
      Proyección de la biomasa desovante hacia el año biológico 2023/24,
      considerando reclutamiento actualizado en Hito 2 y tres escenarios de reclutamientos promedio
      ($R_{actual/proy}$) y mortalidad por pesca actual y proyectada ($F_{supuesto/proy}$) igual 
      al $F_{RMS}$. Biomasa desovante objetivo ($BD_{RMS}$), actualizada en Hito 2 y 
      proyectada ($BD_{actual/proy}$) la probabilidad 
      de sobre-explotación y colapso."
  #-----------
  # TABLA 
  #-----------
  kbl(probEstatus, booktabs = T,format = "latex",position="h!",align="c",escape = FALSE,
      caption = paste(label_tb,caption1,especie,region,sep=" ")) %>% 
    kable_styling(latex_options = c("striped", "condensed"),
                  full_width = FALSE,font_size=10)%>% 
    add_header_above(c(" " = 1,
                       "R actual" = 1,
                       "Escenarios de reclutamiento" = 3),line = F)%>% 
    add_header_above(c(" " = 1,
                       "Estatus actual 2022/23" = 1,"Estatus proyectado 2023/24" = 3),line = F)%>%
    footnote(number = footnote2,threeparttable = T,escape = FALSE,fixed_small_size = TRUE)
}



tb8_H2<-function(label_tb,Rdata_proy_Hito,escRecl,especie,region){
  
  load(Rdata_proy_Hito)
  
  Capturaproy<-rbind('Año biológico proyectada 2022/23'=c(round((cs0[1])/1000,0),"-","-","-"),
                     'Primer Semestre 2023'=c(round((cs0[1]*0.7)/1000,0),"-","-","-"),
                     'Año biológico proyectada 2023/24'=c("-",round(c(cs1[1],cs2[1],cs3[1])/1000,0)),
                     'Segundo Semestre 2023'=c("-",round(c(cs1[1]*0.3,cs2[1]*0.3,cs3[1]*0.3)/1000,0)),
                     'Año calendario 2023'=c("-",round(c((cs0[1]*0.7)+(cs1[1]*0.3),
                                                   (cs0[1]*0.7)+(cs2[1]*0.3),
                                                   (cs0[1]*0.7)+(cs3[1]*0.3))/1000,0)))
  
  colnames(Capturaproy)<-c("$R_{2023}$",escRecl)
  
  caption1 <-"\\label{Tab26m}Proyección de la captura estimada para los años biológicos 
      2022/23 y 2023/24 considerando una mortalidad por pesca igual al $F_{RMS}$, 
      con sus respectivos escenarios de reclutamientos. El cálculo de la Captura 
      para el año calendario 2023 se obtiene como el promedio ponderado según la 
      estacionalidad semestral de la pesquería que a la fecha se asume 70\\% para el primer 
      semestre y 30\\% para el segundo semestre."
  
  #-----------
  # TABLA 
  #-----------
  kbl(Capturaproy, booktabs = T,format = "latex",position="h!",align="c",escape = FALSE,
      caption = paste(label_tb,caption1,especie,region,sep=" ")) %>% 
    kable_styling(latex_options = c("striped", "condensed"),
                  full_width = FALSE,font_size=10)%>% 
    add_header_above(c("Captura" = 1,
                       "R actual" = 1,
                       "Escenarios de reclutamiento" = 3),line = F)
}

#  Considera el descuento del desembarque  segundo  semestre año previo y  remanente  de cuota año previo

tb8b_H2<-function(label_tb,Rdata_proy_Hito,escRecl,especie,region){
  
  load(Rdata_proy_Hito)
  
  Capturaproy<-rbind('Año biológico proyectada 2022/23'=c(round((cs0[1])/1000,0),"-","-","-"),
                     'Desembarque segundo semestre 2022'=c(round(desem2doSem/1000,0),"-","-","-"),
                     'Saldo  cuota 2022'=c(round(remanente/1000,0),"-","-","-"),
                     'Primer Semestre 2023'=c(round((cs0R[1])/1000,0),"-","-","-"),
                     'Año biológico proyectada 2023/24'=c("-",round(c(cs1[1],cs2[1],cs3[1])/1000,0)),
                     'Segundo Semestre 2023'=c("-",round(c(cs1[1]*0.3,cs2[1]*0.3,cs3[1]*0.3)/1000,0)),
                     'Año calendario 2023'=c("-",round(c((cs0R[1])+(cs1[1]*0.3),
                                                         (cs0R[1])+(cs2[1]*0.3),
                                                         (cs0R[1])+(cs3[1]*0.3))/1000,0)))
  
  colnames(Capturaproy)<-c("$R_{2023}$",escRecl)

  
  caption1 <-"\\label{Tab26m}Captura año 2022/23 actualizada en Hito 2 (marzo 2023) y proyectada hacia
                 el año 2023/24 considerando  una mortalidad por pesca al $F_{RMS}$, con sus
                 respectivos escenarios de reclutamientos. El cálculo de la captura para el año 
                 calendario 2023 descontando el desembarque del segundo semestre 2022 y el saldo de
                 cuota 2022 a la captura del año biológico 2022/23 y el 30\\% de la captura proyectada
                 hacia el año 2022/23. "
  
  #-----------
  # TABLA 
  #-----------
  kbl(Capturaproy, booktabs = T,format = "latex",position="h!",align="c",escape = FALSE,
      caption = paste(label_tb,caption1,especie,region,sep=" ")) %>% 
    kable_styling(latex_options = c("striped", "condensed"),
                  full_width = FALSE,font_size=10)%>% 
    add_header_above(c("Captura" = 1,
                       "R actual" = 1,
                       "Escenarios de reclutamiento" = 3),line = F)
}


tb9_H2<-function(label_tb,Rdata_proy_Hito,escRecl,especie,region){
  
  load(Rdata_proy_Hito)
  
  percentiles<-c("10\\%","20\\%","30\\%","40\\%","50\\%")
  
  tCBA<-t(data.frame(CBAp_marzo,CBApstd_marzo,CBA_marzo))
  colnames(tCBA)  <-escRecl
  row.names(tCBA) <-c('mean','sd',percentiles)
  
  caption1 <-"\\label{Tab27m}Primera revisión de la CBA año calendario 2023 calculada bajo 
      $F_{RMS}$ con sus respectivos percentiles de captura entre 10\\% y 50\\% y 
      tres escenarios de reclutamientos."
  #-----------
  # TABLA 
  #-----------
  kbl(round(tCBA,0), booktabs = T,format = "latex",position="h!",align="c",escape = FALSE,
      format.args=list(big.mark = '.'),
      caption = paste(label_tb,caption1,especie,region,sep=" ")) %>% 
    kable_styling(latex_options = c("striped", "condensed"),
                  full_width = FALSE,font_size=10)%>% 
    add_header_above(c(" " = 1,
                       "Escenarios de reclutamiento" = 3),line = F)%>% 
    add_header_above(c(" " = 1,
                       "CBA 2023" = 3))
  
}


tb10_anch_H2<-function(label_tb,Rdata_proy_Hito,escRecl,especie,region){
  
  load(Rdata_proy_Hito)
  percentiles<-c("10\\%","20\\%","30\\%","40\\%","50\\%")
  
  tCBAd<-t(data.frame(CBApd_marzo,CBApdstd_marzo,CBAd_marzo))
  colnames(tCBAd)<-escRecl
  row.names(tCBAd)<-c('mean','sd',percentiles)
  
  caption1 <-"\\label{Tab27m_b}Primera revisión de la CBA año calendario 2023 considerando el
              descuento del 2\\% de descarte."
  
  #-----------
  # TABLA 
  #-----------
  kbl(round(tCBAd,0), booktabs = T,format = "latex",position="h!",escape = FALSE,align="c",
      format.args=list(big.mark = '.'),
      caption = paste(label_tb,caption1,especie,region,sep=" ")) %>% 
    kable_styling(latex_options = c("striped", "condensed"),
                  full_width = FALSE,font_size=10)%>% 
    add_header_above(c(" " = 1,
                       "Escenarios de reclutamiento" = 3),line = F)%>% 
    add_header_above(c(" " = 1,
                       "CBA 2023 - 2\\% descarte" = 3))
  
}


tb11_H2<-function(label_tb,Rdata_proy_Hito,escRecl,especie,region){
  
  load(Rdata_proy_Hito)
  percentiles<-c("10\\%","20\\%","30\\%","40\\%","50\\%")
  
  buffer2<-data.frame(t(buffer))
  colnames(buffer2)<-escRecl
  row.names(buffer2)<-percentiles
  
  caption1 <-"\\label{Tab28m}Nivel de resguardo de cada percentil de Captura al RMS."
  
  #----------
  # TABLA 
  #-----------
  kbl(buffer2, booktabs = T,format = "latex",position="h!",align="c",escape = FALSE,
      format.args=list(decimal.mark = ','),
      caption = paste(label_tb,caption1,especie,region,sep=" ")) %>% 
    kable_styling(latex_options = c("striped", "condensed"),
                  full_width = FALSE,font_size=10)%>% 
    add_header_above(c(" " = 1,
                       "Escenarios de reclutamiento" = 3),line = F)%>% 
    add_header_above(c(" " = 1,
                       "Resguardo" = 3))
  
  
}

# considerando el descuento del desembarque 2do semestre y saldo cuota año previo

tb9B_H2<-function(label_tb,Rdata_proy_Hito,escRecl,especie,region){
  
  load(Rdata_proy_Hito)
  
  percentiles<-c("10\\%","20\\%","30\\%","40\\%","50\\%")
  
  tCBA<-t(data.frame(CBAp_marzoD2SR,CBApstd_marzoD2SR,CBA_marzoD2SR))
  colnames(tCBA)  <-escRecl
  row.names(tCBA) <-c('mean','sd',percentiles)
  
  caption1 <-"\\label{Tab27m1_c2} Primera revisión de la CBA año calendario 2023 calculada bajo 
      $F_{RMS}$ con sus respectivos percentiles de captura entre 10\\% y 50\\% y 
      tres escenarios de reclutamientos."
  #-----------
  # TABLA 
  #-----------
  kbl(round(tCBA,0), booktabs = T,format = "latex",position="h!",align="c",escape = FALSE,
      format.args=list(big.mark = '.'),
      caption = paste(label_tb,caption1,especie,region,sep=" ")) %>% 
    kable_styling(latex_options = c("striped", "condensed"),
                  full_width = FALSE,font_size=10)%>% 
    add_header_above(c(" " = 1,
                       "Escenarios de reclutamiento" = 3),line = F)%>% 
    add_header_above(c(" " = 1,
                       "CBA 2023" = 3))
  
}


tb10B_anch_H2<-function(label_tb,Rdata_proy_Hito,escRecl,especie,region){
  
  load(Rdata_proy_Hito)
  percentiles<-c("10\\%","20\\%","30\\%","40\\%","50\\%")
  
  tCBAd<-t(data.frame(CBApd_marzoD2SdR,CBApdstd_marzoD2SdR,CBAd_marzoD2SdR))
  colnames(tCBAd)<-escRecl
  row.names(tCBAd)<-c('mean','sd',percentiles)
  
  caption1 <-"\\label{Tab27m2_c2} Primera revisión de la CBA año calendario 2023 considerando el
              descuento del 2\\% de descarte."
  
  #-----------
  # TABLA 
  #-----------
  kbl(round(tCBAd,0), booktabs = T,format = "latex",position="h!",escape = FALSE,align="c",
      format.args=list(big.mark = '.'),
      caption = paste(label_tb,caption1,especie,region,sep=" ")) %>% 
    kable_styling(latex_options = c("striped", "condensed"),
                  full_width = FALSE,font_size=10)%>% 
    add_header_above(c(" " = 1,
                       "Escenarios de reclutamiento" = 3),line = F)%>% 
    add_header_above(c(" " = 1,
                       "CBA 2023 - 2\\% descarte" = 3))
  
}


tb11B_H2<-function(label_tb,Rdata_proy_Hito,escRecl,especie,region){
  
  load(Rdata_proy_Hito)
  percentiles<-c("10\\%","20\\%","30\\%","40\\%","50\\%")
  
  buffer2<-data.frame(t(bufferD2SR))
  colnames(buffer2)<-escRecl
  row.names(buffer2)<-percentiles
  
  caption1 <-"\\label{Tab28m_c2} Nivel de resguardo de cada percentil de Captura al RMS."
  
  #----------
  # TABLA 
  #-----------
  kbl(buffer2, booktabs = T,format = "latex",position="h!",align="c",escape = FALSE,
      format.args=list(decimal.mark = ','),
      caption = paste(label_tb,caption1,especie,region,sep=" ")) %>% 
    kable_styling(latex_options = c("striped", "condensed"),
                  full_width = FALSE,font_size=10)%>% 
    add_header_above(c(" " = 1,
                       "Escenarios de reclutamiento" = 3),line = F)%>% 
    add_header_above(c(" " = 1,
                       "Resguardo" = 3))
  
  
}


tb12_anch_H2<-function(label_tb,Rdata_proy_Hito1,Rdata_proy_Hito2,escRecl,especie,region){
  
  percentiles<-c("10\\%","20\\%","30\\%","40\\%","50\\%")
  load(Rdata_proy_Hito1)
  CBAd_sept
  
  
  load(Rdata_proy_Hito2)
  CBAd_marzo
  
  
  tabladiferencias<--round(1-(CBAd_marzo/CBAd_sept),2)*100
  colnames(tabladiferencias)<-percentiles
  row.names(tabladiferencias)<-escRecl
  
  
  
  caption1 <-"\\label{Tab27m_c1} Diferencias porcentuales entre 
  la CBA 2023 proyectada en Hito 1 y  CBA 2023 actualizada en el Hito 2 (Criterio 1)."
  
  #-----------
  # TABLA 
  #-----------
  kbl(round(t(tabladiferencias),0), booktabs = T,format = "latex",position="h!",escape = FALSE,align="c",
      format.args=list(big.mark = '.'),
      caption = paste(label_tb,caption1,especie,region,sep=" ")) %>% 
    kable_styling(latex_options = c("striped", "condensed"),
                  full_width = FALSE,font_size=10)%>% 
    add_header_above(c(" " = 1,
                       "Escenarios de reclutamiento" = 3),line = F)%>% 
    add_header_above(c(" " = 1,
                       "CBA 2023 - 2\\% descarte" = 3))
  
}

tb13_anch_H2<-function(label_tb,Rdata_proy_Hito1,Rdata_proy_Hito2,escRecl,especie,region){
  
  percentiles<-c("10\\%","20\\%","30\\%","40\\%","50\\%")
  load(Rdata_proy_Hito1)
  CBAd_sept
  
  
  load(Rdata_proy_Hito2)
  CBAd_marzoD2SdR
  
  
  tabladiferencias<--round(1-(CBAd_marzoD2SdR/CBAd_sept),2)*100
  colnames(tabladiferencias)<-percentiles
  row.names(tabladiferencias)<-escRecl
  
  
  
  caption1 <-"\\label{Tab27m3_c2} Diferencias porcentuales entre 
  la CBA 2023 proyectada en Hito 1 y  CBA 2023 actualizada en el Hito 2 (Criterio 2)."
  
  #-----------
  # TABLA 
  #-----------
  kbl(round(t(tabladiferencias),0), booktabs = T,format = "latex",position="h!",escape = FALSE,align="c",
      format.args=list(big.mark = '.'),
      caption = paste(label_tb,caption1,especie,region,sep=" ")) %>% 
    kable_styling(latex_options = c("striped", "condensed"),
                  full_width = FALSE,font_size=10)%>% 
    add_header_above(c(" " = 1,
                       "Escenarios de reclutamiento" = 3),line = F)%>% 
    add_header_above(c(" " = 1,
                       "CBA 2023 - 2\\% descarte" = 3))
  
}


# escenarios de supuesto de captura  Hito 2

tb14<-function(label_tb,comp,especie){
  #attach(pbrs)  
  
  Tabla4.1<-rbind("Año biológico"=comp$yearStatus,
                  "$F_{2022/23}$ $(año^{-1})$"=formatC(round(comp$Flast,2), decimal.mark = ","),
                  "$BD_{2022/23}$ (mil t)"=round(comp$BDlast/10^3,0),
                  "$F_{2022/23}/F_{RMS}$ "=formatC(round(comp$F_Frms,2), decimal.mark = ","),
                  "$BD_{2022/23}/BD_{RMS}$ "=formatC(round(comp$BD_BDrms,2), decimal.mark = ","),
                  "$p(BD_{2022/23}<BD_{RMS})^1$"=formatC(round(comp$pa,3), decimal.mark = ","),
                  "$p(F_{2022/23}>F_{RMS})^2$"=formatC(round(comp$pb,3), decimal.mark = ","),
                  "$p(sobre-explotación)^3$"=formatC(round(comp$pc,3), decimal.mark = ","),
                  "$p(agotado/colapsado)^4$"=formatC(round(comp$pd,3), decimal.mark = ","),
                  "$p(sobrepesca)^5$"=formatC(round(comp$pe,3), decimal.mark = ","))
  
  colnames(Tabla4.1)<-comp$Hitos
  
  footnote<-c("Probabilidad que $BD$ del año más reciente sea menor a $BD_{RMS}$ según el diagrama de fase",
              "Probabilidad que $F$ del año más reciente sea mayor a $F_{RMS}$ según el diagrama de fase",
              "Probabilidad de estar en sobreexplotación = $p(0,5<BD_{last}/BD_{RMS}<0,9)$",
              "Probabilidad de estar en colapso =$p(BD_{last}/BD_{RMS}<0,5)$",
              "Probabilidad de estar en sobrepesca = $p(F_{last}/F_{RMS}>1,1)$")
  
  # TABLA 
  
  caption1<-"\\ Comparación entre el caso base y escenarios alternativos de 
  captura del año biológico actuallas estimaciones de mortalidad por pesca, 
  biomasa desovante  del año biológico actual y 
probabilidades de estar bajo $BD_{RMS}$ y sobre $F_{RMS}$ y en sobreexplotación, 
colapsado o sobrepesca. "
  
  kbl(Tabla4.1, booktabs = T,format = "latex",position="h!",align="c",escape = FALSE,
      caption = paste(label_tb,caption1,especie,region,sep=" ")) %>% 
    kable_styling(latex_options = c("striped"),
                  full_width = FALSE,font_size=9) %>%
    footnote(number = footnote,threeparttable = T,escape = FALSE,fixed_small_size = TRUE) %>%
    column_spec(1, width = "25em") %>% 
    add_header_above(c(" " = 1,"Escenarios" = 5))
  
}

#comparación CBA criterio 1
tb10_anch_C1<-function(label_tb,escRecl,especie,region){
  
  
  percentiles<-c("10\\%","20\\%","30\\%","40\\%","50\\%")
  escenarios <-c(rep("S1",5),rep("S2",5),rep("S3",5),rep("S4",5),rep("Caso base",5))
  
  load(Rdata_proy_Hito2)
  CBAd_marzoS0 <- t(CBAd_marzo)
  load(Rdata_proy_Hito2S1)
  CBAd_marzoS1 <- t(CBAd_marzo)
  load(Rdata_proy_Hito2S2)
  CBAd_marzoS2 <- t(CBAd_marzo)
  load(Rdata_proy_Hito2S3)
  CBAd_marzoS3 <- t(CBAd_marzo)
  load(Rdata_proy_Hito2S4)
  CBAd_marzoS4 <- t(CBAd_marzo)
  
  
  tCBAd<-data.frame(escenarios,percentiles,
                    format(round(rbind(CBAd_marzoS1,CBAd_marzoS2,
                                       CBAd_marzoS3,CBAd_marzoS4,CBAd_marzoS0),0),big.mark = '.'))
  
  
  colnames(tCBAd)<-c("Escenarios","Percentiles",escRecl)
  #row.names(tCBAd)<-escenarios
  
  caption1 <-"Estimación de CBA 2023 de los cuatro escenarios de supuesto de captura 2022/23
              considerando el 
              descuento del 2\\% de descarte. Criterio 1."
  
  
  #-----------
  # TABLA 
  #-----------
  kbl(tCBAd, booktabs = T,format = "latex",position="h!",escape = FALSE,align="c",
      caption = paste(label_tb,caption1,especie,region,sep=" ")) %>%
    kable_styling(latex_options = c("striped", "condensed"),
                  full_width = FALSE,font_size=10)%>%
    add_header_above(c(" " = 2,
                       "Escenarios de reclutamiento" = 3),line = F)%>%
    add_header_above(c(" " = 2,
                       "CBA 2023 - 2\\% descarte" = 3))
  
}

#comparación CBA criterio 2
tb10_anch_C2<-function(label_tb,escRecl,especie,region){
  
  
  percentiles<-c("10\\%","20\\%","30\\%","40\\%","50\\%")
  escenarios <-c(rep("S1",5),rep("S2",5),rep("S3",5),rep("S4",5),rep("Caso base",5))
  
  load(Rdata_proy_Hito2)
  CBAd_marzoS0 <- t(CBAd_marzoD2SdR)
  load(Rdata_proy_Hito2S1)
  CBAd_marzoS1 <- t(CBAd_marzoD2SdR)
  load(Rdata_proy_Hito2S2)
  CBAd_marzoS2 <- t(CBAd_marzoD2SdR)
  load(Rdata_proy_Hito2S3)
  CBAd_marzoS3 <- t(CBAd_marzoD2SdR)
  load(Rdata_proy_Hito2S4)
  CBAd_marzoS4 <- t(CBAd_marzoD2SdR)
  
  
  tCBAd<-data.frame(escenarios,percentiles,
                    format(round(rbind(CBAd_marzoS1,CBAd_marzoS2,
                                       CBAd_marzoS3,CBAd_marzoS4,CBAd_marzoS0),0),big.mark = '.'))
  
  
  colnames(tCBAd)<-c("Escenarios","Percentiles",escRecl)
  #row.names(tCBAd)<-escenarios
  
  caption1 <-"Estimación de CBA 2023 de los cuatro escenarios de supuesto de captura 2022/23
              considerando el 
              descuento del 2\\% de descarte. Criterio 2."
  
  
  #-----------
  # TABLA 
  #-----------
  kbl(tCBAd, booktabs = T,format = "latex",position="h!",escape = FALSE,align="c",
      caption = paste(label_tb,caption1,especie,region,sep=" ")) %>%
    kable_styling(latex_options = c("striped", "condensed"),
                  full_width = FALSE,font_size=10)%>%
    add_header_above(c(" " = 2,
                       "Escenarios de reclutamiento" = 3),line = F)%>%
    add_header_above(c(" " = 2,
                       "CBA 2023 - 2\\% descarte" = 3))
  
}

