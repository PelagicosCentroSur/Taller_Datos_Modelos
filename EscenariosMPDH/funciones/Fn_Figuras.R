
# FIGURAS ----

# DESCRIPCION DE DATOS ----
# Índices de abundancia
fig1<-function(archivo.Rdata){
  load(archivo.Rdata)
  
    
  obsC  <- as.data.frame(reclasobs/1000) %>% 
    mutate(year=years) %>% 
    melt(id.vars='year') %>% 
    mutate(type='2.Cruceros de Verano')
  obsP  <- as.data.frame(pelacesobs/1000) %>% 
    mutate(year=years) %>% 
    melt(id.vars='year') %>% 
    mutate(type='3.Cruceros de Otoño')
  obsD  <- as.data.frame(desembarqueobs/1000) %>% 
    mutate(year=years) %>% 
    melt(id.vars='year') %>% 
    mutate(type='1.Desembarques')
  
  Bcru  <-rbind(obsC,obsP,obsD)
  
  p <- ggplot()  +
    geom_bar(data=Bcru, aes(x=year, y =value), 
             stat="identity", fill='gray66', color = 'gray28') + 
    facet_wrap(~type,scale="free",dir = 'v', as.table = TRUE) + 
    labs(x="Años", y="Capturas y biomasas (miles de toneladas)") +  
    theme(panel.background = element_rect(fill ="gray99")) + 
    theme(panel.grid=element_line(color="gray66")) +
    theme(axis.title.y = element_text(margin = margin(r = 20)))
  
  p
  
}

fig2<-function(archivo.Rdata){
  load(archivo.Rdata)
  
  WmF <- as.data.frame(WmedF) %>% 
    mutate(years=years) %>% 
    melt(id.vars='years') %>%             
    mutate(edad = rep(age, each=nyears)) %>% 
    mutate(type='WmedF')
  
  pobsF <- as.data.frame(pfobs) %>% 
    mutate(years=years) %>% 
    melt(id.vars='years') %>%             
    mutate(edad = rep(age, each=nyears)) %>% 
    mutate(type='pobsF')
  
  f1<-ggplot(pobsF, aes(x = years, y = value, fill=edad))+
    geom_bar(aes(fill=variable), colour="grey80", position="stack", stat="identity",size=0.05) +
    scale_fill_hue(l=45) +
    labs(x = '', y = 'Proporción de captura en N° a la edad',fill="Grupos de edad") +
    scale_x_continuous(breaks = seq(from = 1990, to = 2022, by = 5)) +
    scale_fill_discrete(name = "Grupos \nde edad", labels = c('GE 0','GE 1','GE 2','GE 3','GE 4'))+
    ggtitle("FLOTA")+
    theme_bw(base_size=11) + 
    theme(plot.title = element_text(hjust = 0.5),legend.key.size = unit(.4, "cm"))
  
  f2<-ggplot(WmF, aes(x = years, y = value, group=variable,colour=variable))+
    geom_line() +
    geom_point( size=2, shape=21, fill="white") + 
    labs(x = '', y = 'Pesos medios (grs)',fill="",color="Grupos de edad") +
    scale_x_continuous(breaks = seq(from = 1990, to = 2022, by = 5)) +	
    scale_colour_discrete(name = "Grupos \nde edad",  labels = c('GE 0','GE 1','GE 2','GE 3','GE 4'))+
    ggtitle("FLOTA")+
    theme_bw(base_size=11) + 
    theme(plot.title = element_text(hjust = 0.5), legend.key.size = unit(.4, "cm"))
  
  f1 + f2
  
}
#><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
#  Proporción de edad y pesos medios de la flota
#><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>

fig3<-function(archivo.Rdata){
  load(archivo.Rdata)
 
  pF       <- c(pfobs); pF[pF==0]  <-NA
  Wm       <- c(WmedF); Wm[Wm==0]  <-NA     
  
  anos <- rep(years,nage) 
  edad <- gl(nage,nyears,label=age)   
  
  datosProp=data.frame(x=edad,y=anos,tamanio=pF)
  datosWmed=data.frame(x=edad,y=anos,tamanio=Wm )
  
  g1 <- ggplot (datosProp,aes(x,y)) +
    geom_point(aes(size=tamanio),color = 'gray25',
               shape=21, fill="gray85",alpha = 0.7) +
    scale_size_continuous(breaks = seq(0.05,0.65,0.2),range=c(0,6))+
    labs(x = 'Edades', y = 'Años',size="Proporción") +
    ggtitle("Proporción de edad de la Flota")+
    theme_bw(base_size=9) + 
    theme(plot.title = element_text(hjust = 0.5))
  
  g2 <- ggplot (datosWmed,aes(x,y)) + 
    geom_point(aes(size=tamanio),color = 'gray25',
               shape=21, fill="gray85",alpha=0.7) +
    scale_size_continuous(breaks = seq(15,75,20),range=c(0,6))+
    labs(x = 'Edades', y = 'Años',size="Gramos") +
    ggtitle("Pesos medios de la Flota")+
    theme_bw(base_size=9) + 
    theme(plot.title = element_text(hjust = 0.5))
  
  g1 + g2
  
}
#><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
# Proporción de edad de la flota
#><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>

fig4<-function(archivo.Rdata){
  load(archivo.Rdata)

  pobsR <- as.data.frame(pRobs) %>% 
    mutate(years=years) %>% 
    melt(id.vars='years') %>%             
    mutate(edad = rep(age, each=nyears)) %>% 
    mutate(type='pobsR')
  
  pobsP <- as.data.frame(pPobs) %>%
    mutate(years=years) %>% 
    melt(id.vars='years') %>%             
    mutate(edad = rep(age, each=nyears)) %>% 
    mutate(type='pobsP')
  
  f1<-ggplot(pobsR, aes(x = years, y = value, fill=edad))+
    geom_bar(aes(fill=variable), colour="grey80", position="stack", stat="identity",size=0.05) +
    scale_fill_hue(l=45) +
    labs(x = '', y = 'Proporción de captura en N° a la edad',fill="Grupos de edad") +
										   
    scale_x_continuous(breaks = seq(from = 1990, to = 2022, by = 5)) +
    scale_fill_discrete(name = "Grupos \nde edad", labels = c('GE 0','GE 1','GE 2','GE 3','GE 4'))+
    ggtitle("CRUCERO DE VERANO")+
    theme_bw(base_size=11) + 
    theme(plot.title = element_text(hjust = 0.5),legend.position="none")
  
  f2<-ggplot(pobsP, aes(x = years, y = value, fill=edad))+
    geom_bar(aes(fill=variable), colour="grey80", position="stack", stat="identity",size=0.05) +
    scale_fill_hue(l=45) +
    labs(x = '', y = 'Proporción de captura en N° a la edad',fill="Grupos de edad") +
										   
    scale_x_continuous(breaks = seq(from = 1990, to = 2022, by = 5)) +
    scale_fill_discrete(name = "Grupos \nde edad",  labels = c('GE 0','GE 1','GE 2','GE 3','GE 4'))+
    ggtitle("CRUCERO DE OTOÑO")+
												   
																		  
    theme_bw(base_size=11) + 
    theme(plot.title = element_text(hjust = 0.5),legend.key.size = unit(.4, "cm"))
  
  f1 + f2
}
#><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
# Proporción de edad Cruceros
#><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
fig5<-function(archivo.Rdata){
  load(archivo.Rdata)
  
  pR       <- c(pRobs); pR[pR==0]  <-NA
  pP       <- c(pPobs); pP[pP==0]  <-NA 
  
  anos <- rep(years,length(age)) 
  edad <- gl((length(age)),length(years),label=age)   
  
  datosPropR=data.frame(x=edad,y=anos,tamanio=pR)
  datosPropP=data.frame(x=edad,y=anos,tamanio=pP )
  
  g1 <- ggplot (datosPropR,aes(x,y)) +
    geom_point(aes(size=tamanio),color = 'gray25',shape=21, fill="gray85",alpha = 0.7) +
    scale_size_continuous(breaks = seq(0.05,0.65,0.2),range=c(0,6))+
    labs(x = 'Edades', y = 'Años',size="Proporción") +
    ggtitle("Cruceros de Verano")+
    theme_bw(base_size=11) + 
    theme(plot.title = element_text(hjust = 0.5))
  
  g2 <- ggplot (datosPropP,aes(x,y)) + 
    geom_point(aes(size=tamanio),color = 'gray25',shape=21, fill="gray85",alpha=0.7) +
    scale_size_continuous(breaks = seq(0.05,0.65,0.2),range=c(0,6))+
    labs(x = 'Edades', y = 'Años',size="Proporción") +
    ggtitle("Cruceros de otoño")+
    theme_bw(base_size=11) + 
    theme(plot.title = element_text(hjust = 0.5))
  
  g1 + g2
  
}

# DIAGNOSTICO ----
#><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
# 1. Ajuste Índices ----
#><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
fig6<-function(DATA,col_line,type_line,col_fill){
  
  
  cvBcO<-0.3
  cvBcV<-0.3
  cvdes<-0.01
  
BcV <- ggplot(DATA %>% filter(type!='observado', variable=='Crucero_verano'), 
              aes(yrs,value/1000000)) + 
  geom_line(aes(colour=Asesoria,linetype = Asesoria), size=0.8) +
  scale_colour_manual(values=col_line) +
  scale_linetype_manual(values=type_line)+
  geom_point(data = DATA %>% filter(type=='observado', variable=='Crucero_verano'),
             aes(yrs,value/1000000), shape = 19, colour = 'gray30') +
  geom_errorbar(data = DATA %>% filter(type=='observado', variable=='Crucero_verano'),
                aes(ymin = value*exp(-1.96*cvBcO)*10^-6, ymax = value*exp(1.96*cvBcO)*10^-6), color = 'gray30') +
  scale_x_continuous(breaks = seq(from = 1990, to = 2050, by = 4)) +
  labs(x = '', y = 'Biomasas (millones de t)') +
  theme_bw(base_size=11) + 
  ggtitle('Crucero de verano')+
  theme(plot.title = element_text(hjust = 0.5),legend.position="none")

BcP <- ggplot(DATA %>% filter(type!='observado', variable=='Crucero_otoño'), 
              aes(yrs,value/1000000)) + 
  geom_line(aes(colour=Asesoria,linetype = Asesoria), size=0.8) +
  scale_colour_manual(values=col_line,name="Asesoría") +
  scale_linetype_manual(values=type_line,name="Asesoría")+
  geom_point(data = DATA %>% filter(type=='observado', variable=='Crucero_otoño'),
             aes(yrs,value/1000000), shape = 19, colour = 'gray30') +
  geom_errorbar(data = DATA %>% filter(type=='observado', variable=='Crucero_otoño'),
                aes(ymin = value*exp(-1.96*cvBcV)*10^-6, ymax = value*exp(1.96*cvBcV)*10^-6), color = 'gray30') +
  scale_x_continuous(breaks = seq(from = 1990, to = 2050, by = 4)) +
  labs(x = '', y = 'Biomasas (millones de t)') +
  theme_bw(base_size=11) + 
  ggtitle('Crucero de otoño')+
  theme(plot.title = element_text(hjust = 0.5))

BcH <- ggplot(DATA %>% filter(type!='observado', variable=='Crucero_huevos'), 
              aes(yrs,value/1000)) + 
  geom_line(aes(colour=Asesoria,linetype = Asesoria), size=0.8) +
  scale_colour_manual(values=col_line,name="Asesoría") +
  scale_linetype_manual(values=type_line,name="Asesoría")+
  geom_point(data = DATA %>% filter(type=='observado', variable=='Crucero_huevos'),
             aes(yrs,value/1000), shape = 19, colour = 'gray30') +
  scale_x_continuous(breaks = seq(from = 1990, to = 2050, by = 4)) +
  labs(x = '', y = 'Biomasas (miles de t)') +
  theme_bw(base_size=11) + 
  ggtitle('Crucero de huevos')+
  theme(plot.title = element_text(hjust = 0.5),legend.position="none")


d   <- ggplot(DATA %>% filter(type!='observado', variable=='Desembarques'), 
              aes(yrs,value/1000)) +
  geom_line(aes(colour=Asesoria,linetype = Asesoria), size=0.8) +
  scale_colour_manual(values=col_line) +
  scale_linetype_manual(values=type_line)+
  geom_point(data = DATA %>% filter(type=='observado', variable=='Desembarques'),
             aes(yrs,value/1000), shape = 19, colour = 'gray30') +
  geom_errorbar(data = DATA %>% filter(type=='observado', variable=='Desembarques'),
                aes(ymin = value*exp(-1.96*cvdes)*10^-3, ymax = value*exp(1.96*cvdes)*10^-3), color = 'gray30') +
  scale_x_continuous(breaks = seq(from = 1990, to = 2050, by = 5)) +
  labs(x = '', y = 'Capturas (miles de t)') +
  theme_bw(base_size=11) + 
  ggtitle('Desembarques') +
  theme(plot.title = element_text(hjust = 0.5),legend.position="none")

(BcV/BcP|BcH/d) + plot_layout(guides="collect")

}
#><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
# 2. Residuos índices ----
#><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
fig7<-function(archivo.Rdata){
  
  DataRes <- indobs %>%
              mutate(
              Res=(log(indobs$value)-log(indpred$value)),
              Pred=log(indpred$value))
  
r1   <- ggplot(DataRes, aes(yrs,Res)) + 
        geom_bar(stat='identity', 
                 position='dodge') +
  geom_hline(yintercept = 0) + 
  facet_wrap(. ~ variable, ncol = 1) + 
  labs(x= 'Año', y = 'Residuales (escala log)') +
  theme_bw(base_size=9)

r2   <- ggplot(DataRes, aes(Pred,Res)) + 
  geom_point(size = 1.5) + 
  geom_hline(yintercept = 0) + 
  facet_wrap(. ~ variable, ncol = 1) +
  labs(x= 'Predicho (log)', y = 'Residuales') + 
  theme_bw(base_size=8)

r3   <- ggplot(DataRes, aes(Res)) + 
  geom_histogram(fill='white', position  = 'dodge') +
  facet_wrap(. ~ variable, ncol = 1) + 
  labs(x= 'Residuales', y ='Histograma de Residuos (Frecuencia)') +
  theme_bw(base_size=9)

r4   <- ggplot(DataRes, aes(sample = Res)) + 
  stat_qq() + 
  stat_qq_line() + 
  facet_wrap(. ~ variable, ncol = 1) +
  labs(x= 'Sample Quantiles', y ='Theoretical') + 
  theme_bw(base_size=10)

r1.1  <-ggplot(DataRes, aes(yrs,Res)) + 
  geom_smooth()+
  geom_bar(stat='identity', position='dodge') +
  geom_hline(yintercept = 0) + 
  facet_wrap(. ~ variable, ncol = 1) + 
  labs(x= 'Año', y = 'Residuales (escala log)') +
  theme_bw(base_size=9)

r2.2  <-ggplot(DataRes, aes(Pred,Res)) + 
  geom_smooth()+
  geom_point(size = 1.5) + 
  geom_hline(yintercept = 0) + 
  facet_wrap(. ~ variable, ncol = 1) +
  labs(x= 'Predicho (log)', y = 'Residuales') + 
  theme_bw(base_size=8)

r3.3 <- ggplot(DataRes, aes(Res)) + 
  geom_histogram(fill='white', position  = 'dodge') +
  facet_wrap(. ~ variable, ncol = 1) + 
  labs(x= 'Residuales', y ='Histograma de Residuos (Frecuencia)') +
  theme_bw(base_size=9)

r4.4 <- ggplot(DataRes, aes(sample = Res)) + 
  stat_qq() + 
  stat_qq_line() + 
  facet_wrap(. ~ variable, ncol = 1) +
  labs(x= 'Sample Quantiles', y ='Theoretical') + 
  theme_bw(base_size=9)

r1.1 + r2.2 + r4.4

}

#><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
# 3. Ajustes Composición de edad ----
#><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
# Flota ----
fig8 <-function(compEdad,flota,col_line,type_line){
  
  
  
  figf <- ggplot(compEdad %>% filter(type=='observado')) + 
    geom_bar(aes(x = variable, y = value), stat="identity", fill='gray66', color = 'gray28') + 
    facet_wrap(vars(yrs), dir = 'v', as.table = TRUE) + 
    labs(x = 'Edad', y = 'Proporción') +
    geom_line(data = compEdad %>% filter(type=='predicho'),
              aes(x = as.numeric(variable), y = value,colour=Asesoria,linetype =Asesoria)) +
    scale_colour_manual(values=col_line,name="Asesoría") +
    scale_linetype_manual(values=type_line,name="Asesoría")+
    theme(panel.background = element_rect(fill ="gray99")) + 
    theme(panel.grid=element_line(color=NA)) +
    ggtitle(flota, ) + theme(plot.title = element_text(size = 14, face = "bold"))
  figf
  
}


#><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
# 4. Función de figura Residuos composición de edades ----
fig9 <-function(archivo.Rdata, anos,pobs,ppred,title){
  
res  <-pobs-ppred
rng  <-range(res,na.rm=T)
dd   <-dim(res)
est  <-matrix(NA,nrow=dd[1],ncol=dd[2])

for(j in 1:dd[1]){for(k in 1:dd[2]){val<-res[j,k]
if(val>0){est[j,k]<-val/rng[2]}
else{est[j,k]<-val/rng[1]*-1}}}

#par(mfrow=c(1,3),mar=c(5.4,6.7,2,1),cex.axis=1,cex.lab=1.1)
fig1<-image(age,anos,t(est),col=0,yaxt="n",xlab="",ylab="")
ee  <-dim(est)
for(n in 1:ee[1]){for(m in 1:ee[2]){vol<-est[n,m]
if(is.na(vol)==FALSE){
  if(vol>0){points(age[m],anos[n],pch=19,cex=2.82*sqrt(vol),col=1)}
  if(vol<0){points(age[m],anos[n],pch=1,cex=2.82*sqrt(vol*-1),col=1)}
}}}


mtext(title,side=3,cex=1)
															
axis(2,at=anos,labels=anos,las=2)
									  
box()
}
#><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><> ><>
# 5. Comparación con asesorías previas ----

fig10 <-function(Rdata_RetroHist,archivo.Rdata){
  
  
  #archivo.Rdata<-Rdata_H2
  # Rdata_RetroHist<-Rdata_RetroHist
  
  load(archivo.Rdata)
  RTh<-rt
  BDh<-bd
  FTh<-ft
  
  load(Rdata_RetroHist)
  #Retrospectivo tradicional
  library(ggpubr)
  Rt <- ggplot(Rtcomp) + 
    geom_ribbon(data=RTh,aes(ymin=lower, ymax=upper, x=x, fill = "IC"), alpha = 0.2)+
    geom_line(aes(y=Rt_sept21,x=x, colour = year_retros[6]), size=0.5)+
    geom_line(aes(y=Rt_mar22, x=x, colour = year_retros[5]), size=0.5)+
    geom_line(aes(y=Rt_jul22, x=x, colour = year_retros[4]), size=0.5)+
    geom_line(aes(y=Rt_sept22,x=x, colour = year_retros[3]), size=0.5)+
    geom_line(aes(y=Rt_mar23,x=x, colour = year_retros[2]), size=0.5)+
    geom_line(aes(y=Rt_jul23,x=x, colour = year_retros[1]), size=0.5)+
    labs(x = '', y = 'Reclutamientos ',colour='Asesorías')  +
    scale_x_continuous(breaks = seq(from = 1990, to = 2060, by = 4)) +
    scale_colour_manual("",values=c('purple','orange',"green","blue","red","black"))+
    scale_fill_manual("",values=c("grey30"))+
    theme_bw(base_size=9) +
    ggtitle('')+
    theme(plot.title = element_text(hjust = 0.5))
  
  BD <- ggplot(SSBtcomp) + 
    geom_ribbon(data=BDh,aes(ymin=lower, ymax=upper, x=x, fill = "IC"), alpha = 0.2)+
    geom_line(aes(y=SSBt_sept21,x=x, colour = year_retros[6]), size=0.5)+
    geom_line(aes(y=SSBt_mar22, x=x, colour = year_retros[5]), size=0.5)+
    geom_line(aes(y=SSBt_jul22, x=x, colour = year_retros[4]), size=0.5)+
    geom_line(aes(y=SSBt_sept22, x=x, colour = year_retros[3]), size=0.5)+
    geom_line(aes(y=SSBt_mar23, x=x, colour = year_retros[2]), size=0.5)+
    geom_line(aes(y=SSBt_jul23, x=x, colour = year_retros[1]), size=0.5)+
    labs(x = '', y = 'Biomasa desovante (t)',colour='Asesorías')  +
    scale_x_continuous(breaks = seq(from = 1990, to = 2060, by = 4)) +
    scale_colour_manual("",values=c('purple','orange',"green","blue","red","black"))+
    scale_fill_manual("",values=c("grey30"))+
    theme_bw(base_size=9) +
    ggtitle('')+
    theme(plot.title = element_text(hjust = 0.5))
  
  Ft <- ggplot(Ftcomp) + 
    geom_ribbon(data=FTh,aes(ymin=lower, ymax=upper, x=x, fill = "IC"), alpha = 0.2)+
    geom_line(aes(y=Ft_sept21,x=x, colour = year_retros[6]), size=0.5)+
    geom_line(aes(y=Ft_mar22, x=x, colour = year_retros[5]), size=0.5)+
    geom_line(aes(y=Ft_jul22, x=x, colour = year_retros[4]), size=0.5)+
    geom_line(aes(y=Ft_sept22,x=x, colour = year_retros[3]), size=0.5)+
    geom_line(aes(y=Ft_mar23, x=x, colour = year_retros[2]), size=0.5)+
    geom_line(aes(y=Ft_jul23, x=x, colour = year_retros[1]), size=0.5)+
    labs(x = '', y = 'Mortalidad por pesca (1/año)',colour='Asesorías')  +
    scale_x_continuous(breaks = seq(from = 1990, to = 2060, by = 4)) +
    scale_colour_manual("",values=c('purple','orange',"green","blue","red","black"))+
    scale_fill_manual("",values=c("grey30"))+
    theme_bw(base_size=9) +
    ggtitle('')+
    theme(plot.title = element_text(hjust = 0.5))
  
  # Rt/BD/Ft 
  ggarrange(Rt,BD,Ft, ncol=1, nrow=3, common.legend= TRUE, legend= "right")
  
}

# 6. Análisis retrospectivo relativo ----

fig11 <-function(Rdata_retroA,archivo.Rdata){
  
  #archivo.Rdata<-Rdata_H2
  
  load(archivo.Rdata)
  RTa<-rt
  BDa<-bd
  FTa<-ft
  
  load(Rdata_retroA)
  
  
  #Retrospectivo tradicional
  library(ggpubr)
  Rt <- ggplot(Rt_retro) + 
    geom_ribbon(data=RTa,aes(ymin=lower, ymax=upper, x=x, fill = "IC"), alpha = 0.2)+
    geom_line(aes(y=y1, x=x, colour = year_retros[nretros]), size=0.5)+
    geom_line(aes(y=y2, x=x, colour = year_retros[nretros-1]), size=0.5)+
    geom_line(aes(y=y3, x=x, colour = year_retros[nretros-2]), size=0.5)+
    geom_line(aes(y=y4, x=x, colour = year_retros[nretros-3]), size=0.5)+
    geom_line(aes(y=y5, x=x, colour = year_retros[nretros-4]), size=0.5)+
    labs(x = '', y = 'Reclutamiento (miles de millones)',colour='Asesorías')  +
    scale_x_continuous(breaks = seq(from = 1990, to = 2060, by = 5)) +
    scale_colour_manual("",values=c("orange","green","blue","red","black"))+
    scale_fill_manual("",values=c("grey30"))+
    theme_bw(base_size=12) +
    ggtitle('')+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  BD <- ggplot(BD_retro) + 
    geom_ribbon(data=BDa,aes(ymin=lower, ymax=upper, x=x, fill = "IC"), alpha = 0.2)+
    geom_line(aes(y=y1, x=x, colour = year_retros[nretros]), size=0.5)+
    geom_line(aes(y=y2, x=x, colour = year_retros[nretros-1]), size=0.5)+
    geom_line(aes(y=y3, x=x, colour = year_retros[nretros-2]), size=0.5)+
    geom_line(aes(y=y4, x=x, colour = year_retros[nretros-3]), size=0.5)+
    geom_line(aes(y=y5, x=x, colour = year_retros[nretros-4]), size=0.5)+
    labs(x = '', y = 'BD (millones de t)',colour='Asesorías')  +
    scale_x_continuous(breaks = seq(from = 1990, to = 2060, by = 5)) +
    scale_colour_manual("",values=c("orange","green","blue","red","black"))+
    scale_fill_manual("",values=c("grey30"))+
    theme_bw(base_size=12) +
    ggtitle('')+
    theme(plot.title = element_text(hjust = 0.5))
  
  Ft <- ggplot(Ft_retro) + 
    geom_ribbon(data=FTa,aes(ymin=lower, ymax=upper, x=x, fill = "IC"), alpha = 0.2)+
    geom_line(aes(y=y1, x=x, colour = year_retros[nretros]), size=0.5)+
    geom_line(aes(y=y2, x=x, colour = year_retros[nretros-1]), size=0.5)+
    geom_line(aes(y=y3, x=x, colour = year_retros[nretros-2]), size=0.5)+
    geom_line(aes(y=y4, x=x, colour = year_retros[nretros-3]), size=0.5)+
    geom_line(aes(y=y5, x=x, colour = year_retros[nretros-4]), size=0.5)+
    labs(x = '', y = 'F (1/año)',colour='Asesorías')  +
    scale_x_continuous(breaks = seq(from = 1990, to = 2060, by = 5)) +
    scale_colour_manual("",values=c("orange","green","blue","red","black"))+
    scale_fill_manual("",values=c("grey30"))+
    theme_bw(base_size=12) +
    ggtitle('')+
    theme(plot.title = element_text(hjust = 0.5))
  
  #Retrospectivo relativo
  Rtrel <- ggplot(Rt_retroRel) + lims(y=c(-1,1)) +
    geom_line(aes(y=y1, x=x, colour = year_retros[nretros]), size=0.5)+
    geom_line(aes(y=y2, x=x, colour = year_retros[nretros-1]), size=0.5)+
    geom_line(aes(y=y3, x=x, colour = year_retros[nretros-2]), size=0.5)+
    geom_line(aes(y=y4, x=x, colour = year_retros[nretros-3]), size=0.5)+
    geom_line(aes(y=y5, x=x, colour = year_retros[nretros-4]), size=0.5)+
    annotate("text", x=2000, y=0.5,label=paste("Rho =",round(ave.mohn.r,2))) +
    labs(x = '', y = 'Diferencia % último año',colour='Asesorías')  +
    scale_x_continuous(breaks = seq(from = 1990, to = 2060, by = 5)) +
    scale_colour_manual("",values=c("orange","green","blue","red","black"))+
    scale_fill_manual("",values=c("grey30"))+
    theme_bw(base_size=12) +
    ggtitle('')+
    theme(plot.title = element_text(hjust = 0.5))
  
  BDrel <- ggplot(BD_retroRel) + lims(y=c(-1,1)) +
    geom_line(aes(y=y1, x=x, colour = year_retros[nretros]), size=0.5)+
    geom_line(aes(y=y2, x=x, colour = year_retros[nretros-1]), size=0.5)+
    geom_line(aes(y=y3, x=x, colour = year_retros[nretros-2]), size=0.5)+
    geom_line(aes(y=y4, x=x, colour = year_retros[nretros-3]), size=0.5)+
    geom_line(aes(y=y5, x=x, colour = year_retros[nretros-4]), size=0.5)+
    annotate("text", x=2000, y=0.5,label=paste("Rho =",round(ave.mohn.ssb,2))) +
    labs(x = '', y = 'Diferencia % último año',colour='Asesorías')  +
    scale_x_continuous(breaks = seq(from = 1990, to = 2060, by = 5)) +
    scale_colour_manual("",values=c("orange","green","blue","red","black"))+
    scale_fill_manual("",values=c("grey30"))+
    theme_bw(base_size=12) +
    ggtitle('')+
    theme(plot.title = element_text(hjust = 0.5))
  
  Ftrel <- ggplot(Ft_retroRel) + lims(y=c(-1,1)) +
    geom_line(aes(y=y1, x=x, colour = year_retros[nretros]), size=0.5)+
    geom_line(aes(y=y2, x=x, colour = year_retros[nretros-1]), size=0.5)+
    geom_line(aes(y=y3, x=x, colour = year_retros[nretros-2]), size=0.5)+
    geom_line(aes(y=y4, x=x, colour = year_retros[nretros-3]), size=0.5)+
    geom_line(aes(y=y5, x=x, colour = year_retros[nretros-4]), size=0.5)+
    annotate("text", x=2000, y=0.5,label=paste("Rho =",round(ave.mohn.f,2))) +
    labs(x = '', y = 'Diferencia % último año',colour='Asesorías')  +
    scale_x_continuous(breaks = seq(from = 1990, to = 2060, by = 5)) +
    scale_colour_manual("",values=c("orange","green","blue","red","black"))+
    scale_fill_manual("",values=c("grey30"))+
    theme_bw(base_size=12) +
    ggtitle('')+
    theme(plot.title = element_text(hjust = 0.5))
  
  # Rt/BD/Ft |Rtrel/BDrel/Ftrel
  
  ggarrange(Rt,Rtrel,BD,BDrel,Ft,Ftrel, ncol=2, nrow=3, common.legend= TRUE, legend= "top")
  
  
}
# 7. Perfil de verosimilitud ----

fig12 <-function(Rdata_perfil,asesoria,xlimperfil){

  load(Rdata_perfil)
  
  fig<-ggplot() +
    geom_line(aes(x = TLk2$Ro, y = TLk2$Bio_Reclas, colour=names[2])) +
    geom_line(aes(x = TLk2$Ro, y = TLk2$Bio_Pelaces, colour=names[3])) +
    geom_line(aes(x = TLk2$Ro, y = TLk2$Desembarques,colour=names[4])) +
    geom_line(aes(x = TLk2$Ro, y = TLk2$Bio_Mph,colour=names[5])) +
    geom_line(aes(x = TLk2$Ro, y = TLk2$C.Edad_Flota,colour=names[6])) +
    geom_line(aes(x = TLk2$Ro, y = TLk2$C.Edad_Recl,colour=names[7])) +
    geom_line(aes(x = TLk2$Ro, y = TLk2$C.Edad_Pel,colour=names[8])) +
    geom_line(aes(x = TLk2$Ro, y = TLk2$Total,colour=names[17])) +
    geom_hline(yintercept = 2,colour='black',lty=2) +
    scale_colour_manual("",values=c('cyan','magenta','gray',"orange","green","blue","red","black"))+
    coord_cartesian(xlim = xlimperfil, ylim = c(0, 3))+
     xlab(expression("R"[0])) + 
    ylab("L-min(L)") + 
    ggtitle(asesoria)+
    theme_bw(base_size=10) +
    theme(plot.title = element_text(hjust = 0.5),legend.position="right")  
  #fig<-ggplot(data = datares, aes(x = reorder(names_res,residuos), y = residuos)) +
  #      ylim(-2, 3)+
  #      geom_bar(stat = "identity",fill="azure3") +
  #      coord_flip() + # Barras horizontales
  #      xlab("") + 
  #      ylab("Diferencia porcentual de Ro") + 
  #      theme_bw(base_size=10) +
  #     ggtitle("Anchoveta Centro-sur")+
  #      theme(plot.title = element_text(hjust = 0.5))
  
  fig
}

# 7.1. Perfil de verosimilitud ----

fig12.1 <-function(Rdata_perfil,asesoria,xlimperfil){
  
  load(Rdata_perfil)
  
  fig<-ggplot() +
    geom_line(aes(x = TLk2$Ro, y = TLk2$Bio_Reclas, colour=names[2])) +
    geom_line(aes(x = TLk2$Ro, y = TLk2$Bio_Pelaces, colour=names[3])) +
    #geom_line(aes(x = TLk2$Ro, y = TLk2$Desembarques,colour=names[4])) +
    geom_line(aes(x = TLk2$Ro, y = TLk2$Bio_Mph,colour=names[5])) +
    #geom_line(aes(x = TLk2$Ro, y = TLk2$C.Edad_Flota,colour=names[6])) +
    #geom_line(aes(x = TLk2$Ro, y = TLk2$C.Edad_Recl,colour=names[7])) +
    #geom_line(aes(x = TLk2$Ro, y = TLk2$C.Edad_Pel,colour=names[8])) +
    geom_line(aes(x = TLk2$Ro, y = TLk2$Total,colour=names[17])) +
    geom_hline(yintercept = 2,colour='black',lty=2) +
    #scale_colour_manual("",values=c('cyan','magenta','gray',"orange","green","blue","red","black"))+
    scale_colour_manual("",values=c("green","blue","red","black"))+
    coord_cartesian(xlim = xlimperfil, ylim = c(0, 3))+
    xlab(expression("R"[0])) + 
    ylab("L-min(L)") + 
    ggtitle(asesoria)+
    theme_bw(base_size=10) +
    theme(plot.title = element_text(hjust = 0.5),legend.position="right")  
  #fig<-ggplot(data = datares, aes(x = reorder(names_res,residuos), y = residuos)) +
  #      ylim(-2, 3)+
  #      geom_bar(stat = "identity",fill="azure3") +
  #      coord_flip() + # Barras horizontales
  #      xlab("") + 
  #      ylab("Diferencia porcentual de Ro") + 
  #      theme_bw(base_size=10) +
  #     ggtitle("Anchoveta Centro-sur")+
  #      theme(plot.title = element_text(hjust = 0.5))
  
  fig
}

# 8. Variables poblacionales ----

fig13<-function(Rdata_H1,Rdata_H2,Rdata_H3,HITOact,opt.IND,opt.name_IND,col_line,type_line,col_fill,xtext){
  
  load(Rdata_H1); VarPob_H1<-Var
  load(Rdata_H2); VarPob_H2<-Var
  load(Rdata_H3); VarPob_H3<-Var
  
  
  if(HITOact=="Hito 2"){DATA<-rbind(VarPob_H1,VarPob_H2)}
  if(HITOact=="Hito 3"){DATA<-rbind(VarPob_H1,VarPob_H2,VarPob_H3)}
  
  IND   <-c("Rt","BT","BD","Ft")
  name_IND<-c("Reclutamiento (millones)","Biomasa total (t)","Biomasa desovante (t)","Mortalidad por Pesca (1/año)")
  
  ind<-DATA %>% filter(Hito==HITOact,indicador==IND[opt.IND]) 
  meanind<-mean(ind$value)
  name_meanind<-paste(ind$indicador[1],' promedio',sep='')
  if(IND[opt.IND]=="Ft"){name_meanind<-'Ft mediana'}
  
  ggplot()+
    geom_line(data=DATA %>% filter(indicador==IND[opt.IND]),aes(x=x,y=value,colour=Hito,linetype=Hito))+
    geom_hline(yintercept = meanind,colour='black',lty=2) +
    geom_ribbon(data=DATA %>% filter(Hito==HITOact,indicador==IND[opt.IND]),aes(ymin=lower,ymax=upper,x=x),alpha=0.2)+
    geom_text(aes(x=xtext, y=meanind,label=name_meanind,vjust=-0.5)) +
    guides(color = guide_legend(title = "Asesorías"),
           linetype = guide_legend(title = "Asesorías"))+
    labs(x = '', y = name_IND[opt.name_IND])  +
    scale_x_continuous(breaks = seq(from = 1990, to = 2060, by = 5)) +
    scale_colour_manual("",values=col_line)+
    scale_linetype_manual('',values=type_line)+
    scale_fill_manual("",values=col_fill)+
    theme_bw(base_size=12) +
    ggtitle('')+
    theme(plot.title = element_text(hjust = 0.5),legend.position="top")
  
  
}

# 9. Selectividades ----

fig14<-function(archivo.Rdata){
  
  load(archivo.Rdata)
  age     <- seq(0,4,1)                                       
  nage    <- length(age)   
  
  g1 <- ggplot () +
    lims(y=c(0,1))+
    #lineas
    geom_line(aes(x=age,y=sel_Flota))+
    geom_line(aes(x=age,y=sel_CruV))+
    geom_line(aes(x=age,y=sel_CruO),linetype="dashed")+
    #puntos
    geom_point(aes(x=age,y=sel_Flota,shape="Flota"),size=2.5) +
    geom_point(aes(x=age,y=sel_CruV,shape="Cruceros de Verano"),size=2.5) +
    geom_point(aes(x=age,y=sel_CruO,shape="Cruceros de Otoño"),size=2.5) +
    #parámetros
    labs(x = "Edad (años)", y = "Patrón de explotación",shape="Selectividades") +
    ggtitle("")+
    theme_bw(base_size=11) + 
    theme(plot.title = element_text(hjust = 0.5),legend.justification=c(1.1,0), legend.position=c(1,0.1))
  g1
  
}

# Puntos biológicos de referencia en sobre  series BD y F

fig15<-function(Rdata_H1,Rdata_H2,Rdata_H3,HITOact,opt.IND,opt.name_IND,col_line,type_line,col_fill,x_pbrs){
  
  load(Rdata_H1); VarPob_H1<-Var
  load(Rdata_H2); VarPob_H2<-Var
  load(Rdata_H3); VarPob_H3<-Var
  
  if(HITOact=="Hito 2"){DATA<-rbind(VarPob_H1,VarPob_H2)}
  if(HITOact=="Hito 3"){DATA<-rbind(VarPob_H1,VarPob_H2,VarPob_H3)}
  
  IND   <-c("Rt","BT","BD","Ft")
  name_IND<-c("Reclutamiento (millones)","Biomasa total (t)","Biomasa desovante (t)","Mortalidad por Pesca (1/año)")
  
  if(IND[opt.IND]=="BD"){
    lpbrs<-c(BRMS,BLIM,Bo,Bmed)
    ypbrs<-c(BRMS,BLIM,Bo,Bmed)
    col_pbrs<-c('green3','red','blue','black')
    #x_pbrs<-c(rep(2012,3),2005)
    label_pbr<-c(expression("BD"[RMS]),expression("BD"[LIM]),expression("BD"[0]),expression("BD"[promedio]))
  }
  
  if(IND[opt.IND]=="Ft"){
    lpbrs<-c(FRMS,Fmedian)
    ypbrs<-c(0.15,0.35)
    col_pbrs<-c('green3','black')
    #x_pbrs<-c(2011,2001)
    label_pbr<-c(expression("F"[RMS]),expression("F"[mh]))
  }
  
  ggplot()+
    geom_line(data=DATA %>% filter(indicador==IND[opt.IND]),aes(x=x,y=value,colour=Hito,linetype=Hito))+
    geom_hline(yintercept = lpbrs,colour=col_pbrs)+
    annotate("text", x=x_pbrs, y=ypbrs,label=label_pbr,vjust=-0.4, colour=col_pbrs) +
    geom_ribbon(data=DATA %>% filter(Hito==HITOact,indicador==IND[opt.IND]),
                aes(ymin=lower,ymax=upper,x=x),alpha=0.2)+
    guides(color = guide_legend(title = "Asesorías"),
           linetype = guide_legend(title = "Asesorías"))+
    labs(x = '', y = name_IND[opt.name_IND])  +
    scale_x_continuous(breaks = seq(from = 1960, to = 2060, by = 5)) +
    scale_colour_manual("",values=col_line)+
    scale_linetype_manual('',values=type_line)+
    scale_fill_manual("",values=col_fill)+
    theme_bw(base_size=12) +
    ggtitle('')+
    theme(plot.title = element_text(hjust = 0.5),legend.position="top") 
}


# 10. Puntos biológicos de referencia ----

fig16<-function(archivo.Rdata,Fspr,BDspr){
  
  load(archivo.Rdata)
  
  age<-seq(0,4,1)
  madurez   <- dat$madurezsexual
  
  g1 <- ggplot () +
    lims(y=c(0,1))+
    geom_line(aes(x=age,y=sel_Flota))+
    geom_line(aes(x=age,y=madurez),linetype="dashed")+
    geom_point(aes(x=age,y=sel_Flota,shape="Selectividad de la flota"),size=2.5) +
    geom_point(aes(x=age,y=madurez,shape="Madurez sexual"),size=2.5) +
    labs(x = 'Edad (años)', y = 'Madurez y selectividad',shape="") +
    ggtitle("")+
    theme_bw(base_size=11) + 
    theme(plot.title = element_text(hjust = 0.5),legend.justification=c(1.1,0), legend.position=c(1,0.1))
  
  
  g2 <- ggplot () +
    geom_line(aes(x=Fspr,y=BDspr))+
    geom_hline(yintercept = 0.6,colour=c('gray35'),linetype="dashed") +
    geom_vline(xintercept = FRMS,colour=c('gray35'),linetype="dashed") +
    annotate("text", x=2, y=0.6,label=c(expression("F"[RMS])),vjust=-0.5) +
    labs(x = 'Mortalidad por pesca (F)', y = '%BDPR',shape="") +
    ggtitle("")+
    theme_bw(base_size=11) + 
    theme(plot.title = element_text(hjust = 0.5),legend.justification=c(1.1,0), legend.position=c(1,0.1))
  g1 + g2
  
}


# 11. Indicadores de estatus relativos a BDrms y Frms----

fig17<-function(Rdata_H1,Rdata_H2,Rdata_H3,HITOact,opt.IND,opt.name_IND,col_line,type_line,col_fill,x_pbrs){
  
  load(Rdata_H1); IndStatus_H1<-IndStatus
  load(Rdata_H2); IndStatus_H2<-IndStatus
  load(Rdata_H3); IndStatus_H3<-IndStatus
  
  if(HITOact=="Hito 2"){DATA<-rbind(IndStatus_H1,IndStatus_H2)}
  if(HITOact=="Hito 3"){DATA<-rbind(IndStatus_H1,IndStatus_H2,IndStatus_H3)}
  
  IND   <-c("F_FRMS","BD_BDRMS")
  name_IND<-c("F/F","BD/BD")
  
  if(IND[opt.IND]=="BD_BDRMS"){
    lpbrs<-c(BRMS/BRMS,BLIM/BRMS)
    col_pbrs<-c('green3','red')
    label_pbr<-c(expression("BD"[RMS]),expression("BD"[LIM]))
  }
  
  if(IND[opt.IND]=="F_FRMS"){
    lpbrs<-c(FRMS/FRMS)
    col_pbrs<-c('green3')
    label_pbr<-c(expression("F"[RMS]))
  }
  
  
  ggplot()+
    geom_line(data=DATA %>% filter(indicador==IND[opt.IND]),aes(x=x,y=value,colour=Hito,linetype=Hito))+
    geom_hline(yintercept = lpbrs,colour=col_pbrs)+
    annotate("text", x=x_pbrs, y=lpbrs,label=label_pbr,vjust=-0.5) +
    geom_ribbon(data=DATA %>% filter(Hito==HITOact,indicador==IND[opt.IND]),
                aes(ymin=lower,ymax=upper,x=x),alpha=0.2)+
    guides(color = guide_legend(title = "Asesorías"),
           linetype = guide_legend(title = "Asesorías"))+
    {if(IND[opt.IND]=="BD_BDRMS")labs(x = '', y = expression("BD/BD"[RMS]))}  +
    {if(IND[opt.IND]=="F_FRMS")labs(x = '', y = expression("F/F"[RMS]))}  +
    scale_x_continuous(breaks = seq(from = 1960, to = 2060, by = 4)) +
    scale_colour_manual("",values=col_line)+
    scale_linetype_manual('',values=type_line)+
    scale_fill_manual("",values=col_fill)+
    theme_bw(base_size=11) +
    ggtitle('')+
    theme(plot.title = element_text(hjust = 0.5),legend.position="top")
  
}



fig18<-function(Rdata_H1,Rdata_H2,Rdata_H3,HITOact,opt.IND,
                col_line,type_line,col_fill,ylim,xdens,ydens){
  
  IND   <-c("F_FRMS","BD_BDRMS")
  
  load(Rdata_H1); densb_H1<-densb; densbf_H1<-densbf;icb_H1<-icb;icf_H1<-icf
  IC_densb_h1<-paste("IC95%_",hito," = [",round(icb_H1[1],3),"-",round(icb_H1[2],3),"]",sep=" ")
  IC_densf_h1<-paste("IC95%_",hito,"= [",round(icf_H1[1],3),"-",round(icf_H1[2],3),"]",sep=" ")
  
  load(Rdata_H2); densb_H2<-densb; densbf_H2<-densbf;icb_H2<-icb;icf_H2<-icf
  IC_densb_h2<-paste("IC95%_",hito," = [",round(icb_H2[1],3),"-",round(icb_H2[2],3),"]",sep=" ")
  IC_densf_h2<-paste("IC95%_",hito,"= [",round(icf_H2[1],3),"-",round(icf_H2[2],3),"]",sep=" ")
  
  load(Rdata_H3); densb_H3<-densb; densbf_H3<-densbf;icb_H3<-icb;icf_H3<-icf
  IC_densb_h3<-paste("IC95%_",hito," = [",round(icb_H3[1],3),"-",round(icb_H3[2],3),"]",sep=" ")
  IC_densf_h3<-paste("IC95%_",hito,"= [",round(icf_H3[1],3),"-",round(icf_H3[2],3),"]",sep=" ")
  
  if(HITOact=="Hito 2"){DATA<-rbind(densb_H1,densb_H2)%>% filter(indicador==IND[opt.IND])}
  if(HITOact=="Hito 3"){DATA<-rbind(densb_H1,densb_H2,densb_H3)%>% filter(indicador==IND[opt.IND])} 
  
  if(HITOact=="Hito 2"){DATA2<-rbind(densbf_H1,densbf_H2)%>% filter(indicador==IND[opt.IND])}
  if(HITOact=="Hito 3"){DATA2<-rbind(densbf_H1,densbf_H2,densbf_H3)%>% filter(indicador==IND[opt.IND])}
  #
  
  if(IND[opt.IND]=="BD_BDRMS"){labeldens<-c(IC_densb_h1,IC_densb_h2,IC_densb_h3)}
  if(IND[opt.IND]=="F_FRMS"){labeldens<-c(IC_densf_h1,IC_densf_h2,IC_densf_h3)}
  
  # plot
  ggplot() + 
    lims(y=ylim) +
    geom_polygon(data=DATA,aes(x=x, y=y,fill=Hito),alpha=0.2) +
    geom_line(data=DATA2,aes(x,y,colour=Hito,linetype=Hito))+
    annotate("text",
             xdens,
             ydens,
             colour = col_line,
             size = 2.5,
             label=labeldens) +
    {if(IND[opt.IND]=="F_FRMS")labs(x = expression("F"[last]*"/F"[RMS]),y = 'Densidad de probabilidad')} +
    {if(IND[opt.IND]=="BD_BDRMS")labs(x = expression("BD"[last]*"/BD"[RMS]),y = 'Densidad de probabilidad')} +
    scale_fill_manual("",values=col_fill)+
    scale_colour_manual("",values=col_line)+
    scale_linetype_manual('',values=type_line)+
    theme_bw(base_size=11) + 
    theme(plot.title = element_text(hjust = 0.5),legend.position="top")
  
}


# 12. Estatus ><> ><> ><> ><> ----

# diagrama fase  ----
fig19<-function(archivo.Rdata){
  
  source(here(dir.fun,"Fn_DiagramaFase2.R"))
  
  load(archivo.Rdata)
  
  name<-paste("Asesoría",hito,"-",meshito,years[nyears],sep=" ")
  preliminar=F
  completo=T
  if(hito=="Hito 2"){preliminar=T;completo=F}
  if(hito=="Hito 3"){preliminar=T;completo=F}
  
  yearpoints<-c(yearsbiol[1],yearsbiol[nyears],yearsbiol[nyears-1])
  
  DiagramaFase2(name,
                years[1:nyears-1],
                SpB[1:nyears-1],
                SpBSE[1:nyears-1],
                ln_Fyr[1:nyears-1],
                ln_FSE[1:nyears-1],
                SpB[nyears],
                SpBSE[nyears],
                ln_Fyr[nyears],
                ln_FSE[nyears],
                FRMS,
                BRMS,
                BLIM,
                FLIM,
                color=F,
                dir.1,
                etiqueta=F,
                preliminar,
                completo)
  
  text(c(SpB[1]/BRMS,
         SpB[nyears]/BRMS,
         SpB[nyears-1]/BRMS),
       c(exp(ln_Fyr[1])/FRMS-0.05,
         exp(ln_Fyr[nyears])/FRMS-0.05,
         exp(ln_Fyr[nyears-1])/FRMS+0.05),
       yearpoints,cex=0.9)
  text(1.8,2.5,name,cex=1,col=1)
}


# 13. Proyección ><> ><> ><> ><>  ----

fig20<-function(Rdata_proy_Hito,escRecl,col_escRecl){
  
  load(Rdata_proy_Hito)
  
   ggplot(data=DataProy,aes(y=value, x=years, colour = variable)) + 
   annotate("rect", 
            xmin = yearProy[1]-1, xmax = yearProy[2]+1,
            ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "pink") + 
   geom_vline(xintercept = yearProy[1]-1,linetype ="dashed") +
   geom_line()+
   geom_point()+
   facet_wrap(vars(fct_recode(indicador,
                              "Reclutamiento (millones)"="Rt",
                              "Biomasa desovante (miles de t)"="BD",
                              "BD/BDrms" ="BD_BDrms",
                              "Capturas (miles de t)"="Ct")%>% 
                     # Change factor level order
                     fct_relevel("Reclutamientos")),
              dir = 'h', as.table = TRUE,scales='free_y',ncol=2) +  
   scale_x_continuous(breaks = seq(from = 1960, to = 2062, by = 4)) +
   scale_colour_manual(values=col_escRecl,
                       labels=escRecl,
                       name='Escenarios de\n Reclutamientos')+
   labs(x = '', y = '') +
   theme_bw(base_size=12) + 
   theme(plot.title = element_text(hjust = 0.5),legend.position="top")
 
}




# diagrama fase proyectado ----

fig21<-function(archivo.Rdata,dir.0,carpetaCBA,escRecl,op_escRecl){
  
  source(here("funciones","Fn_DiagramaFase_proy.R"))
  
  load(archivo.Rdata)
  
  escR<-c(rep(escRecl[1],1), 
          rep(escRecl[2],1),
          rep(escRecl[3],1))
  
  pF<-rep(c("Frms*1",
            "Frms*0,9",
            "Frms*1,1"),3)
  
  nesc_escR<-c(11,12,13)
  nesc_pF<-c(11,12,13)
  
  yearpoints<-c(yearsbiol[1],yearsbiol[nyears],yearsbiol[nyears-1])
  

  name<-escR[op_escRecl]
  DiagramaFase_proy(name,
                years[(nyears-1):nyears],
                SpB[(nyears-1):nyears],
                SpBSE[(nyears-1):nyears],
                ln_Fyr[(nyears-1):nyears],
                ln_FSE[(nyears-1):nyears],
                FRMS,
                BRMS,
                BLIM,
                FLIM,
                color=F,
                dir.0,
                etiqueta=F)

  setwd(paste(dir.0,carpetaCBA,sep=""))
  rep     <- reptoRlist(paste(admb_base,nesc_escR[op_escRecl],".rep",sep=""))
  std     <- read.table(paste(admb_base,nesc_escR[op_escRecl],".std",sep=""),header=T,sep="",na="NA",fill=T) 
  SpBp    <- subset(std,name=="BD_p0")$value                              
  SpBSEp  <- subset(std,name=="BD_p0")$std    
  Frmsp    <-c(rep$Fref,rep$Fref)
  
  
  lines(c(rprlast,SpBp/BRMS),
        c(Frprlast,Frmsp/FRMS))
  
  
  points(c(SpBp/BRMS),
         c(Frmsp/FRMS),
         pch=19,col=c(2,2),cex=c(1,1))
  
  text((SpBp/BRMS),
       c(Frmsp/FRMS-0.07),
       yearbiol_proy,cex=c(1,1)) # diferencia a fig21 hito 2
  
  
  text(c(SpB[1]/BRMS,
         SpB[nyears]/BRMS,
         SpB[nyears-1]/BRMS),
       c(exp(ln_Fyr[1])/FRMS-0.05,
         exp(ln_Fyr[nyears])/FRMS-0.05,
         exp(ln_Fyr[nyears-1])/FRMS+0.05),
       yearpoints,cex=1)
  text(1.5,1.9,name,cex=1,col=1)
  
}

# Figura de proyección a la edad ----
fig22_sard<-function(Rdata_proy_Hito,escRecl,col_escRecl,ylim,yearProy){

  yearbiol_proy  <-c(paste(yearProy[1]-1,"/",str_sub(yearProy[1],3,4),sep=""),
                     paste(yearProy[1],"/",str_sub(yearProy[2],3,4),sep=""))
  load(Rdata_proy_Hito)
  
  
  
  par(mfrow=c(1,2),mar=c(4,4,1,1)+0.5)
  #===============================================================================================================
  # plot captura a la talla
  YTPalto     <-reps2a$YTP_p0W_proyectada
  YTPmedio    <-reps3a$YTP_p0W_proyectada
  YTPbajo     <-reps1a$YTP_p0W_proyectada
  
  C1alto     <- C1eryearR2
  C1medio    <- C1eryearR3
  C1bajo     <- C1eryearR1
  
  C1alto2     <- C1eryearR22
  C1medio2    <- C1eryearR32
  C1bajo2     <- C1eryearR12
  
  C1alto2a     <- C1eryearR22a
  C1medio2a    <- C1eryearR32a
  C1bajo2a     <- C1eryearR12a
  
  plot(seq(0,4,1),YTPalto[1,],
       type="h",main=paste("Año",yearbiol_proy[1] ,sep=" "),lwd=20,ylim=c(0,ylim),yaxs= "i",
       ylab="Captura  a la edad (t)",xlab="Edades",
       cex.axis=1,cex.main=1,cex.lab=1,col="red")
   lines(seq(0,4,1),YTPmedio[1,],type="h",col="green",lwd=15)
   lines(seq(0,4,1),YTPbajo[1,],type="h",col="blue",lwd=10)
   
   text(rep(0.35,3),c(YTPalto[1,1],
                      YTPmedio[1,1],
                      YTPbajo[1,1]),
              paste(c(C1alto,
                      C1medio,
                      C1bajo)*100,"%",sep=""))
   
   text(4,ylim-10000,'a)')

   
  plot(seq(0,4,1),YTPalto[2,],type="h",main=paste("Año",yearbiol_proy[2],sep=" "),
       lwd=20,ylim=c(0,ylim),yaxs= "i", ylab="Captura a la edad (t)",
       xlab="Edades",cex.axis=1,cex.main=1,cex.lab=1,col="red")
  lines(seq(0,4,1),YTPmedio[2,],type="h",col="green",lwd=15)
  lines(seq(0,4,1),YTPbajo[2,],type="h",col='blue2',lwd=10)
  
  
  text(rep(0.35,3),c(YTPalto[2,1],
                     YTPmedio[2,1],
                     YTPbajo[2,1]),
       paste(c(C1alto2,C1medio2,C1bajo2)*100,"%",sep=""))
  
  text(rep(1.4,3),c(YTPalto[2,2],
                    YTPmedio[2,2],
                    YTPbajo[2,2]),
       paste(c(C1alto2a,C1medio2a,C1bajo2a)*100,"%",sep=""))
  
  text(4,ylim-10000,'b)')
  
}

# Figura de densidad de probabilidad CBA ----

fig23<-function(Rdata_proy_Hito,escRecl,col_escRecl,ylim1,ylim2,xlim1,
                leg1,leg2,text1,text2,yearProy){
  
  escR<-c(rep(escRecl[1],1), 
          rep(escRecl[2],1),
          rep(escRecl[3],1))
  
  yearbiol_proy  <-c(paste(yearProy[1]-1,"/",str_sub(yearProy[1],3,4),sep=""),
                     paste(yearProy[1],"/",str_sub(yearProy[2],3,4),sep=""))
  load(Rdata_proy_Hito)
  
  
  par(mfrow=c(1,2),mar=c(4,4,1,1)+0.5)
  #===============================================================================
  # plot densidad de probabilidad CBA
  plot(xca,yca,type="l",ylab="Densidad de probabilidad",
       xaxs="i",yaxs= "i",xlab=paste("CBA", yearProy[1] ,"(toneladas)",sep=" "),main="",
       cex.axis=1,cex.main=1,cex.lab=1,xlim=xlim1,ylim=ylim1)
  
  polygon(xxca,yyca,col=gray(0.8,0.7),border="gray80")
  lines(xca,yca,lwd=1,col="blue2",lty=1)
  
  polygon(xxcb,yycb,col=gray(0.8,0.7),border="gray80")
  lines(xcb,ycb,lwd=1,col="red2",lty=1)
  
  polygon(xxcc,yycc,col=gray(0.8,0.7),border="gray80")
  lines(xcc,ycc,lwd=1,col="green2",lty=1)
  
  legend(leg1[1],leg1[2],
         c(paste("CBA",yearProy[1],escR[1],"= [",round(icca[1],0),"-",round(icca[2],0),"]",sep=" "),
           paste("CBA",yearProy[1],escR[2],"= [",round(iccb[1],0),"-",round(iccb[2],0),"]",sep=" "),
           paste("CBA",yearProy[1],escR[3],"= [",round(iccc[1],0),"-",round(iccc[2],0),"]",sep=" ")),
          lty=c(1,1,1),col=c('blue2',"red2","green2"),bty="n",lwd=1,cex=0.8)
  text(text1[1],text1[2],'a)')
  box()
  
  #===============================================================================
  # plot percentiles de captura
  plot(seq(10,50,10),CBA_sept[1,],type="o", 
       ylab=paste("CBA", yearProy[1] ,"(toneladas)",sep=" "),xlab="Percentil de Captura (P*)",
       main="",cex.axis=1,cex.main=1,cex.lab=1,col="blue2",ylim=ylim2)
  lines(seq(10,50,10),CBA_sept[2,],type="o",col="red2")
  lines(seq(10,50,10),CBA_sept[3,],type="o",col="green2")
  
  legend(leg2[1],leg2[2],escR,
         col=c('blue2',"red2","green2"),pch=c(19,19,19),lwd=1,bty="n",cex=0.8)
  text(text2[1],text2[2],'b)')
  
  
}



fig21_H2y3<-function(archivo.Rdata,dir.0,carpetaCBA,escRecl,op_escRecl){
  
  
  source(here("funciones","Fn_DiagramaFase_proy.R"))
  
  load(archivo.Rdata)
  
  escR<-c(rep(escRecl[1],1), 
          rep(escRecl[2],1),
          rep(escRecl[3],1))
  
  pF<-rep(c("Frms*1",
            "Frms*0,9",
            "Frms*1,1"),3)
  
  nesc_escR<-c(11,12,13)
  nesc_pF<-c(11,12,13)
  
  yearpoints<-c(yearsbiol[1],yearsbiol[nyears],yearsbiol[nyears-1])
  
  
  name<-escR[op_escRecl]
  DiagramaFase_proy(name,
                    years[(nyears-1):nyears],
                    SpB[(nyears-1):nyears],
                    SpBSE[(nyears-1):nyears],
                    ln_Fyr[(nyears-1):nyears],
                    ln_FSE[(nyears-1):nyears],
                    FRMS,
                    BRMS,
                    BLIM,
                    FLIM,
                    color=F,
                    dir.0,
                    etiqueta=F)
  
  setwd(paste(dir.0,carpetaCBA,sep=""))
  rep     <- reptoRlist(paste(admb_base,nesc_escR[op_escRecl],".rep",sep=""))
  std     <- read.table(paste(admb_base,nesc_escR[op_escRecl],".std",sep=""),header=T,sep="",na="NA",fill=T) 
  SpBp    <- subset(std,name=="BD_p0")$value[1]                              
  SpBSEp  <- subset(std,name=="BD_p0")$std[1]    
  Frmsp    <-c(rep$Fref) # DIFERENTE EN ANCHOVETA 
  
  
  lines(c(rprlast,SpBp/BRMS),
        c(Frprlast,Frmsp/FRMS))
  
  
  points(c(SpBp/BRMS),
         c(Frmsp/FRMS),
         pch=19,col=c(2,2),cex=c(1,1))
  
  text((SpBp/BRMS),
       c(Frmsp/FRMS+0.07),
       yearbiol_proy[1],cex=c(1,1))
  
  
  text(c(SpB[1]/BRMS,
         SpB[nyears]/BRMS,
         SpB[nyears-1]/BRMS),
       c(exp(ln_Fyr[1])/FRMS-0.05,
         exp(ln_Fyr[nyears])/FRMS-0.05,
         exp(ln_Fyr[nyears-1])/FRMS+0.05),
       yearpoints,cex=1)
  text(1.5,1.9,name,cex=1,col=1)
  
}

# Figura de proyección a la edad ----
fig22_H2y3<-function(Rdata_proy_Hito,escRecl,col_escRecl,ylim,yearProy){
  
  yearbiol_proy  <-c(paste(yearProy[1]-1,"/",str_sub(yearProy[1],3,4),sep=""),
                     paste(yearProy[1],"/",str_sub(yearProy[2],3,4),sep=""))
  load(Rdata_proy_Hito)
  
  
  
  par(mfrow=c(1,2),mar=c(4,4,1,1)+0.5)
  #===============================================================================================================
  # plot captura a la talla
  YTPalto     <-reps2a$YTP_p0W_proyectada
  YTPmedio    <-reps3a$YTP_p0W_proyectada
  YTPbajo     <-reps1a$YTP_p0W_proyectada
  
  C1alto     <- C1eryearR2
  C1medio    <- C1eryearR3
  C1bajo     <- C1eryearR1
  
  # C1alto2     <- C1eryearR12
  # C1medio2    <- C1eryearR22
  # C1bajo2     <- C1eryearR32
  # 
  # C1alto2a     <- C1eryearR12a
  # C1medio2a    <- C1eryearR22a
  # C1bajo2a     <- C1eryearR32a
  # 
  yearbiol_actual  <-c(paste(years[nyears-1],"/",str_sub(years[nyears],3,4),sep=""))
  
  plot(seq(0,4,1),reps1a$YTP_r0W_actual,
       type="h",main=paste("Año",yearbiol_actual ,sep=" "),lwd=20,ylim=c(0,ylim),yaxs= "i",
       ylab="Captura  a la edad (t)",xlab="Edades",
       cex.axis=1,cex.main=1,cex.lab=1,col='gray30')
  
  
  
  text(c(rep(0.35,1),1.35),
       c(reps1a$YTP_r0W_actual[1],
         reps1a$YTP_r0W_actual[2]),
       paste(c(C1eryearR1act,C1eryearR1act2)*100,"%",sep=""))
  
  text(4,ylim-10000,'a)')
  
  
  plot(seq(0,4,1),YTPalto[1,],
       type="h",main=paste("Año",yearbiol_proy[1] ,sep=" "),lwd=20,ylim=c(0,ylim),yaxs= "i",
       ylab="Captura  a la edad (t)",xlab="Edades",
       cex.axis=1,cex.main=1,cex.lab=1,col='red2')
  lines(seq(0,4,1),YTPmedio[1,],type="h",col='green2',lwd=15)
  lines(seq(0,4,1),YTPbajo[1,],type="h",col='blue2',lwd=10)
  
  text(rep(0.35,3),c(YTPalto[1,1],
                     YTPmedio[1,1],
                     YTPbajo[1,1]),
       paste(c(C1alto,
               C1medio,
               C1bajo)*100,"%",sep=""))
  
  text(4,ylim-10000,'b)')
  
}

fig23_HITO2y3<-function(Rdata_proy_Hito,escRecl,col_escRecl,ylim1,ylim2,xlim1,
                        leg1,leg2,text1,text2,yearProy){
  
  escR<-c(rep(escRecl[1],1), 
          rep(escRecl[2],1),
          rep(escRecl[3],1))
  
  yearbiol_proy  <-c(paste(yearProy[1]-1,"/",str_sub(yearProy[1],3,4),sep=""))
  load(Rdata_proy_Hito)
  
  
  par(mfrow=c(1,2),mar=c(4,4,1,1)+0.5)
  #===============================================================================
  # plot densidad de probabilidad CBA
  plot(xca,yca,type="l",ylab="Densidad de probabilidad",
       xaxs="i",yaxs= "i",xlab=paste("CBA", yearProy[1]-1 ,"(toneladas)",sep=" "),main="",
       cex.axis=1,cex.main=1,cex.lab=1,xlim=xlim1,ylim=ylim1)
  
  polygon(xxca,yyca,col=gray(0.8,0.7),border="gray80")
  lines(xca,yca,lwd=1,col=col_escRecl[1],lty=1)
  
  polygon(xxcb,yycb,col=gray(0.8,0.7),border="gray80")
  lines(xcb,ycb,lwd=1,col=col_escRecl[2],lty=1)
  
  polygon(xxcc,yycc,col=gray(0.8,0.7),border="gray80")
  lines(xcc,ycc,lwd=1,col=col_escRecl[3],lty=1)
  
  legend(leg1[1],leg1[2],
         c(paste("CBA",yearProy[1]-1,escR[1],"= [",round(icca[1],0),"-",round(icca[2],0),"]",sep=" "),
           paste("CBA",yearProy[1]-1,escR[2],"= [",round(iccb[1],0),"-",round(iccb[2],0),"]",sep=" "),
           paste("CBA",yearProy[1]-1,escR[3],"= [",round(iccc[1],0),"-",round(iccc[2],0),"]",sep=" ")),
         lty=c(1,1,1),col=col_escRecl,bty="n",lwd=1,cex=0.8)
  text(text1[1],text1[2],'a)')
  box()
  
  #===============================================================================
  # plot percentiles de captura
  plot(seq(10,50,10),CBA_marzo[1,],type="o", 
       ylab=paste("CBA", yearProy[1]-1 ,"(toneladas)",sep=" "),xlab="Percentil de Captura (P*)",
       main="",cex.axis=1,cex.main=1,cex.lab=1,col=col_escRecl[1],ylim=ylim2)
  lines(seq(10,50,10),CBA_marzo[2,],type="o",col=col_escRecl[2])
  lines(seq(10,50,10),CBA_marzo[3,],type="o",col=col_escRecl[3])
  
  legend(leg2[1],leg2[2],escR,
         col=col_escRecl,pch=c(19,19,19),lwd=1,bty="n",cex=0.8)
  text(text2[1],text2[2],'b)')
  
  
}

#comparacion de escenarios de captura
fig24_HITO2y3<-function(escRecl,col_escRecl,ylim1,xlim1,
                        leg1,text1,yearProy){
  
  escR<-c(rep(escRecl[1],1), 
          rep(escRecl[2],1),
          rep(escRecl[3],1))
  
  yearbiol_proy  <-c(paste(yearProy[1]-1,"/",str_sub(yearProy[1],3,4),sep=""))
  load(Rdata_proy_Hito2)
  xccS0<-xcc
  yccS0<-ycc
  xxccS0<-xxcc
  yyccS0<-yycc
  icccS0<-iccc
  
  load(Rdata_proy_Hito2S1)
  xccS1<-xcc
  yccS1<-ycc
  xxccS1<-xxcc
  yyccS1<-yycc
  icccS1<-iccc
  
  load(Rdata_proy_Hito2S2)
  xccS2<-xcc
  yccS2<-ycc
  xxccS2<-xxcc
  yyccS2<-yycc
  icccS2<-iccc
  
  load(Rdata_proy_Hito2S3)
  xccS3<-xcc
  yccS3<-ycc
  xxccS3<-xxcc
  yyccS3<-yycc
  icccS3<-iccc
  
  load(Rdata_proy_Hito2S4)
  xccS4<-xcc
  yccS4<-ycc
  xxccS4<-xxcc
  yyccS4<-yycc
  icccS4<-iccc
  
  
  par(mfrow=c(1,1),mar=c(4,4,1,1)+0.5)
  #===============================================================================
  # plot densidad de probabilidad CBA
  plot(xccS0,yccS0,type="l",ylab="Densidad de probabilidad",
       xaxs="i",yaxs= "i",xlab=paste("CBA", yearProy[1]-1 ,"(toneladas)",sep=" "),main="",
       cex.axis=0.8,cex.main=0.8,cex.lab=0.8,xlim=xlim1,ylim=ylim1)
  
  polygon(xxccS0,yyccS0,col=gray(0.8,0.7),border="gray80")
  lines(xccS0,yccS0,lwd=1,col=1,lty=1)
  
  polygon(xxccS1,yyccS1,col=gray(0.8,0.7),border="gray80")
  lines(xccS1,yccS1,lwd=1,col=2,lty=1)
  
  polygon(xxccS2,yyccS2,col=gray(0.8,0.7),border="gray80")
  lines(xccS2,yccS2,lwd=1,col=3,lty=1)
  
  polygon(xxccS3,yyccS3,col=gray(0.8,0.7),border="gray80")
  lines(xccS3,yccS3,lwd=1,col=4,lty=1)
  
  polygon(xxccS4,yyccS4,col=gray(0.8,0.7),border="gray80")
  lines(xccS4,yccS4,lwd=1,col=5,lty=1)
  
  
  legend(leg1[1],leg1[2],
         c(paste("CBA",yearProy[1]-1,"Caso base","= [",round(icccS0[1],0),"-",round(icccS0[2],0),"]",sep=" "),
           paste("CBA",yearProy[1]-1,"S1","= [",round(icccS1[1],0),"-",round(icccS1[2],0),"]",sep=" "),
           paste("CBA",yearProy[1]-1,"S2","= [",round(icccS2[1],0),"-",round(icccS2[2],0),"]",sep=" "),
           paste("CBA",yearProy[1]-1,"S3","= [",round(icccS3[1],0),"-",round(icccS3[2],0),"]",sep=" "),
           paste("CBA",yearProy[1]-1,"S4","= [",round(icccS4[1],0),"-",round(icccS4[2],0),"]",sep=" ")),
         lty=c(1,1,1),col=c(1,2,3,4,5),bty="n",lwd=1,cex=0.7)
  text(text1[1],text1[2],'a)')
  box()
  
 
  
  
}

