
CBA<-function(dir.0,dir.1,Carpeta,admb,l_opt_proy,l_opRec,l_opt_wmed,l_opt_Str,l_mf,opt_proy,opt_wmed,opt_Str,system){
  
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
  
  opRec<-seq(4,6,1) #//ESCENARIO DE RECLUTAMIENTO PROMEDIO
  mf<-c(1,0.9,1.1)
  
  data_file<-(paste(admb,".dat",sep="")) #// archivos de datos originales
  S<-readLines(data_file,encoding="UTF-8")
  Se<-S
  
  for( i in 1:3){
    for(j in 1:3){
      Se[l_opt_proy]<-opt_proy #opt_proy 
      Se[l_opRec]<-opRec[i] #opRec
      Se[l_opt_wmed]<-opt_wmed #opWed 
      Se[l_opt_Str]<-opt_Str
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
