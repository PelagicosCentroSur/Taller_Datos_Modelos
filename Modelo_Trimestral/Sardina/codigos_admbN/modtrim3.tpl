GLOBALS_SECTION
 #include <admodel.h>
 #include <stdio.h>
 #include <time.h>
 time_t start,finish;
 long hour,minute,second;
 double elapsed_time;
 ofstream mcmc_report("mcmc.csv");

TOP_OF_MAIN_SECTION
 time(&start);
 arrmblsize = 90000000; 
 gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7); 
 gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7); 
 gradient_structure::set_MAX_NVAR_OFFSET(5000); 
 gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000); 


DATA_SECTION
 init_int nyear
 init_int nt  //numero	de tiempos


 number ntime
 !! ntime=nyear*nt;

 init_int ntallas
 init_matrix mdatos(1,ntime,1,15)
 init_vector Tallas(1,ntallas)


 init_int N_ftc
 !!if(N_ftc>0)
 init_ivector nanos_ftc(1,N_ftc)
 init_matrix Ctot(1,N_ftc,1,ntallas)

 init_int N_ftr
 !!if(N_ftr>0)
 init_ivector nanos_ftr(1,N_ftr)
 init_matrix Nrecl(1,N_ftr,1,ntallas)

 init_int N_ftp
 !!if(N_ftp>0)
 init_ivector nanos_ftp(1,N_ftp)
 init_matrix Npel(1,N_ftp,1,ntallas)

 init_vector msex(1,ntallas)
 init_vector Wmed(1,ntallas)

 //!!ad_comm::change_datafile_name("Modtrim.ctl");
 init_number sigmaR
 init_vector dt(1,2)
 init_vector Par_bio(1,7)
 init_vector cv_Par_bio(1,7)
 init_int    fase_Linf 
 init_int    fase_k 
 init_int    fase_Lr 
 init_int    fase_sr
 init_int    fase_beta
 init_int    fase_M 
 init_int    fase_h 

  number log_Linf_prior
  number log_k_prior
  number log_Lr_prior
  number log_sr_prior
  number log_beta_prior
  number log_M_prior
  number log_h_prior
  
  
  !! log_Linf_prior = log(Par_bio(1));
  !! log_k_prior = log(Par_bio(2));
  !! log_Lr_prior = log(Par_bio(3));
  !! log_sr_prior = log(Par_bio(4));
  !! log_beta_prior= log(Par_bio(5));
  !! log_M_prior= log(Par_bio(6));
  !! log_h_prior= log(Par_bio(7));


 init_int    minedad // lag de tiempo S/R

 init_number L50prior
 init_number s1prior
 init_number s2prior
 init_vector cv_parsel(1,3) // CV de las priors Selectividad
 init_int fases_flo1 // fases Sel Flota
 init_int fases_flo2 // 
 init_int fases_flo3 //

 number log_L50prior
 number log_s1prior
 number log_s2prior

 !! log_L50prior = log(L50prior);
 !! log_s1prior = log(s1prior);
 !! log_s2prior = log(s2prior);
 
 init_number L50pel_prior //pelas
 init_number s1pel_prior
 init_vector cv_parselp(1,2) // CV de las priors Selectividad
 init_int fases_pel1 // fases Sel Flota
 init_int fases_pel2 // 
 

 number log_L50pel_prior
 number log_s1pel_prior
 
 !! log_L50pel_prior = log(L50pel_prior);
 !! log_s1pel_prior = log(s1pel_prior);
 
 init_number L50recl_prior //recla
 init_number s1recl_prior
 init_vector cv_parselr(1,2) // CV de las priors Selectividad
 init_int fases_recl1 // fases Sel Flota
 init_int fases_recl2 // 

 number log_L50recl_prior
 number log_s1recl_prior
 
 !! log_L50recl_prior = log(L50recl_prior);
 !! log_s1recl_prior = log(s1recl_prior);


 init_int    nbloques1
 init_vector ybloques1(1,nbloques1)

 init_int    nqbloques
 init_vector yqbloques(1,nqbloques)
 
 init_int    fases_qf //capturabilidad_Flota
 init_int    fases_qrecl //capturabilidad_reclas
 init_int    fases_qpel //capturabilidad_pelas
 init_int    fases_F //_mortalidad_por_pesc
 init_int    fases_devRt
 init_int    fases_devRs

 init_int    npbr
 init_vector pbr(1,npbr)
 init_int ntime_sim //Años_a_proyectar


INITIALIZATION_SECTION

  log_Linf       log_Linf_prior
  log_k          log_k_prior
  log_Lr         log_Lr_prior
  log_sr         log_sr_prior
  log_beta       log_beta_prior
  log_M          log_M_prior
  log_h          log_h_prior
  log_L50        log_L50prior 
  log_sigma1     log_s1prior 
  log_sigma2     log_s2prior
  log_L50pel     log_L50pel_prior 
  log_spel       log_s1pel_prior
  log_L50recl    log_L50recl_prior 
  log_srecl      log_s1recl_prior

PARAMETER_SECTION
 
 init_vector log_L50(1,nbloques1,fases_flo1)  
 init_vector log_sigma1(1,nbloques1,fases_flo2)
 init_vector log_sigma2(1,nbloques1,fases_flo3)

 init_number log_L50pel(fases_pel1)  
 init_number log_spel(fases_pel2)
 init_number log_L50recl(fases_recl1)  
 init_number log_srecl(fases_recl2)

// parametros reclutamientos y mortalidades)
 init_number log_Rmed(1)
// init_bounded_dev_vector log_desv_Rt(1,nyear,-10,10,fases_devRt)
 init_bounded_vector log_desv_Rt(1,nyear,-10,10,fases_devRt)
 init_bounded_dev_vector log_desv_Rs(1,nt,-10,10,fases_devRs)
 init_vector log_F(1,ntime,fases_F) // log  mortalidad por pesca por flota

// capturabilidades
 init_vector log_qflo(1,nqbloques,fases_qf)
 init_number log_q_pel(fases_qpel)
 init_number log_q_recl(fases_qrecl)

// Crecimiento
 init_number log_Linf(fase_Linf)
 init_number log_k(fase_k)
 init_number log_Lr(fase_Lr)
 init_number log_sr(fase_sr)
 init_number log_beta(fase_beta)

// M y h
 init_number log_M(fase_M)
 init_number log_h(fase_h)

//---------------------------------------------------------------------------------
//Defino las variables de estado 
 vector BMflo(1,ntime)
 vector Brec(1,ntime)
 vector pred_CPUE(1,ntime);
 vector pred_Desemb(1,ntime);
 vector pred_Brecl(1,ntime);
 vector pred_Bpel(1,ntime);
 vector likeval(1,15);
 vector Neq(1,ntallas);

 vector Rpred(1,ntime);
 vector Unos_year(1,ntime);
 vector Unos_tallas(1,ntallas);
 vector delta(1,ntallas)
 vector Lesp(1,ntallas)
 vector sigmaL(1,ntallas)
 vector pre(1,ntallas)
 vector fobs_marg(1,ntallas)
 vector fpred_marg(1,ntallas)

 vector BDo(1,ntime);
 vector yrs(1,ntime)
 vector Desemb(1,ntime);
 vector CPUE(1,ntime);
 vector Brecl(1,ntime);
 vector Bpel(1,ntime);
 vector Reclutas(1,ntime)
 vector penalty(1,5)
 vector Lmed_obs(1,N_ftc)
 vector Lmed_pred(1,N_ftc)
 vector Lmed_obs_r(1,N_ftr)
 vector Lmed_pred_r(1,N_ftr)
 vector Lmed_obs_p(1,N_ftp)
 vector Lmed_pred_p(1,N_ftp)
 vector BDo_din(1,ntime)
 vector pri(1,14)
 vector years_ftc(1,N_ftc)
 vector years_ftr(1,N_ftr)
 vector years_ftp(1,N_ftp)

 matrix cv_index(1,4,1,ntime)
 matrix nm(1,3,1,ntime)
 matrix N0_din(1,ntime,1,ntallas)

 matrix S1(1,nbloques1,1,ntallas)
 matrix Sel(1,ntime,1,ntallas)
 vector Sel_recl(1,ntallas)//
 vector Sel_pel(1,ntallas)

 matrix F(1,ntime,1,ntallas)
 matrix Z(1,ntime,1,ntallas)
 matrix S(1,ntime,1,ntallas)


 matrix N(1,ntime,1,ntallas)

 matrix NM(1,ntime,1,ntallas)
 matrix NMD(1,ntime,1,ntallas)
 matrix NDv(1,ntime,1,ntallas)
 matrix Nrec(1,ntime,1,ntallas)
 matrix NVflo(1,ntime,1,ntallas)

 matrix pred_Ctot(1,ntime,1,ntallas)
 matrix pred_Nrecl(1,ntime,1,ntallas);
 matrix pred_Npel(1,ntime,1,ntallas);
 vector No(1,ntallas)


 matrix T(1,ntallas,1,ntallas)
 matrix Diag(1,ntallas,1,ntallas)
 matrix Id(1,ntallas,1,ntallas)

 matrix NMDv(1,ntime,1,ntallas)
 matrix pobs_recl(1,N_ftr,1,ntallas)
 matrix ppred_recl(1,N_ftr,1,ntallas)
 matrix pobs_pel(1,N_ftp,1,ntallas)
 matrix ppred_pel(1,N_ftp,1,ntallas)
 matrix pobs(1,N_ftc,1,ntallas)
 matrix ppred(1,N_ftc,1,ntallas)


 number suma1
 number suma2
 number suma3
 number suma4
 number suma5
 number suma6
 number suma7
 number suma8
 number suma9
 number suma10
  
 number So
 number alfa
 number beta

 number Linf
 number k
 number M
 number Lr
 number sr
 number h

 number BDp
 number Npplus
 number Bp_anch 

 number nm1;
 number cuenta1;
 number alfa_sr;
 number beta_sr;
 number pF

 vector Np(1,ntallas)
 vector Zpbr(1,ntallas)
 vector Fpbr(1,ntallas)
 vector Sp(1,ntallas)
 vector BT(1,ntime) // 

 matrix Bp(1,npbr,1,ntime_sim)
 vector CTPp(1,ntallas)
 matrix Rpp(1,npbr,1,ntime_sim)
 number SSBo

 
 objective_function_value f
  
 sdreport_vector BD(1,ntime) // 
 sdreport_vector Fmort(1,ntime) // 
 sdreport_vector RPRlp(1,ntime) // 
 sdreport_vector RPRdin(1,ntime) // 


PRELIMINARY_CALCS_SECTION


 yrs=column(mdatos,1);
 Desemb=column(mdatos,13);
 CPUE=column(mdatos,11);
 Brecl=column(mdatos,3);
 Bpel=column(mdatos,6);
 

 nm(1)=column(mdatos,5); //nm_r
 nm(2)=column(mdatos,8); //nm_p
 nm(3)=column(mdatos,15); // nm_c
 
 cv_index(1)=column(mdatos,4); //cv_r
 cv_index(2)=column(mdatos,7); //cv_p
 cv_index(3)=column(mdatos,12); // cv_cpue
 cv_index(4)=column(mdatos,14); //cv_y

 Unos_tallas=1;// lo uso en operaciones matriciales con tallas
 Unos_year=1;// lo uso en operaciones matriciales con el año


RUNTIME_SECTION
  convergence_criteria 1.e-1,1.e-01,1.e-03,1e-3,1e-5
  maximum_function_evaluations 100,100,200,3000,3500

PROCEDURE_SECTION
// se listan las funciones que contienen los calculos
 Eval_Trans_talla_talla();
 Eval_selectividad();
 Eval_mortalidades();
 Eval_abundancia();
 Eval_biomasas();
 Eval_capturas_predichas();
 Eval_indices();
 Eval_logverosim();
 Eval_funcion_objetivo();

//-----------------------------------------------------------------
FUNCTION Eval_Trans_talla_talla


 Linf=exp(log_Linf);
 k=exp(log_k);
 Lr=exp(log_Lr);
 sr=exp(log_sr);
 beta=exp(log_beta);
 M=exp(log_M);
 h=exp(log_h);


 int i, j;
 
// matriz de transicion modelo normal

  delta=(Linf-Tallas)*(1-mfexp(-k));// incremento en tallas

  for (i=1;i<=ntallas;i++){
   if(delta(i)<0){delta(i)=0.01;}
  }
  
  Lesp=Tallas+delta; // talla esperada luego del crecimiento
  sigmaL=delta*beta;

 
  for (i=1;i<=ntallas;i++){
    for (j=1;j<=ntallas;j++){
      if(i==j){
         T(i,j)=1.0;}}
   }


  for (i=1;i<=ntallas;i++){

    for (j=1;j<=ntallas;j++){
     if(sigmaL(i)>0){
     T(i,j)=mfexp(-0.5*square((Lesp(i)-Tallas(j))/sigmaL(i)));}}
   }


  for (j=1;j<=ntallas;j++){
  T(j)/=sum(T(j));
  } 

//----------------------------------------------------------------------

FUNCTION Eval_selectividad
 int i,j;
  //<<2<<endl;

 // FLOTA...................

 for (j=1;j<=nbloques1;j++){

 S1(j)=exp(-0.5*square(Tallas-exp(log_L50(j)))/square(exp(log_sigma1(j))));


    for (i=1;i<=ntallas;i++){

      if(Tallas(i)>=exp(log_L50(j))){
      S1(j,i)= exp(-0.5*square(Tallas(i)-exp(log_L50(j)))/square(exp(log_sigma2(j))));
      }

 }}
 

   for (i=1;i<=ntime;i++){
      for (j=1;j<=nbloques1;j++){
              if (yrs(i)>=ybloques1(j)){
                Sel(i)=S1(j);}
       }}


// RECLAS
 Sel_recl=1/(1+exp(-log(19)*(Tallas-exp(log_L50recl))/exp(log_srecl)));

// Pelas
 Sel_pel=1/(1+exp(-log(19)*(Tallas-exp(log_L50pel))/exp(log_spel)));



FUNCTION Eval_mortalidades
  //<<3<<endl;

 Fmort=exp(log_F);
 for (int j=1;j<=ntime;j++){
  F(j)=Sel(j)*Fmort(j);
  }

 Z=F+M;

 S=mfexp(-1.0*Z);


FUNCTION Eval_abundancia
   //<<4<<endl;

 int i, j;

  Lr=Par_bio(3);
  sr=Par_bio(4);

  if (active(log_Lr)){Lr=mfexp(log_Lr);}
  if (active(log_sr)){sr=mfexp(log_sr);}


// genero la composicion de tallas del reclutamiento
  pre=exp(-0.5*square((Tallas-Lr)/sr));
  pre/=sum(pre);


// genero la poblacion en equilibrio virginal de LP;

// matriz identidad
  for (int j=1;j<=ntallas;j++){
  Id(j,j)=1;}

  No=pre*exp(log_Rmed)*inv(Id-(exp(-M)*Id)*T);

  SSBo=sum(elem_prod(No*mfexp(-dt(1)*M),elem_prod(Wmed,msex)));
  alfa_sr=4*h*exp(log_Rmed)/(5*h-1);//
  beta_sr=(1-h)*SSBo/(5*h-1);// Reclutamiento

// -----------------primer año y primer trimestre
  Reclutas(1)=mfexp(log_Rmed+log_desv_Rt(1)+log_desv_Rs(1));
  Rpred(1)=mfexp(log_Rmed);


// Por defecto genero una estructura inicial en torno a Z del primer año y trimestre;

  for (int j=1;j<=ntallas;j++){
  Diag(j,j)=exp(-Z(1,j));}
  
  Neq=pre*exp(log_Rmed)*inv(Id-Diag*T);
  N(1)=Neq;

  NMD(1)=elem_prod(elem_prod(N(1),mfexp(-dt(1)*Z(1))),msex);
  BD(1)=sum(elem_prod(Wmed,NMD(1)));


// --------------------dinamica anual

  N0_din(1)=No;
  BDo_din(1)=sum(elem_prod(elem_prod(N0_din(1)*exp(-dt(1)*M),msex),Wmed));
  int cuenta;
  cuenta=1;

 
  for (i=1;i<=nyear;i++){

      for (j=1;j<=nt;j++){

    if(cuenta>1){
        Reclutas(cuenta)=mfexp(log_Rmed+log_desv_Rt(i)+log_desv_Rs(j));
        Rpred(cuenta)=Reclutas(cuenta);

        if(cuenta>minedad){
        Rpred(cuenta)=(alfa_sr*BD(cuenta-minedad)/(beta_sr+BD(cuenta-minedad)));
        Reclutas(cuenta)=Rpred(cuenta)*mfexp(log_desv_Rt(i) + log_desv_Rs(j));}

        N(cuenta)=(elem_prod(N(cuenta-1),S(cuenta-1)))*T+pre*Reclutas(cuenta);
        NMD(cuenta)=elem_prod(elem_prod(N(cuenta),mfexp(-dt(1)*Z(cuenta))),msex);
        BD(cuenta)=sum(elem_prod(Wmed,NMD(cuenta)));
      
        N0_din(cuenta)=(N0_din(cuenta-1)*exp(-M))*T+pre*Reclutas(cuenta); //------------nueva linea,.. N0_dinamino
        BDo_din(cuenta)=sum(elem_prod(elem_prod(N0_din(cuenta)*exp(-dt(1)*M),msex),Wmed));
   }
     
      cuenta+=1;

  }} //



FUNCTION Eval_biomasas
  //<<5<<endl;
 
 NVflo=elem_prod(elem_prod(N,mfexp(-dt(2)*(Z))),Sel);


// vectores de biomasas derivadas
 BMflo=Wmed*trans(NVflo);
 BT=Wmed*trans(N);

 //BDo=sum(elem_prod(No(nedades),Wmed));
 RPRlp=BD/SSBo;
 RPRdin=elem_div(BD,BDo_din);

FUNCTION Eval_capturas_predichas
  //<<6<<endl;

// matrices de capturas predichas por edad y año
 pred_Ctot=elem_prod(elem_div(F,Z),elem_prod(1.-S,N));

// vectores de desembarques predichos por año
 pred_Desemb=Wmed*trans(pred_Ctot);


 for (int i=1;i<=ntime;i++){

 pred_Nrecl(i)=elem_prod(N(i),Sel_recl); // N del reclas
 pred_Npel(i)=elem_prod(N(i),Sel_pel); // N del pelaces

 pred_Brecl(i)=exp(log_q_recl)*sum(elem_prod(Wmed,pred_Nrecl(i))); // biomasa del reclas
 pred_Bpel(i)=exp(log_q_pel)*sum(elem_prod(Wmed,pred_Npel(i))); // biomasa del pelaces

 }

 if(N_ftc>0){
 for (int i=1;i<=N_ftc;i++){
 pobs(i)=Ctot(i)/sum(Ctot(i));
 ppred(i)=pred_Ctot(nanos_ftc(i))/sum(pred_Ctot(nanos_ftc(i)));
 }
 Lmed_obs=Tallas*trans(pobs);
 Lmed_pred=Tallas*trans(ppred);
 }

 if(N_ftr>0){
 for (int i=1;i<=N_ftr;i++){
 pobs_recl(i)=Nrecl(i)/sum(Nrecl(i));
 ppred_recl(i)= pred_Nrecl(nanos_ftr(i))/sum( pred_Nrecl(nanos_ftr(i)));
 }
 Lmed_pred_r=Tallas*trans(ppred_recl);
 Lmed_obs_r=Tallas*trans(pobs_recl);
 }

 if(N_ftp>0){
 for (int i=1;i<=N_ftp;i++){
 pobs_pel(i)=Npel(i)/sum(Npel(i));
 ppred_pel(i)= pred_Npel(nanos_ftr(i))/sum( pred_Npel(nanos_ftr(i)));
 }
 Lmed_pred_p=Tallas*trans(ppred_pel);
 Lmed_obs_p=Tallas*trans(pobs_pel);
 }

FUNCTION Eval_indices



   for (int i=1;i<=ntime;i++){
      for (int j=1;j<=nqbloques;j++){
              if (yrs(i)>=yqbloques(j)){
                 pred_CPUE(i)=exp(log_qflo(j))*BMflo(i);}
       }
   }


FUNCTION Eval_logverosim


 int i;

 suma1=0; suma2=0; suma3=0; 

 for (i=1;i<=ntime;i++)
 {

  if (Brecl(i)>0){
    suma1+=square(log(Brecl(i)/pred_Brecl(i))*1/cv_index(1,i));}

  if (Bpel(i)>0){
    suma2+=square(log(Bpel(i)/pred_Bpel(i))*1/cv_index(2,i));}

  if (CPUE(i)>0){
    suma3+=square(log(CPUE(i)/pred_CPUE(i))*1/cv_index(3,i));}

 }


FUNCTION Eval_funcion_objetivo
  //<<8<<endl;

 suma4=0;suma5=0;suma6=0; penalty=0;

 likeval(1)=0.5*suma1;//reclas
 likeval(2)=0.5*suma2;// pelaces
 likeval(3)=0.5*suma3;// cpue
 likeval(4)=0.5*norm2(elem_div(log(elem_div(Desemb,pred_Desemb)),cv_index(4)));// desemb

 if(N_ftc>0){
 for (int i=1;i<=N_ftc;i++){
 suma4+=-nm(3,nanos_ftc(i))*sum(elem_prod(pobs(i),log(ppred(i)+1e-10)));
 }}

 if(N_ftr>0){
 for (int i=1;i<=N_ftr;i++){
 suma5+=-nm(1,nanos_ftr(i))*sum(elem_prod(pobs_recl(i),log(ppred_recl(i)+1e-10)));
 }}

 if(N_ftp>0){
 for (int i=1;i<=N_ftp;i++){
 suma6+=-nm(2,nanos_ftp(i))*sum(elem_prod(pobs_pel(i),log(ppred_pel(i)+1e-10)));
 }}
 
 likeval(5)=suma4;//
 likeval(6)=suma5;//
 likeval(7)=suma6;//

// Reclutas
   if(active(log_desv_Rt)){
   likeval(8)=1./(2*square(sigmaR))*(norm2(log_desv_Rt)+norm2(log_desv_Rs));}
 

 // distribuciones a priori
 
  pri(1)=0.5*square((log_Linf_prior-log_Linf)/cv_Par_bio(1));
  pri(2)=0.5*square((log_k_prior-log_k)/cv_Par_bio(2));
  pri(3)=0.5*square((log_Lr_prior-log_Lr)/cv_Par_bio(3));
  pri(4)=0.5*square((log_sr_prior-log_sr)/cv_Par_bio(4));
  pri(5)=0.5*square((log_beta_prior-log_beta)/cv_Par_bio(5));
  pri(6)=0.5*square((log_M_prior-log_M)/cv_Par_bio(6));
  pri(7)=0.5*square((log_h_prior-log_h)/cv_Par_bio(7));

  pri(8)=0.5*norm2((log_L50prior-log_L50)/cv_parsel(1));
  pri(9)=0.5*norm2((log_s1prior-log_sigma1)/cv_parsel(2));
  pri(10)=0.5*norm2((log_s2prior-log_sigma2)/cv_parsel(3));

  pri(11)=0.5*square((log_L50pel-log_L50pel_prior)/cv_parselp(1));
  pri(12)=0.5*square((log_spel-log_s1pel_prior)/cv_parselp(2));
  pri(13)=0.5*square((log_L50recl-log_L50recl_prior)/cv_parselr(1));
  pri(14)=0.5*square((log_srecl-log_s1recl_prior)/cv_parselr(2));

   if (active(log_F)){
   pF=1000*norm2(log_F-mean(log_F));}

   f=sum(likeval)+sum(penalty)+pF+sum(pri);

   if(last_phase){
   f=sum(likeval)+sum(pri);}
 
REPORT_SECTION
   //<<9<<endl;

 report << "Yrs" << endl;
 report << yrs << endl;
 report << "CPUE_obs_pred" << endl;
 report << CPUE << endl;
 report << pred_CPUE << endl;
 report << "RECLAS_obs_pred" << endl;
 report << Brecl << endl;
 report << pred_Brecl << endl;
 report << "PELACES_obs_pred" << endl;
 report << Bpel << endl;
 report << pred_Bpel << endl;
 report << "Desemb_obs_pred" << endl;
 report << Desemb << endl;
 report << pred_Desemb << endl;
 report << "Lmed_obs_pred_FLOTA" << endl;
 for (int i=1;i<=N_ftc;i++){
 years_ftc(i)=yrs(nanos_ftc(i));
 }
 report << years_ftc << endl;
 report << Lmed_obs << endl;
 report << Lmed_pred << endl;
 report << "Lmed_obs_pred_RECLAS" << endl;
 for (int i=1;i<=N_ftr;i++){
 years_ftr(i)=yrs(nanos_ftr(i));
 }
 report << years_ftr << endl;
 report << Lmed_obs_r << endl;
 report << Lmed_pred_r << endl;
 report << "Lmed_obs_pred_PELACES" << endl;
 for (int i=1;i<=N_ftp;i++){
 years_ftp(i)=yrs(nanos_ftp(i));
 }
 report << years_ftp << endl;
 report << Lmed_obs_p << endl;
 report << Lmed_pred_p << endl;
 report << "Biomasa_desovante" << endl;
 report << BD << endl;
 report << "Biomasa_total" << endl;
 report << BT << endl;
 report << "Biomasa_explotable" << endl;
 report << BMflo << endl;
 report << "Reclutamiento" << endl;
 report << Reclutas<< endl;
 report << Rpred<< endl;
 report << "F " << endl;
 report << Fmort << endl;
 report << "log_anomalia_R_t " << endl;
 report << log_desv_Rt << endl;
 report << "log_anomalia_R_s " << endl;
 report << log_desv_Rs << endl;

 report<<"Tallas"<<endl;
 report<<Tallas<<endl;

 report<<"Frec_tallas_marg_FLOTA"<<endl;

 report<<colsum(pobs)<<endl;
 report<<colsum(ppred)<<endl;

 report<<"Frec_tallas_marg_RECLAS"<<endl;

 report<<colsum(pobs_recl)<<endl;
 report<<colsum(ppred_recl)<<endl;

 report<<"Frec_tallas_marg_PELACES"<<endl;

 report<<colsum(pobs_pel)<<endl;
 report<<colsum(ppred_pel)<<endl;


 report<<"N_talla"<<endl;
 report<<N<<endl;

 report<<"Sel_Flota"<<endl;
 report<<S1<<endl;

 report<<"Sel_Reclas"<<endl;
 report<<Sel_recl<<endl;

 report<<"Sel_Pelaces"<<endl;
 report<<Sel_pel<<endl;

 report << "Prop_obs_FLOTA" << endl;
 report << pobs<< endl;
 report << "Prop_pred_FLOTA" << endl;
 report << ppred<< endl;

 report << "Prop_obs_RECLAS" << endl;
 report << pobs_recl<< endl;
 report << "Prop_pred_RECLAS" << endl;
 report << ppred_recl<< endl;

 report << "Prop_obs_PELACES" << endl;
 report << pobs_pel<< endl;
 report << "Prop_pred_PELACES" << endl;
 report << ppred_pel<< endl;

 report << "BD0_dinam" << endl;
 report << BDo_din << endl;
 report << "BD0_LP" << endl;
 report << SSBo << endl;

 report << "RPR_LP_Din" << endl;
 report << RPRlp << endl;
 report << RPRdin << endl;

 report << "Lesp_sd" << endl;
 report << Lesp << endl;
 report << sigmaL << endl;
 report << "DistrL_recl" << endl;
 report << Tallas<< endl;
 report << pre<< endl;
 report << "Mtrans" << endl;
 report << T << endl;
 report << "Proy_No" << endl;
  No=pre*exp(log_Rmed);
  for (int j=2;j<=20;j++){
  No=(No*exp(-1.*M))*T; // proyeccion poblacional virgen
  report <<No<<endl;
  }



//-------------------------------------------------------------------
// ESTIMA nm y CV

  suma1=0; suma2=0;nm1=1;cuenta1=0;

  for (int i=1;i<=N_ftc;i++){ //

   if (sum(pobs(i))>0){
      suma1=sum(elem_prod(ppred(i),1-ppred(i)));
      suma2=norm2(pobs(i)-ppred(i));
      nm1=nm1*suma1/suma2;
      cuenta1+=1;
   }}

 report << "nm_f" <<endl;
 report <<pow(nm1,1/cuenta1)<< endl;

 report << "Componentes de verosimilitud" << endl;
 report << "Brecl   Bpela   CPUE  Desem  pflo  precl   ppel  devR"<< endl;
 report << likeval << endl;
 report << "Verosimilitud total" << endl;
 report << sum(likeval) << endl;
 report << "-----------------------------------------------" << endl;
  report << "Parametros biologicos y log_prior" << endl;
  report<<"Linf="<<" "<<exp(log_Linf)<<endl;
  report<<"k   ="<<" "<<exp(log_k)<<endl;
  report<<"Lr  ="<<" "<<exp(log_Lr)<<endl;
  report<<"sr  ="<<" "<<exp(log_sr)<<endl;
  report<<"beta="<<" "<<exp(log_beta)<<endl;
  report<<"M   ="<<" "<<exp(log_M)<<endl;
  report<<"h   ="<<" "<<exp(log_h)<<endl;


FINAL_SECTION

 time(&finish);
 elapsed_time=difftime(finish,start);
 hour=long(elapsed_time)/3600;
 minute=long(elapsed_time)%3600/60;
 second=(long(elapsed_time)%3600)%60;
 cout<<endl<<endl<<"*********************************************"<<endl;
 cout<<"--Start time:  "<<ctime(&start)<<endl;
 cout<<"--Finish time: "<<ctime(&finish)<<endl;
 cout<<"--Runtime: ";
 cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
 cout<<"*********************************************"<<endl;


