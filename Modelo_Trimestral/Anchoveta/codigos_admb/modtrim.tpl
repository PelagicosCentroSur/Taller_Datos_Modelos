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

 init_int nedades
 init_int ntallas
 init_matrix mdatos(1,ntime,1,15)
 init_vector Tallas(1,ntallas)
 //!!cout<<nyear<<endl;exit(1);
 
 init_matrix Ctot(1,ntime,1,ntallas)
 init_matrix Nrecl(1,ntime,1,ntallas)
 init_matrix Npel(1,ntime,1,ntallas)

 init_vector msex(1,ntallas)
 init_vector Wmed(1,ntallas)
// !!cout<<msex<<endl;exit(1);

 !!ad_comm::change_datafile_name("Modtrim.ctl");
 init_number sigmaR
 init_vector dt(1,2)
 init_vector Par_bio(1,7)
 init_vector cv_Par_bio(1,7)
//!!cout<<cv_Par_bio<<endl;exit(1);

// !!cout<<dt<<endl;
// !!cout<<cv_Par_bio<<endl;exit(1);

  number log_Lr_prior
  number log_sr_prior
  number log_b_prior
  number log_beta_prior
  
  
  !! log_Lr_prior = log(Par_bio(3));
  !! log_sr_prior = log(Par_bio(4));
  !! log_beta_prior= log(Par_bio(5));
  !! log_b_prior = log(bprior);


 init_int    minedad
 init_number bprior //Hiperestabilidad


 init_number L50prior
 init_number s1prior
 init_number s2prior
 init_int opt_sel2 // opcion domo en flota

 number log_L50prior
 number log_s1prior
 number log_s2prior

 !! log_L50prior = log(L50prior);
 !! log_s1prior = log(s1prior);
 !! log_s2prior = log(s2prior);
 
 init_number L50pel_prior //pelas
 init_number s1pel_prior

 number log_L50pel_prior
 number log_s1pel_prior
 
 !! log_L50pel_prior = log(L50pel_prior);
 !! log_s1pel_prior = log(s1pel_prior);
 
 init_number L50recl_prior //recla
 init_number s1recl_prior

 number log_L50recl_prior
 number log_s1recl_prior
 
 !! log_L50recl_prior = log(L50recl_prior);
 !! log_s1recl_prior = log(s1recl_prior);

 init_int    nbloques1
 init_vector ybloques1(1,nbloques1)

 init_int    nqbloques
 init_vector yqbloques(1,nqbloques)
 

 init_int    opt_qf //capturabilidad_Flota
 init_int    opt1_fase // selectividad flota

 init_int opt1_fase_srecl // selectividad del reclas
 init_int  opt_qrecl //capturabilidad_reclas

 init_int opt1_fase_spel // selectividad pelas
 init_int  opt_qpel //capturabilidad_pelas

 init_int    opt_Lr 
 init_int    opt_sr
 init_int    opt_beta
 init_int    opt_F //_mortalidad_por_pesc
 init_int    opt_devRt
 init_int    opt_devs
 init_int  opt_devNo//Condicion inicial (Si no estima significa poblaciónen equilibrio)
 init_int    opt_bpow //hiperestabilidad
 init_int    npbr
 init_vector pbr(1,npbr)
 init_int ntime_sim //Años_a_proyectar


INITIALIZATION_SECTION

  log_Lr         log_Lr_prior
  log_sr         log_sr_prior
  log_L50        log_L50prior 
  log_sigma1     log_s1prior 
  log_sigma2     log_s2prior
  log_b          log_b_prior 
  log_beta       log_beta_prior 
  log_L50pel     log_L50pel_prior 
  log_spel       log_s1pel_prior
  log_L50recl    log_L50recl_prior 
  log_srecl      log_s1recl_prior
  log_b          0

PARAMETER_SECTION


// selectividad paramétrica a la talla común
// init_bounded_vector log_L50f(1,nbloques1,-5,8,opt1_fase)  
 
 init_vector log_L50(1,nbloques1,opt1_fase)  
 init_vector log_sigma1(1,nbloques1,opt1_fase)
 init_vector log_sigma2(1,nbloques1,opt_sel2)

 init_number log_L50pel(opt1_fase_spel)  
 init_number log_spel(opt1_fase_spel)
 init_number log_L50recl(opt1_fase_srecl)  
 init_number log_srecl(opt1_fase_srecl)



// parametros reclutamientos y mortalidades)
 init_number log_Rmed(1)
 init_bounded_dev_vector log_desv_Rt(1,nyear,-10,10,opt_devRt)
 init_bounded_dev_vector log_desv_Rs(1,nt,-10,10,opt_devs)
 //init_bounded_vector log_desv_No(1,nedades,-10,10,opt_devNo)
 init_bounded_vector log_F(1,ntime,-20,2,opt_F) // log  mortalidad por pesca por flota

// capturabilidades
 init_vector log_qflo(1,nqbloques,opt_qf)
 init_number log_q_pel(opt_qpel)
 init_number log_q_recl(opt_qrecl)


 init_number log_b(opt_bpow)

// Crecimiento
 init_number log_Lr(opt_Lr)
 init_number log_sr(opt_sr)
 init_number log_beta(opt_beta)

//---------------------------------------------------------------------------------
//Defino las variables de estado 
 vector BMflo(1,ntime)
 vector Brec(1,ntime)
 vector pred_CPUE(1,ntime);
 vector pred_Desemb(1,ntime);
 vector pred_Brecl(1,ntime);
 vector pred_Bpel(1,ntime);
 matrix pred_Nrecl(1,ntime,1,ntallas);
 matrix pred_Npel(1,ntime,1,ntallas);
 vector likeval(1,15);
 vector Neq(1,ntallas);

 vector Rpred(1,ntime);
 vector Unos_edad(1,nedades);
 vector Unos_year(1,ntime);
 vector Unos_tallas(1,ntallas);
 vector delta(1,ntallas)
 vector Lesp(1,ntallas)
 vector sigmaL(1,ntallas)
 vector pre(1,ntallas)
 vector fobs_marg(1,ntallas)
 vector fpred_marg(1,ntallas)
 vector G(1,ntallas)

 vector mu_edad(1,nedades)
 vector sigma_edad(1,nedades)
 vector BDo(1,ntime);
 matrix No(1,nedades,1,ntallas)
 vector prior(1,7)
 vector yrs(1,ntime)
 vector Desemb(1,ntime);
 vector CPUE(1,ntime);
 vector Brecl(1,ntime);
 vector Bpel(1,ntime);
 vector Lmed_obs(1,ntime)
 vector Lmed_pred_r(1,ntime)
 vector Lmed_obs_r(1,ntime)
 vector Lmed_pred_p(1,ntime)
 vector Lmed_obs_p(1,ntime)
 vector Lmed_pred(1,ntime)
 vector edades(1,nedades)
 vector Reclutas(1,ntime)
 vector penalty(1,5)
 matrix N0_din(1,ntime,1,ntallas)
 vector BDo_din(1,ntime)

 matrix cv_index(1,4,1,ntime)
 matrix nm(1,3,1,ntime)

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

 matrix pobs(1,ntime,1,ntallas)
 matrix ppred(1,ntime,1,ntallas)

 matrix T(1,ntallas,1,ntallas)

 matrix Nv(1,ntime,1,nedades)
 matrix NMDv(1,ntime,1,ntallas)
 matrix pobs_recl(1,ntime,1,ntallas)
 matrix ppred_recl(1,ntime,1,ntallas)
 matrix pobs_pel(1,ntime,1,ntallas)
 matrix ppred_pel(1,ntime,1,ntallas)


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
 number Linfh
 number M
 number Lr
 number sr
 number Lm
 number Rm
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

 matrix Bp(1,npbr,1,ntime_sim)
 vector CTPp(1,ntallas)
 matrix Rpp(1,npbr,1,ntime_sim)

 
 objective_function_value f
  
 sdreport_vector BD(1,ntime) // 
 sdreport_vector BT(1,ntime) // 
 sdreport_vector RPRlp(1,ntime) // 
 sdreport_number SSBo
 sdreport_matrix Yp(1,npbr,1,ntime_sim)


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


 Linf=Par_bio(1);
 k=Par_bio(2);
 M=Par_bio(6);
 h=Par_bio(7);

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

 if(last_phase()){Eval_CTP();}


//-----------------------------------------------------------------
FUNCTION Eval_Trans_talla_talla

  Linf=Par_bio(1);
  k=Par_bio(2);
  beta=Par_bio(5);


//  if(active(log_k)){k=mfexp(log_k);}
  if(active(log_beta)){beta=mfexp(log_beta);}

 int i, j;
 
// matriz de transicion modelo normal

  delta=(Linf-Tallas)*(1-mfexp(-k));// incremento en tallas
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

 F=elem_prod(Sel,outer_prod(mfexp(log_F),Unos_tallas));

 Z=F+M;

 S=mfexp(-1.0*Z);


FUNCTION Eval_abundancia
 int i, j;

  Lr=Par_bio(3);
  sr=Par_bio(4);

  if (active(log_Lr)){Lr=mfexp(log_Lr);}
  if (active(log_sr)){sr=mfexp(log_sr);}


// genero la composicion de tallas del reclutamiento
  pre=exp(-0.5*square((Tallas-Lr)/sr));
  pre/=sum(pre);


// genero la poblacion en equilibrio virginal de LP;
  No(1)=pre*exp(log_Rmed);
  for (int j=2;j<=nedades;j++){
  No(j)=(No(j-1)*exp(-1.*M))*T+pre*exp(log_Rmed); // proyeccion poblacional virgen
  }

  SSBo=sum(elem_prod(No(nedades)*mfexp(-dt(1)*M),elem_prod(Wmed,msex)));
  alfa_sr=4*h*exp(log_Rmed)/(5*h-1);//
  beta_sr=(1-h)*SSBo/(5*h-1);// Reclutamiento



// -----------------primer año y primer trimestre
  Reclutas(1)=mfexp(log_Rmed+log_desv_Rt(1)+log_desv_Rs(1));
  Rpred(1)=mfexp(log_Rmed);


// Por defecto genero una estructura inicial en torno a Z del primer año y trimestre;
  Neq=pre*Rpred(1);
  for (j=1;j<=nedades;j++)
   {Neq=elem_prod(Neq,exp(-1.*Z(1)))*T+pre*Rpred(1);}
    N(1)=Neq;

  NMD(1)=elem_prod(elem_prod(N(1),mfexp(-dt(1)*Z(1))),msex);

// En el caso que el año 1  y trimestre 1 sea virginal

  if(opt_devNo<0){
  N(1)=No(nedades); 
  NMD(1)=elem_prod(N(1)*mfexp(-dt(1)*M),msex);}
  
  BD(1)=sum(elem_prod(Wmed,NMD(1)));


// --------------------dinamica anual

  N0_din(1)=No(nedades);
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
 
 NVflo=elem_prod(elem_prod(N,mfexp(-dt(2)*(Z))),Sel);


// vectores de biomasas derivadas
 BMflo=Wmed*trans(NVflo);
 BT=Wmed*trans(N);

 //BDo=sum(elem_prod(No(nedades),Wmed));
 RPRlp=BD/SSBo;


FUNCTION Eval_capturas_predichas

// matrices de capturas predichas por edad y año
 pred_Ctot=elem_prod(elem_div(F,Z),elem_prod(1.-S,N));

// vectores de desembarques predichos por año
 pred_Desemb=Wmed*trans(pred_Ctot);

// matrices de proporcion de capturas por talla y año
 pobs=elem_div(Ctot,outer_prod(rowsum(Ctot+1e-10),Unos_tallas));
 ppred=elem_div(pred_Ctot,outer_prod(rowsum(pred_Ctot+1e-10),Unos_tallas));

 Lmed_pred=Tallas*trans(ppred);
 Lmed_obs=Tallas*trans(pobs);


FUNCTION Eval_indices



   for (int i=1;i<=ntime;i++){
      for (int j=1;j<=nqbloques;j++){
              if (yrs(i)>=yqbloques(j)){
                 pred_CPUE(i)=exp(log_qflo(j))*pow(BMflo(i),exp(log_b));}
       }
   }


 for (int i=1;i<=ntime;i++){

 pred_Nrecl(i)=elem_prod(N(i),Sel_recl); // N del reclas
 pred_Npel(i)=elem_prod(N(i),Sel_pel); // N del pelaces

 pred_Brecl(i)=exp(log_q_recl)*sum(elem_prod(Wmed,pred_Nrecl(i))); // biomasa del reclas
 pred_Bpel(i)=exp(log_q_pel)*sum(elem_prod(Wmed,pred_Npel(i))); // biomasa del pelaces

 }

 pobs_recl=elem_div(Nrecl,outer_prod(rowsum(Nrecl+1e-10),Unos_tallas));
 ppred_recl=elem_div(pred_Nrecl,outer_prod(rowsum(pred_Nrecl+1e-10),Unos_tallas));

 pobs_pel=elem_div(Npel,outer_prod(rowsum(Npel+1e-10),Unos_tallas));
 ppred_pel=elem_div(pred_Npel,outer_prod(rowsum(pred_Npel+1e-10),Unos_tallas));



 Lmed_pred_r=Tallas*trans(ppred_recl);
 Lmed_obs_r=Tallas*trans(pobs_recl);


 Lmed_pred_p=Tallas*trans(ppred_pel);
 Lmed_obs_p=Tallas*trans(pobs_pel);




FUNCTION Eval_logverosim
// esta funcion evalua el nucleo de las -log-verosimilitudes marginales para
// series con datos 0.
 int i;

 suma1=0; suma2=0; suma3=0; penalty=0;

 for (i=1;i<=ntime;i++)
 {

  if (Brecl(i)>0){
    suma2+=square(log(Brecl(i)/pred_Brecl(i))*1/cv_index(1,i));}

  if (Bpel(i)>0){
    suma3+=square(log(Bpel(i)/pred_Bpel(i))*1/cv_index(2,i));}

  if (CPUE(i)>0){
    suma1+=square(log(CPUE(i)/pred_CPUE(i))*1/cv_index(3,i));}

 }


FUNCTION Eval_funcion_objetivo // punto para resolver

 suma4=0;suma5=0;suma6=0; penalty=0;

 likeval(1)=0.5*suma1;//CPUE
 likeval(2)=0.5*suma2;// Breclas
 likeval(3)=0.5*suma3;// Bpelaces



 likeval(4)=0.5*norm2(elem_div(log(elem_div(Desemb,pred_Desemb)),cv_index(4)));// desemb
 likeval(5)=-1.*mean(nm(1))*sum(elem_prod(pobs_recl,log(ppred_recl)));
 likeval(6)=-1.*mean(nm(2))*sum(elem_prod(pobs_pel,log(ppred_pel)));


 for (int i=1;i<=ntime;i++){
 suma4+=-nm(3,i)*sum(elem_prod(pobs(i),log(ppred(i)+1e-10)));
 }

 for (int i=1;i<=ntime;i++){
 suma5+=-nm(1,i)*sum(elem_prod(pobs_recl(i),log(ppred_recl(i)+1e-10)));
 }

 for (int i=1;i<=ntime;i++){
 suma6+=-nm(2,i)*sum(elem_prod(pobs_pel(i),log(ppred_pel(i)+1e-10)));
 }


 likeval(7)=suma4;//
 likeval(8)=suma5;//
 likeval(9)=suma6;//



 //lognormal Ninicial y Reclutas
   if(active(log_desv_Rt)){
   likeval(10)=1./(2*square(sigmaR))*(norm2(log_desv_Rt)+norm2(log_desv_Rs));}

  /*
  if(active(log_desv_No)){
   likeval(9)=1./(2*square(sigmaR))*norm2(log_desv_No);}
  */
  
  if(active(log_Lr)){
   likeval(11)=1./(2*square(cv_Par_bio(3)))*square(log_Lr-log_Lr_prior);}

  if (active(log_F)){
   pF=1000*norm2(log_F-mean(log_F));}


  penalty=500*square(log_desv_Rt(nyear));


  f=(sum(likeval)+sum(penalty)+pF);

   if(last_phase){
   f=(sum(likeval)+sum(penalty));}
 


FUNCTION  Eval_CTP

//-----------------------------------------------------------------


  for (int i=1;i<=npbr;i++){ // ciclo de PBR

  Np=N(ntime);
  Sp=S(ntime);


  Fpbr=F(ntime)*pbr(i);//
  Zpbr=Fpbr+M;

  CTPp=elem_prod(elem_div(Fpbr,Zpbr),elem_prod(1.-exp(-1.*Zpbr),Np));

  for (int j=1;j<=ntime_sim;j++){ // ciclo de años


  if(j<=minedad){
  Np=(elem_prod(Np,Sp))*T+pre*(alfa_sr*BD(ntime-minedad+1)/(beta_sr+BD(ntime-minedad+1)));} //

  if(j>minedad){
  Np=(elem_prod(Np,Sp))*T+pre*(alfa_sr*Bp(i,j-minedad)/(beta_sr+Bp(i,j-minedad)));} //

  Rpp(i,j)=(alfa_sr*BD(ntime-minedad+1)/(beta_sr+BD(ntime-minedad+1)));

  Bp(i,j)=sum(elem_prod(elem_prod(Np,exp(-dt(1)*Zpbr)),elem_prod(msex,Wmed)));
  CTPp=elem_prod(elem_div(Fpbr,Zpbr),elem_prod(1.-exp(-1.*Zpbr),Np));
  Yp(i,j)=sum(elem_prod(CTPp,Wmed));
  Sp=exp(-1.*Zpbr);
  }}


REPORT_SECTION

 report << "Años" << endl;
 report << yrs << endl;
 report << "CPUE_obs & pred" << endl;
 report << CPUE << endl;
 report << pred_CPUE << endl;
 report << "RECLAS_obs & pred" << endl;
 report << Brecl << endl;
 report << pred_Brecl << endl;
 report << "PELACES_obs & pred" << endl;
 report << Bpel << endl;
 report << pred_Bpel << endl;
 report << "Desemb_obs & pred" << endl;
 report << Desemb << endl;
 report << pred_Desemb << endl;
 report << "Lmed_obs & pred FLOTA" << endl;
 report << Lmed_obs << endl;
 report << Lmed_pred << endl;
 report << "Lmed_obs & pred RECLAS" << endl;
 report << Lmed_obs_r << endl;
 report << Lmed_pred_r << endl;
 report << "Lmed_obs & pred PELACES" << endl;
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
 report << exp(log_F) << endl;

 report << "log_anomalia_R_t " << endl;
 report << log_desv_Rt << endl;
 report << "log_anomalia_R_s " << endl;
 report << log_desv_Rs << endl;



 report<<"Tallas"<<endl;
 report<<Tallas<<endl;

 report<<"Frecuencias de tallas marginales FLOTA"<<endl;


 fobs_marg=0;fpred_marg=0;
 for (int j=1;j<=ntime;j++){

    if(sum(pobs(j))>0){

    fobs_marg+=pobs(j); 
    fpred_marg+=ppred(j);
    }}

 report<<fobs_marg<<endl;
 report<<fpred_marg<<endl;


 report<<"Frecuencias de tallas marginales RECLAS"<<endl;


 fobs_marg=0;fpred_marg=0;
 for (int j=1;j<=ntime;j++){

    if(sum(pobs_recl(j))>0){

    fobs_marg+=pobs_recl(j); 
    fpred_marg+=ppred_recl(j);
    }}

 report<<fobs_marg<<endl;
 report<<fpred_marg<<endl;

 report<<"Frecuencias de tallas marginales PELACES"<<endl;

 fobs_marg=0;fpred_marg=0;
 for (int j=1;j<=ntime;j++){

    if(sum(pobs_pel(j))>0){

    fobs_marg+=pobs_pel(j); 
    fpred_marg+=ppred_pel(j);
    }}

 report<<fobs_marg<<endl;
 report<<fpred_marg<<endl;




 report<<"Abundancia a la talla"<<endl;
 report<<N<<endl;

 report<<"Selectividad a la talla Flota"<<endl;
 report<<Sel<<endl;

 report<<"Selectividad a la talla Reclas"<<endl;
 report<<Sel_recl<<endl;

 report<<"Selectividad a la talla Pelaces"<<endl;
 report<<Sel_pel<<endl;

 report << "Prop_obs FLOTA" << endl;
 report << pobs<< endl;
 report << "Prop_pred FLOTA" << endl;
 report << ppred<< endl;

 report << "Prop_obs RECLAS" << endl;
 report << pobs_recl<< endl;
 report << "Prop_pred RECLAS" << endl;
 report << ppred_recl<< endl;

 report << "Prop_obs PELACES" << endl;
 report << pobs_pel<< endl;
 report << "Prop_pred PELACES" << endl;
 report << ppred_pel<< endl;


 report << "BD virgen anual" << endl;
 report << BDo_din << endl;
 report << "BD virgen largo plazo" << endl;
 report << SSBo << endl;

 report << "Reducción del stock (BD/BDo) de LP y Dinámica" << endl;
 report << RPRlp << endl;
 report << elem_div(BD,BDo_din) << endl;


 report << " " << endl;

 report << "Crecimiento esperado para cada talla y su desviación" << endl;
 report << Lesp << endl;
 report << sigmaL << endl;


 report << "-----------------------------------------------" << endl;
 report << "Tallas y función de reclutamiento a la talla" << endl;
 report << Tallas<< endl;
 report << pre<< endl;
 
 report << "-----------------------------------------------" << endl;

 report << "P(talla-talla)" << endl;
 report << T << endl;

 report << "-----------------------------------------------" << endl;
 report << "bCPUE  Lr  Sr  beta  h " << endl;
 report << exp(log_b)<<" "<<exp(log_Lr)<<" "<<exp(log_sr)<<" "<<exp(log_beta)<<" "<<h<< endl;
 report << "-----------------------------------------------" << endl;
 report << "Proyeccion de la abundancia virginal" << endl;

 for (int j=2;j<=nedades;j++){
 No(j)=(No(j-1)*exp(-1.*M))*T; // proyeccion poblacional virgen
 }

 report << No << endl;

//-------------------------------------------------------------------
// ESTIMA nm y CV

  suma1=0; suma2=0;nm1=1;cuenta1=0;

  for (int i=1;i<=ntime;i++){ //

   if (sum(pobs(i))>0){
      suma1=sum(elem_prod(ppred(i),1-ppred(i)));
      suma2=norm2(pobs(i)-ppred(i));
      nm1=nm1*suma1/suma2;
      cuenta1+=1;
   }}

 report << "Tamaño muestra ideal flota " <<endl;
 report <<pow(nm1,1/cuenta1)<< endl;

 
 report << "-----------------------------------------------" << endl;
 report << "Biomasa desovante proyectada para cada multiplo de Flast" << endl;
 report << Bp << endl;
 report << "Capturas proyectadas para cada Fpbr" << endl;
 report << Yp << endl;
 report << "-----------------------------------------------" << endl;
 report << "Componentes de verosimilitud" << endl;
 report << "CPUE Brecl Bpel Des pobsrecl pobspel  Rt No Lr F"<< endl;
 report << likeval << endl;
 report << "Verosimilitud total" << endl;
 report << sum(likeval) << endl;






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


