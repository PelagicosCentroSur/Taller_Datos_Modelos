#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
 #include <admodel.h>
 #include <stdio.h>
 #include <time.h>
 time_t start,finish;
 long hour,minute,second;
 double elapsed_time;
 ofstream mcmc_report("mcmc.csv");
#ifdef DEBUG
  #include <chrono>
#endif
#include <admodel.h>
#ifdef USE_ADMB_CONTRIBS
#include <contrib.h>

#endif
  extern "C"  {
    void ad_boundf(int i);
  }
#include <modtrim3.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  adstring tmpstring;
  tmpstring=adprogram_name + adstring(".dat");
  if (argc > 1)
  {
    int on=0;
    if ( (on=option_match(argc,argv,"-ind"))>-1)
    {
      if (on>argc-2 || argv[on+1][0] == '-')
      {
        cerr << "Invalid input data command line option"
                " -- ignored" << endl;
      }
      else
      {
        tmpstring = adstring(argv[on+1]);
      }
    }
  }
  global_datafile = new cifstream(tmpstring);
  if (!global_datafile)
  {
    cerr << "Error: Unable to allocate global_datafile in model_data constructor.";
    ad_exit(1);
  }
  if (!(*global_datafile))
  {
    delete global_datafile;
    global_datafile=NULL;
  }
  nyear.allocate("nyear");
  nt.allocate("nt");
 ntime=nyear*nt;
  ntallas.allocate("ntallas");
  mdatos.allocate(1,ntime,1,15,"mdatos");
  Tallas.allocate(1,ntallas,"Tallas");
  N_ftc.allocate("N_ftc");
if(N_ftc>0)
  nanos_ftc.allocate(1,N_ftc,"nanos_ftc");
  Ctot.allocate(1,N_ftc,1,ntallas,"Ctot");
  N_ftr.allocate("N_ftr");
if(N_ftr>0)
  nanos_ftr.allocate(1,N_ftr,"nanos_ftr");
  Nrecl.allocate(1,N_ftr,1,ntallas,"Nrecl");
  N_ftp.allocate("N_ftp");
if(N_ftp>0)
  nanos_ftp.allocate(1,N_ftp,"nanos_ftp");
  Npel.allocate(1,N_ftp,1,ntallas,"Npel");
  msex.allocate(1,ntallas,"msex");
  Wmed.allocate(1,ntallas,"Wmed");
  sigmaR.allocate("sigmaR");
  dt.allocate(1,2,"dt");
  Par_bio.allocate(1,7,"Par_bio");
  cv_Par_bio.allocate(1,7,"cv_Par_bio");
  fase_Linf.allocate("fase_Linf");
  fase_k.allocate("fase_k");
  fase_Lr.allocate("fase_Lr");
  fase_sr.allocate("fase_sr");
  fase_beta.allocate("fase_beta");
  fase_M.allocate("fase_M");
  fase_h.allocate("fase_h");
 log_Linf_prior = log(Par_bio(1));
 log_k_prior = log(Par_bio(2));
 log_Lr_prior = log(Par_bio(3));
 log_sr_prior = log(Par_bio(4));
 log_beta_prior= log(Par_bio(5));
 log_M_prior= log(Par_bio(6));
 log_h_prior= log(Par_bio(7));
  minedad.allocate("minedad");
  L50prior.allocate("L50prior");
  s1prior.allocate("s1prior");
  s2prior.allocate("s2prior");
  cv_parsel.allocate(1,3,"cv_parsel");
  fases_flo1.allocate("fases_flo1");
  fases_flo2.allocate("fases_flo2");
  fases_flo3.allocate("fases_flo3");
 log_L50prior = log(L50prior);
 log_s1prior = log(s1prior);
 log_s2prior = log(s2prior);
  L50pel_prior.allocate("L50pel_prior");
  s1pel_prior.allocate("s1pel_prior");
  cv_parselp.allocate(1,2,"cv_parselp");
  fases_pel1.allocate("fases_pel1");
  fases_pel2.allocate("fases_pel2");
 log_L50pel_prior = log(L50pel_prior);
 log_s1pel_prior = log(s1pel_prior);
  L50recl_prior.allocate("L50recl_prior");
  s1recl_prior.allocate("s1recl_prior");
  cv_parselr.allocate(1,2,"cv_parselr");
  fases_recl1.allocate("fases_recl1");
  fases_recl2.allocate("fases_recl2");
 log_L50recl_prior = log(L50recl_prior);
 log_s1recl_prior = log(s1recl_prior);
  nbloques1.allocate("nbloques1");
  ybloques1.allocate(1,nbloques1,"ybloques1");
  nqbloques.allocate("nqbloques");
  yqbloques.allocate(1,nqbloques,"yqbloques");
  fases_qf.allocate("fases_qf");
  fases_qrecl.allocate("fases_qrecl");
  fases_qpel.allocate("fases_qpel");
  fases_F.allocate("fases_F");
  fases_devRt.allocate("fases_devRt");
  fases_devRs.allocate("fases_devRs");
  npbr.allocate("npbr");
  pbr.allocate(1,npbr,"pbr");
  ntime_sim.allocate("ntime_sim");
}

void model_parameters::initializationfunction(void)
{
  log_Linf.set_initial_value(log_Linf_prior);
  log_k.set_initial_value(log_k_prior);
  log_Lr.set_initial_value(log_Lr_prior);
  log_sr.set_initial_value(log_sr_prior);
  log_beta.set_initial_value(log_beta_prior);
  log_M.set_initial_value(log_M_prior);
  log_h.set_initial_value(log_h_prior);
  log_L50.set_initial_value(log_L50prior);
  log_sigma1.set_initial_value(log_s1prior);
  log_sigma2.set_initial_value(log_s2prior);
  log_L50pel.set_initial_value(log_L50pel_prior);
  log_spel.set_initial_value(log_s1pel_prior);
  log_L50recl.set_initial_value(log_L50recl_prior);
  log_srecl.set_initial_value(log_s1recl_prior);
  if (global_datafile)
  {
    delete global_datafile;
    global_datafile = NULL;
  }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  log_L50.allocate(1,nbloques1,fases_flo1,"log_L50");
  log_sigma1.allocate(1,nbloques1,fases_flo2,"log_sigma1");
  log_sigma2.allocate(1,nbloques1,fases_flo3,"log_sigma2");
  log_L50pel.allocate(fases_pel1,"log_L50pel");
  log_spel.allocate(fases_pel2,"log_spel");
  log_L50recl.allocate(fases_recl1,"log_L50recl");
  log_srecl.allocate(fases_recl2,"log_srecl");
  log_Rmed.allocate(1,"log_Rmed");
  log_desv_Rt.allocate(1,nyear,-10,10,fases_devRt,"log_desv_Rt");
  log_desv_Rs.allocate(1,nt,-10,10,fases_devRs,"log_desv_Rs");
  log_F.allocate(1,ntime,fases_F,"log_F");
  log_qflo.allocate(1,nqbloques,fases_qf,"log_qflo");
  log_q_pel.allocate(fases_qpel,"log_q_pel");
  log_q_recl.allocate(fases_qrecl,"log_q_recl");
  log_Linf.allocate(fase_Linf,"log_Linf");
  log_k.allocate(fase_k,"log_k");
  log_Lr.allocate(fase_Lr,"log_Lr");
  log_sr.allocate(fase_sr,"log_sr");
  log_beta.allocate(fase_beta,"log_beta");
  log_M.allocate(fase_M,"log_M");
  log_h.allocate(fase_h,"log_h");
  BMflo.allocate(1,ntime,"BMflo");
  #ifndef NO_AD_INITIALIZE
    BMflo.initialize();
  #endif
  Brec.allocate(1,ntime,"Brec");
  #ifndef NO_AD_INITIALIZE
    Brec.initialize();
  #endif
  pred_CPUE.allocate(1,ntime,"pred_CPUE");
  #ifndef NO_AD_INITIALIZE
    pred_CPUE.initialize();
  #endif
  pred_Desemb.allocate(1,ntime,"pred_Desemb");
  #ifndef NO_AD_INITIALIZE
    pred_Desemb.initialize();
  #endif
  pred_Brecl.allocate(1,ntime,"pred_Brecl");
  #ifndef NO_AD_INITIALIZE
    pred_Brecl.initialize();
  #endif
  pred_Bpel.allocate(1,ntime,"pred_Bpel");
  #ifndef NO_AD_INITIALIZE
    pred_Bpel.initialize();
  #endif
  likeval.allocate(1,15,"likeval");
  #ifndef NO_AD_INITIALIZE
    likeval.initialize();
  #endif
  Neq.allocate(1,ntallas,"Neq");
  #ifndef NO_AD_INITIALIZE
    Neq.initialize();
  #endif
  Rpred.allocate(1,ntime,"Rpred");
  #ifndef NO_AD_INITIALIZE
    Rpred.initialize();
  #endif
  Unos_year.allocate(1,ntime,"Unos_year");
  #ifndef NO_AD_INITIALIZE
    Unos_year.initialize();
  #endif
  Unos_tallas.allocate(1,ntallas,"Unos_tallas");
  #ifndef NO_AD_INITIALIZE
    Unos_tallas.initialize();
  #endif
  delta.allocate(1,ntallas,"delta");
  #ifndef NO_AD_INITIALIZE
    delta.initialize();
  #endif
  Lesp.allocate(1,ntallas,"Lesp");
  #ifndef NO_AD_INITIALIZE
    Lesp.initialize();
  #endif
  sigmaL.allocate(1,ntallas,"sigmaL");
  #ifndef NO_AD_INITIALIZE
    sigmaL.initialize();
  #endif
  pre.allocate(1,ntallas,"pre");
  #ifndef NO_AD_INITIALIZE
    pre.initialize();
  #endif
  fobs_marg.allocate(1,ntallas,"fobs_marg");
  #ifndef NO_AD_INITIALIZE
    fobs_marg.initialize();
  #endif
  fpred_marg.allocate(1,ntallas,"fpred_marg");
  #ifndef NO_AD_INITIALIZE
    fpred_marg.initialize();
  #endif
  BDo.allocate(1,ntime,"BDo");
  #ifndef NO_AD_INITIALIZE
    BDo.initialize();
  #endif
  yrs.allocate(1,ntime,"yrs");
  #ifndef NO_AD_INITIALIZE
    yrs.initialize();
  #endif
  Desemb.allocate(1,ntime,"Desemb");
  #ifndef NO_AD_INITIALIZE
    Desemb.initialize();
  #endif
  CPUE.allocate(1,ntime,"CPUE");
  #ifndef NO_AD_INITIALIZE
    CPUE.initialize();
  #endif
  Brecl.allocate(1,ntime,"Brecl");
  #ifndef NO_AD_INITIALIZE
    Brecl.initialize();
  #endif
  Bpel.allocate(1,ntime,"Bpel");
  #ifndef NO_AD_INITIALIZE
    Bpel.initialize();
  #endif
  Reclutas.allocate(1,ntime,"Reclutas");
  #ifndef NO_AD_INITIALIZE
    Reclutas.initialize();
  #endif
  penalty.allocate(1,5,"penalty");
  #ifndef NO_AD_INITIALIZE
    penalty.initialize();
  #endif
  Lmed_obs.allocate(1,N_ftc,"Lmed_obs");
  #ifndef NO_AD_INITIALIZE
    Lmed_obs.initialize();
  #endif
  Lmed_pred.allocate(1,N_ftc,"Lmed_pred");
  #ifndef NO_AD_INITIALIZE
    Lmed_pred.initialize();
  #endif
  Lmed_obs_r.allocate(1,N_ftr,"Lmed_obs_r");
  #ifndef NO_AD_INITIALIZE
    Lmed_obs_r.initialize();
  #endif
  Lmed_pred_r.allocate(1,N_ftr,"Lmed_pred_r");
  #ifndef NO_AD_INITIALIZE
    Lmed_pred_r.initialize();
  #endif
  Lmed_obs_p.allocate(1,N_ftp,"Lmed_obs_p");
  #ifndef NO_AD_INITIALIZE
    Lmed_obs_p.initialize();
  #endif
  Lmed_pred_p.allocate(1,N_ftp,"Lmed_pred_p");
  #ifndef NO_AD_INITIALIZE
    Lmed_pred_p.initialize();
  #endif
  BDo_din.allocate(1,ntime,"BDo_din");
  #ifndef NO_AD_INITIALIZE
    BDo_din.initialize();
  #endif
  pri.allocate(1,14,"pri");
  #ifndef NO_AD_INITIALIZE
    pri.initialize();
  #endif
  years_ftc.allocate(1,N_ftc,"years_ftc");
  #ifndef NO_AD_INITIALIZE
    years_ftc.initialize();
  #endif
  years_ftr.allocate(1,N_ftr,"years_ftr");
  #ifndef NO_AD_INITIALIZE
    years_ftr.initialize();
  #endif
  years_ftp.allocate(1,N_ftp,"years_ftp");
  #ifndef NO_AD_INITIALIZE
    years_ftp.initialize();
  #endif
  cv_index.allocate(1,4,1,ntime,"cv_index");
  #ifndef NO_AD_INITIALIZE
    cv_index.initialize();
  #endif
  nm.allocate(1,3,1,ntime,"nm");
  #ifndef NO_AD_INITIALIZE
    nm.initialize();
  #endif
  N0_din.allocate(1,ntime,1,ntallas,"N0_din");
  #ifndef NO_AD_INITIALIZE
    N0_din.initialize();
  #endif
  S1.allocate(1,nbloques1,1,ntallas,"S1");
  #ifndef NO_AD_INITIALIZE
    S1.initialize();
  #endif
  Sel.allocate(1,ntime,1,ntallas,"Sel");
  #ifndef NO_AD_INITIALIZE
    Sel.initialize();
  #endif
  Sel_recl.allocate(1,ntallas,"Sel_recl");
  #ifndef NO_AD_INITIALIZE
    Sel_recl.initialize();
  #endif
  Sel_pel.allocate(1,ntallas,"Sel_pel");
  #ifndef NO_AD_INITIALIZE
    Sel_pel.initialize();
  #endif
  F.allocate(1,ntime,1,ntallas,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  Z.allocate(1,ntime,1,ntallas,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  S.allocate(1,ntime,1,ntallas,"S");
  #ifndef NO_AD_INITIALIZE
    S.initialize();
  #endif
  N.allocate(1,ntime,1,ntallas,"N");
  #ifndef NO_AD_INITIALIZE
    N.initialize();
  #endif
  NM.allocate(1,ntime,1,ntallas,"NM");
  #ifndef NO_AD_INITIALIZE
    NM.initialize();
  #endif
  NMD.allocate(1,ntime,1,ntallas,"NMD");
  #ifndef NO_AD_INITIALIZE
    NMD.initialize();
  #endif
  NDv.allocate(1,ntime,1,ntallas,"NDv");
  #ifndef NO_AD_INITIALIZE
    NDv.initialize();
  #endif
  Nrec.allocate(1,ntime,1,ntallas,"Nrec");
  #ifndef NO_AD_INITIALIZE
    Nrec.initialize();
  #endif
  NVflo.allocate(1,ntime,1,ntallas,"NVflo");
  #ifndef NO_AD_INITIALIZE
    NVflo.initialize();
  #endif
  pred_Ctot.allocate(1,ntime,1,ntallas,"pred_Ctot");
  #ifndef NO_AD_INITIALIZE
    pred_Ctot.initialize();
  #endif
  pred_Nrecl.allocate(1,ntime,1,ntallas,"pred_Nrecl");
  #ifndef NO_AD_INITIALIZE
    pred_Nrecl.initialize();
  #endif
  pred_Npel.allocate(1,ntime,1,ntallas,"pred_Npel");
  #ifndef NO_AD_INITIALIZE
    pred_Npel.initialize();
  #endif
  No.allocate(1,ntallas,"No");
  #ifndef NO_AD_INITIALIZE
    No.initialize();
  #endif
  T.allocate(1,ntallas,1,ntallas,"T");
  #ifndef NO_AD_INITIALIZE
    T.initialize();
  #endif
  Diag.allocate(1,ntallas,1,ntallas,"Diag");
  #ifndef NO_AD_INITIALIZE
    Diag.initialize();
  #endif
  Id.allocate(1,ntallas,1,ntallas,"Id");
  #ifndef NO_AD_INITIALIZE
    Id.initialize();
  #endif
  NMDv.allocate(1,ntime,1,ntallas,"NMDv");
  #ifndef NO_AD_INITIALIZE
    NMDv.initialize();
  #endif
  pobs_recl.allocate(1,N_ftr,1,ntallas,"pobs_recl");
  #ifndef NO_AD_INITIALIZE
    pobs_recl.initialize();
  #endif
  ppred_recl.allocate(1,N_ftr,1,ntallas,"ppred_recl");
  #ifndef NO_AD_INITIALIZE
    ppred_recl.initialize();
  #endif
  pobs_pel.allocate(1,N_ftp,1,ntallas,"pobs_pel");
  #ifndef NO_AD_INITIALIZE
    pobs_pel.initialize();
  #endif
  ppred_pel.allocate(1,N_ftp,1,ntallas,"ppred_pel");
  #ifndef NO_AD_INITIALIZE
    ppred_pel.initialize();
  #endif
  pobs.allocate(1,N_ftc,1,ntallas,"pobs");
  #ifndef NO_AD_INITIALIZE
    pobs.initialize();
  #endif
  ppred.allocate(1,N_ftc,1,ntallas,"ppred");
  #ifndef NO_AD_INITIALIZE
    ppred.initialize();
  #endif
  suma1.allocate("suma1");
  #ifndef NO_AD_INITIALIZE
  suma1.initialize();
  #endif
  suma2.allocate("suma2");
  #ifndef NO_AD_INITIALIZE
  suma2.initialize();
  #endif
  suma3.allocate("suma3");
  #ifndef NO_AD_INITIALIZE
  suma3.initialize();
  #endif
  suma4.allocate("suma4");
  #ifndef NO_AD_INITIALIZE
  suma4.initialize();
  #endif
  suma5.allocate("suma5");
  #ifndef NO_AD_INITIALIZE
  suma5.initialize();
  #endif
  suma6.allocate("suma6");
  #ifndef NO_AD_INITIALIZE
  suma6.initialize();
  #endif
  suma7.allocate("suma7");
  #ifndef NO_AD_INITIALIZE
  suma7.initialize();
  #endif
  suma8.allocate("suma8");
  #ifndef NO_AD_INITIALIZE
  suma8.initialize();
  #endif
  suma9.allocate("suma9");
  #ifndef NO_AD_INITIALIZE
  suma9.initialize();
  #endif
  suma10.allocate("suma10");
  #ifndef NO_AD_INITIALIZE
  suma10.initialize();
  #endif
  So.allocate("So");
  #ifndef NO_AD_INITIALIZE
  So.initialize();
  #endif
  alfa.allocate("alfa");
  #ifndef NO_AD_INITIALIZE
  alfa.initialize();
  #endif
  beta.allocate("beta");
  #ifndef NO_AD_INITIALIZE
  beta.initialize();
  #endif
  Linf.allocate("Linf");
  #ifndef NO_AD_INITIALIZE
  Linf.initialize();
  #endif
  k.allocate("k");
  #ifndef NO_AD_INITIALIZE
  k.initialize();
  #endif
  M.allocate("M");
  #ifndef NO_AD_INITIALIZE
  M.initialize();
  #endif
  Lr.allocate("Lr");
  #ifndef NO_AD_INITIALIZE
  Lr.initialize();
  #endif
  sr.allocate("sr");
  #ifndef NO_AD_INITIALIZE
  sr.initialize();
  #endif
  h.allocate("h");
  #ifndef NO_AD_INITIALIZE
  h.initialize();
  #endif
  BDp.allocate("BDp");
  #ifndef NO_AD_INITIALIZE
  BDp.initialize();
  #endif
  Npplus.allocate("Npplus");
  #ifndef NO_AD_INITIALIZE
  Npplus.initialize();
  #endif
  Bp_anch.allocate("Bp_anch");
  #ifndef NO_AD_INITIALIZE
  Bp_anch.initialize();
  #endif
  nm1.allocate("nm1");
  #ifndef NO_AD_INITIALIZE
  nm1.initialize();
  #endif
  cuenta1.allocate("cuenta1");
  #ifndef NO_AD_INITIALIZE
  cuenta1.initialize();
  #endif
  alfa_sr.allocate("alfa_sr");
  #ifndef NO_AD_INITIALIZE
  alfa_sr.initialize();
  #endif
  beta_sr.allocate("beta_sr");
  #ifndef NO_AD_INITIALIZE
  beta_sr.initialize();
  #endif
  pF.allocate("pF");
  #ifndef NO_AD_INITIALIZE
  pF.initialize();
  #endif
  Np.allocate(1,ntallas,"Np");
  #ifndef NO_AD_INITIALIZE
    Np.initialize();
  #endif
  Zpbr.allocate(1,ntallas,"Zpbr");
  #ifndef NO_AD_INITIALIZE
    Zpbr.initialize();
  #endif
  Fpbr.allocate(1,ntallas,"Fpbr");
  #ifndef NO_AD_INITIALIZE
    Fpbr.initialize();
  #endif
  Sp.allocate(1,ntallas,"Sp");
  #ifndef NO_AD_INITIALIZE
    Sp.initialize();
  #endif
  BT.allocate(1,ntime,"BT");
  #ifndef NO_AD_INITIALIZE
    BT.initialize();
  #endif
  Bp.allocate(1,npbr,1,ntime_sim,"Bp");
  #ifndef NO_AD_INITIALIZE
    Bp.initialize();
  #endif
  CTPp.allocate(1,ntallas,"CTPp");
  #ifndef NO_AD_INITIALIZE
    CTPp.initialize();
  #endif
  Rpp.allocate(1,npbr,1,ntime_sim,"Rpp");
  #ifndef NO_AD_INITIALIZE
    Rpp.initialize();
  #endif
  SSBo.allocate("SSBo");
  #ifndef NO_AD_INITIALIZE
  SSBo.initialize();
  #endif
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  BD.allocate(1,ntime,"BD");
  Fmort.allocate(1,ntime,"Fmort");
  RPRlp.allocate(1,ntime,"RPRlp");
  RPRdin.allocate(1,ntime,"RPRdin");
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
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
}

void model_parameters::set_runtime(void)
{
  dvector temp("{1.e-1,1.e-01,1.e-03,1e-3,1e-5}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
  dvector temp1("{100,100,200,3000,3500}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
}

void model_parameters::userfunction(void)
{
  f =0.0;
 Eval_Trans_talla_talla();
 Eval_selectividad();
 Eval_mortalidades();
 Eval_abundancia();
 Eval_biomasas();
 Eval_capturas_predichas();
 Eval_indices();
 Eval_logverosim();
 Eval_funcion_objetivo();
}

void model_parameters::Eval_Trans_talla_talla(void)
{
 Linf=exp(log_Linf);
 k=exp(log_k);
 Lr=exp(log_Lr);
 sr=exp(log_sr);
 beta=exp(log_beta);
 M=exp(log_M);
 h=exp(log_h);
 int i, j;
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
}

void model_parameters::Eval_selectividad(void)
{
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
 Sel_recl=1/(1+exp(-log(19)*(Tallas-exp(log_L50recl))/exp(log_srecl)));
 Sel_pel=1/(1+exp(-log(19)*(Tallas-exp(log_L50pel))/exp(log_spel)));
}

void model_parameters::Eval_mortalidades(void)
{
  //<<3<<endl;
 Fmort=exp(log_F);
 for (int j=1;j<=ntime;j++){
  F(j)=Sel(j)*Fmort(j);
  }
 Z=F+M;
 S=mfexp(-1.0*Z);
}

void model_parameters::Eval_abundancia(void)
{
   //<<4<<endl;
 int i, j;
  Lr=Par_bio(3);
  sr=Par_bio(4);
  if (active(log_Lr)){Lr=mfexp(log_Lr);}
  if (active(log_sr)){sr=mfexp(log_sr);}
  pre=exp(-0.5*square((Tallas-Lr)/sr));
  pre/=sum(pre);
  for (int j=1;j<=ntallas;j++){
  Id(j,j)=1;}
  No=pre*exp(log_Rmed)*inv(Id-(exp(-M)*Id)*T);
  SSBo=sum(elem_prod(No*mfexp(-dt(1)*M),elem_prod(Wmed,msex)));
  alfa_sr=4*h*exp(log_Rmed)/(5*h-1);//
  beta_sr=(1-h)*SSBo/(5*h-1);// Reclutamiento
  Reclutas(1)=mfexp(log_Rmed+log_desv_Rt(1)+log_desv_Rs(1));
  Rpred(1)=mfexp(log_Rmed);
  for (int j=1;j<=ntallas;j++){
  Diag(j,j)=exp(-Z(1,j));}
  Neq=pre*exp(log_Rmed)*inv(Id-Diag*T);
  N(1)=Neq;
  NMD(1)=elem_prod(elem_prod(N(1),mfexp(-dt(1)*Z(1))),msex);
  BD(1)=sum(elem_prod(Wmed,NMD(1)));
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
}

void model_parameters::Eval_biomasas(void)
{
  //<<5<<endl;
 NVflo=elem_prod(elem_prod(N,mfexp(-dt(2)*(Z))),Sel);
 BMflo=Wmed*trans(NVflo);
 BT=Wmed*trans(N);
 //BDo=sum(elem_prod(No(nedades),Wmed));
 RPRlp=BD/SSBo;
 RPRdin=elem_div(BD,BDo_din);
}

void model_parameters::Eval_capturas_predichas(void)
{
  //<<6<<endl;
 pred_Ctot=elem_prod(elem_div(F,Z),elem_prod(1.-S,N));
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
}

void model_parameters::Eval_indices(void)
{
   for (int i=1;i<=ntime;i++){
      for (int j=1;j<=nqbloques;j++){
              if (yrs(i)>=yqbloques(j)){
                 pred_CPUE(i)=exp(log_qflo(j))*BMflo(i);}
       }
   }
}

void model_parameters::Eval_logverosim(void)
{
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
}

void model_parameters::Eval_funcion_objetivo(void)
{
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
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
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
}

void model_parameters::final_calcs()
{
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
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
 time(&start);
 arrmblsize = 90000000; 
 gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7); 
 gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7); 
 gradient_structure::set_MAX_NVAR_OFFSET(5000); 
 gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000); 
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
  auto start = std::chrono::high_resolution_clock::now();
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  std::cout << endl << argv[0] << " elapsed time is " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << " microseconds." << endl;
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
