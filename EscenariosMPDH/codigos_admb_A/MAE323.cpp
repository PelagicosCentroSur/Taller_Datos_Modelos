#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
  #include  <admodel.h>
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
#include <MAE323.htp>

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
  nanos.allocate("nanos");
  nedades.allocate("nedades");
  ntallas.allocate("ntallas");
  edades.allocate(1,nedades,"edades");
  Tallas.allocate(1,ntallas,"Tallas");
  msex.allocate(1,nedades,"msex");
  matdat.allocate(1,nanos,1,9,"matdat");
  Ctot.allocate(1,nanos,1,nedades,"Ctot");
  Ccru_a.allocate(1,nanos,1,nedades,"Ccru_a");
  Ccru_pel.allocate(1,nanos,1,nedades,"Ccru_pel");
  Ccru_l.allocate(1,nanos,1,ntallas,"Ccru_l");
  Wmed.allocate(1,nanos,1,nedades,"Wmed");
  Wmedp_3.allocate(1,5,"Wmedp_3");
  Win.allocate(1,nanos,1,nedades,"Win");
  Winip_3.allocate(1,5,"Winip_3");
  error_edad.allocate(1,nedades,1,nedades,"error_edad");
  sigmaR.allocate("sigmaR");
  cvpriorq_reclas.allocate("cvpriorq_reclas");
  cvpriorq_pelaces.allocate("cvpriorq_pelaces");
  log_priorRo.allocate("log_priorRo");
  nmus.allocate(1,4,"nmus");
  dt.allocate(1,4,"dt");
  opt1_fase.allocate("opt1_fase");
  opt2_fase.allocate("opt2_fase");
  opt_Scru1.allocate("opt_Scru1");
  opt_Scru2.allocate("opt_Scru2");
  lam.allocate(1,4,"lam");
  opt_qrecl.allocate("opt_qrecl");
  opt_qpela.allocate("opt_qpela");
  opt_qmph.allocate("opt_qmph");
  pars_Bio.allocate(1,5,"pars_Bio");
  opt_Lo.allocate("opt_Lo");
  opt_cv.allocate("opt_cv");
  erredad.allocate("erredad");
  opt_M.allocate("opt_M");
  opt_Ro.allocate("opt_Ro");
  opt_devR.allocate("opt_devR");
  opt_devNo.allocate("opt_devNo");
  opt_F.allocate("opt_F");
  opt_Fspr.allocate("opt_Fspr");
  npbr.allocate("npbr");
  ratio.allocate(1,npbr,"ratio");
  nproy.allocate("nproy");
  opProy.allocate("opProy");
  oprec.allocate("oprec");
  prop.allocate(1,2,"prop");
  mF.allocate("mF");
  Desembactual.allocate("Desembactual");
  prop_est.allocate(1,5,"prop_est");
  Fmedian1.allocate("Fmedian1");
  FrmsH2.allocate("FrmsH2");
  Desemb2doSem.allocate("Desemb2doSem");
  Descarte.allocate("Descarte");
  Remanente.allocate("Remanente");
}

void model_parameters::initializationfunction(void)
{
  log_Ro.set_initial_value(12.54);
  log_Lo.set_initial_value(2);
  log_cv_edad.set_initial_value(-2.52);
  log_M.set_initial_value(0);
  log_qrecl.set_initial_value(0);
  log_qpela.set_initial_value(0);
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
  A50f.allocate(-1,2,opt1_fase,"A50f");
  log_rangof.allocate(-4,0,opt1_fase,"log_rangof");
  A50c.allocate(-1,2,opt_Scru1,"A50c");
  log_rangoc.allocate(-4,0.6,opt_Scru1,"log_rangoc");
  A50pela.allocate(-1,2,opt1_fase,"A50pela");
  log_rangopela.allocate(-4,0.6,opt1_fase,"log_rangopela");
  log_Sflo.allocate(1,nedades,-11,0,opt2_fase,"log_Sflo");
  log_Scru.allocate(1,nedades,-11,0,opt_Scru2,"log_Scru");
  log_Scru_pela.allocate(1,nedades,-11,0,opt_Scru2,"log_Scru_pela");
  log_Ro.allocate(5,20,opt_Ro,"log_Ro");
  log_desv_No.allocate(1,nedades-1,-10,10,opt_devNo,"log_desv_No");
  log_desv_Rt.allocate(1,nanos,-10,10,opt_devR,"log_desv_Rt");
  log_Ft.allocate(1,nanos,-6,1.6,opt_F,"log_Ft");
  log_qrecl.allocate(opt_qrecl,"log_qrecl");
  log_qpela.allocate(opt_qpela,"log_qpela");
  log_qmph.allocate(opt_qmph,"log_qmph");
  log_Lo.allocate(1,2.4,opt_Lo,"log_Lo");
  log_cv_edad.allocate(-4,.69,opt_cv,"log_cv_edad");
  log_M.allocate(-0.3,0.4,opt_M,"log_M");
  log_Fref.allocate(1,npbr,0.01,2.,opt_Fspr,"log_Fref");
  anos.allocate(1,nanos,"anos");
  #ifndef NO_AD_INITIALIZE
    anos.initialize();
  #endif
  Unos_edad.allocate(1,nedades,"Unos_edad");
  #ifndef NO_AD_INITIALIZE
    Unos_edad.initialize();
  #endif
  Unos_anos.allocate(1,nanos,"Unos_anos");
  #ifndef NO_AD_INITIALIZE
    Unos_anos.initialize();
  #endif
  Unos_tallas.allocate(1,ntallas,"Unos_tallas");
  #ifndef NO_AD_INITIALIZE
    Unos_tallas.initialize();
  #endif
  Sel_f.allocate(1,nanos,1,nedades,"Sel_f");
  #ifndef NO_AD_INITIALIZE
    Sel_f.initialize();
  #endif
  Scru.allocate(1,nanos,1,nedades,"Scru");
  #ifndef NO_AD_INITIALIZE
    Scru.initialize();
  #endif
  Scru_pela.allocate(1,nanos,1,nedades,"Scru_pela");
  #ifndef NO_AD_INITIALIZE
    Scru_pela.initialize();
  #endif
  pSf.allocate(1,nedades-2,"pSf");
  #ifndef NO_AD_INITIALIZE
    pSf.initialize();
  #endif
  pSc.allocate(1,nedades-2,"pSc");
  #ifndef NO_AD_INITIALIZE
    pSc.initialize();
  #endif
  pStotf.allocate(1,nanos,"pStotf");
  #ifndef NO_AD_INITIALIZE
    pStotf.initialize();
  #endif
  pStotc.allocate(1,nanos,"pStotc");
  #ifndef NO_AD_INITIALIZE
    pStotc.initialize();
  #endif
  Linf.allocate("Linf");
  #ifndef NO_AD_INITIALIZE
  Linf.initialize();
  #endif
  k.allocate("k");
  #ifndef NO_AD_INITIALIZE
  k.initialize();
  #endif
  Lo.allocate("Lo");
  #ifndef NO_AD_INITIALIZE
  Lo.initialize();
  #endif
  cv_edad.allocate("cv_edad");
  #ifndef NO_AD_INITIALIZE
  cv_edad.initialize();
  #endif
  mu_edad.allocate(1,nedades,"mu_edad");
  #ifndef NO_AD_INITIALIZE
    mu_edad.initialize();
  #endif
  sigma_edad.allocate(1,nedades,"sigma_edad");
  #ifndef NO_AD_INITIALIZE
    sigma_edad.initialize();
  #endif
  Prob_talla.allocate(1,nedades,1,ntallas,"Prob_talla");
  #ifndef NO_AD_INITIALIZE
    Prob_talla.initialize();
  #endif
  P1.allocate(1,nedades,1,ntallas,"P1");
  #ifndef NO_AD_INITIALIZE
    P1.initialize();
  #endif
  P2.allocate(1,nedades,1,ntallas,"P2");
  #ifndef NO_AD_INITIALIZE
    P2.initialize();
  #endif
  P3.allocate(1,nedades,1,ntallas,"P3");
  #ifndef NO_AD_INITIALIZE
    P3.initialize();
  #endif
  M.allocate("M");
  #ifndef NO_AD_INITIALIZE
  M.initialize();
  #endif
  Ftot.allocate(1,nanos,1,nedades,"Ftot");
  #ifndef NO_AD_INITIALIZE
    Ftot.initialize();
  #endif
  Z.allocate(1,nanos,1,nedades,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  S.allocate(1,nanos,1,nedades,"S");
  #ifndef NO_AD_INITIALIZE
    S.initialize();
  #endif
  Neq.allocate(1,nedades,"Neq");
  #ifndef NO_AD_INITIALIZE
    Neq.initialize();
  #endif
  SSBo.allocate("SSBo");
  #ifndef NO_AD_INITIALIZE
  SSBo.initialize();
  #endif
  N.allocate(1,nanos,1,nedades,"N");
  #ifndef NO_AD_INITIALIZE
    N.initialize();
  #endif
  NM.allocate(1,nanos,1,nedades,"NM");
  #ifndef NO_AD_INITIALIZE
    NM.initialize();
  #endif
  NVflo.allocate(1,nanos,1,nedades,"NVflo");
  #ifndef NO_AD_INITIALIZE
    NVflo.initialize();
  #endif
  NVcru.allocate(1,nanos,1,nedades,"NVcru");
  #ifndef NO_AD_INITIALIZE
    NVcru.initialize();
  #endif
  NVpel.allocate(1,nanos,1,nedades,"NVpel");
  #ifndef NO_AD_INITIALIZE
    NVpel.initialize();
  #endif
  NVpel_l.allocate(1,nanos,1,ntallas,"NVpel_l");
  #ifndef NO_AD_INITIALIZE
    NVpel_l.initialize();
  #endif
  NMD.allocate(1,nanos,1,nedades,"NMD");
  #ifndef NO_AD_INITIALIZE
    NMD.initialize();
  #endif
  Reclutas.allocate(1,nanos,"Reclutas");
  SSB.allocate(1,nanos,"SSB");
  BD.allocate(1,nanos,"BD");
  #ifndef NO_AD_INITIALIZE
    BD.initialize();
  #endif
  BT.allocate(1,nanos,"BT");
  BMflo.allocate(1,nanos,"BMflo");
  #ifndef NO_AD_INITIALIZE
    BMflo.initialize();
  #endif
  BMpel.allocate(1,nanos,"BMpel");
  #ifndef NO_AD_INITIALIZE
    BMpel.initialize();
  #endif
  Bcru.allocate(1,nanos,"Bcru");
  #ifndef NO_AD_INITIALIZE
    Bcru.initialize();
  #endif
  pred_Ctot.allocate(1,nanos,1,nedades,"pred_Ctot");
  #ifndef NO_AD_INITIALIZE
    pred_Ctot.initialize();
  #endif
  pobs_f.allocate(1,nanos,1,nedades,"pobs_f");
  #ifndef NO_AD_INITIALIZE
    pobs_f.initialize();
  #endif
  ppred_f.allocate(1,nanos,1,nedades,"ppred_f");
  #ifndef NO_AD_INITIALIZE
    ppred_f.initialize();
  #endif
  pobs_crua.allocate(1,nanos,1,nedades,"pobs_crua");
  #ifndef NO_AD_INITIALIZE
    pobs_crua.initialize();
  #endif
  ppred_crua.allocate(1,nanos,1,nedades,"ppred_crua");
  #ifndef NO_AD_INITIALIZE
    ppred_crua.initialize();
  #endif
  pobs_pel.allocate(1,nanos,1,nedades,"pobs_pel");
  #ifndef NO_AD_INITIALIZE
    pobs_pel.initialize();
  #endif
  ppred_pel.allocate(1,nanos,1,nedades,"ppred_pel");
  #ifndef NO_AD_INITIALIZE
    ppred_pel.initialize();
  #endif
  pobs_crul.allocate(1,nanos,1,ntallas,"pobs_crul");
  #ifndef NO_AD_INITIALIZE
    pobs_crul.initialize();
  #endif
  ppred_crul.allocate(1,nanos,1,ntallas,"ppred_crul");
  #ifndef NO_AD_INITIALIZE
    ppred_crul.initialize();
  #endif
  qrecl.allocate("qrecl");
  #ifndef NO_AD_INITIALIZE
  qrecl.initialize();
  #endif
  qpela.allocate("qpela");
  #ifndef NO_AD_INITIALIZE
  qpela.initialize();
  #endif
  Reclas.allocate(1,nanos,"Reclas");
  #ifndef NO_AD_INITIALIZE
    Reclas.initialize();
  #endif
  Reclas_pred.allocate(1,nanos,"Reclas_pred");
  #ifndef NO_AD_INITIALIZE
    Reclas_pred.initialize();
  #endif
  Pelaces.allocate(1,nanos,"Pelaces");
  #ifndef NO_AD_INITIALIZE
    Pelaces.initialize();
  #endif
  Pelaces_pred.allocate(1,nanos,"Pelaces_pred");
  #ifndef NO_AD_INITIALIZE
    Pelaces_pred.initialize();
  #endif
  MPH.allocate(1,nanos,"MPH");
  #ifndef NO_AD_INITIALIZE
    MPH.initialize();
  #endif
  MPH_pred.allocate(1,nanos,"MPH_pred");
  #ifndef NO_AD_INITIALIZE
    MPH_pred.initialize();
  #endif
  Desemb.allocate(1,nanos,"Desemb");
  #ifndef NO_AD_INITIALIZE
    Desemb.initialize();
  #endif
  Desemb_pred.allocate(1,nanos,"Desemb_pred");
  #ifndef NO_AD_INITIALIZE
    Desemb_pred.initialize();
  #endif
  Nv.allocate(1,nanos,1,nedades,"Nv");
  #ifndef NO_AD_INITIALIZE
    Nv.initialize();
  #endif
  NDv.allocate(1,nanos,1,nedades,"NDv");
  #ifndef NO_AD_INITIALIZE
    NDv.initialize();
  #endif
  NMDv.allocate(1,nanos,1,nedades,"NMDv");
  #ifndef NO_AD_INITIALIZE
    NMDv.initialize();
  #endif
  BDo.allocate(1,nanos,"BDo");
  #ifndef NO_AD_INITIALIZE
    BDo.initialize();
  #endif
  RPRdin.allocate(1,nanos,"RPRdin");
  #ifndef NO_AD_INITIALIZE
    RPRdin.initialize();
  #endif
  RPRequ.allocate(1,nanos,"RPRequ");
  #ifndef NO_AD_INITIALIZE
    RPRequ.initialize();
  #endif
  RPRequ2.allocate(1,nanos,"RPRequ2");
  #ifndef NO_AD_INITIALIZE
    RPRequ2.initialize();
  #endif
  RPRequ3.allocate(1,nanos,"RPRequ3");
  Frpr.allocate(1,nanos,"Frpr");
  Fspr.allocate(1,nedades,"Fspr");
  #ifndef NO_AD_INITIALIZE
    Fspr.initialize();
  #endif
  Zspr.allocate(1,nedades,"Zspr");
  #ifndef NO_AD_INITIALIZE
    Zspr.initialize();
  #endif
  Nspro.allocate(1,nedades,"Nspro");
  #ifndef NO_AD_INITIALIZE
    Nspro.initialize();
  #endif
  Nspr.allocate(1,nedades,"Nspr");
  #ifndef NO_AD_INITIALIZE
    Nspr.initialize();
  #endif
  Nmed.allocate(1,nedades,"Nmed");
  #ifndef NO_AD_INITIALIZE
    Nmed.initialize();
  #endif
  Bspro.allocate("Bspro");
  #ifndef NO_AD_INITIALIZE
  Bspro.initialize();
  #endif
  Bspr.allocate("Bspr");
  #ifndef NO_AD_INITIALIZE
  Bspr.initialize();
  #endif
  ratio_spr.allocate(1,npbr,"ratio_spr");
  #ifndef NO_AD_INITIALIZE
    ratio_spr.initialize();
  #endif
  Fmedian.allocate("Fmedian");
  #ifndef NO_AD_INITIALIZE
  Fmedian.initialize();
  #endif
  Fmed.allocate(1,nedades,"Fmed");
  #ifndef NO_AD_INITIALIZE
    Fmed.initialize();
  #endif
  Zmed.allocate(1,nedades,"Zmed");
  #ifndef NO_AD_INITIALIZE
    Zmed.initialize();
  #endif
  Bsprmed.allocate("Bsprmed");
  #ifndef NO_AD_INITIALIZE
  Bsprmed.initialize();
  #endif
  ratio_Fmed.allocate("ratio_Fmed");
  #ifndef NO_AD_INITIALIZE
  ratio_Fmed.initialize();
  #endif
  Bmed.allocate("Bmed");
  #ifndef NO_AD_INITIALIZE
  Bmed.initialize();
  #endif
  Bo.allocate("Bo");
  #ifndef NO_AD_INITIALIZE
  Bo.initialize();
  #endif
  Brms.allocate(1,npbr,"Brms");
  Blim.allocate("Blim");
  #ifndef NO_AD_INITIALIZE
  Blim.initialize();
  #endif
  cvar.allocate(1,4,1,nanos,"cvar");
  #ifndef NO_AD_INITIALIZE
    cvar.initialize();
  #endif
  likeval.allocate(1,15,"likeval");
  #ifndef NO_AD_INITIALIZE
    likeval.initialize();
  #endif
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  Fref_r0.allocate("Fref_r0");
  #ifndef NO_AD_INITIALIZE
  Fref_r0.initialize();
  #endif
  Frms_r0.allocate(1,nedades,"Frms_r0");
  #ifndef NO_AD_INITIALIZE
    Frms_r0.initialize();
  #endif
  Zrms_r0.allocate(1,nedades,"Zrms_r0");
  #ifndef NO_AD_INITIALIZE
    Zrms_r0.initialize();
  #endif
  CTP_r0.allocate(1,nedades,"CTP_r0");
  #ifndef NO_AD_INITIALIZE
    CTP_r0.initialize();
  #endif
  YTP_r0W.allocate(1,nedades,"YTP_r0W");
  #ifndef NO_AD_INITIALIZE
    YTP_r0W.initialize();
  #endif
  YTP_r0.allocate("YTP_r0");
  YTP_r01ersem.allocate("YTP_r01ersem");
  YTP_r01ersemR.allocate("YTP_r01ersemR");
  YTP_r0d.allocate("YTP_r0d");
  YTP_r0d1ersem.allocate("YTP_r0d1ersem");
  YTP_r0d1ersemR.allocate("YTP_r0d1ersemR");
  NMD_r0.allocate(1,nedades,"NMD_r0");
  #ifndef NO_AD_INITIALIZE
    NMD_r0.initialize();
  #endif
  BD_r0.allocate("BD_r0");
  #ifndef NO_AD_INITIALIZE
  BD_r0.initialize();
  #endif
  RPR_r0.allocate("RPR_r0");
  #ifndef NO_AD_INITIALIZE
  RPR_r0.initialize();
  #endif
  Np.allocate(1,nedades,"Np");
  #ifndef NO_AD_INITIALIZE
    Np.initialize();
  #endif
  Sp.allocate(1,nedades,"Sp");
  #ifndef NO_AD_INITIALIZE
    Sp.initialize();
  #endif
  Nvp.allocate(1,nedades,"Nvp");
  #ifndef NO_AD_INITIALIZE
    Nvp.initialize();
  #endif
  RPRp.allocate(1,nproy,"RPRp");
  #ifndef NO_AD_INITIALIZE
    RPRp.initialize();
  #endif
  Npp.allocate(1,nedades,"Npp");
  #ifndef NO_AD_INITIALIZE
    Npp.initialize();
  #endif
  Wmedp.allocate(1,nedades,"Wmedp");
  #ifndef NO_AD_INITIALIZE
    Wmedp.initialize();
  #endif
  Winp.allocate(1,nedades,"Winp");
  #ifndef NO_AD_INITIALIZE
    Winp.initialize();
  #endif
  Fref_p0.allocate("Fref_p0");
  #ifndef NO_AD_INITIALIZE
  Fref_p0.initialize();
  #endif
  Frms_p0.allocate(1,nedades,"Frms_p0");
  #ifndef NO_AD_INITIALIZE
    Frms_p0.initialize();
  #endif
  Zrms_p0.allocate(1,nedades,"Zrms_p0");
  #ifndef NO_AD_INITIALIZE
    Zrms_p0.initialize();
  #endif
  CTP_p0.allocate(1,nproy,1,nedades,"CTP_p0");
  #ifndef NO_AD_INITIALIZE
    CTP_p0.initialize();
  #endif
  YTP_p0W.allocate(1,nproy,1,nedades,"YTP_p0W");
  #ifndef NO_AD_INITIALIZE
    YTP_p0W.initialize();
  #endif
  YTP_p0.allocate(1,nproy,"YTP_p0");
  YTP_p0d.allocate(1,nproy,"YTP_p0d");
  BD_p0.allocate(1,nproy,"BD_p0");
  RPR_p0.allocate(1,nproy,"RPR_p0");
  NVrecl_p0.allocate(1,nedades,"NVrecl_p0");
  #ifndef NO_AD_INITIALIZE
    NVrecl_p0.initialize();
  #endif
  NVpel_p0.allocate(1,nedades,"NVpel_p0");
  #ifndef NO_AD_INITIALIZE
    NVpel_p0.initialize();
  #endif
  Brecl_p0.allocate(1,nproy,"Brecl_p0");
  #ifndef NO_AD_INITIALIZE
    Brecl_p0.initialize();
  #endif
  Bpel_p0.allocate(1,nproy,"Bpel_p0");
  #ifndef NO_AD_INITIALIZE
    Bpel_p0.initialize();
  #endif
  CBA_c0.allocate("CBA_c0");
  CBA_c0D.allocate("CBA_c0D");
  CBA_c0R.allocate("CBA_c0R");
  CBA_c0d.allocate("CBA_c0d");
  CBA_c0dD.allocate("CBA_c0dD");
  CBA_c0dR.allocate("CBA_c0dR");
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
  nm1.allocate("nm1");
  #ifndef NO_AD_INITIALIZE
  nm1.initialize();
  #endif
  nm2.allocate("nm2");
  #ifndef NO_AD_INITIALIZE
  nm2.initialize();
  #endif
  nm3.allocate("nm3");
  #ifndef NO_AD_INITIALIZE
  nm3.initialize();
  #endif
  nm4.allocate("nm4");
  #ifndef NO_AD_INITIALIZE
  nm4.initialize();
  #endif
  alfa.allocate("alfa");
  #ifndef NO_AD_INITIALIZE
  alfa.initialize();
  #endif
  beta.allocate("beta");
  #ifndef NO_AD_INITIALIZE
  beta.initialize();
  #endif
  cuenta1.allocate("cuenta1");
  #ifndef NO_AD_INITIALIZE
  cuenta1.initialize();
  #endif
  cuenta2.allocate("cuenta2");
  #ifndef NO_AD_INITIALIZE
  cuenta2.initialize();
  #endif
  cuenta3.allocate("cuenta3");
  #ifndef NO_AD_INITIALIZE
  cuenta3.initialize();
  #endif
  cuenta4.allocate("cuenta4");
  #ifndef NO_AD_INITIALIZE
  cuenta4.initialize();
  #endif
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
  anos=column(matdat,1);// asigno la 1 columna de indices a "anos"
  Reclas=column(matdat,2);
  cvar(1)=column(matdat,3);
  Pelaces=column(matdat,4);
  cvar(2)=column(matdat,5);
  MPH=column(matdat,6);
  cvar(3)=column(matdat,7);
  Desemb=column(matdat,8);
  cvar(4)=column(matdat,9);
  Unos_edad=1;;// lo uso en  operaciones matriciales con la edad
  Unos_anos=1;// lo uso en operaciones matriciales con el a?o
  Unos_tallas=1;// lo uso en operaciones matriciales con el a?o
  reporte_mcmc=0;
}

void model_parameters::set_runtime(void)
{
  dvector temp1("{200,1000,5000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
  dvector temp("{1e-3,1e-5,1e-6}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
}

void model_parameters::userfunction(void)
{
  f =0.0;
  if(opt1_fase>0)
  {
      Eval_selectividad_logis();
  }
  if(opt2_fase>0)
  {
      Eval_selectividad_nopar();
  }
  Eval_prob_talla_edad();
  Eval_mortalidades();
  Eval_abundancia();
  Eval_biomasas();
  Eval_capturas_predichas();
  Eval_indices();
  Eval_PBR();
  Eval_Estatus();
  Eval_logverosim();
  Eval_funcion_objetivo();
  Eval_CTP();
  Eval_mcmc();
}

void model_parameters::Eval_selectividad_logis(void)
{
   Sel_f     = outer_prod(Unos_anos,(elem_div(Unos_edad,(1+exp(-1.0*log(19)*(edades-A50f)/exp(log_rangof))))));
   Scru_pela = outer_prod(Unos_anos,(elem_div(Unos_edad,(1+exp(-1.0*log(19)*(edades-A50pela)/exp(log_rangopela))))));
   if (opt_Scru1>0)// evaluo si el indice es >0
    {
       Scru  = outer_prod(Unos_anos,(elem_div(Unos_edad,(1+exp(-1.0*log(19)*(edades-A50c)/exp(log_rangoc))))));
    }    
   else
    {
       Scru  = 1.;
    }
}

void model_parameters::Eval_selectividad_nopar(void)
{
  int i;
  for (i=1;i<=nanos;i++)
   {
          Sel_f(i)(1,nedades)=mfexp(log_Sflo);
      if (opt_Scru2>0)
      {
          Scru(i)(1,nedades)=mfexp(log_Scru);
          Scru_pela(i)(1,nedades)=mfexp(log_Scru_pela);
      }
      else
      {
          Scru=1.;
      }
      if (Sel_f(i)(nedades)>Sel_f(i)(nedades-1))
       {
           Sel_f(i)(nedades)=Sel_f(i)(nedades-1);
       }
      if (Scru(i)(nedades)>Scru(i)(nedades-1))
       {
           Scru(i)(nedades)=Scru(i)(nedades-1);
       }
   }
}

void model_parameters::Eval_prob_talla_edad(void)
{
  Linf    = pars_Bio(1);
  k       = pars_Bio(2);
  Lo      = pars_Bio(3);
  cv_edad = pars_Bio(4);
  if(active(log_Lo))
  {
   Lo=exp(log_Lo);
  }
  			if(active(log_cv_edad))
  				{
     				cv_edad=exp(log_cv_edad);
  				}
  int i, j;
  mu_edad(1)=Lo;
  for (i=2;i<=nedades;i++)
   {
      mu_edad(i)=Linf*(1-exp(-k))+exp(-k)*mu_edad(i-1);
   }
      sigma_edad=cv_edad*mu_edad;
  for (i=1;i<=nedades;i++)
   {
       P1(i)=(Tallas-mu_edad(i))/sigma_edad(i); //  standarizated deviation of the length (respect to expected length-at-age)
     for (j=1;j<=ntallas;j++) // cumulative pdf (normal distribution)
     {
       P2(i,j)=cumd_norm(P1(i,j));
     }
   } 
  for (i=1;i<=nedades;i++)// estimation of probabilities
  {
     for (j=2;j<=ntallas;j++)
     {
       P3(i,j)=P2(i,j)-P2(i,j-1);
     }
  } 
 /*
  P1=elem_div(outer_prod(Unos_edad,Unos_tallas),sqrt(2*3.1416)*outer_prod(sigma_edad,Unos_tallas));
  P2=mfexp(elem_div(-square(outer_prod(Unos_edad,Tallas)-outer_prod(mu_edad,Unos_tallas)),2*square(outer_prod(sigma_edad,Unos_tallas))));
  P3=elem_prod(P1,P2);
 */
  Prob_talla = elem_div(P3+1e-16,outer_prod(rowsum(P3+1e-16),Unos_tallas));// normalizo para que la suma sobre las edades sea 1.0
}

void model_parameters::Eval_mortalidades(void)
{
  if(opt_M>0)
  {
  M = exp(log_M);
  }
  else
  {
  M = pars_Bio(5);
  }
  Ftot = elem_prod(Sel_f,outer_prod(mfexp(log_Ft),Unos_edad));
  Z    = Ftot+M;
  S    = mfexp(-1.0*Z);
}

void model_parameters::Eval_abundancia(void)
{
 int i, j;
  if(opt_Ro<0)
  {
   log_Ro=log_priorRo;
  }
  for (i=1;i<=nanos;i++)
  {
    N(i,1) = mfexp(log_Ro+log_desv_Rt(i)+0.5*square(sigmaR)); 
  }
    Neq(1) = mfexp(log_Ro+0.5*square(sigmaR));
  for (i=2;i<=nedades;i++)
  { 
    Neq(i) = Neq(i-1)*exp(-1*M);
  }
    Neq(nedades)=Neq(nedades)/(1-exp(-1*M));
    SSBo=sum(elem_prod(Neq*exp(-dt(4)*M),elem_prod(msex,colsum(Win)/nanos)));
  for (i=2;i<=nedades;i++)
  {
  N(1)(i)=Neq(i)*exp(log_desv_No(i-1)+0.5*square(sigmaR));
  }
  for (i=2;i<=nanos;i++)
  {
      N(i)(2,nedades)=++elem_prod(N(i-1)(1,nedades-1),S(i-1)(1,nedades-1));
      N(i,nedades)+=N(i-1,nedades)*S(i-1,nedades); 
  }
}

void model_parameters::Eval_biomasas(void)
{
  NVcru    = elem_prod(elem_prod(N,mfexp(-dt(1)*Z)),Scru);// Crucero Reclas
  NVpel    = elem_prod(elem_prod(N,mfexp(-dt(2)*Z)),Scru_pela);// Pelaces
  if(erredad==1)
  {
  NVcru    = NVcru*error_edad;
  NVpel    = NVpel*error_edad;
  }
  NMD      = elem_prod(elem_prod(N,mfexp(-dt(3)*Z)),outer_prod(Unos_anos,msex));// desovante y MPH
  NVflo    = elem_prod(elem_prod(N,mfexp(-dt(4)*Z)),Sel_f);// explotable
  Reclutas = column(N,1);
  BD       = rowsum(elem_prod(NMD,Win));      // Desovante
  BT       = rowsum(elem_prod(N,Win));        // Total inicios de a?o biol
  BMflo    = rowsum(elem_prod(NVflo,Win));    // Biomasa explotable
  BMpel    = rowsum(elem_prod(NVpel,Win));    // pelaces
  Bcru     = rowsum(elem_prod(NVcru,Wmed));   // Reclas, mitad a?o biol
}

void model_parameters::Eval_capturas_predichas(void)
{
  pred_Ctot    = (elem_prod(elem_div(Ftot,Z),elem_prod(1.-S,N)));
  if(erredad==1)
  {
  pred_Ctot    = pred_Ctot*error_edad;
  }
  Desemb_pred  = rowsum(elem_prod(pred_Ctot,Wmed));
  pobs_f       = elem_div(Ctot,outer_prod(rowsum(Ctot+1e-10),Unos_edad));
  ppred_f      = elem_div(pred_Ctot,outer_prod(rowsum(pred_Ctot),Unos_edad));
  pobs_crua    = elem_div(Ccru_a,outer_prod(rowsum(Ccru_a+1e-10),Unos_edad));
  ppred_crua   = elem_div(NVcru,outer_prod(rowsum(NVcru),Unos_edad));
  pobs_pel     = elem_div(Ccru_pel,outer_prod(rowsum(Ccru_pel+1e-10),Unos_edad));
  ppred_pel    = elem_div(NVpel,outer_prod(rowsum(NVpel),Unos_edad));
  pobs_crul    = elem_div(Ccru_l,outer_prod(rowsum(Ccru_l+1e-10),Unos_tallas));
  ppred_crul   = elem_div(NVpel*Prob_talla,outer_prod(rowsum(NVpel),Unos_tallas));
}

void model_parameters::Eval_indices(void)
{
 Reclas_pred   = exp(log_qrecl)*Bcru;
 Pelaces_pred  = exp(log_qpela)*BMpel;
 MPH_pred      = exp(log_qmph)*BD;
 qrecl         = exp(log_qrecl);
 qpela         = exp(log_qpela);
 //===============================================================================
}

void model_parameters::Eval_PBR(void)
{
  //=============================================================================== 
  if(opt_Ro<0)
  {
    log_Ro=log_priorRo;
  }
  // Frms proxy (60%SPR y otros) y xx%SPR de Fmediana hist?rica
  for(int i=1;i<=npbr;i++){
    Fspr  = Sel_f(nanos)*log_Fref(i);
    Zspr  = Fspr+M;
    //Fmedian = Fmedian_ext;
    Fmedian = Fmedian1; 
    Fmed    = Sel_f(nanos)*Fmedian;
    Zmed    = Fmed+M;
    Nspro(1)=mfexp(log_Ro+0.5*square(sigmaR)); 
    Nspr(1)=mfexp(log_Ro+0.5*square(sigmaR)); 
    Nmed(1)=mfexp(log_Ro+0.5*square(sigmaR)); 
    for (int j=2;j<=nedades;j++)
    { 
      Nspro(j) = Nspro(j-1)*exp(-1*M);
      Nspr(j)  = Nspr(j-1)*exp(-Zspr(j-1));
      Nmed(j)  = Nmed(j-1)*exp(-Zmed(j-1));
    }
    Nspro(nedades) = Nspro(nedades)/(1-exp(-1*M));
    Nspr(nedades)  = Nspr(nedades)/(1-exp(-Zspr(nedades)));
    Nmed(nedades)  = Nmed(nedades)/(1-exp(-Zmed(nedades)));
    Bspro   = sum(elem_prod(Nspro*exp(-dt(3)*M),elem_prod(msex,colsum(Win)/nanos)));
    Bspr    = sum(elem_prod(elem_prod(elem_prod(Nspr,mfexp(-dt(3)*Zspr)),msex),colsum(Win)/nanos));
    Bsprmed = sum(elem_prod(elem_prod(elem_prod(Nmed,mfexp(-dt(3)*Zmed)),msex),colsum(Win)/nanos));
    ratio_spr(i) = Bspr/Bspro;	
    ratio_Fmed   = Bsprmed/Bspro;// xx%SPR de Fmediana
    // Bo y Brms proxy  seg?n metodolog?a Taller PBRs 2014
    //Bmed    = mean(BD(1,12));
    Bmed=mean(BD(1,nanos));
    Bo      = Bmed/(ratio_Fmed-0.05);
    Brms(i) = Bo*(ratio_spr(i)-0.05);
  }    	
}

void model_parameters::Eval_Estatus(void)
{
  SSB=BD(1,nanos);   // variables de inter?s para mcmc 
  Nv    = N;// solo para empezar los calculos
 for (int i=2;i<=nanos;i++)
  {
      Nv(i)(2,nedades)=++Nv(i-1)(1,nedades-1)*exp(-1.0*M);
      Nv(i)(nedades)+=Nv(i-1)(nedades)*exp(-1.0*M);
  }
  NDv  = elem_prod(Nv*exp(-dt(3)*M),outer_prod(Unos_anos,msex));
  BDo  = rowsum(elem_prod(NDv,Win));
  // INDICADORES DE REDUCCI?N DEL STOCK
  RPRdin =  elem_div(BD,BDo);                       // RPR BDspr_t, din?mico
  RPRequ =  BD/Bspro;                               // RPR con BDspro
  RPRequ2 = BD/Bo;                                 // RPR con Bo proxy
  RPRequ3 = BD/Brms(1);                            // Raz?n para diagrama de fase
  Frpr    = exp(log_Ft)/log_Fref(1);
}

void model_parameters::Eval_logverosim(void)
{
  int i;
  suma1=0;
  suma2=0;
  suma3=0;
  suma4=0;
  for (i=1;i<=nanos;i++)
  {
      if (Reclas(i)>0)
      {
         suma1   += square((log(Reclas(i))-log(Reclas_pred(i)))/cvar(1,i));
      }
      if (Pelaces(i)>0)
      {
         suma2   += square((log(Pelaces(i))-log(Pelaces_pred(i)))/cvar(2,i));
      }
      if (MPH(i)>0)
      {
         suma3   += square((log(MPH(i))-log(MPH_pred(i)))/cvar(3,i));
      }
      if (Desemb(i)>0)
      {
         suma4   += square((log(Desemb(i))-log(Desemb_pred(i)))/cvar(4,i));}
  }
}

void model_parameters::Eval_funcion_objetivo(void)
{
  likeval(1)   = 0.5*suma1;//Reclas
  likeval(2)   = 0.5*suma2;//pelaces
  likeval(3)   = 0.5*suma4;//Desemb
  likeval(4)   = 0.5*suma3;//MPH
  likeval(5)   = -nmus(1)*sum(elem_prod(pobs_f,log(ppred_f)));
  likeval(6)   = -nmus(2)*sum(rowsum(elem_prod(pobs_crua,log(ppred_crua))));
  likeval(7)   = -nmus(3)*sum(rowsum(elem_prod(pobs_pel,log(ppred_pel))));
  likeval(8)   = -nmus(4)*sum(rowsum(elem_prod(pobs_crul,log(ppred_crul))));
  likeval(9)   = 1./(2*square(sigmaR))*norm2(log_desv_Rt);
  likeval(10)  = 1./(2*square(cvpriorq_reclas))*square(log_qrecl);
  likeval(11)  = 1./(2*square(cvpriorq_pelaces))*square(log_qpela);
  //likeval(12) = 1000*(square(log_Ft(2)-mean(log_Ft))+square(log_Ft(3)-mean(log_Ft)));  // S12
  if(active(log_Fref)){
  likeval(12) = 1000*norm2(ratio_spr-ratio);}
  if (opt2_fase>0)
  {
    for (int i=1;i<=nedades-2;i++)
    {
       pSf(i) = lam(1)*sum(square(log(Sel_f(i))-2*log(Sel_f(i+1))+log(Sel_f(i+2))));
       pSc(i) = lam(2)*sum(square(log(Scru(i))-2*log(Scru(i+1))+log(Scru(i+2))));
    }
    for (int i=1;i<=nanos;i++)
    {
       pStotf(i)=lam(3)*sum(square(log(Sel_f(i))));
    }
    for (int i=1;i<=nanos;i++)
    {
    pStotc(i)=lam(4)*sum(square(log(Scru(i))));
    }
  }
   if (opt1_fase>0)
   {
   f = sum(likeval);
   }
   if (opt2_fase>0)
   {
   f = sum(likeval)+sum(pSf)+sum(pSc)+sum(pStotf)+sum(pStotc);
   }
   if(mceval_phase())
   {
   Eval_mcmc();
   }
}

void model_parameters::Eval_CTP(void)
{
  // Estimaci?n de CBA A?O BIOL?GICO
  //for (int i=1;i<=npbr;i++) // ciclo de PBR
  //{
   // Fpbr   = Sel_f(nanos)*log_Fref(1)*1;
    //Zpbr   = Fpbr+M;
    //CTP    = elem_prod(elem_div(Fpbr,Zpbr),elem_prod(1.-exp(-1.*Zpbr),N(nanos)));
    //YTP(i) = sum(elem_prod(CTP,Wmed(nanos)));                              
  //}
  //************************************************************************************************
  // Estimaci?n de CBA A?O BIOL?GICO sin proyectar//revisar ?ltimo a?o!!!
    if(FrmsH2==1){ // para cálculo de C_rms segundo hito marzo)
    Fref_r0 = log_Fref(1);}
    if(FrmsH2<1){ //para cálculo de C_actual primer hito julio y septiembre)
    Fref_r0 = exp(log_Ft(nanos));} //log_Fref(1);//aqu? deber?a ir F del ?ltimo a?o
    //Fref_r0 = exp(log_Ft(nanos)); //log_Fref(1);//aqu? deber?a ir F del ?ltimo a?o
  	Frms_r0 = Sel_f(nanos)*Fref_r0;
    Zrms_r0 = Frms_r0+M;
    CTP_r0  = elem_prod(elem_div(Frms_r0,Zrms_r0),elem_prod(1.-exp(-1.*Zrms_r0),N(nanos)));
    YTP_r0W  = elem_prod(CTP_r0,Wmed(nanos));  
    YTP_r0  = sum(elem_prod(CTP_r0,Wmed(nanos)));    
    YTP_r01ersem = YTP_r0 - Desemb2doSem; // Captura en toneladas del primer semestre (enero a junio)
    YTP_r01ersemR = YTP_r0 - (Desemb2doSem+Remanente); // Captura en toneladas del primer semestre (enero a junio)
    // descuento del descarte para CBA
    YTP_r0d  = YTP_r0*Descarte; // Captura en toneladas año biológico actual - el 4% de descarte supuesto
    YTP_r0d1ersem = YTP_r0d - Desemb2doSem ; // Captura en toneladas del primer semestre (enero a junio)
    YTP_r0d1ersemR = YTP_r0d - (Desemb2doSem+Remanente) ; // Captura en toneladas del primer semestre (enero a junio)
  	NMD_r0  = elem_prod(elem_prod(N(nanos),mfexp(-dt(3)*Zrms_r0)),msex);
    BD_r0   = sum(elem_prod(NMD_r0,Win(nanos)));	
    RPR_r0  = BD_r0/Brms(1);
  if(opt_Ro<0){ log_Ro=log_priorRo; }//Esto corre cuando se hace perfil de verosimilitud
    Np     = N(nanos);
    Sp     = S(nanos);
    Nvp    = Nv(nanos);
    RPRp(1)= RPRequ3(nanos);
   for (int j=1;j<=nproy;j++){ // ciclo de 5 a?os
    Np(2,nedades)=++elem_prod(Np(1,nedades-1),Sp(1,nedades-1));
    Np(nedades)+=Np(nedades)*Sp(nedades);
 // Escenarios de reclutamiento promedio
    if(oprec==1){Np(1)=mfexp(mean(log(Reclutas(1,nanos-13))));} // reclutamiento promedio 1997-2009 (primeros 13 a?os)
    if(oprec==2){Np(1)=mfexp(mean(log(Reclutas)));} // reclutamiento promedio hist?rico (1997-a?o m?s reciente)
    if(oprec==3){Np(1)=mfexp(mean(log(Reclutas(nanos-12,nanos))));}  // reclutamiento promedio ?ltimos a?os a partir del 2010 (2010-a?o m?s reciente)
    //if(oprec==1){Np(1)=mean(Reclutas(1,nanos-13));} // reclutamiento promedio 1997-2009 (primeros 13 a?os)
    //if(oprec==2){Np(1)=mean(Reclutas);} // reclutamiento promedio hist?rico (1997-a?o m?s reciente)
    //if(oprec==3){Np(1)=mean(Reclutas(nanos-12,nanos));}  // reclutamiento promedio ?ltimos a?os a partir del 2010 (2010-a?o m?s reciente)
    Npp = elem_prod(prop_est,Np);
    Wmedp=Wmedp_3;
	  Winp=Winip_3;
   Fref_p0 = mF*log_Fref(1);
   Frms_p0 = Sel_f(nanos)*Fref_p0;
   Zrms_p0 = Frms_p0+M;
   CTP_p0(j)    = elem_prod(elem_div(Frms_p0,Zrms_p0),elem_prod(1.-exp(-1.*Zrms_p0),Npp));
   YTP_p0W(j)   = elem_prod(CTP_p0(j),Wmedp);
   YTP_p0(j) = sum(elem_prod(CTP_p0(j),Wmedp)); 
   YTP_p0d(j) = (YTP_p0(j)*Descarte); //descuenta descarte
   BD_p0(j)  = sum(elem_prod(elem_prod(elem_prod(Npp,mfexp(-dt(3)*Zrms_p0)),msex),Winp)); 
   RPR_p0(j) = BD_p0(j)/Brms(1);
   //Nap(j)   = Npp;    
   Sp       = exp(-1.*Zrms_p0); 
   NVrecl_p0   = elem_prod(elem_prod(Npp,mfexp(-dt(1)*Zrms_p0)),Scru(nanos));//considerar s?lo mortalidad natural- Crucero Reclas!!!
   NVpel_p0    = elem_prod(elem_prod(Npp,mfexp(-dt(2)*Zrms_p0)),Scru_pela(nanos));
   Brecl_p0(j) = qrecl*sum(elem_prod(NVrecl_p0,Wmedp));
   Bpel_p0(j)  = qpela*sum(elem_prod(NVpel_p0,Winp)); 
	}		   
  //----------------------------------------------------------------
 // regla Fconstante=Frms (0 = Fconst, 1 = regla mixta, r = mismo a?o, p=proyectado)
   // if(opProy==1){CBA_c0=prop(1)*YTP_r0+prop(2)*YTP_p0(1);
  //                CBA_c0d=prop(1)*(YTP_r0*0.98)+prop(2)*(YTP_p0(2)*0.98); } 
    if(opProy==1) // Opci?n 1: 1era y 2da revisi?n (para el mismo a?o)
    {
      // regla Fconstante=Frms (0 = Fconst, 1 = regla mixta, r = mismo a?o, p=proyectado)
      CBA_c0  =  (prop(1)*YTP_r0) + (prop(2)*YTP_p0(1)); // proporción 70/30 sin descontar el descarte 
      CBA_c0D  = YTP_r01ersem + (prop(2)*YTP_p0(1)); // con descuento del desembarque 2do semestre y sin descontar el descarte 
      CBA_c0R = YTP_r01ersemR + (prop(2)*YTP_p0(1)); //  con descuento del desembarque 2do semestre y con descuento del remanente 1er semestre
      CBA_c0d = (prop(1)*YTP_r0d) + (prop(2)*YTP_p0d(1)); //  proporción 70/30 con descuento del descarte
      CBA_c0dD = YTP_r0d1ersem + (prop(2)*YTP_p0d(1)); //  ccon descuento del desembarque 2do semestre y con descuento del descarte
      CBA_c0dR = YTP_r0d1ersemR + (prop(2)*YTP_p0d(1)); //   con descuento del desembarque 2do semestre y con descuento del remanente 1er semestre y con descuento del descarte
    }
    if(opProy==2)
      {
      CBA_c0=(prop(1)*YTP_p0(1))+(prop(2)*YTP_p0(2));  // Opcion 2: CBA inicial (proyeccion de un ano calendario)
      CBA_c0d=(prop(1)*YTP_p0d(1))+(prop(2)*YTP_p0d(2)); 
      } 
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
  report << "years"<<endl;
  report << anos << endl;
  report << "reclasobs" << endl;
  report << Reclas << endl;
  report << "reclaspred" << endl;
  report << Reclas_pred << endl;
  report << "pelacesobs" << endl;
  report << Pelaces << endl;
  report << "pelacespred" << endl;
  report << Pelaces_pred << endl;
  report << "mphobs" << endl;
  report << MPH << endl;
  report << "mphpred" << endl;
  report << MPH_pred << endl;
  report << "desembarqueobs" << endl;
  report << Desemb << endl;
  report << "desembarquepred" << endl;
  report << Desemb_pred << endl;
  report << "Reclutas" << endl;
  report << Reclutas << endl;
  report << "log_desv_Rt" << endl;
  report << log_desv_Rt << endl;
  report << "SSB" << endl;
  report << BD << endl;
  report << "BT" << endl;
  report << BT << endl;
  report << "Ftot" << endl;
  report << exp(log_Ft) << endl;
  report << "BD_Bspro_t" << endl;
  report << RPRdin << endl;
  report << "BD_Bspro" << endl;
  report << RPRequ << endl;
  report << "BD_Bo" << endl;
  report << RPRequ2 << endl;
  report << "BD_Brms" << endl;
  report << RPRequ3 << endl;
  report << "S_f" << endl;
  report << Sel_f << endl;
  report << "Scru_reclas" << endl;
  report << Scru << endl;
  report << "Scru_pelaces" << endl;
  report << Scru_pela << endl;
  report << "pf_obs " << endl;
  report << pobs_f << endl;
  report << "pf_pred " << endl;
  report << ppred_f << endl;
  report << "pobs_RECLAS" << endl;
  report << pobs_crua << endl;
  report << "ppred_RECLAS" << endl;
  report << ppred_crua << endl;
  report << "pobs_PELACES" << endl;
  report << pobs_pel << endl;
  report << "ppred_PELACES" << endl;
  report << ppred_pel << endl;
  report << "pobs_pel_tallas" << endl;
  report << pobs_crul << endl;
  report << "ppred_pel_tallas" << endl;
  report << ppred_crul << endl;
  //----------------------------------------
  // PUNTOS BIOL?GICOS DE REFERENCIA TALLER
  //----------------------------------------
  report << "pSPR Fmed_Fpbrs"<<endl;
  report << ratio_Fmed<<"  "<<ratio_spr<<endl;
  report << "Fs Fmed_Fpbrs"<<endl;
  report << Fmedian <<"  "<<log_Fref<<endl;
  report << "SSBpbr Bo_Bmed_Bpbrs"<<endl;
  report <<  Bo <<"  "<< Bmed << " " << Brms<<endl;
  report << "Ro"<<endl;
  report <<  Nspro(1) << endl;
  report << "SPR SPRFo_SPRFmed_SPRFpbrs"<<endl;
  report << Bspro <<" " << Bsprmed <<" " << Bspr<< endl;
  report << "log_Ro" << endl;
  report << log_Ro << endl;
  report << "likeval ReclasPelacesDesembMPH_pf_preclas_ppelaces_ptallas_desvR_qrecl_qpela" << endl;
  report << likeval << endl;
  report << "q_reclas_q_pelaces" << endl;
  report << exp(log_qrecl)<<" "<<exp(log_qpela)<< endl;
  report << "M"<<endl;
  report << M << endl;
  report << "N" << endl;
  report << N << endl;
  report << "pred_Ctot" << endl;
  report << pred_Ctot << endl;
  report << "F" << endl;
  report << Ftot << endl;
  report << "Np" <<endl;
  report << Np <<endl;
  report << "YTP_r0W_actual" << endl;
  report << YTP_r0W << endl;
  report << "YTP_p0W_proyectada" << endl;
  report << YTP_p0W << endl;
  report << "CTP_r0_actual" << endl;
  report << CTP_r0 << endl;
  report << "CTP_p0_proyectada" << endl;
  report << CTP_p0 << endl;
  report << "Wmedp" << endl;
  report << Wmedp << endl;
 suma1=0; suma2=0;nm1=1;cuenta1=0;cuenta2=0;
  for (int i=1;i<=nanos;i++){ //
   if (sum(pobs_f(i))>0){
      suma1=sum(elem_prod(ppred_f(i),1-ppred_f(i)));
      suma2=norm2(pobs_f(i)-ppred_f(i));
      nm1=nm1*suma1/suma2;
      cuenta1+=1;
   }}
  suma1=0; suma2=0;nm2=1;cuenta2=0;
  for (int i=1;i<=nanos;i++){ //
   if (sum(pobs_crua(i))>0){
      suma1=sum(elem_prod(ppred_crua(i),1-ppred_crua(i)));
      suma2=norm2(pobs_crua(i)-ppred_crua(i));
      nm2=nm2*suma1/suma2;
      cuenta2+=1;
   }}
  suma1=0; suma2=0;nm3=1;cuenta3=0;
  for (int i=1;i<=nanos;i++){ //
   if (sum(pobs_pel(i))>0){
      suma1=sum(elem_prod(ppred_pel(i),1-ppred_pel(i)));
      suma2=norm2(pobs_pel(i)-ppred_pel(i));
      nm3=nm3*suma1/suma2;
      cuenta3+=1;
   }}
  suma1=0; suma2=0;nm4=1; cuenta1=0;cuenta4=0;
  for (int i=1;i<=nanos;i++){ //
   if (sum(pobs_crul(i))>0){
      suma1=sum(elem_prod(ppred_crul(i),1-ppred_crul(i)));
      suma2=norm2(pobs_crul(i)-ppred_crul(i));
      nm4=nm4*suma1/suma2;
      cuenta4+=1;
   }}
  report << "nmprior"<<endl;
  report <<nmus<<endl;
  report << "nm_flota  nm_reclas  nm_pelaces    nm_pelaL" << endl;
  report<<pow(nm1,1/cuenta1)<<" "<<pow(nm2,1/cuenta2)<<" "<<pow(nm3,1/cuenta3)<<" "<<pow(nm3,1/cuenta4)<<endl;
  suma1=0;  suma2=0;  suma3=0;   suma4=0;  cuenta1=0;     cuenta2=0;   cuenta3=0;
  for (int i=1;i<=nanos;i++)
  {
   if (Reclas(i)>0){
    suma1+=square(log(Reclas(i))-log(Reclas_pred(i)));
    cuenta1+=1;}
   if (Pelaces(i)>0){
    suma2+=square(log(Pelaces(i))-log(Pelaces_pred(i)));
    cuenta2+=1;}
   if (MPH(i)>0){
    suma4+=square(log(MPH(i))-log(MPH_pred(i)));
   cuenta4+=1;}
  }
 report << "cv_recla  cv_pelaces  cv_mph" << endl;
 report<<sqrt(suma1/cuenta1)<<" "<<sqrt(suma2/cuenta2)<<" "<<sqrt(suma4/cuenta4)<<endl;
  if(erredad==1){
  report << " ------------------------------------------------" << endl;
  report << "Matriz de error "<< endl;
  report << error_edad << endl;}
  report << " ------------------------------------------------" << endl;
  report << "Talla a la edad & desviaci?n "<< endl;
  report << mu_edad << endl;
  report << sigma_edad << endl;
}

void model_parameters::Eval_mcmc(void)
{
  if(reporte_mcmc == 0)
  mcmc_report<<"f, RPR, RBV, F_last, CTP_60, CTP_obj, Recl_last"<<endl;
  //mcmc_report<<f<<","<<RPRdin(nanos)<<","<<RPRequ(nanos)<<","<<max(Ftot(nanos))<<"," <<YTP(1)<<","<<YTP(2)<<","<<Reclutas(nanos)<<endl;
  reporte_mcmc++;
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

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
 //arrmblsize=60000000; // 
 // gradient_structure::set_GRADSTACK_BUFFER_SIZE(30000000);
 // gradient_structure::set_CMPDIF_BUFFER_SIZE(50000000);
 // gradient_structure::set_MAX_NVAR_OFFSET(100000);
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
