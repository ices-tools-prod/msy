// srmsymc.tpl: ADMB code developed by José De Oliveira at Cefas, UK, April 2010,
//  in collaboration with Chris Darby and Timothy Earl, also at Cefas.
// Features:
//  1. three stock recruit curves are fitted: Ricker, Beverton-Holt and Smooth Hockey-stick (after Mesnil in press)
//  2. data are normalised to the mean before being used to fit stock-recruit curves, to ease estimation
//  3. stock recruit parameters for Ricker and Beverton-Holt are re-parameterised (more orthoganal) to further ease estimation
//  4. program uses mcmc sampling feature of ADMB to sample stock-recruuit parameters from the pdf.
//     program therefore needs to be run first with "-mcmc n -mcsave m" and then with "-mceval" (n=number samples, m=degree of thinning, i.e. save every mth sample)
//  5. for the SHS, the functions get_y() and get_Bssb() are adjusted to give negatives for F>Fcrash because of the discontinuity at 2a 
//  6. various switch parameters are in use: Ropt, senopt, penopt (see comments below)
//
// Modified by TJE 18/10/2012 to calculate just per recruit statistics when Ropt = 0, for stocks with where none of the SR relationships are appropriate

GLOBALS_SECTION
  #include <admodel.h>
  ofstream simpar("simpar.dat");
  ofstream simparypr("simparypr.dat");
  ofstream simpary("simpary.dat");
  ofstream simparssbpr("simparssbpr.dat");
  ofstream simparssb("simparssb.dat");
  ofstream biopar("biopar.dat");
  random_number_generator rng(12135);
  const double pi = 3.141592654;
  const double NaN = 1.0/0.0;
  const int n_spr_vals = 5;
  const double spr_vals[n_spr_vals] = {0.20, 0.25, 0.30, 0.35, 0.40};  //Must be in ascending order
  
DATA_SECTION
  init_adstring stkname;
  init_adstring filname;
  init_int ybeg;
  init_int yend;
  init_int r; //recruitment age
  init_int A; //plusgroup age
  init_int Ropt; //0=None, 1=Ricker, 2=B-H, else=SHS (smooth hockey-stick)
  init_int simopt; //0=no sim, 1=run sim - NOT IN USE FOR THIS VERSION
  init_int senopt; //0=error only in recr, 1=error in recr & steady-state vectors
  init_int penopt; //0=don't constrain SR params, 1=constrain SR params
  init_matrix srdat(ybeg,yend,1,2);
 !! ad_comm::change_datafile_name(filname);
  init_int fno; //nr of fleets
  init_int sno; //fleet for ypr stats (0=total)
  init_number f; //proportion F before spawning
  init_number m; //proportion M before spawning
  init_matrix sdat(r,A,1,fno);
  init_matrix sdatcv(r,A,1,fno);
  init_matrix wdat(r,A,1,fno);
  init_matrix wdatcv(r,A,1,fno);
  init_matrix biodat(r,A,1,3);
  init_matrix biodatcv(r,A,1,3);

  int i;
  int n;
  int outlines;
  number Rav;
  number Bssbav;
  number Bssbmin;
  number Bssbmax;
  number apinit; //1st stock-recruit parameter (initial value)
  number bpinit; //2nd stock-recruit parameter (initial value)
  number gp; //3rd stock-recruit parameter for SHS model only
  number g; 
  vector R(ybeg,yend); //Recruitment
  vector Bssb(ybeg,yend); //SSB
  vector sc(r,A); //total selectivity
  vector sf(r,A); //selectivity of fleet for ypr stats
  vector wf(r,A); //mean weight from fleet for ypr stats
  vector M(r,A); //natural mortality
  vector Q(r,A); //maturity
  vector ws(r,A); //mean weight in stock
  matrix sdatx(r,A,1,fno);
  matrix wdatx(r,A,1,fno);
  matrix epss(r,A,1,fno);
  matrix epsw(r,A,1,fno);
  matrix epsbio(r,A,1,3);
  vector sigws(r,A);
  matrix sigwf(r,A,1,fno);

 LOCAL_CALCS
  outlines = 0;
  n = yend-ybeg-r+1;
  R = column(srdat,1);
  Bssb = column(srdat,2);
  Rav = mean(R);
  Bssbav = mean(Bssb);
  R /= Rav;
  Bssb /= Bssbav;
  Bssbmin = min(Bssb);
  Bssbmax = max(Bssb);
  gp = 0.001;
  g = gp*Bssbav;
  if (Ropt == 1) {
    apinit = mfexp(1-Bssb(ybeg));
    bpinit = Bssb(ybeg);
  } else if (Ropt == 2) {
    apinit = Bssb(ybeg)/2.;
    bpinit = (Bssb(ybeg)+1)/2.;
  } else {
    apinit = 1./(1.-gp/2.+sqrt(1.+gp*gp/4.));
    bpinit = 1.;
  }

  sc = rowsum(sdat);
  if (sno == 0) {
    sf = sc;
    for (i=r;i<=A;i++) wf(i) = sdat(i)*wdat(i)/sum(sdat(i));
  } else {
    sf = column(sdat,sno);
    wf = column(wdat,sno);
  }
  M = column(biodat,1);
  Q = column(biodat,2);
  ws = column(biodat,3);
  cout<<"----------------------------------------------------------"<<endl;
  cout<<"sample "<<outlines<<endl;
  cout<<"sc "<<sc<<endl;
  cout<<"sf "<<sf<<endl;
  cout<<"wf "<<wf<<endl;
  cout<<"M "<<M<<endl;
  cout<<"Q "<<Q<<endl;
  cout<<"ws "<<ws<<endl;
 END_CALCS

INITIALIZATION_SECTION
  ap apinit;
  bp bpinit;

PARAMETER_SECTION
  init_number ap;
  init_number bp;

  sdreport_number a;
  sdreport_number b;
  sdreport_number sigR;
  sdreport_number scor;

  number Fcrash;
  number Fmax;
  number F01;
  vector FXspr(1,n_spr_vals);;
  number FXsprNUM;
  number Fmsy;
  number MSYpr;
  sdreport_number MSY;
  number MSYR;
  number Bmsypr;
  sdreport_number Bmsy;
  number fnFcrash;
  number fnFmax;
  number fnF01;
  vector fnFXspr(1,n_spr_vals);
  //  number fnF40spr;
  number dYdF;
  vector Rh(ybeg+r,yend);
  vector epsR(ybeg+r,yend);
  number pen;
  objective_function_value nll;

PROCEDURE_SECTION
  get_likelihood();
  if (sd_phase()) get_SRpar();
  if (mceval_phase()) get_mcmc_outputs();

FUNCTION get_likelihood
  dvariable penx;
  pen = 0.;
  penx = posfun(ap,1.e-6,pen);
  penx = posfun(bp,1.e-6,pen);
  if (Ropt == 0)
    Rh = R + pow(ap,2) + pow(bp,2)+1; //Something to minimise, just to find ap and bp = 0
  else if (Ropt == 1)
    Rh = ap*elem_prod(Bssb(ybeg,yend-r),mfexp(-bp*(Bssb(ybeg,yend-r)/Bssb(ybeg)-1.))).shift(Rh.indexmin());
  else if (Ropt == 2) {
    Rh = elem_div(Bssb(ybeg,yend-r),bp+ap*(Bssb(ybeg,yend-r)/Bssb(ybeg)-1.)).shift(Rh.indexmin());
    penx = posfun((bp-ap)*Bssb(ybeg)*Bssbav/ap,1.e-6,pen);
  } else {
    Rh = ap*(Bssb(ybeg,yend-r)+sqrt(bp*bp+gp*gp/4.)-sqrt(square(Bssb(ybeg,yend-r)-bp)+gp*gp/4.)).shift(Rh.indexmin());
    penx = posfun(bp-Bssbmin,1.e-6,pen);
    penx = posfun(Bssbmax-bp,1.e-6,pen);
  }
  epsR = log(R(ybeg+r,yend))-log(Rh);
  sigR = norm(epsR)/sqrt(double(n));
  nll = 0.5*n*log(2*pi*sigR*sigR) + norm2(epsR)/(2*sigR*sigR);
  if (penopt == 1) nll += 1000000000.*pen;

FUNCTION get_mcmc_outputs
  int ix = 50;
  dvariable Fx;
  epss.fill_randn(rng);
  epsw.fill_randn(rng);
  epsbio.fill_randn(rng);
  if (outlines == 0) {
    int rv=system("srmsymc2 -noest");  
    ifstream srmsymc_pe("srmsymc_pe.dat");
    srmsymc_pe>>ap>>bp;
    srmsymc_pe.close();
    rv=system("del srmsymc_pe.dat");  
    get_likelihood();
    for (i=0;i<=ix;i++) {
      simparypr<<i*2./double(ix)<<" ";
      simpary<<i*2./double(ix)<<" ";
      simparssbpr<<i*2./double(ix)<<" ";
      simparssb<<i*2./double(ix)<<" ";
    }
  }
  if (outlines==0 || outlines>10) {
    if (senopt==1 && outlines>0) get_pseudodata();
    get_params();
    simpar<<outlines<<" "<<ap<<" "<<bp<<" "<<a<<" "<<b<<" "<<sigR<<" "<<scor<<" "
      <<Fcrash<<" "<<Fmax<<" "<<F01<<" "<<FXspr<<" "<<Fmsy<<" "
      <<MSY<<" "<<Bmsy<<" "<<MSYpr<<" "<<Bmsypr<<" "<<MSYR<<" "
      <<fnFcrash<<" "<<fnFmax<<" "<<fnF01<<" "<<fnFXspr<<" "<<dYdF<<" "
      <<pen<<" "<<nll<<" "<<2.*nll+2.*3.*double(n)/(double(n)-3.-1.)<<endl;
      if (outlines==0) {
        simparypr<<endl;
        simpary<<endl;
        simparssbpr<<endl;
        simparssb<<endl;
      }
    for (i=0;i<=ix;i++) {
      Fx = i*2./double(ix);
      simparypr<<get_ypr(Fx)<<" ";
      simparssbpr<<get_Bssbpr(Fx)<<" ";
      simpary<<get_y(Fx)<<" ";
      simparssb<<get_Bssb(Fx)<<" ";
    }
    simparypr<<endl;
    simpary<<endl;
    simparssbpr<<endl;
    simparssb<<endl;
  }
  outlines++;

FUNCTION get_pseudodata
  int ii;
  sdatx = sdat;
  sdatx = elem_prod(sdatx,1.0+elem_prod(epss,sdatcv));
  for (i=r;i<=A;i++) for (ii=1;ii<=fno;ii++) if (sdatx(i,ii)<0.) sdatx(i,ii)=0.;

  sigwf = sqrt(log(pow(wdatcv,2)+1.0)); //Convert CV to SD
  wdatx = elem_prod(wdat,mfexp(elem_prod(epsw,sigwf)-square(sigwf)/2.));

  sc = rowsum(sdatx);
  if (sno == 0) {
    sf = sc;
    for (i=r;i<=A;i++) wf(i) = sdatx(i)*wdatx(i)/sum(sdatx(i));
  } else {
    sf = column(sdatx,sno);
    wf = column(wdatx,sno);
  }

  M = column(biodat,1);
  Q = column(biodat,2);
  ws = column(biodat,3);

  M = elem_prod(M, 1.0+elem_prod(column(epsbio,1),column(biodatcv,1)));
  for (i=r;i<=A;i++) if (M(i)<0.001) M(i)=0.001;
  for (i=r;i<=A;i++) if (Q(i)==1.) Q(i)=1.-0.0000001;

  Q = log(elem_div(Q,1.-Q))
      + elem_prod(column(epsbio,2),elem_div(column(biodatcv,2),1.-Q));
  Q = elem_div(mfexp(Q),1.+mfexp(Q));
  sigws = sqrt(log(pow(column(biodatcv,3),2)+1.0)); //Convert CV to SD
  ws = elem_prod(ws,mfexp(elem_prod(column(epsbio,3),sigws)-square(sigws)/2.));
  cout<<endl;
  cout<<"----------------------------------------------------------"<<endl;
  cout<<"sample "<<outlines<<endl;
  cout<<"sc "<<sc<<endl;
  cout<<"sf "<<sf<<endl;
  cout<<"wf "<<wf<<endl;
  cout<<"M "<<M<<endl;
  cout<<"Q "<<Q<<endl;
  cout<<"ws "<<ws<<endl;
 //  biopar <<  sc << sf << wf << ws << M << Q <<endl;

FUNCTION get_params
  cout<<"----------------------------------------------------------"<<endl;
  cout<<"Ropt "<<Ropt<<" sample "<<outlines<<endl;
  get_SRpar();
  get_Fcrash();
  get_Fmax_F01();
  get_FXspr();
  get_Fmsy();
  cout<<"pen "<<pen<<" nll "<<nll<<" "<<" AICc "<<2.*nll+2.*3.*double(n)/(double(n)-3.-1.)<<endl;
  cout<<"----------------------------------------------------------"<<endl;
  cout<<endl;

FUNCTION get_SRpar
  scor = epsR(ybeg+r,yend-1)*epsR(ybeg+r+1,yend).shift(epsR.indexmin())
         /(norm(epsR(ybeg+r,yend-1))*norm(epsR(ybeg+r+1,yend)));
  if (Ropt == 0) {
    a = 0;
    b = 0;
  } else if (Ropt == 1) {
    a = ap*mfexp(bp)*Rav/Bssbav;
    b = bp/(Bssb(ybeg)*Bssbav);
  } else if (Ropt == 2) {
    a = Bssb(ybeg)*Rav/ap;
    b = (bp-ap)*Bssb(ybeg)*Bssbav/ap;
  } else {
    a = ap*Rav/Bssbav;
    b = bp*Bssbav;
  }

FUNCTION get_Fcrash
  dvariable Fmn,Fmx;

  Fmn = 0.;
  Fmx = 5.  ;
  for (i=1;i<=30;i++)
  {
    Fcrash = (Fmn+Fmx)/2.;
    if (Ropt == 1) fnFcrash = a-1./get_Bssbpr(Fcrash);
    else if (Ropt == 2) fnFcrash = a/b-1./get_Bssbpr(Fcrash);
    else fnFcrash = a*(1.+b/sqrt(b*b+g*g/4.))-1./get_Bssbpr(Fcrash);
    if (fnFcrash>0.) Fmn = Fcrash;
    else Fmx = Fcrash;
  }
  Fcrash = (Fmn+Fmx)/2.;
  if (Ropt == 1) fnFcrash = a-1./get_Bssbpr(Fcrash);
  else if (Ropt == 2) fnFcrash = a/b-1./get_Bssbpr(Fcrash);
  else fnFcrash = a*(1.+b/sqrt(b*b+g*g/4.))-1./get_Bssbpr(Fcrash);
  cout<<"Fcrash "<<Fcrash<<" fnFcrash "<<fnFcrash<<endl;

FUNCTION get_Fmax_F01
  dvariable dypr_dF_F0,Fmn,Fmx,F0;

  F0 = 0.;
  Fmn = 0.;
  Fmx = 3.;
  for (i=1;i<=30;i++)
  {
    Fmax = (Fmn+Fmx)/2.;
    fnFmax = dypr_dF(Fmax);
    if (fnFmax>0.) Fmn = Fmax;
    else Fmx = Fmax;
  }
  Fmax = (Fmn+Fmx)/2.;
  fnFmax = dypr_dF(Fmax);
  cout<<"Fmax "<<Fmax<<" fnFmax "<<fnFmax<<endl;

  Fmn = 0.;
  Fmx = Fmax;
  dypr_dF_F0 = dypr_dF(F0);
  for (i=1;i<=30;i++)
  {
    F01 = (Fmn+Fmx)/2.;
    fnF01 = dypr_dF(F01)-0.1*dypr_dF_F0;
    if (fnF01>0.) Fmn = F01;
    else Fmx = F01;
  }
  F01 = (Fmn+Fmx)/2.;
  fnF01 = dypr_dF(F01)-0.1*dypr_dF_F0;
  cout<<"F01 "<<F01<<" fnF01 "<<fnF01<<endl;

FUNCTION get_FXspr
  dvariable Bssbpr_F0,Fmn,Fmx,F0;
  

  F0 = 0.;
  Fmn = 0.;
  Fmx = 3.;
  Bssbpr_F0 = get_Bssbpr(F0);
  for (int j=1; j<=n_spr_vals; j++)
  {
    for (i=1;i<=35;i++)
    {
      FXsprNUM = (Fmn+Fmx)/2.;
      fnFXspr(j) = get_Bssbpr(FXsprNUM)-spr_vals[j-1]*Bssbpr_F0;
      if (fnFXspr(j)>0.) Fmn = FXsprNUM;
      else Fmx = FXsprNUM;
    }
    FXspr(j) = (Fmn+Fmx)/2.;
    FXsprNUM = (Fmn+Fmx)/2.;
    fnFXspr(j) = get_Bssbpr(FXsprNUM)-spr_vals[j-1]*Bssbpr_F0;
    cout<<"F" << spr_vals[j-1]*100 << "spr "<<FXspr(j) <<" fnF"<<spr_vals[j-1]*100 << "spr "<<fnFXspr(j)<<endl;
    Fmn = 0.;
    Fmx = FXspr(j);
  }


FUNCTION get_Fmsy
  dvariable Fmn,Fmx,Y1,Y2,F1,F2;

  if (Ropt !=0)
  {
  Fmn = 0.;
  Fmx = Fcrash;
  if (Fmx > 3.) Fmx = 3.;
  for (i=1;i<=40;i++)
  {
    Fmsy = (Fmn+Fmx)/2.;
    F1 = Fmsy + 0.00001;
    F2 = Fmsy - 0.00001;
    Y1 = get_y(F1);
    Y2 = get_y(F2);
    dYdF = (Y1-Y2)/0.00002;
    if (dYdF>0.) Fmn = Fmsy;
    else Fmx = Fmsy;
  }
  Fmsy = (Fmn+Fmx)/2.;
  MSYpr = get_ypr(Fmsy);
  MSY = get_y(Fmsy);
  Bmsypr = get_Bssbpr(Fmsy);
  Bmsy = get_Bssb(Fmsy);
  MSYR = MSY/Bmsy;
  cout<<"Fmsy "<<Fmsy<<" dYdF "<<dYdF<<endl;
  cout<<"MSY "<<MSY<<" Bmsy "<<Bmsy<<endl;
  cout<<"MSYpr "<<MSYpr<<" Bmsypr "<<Bmsypr<<" MSYR "<<MSYR<<endl;
  }

FUNCTION dvariable dypr_dF(dvariable& Fx)
  RETURN_ARRAYS_INCREMENT();
  dvariable dypr_dFx;
  dvar_vector sx = sc;
  dvar_vector esum = get_esum(Fx);
  dvar_vector e1min = get_e1min(Fx);
  dvar_vector Zprop = get_Zprop(Fx,sx);
  dvar_vector sFsum = get_sFsum(Fx);
  dvar_vector sFe = get_sFe(Fx);
  dypr_dFx = elem_div(elem_prod(wf,elem_prod(sf,esum)),sc*Fx+M)
             *(elem_prod(e1min,1.-Zprop-sFsum)+sFe);
  RETURN_ARRAYS_DECREMENT();
  return dypr_dFx;

FUNCTION dvariable get_ypr(dvariable& Fx)
  RETURN_ARRAYS_INCREMENT();
  dvariable  x;
  dvar_vector sx = sf;
  dvar_vector esum = get_esum(Fx);
  dvar_vector e1min = get_e1min(Fx);
  dvar_vector Zprop = get_Zprop(Fx,sx);
  x = elem_prod(esum,wf)*elem_prod(Zprop,e1min);
  RETURN_ARRAYS_DECREMENT();
  return x;

FUNCTION dvariable get_y(dvariable& Fx) //note adjustment of yield calc for Ropt=3
  RETURN_ARRAYS_INCREMENT();
  dvariable x,yprx,Bssbprx;
  yprx = get_ypr(Fx);
  Bssbprx = get_Bssbpr(Fx);
  if (Ropt == 1) x = yprx*log(a*Bssbprx)/(b*Bssbprx);
  else if (Ropt == 2) x = yprx*(a*Bssbprx-b)/Bssbprx;
  else x = yprx*2.*a*(sqrt(b*b+g*g/4.)*(1.-1./(a*Bssbprx))+b)/sqrt(square(2.-1./(a*Bssbprx))); //sqrt(square()) is to ensure yield is negative for F>Fcrash
  RETURN_ARRAYS_DECREMENT();
  return x;

FUNCTION dvariable get_Bssbpr(dvariable& Fx)
  RETURN_ARRAYS_INCREMENT();
  dvariable  x;
  dvar_vector esum = get_esum(Fx);
  dvar_vector essb = get_essb(Fx);
  x = elem_prod(esum,ws)*elem_prod(Q,essb);
  RETURN_ARRAYS_DECREMENT();
  return x;

FUNCTION dvariable get_Bssb(dvariable& Fx) //note adjustment of Bssb calc for Ropt=3
  RETURN_ARRAYS_INCREMENT();
  dvariable  x,Bssbprx;
  Bssbprx = get_Bssbpr(Fx);
  if (Ropt == 1) x = log(a*Bssbprx)/b;
  else if (Ropt == 2) x = a*Bssbprx-b;
  else x = 2.*a*Bssbprx*(sqrt(b*b+g*g/4.)*(1.-1./(a*Bssbprx))+b)/sqrt(square(2.-1./(a*Bssbprx))); //sqrt(square()) is to ensure Bssb is negative for F>Fcrash
  RETURN_ARRAYS_DECREMENT();
  return x;

FUNCTION dvar_vector get_esum(dvariable& Fx)
  RETURN_ARRAYS_INCREMENT();
  dvar_vector x(r,A);
  x(r) = 0.;
  for (int j=r+1;j<=A;j++) x(j) = -sum(sc(r,j-1)*Fx+M(r,j-1));
  x = mfexp(x);
  RETURN_ARRAYS_DECREMENT();
  return x;

FUNCTION dvar_vector get_e1min(dvariable& Fx)
  RETURN_ARRAYS_INCREMENT();
  dvar_vector x(r,A);
  x(r,A-1) = 1.-mfexp(-sc(r,A-1)*Fx-M(r,A-1));
  x(A) = 1.;
  RETURN_ARRAYS_DECREMENT();
  return x;

FUNCTION dvar_vector get_Zprop(dvariable& Fx,dvar_vector& sx)
  RETURN_ARRAYS_INCREMENT();
  dvar_vector x(r,A);
  x = elem_div(sx*Fx,sc*Fx+M);
  RETURN_ARRAYS_DECREMENT();
  return x;

FUNCTION dvar_vector get_sFsum(dvariable& Fx)
  RETURN_ARRAYS_INCREMENT();
  dvar_vector x(r,A);
  x(r) = 0.;
  for (int j=r+1;j<=A;j++) x(j) = sum(sc(r,j-1))*Fx;
  RETURN_ARRAYS_DECREMENT();
  return x;

FUNCTION dvar_vector get_sFe(dvariable& Fx)
  RETURN_ARRAYS_INCREMENT();
  dvar_vector x(r,A);
  x(r,A-1) = elem_prod(sc(r,A-1)*Fx,mfexp(-sc(r,A-1)*Fx-M(r,A-1)));
  x(A) = 0.;
  RETURN_ARRAYS_DECREMENT();
  return x;

FUNCTION dvar_vector get_essb(dvariable& Fx)
  RETURN_ARRAYS_INCREMENT();
  dvar_vector x(r,A);
  x = mfexp(-f*sc*Fx-m*M);
  x(A) /= (1.-mfexp(-sc(A)*Fx-M(A)));
  RETURN_ARRAYS_DECREMENT();
  return x;
