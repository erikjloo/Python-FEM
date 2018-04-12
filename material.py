

//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


MelroMaterial::MelroMaterial 

  ( idx_t rank, const Properties& globdat )
    : Super ( 3, globdat ), melroRank_ ( rank ), globdat_ ( globdat )

{
  poissonP_ = 0.;
  rmTolerance_ = 1.e-10;
  rmMaxIter_ = 10;

  v61_.resize ( 6 );
  v62_.resize ( 6 );
}


MelroMaterial::~MelroMaterial ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void MelroMaterial::configure ( const Properties& props )
{
  Super::configure ( props );

  props.find ( rmTolerance_, "rmTolerance" );
  props.get  ( poissonP_, "poissonP" );
  props.find ( rmMaxIter_, "rmMaxIter" );
  
  Ref<Function> sigmaC = makeFunc_ ( props, "sigmaC" );
  Ref<Function> sigmaT = makeFunc_ ( props, "sigmaT" );

  y_ = newInstance<YieldFunc_> 
       ( young_, poisson_, poissonP_, sigmaC, sigmaT );

  y_-> setRmSettings ( rmTolerance_, rmMaxIter_ );

  stateString_ = "PLANE_STRAIN";

  if ( melroRank_ == 2 )
  {
    props.find( stateString_, STATE_PROP );

    if      ( stateString_ == "PLANE_STRAIN" )
    {
      state_ = PlaneStrain;
    }
    else if ( stateString_ == "PLANE_STRESS" )
    {
      state_ = PlaneStress;
    }
    else if ( stateString_ == "AXISYMMETRIC" )
    {
      state_ = AxiSymmetric;
    }
  }

  historyNames_.resize ( 9 );
  historyNames_[0] = "epsp_xx";
  historyNames_[1] = "epsp_yy";
  historyNames_[2] = "epsp_zz";
  historyNames_[3] = "epsp_xy";
  historyNames_[4] = "epsp_yz";
  historyNames_[5] = "epsp_zy";
  historyNames_[6] = "epspeq";
  historyNames_[7] = "diss";
  historyNames_[8] = "loading";
}

//-----------------------------------------------------------------------
//   makeFunc_
//-----------------------------------------------------------------------

Ref<Function> MelroMaterial::makeFunc_ 

  ( const Properties& props,
    const String&     name ) const

{
  String args = "x";
  props.find ( args, "args" );

  Ref<Function> func = FuncUtils::newFunc ( args, name, props, globdat_ );

  FuncUtils::resolve ( *func, globdat_ );

  return func;
}

//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void MelroMaterial::getConfig ( const Properties& conf ) const
{
  Super::getConfig ( conf );

  conf.set ( "rmMaxIter", rmMaxIter_ );
  conf.set ( "rmTolerance", rmTolerance_ );
  conf.set ( "poissonP", poissonP_ );

  FuncUtils::getConfig ( conf, y_->getSigmaCFunc(), "sigmaC" );
  FuncUtils::getConfig ( conf, y_->getSigmaTFunc(), "sigmaT" );

  if ( melroRank_ == 2 )
  {
    // note that state is not set, because Super has rank 3

    conf.set ( STATE_PROP, stateString_ );
  }
}

//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void MelroMaterial::update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      idx_t                 ipoint )

{
  Tuple<double,6,6>  dmat;
  Vec6   eps   ( fill3DStrain ( strain ) ); 

  for ( idx_t i = 0; i < 6; ++i ) 
  {
    eps[i] -= thermStrain_[i];
  }

  bool   loading = false;
  double dgam;

  Vec6   sig;

  double epspeq0 = preHist_[ipoint].epspeq;
  Vec6   sigtr = tmatmul ( stiffMat_, (eps-preHist_[ipoint].epsp) );

  StressInvariants inv ( sigtr );

  if ( y_->isPlastic ( inv, epspeq0 ) )
  {
    try
    {
      // classical return mapping scheme starting with dgam = 0.

      dgam = y_->findRoot ( 0. );
    }
    catch ( const Exception& ex )
    {
      handleException_ ( ex, strain, ipoint );

      // more robust return mapping scheme

      double gmin, gmax, fmin, fmax;

      y_->findBounds ( gmin, gmax, fmin, fmax );

      dgam = y_->findRoot ( gmin, gmax, fmin, fmax );
    }

    // compute stress

    double ptr = inv.getFirstInvariant() / 3.;
    Vec6   sigDtr = deviatoric ( sigtr, ptr );

    double p    = ptr    / y_->getZetaP();
    Vec6   sigD = sigDtr / y_->getZetaS();

    sig     = sigD;
    sig[0] += p;
    sig[1] += p;
    sig[2] += p;
    
    // compute plastic strain increment 

    double alphaP23 = y_->getAlpha() * p * 2./3.;

    Vec6 mvec;

    mvec[0]  = 3. * sigD[0] + alphaP23;
    mvec[1]  = 3. * sigD[1] + alphaP23;
    mvec[2]  = 3. * sigD[2] + alphaP23;
    mvec[3]  = 6. * sigD[3];  // engng strain factor 2
    mvec[4]  = 6. * sigD[4];  
    mvec[5]  = 6. * sigD[5]; 

    Vec6 depsp = dgam * mvec;

    // consistent tangent

    double betam, pbm, hfac, beta, phi, rho, chi, psi, xi, omega;

    y_->getTangentParameters 
      ( betam, pbm, hfac, beta, phi, rho, chi, psi, xi, omega, dgam );

    Tuple<double,6,6> dmde = aI4_plus_bII_ ( betam, pbm );

    Vec6 hvec = hfac * matmul ( dmde, mvec );

    dmat = aI4_plus_bII_ ( beta, phi-beta/3. );

    for ( idx_t i = 0; i < 6; ++i )
    {
      for ( idx_t j = 0; j < 6; ++j )
      {
        dmat(i,j) -= chi * sigDtr[i] * sigDtr[j];
        dmat(i,j) -= omega * sigDtr[i] * hvec[j];
      }
      for ( idx_t j = 0; j < 3; ++j )
      {
        dmat(i,j) -= rho * sigDtr[i];
        dmat(j,i) -= psi * sigDtr[i];
        dmat(j,i) -= xi  * hvec[i];
      }
    }

    // update history

    double dG = dot ( sig, depsp );
    newHist_[ipoint].epsp = preHist_[ipoint].epsp + depsp;
    newHist_[ipoint].epspeq = y_->getEpspeq();
    newHist_[ipoint].dissipation = preHist_[ipoint].dissipation + dG;
    loading = true;
  }
  else
  {
    sig = sigtr;
    for ( idx_t i = 0; i < 6; ++i )
    {
      for ( idx_t j = 0; j < 6; ++j )
      {
        dmat(i,j) = stiffMat_(i,j);
      }
    }
    newHist_[ipoint].epsp = preHist_[ipoint].epsp;
    newHist_[ipoint].epspeq = preHist_[ipoint].epspeq;
    newHist_[ipoint].dissipation = preHist_[ipoint].dissipation;

  }
  newHist_[ipoint].loading = loading;
  reduce3DVector ( stress, sig );

  latestHist_ = &newHist_;

  // use initial stiffness for oscillating point in desparateMode

  if ( desperateMode_ )
  {
    if ( ! useSecant_[ipoint] )
    {
      bool switches = ( loading != preHist_[ipoint].loading );

      if ( switches != hasSwitched_[ipoint] )
      {
        useSecant_[ipoint] = true;
      }
    }
    if ( useSecant_[ipoint] )
    {
      for ( idx_t i = 0; i < 6; ++i )
      {
        for ( idx_t j = 0; j < 6; ++j )
        {
          dmat(i,j) = stiffMat_(i,j);
        }
      }
    }
  }
  reduce3DMatrix ( stiff, dmat );
}

//-----------------------------------------------------------------------
//   handleException_
//-----------------------------------------------------------------------

void MelroMaterial::handleException_ 

  ( const Exception&      ex,
    const Vector&         strain,
    const idx_t           ipoint )

{
  Ref<PrintWriter> out = newInstance<PrintWriter>
    ( &System::debug("melro") );

  out->nformat.setFractionDigits ( 10 );

  if ( ex.what().equals ( "No convergence" ) )
  {
    *out << "No convergence in return mapping algorithm " << endl;
  }
  else if ( ex.what().equals ( "nan" ) )
  {
    *out << "NaN detected in return mapping algorithm " << endl;
  }
  else if ( ex.what().equals ( "Negative dgam" ) )
  {
    *out << "Negative increment found" << endl;
  }
  else
  {
    *out << "Caught unkown exception: " << ex.what() << endl;
  }
  *out << "epsp " << preHist_[ipoint].epsp << endl;
  *out << "epspeq0 " << preHist_[ipoint].epspeq << endl;
  *out << "strain " << strain << endl << endl;;
}

//-----------------------------------------------------------------------
//   getDissipationStress
//-----------------------------------------------------------------------

void MelroMaterial::getDissipationStress

    ( const Vector&   sstar,
      const Vector&   strain,
      const idx_t     ipoint ) const

{
  Vec6   eps   ( fill3DStrain ( strain ) ); 
  Vec6   sstar6;
  
  // NB, use preHist. After commit, this contains the new values,
  // after cancel the old ones. In both cases the starting values
  // for the next solve
  
  double epspeq = preHist_[ipoint].epspeq;
  Vec6   sig = tmatmul ( stiffMat_, (eps-preHist_[ipoint].epsp) );

  StressInvariants inv ( sig );

  if ( y_->isPlastic ( inv, epspeq ) )
  {
    double p = inv.getFirstInvariant() /3.;
    Vec6   sigD = deviatoric ( sig, p );

    double beta, phi, rho, chi, psi;
    y_->getTangentParameters ( beta, phi, rho, chi, psi );

    Tuple<double,6,6>  dmat = aI4_plus_bII_ ( beta, phi-beta/3. );

    for ( idx_t i = 0; i < 6; ++i )
    {
      for ( idx_t j = 0; j < 6; ++j )
      {
        dmat(i,j) -= chi * sigD[i] * sigD[j];
      }
      for ( idx_t j = 0; j < 3; ++j )
      {
        dmat(i,j) -= rho * sigD[i];
        dmat(j,i) -= psi * sigD[i];
      }
    }

    Vec6 epstmp = 2.*preHist_[ipoint].epsp-eps;
    sstar6 = sig + matmul ( dmat.transpose(), epstmp );
  }
  else
  {
    Tuple<double,6,6>  dmat;

    for ( idx_t i = 0; i < 6; ++i )
    {
      for ( idx_t j = 0; j < 6; ++j )
      {
        dmat(i,j) = stiffMat_(i,j);
      }
    }
    sstar6 = matmul ( dmat, preHist_[ipoint].epsp );
  }
  reduce3DVector ( sstar, sstar6 );
}

//-----------------------------------------------------------------------
//   getHistory
//-----------------------------------------------------------------------

void MelroMaterial::getHistory

  ( const Vector&  hvals,
    const idx_t    mpoint ) const

{
  (*latestHist_)[mpoint].toVector ( hvals );
}

//-----------------------------------------------------------------------
//   setHistory
//-----------------------------------------------------------------------

void MelroMaterial::setHistory

  ( const Vec6&    epsp,
    const double   epspeq,
    const idx_t    ipoint )

{
  preHist_[ipoint].epspeq = epspeq;
  preHist_[ipoint].epsp   = epsp;
}

//-----------------------------------------------------------------------
//   aI4_plus_bII_
//-----------------------------------------------------------------------

Tuple<double,6,6> MelroMaterial::aI4_plus_bII_ 

  ( const double        a,
    const double        b ) const

{
  // make 6x6 matrix with a*I_4^s+b*II in Voigt notation

  Tuple<double,6,6> ret;
  ret = 0.;

  ret(0,0) += a+b;  ret(0,1) += b;    ret(0,2) += b;
  ret(1,0) += b;    ret(1,1) += a+b;  ret(1,2) += b;
  ret(2,0) += b;    ret(2,1) += b;    ret(2,2) += a+b;
  ret(3,3) += .5*a;
  ret(4,4) += .5*a;
  ret(5,5) += .5*a;

  return ret;
}

//-----------------------------------------------------------------------
//   commit
//-----------------------------------------------------------------------

void MelroMaterial::commit () 

{
  newHist_.swap ( preHist_ );

  latestHist_ = &preHist_;
}

//-----------------------------------------------------------------------
//   clone
//-----------------------------------------------------------------------

Ref<Material> MelroMaterial::clone () const

{
  // use default copy constructor

  return newInstance<MelroMaterial> ( *this );
}

//-----------------------------------------------------------------------
//   deviatoric
//-----------------------------------------------------------------------

Vec6 MelroMaterial::deviatoric

  ( const Vec6&   full,
    const double  p ) const

{
  Vec6 ret = full;

  ret[0] -= p;
  ret[1] -= p;
  ret[2] -= p;

  return ret;
}

//-----------------------------------------------------------------------
//   giveHistory
//-----------------------------------------------------------------------

double MelroMaterial::giveHistory ( const idx_t ip ) const
{
  return (*latestHist_)[ip].epspeq;
}

//-----------------------------------------------------------------------
//   allocPoints
//-----------------------------------------------------------------------

void MelroMaterial::allocPoints

  ( idx_t count )

{
  for ( idx_t i = 0; i < count; ++i )
  {
    preHist_.pushBack ( Hist_() );
    newHist_.pushBack ( Hist_() );
  }
}

//-----------------------------------------------------------------------
//   giveDissipation
//-----------------------------------------------------------------------

double MelroMaterial::giveDissipation ( const idx_t ipoint ) const
{   
  return (*latestHist_)[ipoint].dissipation;
}

//-----------------------------------------------------------------------
//   tmatmul
//-----------------------------------------------------------------------

Vec6 MelroMaterial::tmatmul

  ( const Matrix& mat,
    const Vec6&   t6 ) const

{
  t6_to_v6 ( v61_, t6 );

  matmul ( v62_, mat, v61_ );

  return v6_to_t6 ( v62_ );

}

//-----------------------------------------------------------------------
//   Hist_ constructor
//-----------------------------------------------------------------------

MelroMaterial::Hist_::Hist_ () : epspeq ( 0. ), dissipation ( 0. )
{
  epsp = 0.;
  loading = false;
}

//-----------------------------------------------------------------------
//   Hist_ print function
//-----------------------------------------------------------------------

void MelroMaterial::Hist_::print () const

{
  System::out() << "epsp " << epsp << ", epspeq " << epspeq << endl;
}

// -------------------------------------------------------------------
//  Hist_::toVector
// -------------------------------------------------------------------

inline void MelroMaterial::Hist_::toVector

 ( const Vector&  vec ) const

{
  vec[0] = epsp[0];
  vec[1] = epsp[1];
  vec[2] = epsp[2];
  vec[3] = epsp[3];
  vec[4] = epsp[4];
  vec[5] = epsp[5];
  vec[6] = epspeq;
  vec[7] = dissipation;
  vec[8] = loading;
}

//-----------------------------------------------------------------------
//   YieldFunc_ constructor
//-----------------------------------------------------------------------

MelroMaterial::YieldFunc_::YieldFunc_

  ( const double        young,
    const double        poisson,
    const double        poissonP,
    const Ref<Function> sigC,
    const Ref<Function> sigT )

  : sigmaC_ ( sigC ), sigmaT_ ( sigT ), rmTol_(0.), maxIter_(0)

{
  G_ = young / 2. / ( 1. + poisson );
  K_ = young / 3. / ( 1. - 2. * poisson );
  alpha_ = ( 4.5 - 9. * poissonP ) / ( 1. + poissonP );
  alpha2_ = alpha_ * alpha_;
  Ka_ = K_ * alpha_;
  nuPFac_ = sqrt ( 1. / ( 1. + 2. * poissonP * poissonP ) );

  sigC_ = sigmaC_->eval ( 0. );
  sigT_ = sigmaT_->eval ( 0. );
}

//-----------------------------------------------------------------------
//   YieldFunc_::setRmSettings
//-----------------------------------------------------------------------

void MelroMaterial::YieldFunc_::setRmSettings 

  ( const double rmTolerance, 
    const idx_t  rmMaxIter )

{
  rmTol_ = rmTolerance;
  maxIter_ = rmMaxIter;
}

//-----------------------------------------------------------------------
//   YieldFunc_::isPlastic
//-----------------------------------------------------------------------

bool MelroMaterial::YieldFunc_::isPlastic

  ( const Invariants&  inv,
    const double       epspeq0 )

{
  // update all variables that are constant inside the return 
  // mapping scheme

  epspeq0_ = epspeq0;

  j2tr_ = inv.getJ2();
  i1tr_ = inv.getFirstInvariant();
  j2fac_ = 18.*j2tr_;
  i1fac_ = 4.*alpha2_*i1tr_*i1tr_ / 27.;

  // evaluate failure criterion for the trial state

  double crit = evalTrial_();

  return (crit/sigC_/sigT_) > -rmTol_;
}

//-----------------------------------------------------------------------
//   YieldFunc_::findRoot
//-----------------------------------------------------------------------

double MelroMaterial::YieldFunc_::findRoot 

  ( const double dgam0 )

{
  double dgam = dgam0;
  double oldddgam = -1.;
  double oldcrit = evalTrial_();

  for ( idx_t irm = 0; irm < maxIter_; )
  {
    double crit = eval_( dgam );

    double normalized = jem::numeric::abs(crit) / ( sigC_*sigT_  );

    if ( normalized < rmTol_ )
    {
      // System::out() << " return mapping converged in " << irm 
      // << " iterations " << normalized << endl;
      break;
    }
    if ( ++irm == maxIter_ )
    {
      throw Exception ( JEM_FUNC, "No convergence" );
    }
    if ( jem::Float::isNaN(dcdg_) )
    {
      throw Exception ( JEM_FUNC, "nan" );
    }
    // System::out() << " return mapping iteration " << irm << " dgam " <<
    // dgam << ", criterion: " << crit << ", dcdg " << dcdg << endl;

    double ddgam = crit / dcdg_;

    if ( ddgam * oldddgam < 0. )
    {
      // divergence detection
      if ( oldcrit * crit < 0. && jem::numeric::abs(ddgam) > jem::numeric::abs(oldddgam) )
      {
        // there might be an inflection point around the root
        // use linear interpolation rather than linearization

        // System::out() << "divergence prevention" << endl;
        ddgam = -oldddgam * crit / ( crit - oldcrit );
      }
    }

    dgam -= ddgam;

    oldddgam = ddgam;
    oldcrit = crit;
  }

  if ( dgam < -1.e-12 )
  {
    throw Exception ( JEM_FUNC, "Negative dgam" );
  }
  return dgam;
}

//-----------------------------------------------------------------------
//   YieldFunc_::findBounds
//-----------------------------------------------------------------------

void  MelroMaterial::YieldFunc_::findBounds

  (       double&      gmin,
          double&      gmax,
          double&      fmin,
          double&      fmax )

{
  gmin = 0.;
  gmax = -1.;
  fmax = -1.;

  fmin = evalTrial_();

  JEM_ASSERT ( fmin > 0. );

  double dgam = 1.e-16;

  while ( gmax < 0. )
  {
    double crit = eval_( dgam );

    if ( crit > 0. )
    {
      gmin = dgam;
      fmin = crit;
    }
    else
    {
      gmax = dgam;
      fmax = crit;
    }
    dgam *= 10.;
  }
}

//-----------------------------------------------------------------------
//   YieldFunc_::improveBounds_
//-----------------------------------------------------------------------

void  MelroMaterial::YieldFunc_::improveBounds_

  (       double&      gmin,
          double&      gmax,
          double&      fmin,
          double&      fmax )

{
  idx_t np = 10;
  double dg = ( gmax - gmin ) / double(np);
  double dgam = gmin;

  while ( dgam < gmax )
  {
    dgam += dg;
    double crit = eval_( dgam );

    if ( crit > 0. )
    { 
      gmin = dgam;
      fmin = crit;
    }
    else
    {
      gmax = dgam;
      fmax = crit;
      break;
    }
  }
}

//-----------------------------------------------------------------------
//   YieldFunc_::findRoot (with bounds)
//-----------------------------------------------------------------------

double MelroMaterial::YieldFunc_::findRoot

  (       double       gmin,
          double       gmax,
          double       fmin,
          double       fmax )

{
  double dgam;
  try 
  {
    double dgam0 = estimateRoot_ ( gmin, gmax, fmin, fmax );

    dgam = findRoot ( dgam0 );
  }
  catch ( const Exception& ex )
  {
    improveBounds_ ( gmin, gmax, fmin, fmax );

    dgam = findRoot ( gmin, gmax, fmin, fmax );
  }
  return dgam;
}

//-----------------------------------------------------------------------
//   YieldFunc_::estimateRoot_
//-----------------------------------------------------------------------

double MelroMaterial::YieldFunc_::estimateRoot_

  ( const double       gmin,
    const double       gmax,
    const double       fmin,
    const double       fmax ) const

{
  // linear interpolation

  // fmin  o___ 
  //           ---___
  // phi=0 |----------o-----------|------
  //                      ---___
  // fmax  -                    --o
  //     gmin        ret        gmax

  return gmin + fmin * (gmax-gmin) / (fmin-fmax);
}

    
//-----------------------------------------------------------------------
//   YieldFunc_::evalTrial_
//-----------------------------------------------------------------------

double MelroMaterial::YieldFunc_::evalTrial_ () const

{
  // fast version of eval, without computing all dgam dependent variables

  double sigC = sigmaC_->eval ( epspeq0_ );
  double sigT = sigmaT_->eval ( epspeq0_ );
  double sigct = sigC - sigT;

  return 6.*j2tr_ + 2.*i1tr_*sigct - 2.*sigC*sigT;
}

//-----------------------------------------------------------------------
//   YieldFunc_::eval_
//-----------------------------------------------------------------------

double MelroMaterial::YieldFunc_::eval_

  ( const double dgam )

{
  setDGam_ ( dgam );

  return 6.*j2tr_/zS2_ + 2.*i1tr_*sigct_/zetaP_ - 2.*sigC_*sigT_;
}

//-----------------------------------------------------------------------
//   YieldFunc_::resetDGam
//-----------------------------------------------------------------------

void MelroMaterial::YieldFunc_::resetDGam ()

{
  // update all variables that are a function of dgam, using dgam=0.

  zetaS_ = 1.;
  zetaP_ = 1.;

  zS2_ = 1.;
  zP2_ = 1.;

  double j2fac = 18.*j2tr_;
  double i1fac = 4.*alpha2_*i1tr_*i1tr_ / 27.;
  double sqrtA = sqrt ( j2fac + i1fac );

  depspeq_  = 0.;
  epspeq_ = epspeq0_;

  sigC_ = sigmaC_->eval ( epspeq_ );
  sigT_ = sigmaT_->eval ( epspeq_ );
  sigct_ = ( sigC_ - sigT_ );

  HC_ = sigmaC_->deriv ( epspeq_ );
  HT_ = sigmaT_->deriv ( epspeq_ );

  gjzs3_ = G_ * j2tr_;
  kaizp2_ = Ka_ * i1tr_;

  double depedg = nuPFac_ * sqrtA;

  double dsigCe = HC_ * depedg;
  double dsigTe = HT_ * depedg;

  dcdg_ = 2.*i1tr_*(dsigCe-dsigTe) - 4.*kaizp2_*sigct_
          - 72.*gjzs3_ - 2.*(sigC_*dsigTe+sigT_*dsigCe);
}

//-----------------------------------------------------------------------
//   YieldFunc_::setDGam_
//-----------------------------------------------------------------------

void MelroMaterial::YieldFunc_::setDGam_

  ( const double dgam ) 

{
  // update all variables that are a function of dgam

  zetaS_ = 1. + 6. * G_ * dgam;
  zetaP_ = 1. + 2. * Ka_ * dgam;

  zS2_ = zetaS_*zetaS_;
  zP2_ = zetaP_*zetaP_;

  double j2fac = 18.*j2tr_;
  double i1fac = 4.*alpha2_*i1tr_*i1tr_ / 27.;
  double sqrtA = sqrt ( j2fac / zS2_ + i1fac / zP2_ );

  depspeq_  = nuPFac_ * dgam * sqrtA;
  epspeq_ = epspeq0_ + depspeq_;

  sigC_ = sigmaC_->eval ( epspeq_ );
  sigT_ = sigmaT_->eval ( epspeq_ );
  sigct_ = ( sigC_ - sigT_ );

  HC_ = sigmaC_->deriv ( epspeq_ );
  HT_ = sigmaT_->deriv ( epspeq_ );

  gjzs3_ = G_ * j2tr_ / zS2_ / zetaS_;
  kaizp2_ = Ka_ * i1tr_ / zP2_;

  // NB adding an additional I1tr in the final term of depedg
  // (Melro's paper is incorrect here)

  double depedg = nuPFac_ * ( sqrtA - .5 * dgam / sqrtA * 
      ( 216.*gjzs3_ + 16.*alpha2_*kaizp2_*i1tr_/27./zetaP_ ) );

  double dsigCe = HC_ * depedg;
  double dsigTe = HT_ * depedg;

  dcdg_ = 2.*i1tr_*(dsigCe-dsigTe)/zetaP_ - 4.*kaizp2_*sigct_
          - 72.*gjzs3_ - 2.*(sigC_*dsigTe+sigT_*dsigCe);
}

//-----------------------------------------------------------------------
//   YieldFunc_::getTangentParameters
//-----------------------------------------------------------------------

void MelroMaterial::YieldFunc_::getTangentParameters 

  (       double&  betam, 
          double&  pbm, 
          double&  hfac, 
          double&  beta, 
          double&  phi, 
          double&  rho, 
          double&  chi, 
          double&  psi, 
          double&  xi, 
          double&  omega, 
    const double   dgam )

{
  setDGam_ ( dgam );

  betam = 6. * G_ / zetaS_;
  pbm = 2.*Ka_/3./zetaP_ - betam/3.;

  if ( jem::numeric::abs(depspeq_)>1.e-20 )
  {
    hfac = 2. * ( ( HC_ - HT_ ) * i1tr_ / zetaP_ - (sigC_*HT_+sigT_*HC_) ) 
                * dgam*dgam*nuPFac_*nuPFac_/depspeq_;
  }
  else
  {
    hfac =  0.;
  }

  double eta  = -dcdg_;
  beta = 2. * G_ / zetaS_;
  phi  = K_ / zetaP_ * ( 1. - 4. * kaizp2_ * sigct_ / eta );
  rho  = 36. * K_ * G_ * sigct_ / eta / zS2_ / zetaP_;
  chi  = 72. * G_ * G_ / eta / zS2_ / zS2_;
  psi  = 8. * G_ * kaizp2_  / eta / zS2_;
  xi   = 2. * kaizp2_ / eta / 3.;
  omega = 6. * G_ / zS2_ / eta;
}

//-----------------------------------------------------------------------
//   YieldFunc_::getTangentParameters
//-----------------------------------------------------------------------

void MelroMaterial::YieldFunc_::getTangentParameters

  (       double&  beta, 
          double&  phi, 
          double&  rho, 
          double&  chi, 
          double&  psi )

{
  // fast version of getTangentParameters for dgam = 0

  resetDGam();

  double eta = -dcdg_;

  beta = 2. * G_;
  phi  = K_ * ( 1. - 4. * kaizp2_ * sigct_ / eta );
  rho = 36. * K_ * G_ * sigct_ / eta;
  chi = 72. * G_ * G_ / eta;
  psi = 8. * G_ * kaizp2_ / eta;
}

//-----------------------------------------------------------------------
//   YieldFunc_::getPGradientParameters
//-----------------------------------------------------------------------

void MelroMaterial::YieldFunc_::getPGradientParameters

  (       double& phiF,
          double& rhoF,
          double& chiF,
          double& psiF ) const

{
  double eta = -dcdg_;

  phiF = 4. * kaizp2_ * sigct_ / 3. / eta;
  rhoF = 18. * K_ * sigct_ / eta;
  chiF = 36. * G_ / eta;
  psiF = 8. * G_ * alpha_ * i1tr_ / 3. / eta;
}
