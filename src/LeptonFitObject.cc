/*! \file 
 *  \brief Implements class LeptonFitObject
 *  LeptonFitObject works similiar to JetFitObject, but it uses a 1/pt, eta, phi parametrization for the 
 *  leptons, which is e.g. more appropriate for electrons.
 *  Especially the covarianz matrix differs from a common E, theta, phi parametrization.
 */ 

#include "LeptonFitObject.h"
#include <cmath>
#include <cassert>
#include <iostream>

using std::sqrt;
using std::sin;
using std::cos;
using std::cout; 
using std::endl;

// constructor
LeptonFitObject::LeptonFitObject(double Pt, double theta, double phi,  
                           double DPt, double Dtheta, double Dphi, 
                           double m) {
  initCov();                         
//  assert( !isinf(E) );        assert( !isnan(E) );
//  assert( !isinf(theta) );    assert( !isnan(theta) );
//  assert( !isinf(phi) );      assert( !isnan(phi) );
//  assert( !isinf(DE) );       assert( !isnan(DE) );
//  assert( !isinf(Dtheta) );   assert( !isnan(Dtheta) );
//  assert( !isinf(Dphi) );     assert( !isnan(Dphi) );
//  assert( !isinf(m) );        assert( !isnan(m) );
  setMass (m);
  adjustPThetaPhi (m, Pt, theta, phi);
  setParam (0, Pt, true);
  setParam (1, theta, true);
  setParam (2, phi, true);
  setMParam (0, Pt);
  setMParam (1, theta);
  setMParam (2, phi);
  setError (0, DPt);
  setError (1, Dtheta);
  setError (2, Dphi);
  invalidateCache();
//   std::cout << "JetFitObject::JetFitObject: E = " << E << std::endl;
//   std::cout << "JetFitObject::JetFitObject: getParam(0) = " << getParam(0) << std::endl;
//   std::cout << "JetFitObject::JetFitObject: " << *this << std::endl;
//   std::cout << "mpar= " << mpar[0] << ", " << mpar[1] << ", " << mpar[2] << std::endl;
}

// destructor
LeptonFitObject::~LeptonFitObject() {}

LeptonFitObject *LeptonFitObject::copy() const {
  return new LeptonFitObject (*this);
}
    
LeptonFitObject& LeptonFitObject::assign (const BaseFitObject& source) {
  if (const LeptonFitObject *psource = dynamic_cast<const LeptonFitObject *>(&source)) {
    if (psource != this) *this = *psource;
  }
  else {
    assert (0);
  }
  return *this;
}

const char *LeptonFitObject::getParamName (int ilocal) const {
  switch (ilocal) {
    case 0: return "1/Pt";
    case 1: return "theta";
    case 2: return "phi";
  }
  return "undefined";
}

// needed for constructor!
bool LeptonFitObject::setParam (int ilocal, double par_, 
                                    bool measured_, bool fixed_) {
  assert (ilocal >= 0 && ilocal < 3);
  if (measured[ilocal] != measured_ || fixed[ilocal] != fixed_) invalidateCache();
  measured[ilocal] = measured_;
  fixed[ilocal] = fixed_;
  par[ilocal] = par_;
  return true;
}  

bool LeptonFitObject::setParam (int i, double par_ ) {
  invalidateCache();
//   if (i==0) {
//      std::cout << "setParam: par_ = " << par_ << endl;
//      std::cout << "setParam: par[0] = " << par[0] << endl;
//   }   
  bool result = (par_-par[i])*(par_-par[i]) > eps2*cov[i][i]; 
  switch (i) {
    // Pt: positive
    case 0: par[0] = (par_ >= 0) ? par_ : -par_;
            //std::cout << "setParam: par[0] = " << par[0] << endl;
            break;
    // theta: between 0 and pi
    case 1: par[1] = par_;
            break;          
    // phi: any value
    case 2: par[2] = par_;
            break;          
    default: std::cerr << "LeptonFitObject::setParam: Illegal i=" << i << std::endl;
  }
  return result;
} 
 
bool LeptonFitObject::updateParams (double p[], int idim) {
  invalidateCache();
  
  int iPt = getGlobalParNum(0);
  int ith = getGlobalParNum(1);
  int iph = getGlobalParNum(2);
  assert (iPt >= 0 && iPt  < idim);
  assert (ith >= 0 && ith < idim);
  assert (iph >= 0 && iph < idim);
  
  double pt  = p[iPt];
  double th = p[ith];
  double ph = p[iph];
//  assert( !isinf(e) );    assert( !isnan(e) );
//  assert( !isinf(th) );   assert( !isnan(th) );
//  assert( !isinf(ph) );   assert( !isnan(ph) );
  
  if (pt<0) {
    // cout << "JetFitObject::updateParams: mirrored E!\n";
    pt  = -pt;
    th = M_PI-th;
    ph = M_PI+ph;
  }
  // double massPlusEpsilon = 1/(mass*(1.0000001));
  // if (pt < massPlusEpsilon) pt = massPlusEpsilon; 
  
  bool result = ((pt -par[0])*(pt -par[0]) > eps2*cov[0][0]) ||
                ((th-par[1])*(th-par[1]) > eps2*cov[1][1]) ||
                ((ph-par[2])*(ph-par[2]) > eps2*cov[2][2]);
                
//   if (result)
//     cout << "JetFitObject::updateParams " << getName() 
//          << "  e: " <<  par[0] << "->" << e            
//          << "  th: " <<  par[1] << "->" << th            
//          << "  ph: " <<  par[2] << "->" << ph
//          << endl;           
  
  par[0] = pt;
  par[1] = th;
  par[2] = ph;
  p[iPt]  = par[0];         
  p[ith] = par[1];         
  p[iph] = par[2];         
  return result;
}  

// these depend on actual parametrisation!
double LeptonFitObject::getPx() const {
  if (!cachevalid) updateCache();
  return px;
}
double LeptonFitObject::getPy() const {
  if (!cachevalid) updateCache();
  return py;
}
double LeptonFitObject::getPz() const {
  if (!cachevalid) updateCache();
  return pz;
}
double LeptonFitObject::getE() const {
 if (!cachevalid) updateCache();
  return e;
}
double LeptonFitObject::getP() const {
  if (!cachevalid) updateCache();
  return p;
}
double LeptonFitObject::getP2() const {
  if (!cachevalid) updateCache();
  return p2;
}
double LeptonFitObject::getPt() const {
  return par[0];
}
  
double LeptonFitObject::getPt2() const {
  if (!cachevalid) updateCache();
  return pt*pt;
}

double LeptonFitObject::getDPx(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dpxdPt;
      //case 1: return dpxdtheta; 
    case 1: return 0;
    case 2: return dpxdphi;
  }
  return 0; 
}

double LeptonFitObject::getDPy(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dpydPt;
      // case 1: return dpydtheta;
    case 1: return 0;
    case 2: return dpydphi;
  }
  return 0; 
}

double LeptonFitObject::getDPz(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dpzdPt; 
    case 1: return dpzdtheta;
    case 2: return 0;
  }
  return 0; 
}

double LeptonFitObject::getDE(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR); 
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dEdPt;
    case 1: return dEdtheta;
    case 2: return 0;
  }
  return 0; 
}
 
void   LeptonFitObject::addToDerivatives (double der[], int idim, 
                                       double efact, double pxfact, 
                                       double pyfact, double pzfact) const {
  int i_Pt     = globalParNum[0];
  int i_theta = globalParNum[1];
  int i_phi   = globalParNum[2];
  assert (i_Pt     >= 0 && i_Pt     < idim);
  assert (i_theta >= 0 && i_theta < idim);
  assert (i_phi   >= 0 && i_phi   < idim);
  
  if (!cachevalid) updateCache();
  // for numerical accuracy, add up derivatives first,
  // then add them to global vector
  double der_Pt = 0;
  double der_theta = 0;
  double der_phi = 0;
  
  if (pxfact != 0) {
    der_Pt     += pxfact*dpxdPt; 
    //der_theta += pxfact*dpxdtheta;
    der_phi   += pxfact*dpxdphi;
  }
  if (pyfact != 0) {
    der_Pt     += pyfact*dpydPt; 
    //der_theta += pyfact*dpydtheta;
    der_phi   += pyfact*dpydphi;
  }
  if (pzfact != 0) {
    der_Pt     += pzfact*dpzdPt;
    der_theta += pzfact*dpzdtheta;
  }
   if (efact != 0) {
    der_Pt     += efact*dEdPt;
    der_theta += efact*dEdtheta;
   }
  
   der[i_Pt]     += der_Pt;
   der[i_theta] += der_theta;
   der[i_phi]   += der_phi;
}
    
void   LeptonFitObject::addTo2ndDerivatives (double der2[], int idim, 
                                          double efact, double pxfact, 
                                          double pyfact, double pzfact) const {
  int i_Pt  = globalParNum[0];
  int i_th = globalParNum[1];
  int i_ph= globalParNum[2];
  assert (i_Pt  >= 0 && i_Pt  < idim);
  assert (i_th >= 0 && i_th < idim);
  assert (i_ph >= 0 && i_ph < idim);
  
  if (!cachevalid) updateCache();
  // for numerical accuracy, add up derivatives first,
  // then add them to global vector
  double der_PtPt   = 0;
  double der_Ptth  = 0;
  double der_Ptph  = 0;
  double der_thth = 0;
  double der_thph = 0;
  double der_phph = 0;
  
    
  if (pxfact != 0) {
    der_PtPt += pxfact*(2*cphi/(pt*pt*pt)); 
    der_Ptph  += pxfact*(sphi/(pt*pt));
    der_phph -= pxfact*(cphi/pt); 
  }
  if (pyfact != 0) { 
    der_PtPt += pyfact*(2*sphi/(pt*pt*pt)); 
    der_Ptph  -= pyfact*(cphi/(pt*pt));
    der_phph -= pyfact*(sphi/pt);
  }
  if (pzfact != 0) { 
    der_PtPt += pyfact*(2*cottheta/(pt*pt*pt)); 
    der_Ptth  += pzfact/(pt*pt*stheta*stheta);
    der_thth += pzfact*(2*cottheta/(pt*stheta*stheta));
  } 

  if (efact != 0) {
    der_PtPt   += efact*(3/(pt*pt*pt*pt*e*stheta*stheta)-1/(pt*pt*pt*pt*pt*pt*e*e*e*stheta*stheta*stheta*stheta));
    der_Ptth  += efact*(2*cottheta/(pt*pt*pt*e*stheta*stheta)-cottheta/(pt*pt*pt*pt*pt*e*e*e*stheta*stheta*stheta*stheta));
    der_thth += efact*((1/(pt*pt*e*stheta*stheta*stheta*stheta))+(2*cottheta*cottheta/(pt*pt*e*stheta*stheta))-(cottheta*cottheta/(pt*pt*pt*pt*e*e*e*stheta*stheta*stheta*stheta)));
  }

  der2[idim*i_Pt+i_Pt]   += der_PtPt;
  der2[idim*i_Pt+i_th]  += der_Ptth;
  der2[idim*i_Pt+i_ph]  += der_Ptph;
  der2[idim*i_th+i_Pt]  += der_Ptth;
  der2[idim*i_th+i_th] += der_thth;
  der2[idim*i_th+i_ph] += der_thph;
  der2[idim*i_ph+i_Pt]  += der_Ptph;
  der2[idim*i_ph+i_th] += der_thph;
  der2[idim*i_ph+i_ph] += der_phph;
}
    
void   LeptonFitObject::addTo2ndDerivatives (double M[], int idim,  double lambda, double der[]) const {
  addTo2ndDerivatives (M, idim, lambda*der[0], lambda*der[1], lambda*der[2], lambda*der[3]);
}

void   LeptonFitObject::addTo1stDerivatives (double M[], int idim, double der[], int kglobal) const {
  assert (kglobal >= 0 && kglobal < idim);
  int i_Pt  = globalParNum[0];
  int i_th = globalParNum[1];
  int i_ph= globalParNum[2];
  assert (i_Pt  >= 0 && i_Pt  < idim);
  assert (i_th >= 0 && i_th < idim);
  assert (i_ph >= 0 && i_ph < idim);
  
  if (!cachevalid) updateCache();
  
  double dPt  = der[0]*dEdPt + der[1]*dpxdPt + der[2]*dpydPt + der[3]*dpzdPt;
  double dth = der[0]*dEdtheta                               + der[3]*dpzdtheta;
  double dph =               + der[1]*dpxdphi + der[2]*dpydphi;
  
  M[idim*kglobal + i_Pt]  += dPt;  
  M[idim*kglobal + i_th] += dth; 
  M[idim*kglobal + i_ph] += dph; 
  M[idim*i_Pt  + kglobal] += dPt;  
  M[idim*i_th + kglobal] += dth; 
  M[idim*i_ph + kglobal] += dph; 
}

         
void LeptonFitObject::initCov() {
  for (int i = 0; i < NPAR; ++i) {
    for (int j = 0; j < NPAR; ++j) {
      cov[i][j] = static_cast<double>(i == j);
    }
  }    
}

void LeptonFitObject::invalidateCache() const {
  cachevalid = false;
}

void LeptonFitObject::addToGlobalChi2DerVector (double *y, int idim, 
                                             double lambda, double der[]) const {
  int i_Pt  = globalParNum[0];
  int i_th = globalParNum[1];
  int i_ph = globalParNum[2];
  assert (i_Pt  >= 0 && i_Pt  < idim);
  assert (i_th >= 0 && i_th < idim);
  assert (i_ph >= 0 && i_ph < idim);
  
  if (!cachevalid) updateCache();
  
  y[i_Pt]  += lambda*(der[0]*dEdPt + der[1]*dpxdPt     + der[2]*dpydPt     + der[3]*dpzdPt);
  y[i_th] += lambda*(der[0]*dEdtheta                                       + der[3]*dpzdtheta);
  y[i_ph] += lambda*(       + der[1]*dpxdphi        + der[2]*dpydphi);
}

void LeptonFitObject::updateCache() const {
  // std::cout << "LeptonFitObject::updateCache" << std::endl;
  chi2 = calcChi2 ();
  
  double pt    = par[0];
  double theta = par[1];
  double phi   = par[2];

  ctheta = cos(theta);
  stheta = sin(theta);
  cphi   = cos(phi);
  sphi   = sin(phi);
  cottheta = ctheta/stheta;
 
  if(mass > 0){
  p = 1/(pt*stheta);
  p2 = std::abs((p*p));
  e2 = std::abs(p2+mass*mass);
  e = std::sqrt(e2);
  dpdPt = -1/((pt*pt)*stheta);
  assert(p!=0);
  }
  else{
    p = 1/(pt*stheta);
    p2 = std::abs((p*p));
    e2 = std::abs(p2);
    e = std::sqrt(e2);
    dpdPt = -1/((pt*pt)*stheta);
  }

  px = cphi/pt;
  py = sphi/pt;
  pz = cottheta/pt;
 
  dpxdPt = -cphi/(pt*pt);
  dpydPt = -sphi/(pt*pt);
  dpzdPt = -cottheta/(pt*pt);
  //dpxdtheta = pz*cphi;
  //dpydtheta = pz*sphi;
  dpxdphi = -sphi/pt;
  dpydphi = cphi/pt;
  dpzdtheta = -1/(pt*stheta*stheta);
  dEdPt = -1/(pt*pt*pt*e*stheta*stheta);
  dEdtheta = -cottheta/(pt*pt*e*stheta*stheta);
 
  cachevalid = true;
}

double LeptonFitObject::getError2 (double der[]) const {
  if (!cachevalid) updateCache();
  double cov4[4][4]; // covariance of E, px, py, px
  cov4[0][0] = dEdPt*(dEdPt*cov[0][0] + 2*dEdtheta*cov[0][1]) + dEdtheta*dEdtheta*cov[1][1];             // E, E
  cov4[0][1] = cov4[1][0] =                          // E, px 
    dpxdPt*(dEdPt*cov[0][0] + dEdtheta*cov[0][1]) + dpxdphi*(dEdPt*cov[2][0] + dEdtheta*cov[2][1]); 
  cov4[0][2] = cov4[2][0] =                          // E, py
    dpydPt*(dEdPt*cov[0][0] + dEdtheta*cov[0][1]) + dpydphi*(dEdPt*cov[2][0] + dEdtheta*cov[2][1]); 
  cov4[0][3] = cov4[3][0] =                          // E, pz
    dpzdPt*(dEdPt*cov[0][0] + dEdtheta*cov[0][1]) + dpzdtheta*(dEdPt*cov[1][0] + dEdtheta*cov[1][1]); 
  cov4[1][1] =                                       // px, px
    dpxdPt*(dpxdPt*cov[0][0] + 2*dpxdphi*cov[0][2]) + dpxdphi*dpxdphi*cov[2][2]; 
  cov4[1][2] = cov4[2][1] =                         // px, py
    dpydPt*(dpxdPt*cov[0][0] + dpxdphi*cov[0][2]) + dpydphi*(dpxdPt*cov[2][0] + dpxdphi*cov[2][2]); 
  cov4[1][3] = cov4[3][1] =                         // px, pz
    dpzdPt*(dpxdPt*cov[0][0] + dpxdphi*cov[0][2]) + dpzdtheta*(dpxdPt*cov[1][0] + dpxdphi*cov[2][1]); 
  cov4[2][2] =                                       // py, py
    dpydPt*(dpydPt*cov[0][0] + 2*dpydphi*cov[0][2]) + dpydphi*dpydphi*cov[2][2]; 
  cov4[2][3] = cov4[3][2] =                          // py, pz
    dpzdPt*(dpydPt*cov[0][0] + dpydphi*cov[0][2]) + dpzdtheta*(dpydPt*cov[1][0] + dpydphi*cov[2][1]); 
  cov4[3][3] =                                      // pz, pz
       dpzdPt*(dpzdPt*cov[0][0] + 2*dpzdtheta*cov[0][1]) + dpzdtheta*dpzdtheta*cov[1][1]; 
  return der[0]*(der[0]*cov4[0][0] + 2*der[1]*cov4[0][1] + 2*der[2]*cov4[0][2] + 2*der[3]*cov4[0][3])
                             + der[1]*(der[1]*cov4[1][1] + 2*der[2]*cov4[1][2] + 2*der[3]*cov4[1][3])
                                                   + der[2]*(der[2]*cov4[2][2] + 2*der[3]*cov4[2][3])
                                                                          + der[3]*der[3]*cov4[3][3];
}

double LeptonFitObject::getChi2 () const {
  if (!cachevalid) updateCache();
  return chi2;
}

bool LeptonFitObject::adjustPtThetaPhi (double& m, double &Pt, double& theta, double& phi) {
  bool result = false;
  
  if (Pt<0) {
    // cout << "LeptonFitObject::adjustEThetaPhi: mirrored E!\n";
    Pt  = -Pt;
    theta = M_PI-theta;
    phi = M_PI+phi;
    result = true;
  }
  if (theta < -M_PI || theta > M_PI) {
    while (theta < -M_PI) theta += 2*M_PI;
    while (theta >  M_PI) theta -= 2*M_PI;
    result = true;
  }
  
  if (theta<0) {
    // cout << "LeptonFitObject::adjustEThetaPhi: mirrored theta!\n";
    theta = -theta;
    phi = phi > 0 ? phi-M_PI : phi+M_PI;
    result = true;
  }
  else if (theta>M_PI) {
    // cout << "LeptonFitObject::adjustEThetaPhi: mirrored theta!\n";
    theta = 2*M_PI-theta;
    phi = phi > 0 ? phi-M_PI : phi+M_PI;
    result = true;
  }
  if (phi < -M_PI || phi > M_PI) {
    while (phi < -M_PI) phi += 2*M_PI;
    while (phi >  M_PI) phi -= 2*M_PI;
    result = true;
  }

  return result;
}

double LeptonFitObject::calcChi2 () const {
  if (!covinvvalid) calculateCovInv();
  if (!covinvvalid) return -1;
  
  double chi2 = 0;
  static double resid[3];
  //static bool chi2contr[3];
  resid[0] = isParamMeasured(0) && !isParamFixed(0) ? par[0]-mpar[0] : 0;
  resid[1] = isParamMeasured(1) && !isParamFixed(1) ? par[1]-mpar[1] : 0;
  resid[2] = isParamMeasured(2) && !isParamFixed(2) ? par[2]-mpar[2] : 0;
  if      (resid[2] >  M_PI) resid[2] -= 2*M_PI;
  else if (resid[2] < -M_PI) resid[2] += 2*M_PI;
    
  chi2 =     resid[0]*covinv[0][0]*resid[0] 
         + 2*resid[0]*covinv[0][1]*resid[1] 
         + 2*resid[0]*covinv[0][2]*resid[2] 
         +   resid[1]*covinv[1][1]*resid[1] 
         + 2*resid[1]*covinv[1][2]*resid[2] 
         +   resid[2]*covinv[2][2]*resid[2]; 
         
  return chi2;
}

