////////////////////////////////////////////////////////////////
// Class K0Event
//
// Author: Benno List, Jenny Boehme
// Last update: $Date: 2005/01/12 10:11:45 $
//          by: $Author: blist $
// 
// Description: class to generate and fit K0S decays
//               
////////////////////////////////////////////////////////////////
#include "K0Event.h"

#include "JetFitObject.h"
#include "MassConstraint.h"
//#include "cernlib.h"
#include "NeutralParticleTrack.h"
#include "ChargedParticleTrack.h"
#include "ThreeVector.h"
#include "FourVector.h"
#include "VertexConstraint.h"
#include "VertexFitObject.h"
#include "TrackMomentumConstraint.h"

#include <iostream>              // - cout
using std::cout;
using std::endl;
#include <cmath>  
using std::pow;            
using std::atan;            
using std::exp;            
using std::sin;            
using std::cos;            
using std::log;            
using std::abs;     
       
#include <TRandom3.h>

static TRandom *rnd = 0;

// constructor: 
K0Event::K0Event()
 {
  for (int i = 0; i < NFV; ++i) fv[i] = 0;
  for (int i = 0; i < NTFO; ++i) gentrack[i] = rectrack[i] = 0;
  for (int i = 0; i < NVER; ++i) genvert[i] = 0;
}

//destructor: 
K0Event::~K0Event() {
  for (int i = 0; i < NFV; ++i) delete fv[i];
  for (int i = 0; i < NTFO; ++i) {
    delete gentrack[i];
    delete rectrack[i];
  }  
  for (int i = 0; i < NVER; ++i) delete genvert[i];
}


// generate four vectors and tracks
void K0Event::genEvent(){
  /// generate primary vertex
  double randoms[4];
  if (rnd == 0) rnd = new TRandom3();
  for (int irnd = 0; irnd < 4; ++irnd) randoms[irnd] = rnd->Gaus();

  ThreeVector *primvert = genvert[0] = new ThreeVector (-0.4+0.015*randoms[0],
                                                         0.2+0.004*randoms[1],
                                                         1.0+10.5*randoms[2]);
//   ThreeVector *primvert = genvert[0] = new ThreeVector (0.015*randoms[0],
//                                                         0.004*randoms[1],
//                                                         10.5*randoms[2]);
//   ThreeVector *primvert = genvert[0] = new ThreeVector (0, 0, 0);
  /// generate K0S momentum vector and decay length
  if (rnd == 0) rnd = new TRandom3();
  rnd->RndmArray (4, randoms);
  double phi = 6.2831853*randoms[0];
  double ptmin = 1.0;
  double ptmax = 10.;
  int n = 4;           // pt^-4 
  double pt = pow (randoms[1]*(pow (ptmin, -n+1) - pow (ptmax, -n+1)) + pow (ptmax, -n+1), -1./(n-1));
  double etamin = -1.5;
  double etamax = 1.5;
  double eta = randoms[1]*(etamax-etamin) + etamin;
  double theta = 2*atan (exp (-eta));
  
//   cout << "k0s phi = " << phi << ", pt = " << pt << ", eta = " << eta << ", theta = " << theta << endl;
  
  ThreeVector k0smom (pt*cos(phi), pt*sin(phi), pt*cos(theta)/sin(theta)); 
  FourVector *k0s = fv[0] = new FourVector (k0smom, 0.497648);
  
//   cout << "k0s phi = " << k0s->getPhi() << ", pt = " << k0s->getPt() << ", eta = " 
//                << k0s->getEta() << ", theta = " << k0s->getTheta() << endl;
  

  FourVector *pi1 = fv[1] = new FourVector (0.13957, 0, 0, 0);
  FourVector *pi2 = fv[2] = new FourVector (0.13957, 0, 0, 0);

  k0s->decayto (*pi1, *pi2);
  
//   cout << "k0s 4-vector:         " << *k0s << endl;
//   cout << "pi1+p2 4-vector:      " << *pi1+*pi2 << endl;
//   cout << "pi1 4-vector:         " << *pi1 << endl;
//   cout << "pi2 4-vector:         " << *pi2 << endl << endl;
  
  
  double slen = -sin(theta)*log (randoms[3])*2.6842*k0s->getBetaGamma().getMag();
//  double slen = -sin(theta)*log (0.5)*2.6842*k0s->getBetaGamma().getMag();
  
//   cout << "slen = : " << slen << endl << endl;
  
  NeutralParticleTrack *k0strack = new NeutralParticleTrack ("K0 gen", *primvert, k0smom, 0.497648);
  gentrack[0] = k0strack;
  
  ThreeVector *decvert = genvert[1] = new ThreeVector;
  double sstop = k0strack->getArcLength(0) + slen;
  // k0strack->setArcLength(1, sstop);
  k0strack->getTrajectoryPointEx(sstop, *decvert);
  
  ChargedParticleTrack *pi1track;
  ChargedParticleTrack *pi2track;
  gentrack[1] = pi1track = new ChargedParticleTrack ("pi1 gen", *decvert, pi1->getThreeVector(), 0.13957, +1);
  gentrack[2] = pi2track = new ChargedParticleTrack ("pi2 gen", *decvert, pi2->getThreeVector(), 0.13957, -1);
  
//   cout << "k0s track: " << *gentrack[0] << endl;
//   cout << "pi1 track: " << *gentrack[1] << endl;
//   cout << "pi2 track: " << *gentrack[2] << endl << endl;
//   
//   cout << "primary vertex: " << *genvert[0] << endl;
//   cout << "k0s vertex:     " << gentrack[0]->getVertex(0) << endl << endl;
//   
//   cout << "decay vertex:   " << *genvert[1] << endl;
//   cout << "pi1 vertex:     " << gentrack[1]->getVertex(0) << endl;
//   cout << "pi2 vertex:     " << gentrack[2]->getVertex(0) << endl << endl;
//   
//   cout << "k0s 4-vector:         " << *k0s << endl;
//   cout << "k0s track 4-vector:   " << gentrack[0]->getMomentum(0) << endl;
//   cout << "pi1 4-vector:         " << *pi1 << endl;
//   cout << "pi1 track 4-vector:   " << gentrack[1]->getMomentum(0) << endl;
//   cout << "pi2 4-vector:         " << *pi2 << endl;
//   cout << "pi2 track 4-vector:   " << gentrack[2]->getMomentum(0) << endl << endl;
  
  rectrack[1] = createSmearedChargedTrack ("pi1 rec", *pi1track);
  rectrack[2] = createSmearedChargedTrack ("pi2 rec", *pi2track);
  smtrack[1] = new ChargedParticleTrack (* (ChargedParticleTrack*)rectrack[1]);
  smtrack[2] = new ChargedParticleTrack (* (ChargedParticleTrack*)rectrack[2]);
  smtrack[1]->setName ("pi1 sm");
  smtrack[2]->setName ("pi2 sm");
  
  ThreeVector nominalvertex;
  ThreeVector zeromomentum;
  rectrack[0] = new NeutralParticleTrack ("K0 rec", nominalvertex, zeromomentum, 0.497648);
  rectrack[0]->fixVertexParam(0);
  rectrack[0]->fixVertexParam(1, false);
  rectrack[0]->fixParam(5, true);
  
//   cout << "After smearing:\n";
//   cout << "k0s track 4-vector:   " << rectrack[0]->getMomentum(0) << endl;
//   cout << "pi1 track 4-vector:   " << rectrack[1]->getMomentum(0) << endl;
//   cout << "pi2 track 4-vector:   " << rectrack[2]->getMomentum(0) << endl;
  
 
}

ChargedParticleTrack *K0Event::createSmearedChargedTrack (const char *name, const ChargedParticleTrack& in) {
  double kappa = in.getParam(0)*ChargedParticleTrack::parfact[0];
  double phi0  = in.getParam(1)*ChargedParticleTrack::parfact[1];
  double theta = in.getParam(2)*ChargedParticleTrack::parfact[2];
  double dca   = in.getParam(3)*ChargedParticleTrack::parfact[3];
  double z0    = in.getParam(4)*ChargedParticleTrack::parfact[4];

  // pt=cBq/kappa => dpt = cBq/kappa^2 dkappa, dpt/pt = dkappa/kappa = pt * dkappa/cBq
  // => dkappa = dpt/pt^2*cBq = 0.01/GeV*cBq
  double dkappa = 0.01*abs(0.002998*TrackFitObject::bfield);
  double dphi0  = 0.0002;            // 0.2mrad = 100um at 50cm
  double dtheta = 0.08;              // = 2cm / 25cm
  double ddca   = 0.05;              // 0.5mm
  double dz0    = 2.;                // 2cm
  
  double randoms[5];
  if (rnd == 0) rnd = new TRandom3();
  for (int irnd = 0; irnd < 5; ++irnd) randoms[irnd] = rnd->Gaus();
  
  kappa += dkappa*randoms[0];
  phi0  += dphi0 *randoms[1];
  theta += dtheta*randoms[2];
  dca   += ddca  *randoms[3];
  z0    += dz0   *randoms[4];
 
  ChargedParticleTrack *result = new ChargedParticleTrack (name, kappa, phi0, theta, dca, z0, 
                                                           in.getMass(), in.getCharge(), 0);

  result->setError (0, dkappa/ChargedParticleTrack::parfact[0]);
  result->setError (1, dphi0/ChargedParticleTrack::parfact[1]);
  result->setError (2, dtheta/ChargedParticleTrack::parfact[2]);
  result->setError (3, ddca/ChargedParticleTrack::parfact[3]);
  result->setError (4, dz0/ChargedParticleTrack::parfact[4]);
  
  return result;
}

// fit it!
int K0Event::fitEvent (BaseFitter& fitter){
  
  fitter.addFitObject (rectrack[0]);
  fitter.addFitObject (rectrack[1]);
  fitter.addFitObject (rectrack[2]);

  VertexFitObject *k0decvertex = new VertexFitObject ("K0dec", 0, 0, 0);
  fitter.addFitObject (k0decvertex);
  
  rectrack[0]->fixVertexParam (0, true);
  k0decvertex->addTrack (rectrack[0], true, false);
  k0decvertex->addTrack (rectrack[1], false, true);
  k0decvertex->addTrack (rectrack[2], false, true);
  
  k0decvertex->addConstraints (fitter, VertexFitObject::VXYZ|VertexFitObject::PXYZ);
  
  k0decvertex->initForFit();
  
    
  double prob = fitter.fit();
  
  cout << "After fitting:\n";
  cout << "k0s track:            " << *rectrack[0] << endl;
  cout << "pi1 track:            " << *rectrack[1] << endl;
  cout << "pi2 track:            " << *rectrack[2] << endl;
  cout << "k0s track 4-vector:   " << rectrack[0]->getMomentum(0) << endl;
  cout << "pi1 track 4-vector:   " << rectrack[1]->getMomentum(0) << endl;
  cout << "pi2 track 4-vector:   " << rectrack[2]->getMomentum(0) << endl;
  cout << "k0s track vertex:     " << rectrack[0]->getVertex(0) << endl;
  cout << "pi1 track vertex:     " << rectrack[1]->getVertex(0) << endl;
  cout << "pi2 track vertex:     " << rectrack[2]->getVertex(0) << endl;
  cout << "fit probability = " << prob << endl << endl;
   
  return fitter.getError();
}
