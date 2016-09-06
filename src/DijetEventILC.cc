////////////////////////////////////////////////////////////////
// Class DijetEventILC
//
// Author: Jenny List
// Last update: $Date: 2015/11/18 09:58:11 $
//          by: $Author: boehmej $
// 
// Description: class to generate and fit di jet events at ILC
//               
////////////////////////////////////////////////////////////////
#ifdef MARLIN_USE_ROOT
#include "DijetEventILC.h"

#include "JetFitObject.h"
#include "LeptonFitObject.h"

#include <iostream>              // - cout
#include <cmath>            
#include <TRandom3.h>

static TRandom *rnd = 0;

using std::cout;
using std::endl;
using std::abs;
using std::sqrt;
using std::cos;
using std::sin;

// constructor: 
DijetEventILC::DijetEventILC()
: leptonic (false), leptonasjet (false), debug (false),
  pxc (0, 1, 0, 0, 0),
  pyc (0, 0, 1, 0, 0),
  pzc (0, 0, 0, 1, 0),
  ec  (1, 0, 0, 0, 500),
  mc( MassConstraint() )
  {
  for (int i = 0; i < NFV; ++i) fv[i] = 0;
  for (int i = 0; i < NBFO; ++i) bfo[i] = bfosmear[i] = 0;
  mc.setMass (500);
  pxc.setName ("px");
  pyc.setName ("py");
  pzc.setName ("pz");
  ec.setName  ("E");
  mc.setName  ("M=500");
}

//destructor: 
DijetEventILC::~DijetEventILC() {
  for (int i = 0; i < NFV; ++i) delete fv[i];
  for (int i = 0; i < NBFO; ++i) {
    delete bfo[i];
    delete bfosmear[i];
  }  
}


// generate four vectors
void DijetEventILC::genEvent(){

    // reset all constraints
    pxc.resetFOList();
    pyc.resetFOList();
    pzc.resetFOList();
    ec.resetFOList();
    mc.resetFOList();
   
   
  // generate 4-vectors of two jets, like step 0 of TopEvent:
  // 0: top-top-system -> top1 top2
  
  double mj   = 0;
  double Ecm = 500.;
      
  double rw[4];
  if (rnd == 0) rnd = new TRandom3();
  rnd->RndmArray (4, rw);
  
  FourVector *jetpair = fv[0] = new FourVector (Ecm, 0., 0., 0.);
  if (debug) {
    cout << "jetpair: m = " << jetpair->getM() << endl;
  }  
//  double mjet1 = bwrandom (rw[0], mtop, gammatop, mtop-3*gammatop, mtop+3*gammatop);
//  double mjet2 = bwrandom (rw[1], mtop, gammatop, mtop-3*gammatop, mtop+3*gammatop);
  // do something random later
  double mjet1 = mj;
  double mjet2 = mj;
  FourVector *jet1 = fv[1] = new FourVector (mjet1, 0, 0, 0);
  FourVector *jet2 = fv[2] = new FourVector (mjet2, 0, 0, 0);
  
  jetpair->decayto (*jet1, *jet2);
  if (debug) {
    cout << "jet 1: m=" << mjet1 << " = " << jet1->getM() << endl;
    cout << "jet 2: m=" << mjet2 << " = " << jet2->getM() << endl;
  }  
    
  double Eresolhad = 0.35;     // 35% / sqrt (E)
  double Eresolem = 0.10;     // 10% / sqrt (E)
  double thetaResol = 0.1;  // rad
  double phiResol = 0.1;    // rad
  double thetaResolTrack = 0.001;  // rad
  double phiResolTrack = 0.001;    // rad
  
  double Etot=0;
  double pxtot=0;
  double pytot=0;
  double pztot=0;
  
  for (int j = 0; j < 2; ++j) {
    int i = j+1;
    double E = fv[i]->getE();
    double theta = fv[i]->getTheta();
    double phi = fv[i]->getPhi();
    double ptinv = 1/(fv[i]->getPt());
    //double EError = (j==4 && leptonic) ? Eresolem*sqrt(E) : Eresolhad*sqrt(E);
    //double EError = Eresolhad*sqrt(E);  // for jets
    double EError = Eresolhad*sqrt(Ecm/2);  // use a fixed resolution, should be equivalent to jet energy for di-jet events! 
    //double ptinvError = ptinv*ptinv*sqrt(pow(sin(theta)*Eresolem*sqrt(E),2)+pow(E*cos(theta)*thetaResol,2));
    double ptinvError = sqrt(pow(2E-5,2) + pow(1E-3*ptinv/sin(theta),2));
    if (debug) {
      cout << "particle " << j << ": pt = " << 1./ptinv << ", sin(theta) = " << sin(theta) 
                               << ", ptinvError = " << ptinvError << endl;
      cout << "particle " << j << ": E = " << E << ", EError = " << EError << endl;
    }
    if (leptonic && leptonasjet) {
        EError = Eresolem*sqrt(E);
    }
    
    static const char *names[] = {"j1", "j2"};
    // Create fit object with true quantities for later comparisons
    if (!leptonic || leptonasjet) {
      bfo[j] = new JetFitObject (E, theta, phi, EError, thetaResol, phiResol, mj);
      bfo[j]->setName (names[j]);
      if (debug) {
        cout << "true jet " << j << ": E = " << bfo[j]->getParam(0) << " +- " << bfo[j]->getError(0)
             << ", theta = " << bfo[j]->getParam(1) << " +- " << bfo[j]->getError(1)
             << ", phi = " << bfo[j]->getParam(2) << " +- " << bfo[j]->getError(2)
             << endl;
      }
    }  
    else if (leptonic && !leptonasjet) {
      bfo[j] = new LeptonFitObject (ptinv, theta, phi, ptinvError, thetaResolTrack, phiResol, 0);
      bfo[j]->setName (names[j]);
      if (debug) {
        cout << " true Lepton: E = " << bfo[j]->getE() << ", px = " << bfo[j]->getPx() << ", py = " << bfo[j]->getPy() 
             << ", pz = " << bfo[j]->getPz() << endl;
        cout << " true Lepton: ptinv = " << bfo[j]->getParam(0) << " +- " << bfo[j]->getError(0)
             << ", theta = " << bfo[j]->getParam(1) << " +- " << bfo[j]->getError(1)
             << ", phi = " << bfo[j]->getParam(2) << " +- " << bfo[j]->getError(2)
             << endl;
      }       
    }  
    
    double randoms[3];
    if (rnd == 0) rnd = new TRandom3();
    for (int irnd = 0; irnd < 3; ++irnd) randoms[irnd] = rnd->Gaus();
    
    // Create fit object with smeared quantities as fit input
    double ESmear = E + EError*randoms[0];
    double ptinvSmear = ptinv + ptinvError*randoms[0];
    if (debug) {
      cout << "ptinvSmear = " << ptinvSmear << ", ptinv = " << ptinv << ", ptinvError = " << ptinvError 
           << ", randoms[0] = " << randoms[0] << endl;
    }       

    double thetaSmear = theta + thetaResol*randoms[1];
    double phiSmear = phi + phiResol*randoms[2];
    double thetaSmearTrack = theta + thetaResolTrack*randoms[1];
    double phiSmearTrack = phi + phiResolTrack*randoms[2];
    
    
    if (!leptonic || leptonasjet) {
      bfosmear[j] = new JetFitObject (ESmear, thetaSmear, phiSmear, EError, thetaResol, phiResol, mj);
      bfosmear[j]->setName (names[j]);
      bfostart[j] = new JetFitObject (ESmear, thetaSmear, phiSmear, EError, thetaResol, phiResol, mj);
      bfostart[j]->setName (names[j]);
      Etot  += bfosmear[j]->getE();
      pxtot += bfosmear[j]->getPx();
      pytot += bfosmear[j]->getPy();
      pztot += bfosmear[j]->getPz();
      if (debug) {
        cout << "smeared jet " << j << ": E = " << bfosmear[j]->getParam(0) << " +- " << bfosmear[j]->getError(0)
             << ", theta = " << bfosmear[j]->getParam(1) << " +- " << bfosmear[j]->getError(1)
             << ", phi = " << bfosmear[j]->getParam(2) << " +- " << bfosmear[j]->getError(2)
             << endl;
      }
    }
    else if (leptonic && !leptonasjet) {
      bfosmear[j] = new LeptonFitObject (ptinvSmear, thetaSmearTrack, phiSmearTrack, ptinvError, thetaResolTrack, phiResolTrack, 0.);
      bfosmear[j]->setName (names[j]);
      bfostart[j] = new LeptonFitObject (ptinvSmear, thetaSmearTrack, phiSmearTrack, ptinvError, thetaResolTrack, phiResolTrack, 0.);
      bfostart[j]->setName (names[j]);
      Etot  += bfosmear[j]->getE();
      pxtot += bfosmear[j]->getPx();
      pytot += bfosmear[j]->getPy();
      pztot += bfosmear[j]->getPz(); 
      if (debug) {
        cout << "Lepton energy by hand, exact theta: e=sqrt(pow(pt/sintheta,2)+m*m) = " 
             << sqrt(pow(1./ptinvSmear/sin(theta),2)+mj*mj) << endl;
        cout << "Lepton energy by hand, smeared theta: e=sqrt(pow(pt/sintheta,2)+m*m) = " 
             << sqrt(pow(1./ptinvSmear/sin(thetaSmearTrack),2)+mj*mj) << endl;
        cout << "Lepton: E = " << bfosmear[j]->getE() << ", px = " << bfosmear[j]->getPx() << ", py = " << bfosmear[j]->getPy() 
             << ", pz = " << bfosmear[j]->getPz() << endl;
        cout << " smeared Lepton: ptinv = " << bfosmear[j]->getParam(0) << " +- " << bfosmear[j]->getError(0)
             << ", theta = " << bfosmear[j]->getParam(1) << " +- " << bfosmear[j]->getError(1)
             << ", phi = " << bfosmear[j]->getParam(2) << " +- " << bfosmear[j]->getError(2)
             << endl;
      }       
    }
    
    fvsmear[i] = new FourVector (bfosmear[j]->getE(), bfosmear[j]->getPx(), bfosmear[j]->getPy(), bfosmear[j]->getPz());
    if (debug) {
      cout << "jet " << i << ": m = " << fvsmear[i]->getM() << endl;
    }  
    
    
    pxc.addToFOList (*bfosmear[j]);
    pyc.addToFOList (*bfosmear[j]);
    pzc.addToFOList (*bfosmear[j]);
    ec.addToFOList (*bfosmear[j]);
    mc.addToFOList (*bfosmear[j]);
      
  }
  fvsmear[0] = new FourVector (*fvsmear[1]+*fvsmear[2]);
  if (debug) {
    cout << "jet 0: m = " << fvsmear[0]->getM() << endl;
  }  
    
}

// fit it!
int DijetEventILC::fitEvent (BaseFitter& fitter){
  
  if (debug) {
    for (int i = 0; i < 2; ++i) 
      cout << "true four-vector of jet " << i << ": " << *bfo[i] << endl;
    for (int i = 0; i < 2; ++i) 
      cout << "initial four-vector of jet " << i << ": " << *bfosmear[i] << endl;
  }
  
  // reset lists of constraints and fitobjects
  fitter.reset();
  
  if (debug) {
    cout << "DijetEventILC::fitEvent: ==================================================\n";
    cout << "True vectors: \n";
    for (int i = 0; i<2; ++i) {
      cout << bfo[i]->getName() << ": " << *bfo[i] << endl;
    }
    cout << "Start vectors: \n";
    for (int i = 0; i<2; ++i) {
      cout << bfosmear[i]->getName() << ": " << *bfosmear[i] << endl;
    }
    cout << "Total: \n";
    cout << "gen:   " << *fv[0] << ", m=" << fv[0]->getM() << endl;
    cout << "smear: " << *fvsmear[0] << ", m=" << fvsmear[0]->getM() << endl;
    cout << "Jet1: \n";
    cout << "gen:   " << *fv[1] << ", m=" << fv[1]->getM() << endl;
    cout << "smear: " << *fvsmear[1] << ", m=" << fvsmear[1]->getM() << endl;
    cout << "Jet2: \n";
    cout << "gen:   " << *fv[2] << ", m=" << fv[2]->getM() << endl;
    cout << "smear: " << *fvsmear[2] << ", m=" << fvsmear[2]->getM() << endl;
  }
  
   
  for (int i = 0; i < 2; i++) {
    assert (bfosmear[i]);
    fitter.addFitObject (*bfosmear[i]);
  }

    
  fitter.addConstraint (pxc);
  fitter.addConstraint (pyc);
  fitter.addConstraint (pzc);
  fitter.addConstraint (ec);
  //fitter.addConstraint (mc);
  
  double prob = fitter.fit();
  
  if (debug) {
    cout << "fit probability = " << prob << endl;
    for (int i = 0; i < 2; ++i) 
      cout << "final four-vector of jet " << i << ": " << *bfosmear[i] << endl;
    cout << "Constraint ec: "  << ec.getValue()  << endl;
    cout << "Constraint pxc: " << pxc.getValue() << endl;
    cout << "Constraint pyc: " << pyc.getValue() << endl;
    cout << "Constraint pzc: " << pzc.getValue() << endl;
    cout << "Constraint mc: " << mc.getValue() << endl;

    cout << "fit probability = " << prob << ", dof: " << fitter.getDoF() 
         << ", iterations: " << fitter.getIterations() << endl;
  }       

  for (int j = 0; j < 2; ++j) {
    int i = j+1;
    fvfinal[i] = new FourVector (bfosmear[j]->getE(), bfosmear[j]->getPx(), bfosmear[j]->getPy(), bfosmear[j]->getPz());
  }
  
  fvfinal[0] = new FourVector (*fvfinal[1]+*fvfinal[2]);
  
  if (debug) {
    cout << "===============After Fiting ===================================\n";
    cout << "Final vectors: \n";
    for (int i = 0; i<2; ++i) {
      cout << bfosmear[i]->getName() << ": " << *bfosmear[i] << endl;
      if (!leptonic || leptonasjet) {
        cout << "fitted jet " << i << ": E = " << bfosmear[i]->getParam(0) << " +- " << bfosmear[i]->getError(0)
             << ", theta = " << bfosmear[i]->getParam(1) << " +- " << bfosmear[i]->getError(1)
             << ", phi = " << bfosmear[i]->getParam(2) << " +- " << bfosmear[i]->getError(2)
             << endl;
      }
      else if (leptonic && !leptonasjet) {
        cout << " fitted Lepton: ptinv = " << bfosmear[i]->getParam(0) << " +- " << bfosmear[i]->getError(0)
             << ", theta = " << bfosmear[i]->getParam(1) << " +- " << bfosmear[i]->getError(1)
             << ", phi = " << bfosmear[i]->getParam(2) << " +- " << bfosmear[i]->getError(2)
             << endl;
      }       
    }
    cout << "Total: \n";
    cout << "gen:   " << *fv[0] << ", m=" << fv[0]->getM() << endl;
    cout << "final: " << *fvfinal[0] << ", m=" << fvfinal[0]->getM() << endl;
    cout << "Jet1: \n";
    cout << "gen:   " << *fv[1] << ", m=" << fv[1]->getM() << endl;
    cout << "final: " << *fvfinal[1] << ", m=" << fvfinal[1]->getM() << endl;
    cout << "Jet2: \n";
    cout << "gen:   " << *fv[2] << ", m=" << fv[2]->getM() << endl;
    cout << "final: " << *fvfinal[2] << ", m=" << fvfinal[2]->getM() << endl;
    cout << "================================================\n";
  }
  

   
   return fitter.getError();


}
#endif // MARLIN_USE_ROOT
