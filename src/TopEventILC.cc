////////////////////////////////////////////////////////////////
// Class TopEventILC
//
// Author: Benno List, Jenny Boehme
// Last update: $Date: 2008/09/26 09:58:11 $
//          by: $Author: boehmej $
// 
// Description: class to generate and fit top pair events at ILC
//               
////////////////////////////////////////////////////////////////
#ifdef MARLIN_USE_ROOT
#include "TopEventILC.h"

#include "JetFitObject.h"
#include "LeptonFitObject.h"
#include "NeutrinoFitObject.h"
#include "MassConstraint.h"
#include "SoftGaussMassConstraint.h"
#include "OPALFitterGSL.h"

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
//static  double mtop = 174;
//static  double gammatop = 1.4;
//static  double mW   = 80.4;
//static  double gammaW = 2.1;
//static  double mb   = 5.0;
//static  double mj   = 1.0;
//static  double Ecm = 500;

// constructor: 
TopEventILC::TopEventILC()
: leptonic (false), leptonasjet (false), debug (false),
  pxc (0, 1),
  pyc (0, 0, 1),
  pzc (0, 0, 0, 1),
  ec  (1, 0, 0, 0, 500), 
  sw1 (2.1/(2.*sqrt(0.805)), 80.4),  // Thesis Jenny p44
  sw2 (2.1/(2.*sqrt(0.805)), 80.4),
  sw (1.4/sqrt(0.805)),
  w1 (80.4),  
  w2 (80.4),
  w (0)
  {
  for (int i = 0; i < NFV; ++i) fv[i] = 0;
  for (int i = 0; i < NBFO; ++i) bfo[i] = bfosmear[i] = 0;
  pxc.setName ("px=0");
  pyc.setName ("py=0");
  pzc.setName ("pz=0");
  ec.setName  ("E=500");
  sw.setName ("top-equalmass-soft");
  sw1.setName ("w1-mass-soft");
  sw2.setName ("w2-mass-soft");
  w.setName ("top-equalmass");
  w1.setName ("w1-mass");
  w2.setName ("w2-mass");
}

//destructor: 
TopEventILC::~TopEventILC() {
  for (int i = 0; i < NFV; ++i) delete fv[i];
  for (int i = 0; i < NBFO; ++i) {
    delete bfo[i];
    delete bfosmear[i];
  }  
}

// Generate Breit-Wigner Random number
double TopEventILC::bwrandom (double r, double e0, double gamma, double emin, double emax) const {
  double a = atan (2.0*(emax - e0)/gamma);
  double b = atan (2.0*(emin - e0)/gamma);
  return e0 + 0.5*gamma*tan (r*(a - b) + b);
}


// generate four vectors
void TopEventILC::genEvent(){

    // reset all constraints
    pxc.resetFOList();
    pyc.resetFOList();
    pzc.resetFOList();
    ec.resetFOList();
    sw1.resetFOList();
    sw2.resetFOList();
    sw.resetFOList();
    w1.resetFOList();
    w2.resetFOList();
    w.resetFOList();
   
   
  // generate 4-vectors of top-decay:
  // 0: top-top-system -> top1 top2
  // 1: top 1 -> W1 b1
  // 2: top 2 -> W2 b2
  // 3: W1 -> j11 j12
  // 4: W2 -> j21 j22
  // 5: b1
  // 6: j11
  // 7: j12
  // 8: b2
  // 9: j21
  //10: j22 or neutrino
  
  double mtop = 174;
  double gammatop = 1.4;
  double mW   = 80.4;
  double gammaW = 2.1;
  double mb   = 5.0;
  double mj   = 1.0;
  double Ecm = 500;
      
  double rw[4];
  if (rnd == 0) rnd = new TRandom3();
  rnd->RndmArray (4, rw);
  
  FourVector *toppair = fv[0] = new FourVector (Ecm, 0, 0, 0);
  double mtop1 = bwrandom (rw[0], mtop, gammatop, mtop-3*gammatop, mtop+3*gammatop);
  double mtop2 = bwrandom (rw[1], mtop, gammatop, mtop-3*gammatop, mtop+3*gammatop);
  FourVector *top1 = fv[1] = new FourVector (mtop1, 0, 0, 0);
  FourVector *top2 = fv[2] = new FourVector (mtop2, 0, 0, 0);
  
  toppair->decayto (*top1, *top2);
  if (debug) {
    cout << "top 1: m=" << mtop1 << " = " << top1->getM() << endl;
    cout << "top 2: m=" << mtop2 << " = " << top2->getM() << endl;
  }  
  
  double mw1 = bwrandom (rw[2], mW, gammaW, mW-3*gammaW, mW+3*gammaW);
  double mw2 = bwrandom (rw[3], mW, gammaW, mW-3*gammaW, mW+3*gammaW);
  
  FourVector *W1 = fv[3] = new FourVector (mw1, 0, 0, 0);
  FourVector *W2 = fv[4] = new FourVector (mw2, 0, 0, 0);
  FourVector *b1 = fv[5] = new FourVector (mb, 0, 0, 0);
  FourVector *b2 = fv[8] = new FourVector (mb, 0, 0, 0);
  if (debug) {
    cout << "W 1: m=" << mw1 << " = " << W1->getM() << endl;
    cout << "W 2: m=" << mw2 << " = " << W2->getM() << endl;
  }  
  
  top1->decayto (*W1, *b1);
  top2->decayto (*W2, *b2);
  
  FourVector *j11 = fv[6]  = new FourVector (mj, 0, 0, 0);
  FourVector *j12 = fv[7]  = new FourVector (mj, 0, 0, 0);
  if (leptonic) mj = 0; // W2 decays to "massless" particles
  FourVector *j21 = fv[9]  = new FourVector (mj, 0, 0, 0);
  FourVector *j22 = fv[10] = new FourVector (mj, 0, 0, 0);
  
  W1->decayto (*j11, *j12);
  W2->decayto (*j21, *j22);
  
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
  
  for (int j = 0; j < 6; ++j) {
    int i = j+5;
    double E = fv[i]->getE();
    double theta = fv[i]->getTheta();
    double phi = fv[i]->getPhi();
    double ptinv = 1/(fv[i]->getPt());
    //double EError = (j==4 && leptonic) ? Eresolem*sqrt(E) : Eresolhad*sqrt(E);
    double EError = Eresolhad*sqrt(E);  // for jets
    //double ptinvError = ptinv*ptinv*sqrt(pow(sin(theta)*Eresolem*sqrt(E),2)+pow(E*cos(theta)*thetaResol,2));
    double ptinvError = sqrt(pow(2E-5,2) + pow(1E-3*ptinv/sin(theta),2));
    if (debug) {
      cout << "particle " << j << ": pt = " << 1./ptinv << ", sin(theta) = " << sin(theta) 
                               << ", ptinvError = " << ptinvError << endl;
      cout << "particle " << j << ": E = " << E << ", EError = " << EError << endl;
    }
    if (j==4 && leptonic && leptonasjet) {
        EError = Eresolem*sqrt(E);
    }
    
    static const char *names[] = {"b1", "j11", "j12", "b2", "j21", "j22"};
    // Create fit object with true quantities for later comparisons
    if (j < 4 || !leptonic || (j == 4 && leptonasjet)) {
      bfo[j] = new JetFitObject (E, theta, phi, EError, thetaResol, phiResol, 0);
      bfo[j]->setName (names[j]);
      if (debug) {
        cout << "true jet " << j << ": E = " << bfo[j]->getParam(0) << " +- " << bfo[j]->getError(0)
             << ", theta = " << bfo[j]->getParam(1) << " +- " << bfo[j]->getError(1)
             << ", phi = " << bfo[j]->getParam(2) << " +- " << bfo[j]->getError(2)
             << endl;
      }
    }  
    else if (j == 4 && leptonic && !leptonasjet) {
      bfo[4] = new LeptonFitObject (ptinv, theta, phi, ptinvError, thetaResolTrack, phiResol, 0.);
      if (debug) {
        cout << " true Lepton: E = " << bfo[4]->getE() << ", px = " << bfo[4]->getPx() << ", py = " << bfo[4]->getPy() 
             << ", pz = " << bfo[4]->getPz() << endl;
        cout << " true Lepton: ptinv = " << bfo[4]->getParam(0) << " +- " << bfo[4]->getError(0)
             << ", theta = " << bfo[4]->getParam(1) << " +- " << bfo[4]->getError(1)
             << ", phi = " << bfo[4]->getParam(2) << " +- " << bfo[4]->getError(2)
             << endl;
      }       
    }  
    else if (j == 5 && leptonic) {
      bfo[5] = new NeutrinoFitObject (E, theta, phi, 0.01, 0.0001, 0.00001);
      bfo[5]->setName ("n22");
      if (debug) {
        cout << "true Neutrino: E = " << bfo[5]->getE() << ", px = " << bfo[5]->getPx() << ", py = " << bfo[5]->getPy() 
             << ", pz = " << bfo[5]->getPz() << endl;
        cout << "true Neutrino " << j << ": E = " << bfo[5]->getParam(0) << " +- " << bfo[5]->getError(0)
             << ", theta = " << bfo[5]->getParam(1) << " +- " << bfo[5]->getError(1)
             << ", phi = " << bfo[5]->getParam(2) << " +- " << bfo[5]->getError(2)
             << endl;
      }       
    }  
    if (j == 4 && leptonic) bfo[4]->setName ("e22");
    
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
    
    
    if (j < 4 || !leptonic || (j == 4 && leptonasjet)) {
      bfosmear[j] = new JetFitObject (ESmear, thetaSmear, phiSmear, EError, thetaResol, phiResol, 0.);
      bfosmear[j]->setName (names[j]);
      bfostart[j] = new JetFitObject (ESmear, thetaSmear, phiSmear, EError, thetaResol, phiResol, 0.);
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
    else if (j == 4 && leptonic && !leptonasjet) {
      bfosmear[4] = new LeptonFitObject (ptinvSmear, thetaSmearTrack, phiSmearTrack, ptinvError, thetaResolTrack, phiResolTrack, 0.);
      bfostart[4] = new LeptonFitObject (ptinvSmear, thetaSmearTrack, phiSmearTrack, ptinvError, thetaResolTrack, phiResolTrack, 0.);
      Etot  += bfosmear[4]->getE();
      pxtot += bfosmear[4]->getPx();
      pytot += bfosmear[4]->getPy();
      pztot += bfosmear[4]->getPz(); 
      if (debug) {
        cout << "Lepton energy by hand, exact theta: e=sqrt(pow(pt/sintheta,2)+m*m) = " 
             << sqrt(pow(1./ptinvSmear/sin(theta),2)+mj*mj) << endl;
        cout << "Lepton energy by hand, smeared theta: e=sqrt(pow(pt/sintheta,2)+m*m) = " 
             << sqrt(pow(1./ptinvSmear/sin(thetaSmearTrack),2)+mj*mj) << endl;
        cout << "Lepton: E = " << bfosmear[4]->getE() << ", px = " << bfosmear[4]->getPx() << ", py = " << bfosmear[4]->getPy() 
             << ", pz = " << bfosmear[4]->getPz() << endl;
        cout << " smeared Lepton: ptinv = " << bfosmear[4]->getParam(0) << " +- " << bfosmear[4]->getError(0)
             << ", theta = " << bfosmear[4]->getParam(1) << " +- " << bfosmear[4]->getError(1)
             << ", phi = " << bfosmear[4]->getParam(2) << " +- " << bfosmear[4]->getError(2)
             << endl;
      }       
    }
    else if (j == 5 && leptonic) {
      double pxn = -pxtot;
      double pyn = -pytot;
      double pzn = -pztot;
      double pn = sqrt(pxn*pxn+pyn*pyn+pzn*pzn);
//      double en =  sqrt (pxn*pxn+pyn*pyn+pzn*pzn);
      double en =  Ecm - Etot;
      double ptn = sqrt(pxn*pxn+pyn*pyn);
      double theta = acos (pzn/pn); 
      double phi = atan2 (pyn, pxn);
//       if (en <= 0) {
//         cout << "WARNING: negative missing energy = " << en << ", setting to pn = " << pn << ", true pn = " << bfo[5]->getE() << endl;
//         en = pn;
//       }   
      if (debug) {
        cout << "Neutrino: en = " << en << ", theta = " << theta << ", phi = " << phi << endl;
        cout << "Neutrino: pxn = " << pxn << ", pyn = " << pyn << ", pzn = " << pzn << ", pn = " << pn << endl;
        cout << "Neutrino momenta by hand: px = " << ptn*cos(phi) << ", py = " << ptn*sin(phi) << ", pz = " << pn*cos(theta) << endl;
      }   
      bfosmear[5] = new NeutrinoFitObject (pn, theta, phi, 14, 0.32, 0.425);  // adjust such that "pull vs true" has width ~1
      bfostart[5] = new NeutrinoFitObject (pn, theta, phi, 14, 0.32, 0.425);  // adjust such that "pull vs true" has width ~1
    
      bfosmear[5]->setName ("n22");
      bfostart[5]->setName ("n22");
      if (debug) {
        cout << "Neutrino: E = " << bfosmear[5]->getE() << ", px = " << bfosmear[5]->getPx() << ", py = " << bfosmear[5]->getPy() 
             << ", pz = " << bfosmear[5]->getPz() << endl;
        cout << "smeared Neutrino " << j << ": E = " << bfosmear[5]->getParam(0) << " +- " << bfosmear[5]->getError(0)
             << ", theta = " << bfosmear[5]->getParam(1) << " +- " << bfosmear[5]->getError(1)
             << ", phi = " << bfosmear[5]->getParam(2) << " +- " << bfosmear[5]->getError(2)
             << endl;
      }       
    }
    if (j == 4 && leptonic) {
      bfosmear[4]->setName ("e22");
      bfostart[4]->setName ("e22");
    }  
    
    fvsmear[i] = new FourVector (bfosmear[j]->getE(), bfosmear[j]->getPx(), bfosmear[j]->getPy(), bfosmear[j]->getPz());
    
    pxc.addToFOList (*bfosmear[j]);
    pyc.addToFOList (*bfosmear[j]);
    pzc.addToFOList (*bfosmear[j]);
    ec.addToFOList (*bfosmear[j]);
    sw.addToFOList (*bfosmear[j], j<3?1:2);
    w.addToFOList (*bfosmear[j], j<3?1:2);
      
  }
  fvsmear[3] = new FourVector (*fvsmear[6]+*fvsmear[7]);
  fvsmear[4] = new FourVector (*fvsmear[9]+*fvsmear[10]);
  fvsmear[1] = new FourVector (*fvsmear[3]+*fvsmear[5]);
  fvsmear[2] = new FourVector (*fvsmear[4]+*fvsmear[8]);
  fvsmear[0] = new FourVector (*fvsmear[1]+*fvsmear[2]);
  
  sw1.addToFOList (*bfosmear[1]);
  sw1.addToFOList (*bfosmear[2]);
  w1.addToFOList (*bfosmear[1]);
  w1.addToFOList (*bfosmear[2]);
  sw2.addToFOList (*bfosmear[4]);
  sw2.addToFOList (*bfosmear[5]);
  w2.addToFOList (*bfosmear[4]);
  w2.addToFOList (*bfosmear[5]);
    
  if (debug) cout << "finished setting up constraints" << endl;
  
}

// fit it!
int TopEventILC::fitEvent (BaseFitter& fitter){
  
  if (debug) {
    for (int i = 0; i < 6; ++i) 
      cout << "true four-vector of jet " << i << ": " << *bfo[i] << endl;
    for (int i = 0; i < 6; ++i) 
      cout << "initial four-vector of jet " << i << ": " << *bfosmear[i] << endl;
  }
  
  // reset lists of constraints and fitobjects
  fitter.reset();
  
  if (debug) {
    cout << "TopEventILC::fitEvent: ==================================================\n";
    cout << "True vectors: \n";
    for (int i = 0; i<6; ++i) {
      cout << bfo[i]->getName() << ": " << *bfo[i] << endl;
    }
    cout << "Start vectors: \n";
    for (int i = 0; i<6; ++i) {
      cout << bfosmear[i]->getName() << ": " << *bfosmear[i] << endl;
    }
    cout << "Total: \n";
    cout << "gen:   " << *fv[0] << ", m=" << fv[0]->getM() << endl;
    cout << "smear: " << *fvsmear[0] << ", m=" << fvsmear[0]->getM() << endl;
    cout << "Top1: \n";
    cout << "gen:   " << *fv[1] << ", m=" << fv[1]->getM() << endl;
    cout << "smear: " << *fvsmear[1] << ", m=" << fvsmear[1]->getM() << endl;
    cout << "Top2: \n";
    cout << "gen:   " << *fv[2] << ", m=" << fv[2]->getM() << endl;
    cout << "smear: " << *fvsmear[2] << ", m=" << fvsmear[2]->getM() << endl;
    cout << "W1: \n";
    cout << "gen:   " << *fv[3] << ", m=" << fv[3]->getM() << endl;
    cout << "smear: " << *fvsmear[3] << ", m=" << fvsmear[3]->getM() << endl;
    cout << "W2: \n";
    cout << "gen:   " << *fv[4] << ", m=" << fv[4]->getM() << endl;
    cout << "smear: " << *fvsmear[4] << ", m=" << fvsmear[4]->getM() << endl;
  }
  
   
  for (int i = 0; i < 6; i++) {
    assert (bfosmear[i]);
    fitter.addFitObject (*bfosmear[i]);
  }

    
  fitter.addConstraint (pxc);
  fitter.addConstraint (pyc);
  fitter.addConstraint (pzc);
  fitter.addConstraint (ec);
  if ( softmasses && !dynamic_cast<OPALFitterGSL*>(&fitter) ) {
    fitter.addSoftConstraint (sw);
    fitter.addSoftConstraint (sw1);
    fitter.addSoftConstraint (sw2);
  } 
  else {  
    fitter.addConstraint (w);
    fitter.addConstraint (w1);
    fitter.addConstraint (w2);
  }
  
  double prob = fitter.fit();
  
  if (debug) {
    cout << "fit error = " << fitter.getError() << endl;
    cout << "fit probability = " << prob << endl;
    for (int i = 0; i < 6; ++i) 
      cout << "final four-vector of jet " << i << ": " << *bfosmear[i] << endl;
    cout << "Constraint ec: "  << ec.getValue()  << endl;
    cout << "Constraint pxc: " << pxc.getValue() << endl;
    cout << "Constraint pyc: " << pyc.getValue() << endl;
    cout << "Constraint pzc: " << pzc.getValue() << endl;
    cout << "Constraint sw:   " << sw.getValue() << ", top mass: " << sw.getMass() << endl;
    cout << "Constraint sw1:  " << sw1.getValue() << ", W mass: " << sw1.getMass() << endl;
    cout << "Constraint sw2:  " << sw2.getValue() << ", W mass: " << sw2.getMass() << endl;
    cout << "Constraint w:   " << w.getValue() << ", top mass: " << w.getMass() << endl;
    cout << "Constraint w1:  " << w1.getValue() << ", W mass: " << w1.getMass() << endl;
    cout << "Constraint w2:  " << w2.getValue() << ", W mass: " << w2.getMass() << endl;

    cout << "fit probability = " << prob << ", soft top mass: " << sw.getMass() << ", top mass: " << w.getMass() << ", dof: " << fitter.getDoF() 
         << ", iterations: " << fitter.getIterations() << endl;
  }       

  for (int j = 0; j < 6; ++j) {
    int i = j+5;
    fvfinal[i] = new FourVector (bfosmear[j]->getE(), bfosmear[j]->getPx(), bfosmear[j]->getPy(), bfosmear[j]->getPz());
  }
  
  fvfinal[3] = new FourVector (*fvfinal[6]+*fvfinal[7]);
  fvfinal[4] = new FourVector (*fvfinal[9]+*fvfinal[10]);
  fvfinal[1] = new FourVector (*fvfinal[3]+*fvfinal[5]);
  fvfinal[2] = new FourVector (*fvfinal[4]+*fvfinal[8]);
  fvfinal[0] = new FourVector (*fvfinal[1]+*fvfinal[2]);
  
  if (debug) {
    cout << "===============After Fiting ===================================\n";
    cout << "Final vectors: \n";
    for (int i = 0; i<6; ++i) {
      cout << bfosmear[i]->getName() << ": " << *bfosmear[i] << endl;
      if (i < 4 || !leptonic || (i == 4 && leptonasjet)) {
        cout << "fitted jet " << i << ": E = " << bfosmear[i]->getParam(0) << " +- " << bfosmear[i]->getError(0)
             << ", theta = " << bfosmear[i]->getParam(1) << " +- " << bfosmear[i]->getError(1)
             << ", phi = " << bfosmear[i]->getParam(2) << " +- " << bfosmear[i]->getError(2)
             << endl;
      }
      else if (i == 4 && leptonic && !leptonasjet) {
        cout << " fitted Lepton: ptinv = " << bfosmear[4]->getParam(0) << " +- " << bfosmear[4]->getError(0)
             << ", theta = " << bfosmear[4]->getParam(1) << " +- " << bfosmear[4]->getError(1)
             << ", phi = " << bfosmear[4]->getParam(2) << " +- " << bfosmear[4]->getError(2)
             << endl;
      }       
      else if (i == 5 && leptonic) {
        cout << "fitted Neutrino " << i << ": E = " << bfosmear[5]->getParam(0) << " +- " << bfosmear[5]->getError(0)
             << ", theta = " << bfosmear[5]->getParam(1) << " +- " << bfosmear[5]->getError(1)
             << ", phi = " << bfosmear[5]->getParam(2) << " +- " << bfosmear[5]->getError(2)
             << endl;
      }       
    }
    cout << "Total: \n";
    cout << "gen:   " << *fv[0] << ", m=" << fv[0]->getM() << endl;
    cout << "final: " << *fvfinal[0] << ", m=" << fvfinal[0]->getM() << endl;
    cout << "Top1: \n";
    cout << "gen:   " << *fv[1] << ", m=" << fv[1]->getM() << endl;
    cout << "final: " << *fvfinal[1] << ", m=" << fvfinal[1]->getM() << endl;
    cout << "Top2: \n";
    cout << "gen:   " << *fv[2] << ", m=" << fv[2]->getM() << endl;
    cout << "final: " << *fvfinal[2] << ", m=" << fvfinal[2]->getM() << endl;
    cout << "W1: \n";
    cout << "gen:   " << *fv[3] << ", m=" << fv[3]->getM() << endl;
    cout << "final: " << *fvfinal[3] << ", m=" << fvfinal[3]->getM() << endl;
    cout << "W2: \n";
    cout << "gen:   " << *fv[4] << ", m=" << fv[4]->getM() << endl;
    cout << "final: " << *fvfinal[4] << ", m=" << fvfinal[4]->getM() << endl;
    cout << "================================================\n";
  }
  

   
   return fitter.getError();


}
#endif // MARLIN_USE_ROOT
