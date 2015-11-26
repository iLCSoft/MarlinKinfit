////////////////////////////////////////////////////////////////
// Class TopEventILC
//
// Author: Benno List, Jenny Boehme
// Last update: $Date: 2008/02/12 16:43:26 $
//          by: $Author: blist $
// 
// Description: class to generate and fit top pair events at ILC
//               
////////////////////////////////////////////////////////////////
#ifdef MARLIN_USE_ROOT

#ifndef __TOPEVENTILC_H
#define __TOPEVENTILC_H

#include "BaseEvent.h"
#include "JetFitObject.h"
#include "MomentumConstraint.h"
// #include "PxConstraint.h"
// #include "PyConstraint.h"
#include "MassConstraint.h"
#include "SoftGaussMassConstraint.h"

class TopEventILC : public BaseEvent {
  public: 
    TopEventILC();
    virtual ~TopEventILC();
    virtual void genEvent();
    virtual int fitEvent (BaseFitter& fitter);

    double bwrandom (double r, double e0, double gamma, double emin, double emax) const;
    
    MomentumConstraint& getPxConstraint() {return pxc;};
    MomentumConstraint& getPyConstraint() {return pyc;};
    MomentumConstraint& getPzConstraint() {return pzc;};
    MomentumConstraint& getEConstraint()  {return ec;};
    //MassConstraint& getW1Constraint() {return w1;};
    //MassConstraint& getW2Constraint() {return w2;};
    //MassConstraint& getTopConstraint() {return w;};
    SoftGaussMassConstraint& getW1Constraint() {return w1;};
    SoftGaussMassConstraint& getW2Constraint() {return w2;};
    SoftGaussMassConstraint& getTopConstraint() {return w;};
    
    double getW1Mass()  {return w1.getMass();};
    double getW2Mass()  {return w2.getMass();};
    double getTopMass(int flag)  {return w.getMass(flag);};
    double getTop1Mass()  {return fvsmear[1]->getM();};
    double getTop2Mass()  {return fvsmear[2]->getM();};
    
    void setDebug (bool _debug) {debug = _debug;};
    
    ParticleFitObject* getTrueFitObject (int i) {return bfo[i];};
    ParticleFitObject* getStartFitObject (int i) {return bfostart[i];};
    ParticleFitObject* getFittedFitObject (int i) {return bfosmear[i];};
    FourVector* getTrueFourVector (int i) {return fv[i];};
    
    bool leptonic, leptonasjet, debug;
    
  protected:
  
    enum {NFV = 11, NBFO = 6};
    FourVector *fv[NFV];
    FourVector *fvsmear[NFV];
    FourVector *fvfinal[NFV];
    ParticleFitObject *bfo[NBFO];
    ParticleFitObject *bfostart[NBFO];
    ParticleFitObject *bfosmear[NBFO];
    
    MomentumConstraint pxc;
    MomentumConstraint pyc;
    MomentumConstraint pzc;
    MomentumConstraint ec;
    //MassConstraint w1;
    //MassConstraint w2;
    //MassConstraint w;
    SoftGaussMassConstraint w1;
    SoftGaussMassConstraint w2;
    SoftGaussMassConstraint w;
    
    

};


#endif // __TOPEVENTILC_H

#endif // MARLIN_USE_ROOT
