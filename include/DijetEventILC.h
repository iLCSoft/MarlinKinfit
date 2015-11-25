////////////////////////////////////////////////////////////////
// Class DijetEventILC
//
// Author: Jenny Boehme
// Last update: $Date: 2015/11/18 16:43:26 $
//          by: $Author: boehmej $
// 
// Description: class to generate and fit di jet events at ILC
//               
////////////////////////////////////////////////////////////////
#ifdef MARLIN_USE_ROOT

#ifndef __DIJETEVENTILC_H
#define __DIJETEVENTILC_H

#include "BaseEvent.h"
#include "JetFitObject.h"
#include "MomentumConstraint.h"
#include "MassConstraint.h"

class DijetEventILC : public BaseEvent {
  public: 
    DijetEventILC();
    virtual ~DijetEventILC();
    virtual void genEvent();
    virtual int fitEvent (BaseFitter& fitter);
    
    MomentumConstraint& getPxConstraint() {return pxc;};
    MomentumConstraint& getPyConstraint() {return pyc;};
    MomentumConstraint& getPzConstraint() {return pzc;};
    MomentumConstraint& getEConstraint()  {return ec;};
    MassConstraint& getMassConstraint() {return mc;};
        
    void setDebug (bool _debug) {debug = _debug;};
    
    ParticleFitObject* getTrueFitObject (int i) {return bfo[i];};
    ParticleFitObject* getStartFitObject (int i) {return bfostart[i];};
    ParticleFitObject* getFittedFitObject (int i) {return bfosmear[i];};
    FourVector* getTrueFourVector (int i) {return fv[i];};
    
    bool leptonic, leptonasjet, debug;
    
  protected:
  
    enum {NFV = 3, NBFO = 2};
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
    MassConstraint mc;
    

};


#endif // __DIJETEVENTILC_H

#endif // MARLIN_USE_ROOT
