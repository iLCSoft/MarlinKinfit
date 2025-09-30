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
#include "MomentumConstraint.h"
#include "MassConstraint.h"

#include <array>

class DijetEventILC : public BaseEvent {
  public: 
    DijetEventILC();
    DijetEventILC(const DijetEventILC&) = delete;
    DijetEventILC& operator=(const DijetEventILC&) = delete;
    DijetEventILC(DijetEventILC&&) = delete;
    DijetEventILC& operator=(DijetEventILC&&) = delete;

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
    std::array<FourVector*, NFV> fv{};
    std::array<FourVector*, NFV> fvsmear{};
    std::array<FourVector*, NFV> fvfinal{};
    std::array<ParticleFitObject*, NBFO> bfo{};
    std::array<ParticleFitObject*, NBFO> bfostart{};
    std::array<ParticleFitObject*, NBFO> bfosmear{};
    
    MomentumConstraint pxc;
    MomentumConstraint pyc;
    MomentumConstraint pzc;
    MomentumConstraint ec;    
    MassConstraint mc;
    

};


#endif // __DIJETEVENTILC_H

#endif // MARLIN_USE_ROOT
