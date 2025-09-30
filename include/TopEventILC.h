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
#include "MassConstraint.h"
#include "SoftGaussMassConstraint.h"

#include <array>

class TopEventILC : public BaseEvent {
  private:
    constexpr static int NFV = 11;
    constexpr static int  NBFO = 6;

  public: 
    TopEventILC();
    virtual ~TopEventILC();
    virtual void genEvent();
    virtual int fitEvent (BaseFitter& fitter);

    double bwrandom (double r, double e0, double gamma, double emin, double emax) const;
    
    BaseConstraint*  getPxConstraint() {return &pxc;};
    BaseConstraint*  getPyConstraint() {return &pyc;};
    BaseConstraint*  getPzConstraint() {return &pzc;};
    BaseConstraint*  getEConstraint()  {return &ec;};
    BaseConstraint* getW1Constraint() {return softmasses ? (BaseConstraint*) &sw1 : (BaseConstraint*) &w1;};
    BaseConstraint* getW2Constraint() {return softmasses ? (BaseConstraint*) &sw2 : (BaseConstraint*) &w2;};
    BaseConstraint* getTopConstraint() {return softmasses ? (BaseConstraint*) &sw : (BaseConstraint*) &w;};
    
    double getW1Mass()  {return softmasses ? w1.getMass() : sw1.getMass();};
    double getW2Mass()  {return softmasses ? w2.getMass() : sw2.getMass();};
    double getTopMass(int flag)  {return softmasses ? w.getMass(flag) : sw.getMass();};
    double getTop1Mass()  {return fvsmear[1]->getM();};
    double getTop2Mass()  {return fvsmear[2]->getM();};
    
    void setDebug (bool _debug) {debug = _debug;};
    
    ParticleFitObject* getTrueFitObject (int i) {return bfo[i];};
    ParticleFitObject* getStartFitObject (int i) {return bfostart[i];};
    ParticleFitObject* getFittedFitObject (int i) {return bfosmear[i];};
    FourVector* getTrueFourVector (int i) {return fv[i];};
    
    bool softmasses{false}, leptonic{false}, leptonasjet{false}, debug{false};
    
  protected:
  
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
    MassConstraint w1;
    MassConstraint w2;
    MassConstraint w;
    SoftGaussMassConstraint sw1;
    SoftGaussMassConstraint sw2;
    SoftGaussMassConstraint sw;
    
    

};


#endif // __TOPEVENTILC_H

#endif // MARLIN_USE_ROOT
