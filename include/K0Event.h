////////////////////////////////////////////////////////////////
// Class K0Event
//
// Author: Benno List, Jenny Boehme
// Last update: $Date: 2005/01/12 10:11:45 $
//          by: $Author: blist $
// 
// Description: class to generate and fit K0S events
//               
////////////////////////////////////////////////////////////////
#ifndef __K0Event_H
#define __K0Event_H

#include "BaseEvent.h"

class TrackFitObject;
class ThreeVector;
class ChargedParticleTrack;
class NeutralParticleTrack;

class K0Event : public BaseEvent {
  public: 
    K0Event();
    virtual ~K0Event();
    virtual void genEvent();
    virtual int fitEvent (BaseFitter& fitter);

    TrackFitObject* getTrueFitObject (int i) {return gentrack[i];};
    TrackFitObject* getSmearedFitObject (int i) {return rectrack[i];};
    ChargedParticleTrack *createSmearedChargedTrack (const char *name, const ChargedParticleTrack& in);
    
  public:
    enum {NFV = 3, NTFO = 3, NVER=2};
    FourVector *fv[NFV];
    TrackFitObject *gentrack[NTFO];
    TrackFitObject *smtrack[NTFO];
    TrackFitObject *rectrack[NTFO];
    ThreeVector *genvert[NVER];

};


#endif // __K0Event_H
