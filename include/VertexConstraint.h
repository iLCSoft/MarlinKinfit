/*! \file 
 *  \brief Declares class VertexConstraint
 *
 * \b Changelog:
 * - 17.11.04 BL: First version
 * - 6.12.04 BL: Instead of constraining 2 tracks to a common vertex,
 *   the constraint now constrains one track to one vertex
 *
 *  this is only the position part of a vertex fit
 *    (no momentum constraint involved)
 *
 * daniel changes it to a BaseHardConstraint
 */ 

#ifndef __VERTEXCONSTRAINT_H
#define __VERTEXCONSTRAINT_H

#include<vector>
#include<cassert>
#include "BaseHardConstraint.h"
#include "ThreeVector.h"
 
class VertexFitObject;
class TrackParticleFitObject;

//  Class VertexConstraint:
/// Constrains a TrackParticleFitObject to a VertexFitObject
/*
 * 2016 DJeans updated for trackparticlefitobject
 *
 * Author: Benno List, Jenny List
 * $Date: 2008/01/30 09:14:55 $
 * $Author: blist $
 *
 */

class VertexConstraint: public BaseHardConstraint {
  public:
    /// Constructor
    VertexConstraint (const VertexFitObject& vertex_, 
                      const TrackParticleFitObject& track_,
                      int ivertex_,
                      int axis
                     );
    /// Virtual destructor
    virtual ~VertexConstraint();
    
    /// Returns the value of the constraint
    virtual double getValue() const;
    
    /// Get first order derivatives. 
    /// Call this with a predefined array "der" with the necessary number of entries!
    virtual void getDerivatives(int idim, double der[]) const;

    virtual bool secondDerivatives(int, int, double*) const;
    virtual bool firstDerivatives(int, double*) const;

    virtual int getVarBasis() const {return BaseDefs::VARBASIS_VXYZ;}
    
  protected:
      
    const VertexFitObject *vertex; 
    const TrackParticleFitObject *track; 
    int   ivertex;
    ThreeVector factor;    

};

#endif /* #ifndef  __VERTEXCONSTRAINT_H */
