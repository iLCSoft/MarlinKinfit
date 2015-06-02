/*! \file 
 *  \brief Declares class VertexConstraint
 *
 * \b Changelog:
 * - 17.11.04 BL: First version
 * - 6.12.04 BL: Instead of constraining 2 tracks to a common vertex,
 *   the constraint now constrains one track to one vertex
 *
 */ 

#ifndef __VERTEXCONSTRAINT_H
#define __VERTEXCONSTRAINT_H

#include<vector>
#include<cassert>
#include "BaseConstraint.h"
#include "ThreeVector.h"
 
class VertexFitObject;
class TrackFitObject;

//  Class VertexConstraint:
/// Constrains a TrackFitObject to a VertexFitObject
/**
 * 
 *
 * Author: Benno List, Jenny List
 * $Date: 2008/01/30 09:14:55 $
 * $Author: blist $
 *
 */

class VertexConstraint: public BaseConstraint {
  public:
    /// Constructor
    VertexConstraint (const VertexFitObject& vertex_, 
                      const TrackFitObject& track_,
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
    /// Adds first order derivatives to global covariance matrix M
    virtual void add1stDerivativesToMatrix (double *M,      ///< Global covariance matrix, dimension at least idim x idim
                                            int idim        ///< First dimension of array der
                                            ) const
    {assert (false);}
    /// Adds second order derivatives to global covariance matrix M
    virtual void add2ndDerivativesToMatrix (double *M,      ///< Global covariance matrix, dimension at least idim x idim
                                            int idim,       ///< First dimension of array der
                                            double lambda   ///< Lagrange multiplier for this constraint
                                            ) const
    {assert (false);}
    /// Add lambda times derivatives of chi squared to global derivative matrix
    virtual void addToGlobalChi2DerVector (double *y,   ///< Vector of chi2 derivatives
                                           int idim,    ///< Vector size 
                                           double lambda //< The lambda value
                                           ) const
    {assert (false);}
    
    
    /// Accesses position of constraint in global constraint list
    virtual int  getGlobalNum() const {return globalNum;}
    /// Sets position of constraint in global constraint list
    virtual void setGlobalNum(int iglobal) {globalNum = iglobal;}
    
  protected:
      
    /// Position of constraint in global constraint list
    int globalNum;
    
    const VertexFitObject *vertex; 
    const TrackFitObject *track; 
    int   ivertex;
    ThreeVector factor;    

};

#endif /* #ifndef  __VERTEXCONSTRAINT_H */
