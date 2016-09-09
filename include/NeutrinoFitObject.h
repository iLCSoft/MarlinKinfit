
#ifndef __NEUTRINOFITOBJECT_H
#define __NEUTRINOFITOBJECT_H

#include "ParticleFitObject.h"

#include <cmath>


// Class NeutrinoFitObject
/// Class for neutrinos with (E, eta, phi) in kinematic fits
/**
 *
 * Author: Jenny List, Benno List
 * $Date: 2011/03/03 15:03:02 $
 * $Author: blist $
 *
 * \b Changelog:
 * - 30.12.04 BL: addToGlobCov, getDChi2DParam, getDChi2DParam2,
 *            addToGlobalChi2DerMatrix moved up to ParticleFitObject,
 *            getParamName implemented
 */ 
class NeutrinoFitObject : public ParticleFitObject {
  public:
    NeutrinoFitObject(double E, double theta, double phi, 
                      double DE=1, double Dtheta=0.1, double Dphi=0.1);
                 
    /// Copy constructor
    NeutrinoFitObject (const NeutrinoFitObject& rhs              ///< right hand side
                   );
    /// Assignment               
    NeutrinoFitObject& operator= (const NeutrinoFitObject& rhs   ///< right hand side
                             );

    virtual ~NeutrinoFitObject();
    
    /// Return a new copy of itself
    virtual NeutrinoFitObject *copy() const;
    
    /// Assign from anther object, if of same type
    virtual NeutrinoFitObject& assign (const BaseFitObject& source   ///< The source object
                                      );
    
    /// Get name of parameter ilocal
    virtual const char *getParamName (int ilocal     ///< Local parameter number
                                     ) const;

    /// Read values from global vector, readjust vector; return: significant change
    virtual bool   updateParams (double p[],   ///< The parameter vector
                                 int idim      ///< Length of the vector                         
				 );  
  
    // these depend on actual parametrisation!
    virtual double getDPx(int ilocal) const;
    virtual double getDPy(int ilocal) const;
    virtual double getDPz(int ilocal) const;
    virtual double getDE(int ilocal) const;

    virtual double getFirstDerivative_Meta_Local( int iMeta, int ilocal , int metaSet ) const;
    virtual double getSecondDerivative_Meta_Local( int iMeta, int ilocal , int jlocal , int metaSet ) const;

    virtual int getNPar() const {return NPAR;}
      
  protected:
    void updateCache() const;
    
    enum {NPAR=3};
  
    mutable double ctheta, stheta, cphi, sphi,
                   pt, px, py, pz, dptdE, 
                   dpxdE, dpydE, dpxdtheta, dpydtheta,
                   chi2;

};

#endif // __NEUTRINOFITOBJECT_H

