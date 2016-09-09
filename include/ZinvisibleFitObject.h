/*! \file
 *  \brief Declares class ZinvisibleFitObject
 *
 */

#ifndef __ZINVISIBLEFITOBJECT_H
#define __ZINVISIBLEFITOBJECT_H

#include "ParticleFitObject.h"
#include <cmath>

// Class ZinvisibleFitObject
/// Class for Z->neutrinos with (E, eta, phi) in kinematic fits

class ZinvisibleFitObject : public ParticleFitObject {
  public:
    ZinvisibleFitObject(double E, double eta, double phi, 
			double DE=1, double Dtheta=0.1, double Dphi=0.1, double m = 91.1876); 

    /// Copy constructor
    ZinvisibleFitObject (const ZinvisibleFitObject& rhs              ///< right hand side
		       );

    /// Assignment
    ZinvisibleFitObject& operator= (const ZinvisibleFitObject& rhs   ///< right hand side
				  );

    virtual ~ZinvisibleFitObject();
    
    /// Return a new copy of itself
    virtual ZinvisibleFitObject *copy() const;
    
    /// Assign from anther object, if of same type
    virtual ZinvisibleFitObject& assign (const BaseFitObject& source   ///< The source object
                                      );
    
    /// Get name of parameter ilocal
    virtual const char *getParamName (int ilocal     ///< Local parameter number
                                     ) const;
    
    /// Read values from global vector, readjust vector; return: significant change
    virtual bool   updateParams (double p[],   ///< The parameter vector
                                 int idim      ///< Length of the vector                         
                                );  
    
    // these depend on actual parametrisation!
    virtual double getPx() const;
    virtual double getPy() const;
    virtual double getPz() const;
    virtual double getE() const;
    virtual double getPt() const;
    virtual double getP2() const;
    virtual double getPt2() const;
    virtual double getDPx(int ilocal) const;
    virtual double getDPy(int ilocal) const;
    virtual double getDPz(int ilocal) const;
    virtual double getDE(int ilocal) const;
    
    virtual void invalidateCache() const;

    virtual double getFirstDerivative_Meta_Local( int iMeta, int ilocal , int metaSet ) const;
    virtual double getSecondDerivative_Meta_Local( int iMeta, int ilocal , int jlocal , int metaSet ) const;      
    virtual int getNPar() const {return NPAR;}

  protected:
   inline  double getP() const;
    
    void updateCache() const;

    enum {NPAR=3};
  
    mutable bool cachevalid;
    
    mutable double ctheta, stheta, cphi, sphi,
      p2, p, dpdE, pt, px, py, pz, dptdE,
                   dpxdE, dpydE, dpxdtheta, dpydtheta,
                   chi2;

};

#endif // __ZINVISIBLEFITOBJECT_H

