 
#ifndef __SIMPLEPHOTONFITOBJECT_H
#define __SIMPLEPHOTONFITOBJECT_H

#include "ParticleFitObject.h"

// Class SimplePhotonFitObject
/// Class for ISR/Beamstrahlung photons with (p_x,p_y (both fix),  p_z (free)) in kinematic fits
/**
 *
 * Author: Moritz Beckmann
 * $Date: 2010/06/11 20:32:51 $
 * $Author: mbeckman $
 *
 * \b Changelog:
 * - 
 *
 * \b CVS Log messages:
 * - $Log: SimplePhotonFitObject.h,v $
 * - Revision 1.1  2010/06/11 20:32:51  mbeckman
 * - Renamed PhotonFitObjects, cleaned them up for usage
 * -
 * - Revision 1.3  2009/02/23 12:03:18  mbeckman
 * - - PhotonFitObject:     bug fix (1/0), removed dispensable variables
 * - - PhotonFitObjectPxyg: bug fixes (1/0, order of computing variables), modified parametrization
 * -
 * - Revision 1.2  2009/02/18 11:53:42  mbeckman
 * - documentation, debug output
 *
 */ 
class SimplePhotonFitObject : public ParticleFitObject {
  public:
    SimplePhotonFitObject(double px, double py, double pz,   /// initial values for photon (p_x,p_y fix)
                 double Dpz);                          /// sigma for pz
                 
    /// Copy constructor
    SimplePhotonFitObject (const SimplePhotonFitObject& rhs              ///< right hand side
                   );
    
    /// Return a new copy of itself
    virtual SimplePhotonFitObject *copy() const;
    
    /// Assign from anther object, if of same type
    virtual SimplePhotonFitObject& assign (const BaseFitObject& source   ///< The source object
                                 );
    
    /// Assignment               
    SimplePhotonFitObject& operator= (const SimplePhotonFitObject& rhs   ///< right hand side
                             );

    virtual ~SimplePhotonFitObject();
    
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
    
    mutable double pt2, p2, p,
                   dE0, dE1, dE2, 
                   chi2;

    enum {NPAR=3};

};



#endif // __SIMPLEPHOTONFITOBJECT_H
