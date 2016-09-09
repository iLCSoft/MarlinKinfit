 
#ifndef __ISRPHOTONFITOBJECT_H
#define __ISRPHOTONFITOBJECT_H

#include "ParticleFitObject.h"

/// Class ISRPhotonFitObject
/// Documention: arXiv:1006.0436 [hep-ex]

/// Class for ISR photons with (p_x,p_y (both fix),  p_z (free)) in kinematic fits
/// p_z is internally replaced by a parameter p_g
///
/// This class assumes a photon p_z distribution according to dN/d|p_z| = c*|p_z|^(b-1), where the total number of photons should be given by N = \int_{PzMin}^{PzMax} dN/d|p_z| d|p_z|
/// The parameters b, PzMaxB:=PzMax^b and PzMinB:=PzMin^b are required to describe the photon spectrum
/// Only |p_z| values in [PzMin,PzMax[ can be used as start values (assertion in PgFromPz(...)) and will occur in the fit
/// Recommended start value is the missing p_z (fitter will always find a minimum around p_z=0)

/**
 *
 * Author: Moritz Beckmann
 * $Date: 2011/03/16 16:33:24 $
 * $Author: mbeckman $
 *
 * \b Changelog:
 * - 
 *
 * \b CVS Log messages:
 * - $Log: ISRPhotonFitObject.h,v $
 * - Revision 1.2  2011/03/16 16:33:24  mbeckman
 * - Compatibility fixes with ILCSoft svn
 * -
 * - Revision 1.1  2010/06/11 20:32:51  mbeckman
 * - Renamed PhotonFitObjects, cleaned them up for usage
 * -
 * - Revision 1.6  2009/03/26 08:47:16  mbeckman
 * - Bug fix (measured p = 0 instead of start value), extended documentation
 * -
 * - Revision 1.5  2009/02/23 12:03:18  mbeckman
 * - - PhotonFitObject:     bug fix (1/0), removed dispensable variables
 * - - PhotonFitObjectPxyg: bug fixes (1/0, order of computing variables), modified parametrization
 * -
 * - Revision 1.4  2009/02/18 11:53:42  mbeckman
 * - documentation, debug output
 *
 */ 
class ISRPhotonFitObject : public ParticleFitObject {
  public:
    ISRPhotonFitObject(double px, double py, double pz,                   /// initial values for photon (p_x,p_y fix)
                        double b_, double PzMaxB_, double PzMinB_ = 0.);  /// photon spectrum parametrization (see above)
    
    /// Copy constructor
    ISRPhotonFitObject (const ISRPhotonFitObject& rhs              ///< right hand side
                   );
    /// Assignment               
    ISRPhotonFitObject& operator= (const ISRPhotonFitObject& rhs   ///< right hand side
                             );
                             
    virtual ~ISRPhotonFitObject();
    
    /// Return a new copy of itself
    virtual ISRPhotonFitObject *copy() const;
    
    /// Assign from anther object, if of same type
    virtual ISRPhotonFitObject& assign (const BaseFitObject& source   ///< The source object
                                 );
    
    /// Get name of parameter ilocal
    virtual const char *getParamName (int ilocal     ///< Local parameter number
                                      ) const;

    /// Read values from global vector, readjust vector; return: significant change
    virtual bool   updateParams (double p[],   ///< The parameter vector
				 int idim      ///< Length of the vector                         
                                 );  
    
    
    // // these depend on actual parametrisation!
    virtual double getDPx(int ilocal) const;
    virtual double getDPy(int ilocal) const;
    virtual double getDPz(int ilocal) const;
    virtual double getDE(int ilocal) const;

    virtual double getFirstDerivative_Meta_Local( int iMeta, int ilocal , int metaSet ) const;
    virtual double getSecondDerivative_Meta_Local( int iMeta, int ilocal , int jlocal, int metaSet ) const;

    virtual int getNPar() const {return NPAR;}
  
  protected:
    
    enum {NPAR=3}; // well, it's actually 1...Daniel should update

    virtual double PgFromPz(double pz);
    
    void updateCache() const;
  
    mutable bool cachevalid;
    
    mutable double pt2, p2, p, pz,
                   dpx0, dpy0, dpz0, dE0, dpx1, dpy1, dpz1, dE1,
                   dpx2, dpy2, dpz2, dE2, d2pz22, d2E22,
                   chi2,                   
                   b, PzMinB, PzMaxB, dp2zFact;
};

#endif // __ISRPHOTONFITOBJECT_H
