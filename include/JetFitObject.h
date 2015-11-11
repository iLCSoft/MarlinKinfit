/*! \file 
 *  \brief Declares class JetFitObject
 *
 * \b Changelog:
 * - 
 *
 * \b CVS Log messages:
 * - $Log: JetFitObject.h,v $
 * - Revision 1.3  2011/03/03 15:03:02  blist
 * - Latest version, with NewFitterGSL
 * -
 * - Revision 1.2  2009/02/17 12:46:34  blist
 * - Improved version of NewtonFitterGSL, JetFitObject changed
 * -
 * - Revision 1.1  2008/02/12 10:19:05  blist
 * - First version of MarlinKinfit
 * -
 * - Revision 1.9  2008/02/04 17:30:53  blist
 * - NewtonFitter works now!
 * -
 * - Revision 1.8  2008/01/30 21:48:02  blist
 * - Newton Fitter still doesnt work :-(
 * -
 * - Revision 1.7  2008/01/30 09:14:53  blist
 * - Preparations for NewtonFitter
 * -
 * - Revision 1.6  2008/01/29 17:23:00  blist
 * - new addTo2ndDerivatives and setParam
 * -
 * - Revision 1.5  2007/09/17 12:50:15  blist
 * - Some parameters reordered
 * -
 * - Revision 1.4  2007/09/13 08:09:50  blist
 * - Updated 2nd derivatives for px,py,pz,E constraints, improved header documentation
 * -
 *
 */ 
 
#ifndef __JETFITOBJECT_H
#define __JETFITOBJECT_H

#include "ParticleFitObject.h"

// Class JetFitObject
/// Class for jets with (E, eta, phi) in kinematic fits
/**
 *
 * Author: Jenny List, Benno List
 * $Date: 2011/03/03 15:03:02 $
 * $Author: blist $
 *
 * \b Changelog:
 * - 30.12.04 BL: addToGlobCov, getDChi2DParam, getDChi2DParam2, addToGlobalChi2DerMatrix
 *            moved up to ParticleFitObject,
 *            getParamName implemented
 */ 
class JetFitObject : public ParticleFitObject {
  public:
    JetFitObject(double E, double theta, double phi, 
                 double DE, double Dtheta, double Dphi, 
                 double m = 0);
                 
    /// Copy constructor
    JetFitObject (const JetFitObject& rhs              ///< right hand side
                   );
    /// Assignment               
    JetFitObject& operator= (const JetFitObject& rhs   ///< right hand side
                             );

    virtual ~JetFitObject();
    
    /// Return a new copy of itself
    virtual JetFitObject *copy() const;
    
    /// Assign from anther object, if of same type
    virtual JetFitObject& assign (const BaseFitObject& source   ///< The source object
                                 );
    
    /// Get name of parameter ilocal
    virtual const char *getParamName (int ilocal     ///< Local parameter number
                                     ) const;
 
    /// Read values from global vector, readjust vector; return: significant change
    virtual bool   updateParams (double p[],   ///< The parameter vector
                                 int idim      ///< Length of the vector                         
                                );  

    virtual int getNPar() const {return NPAR;}
        
    // these depend on actual parametrisation!
     
    virtual double getDPx(int ilocal) const;
    virtual double getDPy(int ilocal) const;
    virtual double getDPz(int ilocal) const;
    virtual double getDE(int ilocal) const;

    virtual double getCov (int ilocal,    ///< Local parameter number i                                                                                                 
                           int jlocal     ///< Local parameter number j                                                                                             
			   ) const ;

    /// Get error of parameter ilocal
    virtual double getError (int ilocal     ///< Local parameter number
			     ) const;

    
    /// add derivatives to vector der of size idim
    /// pxfact*dpx/dx_i + pyfact*dpy/dx_i + pzfact*dpz/dx_i + efact*dE/dx_i 
//    virtual void   addToDerivatives (double der[],      ///< Derivatives vector, length idim
//                                     int idim,          ///< Length of derivatives vector
//                                     double efact=0,    ///< Factor for dE/dx_i
//                                     double pxfact=0,   ///< Factor for dpx/dx_i
//                                     double pyfact=0,   ///< Factor for dpy/dx_i 
//                                     double pzfact=0    ///< Factor for dpz/dx_i
//                                     ) const;

    // daniel's new method
    // derivatives of intermediate variable wrt local variable
    virtual double getFirstDerivative( int iMeta, int ilocal , int metaSet ) const;
    virtual double getSecondDerivative( int iMeta, int ilocal , int jlocal , int metaSet ) const;

    /// Get chi squared from measured and fitted parameters
    //    virtual double getChi2() const;

  protected:
    
    enum {NPAR=3};
    
    void updateCache() const;

    mutable double ctheta, stheta, cphi, sphi,
                   p2, p, pt, px, py, pz, dpdE, dptdE, 
                   dpxdE, dpydE, dpzdE, dpxdtheta, dpydtheta,
                   chi2;
                   // d2pdE2, d2ptsE2;
                   
    /// Adjust E, theta and phi such that E>=m, 0<=theta<=pi, -pi <= phi < pi; returns true if anything was changed             
    static bool adjustEThetaPhi (double& m, double &E, double& theta, double& phi);
    
    /// Calculate chi2 
    //    double calcChi2 () const;
    
};



#endif // __JETFITOBJECT_H

