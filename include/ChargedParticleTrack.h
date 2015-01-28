/*! \file 
 *  \brief Declares class ChargedParticleTrack
 *
 * \b Changelog:
 * - 17.11.04 BL First version
 *
 * \b CVS Log messages:
 * - $Log: ChargedParticleTrack.h,v $
 * - Revision 1.4  2008/01/30 09:14:53  blist
 * - Preparations for NewtonFitter
 * -
 * - Revision 1.3  2007/09/17 12:50:15  blist
 * - Some parameters reordered
 * -
 * - Revision 1.2  2007/09/14 10:58:42  blist
 * - Better documentation,
 * - added PyConstraint::add1stDerivativesToMatrix,
 * - fixed bug in PzConstraint::add1stDerivativesToMatrix
 * -
 *
 */ 

#ifndef __CHARGEDPARTICLETRACK_H
#define __CHARGEDPARTICLETRACK_H

#include "TrackFitObject.h"

#include <cassert>

// Class ChargedParticleTrack
/// A helix track of a charged particle with start and stop vertex
/**
 *  Parameters:
 *  - 0: kappa (signed: kappa < 0 <=> charge > 0)
 *  - 1: phi0
 *  - 2: theta
 *  - 3: dca
 *  - 4: z0
 *  - 5: mass
 *  - 6: sstart (s value of start vertex)
 *  - 7: sstop (s value of stop vertex)
 */

class ChargedParticleTrack: public TrackFitObject {
  public:
    /// Constructor from (kappa, phi0, theta, dca, z0, mass, sstart, sstop)
    ChargedParticleTrack (const char *name_,
                          double kappa,
                          double phi0,
                          double theta,
                          double dca,
                          double z0,
                          double mass,
                          double charge_,
                          double sstart,
                          double sstop
                          );
    /// Constructor from (kappa, phi0, theta, dca, z0, mass, sstart)
    ChargedParticleTrack (const char *name_,
                          double kappa,
                          double phi0,
                          double theta,
                          double dca,
                          double z0,
                          double mass,
                          double charge_,
                          double sstart
                          );
    /// Constructor from a (start) vertex and a momentum vector
    ChargedParticleTrack (const char *name_,
                          const ThreeVector& vertex,
                          const ThreeVector& momentum,
                          double mass_,
                          double charge_ = 1
                         );
    /// Constructor from an array of parameters, an array of covariances
    ChargedParticleTrack (const char *name_,
                          const float par_[5],
                          const float cov_[15],
                          double mass_,
                          double charge_ = 1,
                          double sstart = 0,
                          double sstop = 0
                         );
    /// Return a new copy of itself
    virtual ChargedParticleTrack *copy() const ;
    
    /// Assign from anther object, if of same type
    virtual ChargedParticleTrack& assign (const BaseFitObject& source   ///< The source object
                                 ) ;
    
    /// Destructor
    virtual ~ChargedParticleTrack();

    /// Get number of parameters of this FitObject
    virtual int getNPar() const;
    
    /// Set value and measured flag of parameter ilocal; return=success
    virtual bool   setParam (int ilocal,         ///< Local parameter number
                             double par_,        ///< New parameter value
                             bool measured_,     ///< New "measured" flag
                             bool fixed_ = false ///< New "fixed" flag
                            );  
    /// Set value of parameter ilocal; return=success
    virtual bool   setParam (int ilocal,         ///< Local parameter number
                             double par_        ///< New parameter value
                            );  

    /// Set parameters such that track passes through a vertex with a given 4-momentum; return=success
    virtual bool setParameters (int ivertex,                  ///< Vertex number: 0=start, 1=stop
                                const ThreeVector& vertex,    ///< Vertex position
                                const FourVector& momentum,   ///< Four-momentum
                                double charge_                ///< Charge (signed, in units of e)
                               );
                                
    /// Get name of parameter ilocal
    virtual const char *getParamName (int ilocal     ///< Local parameter number
                                     ) const;
 
    /// Get point along trajectory into existing 3-vector
    virtual void getTrajectoryPointEx (double s,  
                                       ThreeVector& p 
                                      ) const; 
                                    
    /// Get start (ivertex=0) or stop (ivertex=1) vertex into existing 3-vector
    virtual void getVertexEx (int ivertex, 
                              ThreeVector& p
                             ) const; 

    /// Set start (i=0) or stop (i=1) vertex to a point as close as possible to given point
    virtual void setVertex (int ivertex,         ///< Vertex number: 0=start, 1=stop
                             const TwoVector& p  ///< Vertex position
                           ); 
    
    /// Get derivative of trajectory w.r.t. parameter ilocal into  existing 3-vector
    virtual void getTrajectoryDerivativeEx (double s, 
                                            int ilocal, 
                                            ThreeVector& p
                                           ) const;
    /// Get derivative of vertex w.r.t. parameter ilocal into  existing 3-vector
    virtual void getVertexDerivativeEx (int ivertex, 
                                        int ilocal, 
                                        ThreeVector& p
                                       ) const;

    /// Get momentum along trajectory into existing 4-vector
    virtual void getMomentumAtTrajectoryEx (double s, 
                                            FourVector& p
                                           ) const; 
     /// Get momentum at vertex into existing 4-vector
    virtual void getMomentumEx (int ivertex, 
                                FourVector& p 
                               ) const; 

    /// Get charge in units of e
    virtual double getCharge () const; 

    /// Get momentum derivative along trajectory into existing 4-vector
    virtual void getMomentumDerivativeAtTrajectoryEx (double s, 
                                                      int ilocal, 
                                                      FourVector& p
                                                     ) const; 
    /// Get derivative of momentum w.r.t. parameter ilocal into  existing 4-vector
    virtual void getMomentumDerivativeEx (int ivertex, 
                                          int ilocal, 
                                          FourVector& p
                                         ) const;
    /// Get s (arclength in r/phi) of vertex i
    virtual double getArcLength (int i
                                ) const;
       
    /// Get derivative of chi squared w.r.t. parameter ilocal
    virtual double getDChi2DParam(int ilocal) const {assert (0);};
    /// Get second derivative of chi squared w.r.t. parameters ilocal1 and ilocal2
    virtual double getD2Chi2DParam2(int ilocal1, int ilocal2) const {assert (0);};
    
    /// Add derivatives of chi squared to global covariance matrix
    virtual void addToGlobalChi2DerMatrix (double *M,   ///< Global covariance matrix
                                           int idim     ///< First dimension of global covariance matrix
                                           ) const
    {assert (0);}
    /// Add derivatives of chi squared to global derivative matrix
    virtual void addToGlobalChi2DerVector (double *y,   ///< Vector of chi2 derivatives
                                           int idim     ///< Vector size 
                                           ) const
    {assert (0);}

    /// Get helix that is tangential at a certain arc length s
    virtual JBLHelix getTangentialHelix (double s     ///<  Arc length
                                        );
    /// get Mass
    virtual double getMass() const;
    
    /// Factor between internal and external parameters
    static const double parfact[NPARMAX];
    
    /// Fix parameter(s) pertaining to a vertex, or release it
    virtual bool fixVertexParam (int ivertex,    ///< Vertex number
                                 bool fix=true  ///< fix if true, release if false
                                 ) { return fixParam (6+ivertex, fix); }
  
  protected:
    /// Update the cache values
    void updateCache() const;
    
    /// init covariance matrix from an array
    virtual void initCov(const float cov_[15]);
    
    /// Get smallest s that corresponds to same (x, y)
    double getNormalS (double s) const;
  
    /// Number of parameters
    enum {NPAR = 8};
    
    // Cache variables
    mutable double kappa;         ///< kappa = parfact[0]*par[0]
    mutable double phi0;          ///< phi0  = parfact[1]*par[1]
    mutable double theta;         ///< theta = parfact[2]*par[2]
    mutable double dca;           ///< dca   = parfact[3]*par[3]
    mutable double z0;            ///< z0    = parfact[4]*par[4]
    mutable double mass;          ///< mass  = parfact[5]*par[5]
    mutable double r;             ///< 1/kappa
    mutable double sphi0;         ///< sin(phi0)
    mutable double cphi0;         ///< cos(phi0)
    mutable double dcamir;        ///< dca - r
    mutable double cottheta;      ///< cotan (theta)
    mutable double sintheta;      ///< sin (theta)
    mutable double sin2theta;     ///< sin^2 (theta)
    mutable double cBq;           ///< cB * charge
    mutable double pt;            ///< -cBq/kappa
    mutable double momentum;      ///< -cBq/[kappa*sin(theta)];
    mutable double momderfact;    ///< +cBq/kappa^2;
    mutable double energy;        ///< sqrt ((cB/[kappa*sin(theta)])^2 + mass^2)
    mutable double beta;          ///< momentum/energy;
  
    double charge;                ///< charge (in units of e)


};

#endif // __CHARGEDPARTICLETRACK_H
