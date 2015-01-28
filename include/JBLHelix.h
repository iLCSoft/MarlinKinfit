/*! \file 
 *  \brief Declares class JBLHelix
 *
 * \b Changelog:
 * - 4.1.05 BL: First version
 *
 */ 


#ifndef __JBLHELIX_H
#define __JBLHELIX_H


class ThreeVector;
class TwoVector;

// Class JBLHelix
/// A helix, parametrized by (kappa, phi0, theta, dca, z0)
/**
 * The main purpose of an JBLHelix is that it can be "intersected"
 * with another JBLHelix object, giving estimates for arc length
 * values s1 and s2 of the two helices where the two helices
 * intersect each other in r/phi and come closest in z
 *
 *  Parameters:
 *  - 0: kappa (signed)
 *  - 1: phi0
 *  - 2: theta
 *  - 3: dca
 *  - 4: z0
 */

class JBLHelix {
  public:
    /// Constructor from (kappa, phi0, theta, dca, z0)
    JBLHelix (double kappa,    ///< kappa
              double phi0,     ///< phi0
              double theta,    ///< theta
              double dca,      ///< dca
              double z0        ///< z0
              );
    /// Constructor from an array of values          
    JBLHelix (double par_[]    ///< Parameters  (kappa, phi0, theta, dca, z0)
             );
    /// Destructor
    ~JBLHelix();
    
    /// Get value of parameter i 
    double getPar (int i        ///< Parameter number i (i=0...4)
                  );
    /// Set value of parameter i 
    JBLHelix& setPar (int i,       ///< Parameter number i (i=0...4)
                      double par_  ///< New parameter value
                     );
    
    /// Get arc lengths for point where another helix comes closest
    /**
     * This routine intersects the helices in the (x, y) plane,
     * and then returns the arclengths corresponding to the
     * best intersection point, where the helices come closest in z.
     * If the helices don't intersect in (x, y), the point where
     * they come closest in (x,y) is returned.
     *
     * The return value indicates: 
     * 1: circles touch, just one intersection point
     * 2: circles intersect, 2 intersection points
     * 0: circles dont intersect, point is best approximation
     */
    int getClosestApproach (const JBLHelix& h1,   ///< The second helix
                            double& s0,           ///< This helix's arc length of the best point
                            double& s1,           ///< The 2nd helix's arc length of the best point
                            double& s02nd,        ///< This helix's arc length of the 2nd best point
                            double& s12nd         ///< The 2nd helix's arc length of the2nd best  point
                           ) const;
     /// Return s of point on helix where helix comes closest in (x,y) to given point                      
     double getClosestS (const TwoVector& p            ///< The point
                        ) const;
     /// Return s of point on helix where helix comes closest in (x,y) to given point                      
     double getClosestS (const ThreeVector& p          ///< The point
                        ) const;
    /// Get point along trajectory into existing 3-vector
    void getTrajectoryPointEx (double s,   
                               ThreeVector& p  
                              ) const; 
    /// Get point along trajectory
    ThreeVector getTrajectoryPoint (double s  
                                   ) const; 
    /// Get z value along trajectory
    double getTrajectoryZ (double s  
                          ) const; 
                 
    /// Get center point of trajectory into existing 2-vector
    void getCenterPointEx (TwoVector& p  
                          ) const; 
    /// Get center point of trajectory
    TwoVector getCenterPoint () const; 
    
    /// Get smallest s that corresponds to same (x, y)
    double getNormalS (double s) const;
    
  private:
    enum {NPAR = 5};
    double par[NPAR];  ///< Parameters (kappa, phi0, theta, dca, z0)
};

#endif // __JBLHELIX_H
