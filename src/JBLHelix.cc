/*! \file 
 *  \brief Implements class JBLHelix
 *
 * \b Changelog:
 * - 4.1.05 BL: First version
 *
 */ 
 
#include "JBLHelix.h"
#include "ThreeVector.h"
#include "TwoVector.h"

#undef NDEBUG
#include <cassert>
#include <cmath>

JBLHelix::JBLHelix (double kappa,    
                    double phi0,   
                    double theta,  
                    double dca,    
                    double z0      
                    ) {
  setPar (0, kappa);
  setPar (1, phi0);
  setPar (2, theta);
  setPar (3, dca);
  setPar (4, z0);
}         

JBLHelix::JBLHelix (double par_[]) {
  for (int i = 0; i < NPAR; i++) setPar (i, par_[i]);        
}         

JBLHelix::~JBLHelix() 
{}


double JBLHelix::getPar (int i) {
  assert (i >= 0 && i < NPAR);
  return par[i];          
}              

JBLHelix& JBLHelix::setPar (int i, double par_) {
  assert (i >= 0 && i < NPAR);
  par[i] = par_;
  return *this;
}

int JBLHelix::getClosestApproach (const JBLHelix& h1, 
                                  double& s0,       
                                  double& s1,       
                                  double& s02nd,    
                                  double& s12nd     
                                 ) const {
  if (std::abs(par[0]) > 1E-7 &&
      std::abs(h1.par[0]) > 1E-7) {
    // We really have two helices, i.e. circles in (x, y)
    TwoVector cent0 = getCenterPoint();
    TwoVector cent1 = h1.getCenterPoint();
    // std::cout << "cent0=" << cent0 << ", cent1=" << cent1 << "\n";
    
    double dist2 = (cent0-cent1).getMag2();
    double dist = std::sqrt(dist2);
    double r0 = 1/par[0];
    double r1 = 1/h1.par[0];
    double r0a = std::abs(r0);
    double r1a = std::abs(r1);
    s0 = s02nd = getClosestS (cent1);
    s1 = s12nd = h1.getClosestS (cent0);
    if (dist > r0a+r1a) return 0;
    if (dist == r0a+r1a) return 1;
    double psi0 = std::acos((r0*r0 + dist2 - r1*r1)/(2*r0a*dist));
    double psi1 = std::acos((r1*r1 + dist2 - r0*r0)/(2*r1a*dist));
    double s0try1 =    getNormalS(s0 + r0*psi0);
    double s1try1 = h1.getNormalS(s1 - r1*psi1);
    double dz1 = std::abs(getTrajectoryZ(s0try1) - h1.getTrajectoryZ(s1try1));
    double s0try2 =    getNormalS(s0 - r0*psi0);
    double s1try2 = h1.getNormalS(s1 + r1*psi1);
    double dz2 = std::abs(getTrajectoryZ(s0try2) - h1.getTrajectoryZ(s1try2));
    
    // std::cout << "s0try1=" << s0try1 << "=>" << getTrajectoryPoint (s0try1)
    //           << ", s1try1=" << s1try1 << "=>" << h1.getTrajectoryPoint (s1try1) 
    //           << ". dz1 = " << dz1 << "\n";
    // std::cout << "s0try2=" << s0try2 << "=>" << getTrajectoryPoint (s0try2)
    //           << ", s1try2=" << s1try2 << "=>" << h1.getTrajectoryPoint (s1try2) 
    //           << ". dz2 = " << dz2 << "\n";
    
    if (dz1 < dz2) {
      s0 = s0try1;
      s1 = s1try1;
      s02nd = s0try2;
      s12nd = s1try2;
    }
    else {
      s0 = s0try2;
      s1 = s1try2;
      s02nd = s0try1;
      s12nd = s1try1;
    } 
    return 2; 
  }
  else {
    assert (0);
  }    
}


double JBLHelix::getClosestS (const TwoVector& p) const {
  if (std::abs (par[0]) < 1.e-8) {
    // If helix is almost a straight line,
    // treat it as straight line here.
    // Get point of closest approach to origin
    double sphi0 = std::sin(par[1]);
    double cphi0 = std::cos(par[1]);
    double xc = sphi0*par[3];
    double yc = -cphi0*par[3];
    return (p.getX()-xc)*cphi0 + (p.getY()-yc)*sphi0;
  }
  else {
    // Get center of circle
    double r = 1/par[0];
    double dcamir =  par[3] - r;
    double xc =  dcamir*std::sin(par[1]);
    double yc = -dcamir*std::cos(par[1]);
    double psi = (r > 0) ? 
                 std::atan2(p.getX()-xc, -p.getY()+yc) - par[1] :
                 std::atan2(p.getX()-xc,  p.getY()-yc) + par[1];
    // if (psi < -M_PI) psi += 2*M_PI;
    // else if (psi > M_PI) psi -= 2*M_PI;
    return getNormalS (std::abs(r)*psi);
  }
}

double JBLHelix::getClosestS (const ThreeVector& p) const {
  return getClosestS (TwoVector (p.getX(), p.getY()));
}

void JBLHelix::getTrajectoryPointEx (double s, ThreeVector& p ) const {
  double sphi0 = std::sin(par[1]);
  double cphi0 = std::cos(par[1]);
  double kappas = par[0]*s;
  if (std::abs (kappas) < 1.e-6) {
    double dcamikssq = par[3] - 0.5*kappas*s;
    p.setValues ( sphi0*dcamikssq + cphi0*s,
                 -cphi0*dcamikssq + sphi0*s,
                  par[4] + s*std::cos(par[2])/std::sin(par[2]));
  }
  else {
    double phi0plkappas = par[1] + kappas;
    double r = 1/par[0];
    double dcamir = par[3] - r;
    p.setValues ( dcamir*sphi0 + r*sin(phi0plkappas),
                 -dcamir*cphi0 - r*cos(phi0plkappas),
                 par[4] + s*std::cos(par[2])/std::sin(par[2]));
  }
}
double JBLHelix::getTrajectoryZ (double s) const {
  return par[4] + s*std::cos(par[2])/std::sin(par[2]);
}

ThreeVector JBLHelix::getTrajectoryPoint (double s) const {
  ThreeVector result;
  getTrajectoryPointEx (s, result);
  return result;
}

void JBLHelix::getCenterPointEx (TwoVector& p) const {
  double r = 1/par[0];
  double dcamir =  par[3] - r;
  p.setValues ( dcamir*std::sin(par[1]), -dcamir*std::cos(par[1]));
}
TwoVector JBLHelix::getCenterPoint () const {
  TwoVector result;
  getCenterPointEx (result);
  return result;
}

double JBLHelix::getNormalS (double s) const {
  double kappas = par[0]*s;
  if (kappas >= -M_PI &&kappas < M_PI) return s;
  return std::atan2 (std::sin(kappas), std::cos (kappas))/par[0];
}
