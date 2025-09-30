/*! \file 
 *  \brief Implements class VertexFitObject
 *
 * \b Changelog:
 *
 * \b CVS Log messages:
 * - $Log: VertexFitObject.C,v $
 * - Revision 1.2  2007/09/13 13:33:06  blist
 * - Print methods return os
 * -
 *
 */ 

#include "VertexFitObject.h"
#include "TwoVector.h"
#include "ThreeVector.h"
#include "FourVector.h"
#include "VertexConstraint.h"
#include "MomentumConstraint.h"
#include "TrackParticleFitObject.h"
#include "BaseFitter.h"
#include "JBLHelix.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include <iostream>
using std::cout;
using std::endl;

#undef NDEBUG
#include <cassert>
#include <cmath>
using std::isfinite;

#include <cstring>

static int debug = 0;

VertexFitObject::VertexFitObject(const char *name_,
                                 double x,
                                 double y,
                                 double z
                                )
: tracks (0), constraints (0)
{

  assert( int(NPAR) <= int(BaseDefs::MAXPAR) );

  setParam (0, x, false);
  setParam (1, y, false);
  setParam (2, z, false);
  setMParam (0, x);
  setMParam (1, y);
  setMParam (2, z);
  setName (name_);
  initCov();

}

// We get a warning that BaseFitObject should be explicitly initialized
// here, but I don't want to change this part because, I think everything is
// done properly already and not changing behavior is more important.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wextra"
VertexFitObject::VertexFitObject (const VertexFitObject& rhs)
{
  //  copy (rhs);
  VertexFitObject::assign (rhs);
}
#pragma GCC diagnostic pop


VertexFitObject& VertexFitObject::operator= (const VertexFitObject& rhs) {
  if (this != &rhs) {
    assign (rhs); // calls virtual function assign of derived class
  }
  return *this;
}

VertexFitObject *VertexFitObject::copy() const {
  return new VertexFitObject (*this);
}
    
VertexFitObject& VertexFitObject::assign (const BaseFitObject& source) {
  if (const VertexFitObject *psource = dynamic_cast<const VertexFitObject *>(&source)) {
    if (psource != this) {
      BaseFitObject::assign (source);
      tracks = psource->tracks; // daniel thinks we probably need these too
      constraints = psource->constraints;
    }
  }
  else {
    assert (0);
  }
  return *this;
}

VertexFitObject::~VertexFitObject()
{
  for (CIterator it = constraints.begin(); it != constraints.end(); ++it) {
    delete (*it);
    (*it) = 0;
  }
}

const char *VertexFitObject::getParamName (int ilocal) const {
  switch (ilocal) {
    case 0: return "x";
    case 1: return "y";
    case 2: return "z";
  }
  return "undefined";
}

std::ostream& VertexFitObject::print (std::ostream& os) const {
  os << getName() << ":\n";
  printParams(os);
  return os;
}

ThreeVector VertexFitObject::getVertex () const {
  ThreeVector result;
  getVertexEx (result);
  return result;
}

void VertexFitObject::getVertexEx (ThreeVector& p) const {
  p.setValues (par[0], par[1], par[2]);
} 
                                                    
ThreeVector VertexFitObject::getVertexDerivative (int ilocal) const  {
  ThreeVector result;
  getVertexDerivativeEx (ilocal, result);
  return result;
}

void VertexFitObject::getVertexDerivativeEx (int ilocal, ThreeVector& p) const {
  switch (ilocal) {
    case 0: // d / d par[0] = d / d x
      p.setValues (1, 0, 0);
      break;
    case 1: // d / d par[1] = d / d y
      p.setValues (0, 1, 0);
      break;
    case 2: // d / d par[2] = d / d z
      p.setValues (0, 0, 1);
      break;
    default: // should never happen!
      assert (0);
  }
}

double VertexFitObject::getFirstDerivative_Meta_Local(int iMeta, int ilocal , int metaSet) const {
  // this is d I_iMeta / d p_ilocal
  // I = intermediate variable = XYZ
  // p = local parameter of vertex = xyz
  // in this case they are the same, so it's trivial
  assert ( metaSet==BaseDefs::VARBASIS_VXYZ );
  assert ( iMeta < BaseDefs::nMetaVars[metaSet] );
  assert ( ilocal>=0 && ilocal<NPAR );

  double deriv = iMeta==ilocal ? 1. : 0.;

  return deriv;
}

double VertexFitObject::getSecondDerivative_Meta_Local(int iMeta, int ilocal , int jlocal , int metaSet) const {
  assert ( metaSet==BaseDefs::VARBASIS_VXYZ );
  assert ( iMeta < BaseDefs::nMetaVars[metaSet] );
  assert ( ilocal>=0 && ilocal<NPAR );
  assert ( jlocal>=0 && jlocal<NPAR );
  return 0;
}




void VertexFitObject::addVertexConstraints (BaseFitter& fitter, int axis) {

  // this adds the vertex contraints to the fitter

  //  cout << "hello from addVertexConstraints: axis " << axis << endl;

  for (TIterator it = tracks.begin(); it != tracks.end(); ++it) {
    assert(it->track); 

    int iTrkVtx = it->inbound ? 1 : 0 ;

    // un-fix the vertex parameter of this track
    it->track->releaseVertexParam( iTrkVtx );

    VertexConstraint *con = new VertexConstraint (*this, *(it->track), iTrkVtx, axis);
    con->setName (this->getName());

    //    cout << "vertex constraint value = " << con->getValue() << endl;

    fitter.addConstraint (con);
    constraints.push_back (con);

    it->track->releaseVertexParam (it->inbound ? 1 : 0);
  }
}

void VertexFitObject::addMomentumConstraint (BaseFitter& fitter, int axis) {
  MomentumConstraint *con = new MomentumConstraint (axis);
  fitter.addConstraint (con);
  constraints.push_back (con);
  
  // Add first inbound track
  for (TIterator it = tracks.begin(); it != tracks.end(); ++it) {
    if (it->inbound) {
      assert(it->track); 
      con->addToFOList(*(it->track), 1);
      break;
    }
  }
  // Now add outgoing tracks
  for (TIterator it = tracks.begin(); it != tracks.end(); ++it) {
    assert(it->track); 
    if (!it->inbound) con->addToFOList(*(it->track), 0);
  }
}

void VertexFitObject::addConstraints (BaseFitter& fitter, int mask) {

  if (mask & VX) addVertexConstraints (fitter, 0); else fixParam (0);
  if (mask & VY) addVertexConstraints (fitter, 1); else fixParam (1);
  if (mask & VZ) addVertexConstraints (fitter, 2); else fixParam (2);
  if (mask & PX) addMomentumConstraint (fitter, 1);
  if (mask & PY) addMomentumConstraint (fitter, 2);
  if (mask & PZ) addMomentumConstraint (fitter, 3);
  if (mask & E)  addMomentumConstraint (fitter, 0);

}

void VertexFitObject::addTrack (TrackParticleFitObject *track, bool inbound, bool ismeasured) {
  assert (track);
  tracks.push_back (TrackDescriptor (track, inbound, ismeasured));
}

ThreeVector VertexFitObject::estimatePosition () {
  if (debug) cout << "VertexFitObject::estimatePosition(): starting" << endl;
  ThreeVector position (0, 0, 0);

  ThreeVector commonRefPoint(0,0,0);

  // DANIEL removed for now (to get it to compile without the JBLhelix class)
   int n = 0;
   for (TIterator it0 = tracks.begin(); it0 != tracks.end(); ++it0) {
     if (it0->measured) {
       JBLHelix h0 = it0->track->getJBLHelix(commonRefPoint);
       TIterator it1 = it0;
       for (++it1; it1 != tracks.end(); ++it1) {
         if (it1->measured) {
           if (debug)
             cout << "Intersecting " << it0->track->getName()
                  << " and " << it1->track->getName() << endl;
           JBLHelix h1 = it1->track->getJBLHelix(commonRefPoint);
           double s0, s1, s02nd, s12nd;
           h0.getClosestApproach (h1, s0, s1, s02nd, s12nd);
           position += h0.getTrajectoryPoint (s0);
           position += h1.getTrajectoryPoint (s1);
           if (debug)
             cout << "  point 0: " << h0.getTrajectoryPoint (s0)
                  << ", point 1: " << h1.getTrajectoryPoint (s1) << endl;
           n += 2;
         }
       }
     }
   }
   position *= (1./n);
   if (debug) cout << "Final position estimate: " << position << endl;

  return position;
}

void VertexFitObject::initForFit() {
  if (debug) 
    cout << "VertexFitObject::initForFit(): starting for " << getName() << endl;

  // Estimate and set vertex position
  ThreeVector position = estimatePosition ();
  for (int i = 0; i < 3; i++) par[i] = position.getComponent (i);

  //  cout << "estimated position " << position << " " << getVertex() << endl;

  if (debug) 
    cout << "VertexFitObject now: " << *this << endl;
  if (debug) {
    cout << "TrackFit objects before adjustment: \n";
    for (TIterator it = tracks.begin(); it != tracks.end(); ++it) {
      cout << "   " << it->track->getName() << " = " << *(it->track)
           << (it->measured ? "  measured " : "  not meas.")
           << (it->inbound ? " decays at " : " starts at ")
           << it->track->getVertex(it->inbound ? 1 : 0)
           << endl;
    }
  }  
    
  // For all measured tracks, set s to be close to vertex position,
  // at the same time estimate total 4-momentum and total charge
  FourVector ptot;
  double chargetot = 0;
  int nunm = 0;

  //  cout << " VertexFitObject::initForFit hello 1 " << endl;

  for (TIterator it = tracks.begin(); it != tracks.end(); ++it) {
    if (it->measured) {
      if (it->inbound) {
	//	cout << " VertexFitObject::initForFit hello 2 " << par[0] << " " << par[1] << endl;
        it->track->setVertex (1, TwoVector (par[0], par[1]));
        ptot -= it->track->getMomentum (1);
        chargetot -= it->track->getCharge();
      }
      else {
	//	cout << " VertexFitObject::initForFit hello 3 " << par[0] << " " << par[1] << endl;
        it->track->setVertex (0, TwoVector (par[0], par[1]));
        ptot += it->track->getMomentum (0);
        chargetot += it->track->getCharge();
      }
    }
    else {
      nunm++;
    }
  }
  assert (nunm <= 1);
  // Now, init remaining unmeasured track
  for (TIterator it = tracks.begin(); it != tracks.end(); ++it) {
    if (!it->measured) {
      if (it->inbound) {
        assert (ptot.getE() > 0);
        it->track->setParameters (1, position, ptot, chargetot);
      }
      else {
        assert (ptot.getE() < 0);
        it->track->setParameters (0, position, -ptot, -chargetot);
      }
    }
  }
  if (debug) {
    cout << "TrackFit objects after adjustment: \n";
    for (TIterator it = tracks.begin(); it != tracks.end(); ++it) {
      cout << "   " << it->track->getName() << " = " << *(it->track)
           << (it->measured ? "  measured " : "  not meas.")
           << (it->inbound ? " decays at " : " starts at ")
           << it->track->getVertex(it->inbound ? 1 : 0)
           << endl;
    }
  }  
  return;
}
