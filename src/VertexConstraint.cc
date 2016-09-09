#include "VertexConstraint.h"
#include "TrackParticleFitObject.h"
#include "VertexFitObject.h"

#include<iostream>

#undef NDEBUG
#include<cassert>

using std::cout;
using std::endl;

VertexConstraint::VertexConstraint (const VertexFitObject& vertex_, 
                                    const TrackParticleFitObject& track_,
                                    int ivertex_,
                                    int axis
                                   )
: vertex (&vertex_), track (&track_), ivertex (ivertex_)
{ 

  // ivertex is 0/1 : start or end of track
  // axis is for xyz

  assert (vertex);
  assert (track);
  switch (axis) {
  case 0:
    factor.setValues (1, 0, 0);
    break;
  case 1:
    factor.setValues (0, 1, 0);
    break;
  case 2:
    factor.setValues (0, 0, 1);
    break;
  default:
    cout << "unknown axis " << axis << " allowed values 0->2" << endl;
    assert(0);
  }
}

// destructor
VertexConstraint::~VertexConstraint () {}

// calculate current value of constraint function
double VertexConstraint::getValue() const {

  //   cout << "VertexConstraint: vertex " << vertex->getName()
  //        << " = " << vertex->getVertex()
  //        << ", Track " << track->getName()
  //        << ", vertex " << ivertex 
  //        << " = " << track->getVertex(ivertex)
  //        << "factor = " << factor
  //        << ", result = " << factor*(vertex->getVertex() - track->getVertex (ivertex))
  //        << endl;
  //  

  double value = factor*(vertex->getVertex() - track->getVertex (ivertex));

  //  cout << "VertexConstraint::getValue() : " << vertex->getVertex() << " " << track->getVertex (ivertex) << " " << value << endl;

  return value;
}

// calculate vector/array of derivatives of this contraint 
// w.r.t. to ALL parameters of all fitobjects
void VertexConstraint::getDerivatives(int idim, double der[]) const {

  //  cout << "hello from VertexConstraint::getDerivatives " << endl;

  for (int iglobal = 0; iglobal < idim; iglobal++) der[iglobal] = 0;

  for (int ilocal = 0; ilocal < vertex->getNPar(); ilocal++) {
    if (!vertex->isParamFixed(ilocal)) {
      int iglobal = vertex->getGlobalParNum (ilocal);
      assert (iglobal >= 0 && iglobal < idim);
      der[iglobal] += factor*vertex->getVertexDerivative (ilocal);
    }
  }
  for (int ilocal = 0; ilocal < track->getNPar(); ilocal++) {
    if (!track->isParamFixed(ilocal)) {
      int iglobal = track->getGlobalParNum (ilocal);
      assert (iglobal >= 0 && iglobal < idim);
      der[iglobal] -= factor*track->getVertexDerivative (ivertex, ilocal);
    }
  }

//  cout << "derivative from VertexConstraint : ";
//  for (int iglobal = 0; iglobal < idim; iglobal++) cout <<  der[iglobal] << " ";
//  cout << endl;

  return;
}

// need to be implemented!?

bool VertexConstraint::firstDerivatives(int, double*) const {

  cout << "hello there! should crash!" << endl;

  assert(0);

  return true;
}

bool VertexConstraint::secondDerivatives(int, int, double*) const {

  assert(0);

  return true;
}

