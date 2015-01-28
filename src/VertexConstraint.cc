#include "VertexConstraint.h"
#include "TrackFitObject.h"
#include "VertexFitObject.h"

#include<iostream>
#include<cassert>

using std::cout;
using std::endl;

VertexConstraint::VertexConstraint (const VertexFitObject& vertex_, 
                                    const TrackFitObject& track_,
                                    int ivertex_,
                                    int axis
                                   )
: vertex (&vertex_), track (&track_), ivertex (ivertex_)
{ 
  assert (vertex);
  assert (track);
  switch (axis) {
    case 1:
      factor.setValues (0, 10, 0);
      break;
    case 2:
      factor.setValues (0, 0, 10);
      break;
    default:
      factor.setValues (10, 0, 0);
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
  return factor*(vertex->getVertex() - track->getVertex (ivertex));
}

// calculate vector/array of derivatives of this contraint 
// w.r.t. to ALL parameters of all fitobjects
void VertexConstraint::getDerivatives(int idim, double der[]) const {
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
}
  
  
