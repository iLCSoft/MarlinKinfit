/*! \file 
 *  \brief Implements class ParticleConstraint
 *
 * \b Changelog:
 * - 
 *
 * \b CVS Log messages:
 * - $Log: ParticleConstraint.cc,v $
 * - Revision 1.2  2008/10/17 13:17:16  blist
 * - Avoid variable-size arrays
 * -
 * - Revision 1.1  2008/02/12 10:19:09  blist
 * - First version of MarlinKinfit
 * -
 * - Revision 1.2  2008/02/07 08:21:07  blist
 * - ParticleConstraint.C fixed
 * -
 * - Revision 1.1  2008/02/07 08:18:57  blist
 * - ParticleConstraint,C added
 * -
 */ 

#include "ParticleConstraint.h"
#include "ParticleFitObject.h"
#include <iostream>
#include <cmath>
using namespace std;


// probably these can also be moved to basehardconstraint?

