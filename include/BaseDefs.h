#ifndef BASEDEFS_HH
#define BASEDEFS_HH

class BaseDefs {
 public:

  // define labels for bases (sets of intermediate variables)
  enum { VARBASIS_EPXYZ=0, VARBASIS_VXYZ, VARBASIS_TRKNORMAL, NMETASET };

  // max # of variables in the above bases
  enum {MAXINTERVARS=4};

  // maximum number of parameters for a fit object
  enum {MAXPAR = 10};

  // this is used to store how many variables in each base (should be <= maxInter)
  static const int nMetaVars[NMETASET];

};

#endif


