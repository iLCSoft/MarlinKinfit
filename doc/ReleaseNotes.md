# v00-03

J. List
   - added FourJetZHPairing

F. Gaede
   - made compatible with c++11
   - removed -ansi -pedantic -Wno-long-long


# v00-02

 moved example processors now into separate package MarlinKinfitProcessors,
 added ZinvisibleFitObject,
 debugged covariance matrix on fitted 4-momenta. 

# v00-01-05

 adjusted version number in CMakeLists.txt, removed obsolete PConstraint,  
 added directory not_used and moved there all the Track / Vertex fitting stuff which
 had only recently (r4750) been transfered from H1 version and is now being
 overhauled by Daniel Jeans.  

# v00-01-04

 First version after refactoring of FitObject inheritance tree (D. Jeans)
 Minimized duplication of code in derived classes at the price of a
 overall limit on the maximum number of parameters, to be set in BaseDefs.h/cc
 As a consequence, the logic of the copy constructors / assignment operators
 had to be adjusted. 

# v00-01-03

 last version before refactoring of FitObject inheritance tree

