# v00-07

* 2025-10-06 Thomas Madlener ([PR#6](https://github.com/iLCSoft/MarlinKinfit/pull/6))
  - Remove a CMake policy that is no longer necessary due to the updated cmake configuration

* 2025-10-02 Thomas Madlener ([PR#5](https://github.com/iLCSoft/MarlinKinfit/pull/5))
  - Adopt cmake configuration to use the target based style of configuring
  - Raise the minimum version of cmake to 3.23
    - Keep existing behavior by enabling cmake policy CMP0033
  - Fix dependencies to rely on LCIO instead of Marlin (since Marlin is not truly a dependency, but LCIO is).
  - Make ROOT a required dependency since building without ROOT is not possible
  - Export targets for downstream consumption with pure cmake functionality. The exported target is `MarlinKinfit::MarlinKinfit`.

* 2025-10-01 Thomas Madlener ([PR#4](https://github.com/iLCSoft/MarlinKinfit/pull/4))
  - Add a Key4hep based CI workflow
  - Update clicdp nightlies based workflows
  - **Fix many compiler warnings** and cleanup code slightly
    - deleting all copy & move constructors
    - Removing unused variables
    - Intializing member variables
    - Switching from c-style arrays to std::array
    - Replacing enum constants with constexpr static ints
    - Shadowing warnings

# v00-06-01

* 2022-06-28 Thomas Madlener ([PR#2](https://github.com/iLCSoft/MarlinKinfit/pull/2))
  - Migrate CI to github actions and remove travis CI
  - Make doxygen cmake config work with newer cmakes (>= 3.17)

# v00-06

# v00-05

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

