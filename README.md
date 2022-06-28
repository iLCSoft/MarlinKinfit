# MarlinKinfit
[![linux](https://github.com/iLCSoft/MarlinKinfit/actions/workflows/linux.yml/badge.svg)](https://github.com/iLCSoft/MarlinKinfit/actions/workflows/linux.yml)
[![Coverity Scan Build Status](https://scan.coverity.com/projects/12361/badge.svg)](https://scan.coverity.com/projects/ilcsoft-marlinkinfit)

Kinematic fitting library for MarlinKinfit. Does not include the Processors, which have been moved
to MarlinKinfitKinfitProcessors.

MarlinKinfit is distributed under the [GPLv3 License](http://www.gnu.org/licenses/gpl-3.0.en.html)

[![License](https://www.gnu.org/graphics/gplv3-127x51.png)](https://www.gnu.org/licenses/gpl-3.0.en.html)


## License and Copyright
Copyright (C), MarlinKinfit Authors

MarlinKinfit is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License long with this program.  If not, see <http://www.gnu.org/licenses/>.


# Kinematic fitting library for Marlin.

Currently two examples are available:

- WW5CFit: This processor performs a 5C fit (px, py, pz, E and equal
  mass constraints) on 4 jet events, as one would do for WW or ZZ
  events. The center of mass energy and the name of the input jet
  collection are steerable.

  The provided example steering assumes that you have a standard DBD
  ILD_o1_v5 DST file at 500 GeV. It runs the FastJetProcessor in exclusive
  kt mode to remove the gammagamma->hadrons overlay and then performs
  the fit on the remaining jets.

- TTBarExample: This processor performs a 6C fit (px, py, pz, E and 2 W
  mass constraints) on 6 jet events, as one could do for ttbar
  events. To reduce the combinatorics it is assumed that 2 jets are
  b-tagged. The center of mass energy and the name of the input jet
  collections are steerable.

  The provided example steering assumes that you have an LCIO file which
  already contains two collections of 4 light jets and 2 b-jets per
  event. If you don't have such a file, you have to call digitization,
  particle flow, jet finder and b-tagging before calling TTBarExample.

In addition, TopTester provides an example how to use TopEventILC,
a toy MC of either fully hadronic or semi-leptonic ttbar events
plus a matching kinematic fit. In the toy MC, the generated four-vectors
of the 6 final-state fermions are smeared with the same resolutions as
used in the fit, and the correct jet pairing is cheated. This comes
handy for testing further developments of the fit engines, or of new
types of fit objects, contraints etc.

Further information about the fit engine and the user interface provided
in MarlinKinfit can be found at
https://www.desy.de/~blist/kinfit/doc/html/

and in the LCNotes LC-TOOL-2009-001 and LC-TOOL-2009-004 available from
http://www-flc.desy.de/lcnotes/

26.11.2014 Jenny List
