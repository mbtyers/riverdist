# riverdist 

### River Network Distance Computation and Applications

The 'riverdist' package is intended as a free and readily-available resource for distance calculation along a river network.  This package was written with fisheries research in mind, but could be applied to other fields.  The 'riverdist' package builds upon the functionality of the 'sp' and 'rgdal' packages, which provide the utility of reading GIS shapefiles into the R environment.  What 'riverdist' adds is the ability to treat a linear feature as a connected network, and to calculate travel routes and travel distances along that network.

### Commonly-used functions

* `line2network()` imports a river shapefile, and calculates topologies to create a connected river network.

* `cleanup()` calls a sequence of editing functions on the resulting river network to facilitate performance.  The editing functions are also available by themselves.

* `plot()` when used with a river network object produces a simple map of the network, with segments labeled and differentiated by color or line type.

* `xy2segvert()` converts a set of X-Y coordinates to river network coordinates by "snapping" each point to the nearest river vertex.

* `riverdistance()`, `riverdirection()`, and `upstream()` return the network distance, the travel direction (upstream or downstream), and the directional (upstream) distance between two river locations, respectively.  Options are included for different handling of locations that are flow-connected or flow-separate, as well as net directional distance when locations are flow-separate.

Several automated analyses are built in.  In most cases, there is a direction or directional distance equivalent.

* `homerange()` returns the minimum observed home range for each individual in a data set
* `riverdistanceseq()` and `riverdistanceseqobs()` return different forms of matrices of pairwise network distances between observations of each individual in a dataset
* `riverdistancemat()` returns a matrix of network distances between all observations in a dataset
* `riverdistancetofrom()` returns a matrix of network distances between two datasets, 
* `mouthdistobs()` returns a matrix of distances between each observation and the mouth of the river network, with rows corresponding to unique individual, and columns corresponding to unique survey 

### Installation

The 'riverdist' package is currently available on Github, and can be installed with the following code:
`install.packages("devtools",dependencies=T)`

`install_github("mbtyers/riverdist")`

### Issues

A major dependency of the 'riverdist' package is the 'rgdal' package, which allows importing shapefiles.  For installation on a non-windows machine, please refer to the SystemRequirements given at https://cran.r-project.org/web/packages/rgdal/index.html

'riverdist' passes R CMD check on both Windows 7 and OS X 10.10.3, but package testing using Travis-CI still gives a Build-Error.  This is not a 'riverdist' issue and does not affect the performance of 'riverdist' within R.  This is a testing issue, since I have not been able to configure my .travis.yml file such that the libraries 'rgdal' needs are loaded.

[![Travis-CI Build Status](https://travis-ci.org/mbtyers/riverdist.svg?branch=master)](https://travis-ci.org/mbtyers/riverdist)
