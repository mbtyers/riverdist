# riverdist 

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/riverdist)](https://cran.r-project.org/package=riverdist)

### River Network Distance Computation and Applications

The 'riverdist' package is intended as a free and readily-available resource for distance calculation along a river network.  This package was written with fisheries research in mind, but could be applied to other fields.  The 'riverdist' package builds upon the functionality of the 'sp' and 'rgdal' packages, which provide the utility of reading GIS shapefiles into the R environment.  What 'riverdist' adds is the ability to treat a linear feature as a connected network, and to calculate travel routes and travel distances along that network.

### Commonly-used functions

* `line2network()` imports a river shapefile, and calculates topologies to create a connected river network.

* `cleanup()` calls a sequence of editing functions on the resulting river network to facilitate performance.  The editing functions are also available by themselves.

* `plot()` when used with a river network object produces a simple map of the network, with segments labeled and differentiated by color or line type.

* `xy2segvert()` converts a set of X-Y coordinates to river network coordinates by "snapping" each point to the nearest river vertex.  `pointshp2segvert()` does the same, with the input being a point shapefile.

* `riverdistance()`, `riverdirection()`, and `upstream()` return the network distance, the travel direction (upstream or downstream), and the directional (upstream) distance between two river locations, respectively.  Options are included for different handling of locations that are flow-connected or flow-separate, as well as net directional distance when locations are flow-separate.

Several automated analyses are built in.  In most cases, there is a direction or directional distance equivalent.

* `homerange()` returns the minimum observed home range for each individual in a data set.
* `riverdistanceseq()` and `riverdistanceseqbysurvey()` return different forms of matrices of pairwise network distances between observations of each individual in a dataset.
* `riverdistancemat()` returns a matrix of network distances between all observations in a dataset.
* `riverdistancetofrom()` returns a matrix of network distances between two datasets.
* `mouthdistbysurvey()` returns a matrix of distances between each observation and the mouth of the river network, with rows corresponding to unique individual, and columns corresponding to unique survey.

Summaries and plots are also available at the dataset level, in addition to individuals, which is likely to be much more useful to analysis.

* `makeriverdensity()` calculates a kernel density object which can be plotted with `plot()` to create a kernel density map.  Depending on the usage of `makeriverdensity()`, this may be a sequence of maps.
* `kfunc()` provides plotting of empirical k-functions for each survey event, giving evidence of clustering or dispersal behavior.
* `plotseq()` produces a plot of a distance sequence such as that returned from `mouthdistbysurvey()` providing plots of overall distance or upriver position for each survey event.
* `matbysurveylist()` produces a list of matrices of distances or upstream distances between all survey events, for each individual.  This can be plotted using `plotmatbysurveylist()`, creating a summary plot for all individuals.

### Installation

Version 0.14.0 of the 'riverdist' package is available on CRAN.

The development version is currently available on Github, and can be installed in R with the following code:

`install.packages("devtools",dependencies=T)`

`devtools::install_github("mbtyers/riverdist")`

### Issues

A major dependency of the 'riverdist' package is the 'rgdal' package, which allows importing shapefiles.  For installation on a non-windows machine, please refer to the SystemRequirements given at https://cran.r-project.org/web/packages/rgdal/index.html