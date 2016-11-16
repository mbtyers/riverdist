\name{riverdist-package}
\alias{riverdist-package}
\alias{riverdist}
\docType{package}
\title{
River Network Distance Computation and Applications
}
\description{
Reads river network shape files and computes network distances.
    Also included are a variety of computation and graphical tools designed 
    for fisheries telemetry research, such as minimum home range, kernel density 
    estimation, and clustering analysis using empirical k-functions with 
    a bootstrap envelope.  Tools are also provided for editing the river 
    networks, meaning there is no reliance on external software.
}
\details{
\tabular{ll}{
Package: \tab riverdist\cr
Type: \tab Package\cr
Version: \tab 0.12.2.9000\cr
Date: \tab 2016-08-11\cr
License: \tab GPL-2\cr
}
The riverdist package provides tools for distance calculation along a river network.  The river network is imported from a projected shapefile.  Spatial point data may be imported from a shapefile as well, or directly from coordinates.

Some basic formatting of the river shapefile may be necessary.  If available, formatting in a geographic information system (GIS) prior to importing into R is recommended (projecting, spatial trimming to the study area, and possibly dissolving river segments), but the riverdist package and its dependencies also include tools for accomplishing the necessary formatting within R.
}
\author{
Matt Tyers

Maintainer: Matt Tyers <matt.tyers@alaska.gov>
}


