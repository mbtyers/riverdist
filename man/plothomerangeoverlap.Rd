% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apps_of_riverdistance.R
\name{plothomerangeoverlap}
\alias{plothomerangeoverlap}
\title{Plot Home Range Overlap}
\usage{
plothomerangeoverlap(x, lwd = 3, maxlwd = 10, col = 4, ...)
}
\arguments{
\item{x}{An object returned from \link{homerange}.}

\item{lwd}{Minimum line width to use, defaults to 3.}

\item{maxlwd}{Maximum line width to use, defaults to 10.}

\item{col}{Line color to use, defaults to \code{"blue"}.}

\item{...}{Additional plotting parameters, see \link{plot.rivernetwork}.}
}
\description{
Produces a plot of the overlap of the minimum observed home range for multiple
  observations of each individual fish, with line thickness illustrating the respective number of individuals' homeranges represented.
}
\examples{
data(Gulk, fakefish)
ranges <- with(fakefish, homerange(unique=fish.id, survey=flight, seg=seg, vert=vert, rivers=Gulk))
ranges

# 19 plots will be produced, recommend calling par(mfrow=c(4,5))
plot(ranges)
plot(ranges,cumulative=TRUE,label=TRUE)

homerangeoverlap(ranges)

plothomerangeoverlap(ranges)
with(fakefish, riverpoints(seg=seg, vert=vert, rivers=Gulk))
}
\seealso{
\link{homerange}, \link{plot.homerange}, \link{homerangeoverlap}
}
\author{
Matt Tyers
}
