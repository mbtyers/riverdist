% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apps_of_direction.R
\name{plotseq}
\alias{plotseq}
\title{Plot Sequence of Observations}
\usage{
plotseq(
  seqbysurvey,
  type = "boxplot",
  xlab = "",
  ylab = "",
  main = "",
  cex.axisX = 0.8,
  lowerbound = NULL,
  upperbound = NULL,
  boundtype = "negative",
  surveysareDates = F,
  ...
)
}
\arguments{
\item{seqbysurvey}{A matrix returned from \link{mouthdistbysurvey}, 
\link{riverdistanceseq}, or \link{upstreamseq}.}

\item{type}{The type of plot to generate.  Options are 
\code{"boxplot"},\code{"dotplot"},\code{"boxline"},or \code{"dotline"}. 
Defaults to \code{"boxplot"}.}

\item{xlab}{X-axis label}

\item{ylab}{Y-axis label}

\item{main}{Plot title}

\item{cex.axisX}{Character expansion factor for X-axis labels}

\item{lowerbound}{An optional vector of lower survey bounds}

\item{upperbound}{An optional vector of upper survey bounds}

\item{boundtype}{Method of plotting survey bounds.  Options are
\code{"positive"}, \code{"negative"} (default), and \code{"lines"}.}

\item{surveysareDates}{If surveys are in Date format (see \link{as.Date}), a
value of \code{TRUE} allows the x-coordinates points to be spaced apart
according to date, not equidistantly.  Defaults to \code{FALSE}.  Any formatting of 
the survey variable must be done within the original call to \link{mouthdistbysurvey}, 
\link{riverdistanceseq}, or \link{upstreamseq}.  Dates must already be formatted as dates,
or in the form \code{"YYYY-MM-DD"} or \code{"YYYY/MM/DD"}.}

\item{...}{Additional plotting parameters}
}
\description{
Plots the sequence of observations or movements of each individual (given as 
  segment and vertex).  This function is primarily intended for use with 
  \link{mouthdistbysurvey}, but will also work with \link{riverdistanceseq} and 
  \link{upstreamseq}.
}
\note{
Plots are intended as descriptive only.  Any ANOVA-like inference that 
  is suggested from these plots is strongly discouraged.  The user is instead
  advised to use a mixed-effects model or some other inferential tool that 
  accounts for repeated-measures and/or temporal autocorrelation.
}
\examples{
data(Gulk, fakefish)

x <- mouthdistbysurvey(unique=fakefish$fish.id, survey=fakefish$flight.date, 
    seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk)
    
plotseq(seqbysurvey=x)
plotseq(seqbysurvey=x, type="boxline")
plotseq(seqbysurvey=x, type="dotplot")
plotseq(seqbysurvey=x, type="dotline")

plotseq(seqbysurvey=x, type="dotline", surveysareDates=TRUE)

from_upstreamseq <- upstreamseq(unique=fakefish$fish.id, 
   survey=fakefish$flight, seg=fakefish$seg, vert=fakefish$vert, 
   rivers=Gulk)
plotseq(seqbysurvey=from_upstreamseq)
}
\author{
Matt Tyers
}
