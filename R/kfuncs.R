#' Plotting K-functions for a Set of Surveys
#' @description Plots K-functions for locations in each of a set of surveys.  In
#'   this implementation, this can be interpreted as the proportion of 
#'   additional fish within a given distance.  This will increase as a function 
#'   of distance, and may provide evidence of clustering or dispersion features,
#'   particularly if the envelope is used.
#' @param seg A vector of river locations (segment)
#' @param vert A vector of river locations (vertex)
#' @param survey A vector of survey IDs corresponding to the values of 
#'   \code{seg} and \code{vert}.  Defaults to \code{NULL}.  If this argument is 
#'   used, K-functions will be calculated for each unique survey, and separate 
#'   plots will be produced.
#' @param rivers The river network object to use
#' @param lwd Line width used for plotting.  Defaults to 2.
#' @param envelope Whether to construct and display a 95 percent confidence 
#'   envelope (see note.)  Defaults to \code{TRUE} if \code{survey} is specified, and is automatically \code{FALSE} otherwise.
#' @param envreps Number of bootstrap replicates to use for envelope 
#'   calculation.  Defaults to 1000.
#' @param envcol Color to use for envelope plotting.  Defaults to 
#'   \code{"grey80"}.
#' @param envborder Border color to use for envelope plotting.  Defaults to 
#'   \code{NA}, which will result in no border being plotted.
#' @param maxdist Maximum distance (x-axis value) for plotting.  The default 
#'   value (\code{NULL}) will result in an appropriate value being chosen.
#' @param xlab X-coordinate label for plotting
#' @param ylab Y-coordinate label for plotting
#' @param showN Whether to show the sample size for each survey in each plot 
#'   title.  Defaults to \code{TRUE}.
#' @param whichplots A vector of plots to produce, if multiple plots are 
#'   produced.  For example, specifying \code{whichplot=c(2,3,4)} will result in
#'   only the second, third, and fourth plots of the sequence being produced. 
#'   Accepting the default (\code{NULL}) will result in all plots being 
#'   produced.
#' @param returnoutput Whether to return output instead of producing a plot.  Defaults to \code{FALSE}.
#' @param ... Additional plotting parameters.
#' @note K-function envelopes for each survey are constructed by bootstrapping 
#'   all within-survey distances, that is, the distances between all individuals
#'   within each survey, for all surveys.  This results in a confidence envelope
#'   under the assumption that spacing is independent of survey; therefore a
#'   survey K-function outside the envelope provides evidence of clustering or 
#'   dispersal in that survey that is outside the typical range.  An envelope is not available if only one survey is plotted.
#'   
#'   A K-function above the envelope for a given distance range provides
#'   evidence of a greater number of individuals than expected at that distance
#'   range (clustering); A K-function below the envelope for a given distance
#'   range provides evidence of a smaller number of individuals than expected at
#'   that distance range (dispersal).
#' @note This function is distance-computation intensive, and will be extremely slow-running if a river network is used that does not have segment routes and/or distance lookup tables for fast distance computation.  See \link{buildsegroutes} and/or \link{buildlookup} for more information.
#' @author Matt Tyers
#' @importFrom stats quantile
#' @examples
#' data(Gulk, fakefish)
#' 
#' # # 10 plots will be created - recommend calling
#' # # par(mfrow=c(3,4))
#' 
#' kfunc(seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk, survey=fakefish$flight,
#' envreps=100, maxdist=200000)
#' 
#' # # This shows relatively high amounts of clustering for surveys 1 and 8,
#' # # and relatively high amounts of dispersal in surveys 5 and 6.
#' 
#' 
#' # # plotting the survey locations that led to this calculation, for comparison
#' 
#' # # 10 plots will be created - recommend calling
#' # # par(mfrow=c(3,4))
#' for(i in 1:10) {
#'   plot(x=Gulk, segmentnum=FALSE, color=FALSE, main=i)
#'   riverpoints(seg=fakefish$seg[fakefish$flight==i], 
#'   vert=fakefish$vert[fakefish$flight==i], rivers=Gulk, col=2, pch=15)
#' }
#' @export
kfunc <- function(seg,vert,survey=NULL,rivers,lwd=2,envelope=TRUE,envreps=1000,envcol="grey80",envborder=NA,maxdist=NULL,xlab="Distance",ylab="% within",showN=TRUE,whichplots=NULL,returnoutput=FALSE,...) {
  if(is.null(survey)) survey <- " "
  if(length(unique(survey))==1) envelope <- F
  if(is.null(whichplots)) whichplots <- 1:length(unique(survey))
  dmats <- list()
  i <- 1
  for(surveyi in sort(unique(survey))) {
    dmats[[i]] <- riverdistancemat(seg=seg[survey==surveyi],vert=vert[survey==surveyi],rivers=rivers)
    i <- i+1
  }
  maxes <- NA
  for(i in 1:length(unique(survey))) maxes[i] <- max(unlist(dmats[[i]]))
  maxdist1 <- ifelse(is.null(maxdist),min(maxes),maxdist)
  kdists <- seq(from=0,to=maxdist1,l=100)
  kdistavg <- list()
  i <- 1
  for(surveyi in sort(unique(survey))) {
    kdistavg[[i]] <- NA*kdists
    for(j in 1:length(kdists)) {
      kdistavg[[i]][j] <- (sum(dmats[[i]]<kdists[j])-(dim(dmats[[i]])[1]))/(dim(dmats[[i]])[1])*100/(dim(dmats[[i]])[1]-1)
    }
    i <- i+1
  }
  
  kdistavgavg <- NA*kdists
  dmatall <- NULL
  dim1 <- 0
  for(i in 1:length(kdistavg)) {
    dim1 <- dim1+dim(dmats[[i]])[1]
    dmatall <- c(dmatall,dmats[[i]][upper.tri(dmats[[i]],diag=F)])
  }
  dim2 <- dim1-length(kdistavg)
  for(j in 1:length(kdists)) {
    kdistavgavg[j] <- 100*(sum(dmatall<kdists[j]))/(length(dmatall))
  }
  if(envelope) {
    kdistavgavgboot <- matrix(NA,nrow=envreps,ncol=length(kdists))
  }
  kdistavgbootlower <- kdistavgbootupper <- NA*kdists
  if(returnoutput) {
    output <- list()
    output$lines <- list()
    output$env_low <- list()
    output$env_high <- list()
    output$dists <- list()
  }
  for(i in whichplots) {
    if(!returnoutput){
      pmain <- ifelse(showN,paste0(sort(unique(survey))[i]," (n=",length(seg[survey==sort(unique(survey))[i]]),")"),sort(unique(survey))[i])
      plot(kdists,kdistavg[[i]],col=1,lwd=lwd,xlim=c(0,1.1*max(kdists)),ylim=c(0,max(unlist(kdistavg))),xlab=xlab,ylab=ylab,type='l',main=pmain)
    }
    if(envelope) {
      for(k in 1:envreps) {
        dmatallboot <- sample(dmatall,sum(upper.tri(dmats[[i]])),replace=T)
        for(j in 1:length(kdists)) {
          kdistavgavgboot[k,j] <- 100*(sum(dmatallboot<kdists[j]))/(length(dmatallboot))
        }
      }
      for(j in 1:length(kdists)) {
        kdistavgbootlower[j] <- quantile(kdistavgavgboot[,j],0.025)
        kdistavgbootupper[j] <- quantile(kdistavgavgboot[,j],0.975)
      }
      if(!returnoutput) {
        polygon(c(kdists,kdists[length(kdists):1]),c(kdistavgbootlower,kdistavgbootupper[length(kdists):1]),col=envcol,border=envborder)
        lines(kdists,kdistavg[[i]],col=1,lwd=lwd)
      }
      if(returnoutput) {
        output$lines[[i]] <- kdistavg[[i]]
        output$env_low[[i]] <- kdistavgbootlower
        output$env_high[[i]] <- kdistavgbootupper
      }
    }
    if(!returnoutput) lines(kdists,kdistavgavg,lty=2)
  }
  if(returnoutput) {
    output$dists <- kdists
    return(output)
  }
}