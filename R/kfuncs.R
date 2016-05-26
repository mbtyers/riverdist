# kfunc1 <- function(dmat=NULL,seg=NULL,vert=NULL,rivers,envelope=T,envreps=100,...) {
#   if(is.null(dmat) & !is.null(seg) & !is.null(vert)) {
#     dmat <- riverdistancemat(seg,vert,rivers,...=...)
#   }
#   # print("asdf")
#   kdists <- seq(from=0,to=(max(dmat)/4),l=100)
#   kdistavg <- NA*kdists
#   for(i in 1:length(kdists)) {
#     kdistavg[i] <- (sum(dmat<kdists[i])-(dim(dmat)[1]))/(dim(dmat)[1])
#   }
#   plot(kdists,kdistavg,type='l')
#   
#   if(envelope) {
#     allsegs <- allverts <- NULL
#     for(i in 1:length(rivers$lines)) {
#       allsegs <- c(allsegs,rep(i,(dim(rivers$lines[[i]])[1])))
#       allverts <- c(allverts,(1:(dim(rivers$lines[[i]])[1])))
#     }
#     for(i in 1:envreps) {
#       randi <- sample((1:length(allsegs)),size=dim(dmat)[1],replace=F)
#       rsegs <- allsegs[randi]
#       rverts <- allverts[randi]
#       # plot(Gulk)
#       # riverpoints(rsegs,rverts,Gulk)kdists <- seq(from=0,to=(max(dmat)/4),l=100)
#       kdistavg <- NA*kdists
#       dmat1 <- riverdistancemat(rsegs,rverts,rivers)
#       for(i in 1:length(kdists)) {
#         kdistavg[i] <- (sum(dmat1<kdists[i])-(dim(dmat)[1]))/(dim(dmat)[1])
#       }
#       lines(kdists,kdistavg,col="grey90")
#     }
#   }
# }
# # 
# # kfunc(seg=smallset$seg,vert=smallset$vert,rivers=Gulk)
# # kfunc(seg=fakefish$seg,vert=fakefish$vert,rivers=Gulk,envreps=1)
# 
# 
# kfunc2 <- function(seg,vert,survey=NULL,rivers,envelope=F,envreps=100,scale=F,...) {
#   if(is.null(survey)) survey <- 1
#   dmats <- list()
#   i <- 1
#   for(surveyi in sort(unique(survey))) {
#     dmats[[i]] <- riverdistancemat(seg=seg[survey==surveyi],vert=vert[survey==surveyi],rivers=rivers)
#     i <- i+1
#   }
#   # print("asdf")
#   maxes <- NA
#   for(i in 1:length(unique(survey))) maxes[i] <- max(unlist(dmats[[i]]))
#   kdists <- seq(from=0,to=min(maxes),l=100)
#   kdistavg <- list()
#   i <- 1
#   for(surveyi in sort(unique(survey))) {
#     kdistavg[[i]] <- NA*kdists
#     for(j in 1:length(kdists)) {
#       kdistavg[[i]][j] <- (sum(dmats[[i]]<kdists[j])-(dim(dmats[[i]])[1]))/(dim(dmats[[i]])[1])
#       if(scale) kdistavg[[i]][j] <- 100*kdistavg[[i]][j]/(dim(dmats[[i]])[1]-1)
#     }
#     i <- i+1
#   }
#   ylab <- ifelse(scale,"Percent of points within distance","Number of points within distance")
#   xlab <- "Distance"
#   plot(NA,xlim=c(0,1.1*max(kdists)),ylim=c(0,max(unlist(kdistavg))),xlab=xlab,ylab=ylab)
#   cols <- adjustcolor(rainbow(length(unique(survey)),start=0,end=.7),red.f=.9,green.f=.9,blue.f=.9,alpha.f=.8)
#   for(i in 1:length(unique(survey))) {
#     lines(kdists,kdistavg[[i]],col=cols[i],lwd=2)
#     text(max(kdists),max(kdistavg[[i]]),sort(unique(survey))[i],pos=4,col=cols[i])
#   }
#   
#   if(envelope) {
#     allsegs <- allverts <- NULL
#     for(i in 1:length(rivers$lines)) {
#       allsegs <- c(allsegs,rep(i,(dim(rivers$lines[[i]])[1])))
#       allverts <- c(allverts,(1:(dim(rivers$lines[[i]])[1])))
#     }
#     for(i in 1:envreps) {
#       randi <- sample((1:length(allsegs)),size=dim(dmat)[1],replace=F)
#       rsegs <- allsegs[randi]
#       rverts <- allverts[randi]
#       # plot(Gulk)
#       # riverpoints(rsegs,rverts,Gulk)kdists <- seq(from=0,to=(max(dmat)/4),l=100)
#       kdistavg <- NA*kdists
#       dmat1 <- riverdistancemat(rsegs,rverts,rivers)
#       for(i in 1:length(kdists)) {
#         kdistavg[i] <- (sum(dmat1<kdists[i])-(dim(dmat)[1]))/(dim(dmat)[1])
#       }
#       lines(kdists,kdistavg,col="grey90")
#     }
#   }
#   
#   
# }
# kfunc2(seg=fakefish$seg,vert=fakefish$vert,rivers=Gulk,survey=fakefish$flight,scale=T)
# 
# par(mfrow=c(2,5))
# plot(fakefish_density)
# 
# 
# setwd("~/Kusko_burbot")
# load(file="kusko.Rdata")
# bb11 <- subset(bb_pts,Date>=as.Date("2011-10-24")&Date<=as.Date("2012-08-10"))
# summary(as.factor(bb11$Date))
# 
# par(mfrow=c(1,1))
# kfunc2(seg=bb11$seg,vert=bb11$vert,rivers=kusko,survey=bb11$Date,scale=T)
# kfunc2(seg=bb_pts$seg,vert=bb_pts$vert,rivers=kusko,survey=bb_pts$Date,scale=T)
# kfunc2(seg=ag_pts$seg,vert=ag_pts$vert,rivers=kusko,survey=ag_pts$Date,scale=T)



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
#'   envelope (see note.)  Defaults to \code{TRUE}.
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
#' @param ... Additional plotting parameters.
#' @note K-function envelopes for each survey are constructed by bootstrapping 
#'   all within-survey distances, that is, the distances between all individuals
#'   within each survey, for all surveys.  This results in a confidence envelope
#'   under the assumption that spacing is independent of survey; therefore a
#'   survey K-function outside the envelope provides evidence of clustering or 
#'   dispersal in that survey that is outside the typical range.
#'   
#'   A K-function above the envelope for a given distance range provides
#'   evidence of a greater number of individuals than expected at that distance
#'   range (clustering); A K-function below the envelope for a given distance
#'   range provides evidence of a smaller number of individuals than expected at
#'   that distance range (dispersal).
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
kfunc <- function(seg,vert,survey=NULL,rivers,lwd=2,envelope=TRUE,envreps=1000,envcol="grey80",envborder=NA,maxdist=NULL,xlab="Distance",ylab="% within",showN=TRUE,whichplots=NULL,...) {
  if(is.null(survey)) survey <- " "
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
      # if(scale) kdistavg[[i]][j] <- 100*kdistavg[[i]][j]/(dim(dmats[[i]])[1]-1)
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
#     for(k in 1:envreps) {
#       dmatallboot <- sample(dmatall,length(dmatall),replace=T)
#       for(j in 1:length(kdists)) {
#         kdistavgavgboot[k,j] <- 100*(sum(dmatallboot<kdists[j]))/(length(dmatallboot))
#       }
#     }
  }
  kdistavgbootlower <- kdistavgbootupper <- NA*kdists
#   for(j in 1:length(kdists)) {
#     kdistavgbootlower[j] <- quantile(kdistavgavgboot[,j],0.025)
#     kdistavgbootupper[j] <- quantile(kdistavgavgboot[,j],0.975)
#   }
  for(i in whichplots) {
    pmain <- ifelse(showN,paste0(sort(unique(survey))[i]," (n=",length(seg[survey==sort(unique(survey))[i]]),")"),sort(unique(survey))[i])
    plot(kdists,kdistavg[[i]],col=1,lwd=lwd,xlim=c(0,1.1*max(kdists)),ylim=c(0,max(unlist(kdistavg))),xlab=xlab,ylab=ylab,type='l',main=pmain)
    if(envelope) {
      ######
      for(k in 1:envreps) {
        dmatallboot <- sample(dmatall,sum(upper.tri(dmats[[i]])),replace=T)
        for(j in 1:length(kdists)) {
          kdistavgavgboot[k,j] <- 100*(sum(dmatallboot<kdists[j]))/(length(dmatallboot))
        }
      }
#       dmatallboot <- matrix(sample(dmatall,(envreps*sum(upper.tri(dmats[[i]]))),replace=T),nrow=envreps,ncol=sum(upper.tri(dmats[[i]])))
#       for(j in 1:length(kdists)) {
#         kdistavgavgboot[,j] <- 100*rowMeans(dmatallboot<kdists[j])
#       }
      for(j in 1:length(kdists)) {
        kdistavgbootlower[j] <- quantile(kdistavgavgboot[,j],0.025)
        kdistavgbootupper[j] <- quantile(kdistavgavgboot[,j],0.975)
      }
      ######
      polygon(c(kdists,kdists[length(kdists):1]),c(kdistavgbootlower,kdistavgbootupper[length(kdists):1]),col=envcol,border=envborder)
      lines(kdists,kdistavg[[i]],col=1,lwd=lwd)
    }
    lines(kdists,kdistavgavg,lty=2)
  }
}
# 
# a <- Sys.time()
# par(mfrow=c(3,4))
# kfunc3(seg=fakefish$seg,vert=fakefish$vert,rivers=Gulk,survey=fakefish$flight,envelope=T,envreps=1000,maxdist=200000,lwd=2)
# Sys.time()-a
# 
# a <- Sys.time()
# par(mfrow=c(5,5))
# kfunc3(seg=bb_pts$seg,vert=bb_pts$vert,rivers=kusko,survey=bb_pts$Date,scale=T,envelope=T,envreps=1000)
# Sys.time()-a
# 
# bb_pts2 <- subset(bb_pts,(Date!=as.Date("2011-10-01")&Date!=as.Date("2012-10-01")))
# a <- Sys.time()
# par(mfrow=c(4,5))
# kfunc3(seg=bb_pts2$seg,vert=bb_pts2$vert,rivers=kusko,survey=bb_pts2$Date,scale=T,envelope=T,envreps=1000)
# Sys.time()-a