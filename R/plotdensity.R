#' Calculate Kernel Density Using River Distance
#' @description Uses spatial point data (segment and vertex) to calculate a 
#'   kernel density object to use in the output class plotting method,\link{plot.riverdensity}. Scaled kernel 
#'   density is calculated at approximately regularly-spaced locations, with 
#'   spacing specified by the user.
#'   
#'   If an argument is used in the \code{survey} field, kernel densities will be
#'   calculated for each unique value of \code{survey}, resulting in a separate 
#'   plot for each.
#'   
#'   The purpose of this function is to generate a kernel density object to plot
#'   using plot(), see \link{plot.riverdensity}.
#' @note It is likely that calculation will be very slow.  Use of this function 
#'   with a river network for which segment routes has not yet been calculated 
#'   is not recommended.
#' @param seg A vector of river locations (segment)
#' @param vert A vector of river locations (vertex)
#' @param rivers The river network object to use
#' @param survey A vector of survey IDs corresponding to the values of 
#'   \code{seg} and \code{vert}.  If this argument is used, kernel densities 
#'   will be calculated for each unique survey, and separate plots will be 
#'   produced.
#' @param kernel The type of density kernel to use.  Allowed types are 
#'   \code{"gaussian"} (normal) and \code{"rect"} (rectangular, giving simple 
#'   density).  Defaults to \code{"gaussian"}.
#' @param bw The kernel bandwidth to use.  If \code{kernel} is set to 
#'   \code{"gaussian"}, this provides the standard deviation of the gaussian 
#'   (normal) kernel to use.  If \code{kernel} is set to \code{"rect"}, this 
#'   provides the half-width of the rectangular kernel, or the distance to use 
#'   in simple density.  Accepting the default (\code{NULL}) will result in the 
#'   function determining a value to use, based on the total length of the river
#'   network and the value of the \code{resolution} argument.
#' @param resolution The approximate spacing of the river locations used for
#'   kernel density calculation.  Accepting the default (\code{NULL}) will
#'   result in the function determining a value to use, based on the total
#'   length of the river network.
#' @note This function is distance-computation intensive, and may be slow-running if a river network is used that does not have segment routes and/or distance lookup tables for fast distance computation.  See \link{buildsegroutes} and/or \link{buildlookup} for more information.
#' @return A river density object, see \link{riverdensity-class}.
#' @seealso \link{plot.riverdensity}, \link{plotriverdensitypoints}
#' @author Matt Tyers
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom stats dnorm
#' @examples
#' data(Gulk, fakefish)
#' 
#' Gulk_dens <- makeriverdensity(seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk, 
#'   survey=fakefish$flight.date)
#'   
#' # # 10 plots will be created, recommend calling par(mfrow=c(2,5))
#' plot(x=Gulk_dens)
#' @export
makeriverdensity <- function(seg,vert,rivers,survey=NULL,kernel="gaussian",bw=NULL,resolution=NULL) {
  if(is.null(resolution)) resolution <- sum(rivers$lengths)/500
  if(is.null(bw)) bw <- 10*resolution     
  if(is.null(survey)) survey <- 0
  if(!is.factor(survey)) survey <- as.factor(survey)
  
  # prediction grid
  densverts <- list()
  endptverts <- list()
  for(segi in 1:length(rivers$lines)) {
    segilength <- dim(rivers$lines[[segi]])[1]
    endptlengths <- c(seq(from=0,to=rivers$lengths[segi],by=resolution),(rivers$lengths[segi]))      
    denslengths <- seq(from=(min(resolution,rivers$lengths[segi])/2),to=rivers$lengths[segi],by=resolution)
    if(length(endptlengths)-length(denslengths)==2) {
      denslengths[length(denslengths)+1] <- rivers$lengths[segi]  # this is a hack but shouldn't be too important
    }
    densverts[[segi]] <- endptverts[[segi]] <- NA
    for(endpti in 1:length(endptlengths)) {
      absdiffs <- abs(rivers$cumuldist[[segi]]-endptlengths[endpti])
      endptverts[[segi]][endpti] <- which(absdiffs==min(absdiffs,na.rm=T))[1]
    }
    for(densi in 1:length(denslengths)) {
      absdiffs <- abs(rivers$cumuldist[[segi]]-denslengths[densi])
      densverts[[segi]][densi] <- which(absdiffs==min(absdiffs,na.rm=T))[1]
    }
  }
  
  densities <- list()
  
  del <- 3*bw  
  
  delind <- floor(del/resolution)
  
  isurvey <- 1
  if(interactive()) pb <- txtProgressBar(style=3)          
  # for(surveyi in sort(unique(survey))) {         
  for(surveyi in levels(survey)) {
    densities[[isurvey]] <- list()
    for(segi in 1:length(rivers$lines)) {
      segilength <- dim(rivers$lines[[segi]])[1]
      whichcoarse <- seq(from=1,to=length(densverts[[segi]]),by=delind)
      if(whichcoarse[length(whichcoarse)] != length(densverts[[segi]])) {
        whichcoarse[length(whichcoarse)+1] <- length(densverts[[segi]])
      }
      densvertscoarse <- densverts[[segi]][whichcoarse]
      tofromicoarse <- riverdistancetofrom(seg1=rep(segi,length(densvertscoarse)),vert1=densvertscoarse,seg2=seg[survey==surveyi],vert2=vert[survey==surveyi],rivers=rivers)
      
      densvertsTF <- rep(F,length(densverts[[segi]]))   
      for(i in 1:length(densvertscoarse)) {
        if(min(tofromicoarse[i,])<=del) {
          if(i>1 & i<length(densvertscoarse)) densvertsTF[(whichcoarse[i-1]):(whichcoarse[i+1])] <- T
          if(i==1 & i<length(densvertscoarse)) densvertsTF[(whichcoarse[i]):(whichcoarse[i+1])] <- T
          if(i==length(densvertscoarse) & 1<length(densvertscoarse)) densvertsTF[(whichcoarse[i-1]):(length(densvertsTF))] <- T
          if(i==1 & i==length(densvertscoarse)) densvertsTF[(whichcoarse[i]):(length(densvertsTF))] <- T
        }
      }
      
      if(any(densvertsTF)) {
        tofromi <- riverdistancetofrom(seg1=rep(segi,length(densverts[[segi]][densvertsTF])),vert1=densverts[[segi]][densvertsTF],seg2=seg[survey==surveyi],vert2=vert[survey==surveyi],rivers=rivers)
        if(kernel=="gaussian") {
          densitiesi <- dnorm(x=tofromi,mean=0,sd=bw)
          densitiesi[densitiesi<dnorm(x=2*bw,mean=0,sd=bw)] <- 0
        }
        if(kernel=="rect") densitiesi <- (tofromi<bw)
      }
      densitiesii <- rep(NA,length(densvertsTF))
      if(!any(densvertsTF)) densitiesi <- matrix(0)
      densitiesii[densvertsTF] <- unname(rowSums(densitiesi))
      densitiesii[!densvertsTF] <- 0
      densities[[isurvey]][[segi]] <- densitiesii
      prop_done <- ((isurvey-1)/length(unique(survey))) + (segi/length(rivers$lines)/length(unique(survey)))
      if(interactive()) setTxtProgressBar(pb=pb, value=prop_done)
    }
    isurvey <- isurvey+1
  }
  
  out <- list(densities=densities, endptverts=endptverts, densverts=densverts, pointsegs=seg, pointverts=vert, survey=survey, rivers=rivers)
  class(out) <- "riverdensity"
  return(out)
}


#' Plot Kernel Density Using River Distance
#' @description Produces a kernel density plot from a kernel density object 
#'   created by \link{makeriverdensity}.
#'   
#'   If the kernel density object includes densities from multiple surveys, a 
#'   new plot will be created for each survey.
#'   
#'   Densities can be displayed using either line widths, color, or both.
#'   
#'   The relative densities that are displayed in the plot are calculated 
#'   according to the form (density/maxdensity)^pwr, with the value of pwr set 
#'   by the \code{pwr} argument.  Setting \code{pwr} to a value less than 1 
#'   allows smaller values to be more visible on the plot.
#' @param x A river density object created by \link{makeriverdensity}.
#' @param whichplots A vector of plots to produce, if multiple plots are 
#'   produced.  For example, specifying \code{whichplot=c(2,3,4)} will result in
#'   only the second, third, and fourth plots of the sequence being produced. 
#'   Accepting the default (\code{NULL}) will result in all plots being 
#'   produced.
#' @param points Whether to add the points used for density calculation. 
#'   Defaults to \code{TRUE}.
#' @param bycol Whether to use a color ramp to show densities.  Defaults to 
#'   \code{TRUE}.
#' @param bylwd Whether to use line thickness to show densities.  Defaults to 
#'   \code{TRUE}.
#' @param maxlwd The maximum line width to use if \code{bylwd} is set to 
#'   \code{TRUE}.  Defaults to 10.
#' @param pwr The power to use in the nonlinear transformation calculating the 
#'   relative density values to be displayed (see above.)  Defaults to 0.7.
#' @param scalebyN Whether to display relative density values scaled by sample
#'   size.  Specifying \code{scalebyN=TRUE} will show larger density values
#'   associated with surveys with more points, and may be more appropriate for
#'   displaying total density.  Specifying \code{scalebyN=FALSE} will allow
#'   surveys with smaller sample sizes to be plotted with similar density values
#'   as those with larger sample sizes, and may be more appropriate for
#'   displaying relative density.  Defaults to \code{TRUE}.
#' @param ramp The color ramp used to display densities if \code{bycol} is set 
#'   to \code{TRUE}.  Allowed values are \code{"grey"} (or \code{"gray"}), 
#'   \code{"red"}, \code{"green"}, \code{"blue"}, \code{"heat"}, 
#'   \code{"stoplight"}, and \code{"rainbow"}.  Defaults to \code{"grey"}.
#' @param lwd The line width to use for background lines if \code{bylwd} is set 
#'   to \code{TRUE}, or all lines if \code{bylwd} is set to \code{FALSE}. 
#'   Defaults to 1.
#' @param linecol The line color to use for background lines if \code{bylwd} is 
#'   set to \code{TRUE}.  Defaults to \code{"black"}.
#' @param denscol The line color to use for showing density if \code{bylwd} is 
#'   set to \code{TRUE}.  Defaults to \code{"black"}.
#' @param alpha The opacity value for lines.  This could potentially allow 
#'   multiple density plots to be overlayed with different colors.
#' @param dark A color-saturation adjustment, with values in [0,1].  A value of 
#'   1 uses the true colors, and a value less than 1 will render the colors as 
#'   slightly darker (less saturated), which may be appear better.  Defaults to 
#'   1.
#' @param showN Whether to automatically include the number of points used as 
#'   part of the plot title(s).
#' @param main Plot title(s), either given as a single text string which is 
#'   repeated if multiple plots are produced, or a vector of text strings (one 
#'   for each plot produced).  If multiple plots are produced (resulting from 
#'   multiple surveys), accepting the default (\code{NULL}) will result in each 
#'   unique value of survey being used as the plot title, along with the sample 
#'   size if \code{showN} is set to \code{TRUE}.
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param add Whether to produce a new plot (\code{FALSE}), or add to an 
#'   existing plot (\code{TRUE}).  Defaults to \code{FALSE}.
#' @param scalebar Whether to add a scale bar to plot(s).  Defaults to \code{TRUE}.
#' @param ... Additional plotting parameters.
#' @seealso \link{makeriverdensity}, \link{plotriverdensitypoints}
#' @method plot riverdensity
#' @aliases plotriverdensity
#' @author Matt Tyers
#' @importFrom grDevices grey
#' @importFrom grDevices rainbow
#' @importFrom grDevices heat.colors
#' @importFrom graphics rect
#' @importFrom grDevices adjustcolor
#' @examples
#' data(Gulk, fakefish)
#' 
#' Gulk_dens <- makeriverdensity(seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk, 
#'   survey=fakefish$flight.date)
#'   
#' # # 10 plots will be created, recommend calling par(mfrow=c(2,5))
#' plot(x=Gulk_dens)
#' @export
plot.riverdensity <- function(x,whichplots=NULL,points=TRUE,bycol=TRUE,bylwd=TRUE,maxlwd=10,pwr=0.7,scalebyN=TRUE,ramp="grey",lwd=1,linecol="black",denscol="black",alpha=1,dark=1,showN=TRUE,main=NULL,xlab="",ylab="",add=FALSE,scalebar=TRUE,...) {
  if(class(x)!="riverdensity") stop("Argument x must be an object returned from makeriverdensity().")
  if(dark>1 | dark<0) dark <-1
  if(alpha>1 | alpha<0) alpha <-1
  densities <- x$densities
  endptverts <- x$endptverts
  densverts <- x$densverts
  seg <- x$pointsegs
  vert <- x$pointverts
  survey <- x$survey
  if(!is.factor(survey)) survey <- as.factor(survey)
  rivers <- x$rivers
  if(length(main)==1) main <- rep(main,length(unique(survey)))
  # if(is.null(main) & length(unique(survey))>1) main <- sort(unique(as.character(survey)))
  if(is.null(main) & length(unique(survey))>1) main <- levels(survey)
  
  lines <- rivers$lines
  length <- length(lines)
  xmin <- min(lines[[1]][,1])
  xmax <- max(lines[[1]][,1])
  ymin <- min(lines[[1]][,2])
  ymax <- max(lines[[1]][,2])
  if(length>1) {
    for(j in 2:length) {
      if(min(lines[[j]][,1])<xmin) xmin <- min(lines[[j]][,1])
      if(max(lines[[j]][,1])>xmax) xmax <- max(lines[[j]][,1])
      if(min(lines[[j]][,2])<ymin) ymin <- min(lines[[j]][,2])
      if(max(lines[[j]][,2])>ymax) ymax <- max(lines[[j]][,2])
    }
  }
  
  nsize <- NA
  isurvey <- 1
  # for(surveyi in sort(unique(survey))) {
  for(surveyi in levels(survey)) {
    nsize[isurvey] <- length(seg[survey==surveyi])
    isurvey <- isurvey+1
  }
  
  iisurvey <- 1
  if(is.null(whichplots)) whichplots <- 1:length(unique(survey))
  # whichplotsurvey <- (sort(unique(survey)))[whichplots]
  whichplotsurvey <- (levels(survey))[whichplots]
  
  if(!scalebyN){ 
    isurvey <- 1
    # for(surveyi in sort(unique(survey))) {
    for(surveyi in levels(survey)) {
      for(segi in 1:length(rivers$lines)) {
        densities[[isurvey]][[segi]] <- densities[[isurvey]][[segi]]*max(nsize[whichplots])/nsize[isurvey]
      }
      isurvey <- isurvey+1
    }
  }
  
  for(surveyi in whichplotsurvey) {
    isurvey <- whichplots[iisurvey]
    if(showN) mainforplot <- paste0(main[isurvey],"  (n=",length(seg[survey==surveyi]),")")
    if(!showN) mainforplot <- main[isurvey]
    if(!add) plot(c(xmin,xmax),c(ymin,ymax),col="white",cex.axis=.6,asp=1,xlab=xlab,ylab=ylab,main=mainforplot,...=...)
    
    for(segi in 1:length(rivers$lines)) {
      quants <- (densities[[isurvey]][[segi]]/max(unlist(densities)))^pwr
      if(bycol) {
        if(ramp=="grey" | ramp=="gray") {
          cols <- grey((1-quants)*.8)
          denscol <- 1
          linecol <- grey(.8)
        }
        if(ramp=="red") {
          cols <- rgb(1,(1-quants)*.8,(1-quants)*.8)
          denscol <- 2
          linecol <- rgb(1,.8,.8)
        }
        if(ramp=="green") {
          cols <- rgb((1-quants)*.8,1,(1-quants)*.8)
          denscol <- 3
          linecol <- rgb(.8,1,.8)
        }
        if(ramp=="blue") {
          cols <- rgb((1-quants)*.8,(1-quants)*.8,1)
          denscol <- 4
          linecol <- rgb(.8,.8,1)
        }
        if(ramp=="heat") {
          cols <- heat.colors(1000)[ceiling(900*(1-quants))+1]
          denscol <- heat.colors(1000)[1]
          linecol <- heat.colors(1000)[901]
        }
        if(ramp=="stoplight") {
          cols <- rainbow(1000)[ceiling(300*(1-quants))+1] 
          denscol <- rainbow(1000)[1] 
          linecol <- rainbow(1000)[301]
        }
        if(ramp=="rainbow") {
          cols <- rainbow(1000)[ceiling(700*(1-quants))+1]
          denscol <- rainbow(1000)[1] 
          linecol <- rainbow(1000)[701]
        }
      }
      if(bylwd) {
        lwds <- maxlwd*quants
      }
      if(!bycol) cols <- rep(denscol,(length(endptverts[[segi]])-1))
      if(!bylwd) lwds <- rep(lwd,(length(endptverts[[segi]])-1))
      
      if(alpha<1) cols <- adjustcolor(cols,alpha.f=alpha)
      if(alpha<1) linecol <- adjustcolor(linecol,alpha.f=alpha)
      if(dark<1) cols <- adjustcolor(cols,red.f=dark,green.f=dark,blue.f=dark)
      if(dark<1) linecol <- adjustcolor(linecol,red.f=dark,green.f=dark,blue.f=dark)
      if(dark<1) denscol <- adjustcolor(denscol,red.f=dark,green.f=dark,blue.f=dark)
      
      for(vertsi in 1:(length(endptverts[[segi]])-1)) {
        if(dim(matrix(rivers$lines[[segi]][(endptverts[[segi]][vertsi]):(endptverts[[segi]][vertsi+1]),],ncol=2))[1] > 1) {
          lines(rivers$lines[[segi]][(endptverts[[segi]][vertsi]):(endptverts[[segi]][vertsi+1]),],lwd=lwd,col=linecol,lend=1)
          if(densities[[isurvey]][[segi]][vertsi] > 0) {
            lines(rivers$lines[[segi]][(endptverts[[segi]][vertsi]):(endptverts[[segi]][vertsi+1]),],lwd=lwds[vertsi],col=cols[vertsi],lend=1)
          }
        }
      }
    }
    if(points) riverpoints(seg=seg[survey==surveyi],vert=vert[survey==surveyi],rivers=rivers,pch=21,bg=0,col=denscol)
    if(scalebar) scalebar(rivers)
    iisurvey <- iisurvey+1
  }
}



#' Plot Points Used for Kernel Density
#' @description Plots the points used to calculate a kernel density object 
#'   in \link{makeriverdensity}.
#'   
#'   This function is intended as a visual check that a sufficient resolution was used.
#' @param riverdensity A river density object created by \link{makeriverdensity}.
#' @seealso \link{makeriverdensity}, \link{plot.riverdensity}
#' @author Matt Tyers
#' @examples
#' data(Gulk, fakefish)
#' 
#' Gulk_dens <- makeriverdensity(seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk)
#' 
#' plotriverdensitypoints(riverdensity=Gulk_dens)
#' @export
plotriverdensitypoints <- function(riverdensity) {
  lines <- riverdensity$rivers$lines
  verts <- riverdensity$densverts
  plot(x=riverdensity$rivers,segmentnum=FALSE,scale=FALSE)
  for(segi in 1:length(verts)) {
    riverpoints(seg=rep(segi,length(verts[[segi]])),vert=verts[[segi]],rivers=riverdensity$rivers)
  }
}
