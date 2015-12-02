#' Pythagorean Distance
#' @description Pythagorean distance between two points.  Called internally.
#' @param p1 X-Y coordinates of point 1
#' @param p2 X-Y coordinates of point 2
#' @return Distance (numeric)
#' @author Matt Tyers
#' @examples
#' point1 <- c(1,3)
#' point2 <- c(4,7)
#' 
#' pdist(point1,point2)
#' @export
pdist <- function(p1,p2) {
  dist <- sqrt((p1[1]-p2[1])^2 + (p1[2]-p2[2])^2)
  return(dist)
}

#' Create a River Network Object from a Shapefile
#' @description Uses \link[rgdal]{readOGR} in package 'rgdal' to read a river 
#'   shapefile, and establishes connectivity of segment endpoints based on 
#'   spatial proximity.
#' @param path File path, default is the current working directory
#' @param layer Name of the shapefile, without the file extension
#' @param tolerance Snapping tolerance of segment endpoints to determine 
#'   connectivity.  Default is 100, therefore care should be exercised when 
#'   working with larger units of distance, such as km.
#' @param reproject A valid Proj.4 projection string, if the shapefile is to be 
#'   re-projected.  Re-projection is done using \link[sp]{spTransform} in 
#'   package 'sp'.
#' @return Returns an object of class \code{"rivernetwork"} containing all
#'   spatial and topological information.  See \link{rivernetwork-class}.
#' @note Since distance can only be calculated using projected coordinates, 
#'   \code{line2network()} will generate an error if a non-projected input 
#'   shapefile is detected.  To resolve this, the shapefile can be re-projected 
#'   in a GIS environment, or using \code{reproject=}, shown in the second 
#'   example below.
#' @author Matt Tyers
#' @importFrom rgdal readOGR
#' @importFrom sp spTransform
#' @importFrom sp is.projected
#' @importFrom sp CRS
#' @examples 
#' filepath <- system.file("extdata", package="riverdist")
#' 
#' Gulk_UTM5 <- line2network(path=filepath, layer="Gulk_UTM5")
#' plot(Gulk_UTM5)
#' 
#' # # Re-projecting in Alaska Albers Equal Area projection:
#' 
#' AKalbers <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 
#'     +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
#'     
#' Gulk_AKalbers <- line2network(path=filepath, layer="Gulk_UTM5", reproject=AKalbers)
#' plot(Gulk_AKalbers)
#' 
#' @export
line2network <- function(path=".",layer,tolerance=100,reproject=NULL) {
  sp <- suppressWarnings(rgdal::readOGR(dsn=path,layer=layer,verbose=F)) # reading the shapefile as an sp object
  if(class(sp)!="SpatialLinesDataFrame") stop("Error - specified shapefile is not a linear feature.")
  proj4 <- strsplit(sp@proj4string@projargs,split=" ")
  projected <- sp::is.projected(sp)
  if(is.null(reproject) & !projected) stop("Error - Distances can only be computed from a projected coordinate system.  Use reproject= to specify a Proj.4 projection to use.")
  
  if(!is.null(reproject)) {
    sp <- spTransform(sp,CRS(reproject))
  }
  
  units <- "unknown"
  for(i in 1:length(proj4[[1]])) {
    proj4arg <- strsplit(proj4[[1]][i],split="=")
    if(proj4arg[[1]][1]=="+units") {
      units <- proj4arg[[1]][2]
      cat('\n',"Units:",proj4arg[[1]][2],'\n')
    }
  }
  
  if(length(sp@lines) > 1) {
    sp_line <- NA
    sp_seg <- NA
    #length <- length(sp@lines)
    lines <- list()
    j<-1
    for(i in 1:length(sp@lines)) {
      for(k in 1:length(sp@lines[i][[1]]@Lines)) {
        lines[[j]] <- sp@lines[i][[1]]@Lines[[k]]@coords
        sp_line[j] <- i
        sp_seg[j] <- k
        j<-j+1
      }
    }
  }
  if(length(sp@lines) == 1) {
    lines <- sp@lines[1][[1]]@Lines    # extracting just a list of lines and coordinates
    
    length <- length(lines) # number of line segments
    
    lines.new <- list()
    for(i in 1:length) {
      lines.new[[i]] <- lines[[i]]@coords
    }
    lines <- lines.new 
    sp_line <- rep(1,length)
    sp_seg <- 1:length
  }
  length <- length(lines)
  
  rivID <- 1:length
  lineID <- data.frame(rivID,sp_line,sp_seg)
  
  # defining a connectivity matrix...
  # connection type 1: beginning - beginning
  # connection type 2: beginning - end
  # connection type 3: end - beginning
  # connection type 4: end - end
  connections <- matrix(NA,nrow=length,ncol=length)
  for(i in 1:length) {
    for(j in 1:length) {
      i.max <- dim(lines[[i]])[1]
      j.max <- dim(lines[[j]])[1]
      if(pdist(lines[[i]][1,],lines[[j]][1,])<tolerance & i!=j) {
        connections[i,j] <- 1
      }
      if(pdist(lines[[i]][1,],lines[[j]][j.max,])<tolerance & i!=j) {
        connections[i,j] <- 2
      }
      if(pdist(lines[[i]][i.max,],lines[[j]][1,])<tolerance & i!=j) {
        connections[i,j] <- 3
      }
      if(pdist(lines[[i]][i.max,],lines[[j]][j.max,])<tolerance & i!=j) {
        connections[i,j] <- 4
      }
    }
  }
  
  # making a vector of total segment lengths
  lengths <- rep(NA,length)
  for(i in 1:length) {
    sum<-0
    linelength <- dim(lines[[i]])[1]
    for(j in 1:(linelength-1)) {
      sum <- sum+pdist(lines[[i]][j,],lines[[i]][(j+1),])
    }
    lengths[i]<-sum
  }
  
  names <- rep(NA,length)
  mouth.seg <- NA
  mouth.vert <- NA
  mouth <- list(mouth.seg,mouth.vert)
  names(mouth) <- c("mouth.seg","mouth.vert")
  sequenced <- FALSE
  braided <- NA
  
  out.names <- c("sp","lineID","lines","connections","lengths","names","mouth","sequenced","tolerance","units","braided")
  out <- list(sp,lineID,lines,connections,lengths,names,mouth,sequenced,tolerance,units,braided)
  names(out) <- out.names
  class(out) <- "rivernetwork"
  return(out)
}


#' Plotting a River Network
#' @description S3 plotting method for the \link{rivernetwork-class}.  Produces a map of all river segments of a river network object.
#' @aliases mapriver
#' @param x The river network object to plot
#' @param segmentnum Whether or not to plot segment numbers (defaults to TRUE)
#' @param lwd Line width
#' @param cex Global character expansion factor for plotting
#' @param scale Whether or not to give x- and y-axes the same scale
#' @param xlab Label for X-axis (defaults to "")
#' @param ylab Label for Y-axis (defaults to "")
#' @param ... Additional plotting arguments (see \link[graphics]{par})
#' @param color How to differentiate segments.  If \code{color==TRUE} (default),
#'   segments will be drawn in solid lines with differing colors.  If
#'   \code{color==FALSE}, segments will be drawn in black with differing line
#'   types.
#' @author Matt Tyers
#' @note This function is intended to provide basic visual checks for the user,
#'   not for any real mapping.
#' @examples 
#' data(Gulk)
#' plot(x=Gulk)
#' @method plot rivernetwork
#' @importFrom grDevices rgb
#' @importFrom graphics plot
#' @importFrom graphics par
#' @importFrom graphics text
#' @importFrom stats cor
#' @importFrom graphics axTicks
#' @export
plot.rivernetwork <- function(x,segmentnum=TRUE,lwd=1,cex=.6,scale=TRUE,color=TRUE,xlab="",ylab="",...) {
  rivers <- x
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
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
      # print(c(j,xmin,xmax,ymin,ymax))
    }
  }
  plot(c(xmin,xmax),c(ymin,ymax),col="white",cex.axis=.6,asp=1,xlab=xlab,ylab=ylab,...=...)
  #plot(c(xmin,xmax),c(ymin,ymax),col="white",main=main,xlab="",ylab="",cex.axis=.6,asp=1,xlim=xlim,ylim=ylim)
  midptx <- midpty <- NA
  for(j in 1:length) {
    if(!color) lines(lines[[j]],col=1,lty=j,lwd=lwd)
    if(color) lines(lines[[j]],col=rgb((sin(j)+1)/2.3,(cos(7*j)+1)/2.3,(sin(3*(length-j))+1)/2.3),lty=1,lwd=lwd)
    linelength <- dim(lines[[j]])[1]
    
    xplot <- lines[[j]][,1][lines[[j]][,1]>par("usr")[1] & lines[[j]][,1]<par("usr")[2] & lines[[j]][,2]>par("usr")[3] & lines[[j]][,2]<par("usr")[4]]
    yplot <- lines[[j]][,2][lines[[j]][,1]>par("usr")[1] & lines[[j]][,1]<par("usr")[2] & lines[[j]][,2]>par("usr")[3] & lines[[j]][,2]<par("usr")[4]]
    if(length(xplot)>0) midptx[j] <- mean(xplot)
    if(length(yplot)>0) midpty[j] <- mean(yplot)
    if(length(xplot)==0 | length(yplot)==0) midptx[j] <- midpty[j] <-NA
    if(segmentnum) text(midptx[j],midpty[j],j,cex=cex)
  }
  if(scale) {
    if(length>1) corthing <- cor(midptx,midpty,use="complete.obs")
    if(length==1) corthing <- (lines[[1]][1,1]-lines[[1]][dim(lines[[1]])[1],1])*(lines[[1]][1,2]-lines[[1]][dim(lines[[1]])[1],2])
    if(corthing<=0) {
      if(length(axTicks(1))>5) {
        scalex <- axTicks(1)[c(1,3)]
        labx <- axTicks(1)[2]
      }
      if(length(axTicks(1))<=5) {
        scalex <- axTicks(1)[c(1,2)]
        labx <- mean(axTicks(1)[c(1,2)])
      }
      scaley <- axTicks(2)[1]
    }
    if(corthing>0) {
      if(length(axTicks(1))>5) {
        scalex <- axTicks(1)[c(length(axTicks(1))-2,length(axTicks(1)))]
        labx <- axTicks(1)[length(axTicks(1))-1]
      }
      if(length(axTicks(1))<=5) {
        scalex <- axTicks(1)[c(length(axTicks(1))-1,length(axTicks(1)))]
        labx <- mean(axTicks(1)[c(length(axTicks(1))-1,length(axTicks(1)))])
      }
      scaley <- axTicks(2)[1]
    }
    lines(scalex,rep(scaley,2))
    if(rivers$units=="m") text(labx,scaley,labels=paste((scalex[2]-scalex[1])/1000,"km"),pos=3,cex=cex)
    if(rivers$units!="m") text(labx,scaley,labels=paste((scalex[2]-scalex[1]),rivers$units),pos=3,cex=cex)
  }
}



#' Highlight Segments
#' @description Plots a river network object and displays specified segments in bold, for easy identification.
#' @param seg A vector of segments to highlight
#' @param rivers The river network object to use
#' @param cex The character expansion factor to use for segment labels
#' @param lwd The line width to use for highlighted segments
#' @param add Whether to add the highlighted segments to an existing plot (\code{TRUE}) or call a new plot (\code{FALSE}).  Defaults to \code{FALSE}.
#' @param ... Additional plotting arguments (see \link[graphics]{par})
#' @author Matt Tyers
#' @examples
#' data(Kenai3)
#' plot(Kenai3)
#' highlightseg(seg=c(3,9,27),rivers=Kenai3)
#' @importFrom grDevices rgb
#' @importFrom graphics plot
#' @importFrom graphics par
#' @importFrom graphics text
#' @export
highlightseg <- function(seg,rivers,cex=0.8,lwd=3,add=FALSE,...) {
  length<-length(rivers$lines)
  if(!add) plot(rivers,color=FALSE,segmentnum=FALSE,...=...)
  lines<-rivers$lines
  midptx <- midpty <- NA
  for(j in seg) {
    lines(lines[[j]],col=rgb((sin(j)+1)/2.3,(cos(7*j)+1)/2.3,(sin(3*(length-j))+1)/2.3),lty=1,lwd=lwd)
    linelength <- dim(lines[[j]])[1]
    xplot <- lines[[j]][,1][lines[[j]][,1]>par("usr")[1] & lines[[j]][,1]<par("usr")[2] & lines[[j]][,2]>par("usr")[3] & lines[[j]][,2]<par("usr")[4]]
    yplot <- lines[[j]][,2][lines[[j]][,1]>par("usr")[1] & lines[[j]][,1]<par("usr")[2] & lines[[j]][,2]>par("usr")[3] & lines[[j]][,2]<par("usr")[4]]
    if(length(xplot)>0) midptx[j] <- mean(xplot)
    if(length(yplot)>0) midpty[j] <- mean(yplot)
    if(length(xplot)==0 | length(yplot)==0) midptx[j] <- midpty[j] <-NA
    text(midptx[j],midpty[j],j,cex=cex)
  }
}


#' Check Connectivity of a River Network Object
#' @description Produces a graphical check of the connectivity of a river 
#'   network object.  It produces a \link{plot} of the river network object,
#'   and overlays red dots at non-connected endpoints and green dots at 
#'   connected endpoints.
#' @param rivers The river network object to check
#' @param add Whether call a new plot (\code{FALSE}) or add dots to an existing
#'   plot (\code{TRUE}).  Defaults to \code{FALSE}.
#' @param ... Additional plotting arguments (see \link[graphics]{par})
#' @author Matt Tyers
#' @seealso \link{line2network}
#' @examples 
#' data(Gulk)
#' topologydots(rivers=Gulk)
#' @importFrom graphics plot
#' @importFrom graphics points
#' @export
topologydots <- function(rivers,add=FALSE,...) {
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  if(!add) plot(rivers,color=F,...=...)
  connections <- rivers$connections
  lines <- rivers$lines
  connections[is.na(connections)==T]<-0
  for(i in 1:dim(connections)[2]) {
    low <- 2
    high <- 2
    if(length(connections[i,][connections[i,]==1 | connections[i,]==2])>0) low <- 3
    if(length(connections[i,][connections[i,]==3 | connections[i,]==4])>0) high <- 3
    points(lines[[i]][1,1],lines[[i]][1,2],pch=16,col=low)
    points(lines[[i]][dim(lines[[i]])[1],1],
           lines[[i]][dim(lines[[i]])[1],2],pch=16,col=high)
  }
}


#' Check Which Segments are Connected to a Given Segment.
#' @description Returns which segments are connected to a specified segment within a river network.  It may be useful for error checking.
#' @param seg The segment to check
#' @param rivers The river network object it belongs to
#' @return A vector of segment numbers
#' @author Matt Tyers
#' @examples 
#' data(Gulk)
#' plot(Gulk)
#' whoconnected(seg=4, rivers=Gulk)
#' @export
whoconnected <- function(seg,rivers) {
  connections1 <- !is.na(rivers$connections)
  connected <- which(connections1[seg,])
  return(connected)
}


#' Convert XY Coordinates to River Coordinates
#' @description This function determines the closest vertex in the river network
#'   to each point of XY data and returns a list of river coordinates, defined 
#'   as segment numbers and vertex numbers.
#' @param x A vector of x-coordinates to transform
#' @param y A vector of y-coordinates to transform
#' @param rivers The river network object to use
#' @return A data frame of river coordinates, with segment numbers in \code{$seg} and vertex numbers in \code{$vert}.
#' @author Matt Tyers
#' @note Conversion to river coordinates is only valid if the input XY 
#'   coordinates and river network are in the same projected coordinate system. 
#'   Point data in geographic coordinates can be projected using 
#'   \link[rgdal]{project} in package 'rgdal', and an example is shown below.
#' @examples 
#' data(Gulk,fakefish)
#' head(fakefish)
#' 
#' fakefish.riv <- xy2segvert(x=fakefish$x, y=fakefish$y, rivers=Gulk)
#' head(fakefish.riv)
#' 
#' plot(x=Gulk, xlim=c(862000,882000), ylim=c(6978000,6993000))
#' points(fakefish$x, fakefish$y, pch=16, col=2)
#' riverpoints(seg=fakefish.riv$seg, vert=fakefish.riv$vert, rivers=Gulk, pch=15, col=4)
#' 
#' 
#' ## converting a matrix of points stored in long-lat to Alaska Albers Equal Area:
#' data(line98, Kenai1)
#' head(line98)  # note that coordinates are stored in long-lat, NOT lat-long
#' 
#' library(rgdal)
#' line98albers <- project(line98,proj="+proj=aea +lat_1=55 +lat_2=65 
#'     +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs 
#'     +ellps=GRS80 +towgs84=0,0,0")
#' head(line98albers)
#' 
#' zoomtoseg(seg=c(162,19), rivers=Kenai1)
#' points(line98albers)
#' @export
xy2segvert <- function(x,y,rivers) {
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  seg <- rep(NA,length(x))
  vert <- rep(NA,length(x))
  lines <- rivers$lines
  
  for(i in 1:length(x)) {
    min.dist1 <- 1000000
    for(river.i in 1:length(lines)) {
      for(vert.i in 1:dim(lines[[river.i]])[1]) {
        dist.trial <- pdist(c(x[i],y[i]),lines[[river.i]][vert.i,])
        if(dist.trial<min.dist1) {
          min.dist1 <- dist.trial
          seg[i] <- river.i
          vert[i] <- vert.i
        }
      }
    }
  }
  out <- data.frame(seg,vert)
  return(out)
}

#' Draw Points from River Coordinates
#' @description Adds points to an active plot.  Works like \link[graphics]{points} but with river coordinates (segments and vertices) rather than xy coordinates.
#' @param seg A vector of segments
#' @param vert A vector of vertices
#' @param rivers The river network object to use
#' @param pch Point character
#' @param col Point color
#' @author Matt Tyers
#' @examples
#' data(fakefish,Gulk)
#' 
#' plot(x=Gulk, xlim=c(862000,882000), ylim=c(6978000,6993000))
#' points(x=fakefish$x, y=fakefish$y, pch=16, col=2)
#' riverpoints(seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk, pch=15, col=4)
#' @importFrom graphics points
#' @export
riverpoints <- function(seg,vert,rivers,pch=1,col=1) {
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  lines <- rivers$lines
  if(max(seg,na.rm=T)>length(lines) | min(seg,na.rm=T)<1) stop("Invalid segment numbers specified.")
  
  for(i in 1:length(seg)) {
    if(vert[i]>dim(lines[[seg[i]]])[1] | vert[i]<1) stop("Invalid vertex numbers specified.")
    points(lines[[seg[i]]][vert[i],1],
           lines[[seg[i]]][vert[i],2],pch=pch,col=col)
  }
}

#' Trim a River Network Object to Specified Segments
#' @description Removes line segments from a river network object.  User can specify which segments to remove (\code{trim}) or which segments to keep (\code{trimto}).
#' @param trim Vector of line segments to remove
#' @param trimto Vector of line segments to keep
#' @param rivers The river network object
#' @return A new river network object
#' @seealso \link{line2network}
#' @note Specifying segments in both trim and trimto arguments will result in an error.
#' @author Matt Tyers
#' @examples
#' data(Kenai1)
#' plot(x=Kenai1)
#' 
#' Kenai1.trim <- trimriver(trim=c(46,32,115,174,169,114,124,142,80), rivers=Kenai1)
#' plot(x=Kenai1.trim)
#'  
#' Kenai1.trim.2 <- trimriver(trimto=c(20,57,118,183,45,162,39,98,19), rivers=Kenai1)
#' plot(x=Kenai1.trim.2)
#' @export
trimriver <- function(trim=NULL,trimto=NULL,rivers) {
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  segs <- 1:length(rivers$lines)
  if(!is.null(trim) & !is.null(trimto)) {
    stop("Error - cannot use both trim and trimto arguments")
  }
  if(!is.null(trim)) {
    if(max(trim,na.rm=T)>length(rivers$lines) | min(trim,na.rm=T)<1) stop("Invalid segment numbers specified.")
    segs <- segs[-trim]
  }
  if(!is.null(trimto)) {
    if(max(trimto,na.rm=T)>length(rivers$lines) | min(trimto,na.rm=T)<1) stop("Invalid segment numbers specified.")
    segs <- segs[trimto]
  }
  trimmed.rivers <- rivers
  trimmed.rivers$lines <- trimmed.rivers$lines[segs]
  trimmed.rivers$connections <- trimmed.rivers$connections[segs,segs]
  trimmed.rivers$names <- rivers$names[segs]
  trimmed.rivers$lengths <- rivers$lengths[segs]
  if(length(segs)==0) stop("Error - resulting river network has no remaining line segments")
  if(!is.na(trimmed.rivers$mouth$mouth.seg) & !is.na(trimmed.rivers$mouth$mouth.vert)) {
    if(any(segs==rivers$mouth$mouth.seg)) {
      trimmed.rivers$mouth$mouth.seg <- which(segs==trimmed.rivers$mouth$mouth.seg)
    }
    if(!any(segs==rivers$mouth$mouth.seg)) {
      trimmed.rivers$mouth$mouth.seg <- NA
      trimmed.rivers$mouth$mouth.vert <- NA
    }
  }
  
  # updating sp object
  id <- rivers$lineID
  sp_lines1 <- rivers$sp@lines[unique(id[segs,2])]   # this is a little kludgy - want to combine Lines in one line when sp_line is the same
  j<-1
  for(i in unique(id[segs,2])) {
    sp_lines1[[j]]@Lines <- sp_lines1[[j]]@Lines[id[segs,3][id[segs,2]==i]]
    j<-j+1
  }
  rivID <- NA
  sp_line <- NA
  sp_seg <- NA
  k<-1
  for(i in 1:length(sp_lines1)) {
    for(j in 1:length(sp_lines1[[i]]@Lines)) {
      sp_line[k] <- i
      sp_seg[k] <- j
      k<-k+1
    }
  }
  rivID <- 1:(k-1)
  lineID <- data.frame(rivID,sp_line,sp_seg)
  trimmed.rivers$lineID <- lineID
  trimmed.rivers$sp@lines <- sp_lines1
  if(dim(rivers$sp@data)[1]==max(rivers$lineID[,2]) & dim(rivers$sp@data)[1]>1) {
    trimmed.rivers$sp@data <- rivers$sp@data[unique(rivers$lineID[segs,2]),]
  }
  return(trimmed.rivers)
}


#' Trim a River Network to a Set of X-Y Coordinates
#' @description Removes line segments from a river network object that are not
#'   adjacent to a set of point data, given in X-Y coordinates.
#' @param x Vector of x-coordinates of point data to buffer around
#' @param y Vector of y-coordinates of point data to buffer around
#' @param rivers The river network object to use
#' @param method Three methods are available.  If \code{"snap"} is specified
#'   (the default), only the closest segment to each point is retained.   If
#'   \code{"snaproute"} is specified, segments are also retained that will
#'   maintain total connectivity in the resulting river network.  If
#'   \code{"buffer"} is specified, all segments with endpoints or midpoints
#'   within \code{dist} units of the input locations are retained.
#' @param dist Distance to use for buffering, if \code{method="buffer"}.  If
#'   this is not specified, the maximum spread in the x- and y- direction will
#'   be used.
#' @return A new river network object (see \link{rivernetwork})
#' @seealso \link{line2network}
#' @note If \code{method=="buffer"}, only distances to segment endpoints 
#'   and midpoints are checked, and still only whole segments are removed.
#' @author Matt Tyers
#' @examples
#' data(Koyukuk2)
#' x <- c(139241.0, 139416.1, 124600.1, 122226.8)
#' y <- c(1917577, 1913864, 1898723, 1898792)
#' 
#' plot(x=Koyukuk2)
#' points(x, y, pch=15, col=4)
#' legend(par("usr")[1], par("usr")[4], legend="points to buffer around", pch=15, col=4, cex=.6)
#' 
#' Koyukuk2.buf1 <- trimtopoints(x, y, rivers=Koyukuk2, method="snap")
#' plot(x=Koyukuk2.buf1)
#' points(x, y, pch=15, col=4)
#' 
#' Koyukuk2.buf2 <- trimtopoints(x, y, rivers=Koyukuk2, method="snaproute")
#' plot(x=Koyukuk2.buf2)
#' points(x, y, pch=15, col=4)
#' 
#' Koyukuk2.buf3 <- trimtopoints(x, y, rivers=Koyukuk2, method="buffer", dist=1000)
#' plot(x=Koyukuk2.buf3)
#' points(x, y, pch=15, col=4)
#' @export
trimtopoints <- function(x,y,rivers,method="snap",dist=NULL) {
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  if(is.null(dist)) dist <- max(c(max(x)-min(x)),c(max(y)-min(y)))
  
  if(method=="buffer") {
    keep <- 0
    k <- 1
    for(i in 1:length(x)) {
      for(j in 1:length(rivers$lines)) {
        line.len <- dim(rivers$lines[[j]])[1]
        if((pdist(c(x[i],y[i]),rivers$lines[[j]][1,]) <= dist) | 
           (pdist(c(x[i],y[i]),rivers$lines[[j]][line.len,]) <= dist) |
           (pdist(c(x[i],y[i]),rivers$lines[[j]][floor(line.len/2),]) <= dist)) {
          if(length(keep[keep==j])==0) {
            keep[k] <- j
            k <- k+1
          }
        }
      }
    }
  }
  
  if(method=="bufferroute") {
    keep <- 0
    k <- 1
    for(i in 1:length(x)) {
      for(j in 1:length(rivers$lines)) {
        line.len <- dim(rivers$lines[[j]])[1]
        if((pdist(c(x[i],y[i]),rivers$lines[[j]][1,]) <= dist) | 
           (pdist(c(x[i],y[i]),rivers$lines[[j]][line.len,]) <= dist) |
           (pdist(c(x[i],y[i]),rivers$lines[[j]][floor(line.len/2),]) <= dist)) {
          if(length(keep[keep==j])==0) {
            keep[k] <- j
            k <- k+1
          }
        }
      }
    }
    
    keep.bool <- rep(FALSE,length(rivers$lines))
    keep.bool[keep] <- TRUE
    
    for(i in keep) {
      for(j in keep) {
        if(i!=j) {
          route1 <- detectroute(start=i,end=j,rivers=rivers)
          keep.bool[route1] <- TRUE
        }
      }
    }
    keep <- (1:length(rivers$lines))[keep.bool]
  }
  
  if(method=="snap") {
    segvert <- xy2segvert(x,y,rivers=rivers)
    keep.bool <- rep(FALSE,length(rivers$lines))
    keep.bool[unique(segvert$seg)] <- TRUE
    keep <- (1:length(rivers$lines))[keep.bool]
  }
  
  if(method=="snaproute") {
    segvert <- xy2segvert(x=x,y=y,rivers)
    keep.bool <- rep(FALSE,length(rivers$lines))
    keep.bool[unique(segvert$seg)] <- TRUE
    
    for(i in segvert$seg) {
      for(j in segvert$seg) {
        if(i!=j) {
          route1 <- detectroute(start=i,end=j,rivers=rivers)
          keep.bool[route1] <- TRUE
        }
      }
    }
    keep <- (1:length(rivers$lines))[keep.bool]
  }
  
  rivers1 <- rivers
  rivers1$lines <- rivers1$lines[keep]
  rivers1$connections <- rivers1$connections[keep,keep]
  rivers1$names <- rivers$names[keep]
  rivers1$lengths <- rivers$lengths[keep]
  if(keep[1]==0) stop("Error - resulting river network has no remaining line segments")
  
  if(!is.na(rivers1$mouth$mouth.seg) & !is.na(rivers1$mouth$mouth.vert)) {
    if(any(keep==rivers$mouth$mouth.seg)) {
      rivers1$mouth$mouth.seg <- which(keep==rivers1$mouth$mouth.seg)
    }
    if(!any(keep==rivers$mouth$mouth.seg)) {
      rivers1$mouth$mouth.seg <- NA
      rivers1$mouth$mouth.vert <- NA
    }
  }
  
  # updating sp object
  id <- rivers$lineID
  sp_lines1 <- rivers$sp@lines[unique(id[keep,2])]   # this is a little kludgy - want to combine Lines in one line when sp_line is the same
  j<-1
  for(i in unique(id[keep,2])) {
    sp_lines1[[j]]@Lines <- sp_lines1[[j]]@Lines[id[keep,3][id[keep,2]==i]]
    j<-j+1
  }
  rivID <- NA
  sp_line <- NA
  sp_seg <- NA
  k<-1
  for(i in 1:length(sp_lines1)) {
    for(j in 1:length(sp_lines1[[i]]@Lines)) {
      sp_line[k] <- i
      sp_seg[k] <- j
      k<-k+1
    }
  }
  rivID <- 1:(k-1)
  lineID <- data.frame(rivID,sp_line,sp_seg)
  rivers1$lineID <- lineID
  rivers1$sp@lines <- sp_lines1
  if(dim(rivers$sp@data)[1]==max(rivers$lineID[,2]) & dim(rivers$sp@data)[1]>1) {
    rivers1$sp@data <- rivers$sp@data[unique(rivers$lineID[keep,2]),]
  }
  return(rivers1)
}


#' Remove Unconnected Segments
#' @description Detects and removes segments that are not connected to the river mouth.
#' @param rivers The river network object to use.
#' @author Matt Tyers
#' @note This function is called within \link{cleanup}, which is recommended in most cases.
#' @examples 
#' data(Kenai1)
#' Kenai1.1 <- dissolve(Kenai1)
#' Kenai1.1 <- setmouth(seg=63,vert=40,rivers=Kenai1.1)
#' plot(Kenai1.1)
#' 
#' Kenai1.2 <- removeunconnected(Kenai1.1)
#' plot(Kenai1.2)
#' @export
removeunconnected <- function(rivers) {
  if(is.na(rivers$mouth$mouth.seg)) stop("River mouth must be specified.")
  origin <- rivers$mouth$mouth.seg
  
  takeout <- NULL
  k <- 1
  
  dists <- rep(NA,length(rivers$lines))
  for(i in 1:length(rivers$lines)) {
    dists[i] <- pdist(rivers$lines[[rivers$mouth$mouth.seg]][rivers$mouth$mouth.vert,],
                      c(mean(rivers$lines[[i]][,1]),mean(rivers$lines[[i]][,2])))
  }
  order <- order(dists)
  rivers1 <- rivers
  rivers1$connections <- rivers$connections[order,order]
  rivers1$lines <- rivers$lines[order]
  origin1 <-which(order==origin)
  
  #plot(rivers1)
  
  for(i in 1:length(rivers1$lines)) {
    # catch <- try(detectroute(end=origin1,start=i,rivers=rivers1))
    # if(class(catch)=="try-error") {
    #   takeout[k] <- i
    #   k <- k+1
    # }
    if(is.na(detectroute(end=origin1,start=i,rivers=rivers1,stopiferror=FALSE,algorithm="sequential")[1])) {
      takeout[k] <- i
      k <- k+1
    }
  }
  
  if(!is.null(takeout)) takeout <- order[takeout]
  
  rivers2 <- trimriver(trim=takeout,rivers=rivers)
  return(rivers2)
}