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

#' Total Pythagorean Distance
#' @description Total Pythagorean distance of a sequence of points.  Called internally.
#' @param xy A matrix of X-Y coordinates of the sequence of points.
#' @return Distance (numeric)
#' @author Matt Tyers
#' @examples
#' points <- matrix(c(1:10), nrow=5, ncol=2, byrow=FALSE)
#' 
#' pdisttot(xy=points)
#' @export
pdisttot <- function(xy) {
  n <- dim(xy)[1]
  if(n==1) dist <- 0
  if(n==2) dist <- pdist(xy[1,],xy[2,])
  if(n>2) dist <- sqrt(((xy[1:(n-1),1] - xy[2:n,1])^2) + ((xy[1:(n-1),2] - xy[2:n,2])^2))
  return(sum(dist))
}


#' Add Cumulative Distance to a River Network
#' @description Adds a vector of cumulative distances to a river network.  Called internally.
#' @param rivers The river network object to use.
#' @return Returns an object of class \code{"rivernetwork"} containing all
#'   spatial and topological information.  See \link{rivernetwork-class}.
#' @author Matt Tyers
#' @examples
#' Gulk1 <- addcumuldist(rivers=Gulk)
#' @export
addcumuldist <- function(rivers) {
  cumuldist <- list()
  for(i in 1:length(rivers$lines)) {
    xy <- rivers$lines[[i]]
    n <- dim(xy)[1]
    cumuldist[[i]] <- c(0,cumsum(sqrt(((xy[1:(n-1),1] - xy[2:n,1])^2) + ((xy[1:(n-1),2] - xy[2:n,2])^2))))
  }
  rivers$cumuldist <- cumuldist
  return(rivers)
}


#' Calculate the Connectivity Matrix for a River Network
#' @description Calculates the connectivity matrix for a river network, during import and editing.  Called internally.
#' @param lines A list of coordinate matrices, each corresponding to a line segment.
#' @param tolerance The spatial tolerance for establishing connectivity.
#' @return A matrix with topological information.  See the \code{$connections} element of the \link{rivernetwork-class}.
#' @author Matt Tyers
#' @examples
#' Gulk_connections <- calculateconnections(lines=Gulk$lines, tolerance=Gulk$tolerance)
#' @export
calculateconnections <- function(lines,tolerance) {
  length <- length(lines)
  # defining a connectivity matrix...
  # connection type 1: beginning - beginning
  # connection type 2: beginning - end
  # connection type 3: end - beginning
  # connection type 4: end - end
  # connection type 5: beginning - beginning and end - end
  # connection type 6: beginning - end and end - beginning
  connections <- matrix(NA,nrow=length,ncol=length)
  begmat <- matrix(unlist(sapply(lines,function(xy) xy[1,],simplify=F)),ncol=2,byrow=T)
  endmat <- matrix(unlist(sapply(lines,function(xy) xy[nrow(xy),],simplify=F)),ncol=2,byrow=T)
  pdist2 <- function(p1,p2mat) {
    dist <- sqrt((p1[1]-p2mat[,1])^2 + (p1[2]-p2mat[,2])^2)
    return(dist)
  }  
  
  for(i in 1:length) {
    mat1 <- pdist2(begmat[i,],begmat)
    mat2 <- pdist2(begmat[i,],endmat)
    mat3 <- pdist2(endmat[i,],begmat)
    mat4 <- pdist2(endmat[i,],endmat)
    connections[i,mat1<tolerance] <- 1
    connections[i,mat2<tolerance] <- 2
    connections[i,mat3<tolerance] <- 3
    connections[i,mat4<tolerance] <- 4
    connections[i,(mat1<tolerance & mat4<tolerance)] <- 5
    connections[i,(mat2<tolerance & mat3<tolerance)] <- 6
  }
  diag(connections) <- NA
  return(connections)
}


#' Create a River Network Object from a Shapefile
#' @description Uses \link[sf]{read_sf} in package 'sf' to read a river 
#'   shapefile, and establishes connectivity of segment endpoints based on 
#'   spatial proximity.
#' @param sf Optional input as an \link[sf]{sf} object, if shapefile has 
#'   already been read into the R environment.
#' @param path File path, default is the current working directory.
#' @param layer Name of the shapefile, without the .shp extension.
#' @param tolerance Snapping tolerance of segment endpoints to determine 
#'   connectivity.  Default is 100, therefore care should be exercised when 
#'   working with larger units of distance, such as km.
#' @param reproject A valid projection, if the shapefile is to be 
#'   re-projected.  Re-projection is done using \link[sf]{st_transform} in 
#'   package 'sf'.
#' @param autofix Whether to automatically apply two corrections: removal of 
#'   duplicate segments, and segments with lengths shorter than the connectivity 
#'   tolerance.  Defaults to `TRUE`.
#' @return Returns an object of class \code{"rivernetwork"} containing all
#'   spatial and topological information.  See \link{rivernetwork-class}.
#' @note Since distance can only be calculated using projected coordinates, 
#'   \code{line2network()} will generate an error if a non-projected input 
#'   shapefile is detected.  To resolve this, the shapefile can be re-projected 
#'   in a GIS environment, or using \code{reproject=}, shown in the second 
#'   example below.
#' @author Matt Tyers, Jemma Stachelek
#' @importFrom sf read_sf
#' @importFrom sf as_Spatial
#' @importFrom sf st_zm
#' @importFrom sf st_is_longlat
#' @importFrom sf st_transform
#' @examples 
#' filepath <- system.file("extdata", package="riverdist")
#' 
#' Gulk_UTM5 <- line2network(path=filepath, layer="Gulk_UTM5")
#' plot(Gulk_UTM5)
#' 
#' ## Reading directly from an sf object
#' 
#' sf <- sf::read_sf(dsn = filepath, layer = "Gulk_UTM5")
#' Gulk_UTM5 <- line2network(sf=sf)
#' plot(Gulk_UTM5)
#' 
#' ## Re-projecting in Alaska Albers Equal Area projection:
#' 
#' AKalbers <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 
#'     +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80"
#'     
#' Gulk_AKalbers <- line2network(path=filepath, layer="Gulk_UTM5", reproject=AKalbers)
#' plot(Gulk_AKalbers)
#' 
#' @export
line2network <- function(sf = NULL, path=".", layer = NA, tolerance=100, 
                         reproject=NULL, autofix=TRUE) {
  
  if(is.null(sf)) {
    ## read to sf here
    sf <- sf::read_sf(dsn = path, layer = layer)
  }
  
  if(!inherits(sf, c("sf"))) stop("Invalid input to sf= argument")
  
  if(!all(sf::st_geometry_type(sf) %in% c("LINESTRING", "MULTILINESTRING"))) {
    stop("Invalid input.  Either specified shapefile is not a linear feature, 
         or not all geometry types are LINESTRING or MULTILINESTRING.")
  } 
  
  projargs <- sf::st_crs(sf)$proj4string
  
  if(is.na(projargs)){ 
    stop("Shapefile projection information is missing.")
  }
  
  projected <- !sf::st_is_longlat(sf)
  
  if(is.null(reproject) & !projected) stop("Distances can only be computed from a projected coordinate system.  Use reproject= to specify a projection to use.")
  
  if(!is.null(reproject)) {
    sf <- sf::st_transform(sf, crs=reproject)
  }
  
  units <- sf::st_crs(sf)$units_gdal
  cat('\n',"Units:",units,'\n')
  
  lines <- list()
  j<-1
  for(i in seq_along(sf$geometry)) {
    if(sf::st_geometry_type(sf$geometry[[i]]) == "LINESTRING") {
      lines[[j]] <- sf$geometry[[i]][,1:2]
      j <- j+1
    }
    if(sf::st_geometry_type(sf$geometry[[i]]) == "MULTILINESTRING") {
      for(k in seq_along(sf$geometry[[i]])) {
        lines[[j]] <- sf$geometry[[i]][[k]][,1:2]
        j<-j+1
      }
    }
  }
  
  length <- length(lines)
  
  connections <- calculateconnections(lines=lines, tolerance=tolerance)
  
  if(any(connections %in% 5:6)) braided <- TRUE
  
  # making a vector of total segment lengths
  lengths <- rep(NA,length)
  for(i in 1:length) {
    lengths[i] <- pdisttot(lines[[i]])
  }
  
  names <- rep(NA,length)
  mouth.seg <- NA
  mouth.vert <- NA
  mouth <- list(mouth.seg,mouth.vert)
  names(mouth) <- c("mouth.seg","mouth.vert")
  sequenced <- FALSE
  braided <- NA
  
  cumuldist <- list()
  for(i in 1:length) {
    xy <- lines[[i]]
    n <- dim(xy)[1]
    cumuldist[[i]] <- c(0,cumsum(sqrt(((xy[1:(n-1),1] - xy[2:n,1])^2) + ((xy[1:(n-1),2] - xy[2:n,2])^2))))
  }
  
  out <- list(sf=sf, lines=lines, connections=connections, lengths=lengths, names=names, mouth=mouth,
              sequenced=sequenced, tolerance=tolerance, units=units, braided=braided, cumuldist=cumuldist)
  class(out) <- "rivernetwork"
  
  if(autofix) {
    length1 <- length(out$lengths)
    suppressMessages(out <- removeduplicates(out))
    length2 <- length(out$lengths)
    if(length2<length1) cat('\n',"Removed",length1-length2,"duplicate segments.",'\n')
    suppressMessages(out <- removemicrosegs(out))
    length3 <- length(out$lengths)
    if(length3<length2) cat('\n',"Removed",length2-length3,"segments with lengths shorter than the connectivity tolerance.",'\n')
  }
  
  return(out)
}


#' Convert a Point Shapefile to River Locations
#' @description This function reads a point shapefile and determines the closest
#'   vertex in the river network to each point of XY data, returning a data
#'   frame with river locations, defined as segment numbers and vertex
#'   numbers, along with the data table read from the input shapefile.
#' @param path File path, default is the current working directory.
#' @param layer Name of the shapefile, without the .shp extension.
#' @param rivers The river network object to use.
#' @return A data frame of river locations, with segment numbers in
#'   \code{$seg}, vertex numbers in \code{$vert}, snapping distances in 
#'   \code{$snapdist}, snapped x- and y-coordinates in \code{$snap_x} and \code{$snap_y},
#'   and the remaining columns
#'   corresponding to the data table in the input point shapefile.
#' @author Matt Tyers
#' @note If the input shapefile is detected to be in a different projection than
#'   the river network, the input shapefile will be re-projected before
#'   conversion to river locations.
#' @seealso \link{xy2segvert}
#' @importFrom sf read_sf
#' @importFrom sf st_coordinates
#' @importFrom sf st_drop_geometry
#' @importFrom sf st_transform
#' @examples 
#' filepath <- system.file("extdata", package="riverdist")
#' 
#' fakefish_UTM5 <- pointshp2segvert(path=filepath, layer="fakefish_UTM5", rivers=Gulk)
#' head(fakefish_UTM5)
#' 
#' plot(x=Gulk)
#' points(fakefish_UTM5$x, fakefish_UTM5$y)
#' riverpoints(seg=fakefish_UTM5$seg, vert=fakefish_UTM5$vert, rivers=Gulk, pch=16, col=2)
#' 
#' @export
pointshp2segvert <- function(path=".",layer,rivers) {
  shp <- sf::read_sf(dsn=path,layer=layer)
  if(sf::st_crs(shp)$proj4string != sf::st_crs(rivers$sf)$proj4string) {
    message('\n',"Point projection detected as different from river network.  Re-projecting points before snapping to river network...")
    shp <- sf::st_transform(shp, crs=sf::st_crs(rivers$sf)) ##
  }
  coords <- sf::st_coordinates(shp)
  segvert <- xy2segvert(x=coords[,1],y=coords[,2],rivers=rivers)
  outdf <- cbind(segvert, sf::st_drop_geometry(shp))
  return(outdf)
}


#' Plotting a River Network
#' @description S3 plotting method for the \link{rivernetwork-class}.  
#'   Produces a map of all river segments of a river network object.
#' @aliases mapriver
#' @param x The river network object to plot
#' @param segmentnum Whether or not to plot segment numbers (defaults to TRUE)
#' @param offset Whether to offset segment numbers from lines (defaults to TRUE)
#' @param lwd Line width
#' @param cex Global character expansion factor for plotting
#' @param scale Whether or not to give x- and y-axes the same scale
#' @param color How to differentiate segments.  If \code{color==TRUE} (default),
#'   segments will be drawn in solid lines with differing colors.  If
#'   \code{color==FALSE}, segments will be drawn in the same color with differing line
#'   types.
#' @param empty Creates an empty plot if set to \code{TRUE}.  Suppresses 
#'   differentiation by line type if \code{color==FALSE}, and suppresses segment 
#'   number labels.  Defaults to \code{FALSE}.
#' @param linecol Line color to use if \code{empty} is \code{TRUE} or 
#'   \code{color} is \code{FALSE}.  Defaults to black.
#' @param xlab Label for X-axis (defaults to "")
#' @param ylab Label for Y-axis (defaults to "")
#' @param ... Additional plotting arguments (see \link[graphics]{par})
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
plot.rivernetwork <- function(x,segmentnum=TRUE,offset=TRUE,lwd=1,cex=.6,scale=TRUE,color=TRUE,empty=FALSE,linecol=1,xlab="",ylab="",...) {
  if(!inherits(x, "rivernetwork")) stop("Argument 'x' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  lines <- x$lines
  length <- length(lines)
  allx <- unlist(lapply(lines,FUN='[',TRUE,1))
  ally <- unlist(lapply(lines,FUN='[',TRUE,2))
  xmax <- max(allx, na.rm=T)
  xmin <- min(allx, na.rm=T)
  ymax <- max(ally, na.rm=T)
  ymin <- min(ally, na.rm=T)
  plot(c(xmin,xmax),c(ymin,ymax),col="white",cex.axis=.6,asp=1,xlab=xlab,ylab=ylab,...=...)
  parusr <- par("usr")
  midptx <- midpty <- NA
  if(color&!empty) {
    linecolors <- rgb((sin(1:length)+1)/2.3,(cos(7*1:length)+1)/2.3,(sin(3*(length-(1:length)))+1)/2.3)
    ltys <- rep(1,length)
  }
  if(!color&!empty) {
    linecolors <- rep(linecol,length)
    ltys <- 1:length
  }
  if(empty) {
    linecolors <- rep(linecol,length)
    ltys <- rep(1,length)
    segmentnum <- F
  }
  
  xtext <- ytext <- direction <- rep(NA,length)
  for(j in 1:length) {
    lines(lines[[j]],col=linecolors[j],lty=ltys[j],lwd=lwd)
    
    xplot <- lines[[j]][,1][lines[[j]][,1]>parusr[1] & lines[[j]][,1]<parusr[2] & lines[[j]][,2]>parusr[3] & lines[[j]][,2]<parusr[4]]
    yplot <- lines[[j]][,2][lines[[j]][,1]>parusr[1] & lines[[j]][,1]<parusr[2] & lines[[j]][,2]>parusr[3] & lines[[j]][,2]<parusr[4]]
    if(length(xplot)>0) midptx[j] <- mean(xplot)
    if(length(yplot)>0) midpty[j] <- mean(yplot)
    if(length(xplot)==0 | length(yplot)==0) midptx[j] <- midpty[j] <-NA
    middle <- floor(length(xplot)/2)
    directionj <- 3+(abs(xplot[floor(length(xplot)*.25)]-xplot[floor(length(xplot)*.75)]) < abs(yplot[floor(length(yplot)*.25)]-yplot[floor(length(yplot)*.75)]))
    direction[j] <- ifelse(length(directionj)==0,1,directionj)
    xtext[j] <- ifelse(length(xplot)>0,xplot[middle],NA)
    ytext[j] <- ifelse(length(yplot)>0,yplot[middle],NA)
  }
  if(segmentnum) {
    if(offset) text(xtext,ytext,labels=1:length,pos=direction,cex=cex,col=linecolors)
    else text(xtext,ytext,labels=1:length,cex=cex,col=linecolors)
  }
  if(scale) {
    if(length>1) corthing <- cor(midptx,midpty,use="complete.obs")
    if(length==1) corthing <- (lines[[1]][1,1]-lines[[1]][dim(lines[[1]])[1],1])*(lines[[1]][1,2]-lines[[1]][dim(lines[[1]])[1],2])
    axticks1 <- axTicks(1)
    if(is.na(corthing)) legendleft <- F
    else legendleft <- corthing<=0
    if(legendleft) {   
      if(length(axticks1)>5) {
        scalex <- axticks1[c(1,3)]
        labx <- axticks1[2]
      }
      if(length(axticks1)<=5) {
        scalex <- axticks1[c(1,2)]
        labx <- mean(axticks1[c(1,2)])
      }
      scaley <- axTicks(2)[1]
    }
    if(!legendleft) {  
      if(length(axticks1)>5) {
        scalex <- axticks1[c(length(axticks1)-2,length(axticks1))]
        labx <- axticks1[length(axticks1)-1]
      }
      if(length(axticks1)<=5) {
        scalex <- axticks1[c(length(axticks1)-1,length(axticks1))]
        labx <- mean(axticks1[c(length(axticks1)-1,length(axticks1))])
      }
      scaley <- axTicks(2)[1]
    }
    lines(scalex,rep(scaley,2))
    if(x$units %in% c("m", "metre")) text(labx,scaley,labels=paste((scalex[2]-scalex[1])/1000,"km"),pos=3,cex=cex)
    if(!x$units %in% c("m", "metre")) text(labx,scaley,labels=paste((scalex[2]-scalex[1]),x$units),pos=3,cex=cex)
  }
}


scalebar <- function(rivers,cex=.6) {
  lines <- rivers$lines
  length <- length(lines)
  
  midptx <- midpty <- rep(NA,length)
  parusr <- par("usr")
  for(j in 1:length) {
    xplot <- lines[[j]][,1][lines[[j]][,1]>parusr[1] & lines[[j]][,1]<parusr[2] & lines[[j]][,2]>parusr[3] & lines[[j]][,2]<parusr[4]]
    yplot <- lines[[j]][,2][lines[[j]][,1]>parusr[1] & lines[[j]][,1]<parusr[2] & lines[[j]][,2]>parusr[3] & lines[[j]][,2]<parusr[4]]
    if(length(xplot)>0) midptx[j] <- mean(xplot)
    if(length(yplot)>0) midpty[j] <- mean(yplot)
  }
  
  if(length>1) corthing <- cor(midptx,midpty,use="complete.obs")
  if(length==1) corthing <- (lines[[1]][1,1]-lines[[1]][dim(lines[[1]])[1],1])*(lines[[1]][1,2]-lines[[1]][dim(lines[[1]])[1],2])
  axticks1 <- axTicks(1)
  if(corthing<=0) {
    if(length(axticks1)>5) {
      scalex <- axticks1[c(1,3)]
      labx <- axticks1[2]
    }
    if(length(axticks1)<=5) {
      scalex <- axticks1[c(1,2)]
      labx <- mean(axticks1[c(1,2)])
    }
    scaley <- axTicks(2)[1]
  }
  if(corthing>0) {
    if(length(axticks1)>5) {
      scalex <- axticks1[c(length(axticks1)-2,length(axticks1))]
      labx <- axticks1[length(axticks1)-1]
    }
    if(length(axticks1)<=5) {
      scalex <- axticks1[c(length(axticks1)-1,length(axticks1))]
      labx <- mean(axticks1[c(length(axticks1)-1,length(axticks1))])
    }
    scaley <- axTicks(2)[1]
  }
  lines(scalex,rep(scaley,2))
  if(rivers$units=="m") text(labx,scaley,labels=paste((scalex[2]-scalex[1])/1000,"km"),pos=3,cex=cex)
  if(rivers$units!="m") text(labx,scaley,labels=paste((scalex[2]-scalex[1]),rivers$units),pos=3,cex=cex)
}








#' Highlight Segments
#' @description Plots a river network object and displays specified segments in 
#'   bold, for easy identification.
#' @param seg A vector of segments to highlight
#' @param rivers The river network object to use
#' @param cex The character expansion factor to use for segment labels
#' @param lwd The line width to use for highlighted segments
#' @param add Whether to add the highlighted segments to an existing plot 
#'   (\code{TRUE}) or call a new plot (\code{FALSE}).  Defaults to \code{FALSE}.
#' @param color Whether to display segment labels as the same color as the 
#'   segments.  Defaults to \code{FALSE}.
#' @param ... Additional plotting arguments (see \link[graphics]{par})
#' @author Matt Tyers
#' @examples
#' data(Kenai3)
#' plot(Kenai3)
#' highlightseg(seg=c(10,30,68),rivers=Kenai3)
#' @importFrom grDevices rgb
#' @importFrom graphics plot
#' @importFrom graphics par
#' @importFrom graphics text
#' @export
highlightseg <- function(seg,rivers,cex=0.8,lwd=3,add=FALSE,color=FALSE,...) {
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
    text(midptx[j],midpty[j],j,cex=cex,col=ifelse(color,rgb((sin(j)+1)/2.3,(cos(7*j)+1)/2.3,(sin(3*(length-j))+1)/2.3),1))
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
  if(!inherits(rivers, "rivernetwork")) stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  if(!add) plot(rivers,color=F,...=...)
  connections <- rivers$connections
  lines <- rivers$lines
  connections[is.na(connections)==T]<-0
  low <- high <- rep(2,dim(connections)[2])
  coordslow <- coordshigh <- matrix(NA,nrow=length(low),ncol=2)
  con1 <- rowSums(connections==1)
  con2 <- rowSums(connections==2)
  con3 <- rowSums(connections==3)
  con4 <- rowSums(connections==4)
  con5 <- rowSums(connections==5)
  con6 <- rowSums(connections==6)
  low[con1+con2+con5+con6>0] <- 3
  high[con3+con4+con5+con6>0] <- 3
  for(i in 1:dim(connections)[2]) {
    coordslow[i,] <- lines[[i]][1,]
    coordshigh[i,] <- lines[[i]][dim(lines[[i]])[1],]
  }
  points(coordslow, pch=16, col=low)
  points(coordshigh, pch=16, col=high)
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


#' Convert XY Coordinates to River Locations
#' @description This function determines the closest vertex in the river network
#'   to each point of XY data and returns a list of river locations, defined 
#'   as segment numbers and vertex numbers.
#' @param x A vector of x-coordinates to transform
#' @param y A vector of y-coordinates to transform
#' @param rivers The river network object to use
#' @return A data frame of river locations, with segment numbers in \code{$seg}, 
#'   vertex numbers in \code{$vert}, and the snapping distance for each point in 
#'   \code{$snapdist}.  Two additional columns are \code{$snap_x} and \code{$snap_y},
#'   which give the x- and y-coordinates snapped to the river network.
#' @author Matt Tyers
#' @note Conversion to river locations is only valid if the input XY 
#'   coordinates and river network are in the same projected coordinate system. 
#'   Point data in geographic coordinates can be projected using 
#'   \link[sf]{sf_project} in package 'sf', and an example is shown below.
#' @seealso \link{pointshp2segvert}
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
#' line98albers <- sf::sf_project(pts=line98, to="+proj=aea +lat_1=55 +lat_2=65 
#'     +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs 
#'     +ellps=GRS80")
#' head(line98albers)
#' 
#' zoomtoseg(seg=c(162,19), rivers=Kenai1)
#' points(line98albers)
#' @export
xy2segvert <- function(x,y,rivers) {
  if(!inherits(rivers, "rivernetwork")) stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  if(any(is.na(x))|any(is.na(y))|!is.numeric(x)|!is.numeric(y)) stop("Missing or non-numeric coordinates.")
  
  # lengthlength <- length(unlist(lines))/2
  lengthlength <- length(unlist(rivers$lines))/2
  
  whichseg <- whichvert <- allx <- ally <- rep(NA,lengthlength)
  istart <- 1
  for(i in 1:length(rivers$lines)) {
    linelength <- dim(rivers$lines[[i]])[1]
    whichseg[istart:(istart+linelength-1)] <- i
    whichvert[istart:(istart+linelength-1)] <- 1:linelength
    allx[istart:(istart+linelength-1)] <- rivers$lines[[i]][,1]
    ally[istart:(istart+linelength-1)] <- rivers$lines[[i]][,2]
    istart <- istart+linelength
  }
  
  pdist2 <- function(p1,p2x,p2y) {
    dist <- sqrt((p1[1]-p2x)^2 + (p1[2]-p2y)^2)
    return(dist)
  }
  
  min.i <- mapply(function(x,y) which.min(pdist2(c(x,y),allx,ally))[1],x,y)  # if there's a tie, accept the first vertex
  seg <- whichseg[min.i]
  vert <- whichvert[min.i]
  snap_x <- allx[min.i]
  snap_y <- ally[min.i]
  snapdist <- sqrt((x-snap_x)^2 + (y-snap_y)^2)
  out <- data.frame(seg, vert, snapdist, snap_x, snap_y)
  return(out)
}



segvert2xy <- function(seg, vert, rivers) {
  if(!inherits(rivers, "rivernetwork")) stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  if(any(is.na(seg))|any(is.na(vert))|!is.numeric(seg)|!is.numeric(vert)) stop("Missing or non-numeric coordinates.")
  if(length(seg) != length(vert)) stop("Length mismatch between seg= and vert= arguments")
  
  snap_x <- snap_y <- rep(NA, length(seg))
  for(i in seq_along(seg)) {
    snap_x[i] <- rivers$lines[[seg[i]]][vert[i],1]
    snap_y[i] <- rivers$lines[[seg[i]]][vert[i],2]
  }
  out <- data.frame(snap_x, snap_y)
  return(out)
}




#' Draw Points from River Locations
#' @description Adds points to an active plot.  Works like \link[graphics]{points} 
#'   but with river locations (segments and vertices) rather than xy coordinates.
#' @param seg A vector of segments
#' @param vert A vector of vertices
#' @param rivers The river network object to use
#' @param pch Point character, as a vector or single value
#' @param col Point color, as a vector or single value
#' @param jitter Maximum amount of random noise to add to "jitter" points if 
#'   desired, so points do not overlap one another
#' @param ... Additional arguments for \link{points}
#' @author Matt Tyers
#' @examples
#' data(fakefish,Gulk)
#' 
#' plot(x=Gulk, xlim=c(862000,882000), ylim=c(6978000,6993000))
#' points(x=fakefish$x, y=fakefish$y, pch=16, col=2)
#' riverpoints(seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk, pch=15, col=4)
#' 
#' plot(x=Gulk, empty=TRUE)
#' with(fakefish, riverpoints(seg=seg, vert=vert, rivers=Gulk, 
#'        pch=16, col=flight, jitter=1000))
#' @importFrom graphics points
#' @importFrom stats runif
#' @export
riverpoints <- function(seg,vert,rivers,pch=1,col=1,jitter=0,...) {
  if(!inherits(rivers, "rivernetwork")) stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  lines <- rivers$lines
  if(max(seg,na.rm=T)>length(lines) | min(seg,na.rm=T)<1) stop("Invalid segment numbers specified.")
  
  x <- y <- rep(NA,length(seg))
  for(i in seq_along(seg)) {
    x[i] <- lines[[seg[i]]][vert[i],1]
    y[i] <- lines[[seg[i]]][vert[i],2]
  }
  if(jitter!=0) {
    x <- x+runif(length(x),min=-jitter,max=jitter)
    y <- y+runif(length(x),min=-jitter,max=jitter)
  }
  points(x,y,pch=pch,col=col,...=...)
}





#' Trim a River Network Object to Specified Segments
#' @description Removes line segments from a river network object.  User can 
#'   specify which segments to remove (\code{trim}) or which segments to 
#'   keep (\code{trimto}).
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
  if(!inherits(rivers, "rivernetwork")) stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
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
  trimmed.rivers$connections <- as.matrix(trimmed.rivers$connections[segs,segs])
  if(!is.na(trimmed.rivers$braided)) {
    if(trimmed.rivers$braided & !any(trimmed.rivers$connections %in% 5:6)) trimmed.rivers$braided <- NA
  }
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
      message("River mouth must be redefined - see help(setmouth)")
    }
  }
  
  if(!is.null(rivers$segroutes) | !is.null(rivers$distlookup)) {
    trimmed.rivers$segroutes <- NULL
    trimmed.rivers$distlookup <- NULL
    warning("Segment routes and/or distance lookup must be rebuilt - see help(buildsegroutes).")
  }
  trimmed.rivers <- addcumuldist(trimmed.rivers)
  
  trimmed.rivers <- update_sf(trimmed.rivers)
  
  message("Note: any point data already using the input river network must be re-transformed to river coordinates using xy2segvert() or ptshp2segvert().")
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
  if(!inherits(rivers, "rivernetwork")) stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
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
  rivers1$connections <- as.matrix(rivers1$connections[keep,keep])
  if(!is.na(rivers1$braided)) {
    if(rivers1$braided & !any(rivers1$connections %in% 5:6)) rivers1$braided <- NA
  }
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
  
  if(!is.null(rivers1$segroutes) | !is.null(rivers1$distlookup)) {
    rivers1$segroutes <- NULL
    rivers1$distlookup <- NULL
    warning("Segment routes and/or distance lookup must be rebuilt - see help(buildsegroutes).")
  }
  rivers1 <- addcumuldist(rivers1)
  
  rivers1 <- update_sf(rivers1)
  
  message("Note: any point data already using the input river network must be re-transformed to river coordinates using xy2segvert() or ptshp2segvert().")
  return(rivers1)
}


#' Remove Unconnected Segments
#' @description Detects and removes segments that are not connected to the river mouth.
#' @param rivers The river network object to use.
#' @author Matt Tyers
#' @note This function is called within \link{cleanup}, which is recommended in most cases.
#' @examples 
#' data(Koyukuk2)
#' Koy_subset <- trimriver(trimto=c(30,28,29,3,19,27,4),rivers=Koyukuk2)
#' Koy_subset <- setmouth(seg=1,vert=427,rivers=Koy_subset)
#' plot(Koy_subset)
#' 
#' Koy_subset_trim <- removeunconnected(Koy_subset)
#' plot(Koy_subset_trim)
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
  
  
  for(i in 1:length(rivers1$lines)) {
    if(is.na(detectroute(end=origin1,start=i,rivers=rivers1,stopiferror=FALSE,algorithm="sequential")[1])) {
      takeout[k] <- i
      k <- k+1
    }
  }
  
  if(!is.null(takeout)) takeout <- order[takeout]
  
  suppressMessages(rivers2 <- trimriver(trim=takeout,rivers=rivers))
  message("Note: any point data already using the input river network must be re-transformed to river coordinates using xy2segvert() or ptshp2segvert().")
  return(rivers2)
}



update_sf <- function(rivers, keep_sp=FALSE) {
  if(!inherits(rivers, "rivernetwork")) stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  
  if(!is.null(rivers$sp) & is.null(rivers$sf)) rivers <- sp2sf(rivers, keep_sp=keep_sp)  ### think if i want to keep this
  
  sfold <- rivers$sf
  sfnew <- sf::st_sf(geometry=sf::st_sfc(sf::st_multilinestring(rivers$lines)))
  sf::st_crs(sfnew) <- sf::st_crs(sfold)
  rivers$sf_current <- sfnew
  return(rivers)
}


#' @importFrom methods as
sp2sf <- function(rivers, keep_sp=FALSE) {
  if(!inherits(rivers, "rivernetwork")) stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  
  if(!is.null(rivers$sp)) {
    message('\n',"River network retains old-style dependency on sp package.")
    if(is.null(rivers$sf)) {
      rivers$sf <- as(rivers$sp, "sf")
      message('\n',"Converting to updated sf package...")
    }
    if(!keep_sp) {
      # rivers <- rivers[!names(rivers) %in% c("sp","lineID")]
      rivers$sp <- NULL
      rivers$lineID <- NULL
    }
  }
  return(rivers)
}

