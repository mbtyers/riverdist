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
#' @details 
#' Connection types. One of:
#' \itemize{
#'  \item 1: beginning - beginning
#'  \item 2: beginning - end
#'  \item 3: end - beginning
#'  \item 4: end - end
#'  \item 5: beginning - beginning and end - end
#'  \item 6: beginning - end and end - beginning
#' }
#' @examples
#' Gulk_connections <- calculateconnections(lines=Gulk$lines, tolerance=Gulk$tolerance)
#' @export
calculateconnections <- function(lines,tolerance) {
  length <- length(lines)
  # defining a connectivity matrix...
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
#' @description Uses \link[rgdal]{readOGR} in package 'rgdal' to read a river 
#'   shapefile, and establishes connectivity of segment endpoints based on 
#'   spatial proximity.
#' @param sp SpatialLinesDataFrame object. optional.
#' @param path File path, default is the current working directory.
#' @param layer Name of the shapefile, without the .shp extension.
#' @param tolerance Snapping tolerance of segment endpoints to determine 
#'   connectivity.  Default is 100, therefore care should be exercised when 
#'   working with larger units of distance, such as km.
#' @param reproject A valid Proj.4 projection string, if the shapefile is to be 
#'   re-projected.  Re-projection is done using \link[sp]{spTransform} in 
#'   package 'sp'.
#' @param supplyprojection A valid Proj.4 projection string, if the input shapefile does not have the projection information attached.
#' @return Returns an object of class \code{"rivernetwork"} containing all
#'   spatial and topological information.  See \link{rivernetwork-class}.
#' @note Since distance can only be calculated using projected coordinates, 
#'   \code{line2network()} will generate an error if a non-projected input 
#'   shapefile is detected.  To resolve this, the shapefile can be re-projected 
#'   in a GIS environment, or using \code{reproject=}, shown in the second 
#'   example below.
#' @author Matt Tyers, Joseph Stachelek
#' @importFrom rgdal readOGR
#' @importFrom sp is.projected
#' @importFrom sp CRS
#' @importFrom sp spTransform
#' @examples 
#' filepath <- system.file("extdata", package="riverdist")
#' 
#' Gulk_UTM5 <- line2network(path=filepath, layer="Gulk_UTM5")
#' plot(Gulk_UTM5)
#' 
#' ## Reading directly from a SpatialLinesDataFrame object
#' 
#' sp <- rgdal::readOGR(dsn = filepath, layer = "Gulk_UTM5", verbose = FALSE)
#' Gulk_UTM5 <- line2network(sp)
#' plot(Gulk_UTM5)
#' 
#' ## Re-projecting in Alaska Albers Equal Area projection:
#' 
#' AKalbers <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 
#'     +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
#'     
#' Gulk_AKalbers <- line2network(path=filepath, layer="Gulk_UTM5", reproject=AKalbers)
#' plot(Gulk_AKalbers)
#' 
#' @export
line2network <- function(sp = NA, path=".", layer = NA, tolerance=100, 
                         reproject=NULL, supplyprojection=NULL) {
  
  if(suppressWarnings(is.na(sp))){
    sp <- suppressWarnings(rgdal::readOGR(dsn = path, layer = layer, verbose = F))   }
  
  if(class(sp)!="SpatialLinesDataFrame"){ 
    stop("Specified shapefile is not a linear feature.")
  }
  
  if(is.na(sp@proj4string@projargs) & !is.null(supplyprojection)){ 
    sp@proj4string@projargs <- supplyprojection
  }
    
  if(is.na(sp@proj4string@projargs)){ 
    stop("Shapefile projection information is missing.  Use supplyprojection= to specify a Proj.4 projection to use.  If the input shapefile is in WGS84 geographic (long-lat) coordinates, this will be +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 (in double-quotes).  If so, it must also be reprojected using reproject=.")
  }
    
  proj4 <- strsplit(sp@proj4string@projargs,split=" ")
  projected <- sp::is.projected(sp)
  if(is.null(reproject) & !projected) stop("Distances can only be computed from a projected coordinate system.  Use reproject= to specify a Proj.4 projection to use.")
  
  if(!is.null(reproject)) {
    sp <- sp::spTransform(sp,sp::CRS(reproject))
    proj4 <- strsplit(sp@proj4string@projargs,split=" ")
  }
  
  units <- "unknown"
  for(i in 1:length(proj4[[1]])) {
    if(proj4[[1]][i]!="") {
      proj4arg <- strsplit(proj4[[1]][i],split="=")
      if(proj4arg[[1]][1]=="+units") {
        units <- proj4arg[[1]][2]
        cat('\n',"Units:",proj4arg[[1]][2],'\n')
      }
    }
  }
  
  if(length(sp@lines) > 1) {
    sp_line <- NA
    sp_seg <- NA
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
  
  out.names <- c("sp","lineID","lines","connections","lengths","names","mouth","sequenced","tolerance","units","braided","cumuldist")
  out <- list(sp,lineID,lines,connections,lengths,names,mouth,sequenced,tolerance,units,braided,cumuldist)
  names(out) <- out.names
  class(out) <- "rivernetwork"
  
  length1 <- length(out$lengths)
  suppressMessages(out <- removeduplicates(out))
  length2 <- length(out$lengths)
  if(length2<length1) cat('\n',"Removed",length1-length2,"duplicate segments.",'\n')
  suppressMessages(out <- removemicrosegs(out))
  length3 <- length(out$lengths)
  if(length3<length2) cat('\n',"Removed",length2-length3,"segments with lengths shorter than the connectivity tolerance.",'\n')
  
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
#'   \code{$seg}, vertex numbers in \code{$vert}, snapping distances in \code{$snapdist}, and the remaining columns
#'   corresponding to the data table in the input point shapefile.
#' @author Matt Tyers
#' @note If the input shapefile is detected to be in a different projection than
#'   the river network, the input shapefile will be re-projected before
#'   conversion to river locations.
#' @importFrom rgdal readOGR
#' @importFrom sp proj4string
#' @importFrom sp CRS
#' @importFrom sp spTransform
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
  shp <- rgdal::readOGR(dsn=path,layer=layer,pointDropZ=TRUE)
  if(sp::proj4string(shp) != sp::proj4string(rivers$sp)) {
    cat('\n',"Point projection detected as different from river network.  Re-projecting points before snapping to river network...")
    projection <- sp::CRS(sp::proj4string(rivers$sp))
    shp <- sp::spTransform(shp,projection) ##
  }
  segvert <- xy2segvert(x=shp@coords[,1],y=shp@coords[,2],rivers=rivers)
  outdf <- cbind(segvert,shp@data)
  return(outdf)
}


#' Plotting a River Network
#' @description S3 plotting method for the \link{rivernetwork-class}.  Produces a map of all river segments of a river network object.
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
#' @param empty Creates an empty plot if set to \code{TRUE}.  Suppresses differentiation by line type if \code{color==FALSE}, and suppresses segment number labels.  Defaults to \code{FALSE}.
#' @param linecol Line color to use if \code{empty} is \code{TRUE} or \code{color} is \code{FALSE}.  Defaults to black.
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
  if(class(x)!="rivernetwork") stop("Argument 'x' must be of class 'rivernetwork'.  See help(line2network) for more information.")
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
    if(is.na(corthing)) legendleft <- F   ######
    else legendleft <- corthing<=0   ######
    if(legendleft) {   ######
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
    if(!legendleft) {   ######
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
    if(x$units=="m") text(labx,scaley,labels=paste((scalex[2]-scalex[1])/1000,"km"),pos=3,cex=cex)
    if(x$units!="m") text(labx,scaley,labels=paste((scalex[2]-scalex[1]),x$units),pos=3,cex=cex)
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


# plotrivernetwork_OLD <- function(x,segmentnum=TRUE,offset=TRUE,lwd=1,cex=.6,scale=TRUE,color=TRUE,empty=FALSE,xlab="",ylab="",...) {
#   rivers <- x
#   if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
#   lines <- rivers$lines
#   length <- length(lines)
#   xmin <- min(lines[[1]][,1])
#   xmax <- max(lines[[1]][,1])
#   ymin <- min(lines[[1]][,2])
#   ymax <- max(lines[[1]][,2])
#   if(length>1) {
#     for(j in 2:length) {
#       if(min(lines[[j]][,1])<xmin) xmin <- min(lines[[j]][,1])
#       if(max(lines[[j]][,1])>xmax) xmax <- max(lines[[j]][,1])
#       if(min(lines[[j]][,2])<ymin) ymin <- min(lines[[j]][,2])
#       if(max(lines[[j]][,2])>ymax) ymax <- max(lines[[j]][,2])
#       # print(c(j,xmin,xmax,ymin,ymax))
#     }
#   }
#   plot(c(xmin,xmax),c(ymin,ymax),col="white",cex.axis=.6,asp=1,xlab=xlab,ylab=ylab,...=...)
#   midptx <- midpty <- NA
#   for(j in 1:length) {
#     if(!color&!empty) lines(lines[[j]],col=1,lty=j,lwd=lwd)
#     linecolor <- rgb((sin(j)+1)/2.3,(cos(7*j)+1)/2.3,(sin(3*(length-j))+1)/2.3)
#     if(color&!empty) lines(lines[[j]],col=linecolor,lty=1,lwd=lwd)
#     if(color&empty) lines(lines[[j]],col=1,lty=1,lwd=lwd)
#     linelength <- dim(lines[[j]])[1]
#     
#     xplot <- lines[[j]][,1][lines[[j]][,1]>par("usr")[1] & lines[[j]][,1]<par("usr")[2] & lines[[j]][,2]>par("usr")[3] & lines[[j]][,2]<par("usr")[4]]
#     yplot <- lines[[j]][,2][lines[[j]][,1]>par("usr")[1] & lines[[j]][,1]<par("usr")[2] & lines[[j]][,2]>par("usr")[3] & lines[[j]][,2]<par("usr")[4]]
#     if(length(xplot)>0) midptx[j] <- mean(xplot)
#     if(length(yplot)>0) midpty[j] <- mean(yplot)
#     if(length(xplot)==0 | length(yplot)==0) midptx[j] <- midpty[j] <-NA
#     middle <- floor(length(xplot)/2)
#     direction <- 3+(abs(xplot[floor(length(xplot)*.25)]-xplot[floor(length(xplot)*.75)]) < abs(yplot[floor(length(yplot)*.25)]-yplot[floor(length(yplot)*.75)]))
#     if(!offset) direction <- NULL
#     xtext <- ifelse(length(xplot)>0,xplot[middle],NA)
#     ytext <- ifelse(length(yplot)>0,yplot[middle],NA)
#     if(segmentnum&!empty) text(x=xtext,y=ytext,labels=j,pos=direction,cex=cex,col=ifelse(color,linecolor,1))
#   }
#   if(scale) {
#     if(length>1) corthing <- cor(midptx,midpty,use="complete.obs")
#     if(length==1) corthing <- (lines[[1]][1,1]-lines[[1]][dim(lines[[1]])[1],1])*(lines[[1]][1,2]-lines[[1]][dim(lines[[1]])[1],2])
#     if(corthing<=0) {
#       if(length(axTicks(1))>5) {
#         scalex <- axTicks(1)[c(1,3)]
#         labx <- axTicks(1)[2]
#       }
#       if(length(axTicks(1))<=5) {
#         scalex <- axTicks(1)[c(1,2)]
#         labx <- mean(axTicks(1)[c(1,2)])
#       }
#       scaley <- axTicks(2)[1]
#     }
#     if(corthing>0) {
#       if(length(axTicks(1))>5) {
#         scalex <- axTicks(1)[c(length(axTicks(1))-2,length(axTicks(1)))]
#         labx <- axTicks(1)[length(axTicks(1))-1]
#       }
#       if(length(axTicks(1))<=5) {
#         scalex <- axTicks(1)[c(length(axTicks(1))-1,length(axTicks(1)))]
#         labx <- mean(axTicks(1)[c(length(axTicks(1))-1,length(axTicks(1)))])
#       }
#       scaley <- axTicks(2)[1]
#     }
#     lines(scalex,rep(scaley,2))
#     if(rivers$units=="m") text(labx,scaley,labels=paste((scalex[2]-scalex[1])/1000,"km"),pos=3,cex=cex)
#     if(rivers$units!="m") text(labx,scaley,labels=paste((scalex[2]-scalex[1]),rivers$units),pos=3,cex=cex)
#   }
# }


plotrivernetwork2 <- function(x,segmentnum=TRUE,offset=TRUE,lwd=1,cex=.6,scale=TRUE,color=TRUE,empty=FALSE,xlab="",ylab="",...) {
  if(class(x)!="rivernetwork") stop("Argument 'x' must be of class 'rivernetwork'.  See help(line2network) for more information.")
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
    linecolors <- rep(1,length)
    ltys <- 1:length
  }
  if(empty) {
    ltys <- linecolors <- rep(1,length)
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
    direction[j] <- 3+(abs(xplot[floor(length(xplot)*.25)]-xplot[floor(length(xplot)*.75)]) < abs(yplot[floor(length(yplot)*.25)]-yplot[floor(length(yplot)*.75)]))
    xtext[j] <- ifelse(length(xplot)>0,xplot[middle],NA)
    ytext[j] <- ifelse(length(yplot)>0,yplot[middle],NA)
  }
  if(offset) text(xtext,ytext,labels=1:length,pos=direction,cex=cex,col=linecolors)
  else text(xtext,ytext,labels=1:length,cex=cex,col=linecolors)
  if(scale) {
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
    if(x$units=="m") text(labx,scaley,labels=paste((scalex[2]-scalex[1])/1000,"km"),pos=3,cex=cex)
    if(x$units!="m") text(labx,scaley,labels=paste((scalex[2]-scalex[1]),x$units),pos=3,cex=cex)
  }
}

# plotrivernetwork2(abstreams)
# 
# par(mfrow=c(5,5))
# system.time(plot(abstreams))
# system.time(plotrivernetwork2(abstreams))
# microbenchmark(plot(abstreams))
# microbenchmark(plotrivernetwork2(abstreams))


#' Highlight Segments
#' @description Plots a river network object and displays specified segments in bold, for easy identification.
#' @param seg A vector of segments to highlight
#' @param rivers The river network object to use
#' @param cex The character expansion factor to use for segment labels
#' @param lwd The line width to use for highlighted segments
#' @param add Whether to add the highlighted segments to an existing plot (\code{TRUE}) or call a new plot (\code{FALSE}).  Defaults to \code{FALSE}.
#' @param color Whether to display segment labels as the same color as the segments.  Defaults to \code{FALSE}.
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
  # if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  # if(!add) plot(rivers,color=F,...=...)
  # connections <- rivers$connections
  # lines <- rivers$lines
  # connections[is.na(connections)==T]<-0
  # for(i in 1:dim(connections)[2]) {
  #   low <- 2
  #   high <- 2
  #   if(length(connections[i,][connections[i,]==1 | connections[i,]==2 | connections[i,]==5 | connections[i,]==6])>0) low <- 3
  #   if(length(connections[i,][connections[i,]==3 | connections[i,]==4 | connections[i,]==5 | connections[i,]==6])>0) high <- 3
  #   points(lines[[i]][1,1],lines[[i]][1,2],pch=16,col=low)
  #   points(lines[[i]][dim(lines[[i]])[1],1],
  #          lines[[i]][dim(lines[[i]])[1],2],pch=16,col=high)
  # }
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
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
#' @return A data frame of river locations, with segment numbers in \code{$seg}, vertex numbers in \code{$vert}, and the snapping distance for each point in \code{$snapdist}.
#' @author Matt Tyers
#' @note Conversion to river locations is only valid if the input XY 
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
  if(any(is.na(x))|any(is.na(y))|!is.numeric(x)|!is.numeric(y)) stop("Missing or non-numeric coordinates.")
  
  lengthlength <- length(unlist(lines))/2
  
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
  snapdist <- sqrt((x-allx[min.i])^2 + (y-ally[min.i])^2)
  out <- data.frame(seg,vert,snapdist)
  return(out)
}

# 
# xy2segvert_OLD <- function(x,y,rivers) {
#   if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
#   if(any(is.na(x))|any(is.na(y))|!is.numeric(x)|!is.numeric(y)) stop("Missing or non-numeric coordinates.")
#   seg <- rep(NA,length(x))
#   vert <- rep(NA,length(x))
#   snapdist <- rep(NA,length(x))
#   lines <- rivers$lines
#   
#   pdist2 <- function(p1,p2mat) {
#     dist <- sqrt((p1[1]-p2mat[,1])^2 + (p1[2]-p2mat[,2])^2)
#     return(dist)
#   }
#   
#   for(i in 1:length(x)) {
#     min.dist1 <- 1000000
#     for(river.i in 1:length(lines)) {
#       dists <- pdist2(c(x[i],y[i]),lines[[river.i]])
#       min.i <- min(dists)
#       if(min.i < min.dist1) {
#         min.dist1 <- min.i
#         seg[i] <- river.i
#         vert[i] <- which(dists==min.i)[1]  # if there's a tie, accept the first vertex
#         snapdist[i] <- pdist(c(x[i],y[i]),lines[[river.i]][vert[i],])
#       }
#     }
#   }
#   out <- data.frame(seg,vert,snapdist)
#   return(out)
# }
# 
# rivtest <- Kenai3
# abpts <- xysim(5000,rivtest)
# asdf <- xy2segvert2(x=abpts$x,y=abpts$y,rivers=rivtest)
# 
# system.time(asdf1 <- xy2segvert(x=abpts$x,y=abpts$y,rivers=rivtest))
# system.time(asdf2 <- xy2segvert2(x=abpts$x,y=abpts$y,rivers=rivtest))
# 
# microbenchmark(asdf1 <- xy2segvert(x=abpts$x,y=abpts$y,rivers=rivtest))
# microbenchmark(asdf2 <- xy2segvert2(x=abpts$x,y=abpts$y,rivers=rivtest))
# all.equal(asdf1,asdf2)
# 



#' Draw Points from River Locations
#' @description Adds points to an active plot.  Works like \link[graphics]{points} but with river locations (segments and vertices) rather than xy coordinates.
#' @param seg A vector of segments
#' @param vert A vector of vertices
#' @param rivers The river network object to use
#' @param pch Point character, as a vector or single value
#' @param col Point color, as a vector or single value
#' @param jitter Maximum amount of random noise to add to "jitter" points if desired, so points do not overlap one another
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
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  lines <- rivers$lines
  if(max(seg,na.rm=T)>length(lines) | min(seg,na.rm=T)<1) stop("Invalid segment numbers specified.")
  
  x <- y <- rep(NA,length(seg))
  for(i in 1:length(seg)) {
    x[i] <- lines[[seg[i]]][vert[i],1]
    y[i] <- lines[[seg[i]]][vert[i],2]
  }
  if(jitter!=0) {
    x <- x+runif(length(x),min=-jitter,max=jitter)
    y <- y+runif(length(x),min=-jitter,max=jitter)
  }
  points(x,y,pch=pch,col=col,...=...)
}

# riverpoints_OLD <- function(seg,vert,rivers,pch=1,col=1,jitter=0,...) {
#   if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
#   if(length(pch)==1) pch <- rep(pch,length(seg))
#   if(length(col)==1) col <- rep(col,length(seg))
#   lines <- rivers$lines
#   if(max(seg,na.rm=T)>length(lines) | min(seg,na.rm=T)<1) stop("Invalid segment numbers specified.")
#   
#   for(i in 1:length(seg)) {
#     if(jitter!=0) {
#       jitterx <- runif(length(lines[[seg[i]]][vert[i],1]),min=-jitter,max=jitter)
#       jittery <- runif(length(lines[[seg[i]]][vert[i],2]),min=-jitter,max=jitter)
#     }
#     else {
#       jitterx <- 0
#       jittery <- 0
#     }
#     if(vert[i]>dim(lines[[seg[i]]])[1] | vert[i]<1) stop("Invalid vertex numbers specified.")
#     points(lines[[seg[i]]][vert[i],1] + jitterx,
#            lines[[seg[i]]][vert[i],2] + jittery, pch=pch[i],col=col[i],...=...)
#   }
# }


# 
# plot(rivtest)
# system.time(riverpoints(asdf2$seg,asdf2$vert,rivtest))
# system.time(riverpoints2(asdf2$seg,asdf2$vert,rivtest,pch=16,col=4))
# 
# microbenchmark(riverpoints(asdf2$seg,asdf2$vert,rivtest))
# microbenchmark(riverpoints2(asdf2$seg,asdf2$vert,rivtest))



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
  
  # updating sp object
  id <- rivers$lineID
  sp_lines1 <- rivers$sp@lines[unique(id[segs,2])]   
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
  
  # updating sp object
  id <- rivers$lineID
  sp_lines1 <- rivers$sp@lines[unique(id[keep,2])]  
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
