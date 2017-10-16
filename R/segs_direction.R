#' Map Segments by Name
#' @description Provides a check that river network segments were appropriately
#'   named.
#' @param rivers The river network object to use.  Function checks segment names
#'   contained in the river network object.
#' @param cex Global character expansion factor for plotting
#' @param scale Whether or not to give x- and y-axes the same scale
#' @param ... Additional plotting arguments (see \link[graphics]{par})
#' @author Matt Tyers
#' @examples
#' data(Gulk)
#' str(Gulk)
#' 
#' Gulk$names <- c("Gulkana River","Trib 1","West Fork","Gulkana River","Trib 1",
#'                 "West Fork","Trib 2","West Fork","Twelvemile Creek","Gulkana River",
#'                 "Middle Fork","Gulkana River","Middle Fork","Hungry Hollow")
#' str(Gulk)
#' 
#' mapbyname(rivers=Gulk)
#' @importFrom graphics plot
#' @importFrom graphics axTicks
#' @importFrom graphics text
#' @importFrom stats cor
#' @importFrom graphics legend
#' @importFrom grDevices rgb
#' @export
mapbyname <- function(rivers,scale=TRUE,cex=.6,...) {
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  lines <- rivers$lines
  names <- rivers$names
  if(any(is.na(names))) stop("Segment names must be specified.")
  length <- length(lines)
  names.cols <- as.numeric(as.factor(names)) 
  names.cols <- rgb((sin(names.cols)+1)/2.3,(cos(0.7*names.cols)+1)/2.3,(sin(3*(length-names.cols))+1)/2.3)  
  xmin <- min(lines[[1]][,1])
  xmax <- max(lines[[1]][,1])
  ymin <- min(lines[[1]][,2])
  ymax <- max(lines[[1]][,2])
  for(j in 2:length) {
    if(min(lines[[j]][,1])<xmin) xmin <- min(lines[[j]][,1])
    if(max(lines[[j]][,1])>xmax) xmax <- max(lines[[j]][,1])
    if(min(lines[[j]][,2])<ymin) ymin <- min(lines[[j]][,2])
    if(max(lines[[j]][,2])>ymax) ymax <- max(lines[[j]][,2])
  }
  plot(c(xmin,xmax),c(ymin,ymax),col="white",xlab="",ylab="",cex.axis=.6,asp=1,...=...)
  for(j in 1:length) {
    lines(lines[[j]],col=names.cols[j],lty=1,lwd=3)
  }
  if(scale) {
    scalex <- axTicks(1)[c(1,3)]
    scaley <- axTicks(2)[1]
    lines(scalex,rep(scaley,2))
    text(axTicks(1)[2],scaley,labels=paste((scalex[2]-scalex[1])/1000,"km"),pos=3,cex=cex)
  }
  
  legend(par("usr")[1],par("usr")[4],legend=unique(names),lty=1,lwd=3,col=unique(names.cols))
}

#' Identify Vertex Coordinates of Segment Endpoints
#' @description Identifies the vertex coordinates (row numbers) of the endpoints
#'   of a given segment.  The main purpose is determining which of the endpoints
#'   is the mouth (or lowest point) of the river system.
#' @param seg The segment (number) to check
#' @param rivers The river network object to use
#' @note This function is called within \link{cleanup}, which is recommended in most cases.
#' @author Matt Tyers
#' @examples
#' data(Gulk)
#' 
#' # say we know that segment 1 is the lowest segment in this river network, but we don't know 
#' # which end is the mouth.
#' showends(seg=1, rivers=Gulk)
#' 
#' # this means that the mouth is row 1, so we can specify this:
#' Gulk <- setmouth(seg=1, vert=1, rivers=Gulk)
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @importFrom graphics points
#' @importFrom graphics text
#' @export
showends <- function(seg, rivers) {
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  if(seg>length(rivers$lines) | seg<1) stop("Invalid segment specified.")
  
  len <- dim(rivers$lines[[seg]])[1]
  
  plot(x = rivers)
  lines(rivers$lines[[seg]], 
         lwd = 2, col = 4)
  points(rivers$lines[[seg]][1,1], 
         rivers$lines[[seg]][1,2], 
         pch = 16, col = 3)
  text(rivers$lines[[seg]][1,1], rivers$lines[[seg]][1,2], 
       labels = paste("beginning", 1), pos = 4)
  points(rivers$lines[[seg]][len,1], 
         rivers$lines[[seg]][len,2], 
         pch = 16, col = 2)
  text(rivers$lines[[seg]][len,1], 
       rivers$lines[[seg]][len,2], 
       labels = paste("end", len), pos = 4)
}

#' Specify the Segment and Vertex of the Mouth of a River Network Object.
#' @description Provides a user-friendly way of specifying the segment and
#'   vertex of the mouth (lowest point) of a river network object.
#' @param seg The segment number to store for the mouth
#' @param vert The vertex number to store for the mouth
#' @param rivers The river network object to use
#' @return A new river network object (see \link{rivernetwork})
#' @note The mouth segment and vertex can also be specified using direct
#'   assignment to the \code{$mouth$seg} and \code{$mouth$vert} components of the river network
#'   object.
#' @note This function is called within \link{cleanup}, which is recommended in most cases.
#' @author Matt Tyers
#' @examples
#' data(Gulk)
#' 
#' # say we know that segment 1 is the lowest segment in this river network, but we don't know 
#' # which end is the mouth.
#' showends(seg=1, rivers=Gulk)
#' 
#' # this means that the mouth is row 1, so we can specify this:
#' Gulk <- setmouth(seg=1, vert=1, rivers=Gulk)
#' @seealso \link{line2network}
#' @export
setmouth <- function(seg,vert,rivers) {
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  if(seg>length(rivers$lines) | seg<1) stop("Invalid segment specified.")
  if(vert>dim(rivers$lines[[seg]])[1] | vert<1) stop("Invalid vertex specified.")
  
  rivers1 <- rivers
  rivers1$mouth$mouth.seg <- seg
  rivers1$mouth$mouth.vert <- vert
  return(rivers1)
}

#' Store Vertices in Ascending Sequence
#' @description Rearranges the vertices of a river network object so that 
#'   vertices are stored sequentially moving up river for all segments 
#'   (coordinates [1,] are the bottom of each segment).
#' @param rivers The river network object to use
#' @return A new river network object (see \link{rivernetwork})
#' @author Matt Tyers
#' @note Even without calling \code{sequenceverts}, the vertices will be stored 
#'   sequentially - either moving up river or down for a given segment.  What 
#'   \code{sequenceverts()} adds is a standardized direction.
#'   
#'   Currently, no function in package 'riverdist' requires the vertices to be stored 
#'   sequentially.
#' @examples
#' data(Gulk)
#' Gulk <- setmouth(seg=1, vert=1, rivers=Gulk)
#' str(Gulk)
#' 
#' Gulk.dir <- sequenceverts(rivers=Gulk)
#' str(Gulk.dir)
#' @seealso \link{line2network}
#' @export
sequenceverts <- function(rivers) {
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  tolerance <- rivers$tolerance
  if(is.na(rivers$mouth$mouth.seg) | is.na(rivers$mouth$mouth.vert)) {
    stop("Error - Need to specify segment & vertex of origin")
  }
  length <- length(rivers$lines)
  for(i in 1:length) {
    seg.length <- dim(rivers$lines[[i]])[1]
    
    # define low vertex indices as downriver locations, i.e. vertices go up as you go upriver
    if(riverdistance(startseg=rivers$mouth$mouth.seg,endseg=i,startvert=rivers$mouth$mouth.vert,endvert=1,rivers=rivers) > 
         riverdistance(startseg=rivers$mouth$mouth.seg,endseg=i,startvert=rivers$mouth$mouth.vert,endvert=seg.length,rivers=rivers)) {
      temp <- matrix(NA,nrow=seg.length,ncol=2)
      for(j in 1:seg.length) {
        temp[j,] <- rivers$lines[[i]][(seg.length-j+1),]
      }
      rivers$lines[[i]] <- temp
      rivers$sp@lines[[rivers$lineID[i,2]]]@Lines[[rivers$lineID[i,3]]]@coords <- temp
      if(i==rivers$mouth$mouth.seg) rivers$mouth$mouth.vert<- 1
    }
  } 
  # redefining connectivity matrix
  connections <- matrix(NA,nrow=length,ncol=length)
  lines <- rivers$lines
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
      if(pdist(lines[[i]][1,],lines[[j]][1,])<tolerance & pdist(lines[[i]][i.max,],lines[[j]][j.max,])<tolerance & i!=j) {     
        connections[i,j] <- 5
      }
      if(pdist(lines[[i]][i.max,],lines[[j]][1,])<tolerance & pdist(lines[[i]][1,],lines[[j]][j.max,])<tolerance & i!=j) {
        connections[i,j] <- 6
      }    
    }
  }
  if(any(connections %in% 5:6)) rivers$braided <- TRUE
  rivers$connections <- connections
  rivers$sequenced <- TRUE
  rivers <- addcumuldist(rivers)
  return(rivers) 
}
