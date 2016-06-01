#'Split Segments by Endpoint Proximity
#'@description Detects cases in which segments should be split to establish
#'  appropriate topology, and splits them.  Specifically, it looks for segment
#'  endpoints intersecting (or within a tolerance of) another
#'  segment.  It then splits the intersected segment at the point where the
#'  endpoint of the other segment breaks it.
#'@param rivers The river network object to use
#'@param tolerance The spatial snapping tolerance to use for detecting
#'  intersection.  If a NULL value is used (default), it will default to the
#'  tolerance that was used in river network creation in \link{line2network}.
#'@return A new, updated river network object
#'@seealso \link{line2network}
#' @note This function is called within \link{cleanup}, which is recommended in most cases.
#'@author Matt Tyers
#'@examples
#' data(Koyukuk1)
#' topologydots(rivers=Koyukuk1)
#' # Segments 7, 8, 13, and 16 need to be split so topologies will work.  
#' # Since endpoints are not in the same place, they are not detected as 
#' # being connected.
#' plot(x=Koyukuk1)
#'
#' Koyukuk1split <- splitsegments(rivers=Koyukuk1)
#' topologydots(rivers=Koyukuk1split)
#' plot(x=Koyukuk1split)
#' @importFrom methods new
#'@export
splitsegments <- function(rivers,tolerance=NULL) {
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  lines <- rivers$lines
  if(is.null(tolerance)) tolerance <- rivers$tolerance
  
  if(!is.na(rivers$mouth$mouth.seg) & !is.na(rivers$mouth$mouth.vert)) {
    mouthcoords <- rivers$lines[[rivers$mouth$mouth.seg]][rivers$mouth$mouth.vert,]
  }
  
  # first identifying where breaks should go
  breaks <- list()
  for(riv.i in 1:length(lines)) {
    breaks.i <- 0
    breaks[[riv.i]] <- NA
    for(seg.i in 1:dim(lines[[riv.i]])[1]) {
      for(riv.j in 1:length(lines)) {
        if((pdist(lines[[riv.i]][seg.i,],lines[[riv.j]][1,])<tolerance | 
            pdist(lines[[riv.i]][seg.i,],lines[[riv.j]][dim(lines[[riv.j]])[1],])<tolerance) & riv.i!=riv.j) {
          breaks.i <- breaks.i+1
          breaks[[riv.i]][breaks.i] <- seg.i
        }
      }
    }
    if(breaks.i>1) {
      for(i in 2:breaks.i) {
        if(abs(breaks[[riv.i]][i]-breaks[[riv.i]][i-1])<=1) breaks[[riv.i]][i-1]<-0  # eliminating sequential runs
        if(breaks[[riv.i]][i] <= i) breaks[[riv.i]][i] <- 0
      }
    }
    breaks[[riv.i]][breaks[[riv.i]]==1] <- 0
    breaks[[riv.i]][breaks[[riv.i]]==dim(lines[[riv.i]])[1]] <- 0
    breaks[[riv.i]] <- breaks[[riv.i]][breaks[[riv.i]]>0]
    if(length(breaks[[riv.i]])==0) breaks[[riv.i]] <- NA
  }
  
  # then breaking them into new segments
  newlines <- list()
  new.i <- 1
  for(i in 1:length(lines)) {
    if(is.na(breaks[[i]][1])) {
      newlines[[new.i]] <- lines[[i]]
      new.i <- new.i+1
    } 
    if(!is.na(breaks[[i]][1])) {
      newlines[[new.i]] <- lines[[i]][1:breaks[[i]][1],]
      new.i <- new.i+1
      if(length(breaks[[i]])>1) {
        for(j in 1:(length(breaks[[i]])-1)) {
          newlines[[new.i]] <- lines[[i]][breaks[[i]][j]:breaks[[i]][j+1],]
          new.i <- new.i+1
        }
      }
      newlines[[new.i]] <- lines[[i]][breaks[[i]][length(breaks[[i]])]:dim(lines[[i]])[1],]
      new.i <- new.i+1
    }
  }
  
  # updating connections matrix
  lines <- newlines
  length <- length(lines)
  
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
  
  # updating rivers object
  rivers$connections <- connections
  rivers$lines <- newlines
  rivers$lengths <- lengths
  rivers$names <- rep(NA,length(lines))
  if(!is.na(rivers$mouth$mouth.seg) & !is.na(rivers$mouth$mouth.vert)) {
    for(i in 1:length(rivers$lines)) {
      for(j in 1:(dim(rivers$lines[[i]])[1])) {
        if(all(mouthcoords==rivers$lines[[i]][j,])) {
          rivers$mouth$mouth.seg <- i
          rivers$mouth$mouth.vert <- j
        }
      }
    }
  }
  
  Id <- 0
  rivers$sp@data <- data.frame(Id)
  rivers$sp@lines <- list(rivers$sp@lines[[1]])
  rivers$sp@lines[[1]]@Lines <- list(rivers$sp@lines[[1]]@Lines[[1]])
  for(i in 1:length) {
    rivers$sp@lines[[1]]@Lines[[i]] <- new("Line",coords = rivers$lines[[i]])
  }
  
  rivID <- 1:length
  sp_line <- rep(1,length)
  sp_seg <- 1:length
  rivers$lineID <- data.frame(rivID,sp_line,sp_seg)
  
  if(!is.null(rivers$segroutes)) {
    rivers$segroutes <- NULL
    warning("Segment routes must be rebuilt - see help(buildsegroutes).")
  }
  
  return(rivers)
}