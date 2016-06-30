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
  
  pdist2 <- function(p1,p2mat) {
    dist <- sqrt((p1[1]-p2mat[,1])^2 + (p1[2]-p2mat[,2])^2)
    return(dist)
  }
  
  # first identifying where breaks should go
  breaks <- list()
  
  for(segi in 1:length(lines)) {
    segibreaks <- NULL
    for(segj in 1:length(lines)) {
      distanceses1 <- pdist2(lines[[segj]][1,], lines[[segi]])
      distanceses2 <- pdist2(lines[[segj]][dim(lines[[segj]])[1],], lines[[segi]])
      newbreaks <- which(distanceses1<tolerance | distanceses2<tolerance)
      segibreaks <- c(segibreaks,newbreaks)
    }
    breaks[[segi]] <- sort(unique(segibreaks))
    if(length(breaks[[segi]]>1)) {
      for(ibreaks in 2:length(breaks[[segi]])) {
        if(abs(breaks[[segi]][ibreaks]-breaks[[segi]][ibreaks-1])<=1) breaks[[segi]][ibreaks-1] <- 0 # eliminating sequential runs
        if(breaks[[segi]][ibreaks] <= ibreaks) breaks[[segi]][ibreaks] <- 0
      }
    }
    breaks[[segi]][breaks[[segi]]==1] <- 0
    breaks[[segi]][breaks[[segi]]==dim(lines[[segi]])[1]] <- 0
    breaks[[segi]] <- breaks[[segi]][breaks[[segi]]>0]
    if(length(breaks[[segi]])==0) breaks[[segi]] <- NA
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
      if(pdist(lines[[i]][1,],lines[[j]][1,])<tolerance & pdist(lines[[i]][i.max,],lines[[j]][j.max,])<tolerance & i!=j) {
        connections[i,j] <- 5
      }
      if(pdist(lines[[i]][i.max,],lines[[j]][1,])<tolerance & pdist(lines[[i]][1,],lines[[j]][j.max,])<tolerance & i!=j) {
        connections[i,j] <- 6
      }
    }
  }
  
  # making a vector of total segment lengths
  lengths <- rep(NA,length)
  for(i in 1:length) {
    lengths[i] <- pdisttot(lines[[i]])
  }
  
  # updating rivers object
  rivers$connections <- connections
  if(any(connections %in% 5:6)) rivers$braided <- TRUE
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
    rivers <- buildsegroutes(rivers,lookup=F)
  }
  rivers <- addcumuldist(rivers)
  if(!is.null(rivers$distlookup)) rivers <- buildlookup(rivers)
  
  message("Note: any point data already using the input river network must be re-transformed to river coordinates using xy2segvert() or ptshp2segvert().")
  return(rivers)
}


#'Split a Segment at a Specified Vertex
#'@description Splits a segment at a specified vertex, creating two new segments.
#'@param seg The segment to split
#'@param vert The vertex to split it at
#'@param rivers The river network object to use
#'@return A new, updated river network object
#'@seealso \link{line2network}
#'@author Matt Tyers
#'@examples
#' data(Gulk)
#' plot(x=Gulk)
#' 
#' Gulk2 <- splitsegmentat(seg=1, vert=400, rivers=Gulk)
#' plot(x=Gulk2)
#' @importFrom methods new
#'@export
splitsegmentat <- function(seg, vert, rivers) {
  lines <- rivers$lines
  lines[[length(lines)+1]] <- lines[[seg]][vert:dim(lines[[seg]])[1],]
  lines[[seg]] <- lines[[seg]][1:vert,]
  
  rivers$lines <- lines
  length <- length(lines)
  
  if(!is.na(rivers$mouth$mouth.seg) & !is.na(rivers$mouth$mouth.vert)) {
    mouthcoords <- rivers$lines[[rivers$mouth$mouth.seg]][rivers$mouth$mouth.vert,]
  }
  
  # defining a connectivity matrix...
  # connection type 1: beginning - beginning
  # connection type 2: beginning - end
  # connection type 3: end - beginning
  # connection type 4: end - end
  tolerance <- rivers$tolerance
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
      if(pdist(lines[[i]][1,],lines[[j]][1,])<tolerance & pdist(lines[[i]][i.max,],lines[[j]][j.max,])<tolerance & i!=j) {
        connections[i,j] <- 5
      }
      if(pdist(lines[[i]][i.max,],lines[[j]][1,])<tolerance & pdist(lines[[i]][1,],lines[[j]][j.max,])<tolerance & i!=j) {
        connections[i,j] <- 6
      }
    }
  }
  
  # making a vector of total segment lengths
  lengths <- rep(NA,length)
  for(i in 1:length) {
    lengths[i] <- pdisttot(lines[[i]])
  }
  
  # updating rivers object
  rivers$connections <- connections
  if(any(connections %in% 5:6)) rivers$braided <- TRUE
  rivers$lengths <- lengths
  rivers$names <- rep(NA,length(lines))
  # rivers$cumuldist <- cumuldist
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
    rivers <- buildsegroutes(rivers,lookup=F)
  }
  rivers <- addcumuldist(rivers)
  if(!is.null(rivers$distlookup)) rivers <- buildlookup(rivers)
  
  message("Note: any point data already using the input river network must be re-transformed to river coordinates using xy2segvert() or ptshp2segvert().")
  return(rivers)
}
