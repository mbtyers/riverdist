#' Dissolve
#' @description Acts like a spatial dissolve within a GIS environment.  Simplifies a river network object by combining "runs" of segments with no other connections.
#' @param rivers The river network object to use
#' @return A new river network object with segments combined
#' @seealso line2network
#' @note This function is called within \link{cleanup}, which is recommended in most cases.
#' @author Matt Tyers
#' @examples
#' data(Kenai2)
#' plot(x=Kenai2)
#' 
#' Kenai2dissolve <- dissolve(rivers=Kenai2)
#' plot(x=Kenai2dissolve)
#' @importFrom methods new
#' @export
dissolve <- function(rivers) {
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  tolerance <- rivers$tolerance
  lines <- rivers$lines
  connections <- rivers$connections
  length <- length(lines)
  if(!is.na(rivers$mouth$mouth.seg) & !is.na(rivers$mouth$mouth.vert)) {
    mouthcoords <- rivers$lines[[rivers$mouth$mouth.seg]][rivers$mouth$mouth.vert,]
  }
  
  # calculating a new connectivity matrix to capture beginning-beginning/end-end and beginning-end/end-beginning connections (special braided case)
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
  
  #internal functions
  n.top <- function(seg,connections) {
    return(length(connections[seg,][(connections[seg,]==1 | connections[seg,]==2 | connections[seg,]==5 | connections[seg,]==6) & is.na(connections[seg,])==F]))
  }
  n.bot <- function(seg,connections) {
    return(length(connections[seg,][(connections[seg,]==3 | connections[seg,]==4 | connections[seg,]==5 | connections[seg,]==6) & is.na(connections[seg,])==F]))
  }
  who.top <- function(seg,connections) {
    return(which((connections[seg,]==1 | connections[seg,]==2 | connections[seg,]==5 | connections[seg,]==6) & is.na(connections[seg,])==F))
  }
  who.bot <- function(seg,connections) {
    return(which((connections[seg,]==3 | connections[seg,]==4 | connections[seg,]==5 | connections[seg,]==6) & is.na(connections[seg,])==F))
  }
  
  # the first big algorithm, it detects "runs" of segments. this is what we ultimately want to combine.
  runs <- list(0)
  avail <- 1:length
  run.i <- 0
  
  for(seg in 1:length) {
    if(n.top(seg,connections) != 1) {
      run.i <- run.i+1
      run.done <- F
      runs[[run.i]] <- seg
      i <- 1
      nextone<-"bot"
      #take the first one out of avail
      while(run.done==F) {
        if(nextone=="top" & n.top(runs[[run.i]][[i]],connections)!=1) run.done<-T
        if(nextone=="bot" & n.bot(runs[[run.i]][[i]],connections)!=1) run.done<-T
        if(nextone=="top" & n.top(runs[[run.i]][i],connections)==1) {
          i <- i+1
          runs[[run.i]][i] <- who.top(runs[[run.i]][i-1],connections)
          # take the next one out of available
          if(n.top(runs[[run.i]][i],connections)==1) {
            if(who.top(runs[[run.i]][i],connections)==runs[[run.i]][i-1]) nextnext <- "bot"
          }
          if(n.bot(runs[[run.i]][i],connections)==1) {
            if(who.bot(runs[[run.i]][i],connections)==runs[[run.i]][i-1]) nextnext <- "top"
          }
        }
        if(nextone=="bot" & n.bot(runs[[run.i]][i],connections)==1) {
          i <- i+1
          runs[[run.i]][i] <- who.bot(runs[[run.i]][i-1],connections)
          # take the next one out of available
          if(n.top(runs[[run.i]][i],connections)==1) {
            if(who.top(runs[[run.i]][i],connections)==runs[[run.i]][i-1]) nextnext <- "bot"
          }
          if(n.bot(runs[[run.i]][i],connections)==1) {
            if(who.bot(runs[[run.i]][i],connections)==runs[[run.i]][i-1]) nextnext <- "top"
          }
        }
        if(!run.done) nextone <- nextnext
      }
    }
  }
  
  # the second big algorithm, it uses the information from the "runs" detected and combines where applicable
  newlines <- list()
  for(run.i in 1:length(runs)) {
    newlines[[run.i]] <- lines[[runs[[run.i]][1]]]
    if(length(runs[[run.i]]) > 1) {
      for(i in 2:length(runs[[run.i]])) {
        if(connections[runs[[run.i]][i-1],runs[[run.i]][i]]==1 | connections[runs[[run.i]][i-1],runs[[run.i]][i]]==3) {
          newlines[[run.i]] <- rbind(newlines[[run.i]],lines[[runs[[run.i]][i]]])
        }
        if(connections[runs[[run.i]][i-1],runs[[run.i]][i]]==2 | connections[runs[[run.i]][i-1],runs[[run.i]][i]]==4) {
          nextline <- NA*lines[[runs[[run.i]][i]]]
          for(j in 1:dim(nextline)[1]) {
            nextline[j,] <- lines[[runs[[run.i]][i]]][(dim(nextline)[1]-j+1),]
          }
          newlines[[run.i]] <- rbind(newlines[[run.i]],nextline)
        }
      }
    }
  }
  
  # updating the connectivity matrix with the new segments
  length <- length(newlines)
  connections <- matrix(NA,nrow=length,ncol=length)
  for(i in 1:length) {
    for(j in 1:length) {
      i.max <- dim(newlines[[i]])[1]
      j.max <- dim(newlines[[j]])[1]
      if(pdist(newlines[[i]][1,],newlines[[j]][1,])<tolerance & i!=j) {
        connections[i,j] <- 1
      }
      if(pdist(newlines[[i]][1,],newlines[[j]][j.max,])<tolerance & i!=j) {
        connections[i,j] <- 2
      }
      if(pdist(newlines[[i]][i.max,],newlines[[j]][1,])<tolerance & i!=j) {
        connections[i,j] <- 3
      }
      if(pdist(newlines[[i]][i.max,],newlines[[j]][j.max,])<tolerance & i!=j) {
        connections[i,j] <- 4
      }
    }
  }
  
  # updating lengths
  lengths <- rep(NA,length)
  for(i in 1:length) {
    # sum<-0
    # linelength <- dim(newlines[[i]])[1]
    # for(j in 1:(linelength-1)) {
    #   sum <- sum+pdist(newlines[[i]][j,],newlines[[i]][(j+1),])
    # }
    # lengths[i]<-sum
    lengths[i] <- pdisttot(newlines[[i]])
  }
  
  # cumuldist <- list()
  # for(i in 1:length(newlines)) {
  #   xy <- newlines[[i]]
  #   n <- dim(xy)[1]
  #   cumuldist[[i]] <- c(0,cumsum(sqrt(((xy[1:(n-1),1] - xy[2:n,1])^2) + ((xy[1:(n-1),2] - xy[2:n,2])^2))))
  # }
  
  #updating rivers object
  rivers1 <- rivers
  #rivers1$mouth$mouth.seg <- NA
  #rivers1$mouth$mouth.vert <- NA
  rivers1$sequenced <- F
  rivers1$names <- rep(NA,length)
  rivers1$lines <- newlines
  # rivers1$cumuldist <- cumuldist
  if(!is.na(rivers$mouth$mouth.seg) & !is.na(rivers$mouth$mouth.vert)) {
    for(i in 1:length(rivers1$lines)) {
      for(j in 1:(dim(rivers1$lines[[i]])[1])) {
        if(all(mouthcoords==rivers1$lines[[i]][j,])) {
          rivers1$mouth$mouth.seg <- i
          rivers1$mouth$mouth.vert <- j
        }
      }
    }
  }
  rivers1$connections <- connections
  rivers1$lengths <- lengths
  
  Id <- 0
  rivers1$sp@data <- data.frame(Id)
  rivers1$sp@lines <- list(rivers$sp@lines[[1]])
  rivers1$sp@lines[[1]]@Lines <- list(rivers1$sp@lines[[1]]@Lines[[1]])
  for(i in 1:length) {
    rivers1$sp@lines[[1]]@Lines[[i]] <- new("Line",coords = rivers1$lines[[i]])
  }
  
  rivID <- 1:length
  sp_line <- rep(1,length)
  sp_seg <- 1:length
  rivers1$lineID <- data.frame(rivID,sp_line,sp_seg)
  
  if(!is.null(rivers1$segroutes)) {
    # rivers1$segroutes <- NULL
    # warning("Segment routes must be rebuilt - see help(buildsegroutes).")
    rivers1 <- buildsegroutes(rivers1,lookup=F)
  }
  rivers1 <- addcumuldist(rivers1)
  if(!is.null(rivers1$distlookup)) rivers1 <- buildlookup(rivers1)
  
  return(rivers1)
}


#' Zoom to segment
#' @description Calls \link{plot.rivernetwork} and automatically zooms to a specified
#'   segment or vector of segments.  Not intended for any real mapping - just
#'   investigating and error checking.
#' @param seg A segment or vector of segments to zoom to
#' @param rivers The river network object to use
#' @param ... Additional plotting arguments (see \link[graphics]{par})
#' @author Matt Tyers
#' @examples
#' data(Kenai3)
#' plot(x=Kenai3)
#' 
#' # checking out a particularly messy region...
#' zoomtoseg(c(109,55), rivers=Kenai3)
#' @importFrom graphics plot
#' @export
zoomtoseg <- function(seg,rivers,...) {
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  if(max(seg,na.rm=T)>length(rivers$lines) | min(seg,na.rm=T)<1) stop("Invalid segment numbers specified.")
  coords <- rivers$lines[[seg[1]]]
  if(length(seg)>1) {
    for(i in 2:length(seg)) {
      coords <- rbind(coords,rivers$lines[[seg[i]]])
    }
  }
  x <- coords[,1]
  y <- coords[,2]
  minx <- min(x)
  maxx <- max(x)
  miny <- min(y)
  maxy <- max(y)
  plot(x=rivers,xlim=c(minx,maxx),ylim=c(miny,maxy),...=...)
}