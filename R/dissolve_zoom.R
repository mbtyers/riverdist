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
  run.i <- 0
  used <- rep(F,length)
  
  for(seg in 1:length) {
    if(((n.top(seg,connections) != 1) | (n.bot(seg,connections) != 1)) & !used[seg]) {
      used[seg] <- T 
      run.i <- run.i+1
      run.done <- F
      runs[[run.i]] <- seg
      i <- 1
      if(n.top(seg,connections) != 1) nextone<-"bot"
      if(n.bot(seg,connections) != 1) nextone<-"top"
      #take the first one out of avail 
      while(!run.done) {
        if(nextone=="top" & n.top(runs[[run.i]][[i]],connections)!=1) run.done<-T    
        if(nextone=="bot" & n.bot(runs[[run.i]][[i]],connections)!=1) run.done<-T
        if(nextone=="top" & n.top(runs[[run.i]][i],connections)==1) {
          i <- i+1
          runs[[run.i]][i] <- who.top(runs[[run.i]][i-1],connections)
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
          if(n.top(runs[[run.i]][i],connections)==1) {
            if(who.top(runs[[run.i]][i],connections)==runs[[run.i]][i-1]) nextnext <- "bot"
          }
          if(n.bot(runs[[run.i]][i],connections)==1) {
            if(who.bot(runs[[run.i]][i],connections)==runs[[run.i]][i-1]) nextnext <- "top"
          }
        }
        if(!run.done) nextone <- nextnext
        if(run.done) used[runs[[run.i]][length(runs[[run.i]])]] <- T
      }
    }
  }
  
  # the second big algorithm, it uses the information from the "runs" detected and combines where applicable
  newlines <- list()
  for(run.i in 1:length(runs)) {
    newlines[[run.i]] <- lines[[runs[[run.i]][1]]]
    if(length(runs[[run.i]]) > 1) {
      for(i in 2:length(runs[[run.i]])) {
        if(connections[runs[[run.i]][i-1],runs[[run.i]][i]]==1) {  
          len <- dim(newlines[[run.i]])[1]
          newlines[[run.i]] <- rbind(newlines[[run.i]][len:1,],lines[[runs[[run.i]][i]]])
        }
        if(connections[runs[[run.i]][i-1],runs[[run.i]][i]]==3) {
          newlines[[run.i]] <- rbind(newlines[[run.i]],lines[[runs[[run.i]][i]]])
        }
        if(connections[runs[[run.i]][i-1],runs[[run.i]][i]]==2) {
          len <- dim(newlines[[run.i]])[1]
          newlines[[run.i]] <- rbind(lines[[runs[[run.i]][i]]],newlines[[run.i]])
        }
        if(connections[runs[[run.i]][i-1],runs[[run.i]][i]]==4) {
          len <- dim(newlines[[run.i]])[1]
          newlines[[run.i]] <- rbind(lines[[runs[[run.i]][i]]],newlines[[run.i]][len:1,])
        }
      }
    }
  }
  
  # updating the connectivity matrix with the new segments
  length <- length(newlines)
  connections <- calculateconnections(lines=newlines, tolerance=tolerance)
  
  # updating lengths
  lengths <- rep(NA,length)
  for(i in 1:length) {
    lengths[i] <- pdisttot(newlines[[i]])
  }
  
  #updating rivers object
  rivers1 <- rivers
  rivers1$sequenced <- F
  rivers1$names <- rep(NA,length)
  rivers1$lines <- newlines
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
  if(any(connections %in% 5:6)) rivers1$braided <- TRUE
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
  
  if(!is.na(rivers1$mouth$mouth.seg) & !is.na(rivers1$mouth$mouth.vert)) {
    if(rivers1$mouth$mouth.vert > 1 & rivers1$mouth$mouth.vert < dim(rivers1$lines[[rivers1$mouth$mouth.seg]])[1]) {
      suppressMessages(rivers1 <- splitsegmentat(seg=rivers1$mouth$mouth.seg, vert=rivers1$mouth$mouth.vert, rivers=rivers1))
    }
  } 
  
  if(!is.null(rivers1$segroutes)) {
    rivers1 <- buildsegroutes(rivers1,lookup=F)
  }
  rivers1 <- addcumuldist(rivers1)
  if(!is.null(rivers1$distlookup)) rivers1 <- buildlookup(rivers1)
  
  message("Note: any point data already using the input river network must be re-transformed to river coordinates using xy2segvert() or ptshp2segvert().")
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
#' zoomtoseg(c(110,63), rivers=Kenai3)
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
