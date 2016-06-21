#' Remove Duplicates
#' @description Removes duplicated line segments, which can sometimes exist within a shapefile. 
#' @param rivers The river network object to use
#' @return A new river network object with duplicated segments removed, see \link{rivernetwork}
#' @seealso line2network
#' @author Matt Tyers
#' @examples
#' data(abstreams0)
#' zoomtoseg(seg=c(170,171,157),rivers=abstreams0)
#' 
#' abstreams1 <- removeduplicates(rivers=abstreams0)
#' zoomtoseg(seg=c(166,167,154),rivers=abstreams1)
#' @export
removeduplicates <- function(rivers) {
  trim.i <- 1
  trim <- NULL
  for(i in 1:length(rivers$lines)) {
    for(j in 1:length(rivers$lines)) {
      if(i!=j&length(rivers$lines[[i]])==length(rivers$lines[[j]])) {
        if(all(rivers$lines[[i]]==rivers$lines[[j]])) {
          trim[trim.i] <- max(c(i,j))
          trim.i <- trim.i+1
        }
      }
    }
  }
  trim <- unique(trim)
  rivers1 <- trimriver(trim=trim,rivers=rivers)
  # message("Note: any point data already using the input river network must be re-transformed to river coordinates using xy2segvert() or ptshp2segvert().")
  return(rivers1)
}


#' Interactive Cleanup of a River Network
#' @description This is the recommended function to use for cleanup of a river
#'   network.  It calls all available river network editing functions in
#'   appropriate sequence, detecting which are needed or recommended, and
#'   prompts user input wherever necessary.
#'   
#'   Currently, it automatically calls \link{removeduplicates}, prompts the user
#'   whether to run \link{dissolve}, automatically runs \link{removemicrosegs}
#'   and \link{splitsegments} if needed, provides user prompts for
#'   \link{addverts} and \link{setmouth}, detects if segments are unconnected
#'   and provides user prompts for \link{removeunconnected} or
#'   \link{connectsegs}, automatically runs \link{checkbraidedTF}, and prompts
#'   the user whether to run \link{buildsegroutes} if no braiding is detected.
#' @param rivers The river network object to use
#' @return A new river network object with duplicated segments removed, see
#'   \link{rivernetwork}
#' @seealso line2network
#' @importFrom graphics hist
#' @author Matt Tyers
#' @examples
#' data(abstreams0,Koyukuk0,Kenai1)
#' 
#' # abstreams_fixed <- cleanup(abstreams0)
#' # Koyukuk <- cleanup(Koyukuk0)
#' # Kenai <- cleanup(Kenai1)
#' @importFrom graphics plot
#' @importFrom graphics text
#' @export
cleanup <- function(rivers) {
  if(!interactive()) stop("The cleanup() function can only be used in an interactive environment.")
  cat("Cleanup started, with",length(rivers$lines),"segments.",'\n')
  plot(rivers)
  cat('\n',"Removing duplicate line segments...",'\n')
  rivers1 <- removeduplicates(rivers=rivers)
  cat("Removed",(length(rivers$lines)-length(rivers1$lines)),"duplicated segments.",'\n')
  if((length(rivers$lines)-length(rivers1$lines))>0) plot(rivers1)
  cat('\n',"Checking if dissolve is recommended...",'\n')
  
  tolerance <- rivers1$tolerance
  lines <- rivers1$lines
  connections <- rivers1$connections
  length <- length(lines)
  
  # calculating a new connectivity matrix to capture beginning-beginning/end-end and beginning-end/end-beginning connections (special braided case)
  # for(i in 1:length) {
  #   for(j in 1:length) {
  #     i.max <- dim(lines[[i]])[1]
  #     j.max <- dim(lines[[j]])[1]
  #     if(pdist(lines[[i]][1,],lines[[j]][1,])<tolerance & i!=j) {
  #       connections[i,j] <- 1
  #     }
  #     if(pdist(lines[[i]][1,],lines[[j]][j.max,])<tolerance & i!=j) {
  #       connections[i,j] <- 2
  #     }
  #     if(pdist(lines[[i]][i.max,],lines[[j]][1,])<tolerance & i!=j) {
  #       connections[i,j] <- 3
  #     }
  #     if(pdist(lines[[i]][i.max,],lines[[j]][j.max,])<tolerance & i!=j) {
  #       connections[i,j] <- 4
  #     }
  #     if(pdist(lines[[i]][1,],lines[[j]][1,])<tolerance & pdist(lines[[i]][i.max,],lines[[j]][j.max,])<tolerance & i!=j) {
  #       connections[i,j] <- 5
  #     }
  #     if(pdist(lines[[i]][i.max,],lines[[j]][1,])<tolerance & pdist(lines[[i]][1,],lines[[j]][j.max,])<tolerance & i!=j) {
  #       connections[i,j] <- 6
  #     }
  #   }
  # }
  n.top <- function(seg,connections) {
    return(length(connections[seg,][(connections[seg,]==1 | connections[seg,]==2 | connections[seg,]==5 | connections[seg,]==6) & is.na(connections[seg,])==F]))
  }
  n.bot <- function(seg,connections) {
    return(length(connections[seg,][(connections[seg,]==3 | connections[seg,]==4 | connections[seg,]==5 | connections[seg,]==6) & is.na(connections[seg,])==F]))
  }
  i<-1
  done<-F
  needed<-F
  while(!done) {
    if(n.top(i,connections)==1 | n.bot(i,connections)==i) {
      needed<-T
      done<-T
    }
    if(i==length) done<-T
    i<-i+1
  }
  
  if(!needed) {
    rivers2<-rivers1
  }
  if(needed) {
    cat("Dissolve recommended. This will combine segments that do not need to be split, but will remove any data stored as a table.",'\n')
    yes<-0
    while(!any(yes==c("y","Y","n","N"))) yes <- readline(prompt="Dissolve? (y/n) ")
    if(yes=="Y" | yes=="y") {
      cat("Dissolving...",'\n')
      suppressMessages(rivers2 <- dissolve(rivers=rivers1))
      cat("Simplified from",length(rivers1$lines),"to",length(rivers2$lines),"segments.",'\n')
      plot(rivers2)
    }
    if(yes!="Y" & yes!="y") rivers2 <- rivers1
  }
  
  cat('\n',"Checking for segments with length less than the connectivity tolerance...",'\n')
  displacement <- NA
  for(i in 1:length(rivers2$lines)) {
    dim.i <- dim(rivers2$lines[[i]])[1]
    displacement[i] <- pdist(rivers2$lines[[i]][1,],rivers2$lines[[i]][dim.i,])
  }
  problems <- NULL
  problems <- (1:length(rivers2$lines))[displacement<=rivers2$tolerance]
  cat('\n',length(problems),"microsegments identified.",'\n')
  if(length(problems)>0) {
    cat('\n',"Removing microsegments...",'\n')
    suppressMessages(rivers3 <- removemicrosegs(rivers2) ) 
    plot(rivers3)
    yes<-0
    while(!any(yes==c("y","Y","n","N"))) yes <- readline(prompt="(Re)dissolve? (y/n) ")
    if(yes=="Y" | yes=="y") {
      cat("Dissolving...",'\n')
      suppressMessages(rivers3a <- dissolve(rivers=rivers3))
      plot(rivers3a)
      cat("Simplified from",length(rivers3$lines),"to",length(rivers3a$lines),"segments.",'\n')
      rivers3 <- rivers3a
    }
  }
  if(length(problems)==0) rivers3<-rivers2
  
  cat('\n',"Checking if splitting segments is needed...",'\n')
  
  pdist2 <- function(p1,p2mat) {
    dist <- sqrt((p1[1]-p2mat[,1])^2 + (p1[2]-p2mat[,2])^2)
    return(dist)
  }
  
  riv.i <- 1
  done <- F
  needed <- F
  while(!done) {
    if(n.top(riv.i,rivers3$connections)==0) {
      for(riv.j in 1:length(rivers3$lines)) {
        # for(vert.j in 1:(dim(rivers3$lines[[riv.j]])[1])) {
        #   if(riv.i!=riv.j & pdist(rivers3$lines[[riv.i]][1,],rivers3$lines[[riv.j]][vert.j,]) < rivers3$tolerance) {
        #     needed <- T
        #     done <- T
        #   }
        # }
        if(riv.i!=riv.j) {
          distanceses <- pdist2(rivers3$lines[[riv.i]][1,],rivers3$lines[[riv.j]])
          if(any(distanceses<rivers3$tolerance)) {
            needed <- T
            done <- T
          }
        }
      }
    }
    if(n.bot(riv.i,rivers3$connections)==0) {
      for(riv.j in 1:length(rivers3$lines)) {
        # for(vert.j in 1:(dim(rivers3$lines[[riv.j]])[1])) {
        #   if(riv.i!=riv.j & pdist(rivers3$lines[[riv.i]][dim(rivers3$lines[[riv.i]])[1],],rivers3$lines[[riv.j]][vert.j,]) < rivers3$tolerance) {
        #     needed <- T
        #     done <- T
        #   }
        # }
        if(riv.i!=riv.j) {
          distanceses <- pdist2(rivers3$lines[[riv.i]][dim(rivers3$lines[[riv.i]])[1],],rivers3$lines[[riv.j]])
          if(any(distanceses<rivers3$tolerance)) {
            needed <- T
            done <- T
          }
        }
      }
    }
    if(riv.i==length(rivers3$lines)) done<-T
    riv.i <- riv.i+1
  }
  
  if(!needed) {
    rivers4<-rivers3
  }
  if(needed) {
    cat("Segments must be split for connectedness to be correct.",'\n')
    cat("Splitting segments...",'\n')
    suppressMessages(rivers4 <- splitsegments(rivers=rivers3))
    cat("Identified",(length(rivers4$lines)-length(rivers3$lines)),"new segment breaks.",'\n')
  }
  
  # checking segment lengths
  distances <- NULL
  for(segi in 1:length(rivers4$lines)) {
    if(dim(rivers4$lines[[segi]])[1]>1) {
      for(verti in 2:(dim(rivers4$lines[[segi]])[1])) {
        distances <- c(distances,pdist(rivers4$lines[[segi]][verti,],rivers4$lines[[segi]][(verti-1),]))
      }
    }
  }
  maxlength <- max(distances,na.rm=T)
  cat("Maximum distance between vertices is",round(maxlength),'\n')
  hist(distances)
  yes <- 0
  while(!any(yes==c("y","Y","n","N"))) yes <- readline(prompt="Insert vertices to reduce distances and increase point snapping precision? (y/n) ")
  if(yes=="Y" | yes=="y") {
    mindist <- as.numeric(readline(prompt="Minimum distance to use: "))
    cat("Inserting vertices...",'\n')
    suppressMessages(rivers4 <- addverts(rivers4, mindist=mindist))
  }
  
  if(is.na(rivers4$mouth$mouth.seg) | is.na(rivers4$mouth$mouth.vert)) {
    cat('\n')
    plot(x=rivers4,cex=.8)
    mouthset <- FALSE
    while(!mouthset) {   
      mouthseg <- as.numeric(readline(prompt="Please identify segment number of river mouth: "))
      if(mouthseg>0 & mouthseg<=length(rivers4$lines)) {
        showends(seg=as.numeric(mouthseg),rivers=rivers4)
        mouthvert <- as.numeric(readline(prompt="Please identify vertex number of river mouth: "))
        if(mouthvert==1 | mouthvert==(dim(rivers4$lines[[mouthseg]])[1])) {
          plot(rivers4)
          riverpoints(seg=mouthseg,vert=mouthvert,pch=16,rivers=rivers4)
          text(x=rivers4$lines[[mouthseg]][mouthvert,1],y=rivers4$lines[[mouthseg]][mouthvert,2],labels="mouth",pos=4)
          accept<-0
          while(!any(accept==c("y","Y","n","N"))) accept <- readline(prompt="Accept mouth assignment? (y/n) ")
          if(accept=="Y" | accept=="y") mouthset <- TRUE
        }
      }
    }
    rivers4 <- setmouth(seg=as.numeric(mouthseg),vert=as.numeric(mouthvert),rivers=rivers4)
  }
  
  dealtwith<-F
  while(!dealtwith) {
    cat('\n',"Checking for unconnected segments...",'\n')
    takeout <- NULL
    k <- 1
    
    dists <- rep(NA,length(rivers4$lines))
    for(i in 1:length(rivers4$lines)) {
      dists[i] <- pdist(rivers4$lines[[rivers4$mouth$mouth.seg]][rivers4$mouth$mouth.vert,],
                        c(mean(rivers4$lines[[i]][,1]),mean(rivers4$lines[[i]][,2])))
    }
    order <- order(dists)
    rivers41 <- rivers4
    rivers41$connections <- rivers4$connections[order,order]
    rivers41$lines <- rivers4$lines[order]
    origin <- rivers4$mouth$mouth.seg
    origin41 <-which(order==origin)
    
    for(i in 1:length(rivers41$lines)) {
      if(is.na(detectroute(end=origin41,start=i,rivers=rivers41,stopiferror=FALSE,algorithm="sequential")[1])) {
        takeout[k] <- i
        k <- k+1
      }
    }
    if(!is.null(takeout)) {
      takeout <- order[takeout]
      highlightseg(seg=takeout,rivers=rivers4)
      cat(length(takeout),"unconnected segments detected.")
      whattodo<-0
      while(!any(whattodo==c("r","R","c","C"))) {
        whattodo <- readline(prompt="Please select - (r)emove unconnected segments or (c)onnect segments & check again (r/c) ")
      }
      if(whattodo=="r" | whattodo=="R") {
        suppressMessages(rivers5 <- trimriver(trim=takeout,rivers=rivers4))
        dealtwith<-T
        cat("No unconnected segments detected.",'\n')
      }
      if(whattodo=="c" | whattodo=="C") {
        plot(rivers4)
        connect1 <- as.numeric(readline("Enter the number of the segment you'd like to connect: "))
        connect2 <- as.numeric(readline("Enter the number of the segment you'd like to connect it to: "))
        howtodo <- 0
        while(!any(howtodo==c("c","C","e","E"))) {
          howtodo <- readline(prompt="Connect at (e)ndpoint or (c)losest point? (e/c) ")
        }
        if(any(howtodo==c("c","C"))) closestpt <- T
        if(any(howtodo==c("e","E"))) closestpt <- F
        if(closestpt) cat('\n',"Connecting and calculating new segment splits...",'\n')
        if(!closestpt) cat('\n',"Connecting...",'\n')
        suppressMessages(rivers4 <- connectsegs(connect=connect1,connectto=connect2,nearestvert=closestpt,rivers=rivers4))
      }
    }
    if(is.null(takeout)) {
      rivers5<-rivers4
      dealtwith<-T
      cat("No unconnected segments detected.",'\n')
    }
  }
  
  plot(rivers5)
  yes <- "Y"
  while(yes=="Y" | yes=="y") {
    yes<-0
    while(!any(yes==c("y","Y","n","N"))) yes<-readline(prompt="Remove any additional segments? (y/n) ")   
    if(yes=="Y" | yes=="y") {
      whichones <- readline("Enter segments to remove, separated by commas: ")
      whichones <- as.numeric(unlist(strsplit(whichones, ",")))
      suppressMessages(rivers5a <- trimriver(trim=whichones,rivers=rivers5))
      plot(rivers5a)
      accept<-0
      while(!any(accept==c("y","Y","n","N"))) accept <- readline(prompt="Accept changes? (y/n) ")
      if(accept=="Y" | accept=="y") rivers5 <- rivers5a
      if(accept!="Y" & accept!="y") plot(rivers5)
      if(is.na(rivers5$mouth$mouth.seg)) {
        mouthset <- FALSE
        while(!mouthset) {    
          mouthseg <- as.numeric(readline(prompt="Please identify new segment number of river mouth: "))
          if(mouthseg>0 & mouthseg<=length(rivers5$lines)) {
            showends(seg=as.numeric(mouthseg),rivers=rivers5)
            mouthvert <- as.numeric(readline(prompt="Please identify new vertex number of river mouth: "))
            if(mouthvert==1 | mouthvert==(dim(rivers5$lines[[mouthseg]])[1])) {
              plot(rivers5)
              riverpoints(seg=mouthseg,vert=mouthvert,pch=16,rivers=rivers5)
              text(x=rivers5$lines[[mouthseg]][mouthvert,1],y=rivers5$lines[[mouthseg]][mouthvert,2],labels="mouth",pos=4)
              accept<-0
              while(!any(accept==c("y","Y","n","N"))) accept <- readline(prompt="Accept mouth assignment? (y/n) ")
              if(accept=="Y" | accept=="y") mouthset <- TRUE
            }
          }
        }
        rivers5 <- setmouth(seg=as.numeric(mouthseg),vert=as.numeric(mouthvert),rivers=rivers5)
      }
    }
  }
  
  cat('\n',"Checking for braiding...",'\n')
  rivers5<-checkbraidedTF(rivers5)
  braided <- rivers5$braided
  rivers6<-rivers5
  if(braided) {
    cat("Braiding detected within river network.",'\n')
    yes <- "Y"
    while(yes=="Y" | yes=="y") {
      yes<-0
      while(!any(yes==c("y","Y","n","N"))) yes<-readline(prompt="Remove any additional segments (y/n) ")  
      if(yes=="Y" | yes=="y") {
        whichones <- readline("Enter segments to remove, separated by commas: ")
        whichones <- as.numeric(unlist(strsplit(whichones, ",")))
        suppressMessages(rivers6a <- trimriver(trim=whichones,rivers=rivers6))
        plot(rivers6a)
        accept<-0
        while(!any(accept==c("y","Y","n","N"))) accept <- readline(prompt="Accept changes? (y/n) ")  
        if(accept=="Y" | accept=="y") rivers6 <- rivers6a
        if(accept!="Y" & accept!="y") plot(rivers6)
        if(is.na(rivers5$mouth$mouth.seg)) {
          mouthseg <- readline(prompt="Please identify new segment number of river mouth ")
          showends(seg=as.numeric(mouthseg),rivers=rivers6)
          mouthvert <- readline(prompt="Please identify new vertex number of river mouth ")
          rivers6 <- setmouth(seg=as.numeric(mouthseg),vert=as.numeric(mouthvert),rivers=rivers6)
        }
      }
      checkagain<-0
      while(!any(checkagain==c("y","Y","n","N"))) checkagain<-readline(prompt="Re-check for braiding? (y/n) ")
      if(checkagain=="Y" | checkagain=="y") {
        rivers6 <- checkbraidedTF(rivers6)
        braided <- rivers6$braided
        if(braided) cat("Braiding still detected.",'\n')
        if(!braided) cat("Braiding no longer detected.",'\n')
      }
    }
  }
  if(!braided) {
    cat("No braiding detected within river network.",'\n','\n')
    yes<-0
    while(!any(yes==c("y","Y","n","N"))) yes<-readline(prompt="Build segment routes?  This will save time in route calculations. (y/n) ")
    if(yes=="Y" | yes=="y") rivers6 <- buildsegroutes(rivers6)
  } 
  
  cat('\n',"Cleanup completed, returning a network with",length(rivers6$lines),"segments.",'\n')
  cat('\n',"Recommend saving output to a .Rdata or .rda file.",'\n')
  cat('\n',"Note: any point data already using the input river network must be re-transformed to river coordinates using xy2segvert() or ptshp2segvert().",'\n')
  return(rivers6)
}




#' Connect Segments
#' @description Provides a method to manually connect unconnected segments
#'   within a river network.  The nearest endpoint (or vertex) of the second segment is
#'   added as a new vertex to the first, and the network topology is then updated.
#' @param connect The segment to connect to the network.  Typically, this is the
#'   segment that is disconnected from the rest of the river network.
#' @param connectto The segment to connect it to.  Typically, this is a segment
#'   that is connected to the rest of the river network.
#' @param nearestvert Whether to connect at the nearest vertex and split the
#'   segment (\code{FALSE}), or connect at the nearest endpoint (\code{TRUE}). 
#'   Defaults to \code{TRUE}.
#' @param rivers The river network object to use.
#' @return A new river network object with the specified segments connected (see
#'   \link{rivernetwork})
#' @seealso line2network
#' @note This function is called within \link{cleanup}, which is recommended in
#'   most cases.
#' @author Matt Tyers
#' @examples
#' data(Koyukuk0)
#' plot(Koyukuk0, ylim=c(1930500,1931500), xlim=c(194900,195100))
#' topologydots(Koyukuk0, add=TRUE)
#' 
#' Koyukuk0.1 <- connectsegs(connect=21, connectto=20, rivers=Koyukuk0)
#' plot(Koyukuk0.1,ylim=c(1930500,1931500), xlim=c(194900,195100))
#' topologydots(Koyukuk0.1, add=TRUE)
#' @export
connectsegs <- function(connect,connectto,nearestvert=F,rivers) {
  if(length(whoconnected(connect,rivers))>0){
    if(any(whoconnected(connect,rivers)==connectto)) stop("Segments are already connected.")
  }
  l1 <- dim(rivers$lines[[connect]])[1]
  l2 <- dim(rivers$lines[[connectto]])[1]
  if(!nearestvert) {
    dists <- pdist(rivers$lines[[connect]][1,],rivers$lines[[connectto]][1,])
    dists[2] <- pdist(rivers$lines[[connect]][1,],rivers$lines[[connectto]][l2,])
    dists[3] <- pdist(rivers$lines[[connect]][l1,],rivers$lines[[connectto]][1,])
    dists[4] <- pdist(rivers$lines[[connect]][l1,],rivers$lines[[connectto]][l2,])
    if(min(dists)==dists[1]) {
      rivers$lines[[connect]] <- rbind(rivers$lines[[connectto]][1,],rivers$lines[[connect]])
    }
    if(min(dists)==dists[2]) {
      rivers$lines[[connect]] <- rbind(rivers$lines[[connectto]][l2,],rivers$lines[[connect]])
    }
    if(min(dists)==dists[3]) {
      rivers$lines[[connect]] <- rbind(rivers$lines[[connect]],rivers$lines[[connectto]][1,])
    }
    if(min(dists)==dists[4]) {
      rivers$lines[[connect]] <- rbind(rivers$lines[[connect]],rivers$lines[[connectto]][l2,])
    }
    
    rivers$sp@lines[[rivers$lineID[connect,2]]]@Lines[[rivers$lineID[connect,3]]]@coords <- rivers$lines[[connect]]
    rivers$lengths[[connect]] <- rivers$lengths[[connect]]+min(dists)
    
    # updating the connectivity matrix 
    length <- length(rivers$lines)
    for(i in 1:length) {
      for(j in 1:length) {
        i.max <- dim(rivers$lines[[i]])[1]
        j.max <- dim(rivers$lines[[j]])[1]
        if(pdist(rivers$lines[[i]][1,],rivers$lines[[j]][1,])<rivers$tolerance & i!=j) {
          rivers$connections[i,j] <- 1
        }
        if(pdist(rivers$lines[[i]][1,],rivers$lines[[j]][j.max,])<rivers$tolerance & i!=j) {
          rivers$connections[i,j] <- 2
        }
        if(pdist(rivers$lines[[i]][i.max,],rivers$lines[[j]][1,])<rivers$tolerance & i!=j) {
          rivers$connections[i,j] <- 3
        }
        if(pdist(rivers$lines[[i]][i.max,],rivers$lines[[j]][j.max,])<rivers$tolerance & i!=j) {
          rivers$connections[i,j] <- 4
        }
        if(pdist(rivers$lines[[i]][1,],rivers$lines[[j]][1,])<rivers$tolerance & pdist(rivers$lines[[i]][i.max,],rivers$lines[[j]][j.max,])<rivers$tolerance & i!=j) {     ##########
          rivers$connections[i,j] <- 5
        }
        if(pdist(rivers$lines[[i]][i.max,],rivers$lines[[j]][1,])<rivers$tolerance & pdist(rivers$lines[[i]][1,],rivers$lines[[j]][j.max,])<rivers$tolerance & i!=j) {
          rivers$connections[i,j] <- 6
        }    ##########
      }
    }
    if(any(rivers$connections %in% 5:6)) rivers$braided <- TRUE
  }
  
  if(nearestvert) {
    dbeg <- dend <- 100000*max(rivers$lengths)
    whichbeg <- whichend <- NA
    for(i in 1:l2) {
      if(pdist(rivers$lines[[connect]][1,],rivers$lines[[connectto]][i,])<=dbeg) {
        dbeg <- pdist(rivers$lines[[connect]][1,],rivers$lines[[connectto]][i,])
        whichbeg <- i
      }
      if(pdist(rivers$lines[[connect]][l1,],rivers$lines[[connectto]][i,])<=dend) {
        dend <- pdist(rivers$lines[[connect]][l1,],rivers$lines[[connectto]][i,])
        whichend <- i
      }
    }
    if(dbeg < dend) {
      rivers$lines[[connect]] <- rbind(rivers$lines[[connectto]][whichbeg,],rivers$lines[[connect]])
    }
    if(dbeg > dend) {
      rivers$lines[[connect]] <- rbind(rivers$lines[[connect]],rivers$lines[[connectto]][whichend,])
    }
    
    rivers$sp@lines[[rivers$lineID[connect,2]]]@Lines[[rivers$lineID[connect,3]]]@coords <- rivers$lines[[connect]]
    rivers$lengths[[connect]] <- rivers$lengths[[connect]]+min(c(dbeg,dend))
    
    # if(!is.null(rivers$segroutes)) {
    #   rivers$segroutes <- NULL
    #   warning("Segment routes must be rebuilt - see help(buildsegroutes).")
    # }
    rivers <- splitsegments(rivers)
  }
  
  if(!is.null(rivers$segroutes)) {
    # rivers$segroutes <- NULL
    # warning("Segment routes must be rebuilt - see help(buildsegroutes).")
    rivers <- buildsegroutes(rivers,lookup=F)
  }
  rivers <- addcumuldist(rivers)
  if(!is.null(rivers$distlookup)) rivers <- buildlookup(rivers)
  
  message("Note: any point data already using the input river network must be re-transformed to river coordinates using xy2segvert() or ptshp2segvert().")
  return(rivers)
}



#' Remove Segments that are Smaller than the Connectivity Tolerance
#' @description Automatically detects and removes segments with total displacement (straight-line distance between endpoints) less than the connectivity tolerance.  These segments do not serve any real purpose, are bypassed in routes, and cannot be dissolved.
#' @param rivers The river network object to use.
#' @return A new river network object with the specified segments connected (see \link{rivernetwork})
#' @seealso line2network
#' @note This function is called within \link{cleanup}, which is recommended in most cases.
#' @author Matt Tyers
#' @examples
#' data(abstreams0)
#' abstreams1 <- removemicrosegs(abstreams0)
#' @export
removemicrosegs <- function(rivers) {
  displacement <- NA
  for(i in 1:length(rivers$lines)) {
    dim.i <- dim(rivers$lines[[i]])[1]
    displacement[i] <- pdist(rivers$lines[[i]][1,],rivers$lines[[i]][dim.i,])
  }
  problems <- (1:length(rivers$lines))[displacement<=rivers$tolerance]
  for(j in problems) {
    connectedto <- whoconnected(seg=j,rivers=rivers)
    for(jj in connectedto) {
      for(jjj in connectedto) {
        if(jj!=jjj & !any(whoconnected(seg=jj,rivers=rivers)==jjj)) {
          rivers <- connectsegs(connect=jj,connectto=jjj,rivers=rivers)
        }
      }
    }
  }
  if(length(problems)>0) rivers <- trimriver(trim=problems,rivers=rivers)
  # message("Note: any point data already using the input river network must be re-transformed to river coordinates using xy2segvert() or ptshp2segvert().")
  return(rivers)
}


#' Add Vertices To Maintain a Minumum Distance Between Vertices
#' @description In certain cases, such as when there is a lake within a river system, there may be long, straight lines in a river network with vertices only at either end.
#' In these cases, any point data along these stretches would be snapped to the vertices at either end.  This function automatically
#' adds equally-spaced vertices along the straight line, according to a specified minimum allowable distance between vertices. 
#' @param rivers The river network object to use.
#' @param mindist The minimum distance to use between vertices.  Defaults to 500.
#' @return A new river network object with the specified segments connected (see \link{rivernetwork})
#' @seealso line2network
#' @note This function is called within \link{cleanup}, which is recommended in most cases.
#' @author Matt Tyers
#' @examples
#' data(Kenai3)
#' Kenai3split <- addverts(Kenai3,mindist=200)
#' 
#' zoomtoseg(seg=c(73,70,71), rivers=Kenai3)
#' points(Kenai3$lines[[71]])        # vertices before adding
#' 
#' zoomtoseg(seg=c(73,70,71), rivers=Kenai3split)
#' points(Kenai3split$lines[[71]])   # vertices after adding
#' @export
addverts <- function(rivers,mindist=500) {

  lines <- rivers$lines 
  rivers1 <- rivers
  for(segi in 1:length(lines)) {
    seginew <- NULL 
    if(dim(lines[[segi]])[1] > 1) {
      for(verti in 2:(dim(lines[[segi]])[1])) {
        if(pdist(lines[[segi]][verti,],lines[[segi]][(verti-1),]) > mindist) {
          mindistx <- mindist*(lines[[segi]][verti,1]-lines[[segi]][(verti-1),1])/pdist(lines[[segi]][verti,],lines[[segi]][(verti-1),])
          mindisty <- mindist*(lines[[segi]][verti,2]-lines[[segi]][(verti-1),2])/pdist(lines[[segi]][verti,],lines[[segi]][(verti-1),])
          if(mindistx!=0) xtoadd <- seq(from=lines[[segi]][(verti-1),1],to=lines[[segi]][(verti),1],by=mindistx)
          if(mindisty!=0) ytoadd <- seq(from=lines[[segi]][(verti-1),2],to=lines[[segi]][(verti),2],by=mindisty)
          if(mindistx==0 & mindisty!=0) xtoadd <- rep(lines[[segi]][(verti),1],length(ytoadd))
          if(mindisty==0 & mindistx!=0) ytoadd <- rep(lines[[segi]][(verti),2],length(xtoadd))
        }
        if(pdist(lines[[segi]][verti,],lines[[segi]][(verti-1),]) <= mindist) {
          xtoadd <- lines[[segi]][(verti-1),1]
          ytoadd <- lines[[segi]][(verti-1),2]
        }
        seginew <- rbind(seginew,cbind(xtoadd,ytoadd))
      }
      seginew <- rbind(seginew,lines[[segi]][dim(lines[[segi]])[1],])
      lines[[segi]] <- unname(seginew)
    }
    
    # updating the sp object!
    rivers1$sp@lines[[rivers$lineID[segi,2]]]@Lines[[rivers$lineID[segi,3]]]@coords <- unname(seginew)
  }
   
  rivers1$lines <- lines
  
  # mouth 
  if(!is.na(rivers$mouth$mouth.seg) & !is.na(rivers$mouth$mouth.vert)) {
    mouthcoords <- rivers$lines[[rivers$mouth$mouth.seg]][rivers$mouth$mouth.vert,]
    for(j in 1:(dim(rivers1$lines[[rivers$mouth$mouth.seg]])[1])) {
      if(all(mouthcoords==rivers1$lines[[rivers$mouth$mouth.seg]][j,])) {
        rivers1$mouth$mouth.vert <- j
      }
    }
  }
  
  if(!is.null(rivers1$segroutes)) {
    # rivers1$segroutes <- NULL
    # warning("Segment routes must be rebuilt - see help(buildsegroutes).")
    rivers1 <- buildsegroutes(rivers1,lookup=F)
  }
  rivers1 <- addcumuldist(rivers1)
  if(!is.null(rivers1$distlookup)) rivers1 <- buildlookup(rivers1)
  
  message("Note: any point data already using the input river network must be re-transformed to river coordinates using xy2segvert() or ptshp2segvert().")
  return(rivers1)
}