#' Detect Route
#' @description Called internally within \link{riverdistance}.  Detects the 
#'   sequential route from one river network segment to another.
#' @param start Segment number of the start of the route
#' @param end Segment number of the end of the route
#' @param rivers The river network object to use
#' @param verbose Whether or not to print all routes being considered (used for 
#'   error checking).  Defaults to FALSE.
#' @param stopiferror Whether or not to exit with an error if a route cannot be
#'   found.  If this is set to \code{FALSE} and a route cannot be found,
#'   \code{detectroute()} will return \code{NA}.  Defaults to \code{TRUE}.
#' @param algorithm Which route detection algorithm to use.  If set to \code{NULL} (the default), the function will automatically make a selection.  Choices are:
#' \itemize{
#' \item Setting \code{algorithm="sequential"} will be quite slow, and may give inaccurate results in the event of braiding.  This algorithm returns the first complete route detected, which may not be the shortest.  This algorithm is not recommended in almost all cases, but is retained as an option for certain checks.  It will not be used unless specified.
#' \item Setting \code{algorithm="Dijkstra"} will be much faster, and will return the shortest route in the event of braiding.  If braiding is present or unknown, this will be the algorithm automatically chosen.
#' \item Setting \code{algorithm="segroutes"} will be the fastest of all, but will only return results in a non-braided network.  This will be the algorithm automatically selected if segment routes are present - see \link{buildsegroutes}.}
#' @return A vector of segment numbers corresponding to the ordered route.
#' @author Matt Tyers
#' @examples
#' data(Gulk)
#' plot(x=Gulk, cex=1)
#' 
#' detectroute(start=6, end=14, rivers=Gulk)
#' 
#' tstart <- Sys.time()
#' detectroute(start=120, end=111, rivers=abstreams, algorithm="sequential")
#' tend <- Sys.time()
#' tend - tstart
#' 
#' data(abstreams)
#' tstart <- Sys.time()
#' detectroute(start=120, end=111, rivers=abstreams, algorithm="Dijkstra")
#' tend <- Sys.time()
#' tend - tstart
#' 
#' tstart <- Sys.time()
#' detectroute(start=120, end=111, rivers=abstreams, algorithm="segroutes")
#' tend <- Sys.time()
#' tend - tstart
#' @export
detectroute <- function(start,end,rivers,verbose=FALSE,stopiferror=TRUE,algorithm=NULL) {
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  
  if(start==end) return(start)
  
  if(is.null(algorithm)) {
    algorithm <- ifelse(is.null(rivers$segroutes), "Dijkstra", "segroutes")
  }
  if(algorithm=="segroutes" & is.null(rivers$segroutes)) { 
    algorithm <- "Dijkstra"
  }
  
  if(!any(algorithm==c("sequential","Dijkstra","segroutes"))) stop("Invalid algorithm specified.")
  if(verbose) cat("Using",algorithm,"algorithm...",'\n')

  if(max(c(start,end),na.rm=T)>length(rivers$lines) | min(c(start,end),na.rm=T)<1) {
    stop("Invalid segments specified.")
  }
  
  if(algorithm=="segroutes") {
    segroutes <- rivers$segroutes
    connections <- rivers$connections
    if(any(segroutes[[start]]==end)) {
      route <- segroutes[[start]][which(segroutes[[start]]==start):which(segroutes[[start]]==end)]
    }
    if(any(segroutes[[end]]==start)) {
      route <- segroutes[[end]][which(segroutes[[end]]==start):which(segroutes[[end]]==end)]
    }
    if(all(segroutes[[end]]!=start) & all(segroutes[[start]]!=end)) {
      
      minl <- min(length(segroutes[[start]]), length(segroutes[[end]]))
      ijunct <- min(which(segroutes[[start]][1:minl] != segroutes[[end]][1:minl]))
      
      route <- c(segroutes[[start]][which(segroutes[[start]]==start):ijunct], segroutes[[end]][ijunct:which(segroutes[[end]]==end)])
      if(is.na(connections[segroutes[[start]][ijunct],segroutes[[end]][ijunct]])) {
        route <- c(segroutes[[start]][which(segroutes[[start]]==start):ijunct], segroutes[[start]][1], segroutes[[end]][ijunct:which(segroutes[[end]]==end)])
      }
    }
    return(route)
  }
  
  if(algorithm=="Dijkstra") {
    connections <- rivers$connections
    length <- length(rivers$lines)
    lengths <- rivers$lengths
    max1 <- 2*sum(lengths)
    
    segs <- 1:length
    dists <- rep(max1,length)
    visited <- rep(FALSE,length)
    minroutes <- list()
    minroutes_vec <- rep(NA,length)  
    
    connected <- !is.na(connections)
    dists[start] <- lengths[start]
    current <- start
    minroutes[[start]] <- start
    
    found <- FALSE
    while(!found) {
      neighbors <- (1:length)[connected[current,]&!visited] 
      if(length(neighbors)>0) {
        dists.tentative <- sum(lengths[minroutes[[current]]]) + lengths[neighbors]  
        check1 <- dists.tentative <= dists[neighbors]  
        if(any(check1)) {
          neighborscheck1 <- neighbors[check1]
          dists[neighborscheck1] <- dists.tentative[check1] 
          for(neighbor in neighborscheck1) {
            minroutes[[neighbor]] <- c(minroutes[[current]],neighbor)
          } 
          minroutes_vec[neighborscheck1] <- neighborscheck1 
        } 
      }
      
      if(is.null(minroutes[[current]][1])) {  
        if(stopiferror) stop("No route detected.")
        found <- TRUE
        minroutes[[end]] <- NA
      }
      
      if(verbose) print(minroutes[[current]])
      visited[current] <- TRUE
      if(current==end) found <- TRUE
      
      if(!found) {
        thing1 <- which.min(dists+(visited*max1))
        if(current == thing1) { 
          dists[current] <- max1
          connected[,current] <- connected[current,] <- F
          current <- minroutes[[current]][length(minroutes[[current]])-1]  
        }
        else current <- thing1 
      }
      
      if(!any(connected[current,minroutes_vec[!is.na(minroutes_vec)]]) & !(current %in% minroutes_vec)) {  
        if(stopiferror) stop("No route detected.")
        found <- TRUE
        minroutes[[end]] <- NA
      }   
      if(all(visited)&!found) {
        if(stopiferror) stop("No route detected.")
        found <- TRUE
        minroutes[[end]] <- NA
      }
    }
    
    return(minroutes[[end]])
  }

  if(algorithm=="sequential") {
    connections <- rivers$connections
    length <- length(rivers$lines)
    if(start==end) path <- start
    if(start != end) {
      is.connected <- is.na(connections)==F
      path <- start             # vector of the current path being considered
      option <- 1            # vector of the option number (of available at each level) the current path is
      used <- rep(F,length)  # which segment numbers have been used
      
      used <- matrix(F,nrow=length,ncol=length)
      used[1,start]<-T
      
      level<-1
      found<-F
      while(found==F) {
        connected.to <- (1:length)[is.connected[,path[level]]]  # which are connected to the current terminal node
        connected.to.new <- connected.to
        for(i in 1:length(connected.to)) {  # checking if target has been reached, defining new possible connections        
          if(length(connected.to)==0) {
            if(stopiferror) {
              cat('\n',"Unable to calculate route between segments",start,"and",end,'\n')
              stop("Gaps exist between specified segments.") 
            }
            if(!stopiferror) {
              path <- NA
              found <- T
            }
          }
          if(!found) {
            if(connected.to[i]==end) {
              found <- T
              path[level+1]<-end
            }
            if(used[level,][connected.to[i]]==T) connected.to.new <- connected.to.new[connected.to.new!=connected.to[i]]
          }
        }
        if(!found & length(connected.to.new)>0) {  # if dead end has not been reached
          used[level,][connected.to.new] <- T
          level <- level+1
          if(found==F) path[level]<-connected.to.new[1]
          option[level]<-1
          used[level,]<-used[level-1,]  
        }
        if(!found & length(connected.to.new)==0) { #if dead end has been reached
          #cat("oops dead end ... ")
          good.to.go <- F
          while(good.to.go==F) {
            option[level] <- option[level]+1
            connected.to.prev <- (1:length)[is.connected[,path[level-1]]]  # which are connected to the previous terminal node
            connected.to.new.prev <- connected.to.prev
            for(i in 1:length(connected.to.prev)) {  # seeing if a different option exists for the previous terminal node
              if(level>2) {
                if(used[level-2,][connected.to.prev[i]]==T) connected.to.new.prev <- connected.to.new.prev[connected.to.new.prev!=connected.to.prev[i]]
              }
            } 
            if(length(connected.to.new.prev)<option[level]) {  # if not, back out another level
              level <- level-1
              path <- path[1:level]
              option <- option[1:level]
            }
            if(length(option[level])==0) {
              if(stopiferror) {
                cat('\n',"Unable to calculate route between segments",start,"and",end,'\n')
                stop("Gaps exist between specified segments.") 
              }
              if(!stopiferror) {
                path <- NA
                found <- T
              } 
            }
            if(level>0) {   
              if(length(connected.to.new.prev)>=option[level]) {  # if so, use the different option
                path[level] <- connected.to.new.prev[option[level]]
                good.to.go <-T
              }
            }
            if(level<1) { 
              if(stopiferror) {
                cat('\n',"Unable to calculate route between segments",start,"and",end,'\n')
                stop("Gaps exist between specified segments.") 
              }
              if(!stopiferror) {
                path <- NA
                found <- T
                good.to.go <- T
              }
            }
          }
        }
        if(verbose) cat(path,'\n')
      }
    }
    return(path)
  }
}


#' Build Segment Routes
#' @description Adds the travel routes from the mouth (lowest point) of a river 
#'   network to each segment, and (optionally) distance lookup tables.  This
#'   greatly reduces the time needed to detect routes, making distance
#'   calculation much more efficient, particularly in the case of multiple
#'   distance calculations.
#' @param rivers The river network object to use
#' @param lookup Whether to build lookup tables as well.  This may take
#'   some time, but will result in even faster distance computation in analyses
#'   (see \link{buildlookup}).  Because of the object size returned, this may
#'   not be advisable in a large river network (more than a few hundred
#'   segments).  Accepts \code{TRUE} or \code{FALSE}, and defaults to
#'   \code{NULL}.  If the default value is accepted, lookup tables will be built
#'   if the river network has 400 segments or fewer.
#' @param verbose Whether or not to print the segment number the function is 
#'   currently building a route for (used for error checking).  Defaults to 
#'   FALSE.
#' @return A rivernetwork object, with a new list element, \code{$segroutes}, 
#'   which gives the route from the mouth to each rivernetwork segment. 
#'   Optionally, it may add \code{$distlookup}, distance lookup tables for even
#'   faster distance computation. (See \link{rivernetwork}.)
#' @note In the event of braiding (multiple channels), it is likely that there 
#'   will be differences in the routes detected.  If this is the case, building 
#'   routes will likely result in a shorter and more efficient route. 
#'   Regardless, extreme caution is always advised in the event of braiding.
#' @note The mouth segment and vertex must be specified (see \link{setmouth}).
#' @author Matt Tyers
#' @note This function is called within \link{cleanup}, which is recommended in
#'   most cases.
#' @examples
#' data(abstreams)
#' plot(x=abstreams)
#' abstreams1 <- abstreams
#' abstreams1$segroutes <- NULL #taking out the $segroutes component
#' 
#' # before
#' tstart <- Sys.time()
#' detectroute(start=120, end=111, rivers=abstreams1)
#' Sys.time() - tstart
#' 
#' # after
#' tstart <- Sys.time()
#' detectroute(start=120, end=111, rivers=abstreams)
#' Sys.time() - tstart
#' @importFrom graphics plot
#' @export
buildsegroutes1 <- function(rivers,lookup=NULL,verbose=FALSE) {
  if(is.na(rivers$mouth$mouth.seg) | is.na(rivers$mouth$mouth.vert)) stop("need to supply the segment and vertex of origin")
  
  if(length(rivers$lines)==1) {
    rivers$segroutes <- list(1)
    if(is.null(lookup)) lookup <- T
    if(lookup) rivers <- buildlookup(rivers)
    return(rivers)
  } else{
  
  if(is.na(rivers$braided)) rivers <- checkbraidedTF(rivers=rivers)
  if(rivers$braided) stop("Building segment routes is inappropriate in a braided river network.")
  
  routes <- list()
  for(i in 1:length(rivers$lines)) {
    routes[[i]] <- NA
  }
  #internal functions
  n.top <- function(seg,connections) {
    return(length(connections[seg,][(connections[seg,]==1 | connections[seg,]==2 | connections[seg,]==5 | connections[seg,]==6) & is.na(connections[seg,])==F]))   
  }
  n.bot <- function(seg,connections) {
    return(length(connections[seg,][(connections[seg,]==3 | connections[seg,]==4 | connections[seg,]==5 | connections[seg,]==6) & is.na(connections[seg,])==F]))
  }
  dists <- rep(NA,length(rivers$lines))
  for(i in 1:length(rivers$lines)) {
    dists[i] <- pdist(rivers$lines[[rivers$mouth$mouth.seg]][rivers$mouth$mouth.vert,],
                      c(mean(rivers$lines[[i]][,1]),mean(rivers$lines[[i]][,2])))
  }
  order <- order(dists)
  rivers1 <- rivers
  rivers1$connections <- rivers$connections[order, order]
  rivers1$lines <- rivers$lines[order]
  rivers1$lengths <- rivers$lengths[order]
  if(verbose) plot(x=rivers1,main="segments used for building routes")
  if(verbose) print("building routes for end segments...")
  for(i in 1:length(rivers$lines)) {
    if((n.top(i,rivers1$connections)==0 & n.bot(i,rivers1$connections)!=0) | 
       (n.top(i,rivers1$connections)!=0 & n.bot(i,rivers1$connections)==0)) {
      if(verbose) print(i)
      routes1a <- detectroute(start=i,end=order(order)[rivers$mouth$mouth.seg],rivers=rivers1,algorithm="Dijkstra")  ###############
      routes1 <- routes1a[length(routes1a):1]
      for(j in 1:length(routes1)) {
        routes[[routes1[j]]] <- routes1[1:j]
      }
    }
  }
  if(verbose) print("building routes for middle segments...")
  for(i in 1:length(rivers$lines)) {
    if(is.na(routes[[i]][1])) {
      if(verbose) print(i)
      routes_ia <- detectroute(start=i,end=order(order)[rivers$mouth$mouth.seg],rivers=rivers1,algorithm="Dijkstra")  ###############
      routes[[i]] <- routes_ia[length(routes_ia):1]
    }
  }
  realroutes <- list()
  for(i in 1:length(rivers$lines)) {
    realroutes[[order[i]]] <- order[routes[[i]]]
  }
  
  rivers$segroutes <- realroutes
  
  if(is.null(lookup) & length(rivers$lines)<=400) lookup <- T
  if(is.null(lookup) & length(rivers$lines)>400) lookup <- F
  if(lookup) rivers <- buildlookup(rivers) 

  return(rivers) }
}


#' River Distance
#' @description Calculates the total river network distance between two points 
#'   on the river network, given in river locations (segment and vertex).
#' @param startseg Segment number of the start of the route
#' @param endseg Segment number of the end of the route
#' @param startvert Vertex number of the start of the route
#' @param endvert Vertex number of the end of the route
#' @param rivers The river network object to use
#' @param path (optional) The vector-format route of segment numbers can also be
#'   supplied instead of the starting and ending segments.
#' @param map Whether or not to draw a sanity-check map, showing the calculated 
#'   route in entirety.  Defaults to FALSE.
#' @param add If \code{map==TRUE}, whether to add the route drawing to an
#'   existing plot (\code{add=TRUE}) or produce a new plot (\code{add=FALSE}).
#' @param stopiferror Whether or not to exit with an error if a route cannot be 
#'   found.  If this is set to \code{FALSE} and a route cannot be found, 
#'   \code{riverdistance()} will return \code{NA}.  Defaults to \code{TRUE}. 
#'   See \link{detectroute}.
#' @param algorithm Which route detection algorithm to use (\code{"Dijkstra"},
#'   \code{"sequential"}, or \code{"segroutes"}).  If left as \code{NULL} (the
#'   default), the function will automatically make a selection.  See
#'   \link{detectroute} for more details.  
#' @return Total route distance, in the units of the coordinate system used 
#'   (this will likely be meters).
#' @note If a distance lookup table (\code{$distlookup}) is present in the river network object, accepting \code{NULL} will bypass route detection and return distance automatically, the fastest algorithm of all.  This is done automatically in \link{buildsegroutes}, but can be called directly using \link{buildlookup}.
#' @author Matt Tyers
#' @note Building routes from the river mouth to each river network segment and/or distance lookup tables will
#'   greatly reduce computation time (see \link{buildsegroutes}).
#' @examples
#' data(Gulk)
#' riverdistance(startseg=6, endseg=14, startvert=100, endvert=200, rivers=Gulk)
#' riverdistance(startvert=100, endvert=200, path=c(6,3,4,10,11,14), rivers=Gulk)
#' riverdistance(startseg=6, endseg=14, startvert=100, endvert=200, rivers=Gulk, map=TRUE)
#'       
#' # speed comparison: 
#' 
#' data(abstreams)
#' 
#' tstart <- Sys.time()
#' riverdistance(startseg=120, startvert=10, endseg=131, endvert=10, rivers=abstreams, 
#'               algorithm="sequential")
#' Sys.time()- tstart
#' 
#' tstart <- Sys.time()
#' riverdistance(startseg=120, startvert=10, endseg=131, endvert=10, rivers=abstreams, 
#'               algorithm="Dijkstra")
#' Sys.time()- tstart
#' 
#' tstart <- Sys.time()
#' riverdistance(startseg=120, startvert=10, endseg=131, endvert=10, rivers=abstreams)
#' 
#' # Note: it is not necessary to specify the algorithm here: the distance function
#' # will automatically select the fastest algorithm unless otherwise specified.
#' Sys.time()- tstart
#' 
#' @export
riverdistance <- function(startseg=NULL,endseg=NULL,startvert,endvert,rivers,path=NULL,map=FALSE,add=FALSE,stopiferror=TRUE,algorithm=NULL) {
  if(is.null(rivers$cumuldist)) {
    rivers <- addcumuldist(rivers)
    warning("River network does not have cumulative distances - recommend adding them with addcumuldist()")
  }
  if(!is.null(rivers$distlookup) & !map & is.null(algorithm)) {
    cumuldist <- rivers$cumuldist
    lengths <- rivers$lengths
    x<-rivers$distlookup
    if(!is.null(path)) {
      startseg <- path[1]
      endseg <- path[length(path)]
    }
    if(startseg<1 | endseg<1 | startseg>length(lengths) | endseg>length(lengths)) stop("Invalid segment specified")
    if(startvert<1 | endvert<1 | startvert>length(cumuldist[[startseg]]) | endvert>length(cumuldist[[endseg]])) stop("Invalid vertex specified")
    if(startseg==endseg) dist <- abs(cumuldist[[startseg]][startvert] - cumuldist[[startseg]][endvert])
    else {
      if(!is.na(x$starttop[startseg,endseg])) { 
        if(x$starttop[startseg,endseg] & x$endtop[startseg,endseg]) dist <- x$middist[startseg,endseg] + cumuldist[[startseg]][startvert] + cumuldist[[endseg]][endvert]
        if(x$starttop[startseg,endseg] & !x$endtop[startseg,endseg]) dist <- x$middist[startseg,endseg] + cumuldist[[startseg]][startvert] + lengths[endseg] - cumuldist[[endseg]][endvert]
        if(!x$starttop[startseg,endseg] & x$endtop[startseg,endseg]) dist <- x$middist[startseg,endseg] + lengths[startseg] - cumuldist[[startseg]][startvert] + cumuldist[[endseg]][endvert]
        if(!x$starttop[startseg,endseg] & !x$endtop[startseg,endseg]) dist <- x$middist[startseg,endseg] + lengths[startseg] - cumuldist[[startseg]][startvert] + lengths[endseg] - cumuldist[[endseg]][endvert]
      }
    }
    if(!is.na(rivers$braided)) if(rivers$braided) {    
      if(!is.na(rivers$connections[startseg,endseg])) {
        if(rivers$connections[startseg,endseg]==5) {
          d1 <- cumuldist[[startseg]][startvert] + cumuldist[[endseg]][endvert]
          d2 <- (lengths[startseg] - cumuldist[[startseg]][startvert]) + (lengths[endseg] - cumuldist[[endseg]][endvert])
          dist <- min(d1,d2)
        }
        if(rivers$connections[startseg,endseg]==6) {
          d1 <- cumuldist[[startseg]][startvert] + (lengths[endseg] - cumuldist[[endseg]][endvert])
          d2 <- (lengths[startseg] - cumuldist[[startseg]][startvert]) + cumuldist[[endseg]][endvert]
          dist <- min(d1,d2)
        }
      }
    }
    return(dist)
  }
  
  connections <- rivers$connections
  seg.lengths <- rivers$lengths
  lines <- rivers$lines
  cumuldist <- rivers$cumuldist
  
  if(is.null(path)) path <- detectroute(start=startseg,end=endseg,rivers=rivers,stopiferror=stopiferror,algorithm=algorithm)
  
  if(is.na(path[1])) return(NA)
  startseg <- path[1]
  endseg <- path[length(path)]
  
  if(startvert>dim(lines[[startseg]])[1] | startvert<1 | endvert>dim(lines[[endseg]])[1] | endvert<1) {
    stop("Invalid vertex specified.")
  } 
  
  if(map) {
    if(!add) plot(x=rivers)
    riverpoints(seg=c(startseg,endseg),vert=c(startvert,endvert),rivers=rivers,pch=15,col=4)
    if(length(path)>1) {
      # distance on the partial segment the beginning is on
      if(any(connections[startseg,path[2]]==c(1,2))) {     
        if(startvert!=1) { 
          lines(lines[[startseg]][1:startvert,],lwd=3,col=4)
        }
      }
      if(any(connections[startseg,path[2]]==c(3,4))) {
        linelength <- dim(lines[[startseg]])[1]
        if(linelength!=startvert) {
          lines(lines[[startseg]][linelength:startvert,],lwd=3,col=4)
        }
      }
      # distance on the full segments between the beginning and end segments
      if(length(path)>2) {
        for(i in path[2:(length(path)-1)]) lines(lines[[i]],col=4,lty=1,lwd=3)
      }
      #distance on the partial segment the end is on
      if(any(connections[path[length(path)-1],endseg]==c(1,3))) {
        if(endvert!=1) { 
          lines(lines[[endseg]][1:endvert,],lwd=3,col=4)
        }
      }
      if(any(connections[path[length(path)-1],endseg]==c(2,4))) {
        linelength <- dim(lines[[endseg]])[1]
        if(linelength!=endvert) {
          lines(lines[[endseg]][linelength:endvert,],lwd=3,col=4)
        }
      }
      # special braided case
      if(length(path)==2 & connections[path[1],path[2]]==5) {
        d1 <- cumuldist[[startseg]][startvert] + cumuldist[[endseg]][endvert]
        d2 <- (seg.lengths[startseg] - cumuldist[[startseg]][startvert]) + (seg.lengths[endseg] - cumuldist[[endseg]][endvert])
        if(d1<d2) {
          lines(lines[[startseg]][1:startvert,],lwd=3,col=4)
          lines(lines[[endseg]][1:endvert,],lwd=3,col=4)
        } else{
          lines(lines[[startseg]][(dim(lines[[startseg]])[1]):startvert,],lwd=3,col=4)
          lines(lines[[endseg]][(dim(lines[[endseg]])[1]):endvert,],lwd=3,col=4)
        }
      }
      if(length(path)==2 & connections[path[1],path[2]]==6) {
        d1 <- cumuldist[[startseg]][startvert] + (seg.lengths[endseg] - cumuldist[[endseg]][endvert])
        d2 <- (seg.lengths[startseg] - cumuldist[[startseg]][startvert]) + cumuldist[[endseg]][endvert]
        if(d1<d2) {
          lines(lines[[startseg]][1:startvert,],lwd=3,col=4)
          lines(lines[[endseg]][(dim(lines[[endseg]])[1]):endvert,],lwd=3,col=4)
        } else{
          lines(lines[[startseg]][(dim(lines[[startseg]])[1]):startvert,],lwd=3,col=4)
          lines(lines[[endseg]][1:endvert,],lwd=3,col=4)
        }
      }  
    }
    # if the beginning and end are on the same segment
    if(length(path)==1) {
      min <- min(startvert,endvert)
      max <- max(startvert,endvert)
      if(min!=max) {
        lines(lines[[path[1]]][min:max,],lwd=3,col=4)
      }
    }
  }
  
  route.dist <- 0
  
  if(length(path)>1) {
    # distance on the partial segment the beginning is on
    if(any(connections[startseg,path[2]]==c(1,2))) {
      if(startvert!=1) { 
        # route.dist <- route.dist + pdisttot(lines[[startseg]][1:startvert,])
        route.dist <- route.dist + cumuldist[[startseg]][startvert]    
      }
    }
    if(any(connections[startseg,path[2]]==c(3,4))) {
      linelength <- dim(lines[[startseg]])[1]
      if(linelength!=startvert) {
        route.dist <- route.dist + seg.lengths[startseg] - cumuldist[[startseg]][startvert]
      }
    } 
    # distance on the full segments between the beginning and end segments
    if(length(path)>2) {
      route.dist <- route.dist+sum(seg.lengths[path[2:(length(path)-1)]])
    } 
    #distance on the partial segment the end is on
    if(any(connections[path[length(path)-1],endseg]==c(1,3))) {
      if(endvert!=1) { 
        route.dist <- route.dist + cumuldist[[endseg]][endvert]
      }
    }
    if(any(connections[path[length(path)-1],endseg]==c(2,4))) {
      linelength <- dim(lines[[endseg]])[1]
      if(linelength!=endvert) {
        route.dist <- route.dist + seg.lengths[endseg] - cumuldist[[endseg]][endvert]
      }
    } 
    # special braided case 
    if(length(path)==2 & connections[path[1],path[2]]==5) {
      d1 <- cumuldist[[startseg]][startvert] + cumuldist[[endseg]][endvert]
      d2 <- (seg.lengths[startseg] - cumuldist[[startseg]][startvert]) + (seg.lengths[endseg] - cumuldist[[endseg]][endvert])
      route.dist <- route.dist + min(d1,d2)
    }
    if(length(path)==2 & connections[path[1],path[2]]==6) {
      d1 <- cumuldist[[startseg]][startvert] + (seg.lengths[endseg] - cumuldist[[endseg]][endvert])
      d2 <- (seg.lengths[startseg] - cumuldist[[startseg]][startvert]) + cumuldist[[endseg]][endvert]
      route.dist <- route.dist + min(d1,d2)
    }  
  }
  # if the beginning and end are on the same segment
  if(length(path)==1) {
    min <- min(startvert,endvert)
    max <- max(startvert,endvert)
    if(startvert!=endvert) {
      route.dist <- route.dist + cumuldist[[startseg]][max] - cumuldist[[startseg]][min]
    }
  } 
  
  return(route.dist)
}


#' Build Lookup Tables for Fast Distance Computation
#' @description Adds lookup tables for distance computation, dramatically
#'   reducing computation time.  It may take some time to calculate,
#'   particularly in a braided network.
#' @param rivers The river network object to use
#' @return A rivernetwork object, with a new list element, \code{$distlookup}, a
#'   list of three matrices.  Element \code{[i,j]} of each matrix corresponds to
#'   the route between segment \code{i} and \code{j}.  The
#'   \code{distlookup$middist} matrix gives the total distance of the "middle"
#'   of each route (between the starting and ending segments"), and the
#'   \code{distlookup$starttop} and \code{distlookup$endtop} matrices have value
#'   \code{TRUE}, \code{FALSE}, or \code{NA} if the segments at the beginning or
#'   end of the route are connected to the resto of the route at the top of the
#'   coordinate matrix, bottom of the coordinate matrix, or if the route is
#'   contained to just one segment, respectively. (See \link{rivernetwork}.)
#' @note This will add three n by n matrices to the river network object, which
#'   will be very large if the river network has many segments.
#' @author Matt Tyers
#' @note This function is called within \link{cleanup}, which is recommended in
#'   most cases.  It is also called within \link{buildsegroutes}, and will add
#'   lookup tables by default if there are fewer than 400 segments in the river
#'   network.
#' @note This function can still be called in the presence of a braided network, but all resulting distances used in subsequent analyses will be the shortest route.
#' @note If segment routes (\code{$segroutes}) are not present, this function may take a very long time to run.
#' @examples
#' data(abstreams)
#' 
#' abstreams1 <- buildlookup(abstreams)
#' @export
buildlookup <- function(rivers) {
  if(is.null(rivers$segroutes) & interactive()) pb <- txtProgressBar(style=3)
  if(!is.na(rivers$braided)) {
    if(rivers$braided) message("Braiding detected in river network - all distances used in analysis will be shortest-path distances, not necessarily travel distances.")
  } else {
    message("Braiding may be present in river network.  If so, all distances used in analysis will be shortest-path distances, not necessarily travel distances.")
  }
  length <- length(rivers$lines)
  connections <- rivers$connections
  middist <- matrix(0,nrow=length,ncol=length)
  starttop <- endtop <- matrix(NA,nrow=length,ncol=length)
  for(i in 1:length) {
    for(j in 1:length) {
      theroute <- detectroute(start=i, end=j, rivers=rivers)
      routelength <- length(theroute)
      if(length(theroute) > 2) middist[i,j] <- sum(rivers$lengths[theroute[2:(routelength-1)]])
      if(length(theroute) > 1) {
        if(any(connections[i, theroute[2]] == c(1,2))) starttop[i,j] <- T
        if(any(connections[i, theroute[2]] == c(3,4))) starttop[i,j] <- F
        if(any(connections[theroute[routelength-1], j] == c(1,3))) endtop[i,j] <- T
        if(any(connections[theroute[routelength-1], j] == c(2,4))) endtop[i,j] <- F
      }
      if(interactive() & is.null(rivers$segroutes)) setTxtProgressBar(pb=pb, value=(i/length + j/length/length))
    }
  }
  distlookup <- list(middist=middist, starttop=starttop, endtop=endtop)
  rivers$distlookup <- distlookup
  
  
  return(rivers)
}
