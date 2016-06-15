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
#' \item Setting \code{algorithm="Dijkstra"} may be slow, but the shortest route will be detected in the event of braiding.  If braiding is present or unknown, this will be the algorithm automatically chosen.
#' \item Setting \code{algorithm="sequential"} will be much faster, but may give inaccurate results in the event of braiding.  This algorithm returns the first complete route detected, which may not be the shortest.  This is the algorithm automatically chosen if the network is known not to be braided, but lacks segment routes.
#' \item Setting \code{algorithm="segroutes"} will be the fastest of all, but will only return results in a non-braided network.  This will be the algorithm automatically selected if segment routes are present.}
#' @return A vector of segment numbers corresponding to the ordered route.
#' @author Matt Tyers
#' @examples
#' data(Gulk)
#' plot(x=Gulk, cex=1)
#' 
#' detectroute(start=6, end=14, rivers=Gulk)
#' 
#' data(abstreams)
#' tstart <- Sys.time()
#' detectroute(start=120, end=111, rivers=abstreams, algorithm="Dijkstra")
#' tend <- Sys.time()
#' tend - tstart
#' 
#' tstart <- Sys.time()
#' detectroute(start=120, end=111, rivers=abstreams, algorithm="sequential")
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
  # connections <- rivers$connections
  # length <- length(rivers$lines)
  
  if(start==end) return(start)
  
  if(is.null(algorithm)) {
    if(length(rivers$segroutes)>0) algorithm <- "segroutes"
    if(!is.na(rivers$braided)) {
      if(rivers$braided) algorithm <- "Dijkstra"
    }
    if(is.na(rivers$braided)) {
      algorithm <- "Dijkstra"
    }
    if(is.null(algorithm)) algorithm <- "sequential"
  }
  if(algorithm=="segroutes" & is.null(rivers$segroutes)) { #(length(rivers$segroutes)==0)) {
    if(!is.na(rivers$braided)) {
      algorithm <- ifelse(rivers$braided,"Dijkstra","sequential")
    }
    if(is.na(rivers$braided)) {
      algorithm <- "Dijkstra"
    }
  }
  
  if(!any(algorithm==c("sequential","Dijkstra","segroutes"))) stop("Invalid algorithm specified.")
  if(verbose) cat("Using",algorithm,"algorithm...",'\n')

  if(max(c(start,end),na.rm=T)>length(rivers$lines) | min(c(start,end),na.rm=T)<1) {
    stop("Invalid segments specified.")
  }
  
  if(algorithm=="segroutes") {
    segroutes <- rivers$segroutes
    if(any(segroutes[[start]]==end)) {
      route <- segroutes[[start]][which(segroutes[[start]]==start):which(segroutes[[start]]==end)]
    }
    if(any(segroutes[[end]]==start)) {
      route <- segroutes[[end]][which(segroutes[[end]]==start):which(segroutes[[end]]==end)]
    }
    if(all(segroutes[[end]]!=start) & all(segroutes[[start]]!=end)) {
      # junct.end <- NA
      # k <- 1
      # while(is.na(junct.end)) {
      #   if(segroutes[[start]][k] != segroutes[[end]][k]) {
      #     junct.start <- segroutes[[start]][k]
      #     junct.end <- segroutes[[end]][k]
      #   }
      #   k <- k+1
      # }
      # part1 <- segroutes[[start]][which(segroutes[[start]]==start):which(segroutes[[start]]==junct.start)]
      # part2 <- segroutes[[end]][which(segroutes[[end]]==junct.end):which(segroutes[[end]]==end)]
      
      # ijunct <- suppressWarnings(min(which(segroutes[[start]] != segroutes[[end]])))
      
      minl <- min(length(segroutes[[start]]), length(segroutes[[end]]))
      ijunct <- min(which(segroutes[[start]][1:minl] != segroutes[[end]][1:minl]))
      
      # part1 <- segroutes[[start]][which(segroutes[[start]]==start):ijunct]
      # part2 <- segroutes[[end]][ijunct:which(segroutes[[end]]==end)]
      # route <- c(part1,part2)
      route <- c(segroutes[[start]][which(segroutes[[start]]==start):ijunct], segroutes[[end]][ijunct:which(segroutes[[end]]==end)])
    }
    return(route)
  }
  
  if(algorithm=="Dijkstra") {
    tolerance <- rivers$tolerance
    lines <- rivers$lines
    connections <- rivers$connections
    length <- length(lines)
    lengths <- rivers$lengths
    max <- 2*sum(lengths)
    
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
    
    segs <- 1:length
    dists <- rep(max,length)
    visited <- rep(FALSE,length)
    minroutes <- list()
    
    connected <- !is.na(connections)
    dists[start] <- lengths[start]
    current <- start
    minroutes[[start]] <- start
    
    found <- FALSE
    while(!found) {
      neighbors <- (1:length)[connected[current,]&!visited]
      for(neighbor in neighbors) {
        route.tentative <- c(minroutes[[current]],neighbor)
        dist.tentative <- sum(lengths[route.tentative])
        if(dist.tentative <= dists[neighbor]) {  # changed to LEQ
          dists[neighbor] <- dist.tentative
          minroutes[[neighbor]] <- route.tentative
        }
      }
      if(is.null(minroutes[[current]][1])) {  # this is a hack
        if(stopiferror) stop("No route detected.")
        found <- TRUE
        minroutes[[end]] <- NA
      }
      if(verbose) print(minroutes[[current]])
      visited[current] <- TRUE
      if(current==end) found <- TRUE
      if(!found) current <- which(dists==min(dists[!visited]))[1]
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
            if(level>0) {    ### this is a bit of a hack
              if(length(connected.to.new.prev)>=option[level]) {  # if so, use the different option
                path[level] <- connected.to.new.prev[option[level]]
                good.to.go <-T
              }
            }
            if(level<1) { ## here comes another hack
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
# microbenchmark(detectroute(start=120, end=111, rivers=abstreams),times=10000)
# microbenchmark(detectroute(start=120, end=111, rivers=abstreams, algorithm="segroutes"),times=10000)
# microbenchmark(detectroute2(start=120, end=111, rivers=abstreams),times=10000)
# microbenchmark(detectroute2a(start=120, end=111, rivers=abstreams),times=10000)
# 
# 
# 
# microbenchmark(detectroute(start=120, end=158, rivers=abstreams),times=10000)
# microbenchmark(detectroute(start=120, end=158, rivers=abstreams, algorithm="segroutes"),times=10000)
# microbenchmark(detectroute2(start=120, end=158, rivers=abstreams),times=10000)
# microbenchmark(detectroute2a(start=120, end=158, rivers=abstreams),times=10000)


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
#' tend <- Sys.time()
#' tend - tstart
#' 
#' # after
#' tstart <- Sys.time()
#' detectroute(start=120, end=111, rivers=abstreams)
#' tend <- Sys.time()
#' tend - tstart
#' @importFrom graphics plot
#' @export
buildsegroutes <- function(rivers,lookup=NULL,verbose=FALSE) {
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
    return(length(connections[seg,][(connections[seg,]==1 | connections[seg,]==2) & is.na(connections[seg,])==F]))
  }
  n.bot <- function(seg,connections) {
    return(length(connections[seg,][(connections[seg,]==3 | connections[seg,]==4) & is.na(connections[seg,])==F]))
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
      routes1a <- detectroute(start=i,end=order(order)[rivers$mouth$mouth.seg],rivers=rivers1,algorithm="sequential")
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
      routes_ia <- detectroute(start=i,end=order(order)[rivers$mouth$mouth.seg],rivers=rivers1,algorithm="sequential")
      routes[[i]] <- routes_ia[length(routes_ia):1]
    }
  }
  realroutes <- list()
  for(i in 1:length(rivers$lines)) {
    realroutes[[order[i]]] <- order[routes[[i]]]
  }
  
  rivers$segroutes <- realroutes
  
  if(is.null(lookup) & length(rivers$lines)<=400) lookup <- T
  if(lookup) rivers <- buildlookup(rivers) 

  return(rivers) }
}


# riverdistance_old <- function(startseg=NULL,endseg=NULL,startvert,endvert,rivers,path=NULL,map=FALSE,add=FALSE,stopiferror=TRUE,algorithm=NULL) {
#   # if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
#   
#   connections <- rivers$connections
#   seg.lengths <- rivers$lengths
#   lines <- rivers$lines
#   
#   if(is.null(path)) path <- detectroute(start=startseg,end=endseg,rivers=rivers,stopiferror=stopiferror,algorithm=algorithm)
#   
#   if(is.na(path[1])) return(NA)
#   startseg <- path[1]
#   endseg <- path[length(path)]
#   
#   # if(max(c(startseg,endseg),na.rm=T)>length(rivers$lines) | min(c(startseg,endseg),na.rm=T)<1) {
#   #   stop("Invalid segments specified.")
#   # } 
#   if(startvert>dim(lines[[startseg]])[1] | startvert<1 | endvert>dim(lines[[endseg]])[1] | endvert<1) {
#     stop("Invalid vertex specified.")
#   } 
#   
#   if(map) {
#     if(!add) plot(x=rivers)
#     riverpoints(seg=c(startseg,endseg),vert=c(startvert,endvert),rivers=rivers,pch=15,col=4)
#     if(length(path)>1) {
#       # distance on the partial segment the beginning is on
#       if(connections[startseg,path[2]]<=2) {
#         if(startvert!=1) { 
#           lines(lines[[startseg]][1:startvert,],lwd=3,col=4)
#         }
#       }
#       if(connections[startseg,path[2]]>=3) {
#         linelength <- dim(lines[[startseg]])[1]
#         if(linelength!=startvert) {
#           lines(lines[[startseg]][linelength:startvert,],lwd=3,col=4)
#         }
#       }
#       # distance on the full segments between the beginning and end segments
#       if(length(path)>2) {
#         for(i in path[2:(length(path)-1)]) lines(lines[[i]],col=4,lty=1,lwd=3)
#       }
#       #distance on the partial segment the end is on
#       if(any(connections[path[length(path)-1],endseg]==c(1,3))) {
#         if(endvert!=1) { 
#           lines(lines[[endseg]][1:endvert,],lwd=3,col=4)
#         }
#       }
#       if(any(connections[path[length(path)-1],endseg]==c(2,4))) {
#         linelength <- dim(lines[[endseg]])[1]
#         if(linelength!=endvert) {
#           lines(lines[[endseg]][linelength:endvert,],lwd=3,col=4)
#         }
#       }
#     }
#     # if the beginning and end are on the same segment
#     if(length(path)==1) {
#       min <- min(startvert,endvert)
#       max <- max(startvert,endvert)
#       if(min!=max) {
#         lines(lines[[path[1]]][min:max])
#       }
#     }
#   }
# 
#   route.dist <- 0
#   
#   if(length(path)>1) {
#     # distance on the partial segment the beginning is on
#     if(any(connections[startseg,path[2]]==c(1,2))) {
#       if(startvert!=1) { 
#         route.dist <- route.dist + pdisttot(lines[[startseg]][1:startvert,])
#       }
#     }
#     if(any(connections[startseg,path[2]]==c(3,4))) {
#       linelength <- dim(lines[[startseg]])[1]
#       if(linelength!=startvert) {
#         route.dist <- route.dist + pdisttot(lines[[startseg]][linelength:startvert,])
#       }
#     } 
#     # distance on the full segments between the beginning and end segments
#     if(length(path)>2) {
#       route.dist <- route.dist+sum(seg.lengths[path[2:(length(path)-1)]])
#     } 
#     #distance on the partial segment the end is on
#     if(any(connections[path[length(path)-1],endseg]==c(1,3))) {
#       if(endvert!=1) { 
#         route.dist <- route.dist + pdisttot(lines[[endseg]][1:endvert,])
#       }
#     }
#     if(any(connections[path[length(path)-1],endseg]==c(2,4))) {
#       linelength <- dim(lines[[endseg]])[1]
#       if(linelength!=endvert) {
#         route.dist <- route.dist + pdisttot(lines[[endseg]][linelength:endvert,])
#       }
#     } 
#   }
#   # if the beginning and end are on the same segment
#   if(length(path)==1) {
#     # min <- min(startvert,endvert)
#     # max <- max(startvert,endvert)
#     if(startvert!=endvert) {
#       route.dist <- route.dist + pdisttot(lines[[startseg]][startvert:endvert,])
#     }
#   } 
#   
#   return(route.dist)
# }
# 
# summary(microbenchmark(riverdistance(startseg=120, endseg=111 ,startvert=120, endvert=120, rivers=abstreams),times=1000))[4:5]
# summary(microbenchmark(detectroute(start=120, end=111, rivers=abstreams),times=1000))[4:5]
# summary(microbenchmark(riverdistance2(startseg=120, endseg=111 ,startvert=120, endvert=120, rivers=abstreams),times=1000))[4:5]
# summary(microbenchmark(riverdistance2a(startseg=120, endseg=111 ,startvert=120, endvert=120, rivers=abstreams),times=1000))[4:5]
# 
# 
# 
# summary(microbenchmark(riverdistance(startseg=1, endseg=14 ,startvert=120, endvert=120, rivers=Gulk),times=1000))[4:5]
# summary(microbenchmark(detectroute(start=1, end=14, rivers=Gulk),times=1000))[4:5]
# summary(microbenchmark(riverdistance2(startseg=1, endseg=14 ,startvert=120, endvert=120, rivers=Gulk),times=1000))[4:5]
# summary(microbenchmark(riverdistance2a(startseg=1, endseg=14 ,startvert=120, endvert=120, rivers=Gulk),times=1000))[4:5]
# 
# 
# addcumuldist <- function(rivers) {
#   cumuldist <- list()
#   for(i in 1:length(rivers$lines)) {
#     xy <- rivers$lines[[i]]
#     n <- dim(xy)[1]
#     cumuldist[[i]] <- c(0,cumsum(sqrt(((xy[1:(n-1),1] - xy[2:n,1])^2) + ((xy[1:(n-1),2] - xy[2:n,2])^2))))
#   }
#   rivers$cumuldist <- cumuldist
#   return(rivers)
# }
# Gulk_boom <- addcumuldist(Gulk)

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
#' @author Matt Tyers
#' @note Building routes from the river mouth to each river network segment may
#'   greatly reduce computation time (see \link{buildsegroutes}).
#' @examples
#' data(Gulk)
#' riverdistance(startseg=6, endseg=14, startvert=100, endvert=200, rivers=Gulk)
#' riverdistance(startvert=100, endvert=200, path=c(6,3,4,10,11,14), rivers=Gulk)
#' riverdistance(startseg=6, endseg=14, startvert=100, endvert=200, rivers=Gulk, map=TRUE)
#'       
#' # speed comparison: before and after building routes for each segment...
#' data(abstreams)
#' plot(x=abstreams)
#' abstreams1 <- abstreams
#' abstreams1$segroutes <- NULL #taking out the $segroutes component
#' 
#' # before
#' tstart <- Sys.time()
#' riverdistance(startseg=120, endseg=111, startvert=20, endvert=20, rivers=abstreams1)
#' tend <- Sys.time()
#' tend - tstart
#' 
#' # after
#' tstart <- Sys.time()
#' riverdistance(startseg=120, endseg=111 ,startvert=20, endvert=20, rivers=abstreams)
#' tend <- Sys.time()
#' tend - tstart
#' 
#' @export
riverdistance <- function(startseg=NULL,endseg=NULL,startvert,endvert,rivers,path=NULL,map=FALSE,add=FALSE,stopiferror=TRUE,algorithm=NULL) {
  # if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  
  if(!is.null(rivers$distlookup) & !map) {
    cumuldist <- rivers$cumuldist
    lengths <- rivers$lengths
    x<-rivers$distlookup
    if(startseg==endseg) dist <- abs(cumuldist[[startseg]][startvert] - cumuldist[[startseg]][endvert])
    else {
      if(x$starttop[startseg,endseg] & x$endtop[startseg,endseg]) dist <- x$middist[startseg,endseg] + cumuldist[[startseg]][startvert] + cumuldist[[endseg]][endvert]
      if(x$starttop[startseg,endseg] & !x$endtop[startseg,endseg]) dist <- x$middist[startseg,endseg] + cumuldist[[startseg]][startvert] + lengths[endseg] - cumuldist[[endseg]][endvert]
      if(!x$starttop[startseg,endseg] & x$endtop[startseg,endseg]) dist <- x$middist[startseg,endseg] + lengths[startseg] - cumuldist[[startseg]][startvert] + cumuldist[[endseg]][endvert]
      if(!x$starttop[startseg,endseg] & !x$endtop[startseg,endseg]) dist <- x$middist[startseg,endseg] + lengths[startseg] - cumuldist[[startseg]][startvert] + lengths[endseg] - cumuldist[[endseg]][endvert]
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
  
  # if(max(c(startseg,endseg),na.rm=T)>length(rivers$lines) | min(c(startseg,endseg),na.rm=T)<1) {
  #   stop("Invalid segments specified.")
  # } 
  if(startvert>dim(lines[[startseg]])[1] | startvert<1 | endvert>dim(lines[[endseg]])[1] | endvert<1) {
    stop("Invalid vertex specified.")
  } 
  
  if(map) {
    if(!add) plot(x=rivers)
    riverpoints(seg=c(startseg,endseg),vert=c(startvert,endvert),rivers=rivers,pch=15,col=4)
    if(length(path)>1) {
      # distance on the partial segment the beginning is on
      if(connections[startseg,path[2]]<=2) {
        if(startvert!=1) { 
          lines(lines[[startseg]][1:startvert,],lwd=3,col=4)
        }
      }
      if(connections[startseg,path[2]]>=3) {
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
        # route.dist <- route.dist + pdisttot(lines[[startseg]][linelength:startvert,])
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
        # route.dist <- route.dist + pdisttot(lines[[endseg]][1:endvert,])
        route.dist <- route.dist + cumuldist[[endseg]][endvert]
      }
    }
    if(any(connections[path[length(path)-1],endseg]==c(2,4))) {
      linelength <- dim(lines[[endseg]])[1]
      if(linelength!=endvert) {
        # route.dist <- route.dist + pdisttot(lines[[endseg]][linelength:endvert,])
        route.dist <- route.dist + seg.lengths[endseg] - cumuldist[[endseg]][endvert]
      }
    } 
  }
  # if the beginning and end are on the same segment
  if(length(path)==1) {
    min <- min(startvert,endvert)
    max <- max(startvert,endvert)
    if(startvert!=endvert) {
      # route.dist <- route.dist + pdisttot(lines[[startseg]][startvert:endvert,])
      route.dist <- route.dist + cumuldist[[startseg]][max] - cumuldist[[startseg]][min]
    }
  } 
  
  return(route.dist)
}
# 
# summary(microbenchmark(riverdistance(startseg=120, endseg=111 ,startvert=120, endvert=120, rivers=abstreams),times=1000))[4:5]
# summary(microbenchmark(detectroute(start=120, end=111, rivers=abstreams),times=1000))[4:5]
# summary(microbenchmark(riverdistance2(startseg=120, endseg=111 ,startvert=120, endvert=120, rivers=abstreams),times=1000))[4:5]
# summary(microbenchmark(riverdistance3(startseg=120, endseg=111 ,startvert=120, endvert=120, rivers=abstreams_boom),times=1000))[4:5]
# 
# 
# 
# summary(microbenchmark(riverdistance(startseg=1, endseg=14 ,startvert=120, endvert=120, rivers=Gulk),times=1000))[4:5]
# summary(microbenchmark(detectroute(start=1, end=14, rivers=Gulk),times=1000))[4:5]
# summary(microbenchmark(riverdistance2(startseg=1, endseg=14 ,startvert=120, endvert=120, rivers=Gulk),times=1000))[4:5]
# summary(microbenchmark(riverdistance3(startseg=1, endseg=14 ,startvert=120, endvert=120, rivers=Gulk_boom),times=1000))[4:5]
# 

# buildallroutes <- function(rivers) {
#   length <- length(rivers$lines)
#   allroutes <- list()
#   for(i in 1:length) {
#     allroutes[[i]] <- list()
#     for(j in 1:length) {
#       allroutes[[i]][[j]] <- detectroute(start=i, end=j, rivers=rivers)
#     }
#   }
#   return(allroutes)
# }
# Gulk_allroutes <- buildallroutes(Gulk)
# abstreams_allroutes <- buildallroutes(abstreams)
# 
# library(microbenchmark)
# summary(microbenchmark(riverdistance(startseg=1, endseg=14 ,startvert=120, endvert=120, rivers=Gulk),times=1000))[4:5]
# summary(microbenchmark(riverdistance(startseg=1, endseg=14 ,startvert=120, endvert=120, rivers=Gulk, path=Gulk_allroutes[[1]][[14]]),times=1000))[4:5]
# 
# summary(microbenchmark(riverdistance(startseg=120, endseg=111 ,startvert=120, endvert=120, rivers=abstreams),times=1000))[4:5]
# summary(microbenchmark(riverdistance(startseg=120, endseg=111 ,startvert=120, endvert=120, rivers=abstreams, path=abstreams_allroutes[[120]][[111]]),times=1000))[4:5]
# 
# calcalldistances <- function(rivers) {
#   length <- length(rivers$lines)
#   alldists <- list()
#   for(i in 1:length) {
#     alldists[[i]] <- list()
#     for(j in 1:length) {
#       theroute <- detectroute(start=i, end=j, rivers=rivers)
#       # nverti <- dim(rivers$lines[[i]])[1]
#       # nvertj <- dim(rivers$lines[[j]])[1]
#       # distmat <- matrix(nrow=nverti,ncol=nvertj)
#       # for(verti in 1:nverti) {
#       #   for(vertj in 1:nvertj) {
#       #     distmat[i,j] <- riverdistance(startseg=i, endseg=j, startvert=verti, endvert=vertj, path=theroute, rivers=rivers)
#       #   }
#       # }
#       # alldists[[i]][[j]] <- distmat
#       # print(c(i,j))
#     }
#   }
#   return(alldists)
# }
# system.time(Gulk_alldists <- calcalldistances(Gulk))
# system.time(abstreams_alldists <- calcalldistances(abstreams))
# 
# 
# summary(microbenchmark(riverdistance(startseg=1, endseg=14 ,startvert=120, endvert=120, rivers=Gulk),times=1000))[4:5]
# summary(microbenchmark(Gulk_alldists[[1]][[14]][120,120],times=1000))[4:5]
# 
# calcdisttable <- function(rivers) {
#   length <- length(rivers$lines) 
#   connections <- rivers$connections
#   middistmat <- matrix(0,nrow=length,ncol=length)
#   startheadmat <- endheadmat <- matrix(NA,nrow=length,ncol=length)
#   for(i in 1:length) {
#     for(j in 1:length) {
#       theroute <- detectroute(start=i, end=j, rivers=rivers)
#       routelength <- length(theroute)
#       if(length(theroute) > 2) middistmat[i,j] <- sum(rivers$lengths[theroute[2:(routelength-1)]])
#       if(length(theroute) > 1) {
#         startheadmat[i,j] <- ifelse(any(connections[i, theroute[2]] == c(1,2)), T, F)
#         endheadmat[i,j] <- ifelse(any(connections[theroute[routelength-1], j] == c(1,3)), T, F)
#       }
#     }
#   }
#   return(list(middist=middistmat, starttop=startheadmat, endtop=endheadmat))
# }
# system.time(Gulk_disttable <- calcdisttable(Gulk))
# system.time(abstreams_disttable <- calcdisttable(abstreams))
# 
# riverdistancefromtable <- function(startseg,startvert,endseg,endvert,rivers,x) {
#   if(startseg==endseg) dist <- abs(rivers$cumuldist[start] - rivers$cumuldist[end])
#   else {
#     if(x$starttop[startseg,endseg] & x$endtop[startseg,endseg]) dist <- x$middist[startseg,endseg] + rivers$cumuldist[[startseg]][startvert] + rivers$cumuldist[[endseg]][endvert]
#     if(x$starttop[startseg,endseg] & !x$endtop[startseg,endseg]) dist <- x$middist[startseg,endseg] + rivers$cumuldist[[startseg]][startvert] + rivers$lengths[endseg] - rivers$cumuldist[[endseg]][endvert]
#     if(!x$starttop[startseg,endseg] & x$endtop[startseg,endseg]) dist <- x$middist[startseg,endseg] + rivers$lengths[startseg] - rivers$cumuldist[[startseg]][startvert] + rivers$cumuldist[[endseg]][endvert]
#     if(!x$starttop[startseg,endseg] & !x$endtop[startseg,endseg]) dist <- x$middist[startseg,endseg] + rivers$lengths[startseg] - rivers$cumuldist[[startseg]][startvert] + rivers$lengths[endseg] - rivers$cumuldist[[endseg]][endvert]
#   }
#   return(dist)
# }
# 
# riverdistancefromtable <- function(startseg,startvert,endseg,endvert,rivers,x) {
#   cumuldist <- rivers$cumuldist
#   lengths <- rivers$lengths
#   if(startseg==endseg) dist <- abs(cumuldist[start] - cumuldist[end])
#   else {
#     if(x$starttop[startseg,endseg] & x$endtop[startseg,endseg]) dist <- x$middist[startseg,endseg] + cumuldist[[startseg]][startvert] + cumuldist[[endseg]][endvert]
#     if(x$starttop[startseg,endseg] & !x$endtop[startseg,endseg]) dist <- x$middist[startseg,endseg] + cumuldist[[startseg]][startvert] + lengths[endseg] - cumuldist[[endseg]][endvert]
#     if(!x$starttop[startseg,endseg] & x$endtop[startseg,endseg]) dist <- x$middist[startseg,endseg] + lengths[startseg] - cumuldist[[startseg]][startvert] + cumuldist[[endseg]][endvert]
#     if(!x$starttop[startseg,endseg] & !x$endtop[startseg,endseg]) dist <- x$middist[startseg,endseg] + lengths[startseg] - cumuldist[[startseg]][startvert] + lengths[endseg] - cumuldist[[endseg]][endvert]
#   }
#   return(dist)
# }
# riverdistance(startseg=1, endseg=14 ,startvert=120, endvert=120, rivers=Gulk)
# riverdistancefromtable(startseg=1, endseg=14 ,startvert=120, endvert=120, rivers=Gulk, x=Gulk_disttable)
# riverdistance(startseg=120, endseg=111 ,startvert=120, endvert=120, rivers=abstreams)
# riverdistancefromtable(startseg=120, endseg=111 ,startvert=120, endvert=120, rivers=abstreams,x=abstreams_disttable)
# 
# microbenchmark(riverdistance(startseg=1, endseg=14 ,startvert=120, endvert=120, rivers=Gulk),times=1000)
# microbenchmark(riverdistancefromtable(startseg=1, endseg=14 ,startvert=120, endvert=120, rivers=Gulk, x=Gulk_disttable),times=1000)
# microbenchmark(riverdistance(startseg=120, endseg=111 ,startvert=120, endvert=120, rivers=abstreams),times=1000)
# microbenchmark(riverdistancefromtable(startseg=120, endseg=111 ,startvert=120, endvert=120, rivers=abstreams,x=abstreams_disttable),times=1000)


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
#' @examples
#' data(abstreams)
#' 
#' abstreams1 <- buildlookup(abstreams)
#' @export
buildlookup <- function(rivers) {
  if(is.null(rivers$segroutes) & interactive()) pb <- txtProgressBar(style=3)
  if(!is.na(rivers$braided)) {
    if(rivers$braided) warning("Braiding detected in river network - all distances used in analysis will be shortest-path distances, not necessarily travel distances.")
  } else {
    warning("Braiding may be present in river network.  If so, all distances used in analysis will be shortest-path distances, not necessarily travel distances.")
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
        starttop[i,j] <- ifelse(any(connections[i, theroute[2]] == c(1,2)), T, F)
        endtop[i,j] <- ifelse(any(connections[theroute[routelength-1], j] == c(1,3)), T, F)
      }
      if(interactive() & is.null(rivers$segroutes)) setTxtProgressBar(pb=pb, value=(i/length + j/length/length))
    }
  }
  distlookup <- list(middist=middist, starttop=starttop, endtop=endtop)
  rivers$distlookup <- distlookup
  
  
  return(rivers)
}
