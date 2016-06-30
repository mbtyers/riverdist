#' Check for Braiding in a River Network
#' @description Detects braiding (multiple flow channels between two locations)
#'   within a river network object.  Braiding can either be checked for in the
#'   route between two segments, or in the river network as a whole.
#' @param rivers The river network object to check.
#' @param startseg Starting segment of a route to investigate.  If this and
#'   \code{endseg} are \code{NULL}, the full river network will be checked.
#' @param endseg Starting segment of a route to investigate.  If this and
#'   \code{startseg} are \code{NULL}, the full river network will be checked.
#' @author Matt Tyers
#' @note This function is called within \link{cleanup}, which is recommended in most cases.
#' @examples
#' data(Gulk)
#' plot(x=Gulk)
#' checkbraided(rivers=Gulk)
#' 
#' data(KilleyW)
#' plot(x=KilleyW)
#' checkbraided(rivers=KilleyW)
#' 
#' Kenai3.subset <- trimriver(trimto=c(22,2,70,30,15,98,96,89,52,3), rivers=Kenai3)
#' plot(x=Kenai3.subset)
#' 
#' checkbraided(startseg=1, endseg=7, rivers=Kenai3.subset)
#' checkbraided(startseg=1, endseg=5, rivers=Kenai3.subset)
#' @export
checkbraided <- function(rivers,startseg=NULL,endseg=NULL) {
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  connections <- rivers$connections
  length <- length(rivers$lines)
  
  invertrivers <- rivers
  invertrivers$connections <- connections[length:1,length:1]
  invertrivers$lines <- rivers$lines[length:1]
  
  if(sum(is.null(startseg),is.null(endseg))==1) {
    stop("Error - need to specify both starting and ending segments, or neither")
  }
  
  if(is.null(startseg) & is.null(endseg)) {
    if(interactive()) pb <- txtProgressBar(style=3)
    finished <- FALSE
    braiding <- FALSE
    i <- 1
    j <- 1
    
    while(!finished & !braiding) {
      if(i!=j) {
        route1 <- detectroute(start=i,end=j,rivers=rivers,algorithm="sequential")
        route2 <- detectroute(start=(length-i+1),end=(length-j+1),rivers=invertrivers,algorithm="sequential")
        route2 <- length-route2+1
        if(length(route1) != length(route2)) braiding<-T
        if(length(route1) == length(route2)) if(any(route1 != route2)) braiding<-T
      }
      if(i==length & j==length) finished<-T
      if(i==length) j<-j+1
      if(i<length) i<-i+1
      if(interactive()) setTxtProgressBar(pb=pb, value=i/length)
    }
    if(interactive()) setTxtProgressBar(pb=pb, value=1)
    if(braiding) cat('\n',"Braiding detected in river network.  Distance measurements may be inaccurate.")
    if(finished & !braiding) cat('\n',"No braiding detected in river network.")
  }
  
  if(!is.null(startseg) & !is.null(endseg)) {
    if(max(c(startseg,endseg),na.rm=T)>length(rivers$lines) | min(c(startseg,endseg),na.rm=T)<1) {
      stop("Invalid segments specified.")
    }
    route1 <- detectroute(start=startseg,end=endseg,rivers=rivers,algorithm="sequential")
    route2 <- detectroute(start=(length-startseg+1),end=(length-endseg+1),rivers=invertrivers,algorithm="sequential")
    route2 <- length-route2+1
    if(length(route1) != length(route2)) cat("Braiding detected between segments.  Distance measurements may be inaccurate.")
    if(length(route1) == length(route2)) if(any(route1 != route2)) cat("Braiding detected between segments.  Distance measurements may be inaccurate.")
    if(length(route1) == length(route2)) if(all(route1 == route2)) cat("No braiding detected between segments.")
  }
}


#' Check for Braiding in a River Network
#' @description Detects braiding (multiple flow channels between two locations)
#'   within a river network object, and returns a logical value for specifying braiding within a river network object.
#' @param rivers The river network object to check.
#' @param toreturn Specifying \code{toreturn="rivers"} (the default) will return a river network object with a value of \code{TRUE} or \code{FALSE} assigned to the \code{$braided} element of the river network object.  Specifying \code{toreturn="logical"} will just return \code{TRUE} if braiding is detected or \code{FALSE} if no braiding is detected.  Specifying \code{toreturn="routes"} will return the first two differing routes detected, which may be useful in identifying where the problem lies.
#' @author Matt Tyers
#' @note This function is called within \link{cleanup}, which is recommended in most cases.
#' @examples
#' data(Gulk,KilleyW)
#' Gulk <- setmouth(seg=1, vert=1, rivers=Gulk)
#' plot(x=Gulk)
#' checkbraidedTF(rivers=Gulk, toreturn="logical")
#' 
#' KilleyW <- setmouth(seg=1, vert=288, rivers=KilleyW)
#' plot(x=KilleyW)
#' checkbraidedTF(rivers=KilleyW, toreturn="logical")
#' checkbraidedTF(rivers=KilleyW, toreturn="routes")
#' 
#' KilleyW.1 <- checkbraidedTF(rivers=KilleyW, toreturn="rivers")
#' str(KilleyW.1)
#' @export
checkbraidedTF <- function(rivers,toreturn="rivers") {
  if(toreturn != "rivers" & toreturn != "logical" & toreturn != "routes") stop("Invalid specification of argument 'toreturn'.  See help for more details.")
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  connections <- rivers$connections
  length <- length(rivers$lines)
  mouthseg <- rivers$mouth$mouth.seg
  if(is.na(mouthseg)) stop("Mouth must be specified.")
  lines <- rivers$lines
  tolerance <- rivers$tolerance
  if(length==1) {
    braiding <- F
  } else{
    
  n.top <- function(seg,connections) {
    return(length(connections[seg,][(connections[seg,]==1 | connections[seg,]==2 | connections[seg,]==5 | connections[seg,]==6) & is.na(connections[seg,])==F]))
  }
  n.bot <- function(seg,connections) {
    return(length(connections[seg,][(connections[seg,]==3 | connections[seg,]==4 | connections[seg,]==5 | connections[seg,]==6) & is.na(connections[seg,])==F]))
  }
  
  invertrivers <- rivers
  invertrivers$connections <- connections[length:1,length:1]
  invertrivers$lines <- rivers$lines[length:1]
  
  ntop <- nbot <- NA
  for(i in 1:length) {
    ntop[i] <- n.top(i,connections = connections)
    nbot[i] <- n.bot(i,connections = connections)
  }
  checkthese <- (1:length)[(ntop==0&nbot!=0) | (ntop!=0&nbot==0)]
  braiding <- F 
  finished <- F
  i<-1
  if(interactive()) pb <- txtProgressBar(style=3)
  while(!finished) {
    route1 <- detectroute(start=checkthese[i],end=mouthseg,rivers=rivers,algorithm="sequential")
    route2 <- detectroute(start=(length-checkthese[i]+1),end=(length-mouthseg+1),rivers=invertrivers,algorithm="sequential")
    route2 <- length-route2+1
    if(length(route1) != length(route2)) {
      braiding<-T
      finished<-T
    }
    if(length(route1) == length(route2)) if(any(route1 != route2)) {
      braiding<-T
      finished<-T
    }
    if(i==length(checkthese)) finished <- T
    i<-i+1
    #print(i)
    if(interactive()) setTxtProgressBar(pb=pb, value=i/length)
  }
  }
  
  if(interactive()) setTxtProgressBar(pb=pb, value=1)
  rivers$braided <- braiding
  if(toreturn=="logical") return(braiding)
  if(toreturn=="rivers") return(rivers)
  if(toreturn=="routes") {
    if(braiding) return(list(route1=route1, route2=route2))
  }
}

#' Detect Multiple Routes
#' @description Called internally within \link{riverdistancelist}.  Detects many
#'   sequential routes from one river network segment to another, in the event
#'   of braiding.  Different routes are detected by randomly reordering the
#'   segment numbers of the input river network object, thus changing the
#'   internal hierarchy of segment selection.
#' @param startseg Segment number of the start of the route
#' @param endseg Segment number of the end of the route
#' @param rivers The river network object to use
#' @param reps Number of randomized reorderings to try.  Larger numbers will 
#'   result in a more comprehensive list of routes, but computation can be slow,
#'   particularly in the presence of a complex river network.
#' @return A list of vectors, each describing a route in segment numbers.
#' @note Since this function uses randomization, there is no guarantee that the 
#'   list of routes will be comprehensive.  Larger numbers of \code{reps} can be
#'   tried, but computation can be slow, particularly in the presence of a
#'   complex river network.  It may be advantageous to use \link{trimriver} to
#'   create a smaller, more specific river network object to work with.
#' @author Matt Tyers
#' @examples
#' data(KilleyW)
#' plot(x=KilleyW)
#' 
#' routelist(startseg=1, endseg=16, rivers=KilleyW, reps=1000)
#' @export
routelist <- function(startseg,endseg,rivers,reps=100) {
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  if(max(c(startseg,endseg),na.rm=T)>length(rivers$lines) | min(c(startseg,endseg),na.rm=T)<1) stop("Invalid segments specified.")
  connections <- rivers$connections
  length <- length(rivers$lines)
  routes <- list()
  
  for(i in 1:reps) {
    codex <- sample(1:length,length)
    scrambled <- rivers
    scrambled$segroutes <- NULL
    scrambled$connections <- rivers$connections[codex,codex]
    scrambled$lines <- rivers$lines[codex]
    scrambled.route <- detectroute(start=order(codex)[startseg],end=order(codex)[endseg],rivers=scrambled,algorithm="sequential")
    routes[[i]] <- codex[scrambled.route]
  }
  unique.routes <- unique(routes)
  return(unique.routes)
}


#' Multiple River Distances
#' @description Used to calculate a list of possible river distances, in the
#'   event of braiding.  Calls \link{routelist} to detect a list of routes from
#'   one river location to another, and uses \link{riverdistance} to calculate
#'   the distances along those routes.  Different routes are detected by
#'   randomly reordering the segment numbers of the input river network object,
#'   thus changing the internal hierarchy of segment selection.
#' @param startseg Segment number of the start of the route
#' @param endseg Segment number of the end of the route
#' @param startvert Vertex number of the start of the route
#' @param endvert Vertex number of the end of the route
#' @param rivers The river network object to use
#' @param reps Number of randomized reorderings to try.  Larger numbers will
#'   result in a more comprehensive list of routes, but computation can be slow,
#'   particularly in the presence of a complex river network.
#' @return A list with two objects, \code{$routes} being a list of detected routes in
#'   ascending order by distance, and \code{$distances} being the respective distances
#'   along the routes detected.
#' @note Since this function uses randomization, there is no guarantee that the
#'   list of routes will be comprehensive.  Larger numbers of reps can be tried,
#'   but computation can be slow, particularly in the presence of a complex
#'   river network.  It may be advantageous to use \link{trimriver} to create a
#'   smaller, more specific river network object to work with.
#' @author Matt Tyers
#' @examples
#' data(KilleyW)
#' plot(x=KilleyW)
#' 
#' Killey.dists <- riverdistancelist(startseg=1, endseg=16, startvert=100, endvert=25,
#'    rivers=KilleyW, reps=1000)
#' Killey.dists  # 18 routes are detected.
#' 
#' # mapping the shortest route detected... 
#' riverdistance(startvert=100, endvert=25, path=Killey.dists$routes[[1]], rivers=KilleyW, map=TRUE)
#' 
#' # mapping the shortest longest detected... 
#' riverdistance(startvert=100, endvert=25, path=Killey.dists$routes[[18]], rivers=KilleyW, map=TRUE)
#' @export
riverdistancelist <- function(startseg,endseg,startvert,endvert,rivers,reps=100) {
  routes <- routelist(startseg=startseg,endseg=endseg,reps=reps,rivers=rivers)
  dists <- NA
  for(i in 1:length(routes)) {
    dists[i] <- riverdistance(startvert=startvert,endvert=endvert,path=routes[[i]],rivers=rivers)
  }
  routes <- routes[order(dists)]
  dists <- dists[order(dists)]
  out <- list(routes,dists)  
  names(out) <- c("routes","distances")
  return(out)
}