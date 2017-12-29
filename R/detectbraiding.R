# checkbraided <- function(rivers,startseg=NULL,endseg=NULL,progress=TRUE) {
#   if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
#   connections <- rivers$connections
#   length <- length(rivers$lines)
#   
#   invertrivers <- rivers
#   invertrivers$connections <- connections[length:1,length:1]
#   invertrivers$lines <- rivers$lines[length:1]
#   
#   if(sum(is.null(startseg),is.null(endseg))==1) {
#     stop("Error - need to specify both starting and ending segments, or neither")
#   }
#   
#   if(is.null(startseg) & is.null(endseg)) {
#     if(interactive() & progress) pb <- txtProgressBar(style=3)
#     finished <- FALSE
#     braiding <- FALSE
#     i <- 1
#     j <- 1
#     
#     while(!finished & !braiding) {
#       if(i!=j) {
#         route1 <- detectroute(start=i,end=j,rivers=rivers,algorithm="sequential")
#         route2 <- detectroute(start=(length-i+1),end=(length-j+1),rivers=invertrivers,algorithm="sequential")
#         route2 <- length-route2+1
#         if(length(route1) != length(route2)) braiding<-T
#         if(length(route1) == length(route2)) if(any(route1 != route2)) braiding<-T
#       }
#       if(i==length & j==length) finished<-T
#       if(i==length) j<-j+1
#       if(i<length) i<-i+1
#       if(interactive() & progress) setTxtProgressBar(pb=pb, value=i/length)
#     }
#     if(interactive() & progress) setTxtProgressBar(pb=pb, value=1)
#     if(braiding) cat('\n',"Braiding detected in river network.  Distance measurements may be inaccurate.")
#     if(finished & !braiding) cat('\n',"No braiding detected in river network.")
#   }
#   
#   if(!is.null(startseg) & !is.null(endseg)) {
#     if(max(c(startseg,endseg),na.rm=T)>length(rivers$lines) | min(c(startseg,endseg),na.rm=T)<1) {
#       stop("Invalid segments specified.")
#     }
#     route1 <- detectroute(start=startseg,end=endseg,rivers=rivers,algorithm="sequential")
#     route2 <- detectroute(start=(length-startseg+1),end=(length-endseg+1),rivers=invertrivers,algorithm="sequential")
#     route2 <- length-route2+1
#     if(length(route1) != length(route2)) cat("Braiding detected between segments.  Distance measurements may be inaccurate.")
#     if(length(route1) == length(route2)) if(any(route1 != route2)) cat("Braiding detected between segments.  Distance measurements may be inaccurate.")
#     if(length(route1) == length(route2)) if(all(route1 == route2)) cat("No braiding detected between segments.")
#   }
# }



#' Check for Braiding in a River Network
#' @description Detects braiding (multiple flow channels between two locations)
#'   within a river network object.  Braiding can either be checked for in the
#'   route between two segments, or in the river network as a whole.
#' @param rivers The river network object to check.
#' @param startseg Starting segment of a route to investigate.  If this and
#'   \code{endseg} are \code{NULL}, the full river network will be checked.
#' @param endseg Starting segment of a route to investigate.  If this and
#'   \code{startseg} are \code{NULL}, the full river network will be checked.
#' @param progress Whether to show the progress bar.  Defaults to \code{TRUE}.
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
checkbraided <- function(rivers,startseg=NULL,endseg=NULL,progress=TRUE) {
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  
  if(is.null(startseg) | is.null(endseg)) {
    braiding <- checkbraidedTF(rivers=rivers,progress=progress,toreturn="logical")
    if(braiding) cat('\n',"Braiding detected in river network.  Distance measurements may be inaccurate.")
    if(!braiding) cat('\n',"No braiding detected in river network.")
  }
  else {
    braiding <- F
    theroute <- detectroute(start=startseg, end=endseg, rivers=rivers, stopiferror=F, algorithm="Dijkstra")           ##################################################################
    if(length(theroute)>2) {
      jdone <- F
      j <- 2
      while(!jdone) {
        rivers_except_not <- rivers
        rivers_except_not$connections[theroute[j],] <- NA
        rivers_except_not$connections[,theroute[j]] <- NA
        trialroute <- detectroute(start=startseg, end=endseg, rivers=rivers_except_not, stopiferror=F, algorithm="Dijkstra")           ##################################################################
        if(!is.na(trialroute[1])) {
          braiding <- T
          finished <- T  ###############
          jdone <- T
        }
        j <- j+1
        if(j>=length(theroute)) jdone<-T
      }
    }
    if(braiding) cat("Braiding detected between segments.  Distance measurements may be inaccurate.")
    else cat("No braiding detected between segments.")
  }
  
}


#' Check for Braiding in a River Network
#' @description Detects braiding (multiple flow channels between two locations)
#'   within a river network object, and returns a logical value for specifying braiding within a river network object.
#' @param rivers The river network object to check.
#' @param toreturn Specifying \code{toreturn="rivers"} (the default) will return a river network object with a value of \code{TRUE} or \code{FALSE} assigned to the \code{$braided} element of the river network object.  Specifying \code{toreturn="logical"} will just return \code{TRUE} if braiding is detected or \code{FALSE} if no braiding is detected.  Specifying \code{toreturn="routes"} will return the first two differing routes detected, which may be useful in identifying where the problem lies.
#' @param progress Whether to show the progress bar.  Defaults to \code{TRUE}.
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
checkbraidedTF <- function(rivers,toreturn="rivers",progress=TRUE) {
  # if(toreturn != "rivers" & toreturn != "logical" & toreturn != "routes") stop("Invalid specification of argument 'toreturn'.  See help for more details.")
  if(!(toreturn %in% c("rivers","logical","routes","allroutes"))) stop("Invalid specification of argument 'toreturn'.  See help for more details.")   #################
  if(toreturn=="allroutes") {
    allroutes <- list()
    allroutesi <- 1
  }
  
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  connections <- rivers$connections
  length <- length(rivers$lines)
  mouthseg <- ifelse(!is.na(rivers$mouth$mouth.seg),rivers$mouth$mouth.seg,1)
  # if(is.na(mouthseg)) stop("Mouth must be specified.")
  lines <- rivers$lines
  tolerance <- rivers$tolerance
  if(length==1) {
    braiding <- F
  } else{
    
  ntop <- rowSums(connections==1,na.rm=T) + rowSums(connections==2,na.rm=T) + rowSums(connections==5,na.rm=T) + rowSums(connections==6,na.rm=T)
  nbot <- rowSums(connections==3,na.rm=T) + rowSums(connections==4,na.rm=T) + rowSums(connections==5,na.rm=T) + rowSums(connections==6,na.rm=T)
  
  checkthese <- (1:length)[xor(ntop==0,nbot==0)]
  
  braiding <- F 
  finished <- F
  i<-1
  if(interactive() & progress) pb <- txtProgressBar(style=3)
  while(!finished) {
    theroute <- detectroute(start=checkthese[i], end=mouthseg, rivers=rivers, stopiferror=F, algorithm="Dijkstra")           ##################################################################
    if(length(theroute)>2) {
      jdone <- F
      j <- 2
      while(!jdone) {
        rivers_except_not <- rivers
        rivers_except_not$connections[theroute[j],] <- NA
        rivers_except_not$connections[,theroute[j]] <- NA
        trialroute <- detectroute(start=checkthese[i], end=mouthseg, rivers=rivers_except_not, stopiferror=F, algorithm="Dijkstra")           ##################################################################
        if(!is.na(trialroute[1])) {
          braiding <- T
          if(toreturn!="allroutes") finished <- T  ###############
          else {
            allroutes[[allroutesi]] <- c(setdiff(theroute,trialroute),setdiff(trialroute,theroute))
            allroutesi <- allroutesi+1
          }  #################
          jdone <- T
          route1 <- theroute
          route2 <- trialroute
        }
        j <- j+1
        if(j>=length(theroute)) jdone<-T
      }
    }
    if(i>=length(checkthese)) finished<-T
    i <- i+1
    
    if(interactive() & progress) setTxtProgressBar(pb=pb, value=i/(length(checkthese)))
  }
  }
  
  if(interactive() & progress & length>1) setTxtProgressBar(pb=pb, value=1)
  rivers$braided <- braiding
  if(toreturn=="logical") return(braiding)
  if(toreturn=="rivers") return(rivers)
  if(toreturn=="routes") {
    if(braiding) return(list(route1=route1, route2=route2))
  }
  if(toreturn=="allroutes") {   ##################
    if(braiding) return(allroutes)
  }  ############
}

#' Detect Multiple Routes
#' @description Called internally within \link{riverdistancelist}.  Detects all possible routes from one river network segment to another, in the event
#'   of braiding.  
#' @param startseg Segment number of the start of the route
#' @param endseg Segment number of the end of the route
#' @param rivers The river network object to use
#' @param reps Deprecated.  Was used in a previous version using randomization.
#' @return A list of vectors, each describing a route in segment numbers.
#' @note The previous version of this function returned many possible routes using randomization - this algorithm now computes all possible routes.
#' @author Matt Tyers
#' @examples
#' data(KilleyW)
#' plot(x=KilleyW)
#' 
#' routelist(startseg=1, endseg=16, rivers=KilleyW, reps=1000)
#' @export
routelist <- function(startseg,endseg,rivers,reps=100) {
  routes <- list()
  routes[[1]] <- detectroute(start=startseg,end=endseg,rivers=rivers,algorithm="Dijkstra",stopiferror=F)
  routetoinvestigate <- 1
  routetoadd <- 2
  done <- F
  brokenlist <- list()
  brokenlist[[1]] <- NA
  while(!done) {
    for(i in 1:length(routes[[routetoinvestigate]])) {
      notrivers <- rivers
      notrivers$connections[routes[[routetoinvestigate]][i],] <- NA
      notrivers$connections[,routes[[routetoinvestigate]][i]] <- NA
      if(!is.na(brokenlist[[routetoinvestigate]][1])) {
        notrivers$connections[brokenlist[[routetoinvestigate]],] <- NA
        notrivers$connections[,brokenlist[[routetoinvestigate]]] <- NA
      }
      theroute <- detectroute(start=startseg,end=endseg,rivers=notrivers,algorithm="Dijkstra",stopiferror=F)
      if(!is.na(theroute[1])) {
        anyofthem <- F
        for(j in 1:length(routes)) {
          if(isTRUE(all.equal(routes[[j]],theroute))) anyofthem<-T
        }
        if(!anyofthem) {
          routes[[routetoadd]] <- theroute
          if(is.na(brokenlist[[routetoinvestigate]][1])) brokenlist[[routetoadd]] <- routes[[routetoinvestigate]][i]
          else brokenlist[[routetoadd]] <- c(brokenlist[[routetoinvestigate]],routes[[routetoinvestigate]][i])
          routetoadd <- routetoadd+1
        }
      }
    }
    routetoinvestigate <- routetoinvestigate+1
    done <- routetoinvestigate>=length(routes)
  }
  return(routes)
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
#' @param reps Deprecated.  Was the number of randomized reorderings to try.
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
