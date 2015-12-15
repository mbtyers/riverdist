#' Check Flow-Connectedness
#' @description Checks to see if two segments are flow-connected.  Called internally within \link{riverdirection} and \link{upstream}.
#' @param seg1 First input segment
#' @param seg2 Second input segment
#' @param rivers The river network object to use
#' @param stopiferror Whether or not to exit with an error if a route cannot be
#'   found.  If this is set to \code{FALSE} and a route cannot be found,
#'   the function will return \code{NA} in the appropriate entry.  Defaults to \code{TRUE}.  See \link{detectroute}.
#' @param algorithm Which route detection algorithm to use (\code{"Dijkstra"},
#'   \code{"sequential"}, or \code{"segroutes"}).  If left as \code{NULL} (the
#'   default), the function will automatically make a selection.  See
#'   \link{detectroute} for more details.
#' @return Logical \code{TRUE} if the two segments are flow-connected, \code{FALSE} if they are not
#' @note The river mouth must be specified (see \link{setmouth}).
#' @author Matt Tyers
#' @examples
#' data(Gulk)
#' plot(Gulk)
#'
#' Gulk <- setmouth(seg=1, vert=1, rivers=Gulk)
#'
#' isflowconnected(seg1=13, seg2=14, rivers=Gulk)
#' isflowconnected(seg1=13, seg2=1, rivers=Gulk)
#' @export
isflowconnected <- function(seg1,seg2,rivers,stopiferror=TRUE,algorithm=NULL) {
  connected <- FALSE
  if(is.na(rivers$mouth$mouth.seg)|is.na(rivers$mouth$mouth.vert)) stop("River mouth must be specified.")
  route1 <- detectroute(start=rivers$mouth$mouth.seg,end=seg1,rivers=rivers,stopiferror=stopiferror,algorithm=algorithm)
  route2 <- detectroute(start=rivers$mouth$mouth.seg,end=seg2,rivers=rivers,stopiferror=stopiferror,algorithm=algorithm)
  if(is.na(route1[1]) | is.na(route2[1])) connected <- NA
  if(!is.na(route1[1]) & !is.na(route2[1])) {
    if(any(route1==seg2)) connected <- TRUE
    if(any(route2==seg1)) connected <- TRUE
  }
  return(connected)
}

#' River Direction
#' @description Calculates direction of travel between two points.  Only works
#'   if river mouth (lowest point) has been specified (see \link{setmouth}).
#' @param startseg Segment number of the start of the route
#' @param endseg Segment number of the end of the route
#' @param startvert Vertex number of the start of the route
#' @param endvert Vertex number of the end of the route
#' @param rivers The river network object to use
#' @param stopiferror Whether or not to exit with an error if a route cannot be
#'   found.  If this is set to \code{FALSE} and a route cannot be found,
#'   the function will return \code{NA} in the appropriate entry.  Defaults to \code{TRUE}.  See \link{detectroute}.
#' @param algorithm Which route detection algorithm to use (\code{"Dijkstra"},
#'   \code{"sequential"}, or \code{"segroutes"}).  If left as \code{NULL} (the
#'   default), the function will automatically make a selection.  See
#'   \link{detectroute} for more details.
#' @param flowconnected If \code{TRUE}, only returns direction if the two input segments are flow-connected.  Defaults to \code{FALSE}.
#' @return Direction: "up", "down", or "0" (character).  Returns NA if \code{flowconnected==TRUE} and the two segments are not flow-connected.
#' @note Building routes from the river mouth to each river network segment may
#'   greatly reduce computation time (see \link{buildsegroutes}).
#' @author Matt Tyers
#' @examples
#' data(Gulk)
#'
#' # Mouth must be specified
#' Gulk$mouth$mouth.seg <- 1
#' Gulk$mouth$mouth.vert <- 1
#'
#' plot(x=Gulk)
#' riverdirection(startseg=6, endseg=3, startvert=40, endvert=40, rivers=Gulk)
#' @seealso \link{setmouth}
#' @export
riverdirection <- function(startseg,endseg,startvert,endvert,rivers,flowconnected=FALSE,stopiferror=TRUE,algorithm=NULL) {
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  if(max(c(startseg,endseg),na.rm=T)>length(rivers$lines) | min(c(startseg,endseg),na.rm=T)<1) {
    stop("Invalid segments specified.")
  }
  if(startvert>dim(rivers$lines[[startseg]])[1] | startvert<1 | endvert>dim(rivers$lines[[endseg]])[1] | endvert<1) {
    stop("Invalid vertex specified.")
  }
  
  if(is.na(rivers$mouth$mouth.seg) | is.na(rivers$mouth$mouth.vert)) {
    stop("Error - Need to specify segment & vertex of river mouth")
  }
  direction <- "0"
  flowc <- isflowconnected(seg1=startseg,seg2=endseg,rivers=rivers,stopiferror=stopiferror,algorithm=algorithm)
  if(!stopiferror & is.na(flowc)) direction <- NA
  if(!is.na(flowc)) {
    if(flowconnected & !flowc) direction <- NA
    if(!flowconnected | flowc) {
      if(riverdistance(startseg=rivers$mouth$mouth.seg,endseg=startseg,startvert=rivers$mouth$mouth.vert,endvert=startvert,
                       rivers=rivers,stopiferror=stopiferror,algorithm=algorithm) <
         riverdistance(startseg=rivers$mouth$mouth.seg,endseg=endseg,startvert=rivers$mouth$mouth.vert,endvert=endvert,
                       rivers=rivers,stopiferror=stopiferror,algorithm=algorithm)) {
        direction <- "up"
      }
      if(riverdistance(startseg=rivers$mouth$mouth.seg,endseg=startseg,startvert=rivers$mouth$mouth.vert,endvert=startvert,
                       rivers=rivers,stopiferror=stopiferror,algorithm=algorithm) >
         riverdistance(startseg=rivers$mouth$mouth.seg,endseg=endseg,startvert=rivers$mouth$mouth.vert,endvert=endvert,
                       rivers=rivers,stopiferror=stopiferror,algorithm=algorithm)) {
        direction <- "down"
      }
    }
  }
  return(direction)
}

#' Upstream River Distance
#' @description Calculates river network distances as +/-, defined as
#'   upriver/downriver.
#'  
#'   Specifying \code{net=TRUE} will compute net upriver distance (3 river km
#'   down a tributary and then 15 river km up the mainstem will mean 12 rkm net.
#'   Otherwise the function will return 18 rkm upriver travel.)
#'  
#'   The mouth (lowest point) segment and vertex must be specified (see
#'   \link{setmouth}).
#' @param startseg Segment number of the start of the route
#' @param endseg Segment number of the end of the route
#' @param startvert Vertex number of the start of the route
#' @param endvert Vertex number of the end of the route
#' @param rivers The river network object to use
#' @param flowconnected If \code{TRUE}, only returns distance if the two input segments are flow-connected.  Defaults to \code{FALSE}.
#' @param net Whether to calculate net distance (\code{net=TRUE}) or total
#'   distance (\code{net=FALSE})
#' @param stopiferror Whether or not to exit with an error if a route cannot be
#'   found.  If this is set to \code{FALSE} and a route cannot be found,
#'   the function will return \code{NA} in the appropriate entry.  Defaults to \code{TRUE}.  See \link{detectroute}.
#' @param algorithm Which route detection algorithm to use (\code{"Dijkstra"},
#'   \code{"sequential"}, or \code{"segroutes"}).  If left as \code{NULL} (the
#'   default), the function will automatically make a selection.  See
#'   \link{detectroute} for more details.
#' @return Upstream distance (numeric).  Returns NA if \code{flowconnected} has value \code{TRUE} and the two segments are not flow-connected.
#' @note Building routes from the river mouth to each river network segment may
#'   greatly reduce computation time (see \link{buildsegroutes}).
#' @author Matt Tyers
#' @examples
#' data(Gulk)
#'
#' # Mouth must be specified
#' Gulk$mouth$mouth.seg <- 1
#' Gulk$mouth$mouth.vert <- 1
#'
#' plot(x=Gulk)
#' riverpoints(seg=c(6,4), vert=c(140,140), pch=16, col=2, rivers=Gulk)
#' upstream(startseg=6, endseg=4, startvert=140, endvert=40, rivers=Gulk, net=TRUE)
#' upstream(startseg=6, endseg=4, startvert=140, endvert=40, rivers=Gulk, net=FALSE)
#' upstream(startseg=6, endseg=4, startvert=140, endvert=40, rivers=Gulk, flowconnected=TRUE)
#' @seealso \link{setmouth}
#' @export
upstream <- function(startseg,endseg,startvert,endvert,rivers,flowconnected=FALSE,net=FALSE,stopiferror=TRUE,algorithm=NULL) {
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  if(max(c(startseg,endseg),na.rm=T)>length(rivers$lines) | min(c(startseg,endseg),na.rm=T)<1) {
    stop("Invalid segments specified.")
  }
  if(startvert>dim(rivers$lines[[startseg]])[1] | startvert<1 | endvert>dim(rivers$lines[[endseg]])[1] | endvert<1) {
    stop("Invalid vertex specified.")
  }
  
  if(is.na(rivers$mouth$mouth.seg) | is.na(rivers$mouth$mouth.vert)) {
    stop("Error - Need to specify segment & vertex of origin")
  }
  flowc <- isflowconnected(seg1=startseg,seg2=endseg,rivers=rivers,stopiferror=stopiferror,algorithm=algorithm)
  if(!stopiferror & is.na(flowc)) distance <- NA 
  if(!is.na(flowc)) {
    if(flowconnected & !flowc) distance <- NA
    if(!flowconnected | flowc) {
      if(net) {
        distance <- riverdistance(startseg=rivers$mouth$mouth.seg,endseg=endseg,startvert=rivers$mouth$mouth.vert,endvert=endvert,
                                  rivers=rivers,stopiferror=stopiferror,algorithm=algorithm) -
          riverdistance(startseg=rivers$mouth$mouth.seg,endseg=startseg,startvert=rivers$mouth$mouth.vert,endvert=startvert,
                        rivers=rivers,stopiferror=stopiferror,algorithm=algorithm)
      }
      if(!net) {
        rawdist <- riverdistance(startseg=startseg,endseg=endseg,startvert=startvert,endvert=endvert,
                                 rivers=rivers,stopiferror=stopiferror,algorithm=algorithm)
        updown <- riverdirection(startseg=startseg,endseg=endseg,startvert=startvert,endvert=endvert,
                                 rivers=rivers,stopiferror=stopiferror,algorithm=algorithm)
        distance <- NA
        distance <- rawdist*(updown=="up") - rawdist*(updown=="down")
        if(rawdist==0) distance <- 0
      }
    }
  }
  return(distance) 
}


#' River Travel Direction Between Sequential Observations
#' @description Returns a matrix of directions traveled by unique fish between
#'   sequential surveys.  The mouth (lowest point) segment and vertex must be
#'   specified (see \link{setmouth}).
#' @param unique A vector of identifiers for each fish.
#' @param survey A vector of identifiers for each survey.  It is recommended to use a numeric or date format (see \link{as.Date}) to preserve survey order.
#' @param seg A vector of river locations (segment component).
#' @param vert A vector of river locations (vertex component).
#' @param rivers The river network object to use.
#' @param logical A boolean vector that can be used for subsetting - if used,
#'   \code{riverdirectionseq()} will only return pairwise distances in which a
#'   specified condition is met.
#' @param flowconnected If \code{TRUE}, only returns direction if the input segments are flow-connected.  Defaults to \code{FALSE}.
#' @param stopiferror Whether or not to exit with an error if a route cannot be
#'   found.  If this is set to \code{FALSE} and a route cannot be found,
#'   the function will return \code{NA} in the appropriate entry.  Defaults to \code{TRUE}.  See \link{detectroute}.
#' @param algorithm Which route detection algorithm to use (\code{"Dijkstra"},
#'   \code{"sequential"}, or \code{"segroutes"}).  If left as \code{NULL} (the
#'   default), the function will automatically make a selection.  See
#'   \link{detectroute} for more details.
#' @return A data frame of directions (character), with rows defined by unique
#'   fish and columns defined by observation increment (1 to 2, 2 to 3, etc.)  See \link{riverdirection} for additional information.
#' @seealso \link{riverdirection}
#' @note Building routes from the river mouth to each river network segment may
#'   greatly reduce computation time (see \link{buildsegroutes}).
#' @author Matt Tyers
#' @examples
#' data(Gulk, fakefish)
#'
#' # Mouth must be specified
#' Gulk$mouth$mouth.seg <- 1
#' Gulk$mouth$mouth.vert <- 1
#'
#' riverdirectionseq(unique=fakefish$fish.id, survey=fakefish$flight, seg=fakefish$seg,
#'    vert=fakefish$vert, rivers=Gulk)
#'
#' riverdirectionseq(unique=fakefish$fish.id, survey=fakefish$flight.date, seg=fakefish$seg,
#'    vert=fakefish$vert, rivers=Gulk)
#' @export
riverdirectionseq <- function(unique,survey,seg,vert,rivers,logical=NULL,flowconnected=FALSE,stopiferror=TRUE,algorithm=NULL) {
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  if(is.null(logical)) logical <- rep(T,length(unique))
  
  unique<-unique[logical]
  survey<-survey[logical]
  seg<-seg[logical]
  vert<-vert[logical]
  tab <- table(unique,survey)
  if(max(tab)>1) cat("Warning: multiple entries exist for at least one unique/survey combination (first one used)")
  dists <- matrix(NA,nrow=dim(tab)[1],ncol=(dim(tab)[2]-1))
  for(i in 1:dim(tab)[1]) {
    for(j in 1:(dim(tab)[2]-1)) {
      if(tab[i,j]*tab[i,(j+1)]!=0) {
        dists[i,j] <- riverdirection(startseg=seg[unique==sort(unique(unique))[i] & survey==sort(unique(survey))[j]][1],
                                     endseg=seg[unique==sort(unique(unique))[i] & survey==sort(unique(survey))[j+1]][1],
                                     startvert=vert[unique==sort(unique(unique))[i] & survey==sort(unique(survey))[j]][1],
                                     endvert=vert[unique==sort(unique(unique))[i] & survey==sort(unique(survey))[j+1]][1],
                                     rivers=rivers,flowconnected=flowconnected,stopiferror=stopiferror,algorithm=algorithm)
      }
    }
  }
  dists<-as.data.frame(dists)
  row.names(dists) <- row.names(tab)
  col.name<-NA
  for(j in 1:(length(dimnames(tab)$survey)-1)) col.name[j] <- paste(dimnames(tab)$survey[j],"to",dimnames(tab)$survey[j+1])
  names(dists) <- col.name
  dists <- dists[rowSums(is.na(dists)) != ncol(dists),]
  return(dists)
}


#' River Direction Matrix of All Observations of an Individual
#' @description Returns a matrix of travel direction between all observations of
#'   one unique fish.
#' @param indiv The unique identifier of the fish in question.
#' @param ID A vector of identifiers for each fish.
#' @param survey A vector of identifiers for each survey.  It is recommended to
#'   use a numeric or date format (see \link{as.Date}) to preserve survey order.
#' @param seg A vector of river locations (segment component).
#' @param vert A vector of river locations (vertex component).
#' @param rivers The river network object to use.
#' @param full Whether to return the full matrix, with \code{NA} values for
#'   missing data (\code{TRUE}), or a the subset of rows and columns
#'   corresponding to successful observations.  Defaults to \code{TRUE}.
#' @param flowconnected If \code{TRUE}, only returns direction if the input segments are flow-connected.  Defaults to \code{FALSE}.
#' @param stopiferror Whether or not to exit with an error if a route cannot be
#'   found.  If this is set to \code{FALSE} and a route cannot be found, the
#'   function will return \code{NA} in the appropriate entry.  Defaults to
#'   \code{TRUE}.  See \link{detectroute}.
#' @param algorithm Which route detection algorithm to use (\code{"Dijkstra"},
#'   \code{"sequential"}, or \code{"segroutes"}).  If left as \code{NULL} (the
#'   default), the function will automatically make a selection.  See
#'   \link{detectroute} for more details.
#' @return A matrix of directions (character), with rows and columns defined by
#'   survey.  In the resulting matrix, the element with the row identified as
#'   \code{A} and column identified as \code{B} is defined as the direction
#'   traveled from survey A to survey B.  Therefore, it is likely that only the
#'   upper triangle of the matrix will be of interest.
#' @seealso \link{riverdirection}
#' @note Building routes from the river mouth to each river network segment may
#'   greatly reduce computation time (see \link{buildsegroutes}).
#' @author Matt Tyers
#' @examples
#' data(Gulk, fakefish)
#' riverdirectionmatobs(indiv=1, ID=fakefish$fish.id, survey=fakefish$flight,
#'       seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk)
#'      
#' riverdirectionmatobs(indiv=1, ID=fakefish$fish.id, survey=fakefish$flight,
#'       seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk, full=FALSE)
#' @export
riverdirectionmatobs <- function(indiv,ID,survey,seg,vert,rivers,full=TRUE,flowconnected=FALSE,stopiferror=TRUE,algorithm=NULL) {
  surveys <- sort(unique(survey))
  surveys_indiv <- sort(unique(survey[ID==indiv]))
  # IDs <- order(unique(ID))
  
  outmat <- matrix(NA,nrow=length(surveys),ncol=length(surveys))
  for(ii in 1:length(surveys)) {
    for(jj in 1:length(surveys)) {
      outmat[ii,jj] <- ifelse((length(seg[ID==indiv & survey==surveys[ii]])==0) | (length(seg[ID==indiv & survey==surveys[jj]])==0),NA,
                              riverdirection(startseg=seg[ID==indiv & survey==surveys[ii]], endseg=seg[ID==indiv & survey==surveys[jj]],
                                             startvert=vert[ID==indiv & survey==surveys[ii]], endvert=vert[ID==indiv & survey==surveys[jj]],
                                             rivers=rivers,flowconnected=flowconnected,stopiferror=stopiferror,algorithm=algorithm))
    }
  }
  dimnames(outmat)[[1]] <- dimnames(outmat)[[2]] <- as.character(surveys)
  if(!full) {
    if(!all(is.na(outmat))) {
      whichnotna <- NA
      iwhichnotna <- 1
      for(i in 1:dim(outmat)[1]) {
        if(!all(is.na(outmat[,i]))) {
          whichnotna[iwhichnotna] <- i
          iwhichnotna <- iwhichnotna+1
        }
      }
      outmat <- outmat[whichnotna,whichnotna]
    }
    if(all(is.na(outmat))) outmat <- NA
  }
  return(outmat)
}

#' River Direction Matrix
#' @description Returns a matrix of calculated travel direction between every
#'   point and every other point of given river coordinates (segment and
#'   vertex), or of a subset.  The mouth (lowest point) segment and vertex must
#'   be specified (see \link{setmouth}).
#' @param seg A vector of river locations (segment component).
#' @param vert A vector of river locations (vertex component).
#' @param rivers The river network object to use
#' @param logical A boolean vector that can be used for subsetting - if used,
#'   \code{riverdirectionmat()} will only return pairwise distances in which a
#'   specified condition is met.
#' @param ID a vector of observation IDs for aid in interpreting the output
#'   table
#' @param flowconnected If \code{TRUE}, only returns direction if the input segments are flow-connected.  Defaults to \code{FALSE}.
#' @param stopiferror Whether or not to exit with an error if a route cannot be
#'   found.  If this is set to \code{FALSE} and a route cannot be found,
#'   the function will return \code{NA} in the appropriate entry.  Defaults to \code{TRUE}.  See \link{detectroute}.
#' @param algorithm Which route detection algorithm to use (\code{"Dijkstra"},
#'   \code{"sequential"}, or \code{"segroutes"}).  If left as \code{NULL} (the
#'   default), the function will automatically make a selection.  See
#'   \link{detectroute} for more details.
#' @return A matrix of directions (character) with rows and columns labeled by
#'   corresponding values of \code{ID}.  See \link{riverdirection} for additional information.
#' @seealso \link{riverdirection}
#' @note Building routes from the river mouth to each river network segment may
#'   greatly reduce computation time (see \link{buildsegroutes}).
#' @author Matt Tyers
#' @examples
#' data(Gulk, fakefish)
#'
#' # Mouth must be specified
#' Gulk$mouth$mouth.seg <- 1
#' Gulk$mouth$mouth.vert <- 1
#'
#' logi1 <- (fakefish$flight.date==as.Date("2015-11-25"))
#'
#' riverdirectionmat(seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk, logical=logi1)
#' @export
riverdirectionmat <- function(seg,vert,rivers,logical=NULL,ID=NULL,flowconnected=FALSE,stopiferror=TRUE,algorithm=NULL) {
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  if(is.null(logical)) logical <- rep(T,length(unique))
  
  len <- length(vert)
  seg<-seg[logical]
  vert<-vert[logical]
  if(is.null(ID)) ID <- 1:len
  ID <- ID[logical]
  dists <- matrix(NA,nrow=length(vert),ncol=length(vert))
  for(i in 1:length(vert)) {
    for(j in  1:length(vert)) {
      dists[i,j] <- riverdirection(startseg=seg[i],endseg=seg[j],startvert=vert[i],endvert=vert[j],
                                   rivers=rivers,flowconnected=flowconnected,stopiferror=stopiferror,algorithm=algorithm)
    }
  }
  dimnames(dists)[[1]] <- dimnames(dists)[[2]] <- ID
  return(dists)
}

#' Upstream Distance Between Sequential Observations
#' @description Returns a matrix of distance with direction by unique fish
#'   between sequential surveys.  The mouth (lowest point) segment and vertex
#'   must be specified (see \link{setmouth}).
#' @param unique A vector of identifiers for each fish.
#' @param survey A vector of identifiers for each survey.  It is recommended to use a numeric or date format (see \link{as.Date}) to preserve survey order.
#' @param seg A vector of river locations (segment component).
#' @param vert A vector of river locations (vertex component).
#' @param rivers The river network object to use.
#' @param logical A boolean vector that can be used for subsetting - if used,
#'   \code{upstreamseq()} will only return pairwise distances in which a
#'   specified condition is met.
#' @param flowconnected If \code{TRUE}, only returns distance if the input segments are flow-connected.  Defaults to \code{FALSE}.
#' @param net Whether to calculate net upstream distance (net=TRUE) or total
#'   distance (net=FALSE, default).
#' @param stopiferror Whether or not to exit with an error if a route cannot be
#'   found.  If this is set to \code{FALSE} and a route cannot be found,
#'   the function will return \code{NA} in the appropriate entry.  Defaults to \code{TRUE}.  See \link{detectroute}.
#' @param algorithm Which route detection algorithm to use (\code{"Dijkstra"},
#'   \code{"sequential"}, or \code{"segroutes"}).  If left as \code{NULL} (the
#'   default), the function will automatically make a selection.  See
#'   \link{detectroute} for more details.
#' @return A data frame of upstream distances (numeric), with rows defined by
#'   unique fish and columns defined by observation increment (1 to 2, 2 to 3,
#'   etc.)  See \link{upstream} for additional information.
#' @seealso \link{upstream}
#' @author Matt Tyers
#' @note Returns either net upstream distance (net=TRUE) or total distance
#'   (net=FALSE, default).  See \link{upstream}.
#' @note Building routes from the river mouth to each river network segment may
#'   greatly reduce computation time (see \link{buildsegroutes}).
#' @examples
#' data(Gulk, fakefish)
#'
#' # Mouth must be specified
#' Gulk$mouth$mouth.seg <- 1
#' Gulk$mouth$mouth.vert <- 1
#'
#' upstreamseq(unique=fakefish$fish.id, survey=fakefish$flight, seg=fakefish$seg,
#'       vert=fakefish$vert, rivers=Gulk)
#'
#' upstreamseq(unique=fakefish$fish.id, survey=fakefish$flight.date, seg=fakefish$seg,
#'       vert=fakefish$vert, rivers=Gulk)
#' @export
upstreamseq <- function(unique,survey,seg,vert,rivers,logical=NULL,flowconnected=FALSE,net=FALSE,stopiferror=TRUE,algorithm=NULL) {
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  if(is.null(logical)) logical <- rep(T,length(unique))
  
  unique<-unique[logical]
  survey<-survey[logical]
  seg<-seg[logical]
  vert<-vert[logical]
  tab <- table(unique,survey)
  if(max(tab)>1) cat("Warning: multiple entries exist for at least one unique/survey combination (first one used)")
  dists <- matrix(NA,nrow=dim(tab)[1],ncol=(dim(tab)[2]-1))
  for(i in 1:dim(tab)[1]) {
    for(j in 1:(dim(tab)[2]-1)) {
      if(tab[i,j]*tab[i,(j+1)]!=0) {
        dists[i,j] <- upstream(startseg=seg[unique==sort(unique(unique))[i] & survey==sort(unique(survey))[j]][1],
                               endseg=seg[unique==sort(unique(unique))[i] & survey==sort(unique(survey))[j+1]][1],
                               startvert=vert[unique==sort(unique(unique))[i] & survey==sort(unique(survey))[j]][1],
                               endvert=vert[unique==sort(unique(unique))[i] & survey==sort(unique(survey))[j+1]][1],
                               rivers=rivers,net=net,flowconnected=flowconnected,stopiferror=stopiferror,algorithm=algorithm)
      }
    }
  }
  dists<-as.data.frame(dists)
  row.names(dists) <- row.names(tab)
  col.name<-NA
  for(j in 1:(length(dimnames(tab)$survey)-1)) col.name[j] <- paste(dimnames(tab)$survey[j],"to",dimnames(tab)$survey[j+1])
  names(dists) <- col.name
  dists <- dists[rowSums(is.na(dists)) != ncol(dists),]
  return(dists)
}


#' Upstream Distance Matrix of All Observations of an Individual
#' @description Returns a matrix of upstream travel distance between all observations of
#'   one unique fish.
#' @param indiv The unique identifier of the fish in question.
#' @param ID A vector of identifiers for each fish.
#' @param survey A vector of identifiers for each survey.  It is recommended to
#'   use a numeric or date format (see \link{as.Date}) to preserve survey order.
#' @param seg A vector of river locations (segment component).
#' @param vert A vector of river locations (vertex component).
#' @param rivers The river network object to use.
#' @param full Whether to return the full matrix, with \code{NA} values for
#'   missing data (\code{TRUE}), or a the subset of rows and columns
#'   corresponding to successful observations.  Defaults to \code{TRUE}.
#' @param flowconnected If \code{TRUE}, only returns direction if the input segments are flow-connected.  Defaults to \code{FALSE}.
#' @param net Whether to calculate net upstream distance (net=TRUE) or total
#'   distance (net=FALSE, default).
#' @param stopiferror Whether or not to exit with an error if a route cannot be
#'   found.  If this is set to \code{FALSE} and a route cannot be found, the
#'   function will return \code{NA} in the appropriate entry.  Defaults to
#'   \code{TRUE}.  See \link{detectroute}.
#' @param algorithm Which route detection algorithm to use (\code{"Dijkstra"},
#'   \code{"sequential"}, or \code{"segroutes"}).  If left as \code{NULL} (the
#'   default), the function will automatically make a selection.  See
#'   \link{detectroute} for more details.
#' @return A matrix of upstream distances (numeric), with rows and columns defined by
#'   survey.  In the resulting matrix, the element with the row identified as
#'   \code{A} and column identified as \code{B} is defined as the upstream distance
#'   traveled from survey A to survey B.  Therefore, it is likely that only the
#'   upper triangle of the matrix will be of interest.
#' @seealso \link{upstream}
#' @note Building routes from the river mouth to each river network segment may
#'   greatly reduce computation time (see \link{buildsegroutes}).
#' @author Matt Tyers
#' @examples
#' data(Gulk, fakefish)
#' upstreammatobs(indiv=1, ID=fakefish$fish.id, survey=fakefish$flight,
#'       seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk)
#'      
#' upstreammatobs(indiv=1, ID=fakefish$fish.id, survey=fakefish$flight,
#'       seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk, full=FALSE)
#' @export
upstreammatobs <- function(indiv,ID,survey,seg,vert,rivers,full=TRUE,flowconnected=FALSE,net=FALSE,stopiferror=TRUE,algorithm=NULL) {
  surveys <- sort(unique(survey))
  surveys_indiv <- sort(unique(survey[ID==indiv]))
  # IDs <- order(unique(ID))
  
  outmat <- matrix(NA,nrow=length(surveys),ncol=length(surveys))
  for(ii in 1:length(surveys)) {
    for(jj in 1:length(surveys)) {
      outmat[ii,jj] <- ifelse((length(seg[ID==indiv & survey==surveys[ii]])==0) | (length(seg[ID==indiv & survey==surveys[jj]])==0),NA,
                              upstream(startseg=seg[ID==indiv & survey==surveys[ii]], endseg=seg[ID==indiv & survey==surveys[jj]],
                                       startvert=vert[ID==indiv & survey==surveys[ii]], endvert=vert[ID==indiv & survey==surveys[jj]],
                                       rivers=rivers,flowconnected=flowconnected,net=net,stopiferror=stopiferror,algorithm=algorithm))
    }
  }
  dimnames(outmat)[[1]] <- dimnames(outmat)[[2]] <- as.character(surveys)
  if(!full) {
    if(!all(is.na(outmat))) {
      whichnotna <- NA
      iwhichnotna <- 1
      for(i in 1:dim(outmat)[1]) {
        if(!all(is.na(outmat[,i]))) {
          whichnotna[iwhichnotna] <- i
          iwhichnotna <- iwhichnotna+1
        }
      }
      outmat <- outmat[whichnotna,whichnotna]
    }
    if(all(is.na(outmat))) outmat <- NA
  }
  return(outmat)
}

#' Upstream Distance Matrix
#' @description Returns a matrix of upstream distance between every point and
#'   every other point of given river coordinates (segment and vertex), or of a
#'   subset.  The mouth (lowest point) segment and vertex must be specified
#'   (see \link{setmouth}).
#' @param seg A vector of river locations (segment component).
#' @param vert A vector of river locations (vertex component).
#' @param rivers The river network object to use.
#' @param logical A boolean vector that can be used for subsetting - if used,
#'   riverdirectionseq() will only return pairwise distances in which a
#'   specified condition is met.
#' @param ID a vector of observation IDs for aid in interpreting the output
#'   table
#' @param flowconnected If \code{TRUE}, only returns distance if the input segments are flow-connected.  Defaults to \code{FALSE}.
#' @param net Whether to calculate net upstream distance (net=TRUE) or total
#'   distance (net=FALSE, default).  See \link{upstream}.
#' @param stopiferror Whether or not to exit with an error if a route cannot be
#'   found.  If this is set to \code{FALSE} and a route cannot be found,
#'   the function will return \code{NA} in the appropriate entry.  Defaults to \code{TRUE}.  See \link{detectroute}.
#' @param algorithm Which route detection algorithm to use (\code{"Dijkstra"},
#'   \code{"sequential"}, or \code{"segroutes"}).  If left as \code{NULL} (the
#'   default), the function will automatically make a selection.  See
#'   \link{detectroute} for more details.
#' @return A matrix of upstream distances (numeric) with rows and columns
#'   labeled by corresponding values of \code{ID}.  See \link{upstream} for additional information.
#' @seealso \link{upstream}
#' @note Building routes from the river mouth to each river network segment may
#'   greatly reduce computation time (see \link{buildsegroutes}).
#' @author Matt Tyers
#' @examples
#' data(Gulk, fakefish)
#'
#' # Mouth must be specified
#' Gulk$mouth$mouth.seg <- 1
#' Gulk$mouth$mouth.vert <- 1
#'
#' logi1 <- (fakefish$flight.date==as.Date("2015-11-25"))
#'
#' upstreammat(seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk, logical=logi1)
#' @export
upstreammat <- function(seg,vert,rivers,logical=NULL,ID=NULL,flowconnected=FALSE,net=FALSE,stopiferror=TRUE,algorithm=NULL) {
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  if(is.null(logical)) logical <- rep(T,length(unique))
  
  len <- length(vert)
  seg<-seg[logical]
  vert<-vert[logical]
  if(is.null(ID)) ID <- 1:len
  ID <- ID[logical]
  dists <- matrix(NA,nrow=length(vert),ncol=length(vert))
  for(i in 1:length(vert)) {
    for(j in  1:length(vert)) {
      dists[i,j] <- upstream(startseg=seg[i],endseg=seg[j],startvert=vert[i],endvert=vert[j],
                             rivers=rivers,net=net,flowconnected=flowconnected,stopiferror=stopiferror,algorithm=algorithm)
    }
  }
  dimnames(dists)[[1]] <- dimnames(dists)[[2]] <- ID
  return(dists)
}

#' Distance From Mouth
#' @description Calculates distance from a river location (given as segment and
#'   vertex) and the specified mouth of the river network.  The mouth must first
#'   be specified (see \link{setmouth}).
#' @param seg Segment of the point in question
#' @param vert Vertex of the point in question
#' @param rivers The river network object to use
#' @param stopiferror Whether or not to exit with an error if a route cannot be
#'   found.  If this is set to \code{FALSE} and a route cannot be found,
#'   the function will return \code{NA} in the appropriate entry.  Defaults to \code{TRUE}.  See \link{detectroute}.
#' @param algorithm Which route detection algorithm to use (\code{"Dijkstra"},
#'   \code{"sequential"}, or \code{"segroutes"}).  If left as \code{NULL} (the
#'   default), the function will automatically make a selection.  See
#'   \link{detectroute} for more details.
#' @return Distance (numeric)
#' @note Building routes from the river mouth to each river network segment may greatly reduce computation time (see \link{buildsegroutes}).
#' @author Matt Tyers
#' @examples
#' data(Gulk)
#'
#' # Mouth must be specified
#' Gulk$mouth$mouth.seg <- 1
#' Gulk$mouth$mouth.vert <- 1
#'
#' mouthdist(seg=4, vert=40, rivers=Gulk)
#' @export
mouthdist <- function(seg,vert,rivers,stopiferror=TRUE,algorithm=NULL) {
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  if(seg>length(rivers$lines) | seg<1) {
    stop("Invalid segments specified.")
  }
  if(vert>dim(rivers$lines[[seg]])[1] | vert<1) {
    stop("Invalid vertex specified.")
  }
  
  if(is.na(rivers$mouth$mouth.seg) | is.na(rivers$mouth$mouth.vert)) {
    stop("Error - Need to specify segment & vertex of origin",'\n')
  }
  dists <- rep(NA,length(seg))
  for(i in 1:length(seg)) {
    dists[i] <- riverdistance(startseg=seg[i],endseg=rivers$mouth$mouth.seg,startvert=vert[i],endvert=rivers$mouth$mouth.vert,
                              rivers=rivers,stopiferror=stopiferror,algorithm=algorithm)
  }
  return(dists)
}


#' Distance From Mouth for All Observations of Individuals
#' @description Calculates distance from the mouth of a river network to all 
#'   observations of each individual (given as segment and vertex). and the 
#'   specified mouth of the river network.  The mouth must first be specified 
#'   (see \link{setmouth}).  Returns a matrix of distances, with a row for each 
#'   unique individual and a column for each survey.
#'   
#'   Plotting a matrix returned from \code{mouthdistobs()} can be done using \link{plotseq}. 
#' @param unique A vector of identifiers for each fish.
#' @param survey A vector of identifiers for each survey.  It is recommended to
#'   use a numeric or date format (see \link{as.Date}) to preserve survey order.
#' @param seg A vector of river coordinates (segment)
#' @param vert A vectpr pf rover coordinates (vertex)
#' @param rivers The river network object to use
#' @param logical A boolean vector that can be used for subsetting - if used, 
#'   \code{mouthdistobs()} will only return distances in which a specified
#'   condition is met.
#' @param stopiferror Whether or not to exit with an error if a route cannot be 
#'   found.  If this is set to \code{FALSE} and a route cannot be found, the 
#'   function will return \code{NA} in the appropriate entry.  Defaults to 
#'   \code{TRUE}.  See \link{detectroute}.
#' @param algorithm Which route detection algorithm to use (\code{"Dijkstra"}, 
#'   \code{"sequential"}, or \code{"segroutes"}).  If left as \code{NULL} (the 
#'   default), the function will automatically make a selection.  See 
#'   \link{detectroute} for more details.
#' @return A vector of river network distances (numeric), with each row
#'   corresponding to a unique fish and each column corresponding to a unique
#'   survey.  Values of \code{NA} indicate the individual not being located
#'   during the survey in question.
#' @note Building routes from the river mouth to each river network segment may 
#'   greatly reduce computation time (see \link{buildsegroutes}).
#' @seealso \link{plotseq}
#' @author Matt Tyers
#' @examples
#' data(Gulk, fakefish)
#' 
#' mouthdistobs(unique=fakefish$fish.id, survey=fakefish$flight.date, 
#'     seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk)
#' @export
mouthdistobs <- function(unique,survey,seg,vert,rivers,logical=NULL,stopiferror=TRUE,algorithm=NULL) {
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  if(is.null(logical)) logical <- rep(T,length(unique))
  
  unique<-unique[logical]
  survey<-survey[logical]
  seg<-seg[logical]
  vert<-vert[logical]
  tab <- table(unique,survey)
  if(max(tab)>1) cat("Warning: multiple entries exist for at least one unique/survey combination (first one used)")
  dists <- matrix(NA,nrow=dim(tab)[1],ncol=(dim(tab)[2]))
  for(i in 1:dim(tab)[1]) {
    for(j in 1:(dim(tab)[2])) {
      if(tab[i,j]!=0) {
        dists[i,j] <- mouthdist(seg=seg[unique==sort(unique(unique))[i] & survey==sort(unique(survey))[j]][1],
                                vert=vert[unique==sort(unique(unique))[i] & survey==sort(unique(survey))[j]][1],
                                rivers=rivers,stopiferror=stopiferror,algorithm=algorithm)
      }
    }
  }
  dists<-as.data.frame(dists)
  row.names(dists) <- row.names(tab)
  col.name<-NA
  names(dists) <- dimnames(tab)$survey
  dists <- dists[rowSums(is.na(dists)) != ncol(dists),]
  return(dists)
}


#' Plot Sequence of Observations
#' @description Plots the sequence of observations of each individual (given as 
#'   segment and vertex).  This function is primarily intended for use with 
#'   \link{mouthdistobs}, but will also work with \link{riverdistanceseq} and 
#'   \link{upstreamseq}.
#' @param seqobs A matrix returned from \link{mouthdistobs}, 
#'   \link{riverdistanceseq}, or \link{upstreamseq}.
#' @param type The type of plot to generate.  Options are 
#'   \code{"boxplot"},\code{"dotplot"},\code{"boxline"},or \code{"dotline"}. 
#'   Defaults to \code{"boxplot"}.
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param main Plot title
#' @param cex.axisX Character expansion factor for X-axis labels
#' @param lowerbound An optional vector of lower survey bounds
#' @param upperbound An optional vector of upper survey bounds
#' @param boundtype Method of plotting survey bounds.  Options are
#'   \code{"positive"}, \code{"negative"} (default), and \code{"lines"}.
#' @param surveysareDates If surveys are in Date format (see \link{as.Date}), a
#'   value of \code{TRUE} allows the x-coordinates points to be spaced apart
#'   according to date, not equidistantly.  Defaults to \code{FALSE}.  Any formatting of 
#'   the survey variable must be done within the original call to \link{mouthdistobs}, 
#'   \link{riverdistanceseq}, or \link{upstreamseq}.  Dates must already be formatted as dates,
#'   or in the form \code{"YYYY-MM-DD"} or \code{"YYYY/MM/DD"}.
#' @param ... Additional plotting parameters
#' @note Plots are intended as descriptive only.  Any ANOVA-like inference that 
#'   is suggested from these plots is strongly discouraged.  The user is instead
#'   advised to use a mixed-effects model or some other inferential tool that 
#'   accounts for repeated-measures and/or temporal autocorrelation.
#' @author Matt Tyers
#' @importFrom graphics polygon
#' @importFrom graphics boxplot
#' @importFrom graphics axis
#' @examples
#' data(Gulk, fakefish)
#' 
#' x <- mouthdistobs(unique=fakefish$fish.id, survey=fakefish$flight.date, 
#'     seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk)
#'     
#' # this is all randomly generated data and will not show much pattern
#' plotseq(seqobs=x)
#' plotseq(seqobs=x, type="boxline")
#' plotseq(seqobs=x, type="dotplot")
#' plotseq(seqobs=x, type="dotline")
#' 
#' plotseq(seqobs=x, type="dotline",surveysareDates=TRUE)
#' 
#' from_upstreamseq <- upstreamseq(unique=fakefish$fish.id, 
#'    survey=fakefish$flight, seg=fakefish$seg, vert=fakefish$vert, 
#'    rivers=Gulk)
#' plotseq(seqobs=from_upstreamseq)
#' @export
plotseq <- function(seqobs,type="boxplot",xlab="",ylab="",main="",cex.axisX=.8,lowerbound=NULL,upperbound=NULL,boundtype="negative",surveysareDates=F,...) {
  if(surveysareDates) xplot <- as.Date(names(seqobs))
  if(!surveysareDates) xplot <- 1:(dim(seqobs)[2])
  if(is.numeric(seqobs[1,1])) {
    plot(NA,xlim=c(xplot[1],xplot[length(xplot)]),ylim=c(min(seqobs,na.rm=T),max(seqobs,na.rm=T)),xaxt='n',xlab=xlab,ylab=ylab,...=...)
    if(!is.null(lowerbound)&!is.null(upperbound)) {
      if(boundtype=="negative") {
        polygon(x=c(xplot[1],xplot,xplot[length(xplot)]),y=c(par("usr")[3],lowerbound,par("usr")[3]),col="grey90",border=NA)
        polygon(x=c(xplot[1],xplot,xplot[length(xplot)]),y=c(par("usr")[4],upperbound,par("usr")[4]),col="grey90",border=NA)
        lines(par("usr")[1:2],par("usr")[c(4,4)])
      }
      if(boundtype=="positive") {
        polygon(x=c(xplot,xplot[(length(xplot)):1]),y=c(lowerbound,(upperbound[(length(upperbound)):1])),col="grey90",border=NA)
      }
      if(boundtype=="lines") {
        del <- .4*min(xplot[2:length(xplot)] - xplot[1:(length(xplot)-1)])
        for(i in 1:length(lowerbound)) {
          lines(xplot[i]+c(-1,1)*del,rep(lowerbound[i],2),lwd=2)
          lines(xplot[i]+c(-1,1)*del,rep(upperbound[i],2),lwd=2)
        }
      }
    }
    if(type=="dotline" | type=="boxline") {
      for(i in 1:(dim(seqobs)[1])) {
        lines(xplot[!is.na(seqobs[i,])],seqobs[i,][!is.na(seqobs[i,])],col="grey60",lty=3)
        lines(xplot,seqobs[i,],col="grey30")
      }
    }
    for(i in 1:(dim(seqobs)[2])) {
      if(type=="boxplot" | type=="boxline") boxplot(seqobs[,i],at=xplot[i],add=T,yaxt='n',col=NA)
      if(type=="dotplot") points(jitter(rep(xplot[i],(dim(seqobs)[1])),amount=.1),seqobs[,i])
      if(type=="dotline") points(rep(xplot[i],(dim(seqobs)[1])),seqobs[,i])
    }
    axis(side=1,at=xplot,labels=names(seqobs),cex.axis=cex.axisX,las=2)
  }
  if(is.character(seqobs[1,1])|is.factor(seqobs[1,1])) {
    stop("Plotting methods do not yet exist for matrices returned from riverdirectionseq().")
#     colname <- NULL
#     contents <- NULL
#     for(i in 1:(dim(seqobs)[2])) {
#       colname <- c(colname,rep(names(seqobs)[i],(dim(seqobs)[1])))
#       contents <- c(contents,seqobs[,i])
#     }
#     contents <- as.character(contents)
#     contents[contents=="down"] <- "3 down"
#     contents[contents=="0"] <- "2 0"
#     contents[contents=="up"] <- "1 up"
#     contents <- as.factor(contents)
#     plot(table(colname,contents))
  }
}


#' River Direction Matrix between Two Datasets
#' @description Returns a matrix of directions between each river location in two datasets, with one expressed as rows and the other expressed as columns.
#' @param seg1 First vector of river locations (segment component).  These are expressed as rows in the output matrix.
#' @param vert1 First vector of river locations (vertex component).  These are expressed as rows in the output matrix.
#' @param seg2 Second vector of river locations (segment component).  These are expressed as columns in the output matrix.
#' @param vert2 Second vector of river locations (vertex component).  These are expressed as columns in the output matrix.
#' @param rivers The river network object to use.
#' @param logical1 A boolean vector that can be used for subsetting.  If used,
#'   \code{riverdirectiontofrom} will only return directions in which a
#'   specified condition is met for the first dataset.
#' @param logical2 A boolean vector that can be used for subsetting.  If used,
#'   \code{riverdirectiontofrom} will only return directions in which a
#'   specified condition is met for the second dataset.
#' @param ID1 a vector of observation IDs for the first dataset that will be used as row names in the output matrix.
#' @param ID2 a vector of observation IDs for the second dataset that will be used as column names in the output matrix.
#' @param flowconnected If \code{TRUE}, only returns distance if the input segments are flow-connected.  Defaults to \code{FALSE}.
#' @param stopiferror Whether or not to exit with an error if a route cannot be
#'   found.  If this is set to \code{FALSE} and a route cannot be found,
#'   the function will return \code{NA} in the appropriate entry.  Defaults to \code{TRUE}.  See \link{detectroute}.
#' @param algorithm Which route detection algorithm to use (\code{"Dijkstra"},
#'   \code{"sequential"}, or \code{"segroutes"}).  If left as \code{NULL} (the
#'   default), the function will automatically make a selection.  See
#'   \link{detectroute} for more details.
#' @return A matrix of directions (character) with rows and columns labeled by corresponding values of \code{ID}.  See \link{riverdirection} for additional information.
#' @seealso \link{riverdirection}
#' @note Building routes from the river mouth to each river network segment may greatly reduce computation time (see \link{buildsegroutes}).
#' @author Matt Tyers
#' @examples
#' data(Gulk)
#'
#' streamlocs.seg <- c(1,8,11)
#' streamlocs.vert <- c(50,70,90)
#' streamlocs.ID <- c("A","B","C")
#'
#' fish.seg <- c(1,4,9,12,14)
#' fish.vert <- c(10,11,12,13,14)
#' fish.ID <- c("fish1","fish2","fish3","fish4","fish5")
#'
#' Gulk <- setmouth(seg=1, vert=1, rivers=Gulk)
#'
#' riverdirectiontofrom(seg1=streamlocs.seg, vert1=streamlocs.vert,
#'   seg2=fish.seg, vert2=fish.vert, rivers=Gulk,
#'   ID1=streamlocs.ID, ID2=fish.ID)
#'
#' logi1 <- streamlocs.ID=="B" | streamlocs.ID=="C"
#' logi2 <- fish.ID!="fish3"
#'
#' riverdirectiontofrom(seg1=streamlocs.seg, vert1=streamlocs.vert,
#'   seg2=fish.seg, vert2=fish.vert, rivers=Gulk, logical1=logi1,
#'   logical2=logi2, ID1=streamlocs.ID, ID2=fish.ID)
#' @export
riverdirectiontofrom <- function(seg1,vert1,seg2,vert2,rivers,logical1=NULL,logical2=NULL,ID1=NULL,ID2=NULL,flowconnected=FALSE,stopiferror=TRUE,algorithm=NULL) {
  if(is.null(logical1)) logical1 <- rep(T,length(seg1))
  if(is.null(logical2)) logical2 <- rep(T,length(seg2))
  
  if(length(logical1) != length(seg1)) stop("logical1 must be the same length as its location vectors")
  if(length(logical2) != length(seg2)) stop("logical2 must be the same length as its location vectors")
  
  if(is.null(ID1)) ID1 <- 1:length(seg1)
  if(is.null(ID2)) ID2 <- 1:length(seg2)
  
  seg1 <- seg1[logical1]
  vert1 <- vert1[logical1]
  seg2 <- seg2[logical2]
  vert2 <- vert2[logical2]
  ID1 <- ID1[logical1]
  ID2 <- ID2[logical2]
  
  dists <- matrix(NA,nrow=length(seg1),ncol=length(seg2))
  
  for(i in 1:length(seg1)) {
    for(j in 1:length(seg2)) {
      dists[i,j] <- riverdirection(startseg=seg1[i],startvert=vert1[i],endseg=seg2[j],endvert=vert2[j],
                                   flowconnected=flowconnected,rivers=rivers,stopiferror=stopiferror,algorithm=algorithm)
    }
  }
  
  class(ID1) <- "list"
  class(ID2) <- "list"
  dimnames(dists)[[1]] <- ID1
  dimnames(dists)[[2]] <- ID2
  return(dists)
}


#' Upstream Distance Matrix between Two Datasets
#' @description Returns a matrix of upstream distances between each river location in two datasets, with one expressed as rows and the other expressed as columns.
#' @param seg1 First vector of river locations (segment component).  These are expressed as rows in the output matrix.
#' @param vert1 First vector of river locations (vertex component).  These are expressed as rows in the output matrix.
#' @param seg2 Second vector of river locations (segment component).  These are expressed as columns in the output matrix.
#' @param vert2 Second vector of river locations (vertex component).  These are expressed as columns in the output matrix.
#' @param rivers The river network object to use.
#' @param logical1 A boolean vector that can be used for subsetting.  If used,
#'   \code{upstreamtofrom} will only return upstream distances in which a
#'   specified condition is met for the first dataset.
#' @param logical2 A boolean vector that can be used for subsetting.  If used,
#'   \code{upstreamtofrom} will only return upstream distances in which a
#'   specified condition is met for the second dataset.
#' @param ID1 a vector of observation IDs for the first dataset that will be used as row names in the output matrix.
#' @param ID2 a vector of observation IDs for the second dataset that will be used as column names in the output matrix.
#' @param net Whether to calculate net upstream distance (\code{TRUE}) or signed total distance (\code{FALSE}).  See \link{upstream}.
#' @param flowconnected If \code{TRUE}, only returns distance if the input segments are flow-connected.  Defaults to \code{FALSE}.
#' @param stopiferror Whether or not to exit with an error if a route cannot be
#'   found.  If this is set to \code{FALSE} and a route cannot be found,
#'   the function will return \code{NA} in the appropriate entry.  Defaults to \code{TRUE}.  See \link{detectroute}.
#' @param algorithm Which route detection algorithm to use (\code{"Dijkstra"},
#'   \code{"sequential"}, or \code{"segroutes"}).  If left as \code{NULL} (the
#'   default), the function will automatically make a selection.  See
#'   \link{detectroute} for more details.
#' @return A matrix of upstream distances (numeric) with rows and columns labeled by corresponding values of \code{ID}.  See \link{upstream} for additional information.
#' @seealso \link{upstream}
#' @note Building routes from the river mouth to each river network segment may greatly reduce computation time (see \link{buildsegroutes}).
#' @author Matt Tyers
#' @examples
#' data(Gulk)
#'
#' streamlocs.seg <- c(1,8,11)
#' streamlocs.vert <- c(50,70,90)
#' streamlocs.ID <- c("A","B","C")
#'
#' fish.seg <- c(1,4,9,12,14)
#' fish.vert <- c(10,11,12,13,14)
#' fish.ID <- c("fish1","fish2","fish3","fish4","fish5")
#'
#' Gulk <- setmouth(seg=1, vert=1, rivers=Gulk)
#'
#' upstreamtofrom(seg1=streamlocs.seg, vert1=streamlocs.vert,
#'   seg2=fish.seg, vert2=fish.vert, rivers=Gulk,
#'   ID1=streamlocs.ID, ID2=fish.ID)
#'
#' logi1 <- streamlocs.ID=="B" | streamlocs.ID=="C"
#' logi2 <- fish.ID!="fish3"
#'
#' upstreamtofrom(seg1=streamlocs.seg, vert1=streamlocs.vert,
#'   seg2=fish.seg, vert2=fish.vert, rivers=Gulk, logical1=logi1,
#'   logical2=logi2, ID1=streamlocs.ID, ID2=fish.ID)
#' @export
upstreamtofrom <- function(seg1,vert1,seg2,vert2,rivers,logical1=NULL,logical2=NULL,ID1=NULL,ID2=NULL,net=FALSE,flowconnected=FALSE,stopiferror=TRUE,algorithm=NULL) {
  if(is.null(logical1)) logical1 <- rep(T,length(seg1))
  if(is.null(logical2)) logical2 <- rep(T,length(seg2))
  
  if(length(logical1) != length(seg1)) stop("logical1 must be the same length as its location vectors")
  if(length(logical2) != length(seg2)) stop("logical2 must be the same length as its location vectors")
  
  if(is.null(ID1)) ID1 <- 1:length(seg1)
  if(is.null(ID2)) ID2 <- 1:length(seg2)
  
  seg1 <- seg1[logical1]
  vert1 <- vert1[logical1]
  seg2 <- seg2[logical2]
  vert2 <- vert2[logical2]
  ID1 <- ID1[logical1]
  ID2 <- ID2[logical2]
  
  dists <- matrix(NA,nrow=length(seg1),ncol=length(seg2))
  
  for(i in 1:length(seg1)) {
    for(j in 1:length(seg2)) {
      dists[i,j] <- upstream(startseg=seg1[i],startvert=vert1[i],endseg=seg2[j],endvert=vert2[j],
                             flowconnected=flowconnected,rivers=rivers,net=net,stopiferror=stopiferror,algorithm=algorithm)
    }
  }
  
  class(ID1) <- "list"
  class(ID2) <- "list"
  dimnames(dists)[[1]] <- ID1
  dimnames(dists)[[2]] <- ID2
  return(dists)
}