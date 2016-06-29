#' River Distance Between Sequential Observations
#' @description Returns a matrix of distances traveled by unique fish, between
#'   sequential surveys.  A plotting method is also provided for the output; see \link{plotseq}
#' @param unique A vector of identifiers for each fish.
#' @param survey A vector of identifiers for each survey.  It is recommended to use a numeric or date format (see \link{as.Date}) to preserve survey order.
#' @param seg A vector of river locations (segment component).
#' @param vert A vector of river locations (vertex component).
#' @param rivers The river network object to use.
#' @param logical A boolean vector that can be used for subsetting.  If used,
#'   \code{riverdistanceseq()} will only return pairwise distances in which a
#'   specified condition is met.
#' @param stopiferror Whether or not to exit with an error if a route cannot be
#'   found.  If this is set to \code{FALSE} and a route cannot be found,
#'   the function will return \code{NA} in the appropriate entry.  Defaults to \code{TRUE}.  See \link{detectroute}.
#' @param algorithm Which route detection algorithm to use (\code{"Dijkstra"},
#'   \code{"sequential"}, or \code{"segroutes"}).  If left as \code{NULL} (the
#'   default), the function will automatically make a selection.  See
#'   \link{detectroute} for more details.
#' @return A data frame of distances (numeric), with rows defined by unique fish and columns defined by observation increment (1 to 2, 2 to 3, etc.)
#' @seealso \link{riverdistance}, \link{plotseq}
#' @note Building routes from the river mouth to each river network segment and/or distance lookup tables will
#'   greatly reduce computation time (see \link{buildsegroutes}).
#' @author Matt Tyers
#' @examples
#' data(Gulk, fakefish)
#' riverdistanceseq(unique=fakefish$fish.id, survey=fakefish$flight,
#'       seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk)
#'      
#' seqbysurvey <- riverdistanceseq(unique=fakefish$fish.id, survey=fakefish$flight.date,
#'       seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk)
#' seqbysurvey
#' plotseq(seqbysurvey)
#' @export
riverdistanceseq <- function(unique,survey,seg,vert,rivers,logical=NULL,stopiferror=TRUE,algorithm=NULL) {
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
        dists[i,j] <- riverdistance(startseg=seg[unique==sort(unique(unique))[i] & survey==sort(unique(survey))[j]][1],
                                    endseg=seg[unique==sort(unique(unique))[i] & survey==sort(unique(survey))[j+1]][1],
                                    startvert=vert[unique==sort(unique(unique))[i] & survey==sort(unique(survey))[j]][1],
                                    endvert=vert[unique==sort(unique(unique))[i] & survey==sort(unique(survey))[j+1]][1],
                                    rivers=rivers,stopiferror=stopiferror,algorithm=algorithm)
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

#' River Distance Matrix of All Observations of an Individual
#' @description Returns a matrix of network distances between all observations
#'   of one unique fish.
#' @param indiv The unique identifier of the fish in question.
#' @param unique A vector of identifiers for each fish.
#' @param survey A vector of identifiers for each survey.  It is recommended to
#'   use a numeric or date format (see \link{as.Date}) to preserve survey order.
#' @param seg A vector of river locations (segment component).
#' @param vert A vector of river locations (vertex component).
#' @param rivers The river network object to use.
#' @param full Whether to return the full matrix, with \code{NA} values for
#'   missing data (\code{TRUE}), or a the subset of rows and columns
#'   corresponding to successful observations.  Defaults to \code{TRUE}.
#' @param stopiferror Whether or not to exit with an error if a route cannot be
#'   found.  If this is set to \code{FALSE} and a route cannot be found, the
#'   function will return \code{NA} in the appropriate entry.  Defaults to
#'   \code{TRUE}.  See \link{detectroute}.
#' @param algorithm Which route detection algorithm to use (\code{"Dijkstra"},
#'   \code{"sequential"}, or \code{"segroutes"}).  If left as \code{NULL} (the
#'   default), the function will automatically make a selection.  See
#'   \link{detectroute} for more details.
#' @return A matrix of distances (numeric), with rows and columns defined by
#'   survey.
#' @seealso \link{riverdistance}
#' @note Building routes from the river mouth to each river network segment and/or distance lookup tables will
#'   greatly reduce computation time (see \link{buildsegroutes}).
#' @author Matt Tyers
#' @examples
#' data(Gulk, fakefish)
#' riverdistancematbysurvey(indiv=1, unique=fakefish$fish.id, survey=fakefish$flight,
#'       seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk)
#'      
#' riverdistancematbysurvey(indiv=1, unique=fakefish$fish.id, survey=fakefish$flight,
#'       seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk, full=FALSE)
#' @export
riverdistancematbysurvey <- function(indiv,unique,survey,seg,vert,rivers,full=TRUE,stopiferror=TRUE,algorithm=NULL) {
  surveys <- sort(unique(survey))
  surveys_indiv <- sort(unique(survey[unique==indiv]))
  
  outmat <- matrix(NA,nrow=length(surveys),ncol=length(surveys))
  for(ii in 1:length(surveys)) {
    for(jj in 1:length(surveys)) {
      outmat[ii,jj] <- ifelse((length(seg[unique==indiv & survey==surveys[ii]])==0) | (length(seg[unique==indiv & survey==surveys[jj]])==0),NA,
                              riverdistance(startseg=seg[unique==indiv & survey==surveys[ii]], endseg=seg[unique==indiv & survey==surveys[jj]],
                                            startvert=vert[unique==indiv & survey==surveys[ii]], endvert=vert[unique==indiv & survey==surveys[jj]],
                                            rivers=rivers,stopiferror=stopiferror,algorithm=algorithm))
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

#' River Distance Matrix
#' @description Returns a matrix of distances between every point and every
#'   other point of given river locations (segment and vertex), or of a subset.
#' @param seg A vector of river locations (segment component).
#' @param vert A vector of river locations (vertex component).
#' @param rivers The river network object to use.
#' @param logical A boolean vector that can be used for subsetting.  If used,
#'   \code{riverdistancemat} will only return pairwise distances in which a
#'   specified condition is met.
#' @param ID a vector of observation IDs for aid in interpreting the output
#'   table
#' @param stopiferror Whether or not to exit with an error if a route cannot be
#'   found.  If this is set to \code{FALSE} and a route cannot be found,
#'   the function will return \code{NA} in the appropriate entry.  Defaults to \code{TRUE}.  See \link{detectroute}
#' @param algorithm Which route detection algorithm to use (\code{"Dijkstra"},
#'   \code{"sequential"}, or \code{"segroutes"}).  If left as \code{NULL} (the
#'   default), the function will automatically make a selection.  See
#'   \link{detectroute} for more details.
#' @return A matrix of distances (numeric) with rows and columns labeled by corresponding values of \code{ID}.
#' @seealso \link{riverdistance}
#' @note Building routes from the river mouth to each river network segment and/or distance lookup tables will
#'   greatly reduce computation time (see \link{buildsegroutes}).
#' @author Matt Tyers
#' @examples
#' data(Gulk, fakefish)
#'
#' logi1 <- (fakefish$flight.date==as.Date("2015-11-25"))
#'
#' riverdistancemat(seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk, logical=logi1)
#' @export
riverdistancemat <- function(seg,vert,rivers,logical=NULL,ID=NULL,stopiferror=TRUE,algorithm=NULL) {
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
      dists[i,j] <- riverdistance(startseg=seg[i],endseg=seg[j],startvert=vert[i],endvert=vert[j],rivers=rivers,
                                  stopiferror=stopiferror,algorithm=algorithm)
    }
  }
  dimnames(dists)[[1]] <- dimnames(dists)[[2]] <- ID
  return(dists)
}


#' Home Range
#' @description Returns the minimum observed home range for multiple
#'   observations of each individual fish.
#' @param unique A vector of unique identifiers for each fish.
#' @param seg A vector of river locations (segment component).
#' @param vert A vector of river locations (vertex component).
#' @param rivers The river network object to use.
#' @param map Boolean (defaults to FALSE) Whether to produce sanity-check maps
#'   of observed locations and calculated home range for each fish.
#' @return A data frame with two columns: \code{$ID} is a list of unique fish
#'   (as specified by \code{unique=}), and \code{$range} is calculated minimum
#'   home range, in the units of the coordinate system (this will likely be
#'   meters).
#' @param algorithm Which route detection algorithm to use (\code{"Dijkstra"},
#'   \code{"sequential"}, or \code{"segroutes"}).  If left as \code{NULL} (the
#'   default), the function will automatically make a selection.  See
#'   \link{detectroute} for more details.
#' @param main Plot title, if \code{map} is set to \code{TRUE}.  If unspecified, the unique ID will be used for the title.
#' @param ... Additional plotting arguments, if \code{map} is set to \code{TRUE}.
#' @note Building routes from the river mouth to each river network segment and/or distance lookup tables will
#'   greatly reduce computation time (see \link{buildsegroutes}).
#' @author Matt Tyers
#' @examples
#' data(Gulk, fakefish)
#' homerange(unique=fakefish$fish.id, seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk)
#'
#' # mapping shown just for fish 15
#' homerange(unique=fakefish$fish.id[fakefish$fish.id==15], seg=fakefish$seg[fakefish$fish.id==15],
#'           vert=fakefish$vert[fakefish$fish.id==15], rivers=Gulk, map=TRUE)
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @export
homerange <- function(unique,seg,vert,rivers,map=FALSE,algorithm=NULL,main=NULL,...) {
  if(class(rivers)!="rivernetwork") stop("Argument 'rivers' must be of class 'rivernetwork'.  See help(line2network) for more information.")
  if((length(unique)!=length(seg))|(length(seg)!=length(vert))) stop("Input vectors must be the same length.")
  ID <- sort(unique(unique))
  range <- rep(0,length(ID))
  
  alllengths <- function(xy) {       #############
    n <- dim(xy)[1]
    if(n==1) dist <- 0
    if(n==2) dist <- pdist(xy[1,],xy[2,])
    if(n>2) dist <- sqrt(((xy[1:(n-1),1] - xy[2:n,1])^2) + ((xy[1:(n-1),2] - xy[2:n,2])^2))
    return(dist)
  }
  subseglength <- list()
  for(i in 1:length(rivers$lines)) {
    subseglength[[i]] <- alllengths(rivers$lines[[i]])
  }   ##################
      
  for(i in 1:length(ID)) {
    # i <- 7
    if(map==T) {
      if(is.null(main)) {
        main1 <- ID[i]
      } else{
        main1 <- main
      }
      plot(x=rivers,main=main1,color=F,segmentnum=F)#,...=...)  ###############
    }
    
    n.entries <- length(unique[unique==ID[i]])
    if(n.entries>1) {
      # create a list of routes taken for fish i
      seg1 <- seg[unique==ID[i]]
      vert1 <- vert[unique==ID[i]]
      routes <- NA
      routes <- list()
      vert2 <- NA
      vert2 <- list()
      for(j in 1:n.entries) {
        for(k in 1:n.entries) {
          routes[[((j-1)*n.entries+k)]] <- detectroute(start=seg1[j],end=seg1[k],rivers=rivers,algorithm=algorithm)
          vert2[[((j-1)*n.entries+k)]] <- c(vert1[j],vert1[k])
        }
      }
      # seg.rep.max <- rep(0,length(rivers$lines))           ############bet i can remove this line
      seg.rep.max2 <- rep(0,length(rivers$lines))   ##############
      
      
      # calculate amounts of each segment represented in each route
      for(j in 1:length(rivers$lines)) {   # segment j
        linelength <- dim(rivers$lines[[j]])[1]
        
        # j <- 2
        subsegused <- rep(F,length(subseglength[[j]]))   ####################
        
        for(k in 1:length(routes)) {       # route k
          # if(length(subsegused)>length(subseglength[[j]])) print(k)
          # k<-8
          
          # if segment j shows up in route k
          if(length(routes[[k]][routes[[k]]==j])>0) {
            
            #middle
            if(length(routes[[k]])>2 & routes[[k]][1]!=j & routes[[k]][length(routes[[k]])]!=j) {
              seg.rep.max[j] <- rivers$lengths[j]            ############bet i can remove this line
              # if(map) lines(rivers$lines[[j]][,1],rivers$lines[[j]][,2],col=4,lwd=3)
              if(map) lines(rivers$lines[[j]],col=4,lwd=3)
              subsegused[] <- T  ##################
            }
            
            #only
            if(length(routes[[k]])==1) {
              # if(riverdistance(startseg=j,endseg=j,startvert=vert2[[k]][1],endvert=vert2[[k]][2],rivers,algorithm=algorithm) > seg.rep.max[j]) {
                # seg.rep.max[j] <- riverdistance(startseg=j,endseg=j,startvert=vert2[[k]][1],endvert=vert2[[k]][2],rivers=rivers,map=map,add=T,algorithm=algorithm)            ############bet i can remove this line
                # subsegused[(vert2[[k]][1]):(vert2[[k]][2]-1)] <- T   ####################   i feel like this was supposed to be commented out
                if(vert2[[k]][1] != vert2[[k]][2]) {
                  subsegused[min((linelength-1),min(vert2[[k]])):(max(vert2[[k]])-1)] <- T   #################### hack
                  if(map) lines(rivers$lines[[j]][(vert2[[k]][1]:vert2[[k]][2]),,drop=F],col=4,lwd=3)
                }
              # }
            }
            
            #beginning
            if(routes[[k]][1]==j & length(routes[[k]])>1) {
              # connected at beginning
              if(rivers$connections[routes[[k]][1],routes[[k]][[2]]]<=2) {
                # if(riverdistance(startseg=j,endseg=j,startvert=vert2[[k]][1],endvert=1,rivers=rivers,algorithm=algorithm) > seg.rep.max[j]) {
                  # seg.rep.max[j] <- riverdistance(startseg=j,endseg=j,startvert=vert2[[k]][1],endvert=1,rivers=rivers,map=map,add=T,algorithm=algorithm)            ############bet i can remove this line
                  subsegused[1:(vert2[[k]][1]-1)] <- T   #####################
                  if(map) lines(rivers$lines[[j]][(1:vert2[[k]][1]),,drop=F],col=4,lwd=3)
                # }
              }
              # connected at end
              if(any(rivers$connections[routes[[k]][1],routes[[k]][[2]]]==3:4)) {
                # if(riverdistance(startseg=j,endseg=j,startvert=vert2[[k]][1],endvert=linelength,rivers=rivers,algorithm=algorithm) > seg.rep.max[j]) {
                  # seg.rep.max[j] <- riverdistance(startseg=j,endseg=j,startvert=vert2[[k]][1],endvert=linelength,rivers=rivers,map=map,add=T,algorithm=algorithm)            ############bet i can remove this line
                  if(vert2[[k]][1] < linelength) subsegused[(vert2[[k]][1]):(linelength-1)] <- T ###################     and same fix here
                  if(map) lines(rivers$lines[[j]][(vert2[[k]][1]):linelength,,drop=F],col=4,lwd=3)
                # }
              }
              #special braided case ---------------------------------------------------------------------------------------
              if(length(routes[[k]]==2) & rivers$connections[j,routes[[k]][2]]==5) {
                d1 <- rivers$cumuldist[[j]][vert2[[k]][1]] + rivers$cumuldist[[routes[[k]][2]]][vert2[[k]][2]]
                d2 <- (rivers$lengths[j] - rivers$cumuldist[[j]][vert2[[k]][1]]) + (rivers$lengths[routes[[k]][2]] - rivers$cumuldist[[routes[[k]][2]]][vert2[[k]][2]])
                if(d1 <= d2) {
                  # if(rivers$cumuldist[[j]][vert2[[k]][1]] > seg.rep.max[j]) {
                    # seg.rep.max[j] <- riverdistance(startseg=j,endseg=j,startvert=vert2[[k]][1],endvert=1,rivers=rivers,map=map,add=T,algorithm=algorithm)            ############bet i can remove this line
                    subsegused[1:(vert2[[k]][1]-1)] <- T   #####################
                    if(map) lines(rivers$lines[[j]][(vert2[[k]][1]:1),,drop=F],col=4,lwd=3)
                  # } 
                } else {
                  # if(rivers$lengths[j] - rivers$cumuldist[[j]][vert2[[k]][1]] > seg.rep.max[j]) {
                    # seg.rep.max[j] <- riverdistance(startseg=j,endseg=j,startvert=vert2[[k]][1],endvert=linelength,rivers=rivers,map=map,add=T,algorithm=algorithm)            ############bet i can remove this line
                    if(vert2[[k]][1] < linelength) subsegused[(vert2[[k]][1]):(linelength-1)] <- T ###################     and same fix here
                    if(map) lines(rivers$lines[[j]][(vert2[[k]][1]:linelength),,drop=F],col=4,lwd=3)
                  # }
                }
              }
              if(length(routes[[k]]==2) & rivers$connections[j,routes[[k]][2]]==6) {
                d1 <- rivers$cumuldist[[j]][vert2[[k]][1]] + (rivers$lengths[routes[[k]][2]] - rivers$cumuldist[[routes[[k]][2]]][vert2[[k]][2]])
                d2 <- (rivers$lengths[j] - rivers$cumuldist[[j]][vert2[[k]][1]]) + rivers$cumuldist[[routes[[k]][2]]][vert2[[k]][2]]
                if(d1 <= d2) {
                  # if(rivers$cumuldist[[j]][vert2[[k]][1]] > seg.rep.max[j]) {
                    # seg.rep.max[j] <- riverdistance(startseg=j,endseg=j,startvert=vert2[[k]][1],endvert=1,rivers=rivers,map=map,add=T,algorithm=algorithm)            ############bet i can remove this line
                    subsegused[1:(vert2[[k]][1]-1)] <- T   #####################
                    if(map) lines(rivers$lines[[j]][(vert2[[k]][1]:1),,drop=F],col=4,lwd=3)
                  # } 
                } else {
                  # if(rivers$lengths[j] - rivers$cumuldist[[j]][vert2[[k]][1]] > seg.rep.max[j]) {
                    # seg.rep.max[j] <- riverdistance(startseg=j,endseg=j,startvert=vert2[[k]][1],endvert=linelength,rivers=rivers,map=map,add=T,algorithm=algorithm)            ############bet i can remove this line
                    if(vert2[[k]][1] < linelength) subsegused[(vert2[[k]][1]):(linelength-1)] <- T ###################     and same fix here
                    if(map) lines(rivers$lines[[j]][(vert2[[k]][1]:linelength),,drop=F],col=4,lwd=3)
                  # }
                }
              }   # -------------------------------------------------------------------------------------------------------
            }
            
            #end
            if(routes[[k]][length(routes[[k]])]==j & length(routes[[k]])>1) {
              # connected at beginning
              if(rivers$connections[routes[[k]][length(routes[[k]])],routes[[k]][[length(routes[[k]])-1]]]<=2) {
                # if(riverdistance(startseg=j,endseg=j,startvert=vert2[[k]][2],endvert=1,rivers=rivers,algorithm=algorithm) > seg.rep.max[j]) {
                  # seg.rep.max[j] <- riverdistance(startseg=j,endseg=j,startvert=vert2[[k]][2],endvert=1,rivers=rivers,map=map,add=T,algorithm=algorithm)            ############bet i can remove this line
                  subsegused[1:(vert2[[k]][2]-1)] <- T   #####################
                  if(map) lines(rivers$lines[[j]][(1:vert2[[k]][2]),,drop=F],col=4,lwd=3)
                # }
              }
              # connected at end
              if(any(rivers$connections[routes[[k]][length(routes[[k]])],routes[[k]][[length(routes[[k]])-1]]]==3:4)) {
                # if(riverdistance(startseg=j,endseg=j,startvert=vert2[[k]][2],endvert=linelength,rivers=rivers,algorithm=algorithm) > seg.rep.max[j]) {
                  # seg.rep.max[j] <- riverdistance(startseg=j,endseg=j,startvert=vert2[[k]][2],endvert=linelength,rivers=rivers,map=map,add=T,algorithm=algorithm)            ############bet i can remove this line
                  if(vert2[[k]][2] < linelength) subsegused[(vert2[[k]][2]):(linelength-1)] <- T ###################          <-------------- HERE IT IS and fixed
                  if(map) lines(rivers$lines[[j]][(linelength:vert2[[k]][2]),,drop=F],col=4,lwd=3)
                # }
              }
              
              #special braided case ---------------------------------------------------------------------------------------
              if(length(routes[[k]]==2) & rivers$connections[routes[[k]][1],routes[[k]][2]]==5) {
                d1 <- rivers$cumuldist[[routes[[k]][1]]][vert2[[k]][1]] + rivers$cumuldist[[j]][vert2[[k]][2]]
                d2 <- (rivers$lengths[routes[[k]][1]] - rivers$cumuldist[[routes[[k]][1]]][vert2[[k]][1]]) + (rivers$lengths[j] - rivers$cumuldist[[j]][vert2[[k]][2]])
                if(d1 <= d2) {
                  # if(rivers$cumuldist[[j]][vert2[[k]][2]] > seg.rep.max[j]) {
                    # seg.rep.max[j] <- riverdistance(startseg=j,endseg=j,startvert=vert2[[k]][2],endvert=1,rivers=rivers,map=map,add=T,algorithm=algorithm)            ############bet i can remove this line
                    subsegused[1:(vert2[[k]][2]-1)] <- T   #####################
                    if(map) lines(rivers$lines[[j]][(1:vert2[[k]][2]),,drop=F],col=4,lwd=3)
                  # } 
                } else {
                  # if(rivers$lengths[j] - rivers$cumuldist[[j]][vert2[[k]][2]] > seg.rep.max[j]) {
                    # seg.rep.max[j] <- riverdistance(startseg=j,endseg=j,startvert=vert2[[k]][2],endvert=linelength,rivers=rivers,map=map,add=T,algorithm=algorithm)            ############bet i can remove this line
                    if(vert2[[k]][2] < linelength) subsegused[(vert2[[k]][2]):(linelength-1)] <- T ###################    and same fix here
                    if(map) lines(rivers$lines[[j]][(linelength:vert2[[k]][2]),,drop=F],col=4,lwd=3)
                  # }
                }
              }
              if(length(routes[[k]]==2) & rivers$connections[routes[[k]][1],routes[[k]][2]]==6) {
                d1 <- (rivers$lengths[routes[[k]][1]] - rivers$cumuldist[[routes[[k]][1]]][vert2[[k]][1]]) + rivers$cumuldist[[j]][vert2[[k]][2]]
                d2 <- rivers$cumuldist[[routes[[k]][1]]][vert2[[k]][1]] + (rivers$lengths[j] - rivers$cumuldist[[j]][vert2[[k]][2]])
                if(d1 <= d2) {
                  # if(rivers$cumuldist[[j]][vert2[[k]][2]] > seg.rep.max[j]) {
                    # seg.rep.max[j] <- riverdistance(startseg=j,endseg=j,startvert=vert2[[k]][2],endvert=1,rivers=rivers,map=map,add=T,algorithm=algorithm)            ############bet i can remove this line
                    subsegused[1:(vert2[[k]][2]-1)] <- T   #####################
                    if(map) lines(rivers$lines[[j]][(1:vert2[[k]][2]),,drop=F],col=4,lwd=3)
                  # } 
                } else {
                  # if(rivers$lengths[j] - rivers$cumuldist[[j]][vert2[[k]][2]] > seg.rep.max[j]) {
                    # seg.rep.max[j] <- riverdistance(startseg=j,endseg=j,startvert=vert2[[k]][2],endvert=linelength,rivers=rivers,map=map,add=T,algorithm=algorithm)            ############bet i can remove this line
                    if(vert2[[k]][2] < linelength) subsegused[(vert2[[k]][2]):(linelength-1)] <- T ###################    and same fix here
                    if(map) lines(rivers$lines[[j]][(linelength:vert2[[k]][2]),,drop=F],col=4,lwd=3)
                  # }
                }
              }   # -------------------------------------------------------------------------------------------------------   
            }
          }
          # print(c(k,length(subsegused)))
        }
        if(length(subsegused)!=length(subseglength[[j]])) {        ########## this bit should get removed when i fix the bug!
          print(c(i,j,k))
          print(subsegused)
          print(subseglength[[j]])
        }       ########## this bit should get removed when i fix the bug!
        seg.rep.max2[j] <- sum(subsegused*subseglength[[j]])   ################    <--------- here's where the bug gets flagged
      }
      # range[i] <- sum(seg.rep.max)
      range[i] <- sum(seg.rep.max2)
    }
    if(map) riverpoints(seg=seg[unique==ID[i]],vert=vert[unique==ID[i]],rivers=rivers,pch=15,col=4)
  }
  
  thing <- data.frame(ID,range)
  thing2 <- subset(thing,range>0)
  return(thing2)
}

#' River Distance Matrix between Two Datasets
#' @description Returns a matrix of distances between each river location in two datasets, with one expressed as rows and the other expressed as columns.
#' @param seg1 First vector of river locations (segment component).  These are expressed as rows in the output matrix.
#' @param vert1 First vector of river locations (vertex component).  These are expressed as rows in the output matrix.
#' @param seg2 Second vector of river locations (segment component).  These are expressed as columns in the output matrix.
#' @param vert2 Second vector of river locations (vertex component).  These are expressed as columns in the output matrix.
#' @param rivers The river network object to use.
#' @param logical1 A boolean vector that can be used for subsetting.  If used,
#'   \code{riverdistancetofrom} will only return distances in which a
#'   specified condition is met for the first dataset.
#' @param logical2 A boolean vector that can be used for subsetting.  If used,
#'   \code{riverdistancetofrom} will only return distances in which a
#'   specified condition is met for the second dataset.
#' @param ID1 a vector of observation IDs for the first dataset that will be used as row names in the output matrix.
#' @param ID2 a vector of observation IDs for the second dataset that will be used as column names in the output matrix.
#' @param stopiferror Whether or not to exit with an error if a route cannot be
#'   found.  If this is set to \code{FALSE} and a route cannot be found,
#'   the function will return \code{NA} in the appropriate entry.  Defaults to \code{TRUE}.  See \link{detectroute}.
#' @param algorithm Which route detection algorithm to use (\code{"Dijkstra"},
#'   \code{"sequential"}, or \code{"segroutes"}).  If left as \code{NULL} (the
#'   default), the function will automatically make a selection.  See
#'   \link{detectroute} for more details.
#' @return A matrix of distances (numeric) with rows and columns labeled by corresponding values of \code{ID}.
#' @seealso \link{riverdistance}
#' @note Building routes from the river mouth to each river network segment and/or distance lookup tables will
#'   greatly reduce computation time (see \link{buildsegroutes}).
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
#' riverdistancetofrom(seg1=streamlocs.seg, vert1=streamlocs.vert,
#'   seg2=fish.seg, vert2=fish.vert, rivers=Gulk, ID1=streamlocs.ID, ID2=fish.ID)
#'
#' logi1 <- streamlocs.ID=="B" | streamlocs.ID=="C"
#' logi2 <- fish.ID!="fish3"
#'
#' riverdistancetofrom(seg1=streamlocs.seg, vert1=streamlocs.vert,
#'   seg2=fish.seg, vert2=fish.vert, rivers=Gulk, logical1=logi1, logical2=logi2,
#'   ID1=streamlocs.ID, ID2=fish.ID)
#' @export
riverdistancetofrom <- function(seg1,vert1,seg2,vert2,rivers,logical1=NULL,logical2=NULL,ID1=NULL,ID2=NULL,stopiferror=TRUE,algorithm=NULL) {
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
      dists[i,j] <- riverdistance(startseg=seg1[i],startvert=vert1[i],endseg=seg2[j],endvert=vert2[j],rivers=rivers,
                                  stopiferror=stopiferror,algorithm=algorithm)
    }
  }
  
  class(ID1) <- "list"
  class(ID2) <- "list"
  dimnames(dists)[[1]] <- ID1
  dimnames(dists)[[2]] <- ID2
  return(dists)
}