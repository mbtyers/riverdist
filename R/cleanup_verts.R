#' Interactive Cleanup of the Vertices of Individual Segments
#' @description A trial version of a function for deep-cleaning a river network.
#' 
#' Sometimes a shapefile contains errors that are not obvious at an initial check, typically vertices that should not be there.  
#' 
#' This function steps through each segment in sequence, and allows the user to interactively remove vertices.
#' @param rivers The river network object to use
#' @param startwith The segment (number) to start with, defaulting to \code{1}.  
#' @note Stepping through a large and messy river network can be time-consuming.  To resume a cleanup session, use the \code{startwith=} argument and the last returned river network.  For example, if \code{rivers1 <- cleanup_verts(rivers)} were initially called and the user selected "save & close" at segment 100, cleanup can be resumed by calling \code{rivers2 <- cleanup_verts(rivers1, startwith=100)}.
#' @return A new river network object, see
#'   \link{rivernetwork}
#' @seealso line2network
#' @author Matt Tyers
#' @examples
#' data(abstreams0,Koyukuk0,Kenai1)
#' 
#' # abstreams_fixed1 <- cleanup_verts(abstreams0)
#' # Koyukuk <- cleanup(Koyukuk0)
#' # Kenai <- cleanup(Kenai1)
#' @importFrom graphics plot
#' @importFrom graphics text
#' @export
cleanup_verts <- function(rivers, startwith=1) {
  lines <- rivers$lines
  connections <- rivers$connections
  
  for(i in 1:length(rivers$lines)) {
    ll <- nrow(rivers$lines[[i]])
    keep <- c(T,(rivers$lines[[i]][2:ll,1]!=rivers$lines[[i]][1:(ll-1),1] & rivers$lines[[i]][2:ll,2]!=rivers$lines[[i]][1:(ll-1),2]))
    rivers$lines[[i]] <- rivers$lines[[i]][keep,]
  }
  
  done <- F
  segi <- startwith
  
  
  while(!done) {
    theguys <- whoconnected(segi,rivers)
    zoomtoseg(c(segi,theguys),rivers=rivers,segmentnum=F)
    highlightseg(theguys,add=T,rivers=rivers,lwd=3,color=T,cex=.6)
    # highlightseg(segi,add=T,rivers=rivers,lwd=4,color=T)
    lines(rivers$lines[[segi]],lwd=1,col=1)
    topologydots(rivers,add=T)
    points(rivers$lines[[segi]])#,pch="+")
    accept <- 0
    zoom <- F
    while(!accept %in% c("Y","y","S","s")) {
      
      # for(i in 2:(nrow(rivers$lines[[segi]])-1)) {
      #   thesesegs <- blm3[blm1==rivers$lines[[segi]][i,1] & blm2==rivers$lines[[segi]][i,2] & blm3!=segi]
      #   if(length(thesesegs)>0) text(rivers$lines[[segi]][i,1], rivers$lines[[segi]][i,2], labels=paste("coincident with segment(s)",paste(thesesegs,collapse=",")), pos=4)
      # }
      
      while(!(accept %in% c("Y","y","N","n","Z","z","s","S"))) accept <- readline(prompt=paste("Segment",segi,"of",length(rivers$lines),"-- Accept segment? (y)es / (n)o, (z)oom, or (s)ave & close: "))
      if(accept %in% c("Z","z"))  {
        zoomtoseg(segi,rivers=rivers,segmentnum=F)
        highlightseg(theguys,add=T,rivers=rivers,lwd=3,color=T,cex=.6)
        # highlightseg(segi,add=T,rivers=rivers,lwd=4,color=T)
        lines(rivers$lines[[segi]],lwd=1,col=1)
        # highlightseg(segi,add=T,rivers=rivers,lwd=3,color=T)
        
        coords <- rivers$lines[[segi]]
        coordID <- 1:nrow(coords)
        coordspaste <- paste(coords[,1],coords[,2])
        
        dupes <- list()
        dupesi <- 1
        for(i in 1:length(coordspaste)) {
          if(coordspaste[i] %in% coordspaste[-i]) {
            dupes[[dupesi]] <- which(coordspaste==coordspaste[i])
            dupesi <- dupesi+1
          }
        }
        
        if(length(dupes)>0) {
          dupes <- unique(dupes)
          for(dupesi in 1:length(dupes)) {
            dupesilab <- paste(dupes[[dupesi]],collapse=",")
            text(coords[dupes[[dupesi]][1],1],coords[dupes[[dupesi]][1],2],labels=dupesilab,pos=4,col="grey30",cex=.8)
          }
        }
        
        text(coords[!(coordID %in% unlist(dupes)),],labels=coordID[!(coordID %in% unlist(dupes))],col="grey30",pos=4,cex=.8)
        
        topologydots(rivers,add=T)
        points(rivers$lines[[segi]])#,pch="+")
        
        accept <- 0
        zoom <- T
      }
      
      if(accept %in% c("N","n")) {
        
        coords <- rivers$lines[[segi]]
        coordID <- 1:nrow(coords)
        coordspaste <- paste(coords[,1],coords[,2])
        
        dupes <- list()
        dupesi <- 1
        for(i in 1:length(coordspaste)) {
          if(coordspaste[i] %in% coordspaste[-i]) {
            dupes[[dupesi]] <- which(coordspaste==coordspaste[i])
            dupesi <- dupesi+1
          }
        }
        if(length(dupes)>0) {
          dupes <- unique(dupes)
          for(dupesi in 1:length(dupes)) {
            dupesilab <- paste(dupes[[dupesi]],collapse=",")
            text(coords[dupes[[dupesi]][1],1],coords[dupes[[dupesi]][1],2],labels=dupesilab,pos=4,col="grey30",cex=.8)
          }
        }
        
        text(coords[!(coordID %in% unlist(dupes)),],labels=coordID[!(coordID %in% unlist(dupes))],col="grey30",pos=4,cex=.8)
        
        whichones <- readline("Enter verts to remove, separated by commas (or 0 to undo): ")
        whichones <- as.numeric(unlist(strsplit(whichones, ",")))
        
        if(whichones[1]>0) {
          linessegi <- rivers$lines[[segi]][-whichones,]
          
          rivers1 <- rivers
          rivers1$lines[[segi]] <- linessegi
          rivers1$connections <- calculateconnections(rivers1$lines,rivers1$tolerance)
          
          ## plot stuff  ----------- get this right this time
          
          if(zoom) zoomtoseg(segi,rivers=rivers1,segmentnum=F)
          else zoomtoseg(c(segi,theguys),rivers=rivers1,segmentnum=F)
          
          highlightseg(theguys,add=T,rivers=rivers1,lwd=3,color=T)
          lines(rivers1$lines[[segi]],lwd=1,col=1)
          
          text(linessegi,labels=1:nrow(linessegi),col="grey30",pos=4,cex=.8)
          
          topologydots(rivers1,add=T)
          points(linessegi)#,pch="+")
          
          accept <- 0
          while(!(accept %in% c("Y","y","N","n"))) accept <- readline(prompt="Accept changes? (y/n) ")
        }
        if(accept %in% c("Y","y")) rivers$lines[[segi]] <- linessegi
        if(accept %in% c("N","n")) {
          theguys <- whoconnected(segi,rivers)
          if(zoom) zoomtoseg(segi,rivers=rivers,segmentnum=F)
          else zoomtoseg(c(segi,theguys),rivers=rivers,segmentnum=F)
          highlightseg(theguys,add=T,rivers=rivers,lwd=3,color=T,cex=.6)
          # highlightseg(segi,add=T,rivers=rivers,lwd=4,color=T)
          lines(rivers$lines[[segi]],lwd=1,col=1)
          topologydots(rivers,add=T)
          points(rivers$lines[[segi]])#,pch="+")
        }
        if(whichones[1]<=0) accept <- T
      }
    }
    if(segi >= length(rivers$lines) | (accept %in% c("S","s"))) done <- T
    if(accept %in% c("S","s")) cat('\n',"Cleanup session ended and rivernetwork returned.  Run again with cleanup_verts(..., startwith=",segi,")",sep="")
    segi <- segi+1
  }
  
  # names(schafer)
  # "sp"          - done I think
  # "lineID"      - done I think   ---- actually shouldn't be needed
  # "lines"       - done
  # "connections" - done
  # "lengths"     - done
  # "names"       - not needed
  # "mouth"       - not needed (though maybe build it in if i can)
  # "sequenced"   - not needed
  # "tolerance"   - not needed
  # "units"       - not needed
  # "braided"     - not needed I think
  # "cumuldist"   - done
  
  # rivers$lines <- lines
  rivers$connections <- calculateconnections(lines,rivers$tolerance)
  #### finish updating rivers object!!!!
  rivers <- addcumuldist(rivers)
  
  # # updating sp object
  # riverssplines <- rivers$sp@lines
  # for(i in 1:length(rivers$lines)) {
  #   # rivers$sp@lines[[rivers$lineID[i,2]]][[rivers$lineID[i,3]]] <- lines[[i]]
  #   riverssplines[[rivers$lineID[i,2]]]@Lines[[rivers$lineID[i,3]]]@coords <- rivers$lines[[i]]
  # }   ## really hope this works
  # rivers$sp@lines <- riverssplines  #### only run this if it does!!!
  
  rivers <- update_sf(rivers)
  
  # making a vector of total segment lengths
  lengths <- rep(NA,length(rivers$lines))
  for(i in 1:length(rivers$lines)) {
    lengths[i] <- pdisttot(rivers$lines[[i]])
  }
  rivers$lengths <- lengths
  
  return(rivers)
  
}