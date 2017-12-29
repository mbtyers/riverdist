#' Dataset: Gulkana River
#' 
#' A stretch of Gulkana River and tributaries.
#' 
#' @docType data
#' @keywords datasets
#' @name Gulk
#' @usage data(Gulk)
#' @format A river network object, see \link{rivernetwork}
NULL

#' Dataset: Kenai River 1
#' 
#' A first pass at a messy river network object.
#' 
#' @docType data
#' @keywords datasets
#' @name Kenai1
#' @usage data(Kenai1)
#' @format A river network object, see \link{rivernetwork}
#' @seealso \link{Kenai2}, \link{Kenai3}
NULL

#' Dataset: Kenai River 2
#' 
#' A second pass at a messy river network object.  In this iteration of cleanup,
#' several non-connected segments have been removed.
#' 
#' @docType data
#' @keywords datasets
#' @name Kenai2
#' @usage data(Kenai2)
#' @format A river network object, see \link{rivernetwork}
#' @seealso \link{Kenai1}, \link{Kenai3}
NULL

#' Dataset: Kenai River 3
#' 
#' A third pass at a messy river network object.  In this iteration of cleanup,
#' several non-connected segments have been removed, and several series of
#' segments have been dissolved into single segments.
#' 
#' @docType data
#' @keywords datasets
#' @name Kenai3
#' @usage data(Kenai3)
#' @format A river network object, see \link{rivernetwork}
#' @seealso \link{Kenai1}, \link{Kenai2}
NULL

#' Dataset: Fakefish
#' 
#' A set of observations of Fakefish on the Gulkana River and its tributaries.
#' 
#' \itemize{ 
#' \item \code{x}. X-coordinate of observation (Alaska Albers Equal Area). Note that the locations do not align with the river network object. 
#' \item \code{y}. Y-coordinate of observation 
#' \item \code{seg}. River segment (with x- and y-coordinates snapped to river network object) 
#' \item \code{vert}. River vertex
#' \item \code{fish.id}. Numeric identifier for each fish (individual fish were observed more than once) 
#' \item \code{flight}. Numeric identifier for each telemetry flight 
#' \item \code{flight.date}. Date of each telemetry flight 
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name fakefish
#' @usage data(fakefish)
#' @format A data frame
#' @seealso \link{Gulk}
NULL

#' Dataset: Smallset
#' 
#' A small set of observations of fakefish on the Gulkana River and its tributaries.
#' 
#' \itemize{ 
#' \item \code{x}. X-coordinate of observation (Alaska Albers Equal Area). Note that the locations do not align with the river network object. 
#' \item \code{y}. Y-coordinate of observation 
#' \item \code{seg}. River segment 
#' \item \code{vert}. River vertex
#' \item \code{fish.id}. Numeric identifier for each fish (individual fish were observed more than once) 
#' \item \code{flight}. Numeric identifier for each telemetry flight
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name smallset
#' @usage data(smallset)
#' @format A data frame
#' @seealso \link{Gulk}
NULL

#' Dataset: Killey River, West Channel
#' 
#' A messy and braided section of the Kenai River network - actually a subset of \link{Kenai3}.
#' 
#' @docType data
#' @keywords datasets
#' @name KilleyW
#' @usage data(KilleyW)
#' @format A river network object, see \link{rivernetwork}
NULL

#' Dataset: Koyukuk River 1
#' 
#' A first pass at a messy river network object.  The way it was dissolved in
#' ArcGIS makes the endpoints appear disconnected to \link{line2network} and the
#' topologies do not work.
#' 
#' @docType data
#' @keywords datasets
#' @name Koyukuk1
#' @usage data(Koyukuk1)
#' @format A river network object, see \link{rivernetwork}
#' @seealso \link{Koyukuk2}
NULL

#' Dataset: Koyukuk River 2
#' 
#' A second pass at a messy river network object, with topologies fixed from \link{Koyukuk1}.
#' 
#' @docType data
#' @keywords datasets
#' @name Koyukuk2
#' @usage data(Koyukuk2)
#' @format A river network object, see \link{rivernetwork}
#' @seealso \link{Koyukuk1}
NULL

#' Dataset: Koyukuk River 0
#' 
#' An unusably messy river network object, included for the purpose of testing river network editing functions.
#' 
#' @docType data
#' @keywords datasets
#' @name Koyukuk0
#' @usage data(Koyukuk0)
#' @format A river network object, see \link{rivernetwork}
#' @seealso \link{Koyukuk1}, \link{Koyukuk2}
NULL

#' Dataset: Line 98 of Kenai River 1 (Long-Lat)
#' 
#' A matrix of coordinates in longitude-latitude, used to illustrate coordinate
#' transformation.  Coordinates come from arbitrary line number 98 in the Kenai
#' River 1 shapefile, rendered in long-lat.
#' @docType data
#' @keywords datasets
#' @name line98
#' @usage data(line98)
#' @format A matrix of values
NULL

#' Dataset: A-B Streams
#' 
#' A complex river network object, a subset of the streams in the Absaroka-Beartooth Wilderness.
#' 
#' @docType data
#' @keywords datasets
#' @name abstreams
#' @usage data(abstreams)
#' @format A river network object, see \link{rivernetwork}
NULL

#' Dataset: A-B Streams 0
#' 
#' An unusably messy river network object, included for the purpose of testing river network editing functions.
#' 
#' @docType data
#' @keywords datasets
#' @name abstreams0
#' @usage data(abstreams0)
#' @format A river network object, see \link{rivernetwork}
NULL


#' Dataset: Fakefish Density
#' 
#' An object created by \link{riverdensity}, describing the density of Fakefish points in the Gulkana River during ten surveys.
#' 
#' Intended for plotting using \link{plotriverdensity}.
#' 
#' @docType data
#' @keywords datasets
#' @name fakefish_density
#' @usage data(fakefish_density)
#' @format A river density object, see \link{riverdensity}, \link{plotriverdensity}, \link{riverdensity-class}
NULL

#' The "rivernetwork" Class
#'
#' A class that holds spatial coordinates for river networks, as well as network topology and attributes.
#'
#' Created by \link{line2network} from an input line shapefile.  Contains all information for network distance calculation, plotting, etc. in the 'riverdist' package.
#'
#' Plotting methods are described in \link{plot.rivernetwork}.
#'@section Elements: 
#'  \describe{
#'    \item{\code{sp}:}{Object of class \code{"SpatialLinesDataFrame"} from package 'sp'; see \link[sp]{SpatialLinesDataFrame-class}.  This is the original object as read by \link[rgdal]{readOGR}, and is preserved to maintain plotting capability.}
#'    \item{\code{lines}:}{Object of class \code{"list"}.  Each list element is a matrix of XY coordinates of the vertices of a single river segment.}
#'    \item{\code{connections}:}{Object of class \code{"matrix"}, with \code{"numeric"} elements.  Defined as a square matrix, with elements describing the type of connection detected between line segments.
#'      \itemize{
#'      \item A value of 1 in element \code{[i,j]} indicates that the beginning of segment \code{i} is connected to the beginning of segment \code{j}.
#'      \item A value of 2 in element \code{[i,j]} indicates that the beginning of segment \code{i} is connected to the end of segment \code{j}.
#'      \item A value of 3 in element \code{[i,j]} indicates that the end of segment \code{i} is connected to the beginning of segment \code{j}.
#'      \item A value of 4 in element \code{[i,j]} indicates that the end of segment \code{i} is connected to the end of segment \code{j}.
#'      \item A value of 5 in element \code{[i,j]} indicates that segments \code{i} and \code{j} are connected at both beginning and end.
#'      \item A value of 6 in element \code{[i,j]} indicates that the beginning of segment \code{i} is connected to the end of segment \code{j}, and the end of segment \code{i} is connected to the beginning of segment \code{j}.
#'      \item A value of NA in element \code{[i,j]} indicates that segments \code{i} and \code{j} are not connected.}}
#'    \item{\code{lengths}:}{Vector of class \code{"numeric"}.  Defined as the calculated total lengths of each river segment.}
#'    \item{\code{names}:}{Vector of class \code{"character"}.  Defined as the names of each river segment.}
#'    \item{\code{mouth}:}{Object of class \code{"list"}, with two elements.  Element \code{mouth.seg} gives the segment number of the mouth (lowest point) of the river network, and \code{mouth.vert} gives the vertex number.}
#'    \item{\code{sequenced}:}{\code{"logical"}: has value of TRUE if line vertices have been stored in upstream sequence using \link{sequenceverts}.}
#'    \item{\code{tolerance}:}{\code{"numeric"}: the spatial tolerance that was used in determining river segment endpoint connectivity; see \link{line2network}, \link{splitsegments}.}
#'    \item{\code{units}:}{\code{"character"}: the spatial units detected from the input shapefile.}
#'    \item{\code{lineID}:}{Object of class \code{"data.frame"} establishing the relationship between river segments as stored in the \code{sp} and \code{lines} elements, and is used for updating the \code{sp} element during river network editing in \link{dissolve}, \link{splitsegments}, \link{sequenceverts}, \link{trimriver}, and \link{trimtopoints}.
#'    \itemize{
#'      \item \code{rivID} gives the list element number of each river segment in \code{lines}.  This is the same number that is used for segment numbering in river locations.
#'      \item \code{sp_line} gives the corresponding list element in \code{sp@@lines}.
#'      \item \code{sp_seg} gives the corresponding list element in \code{sp@@lines[[]]@@Lines}.
#'      }}
#'    \item{\code{braided}:}{\code{"logical"}: Has value of \code{TRUE} if \link{checkbraidedTF} has detected braiding, \code{FALSE} if no braiding has been detected, and \code{NA} if braiding has not yet been checked.}
#'    \item{\code{cumuldist}:}{List of class \code{"numeric"}: Each element is a vector of cumulative distances along each river segment, beginning with 0.}
#'    \item{\code{segroutes}:}{Object of class \code{"list"}, with each element defined as a vector of class \code{"numeric"}, describing the route from the mouth segment to the specific segment.  This element only exists if \link{buildsegroutes} has been run, and can greatly speed up route and distance calculation.}
#'    \item{\code{distlookup}:}{List of three matrices, of class \code{"numeric"} or \code{"logical"}.  Element \code{[i,j]} of each matrix corresponds to
#'   the route between segment \code{i} and \code{j}.  The
#'   \code{distlookup$middist} matrix gives the total distance of the "middle"
#'   of each route (between the starting and ending segments"), and the
#'   \code{distlookup$starttop} and \code{distlookup$endtop} matrices have value
#'   \code{TRUE}, \code{FALSE}, or \code{NA} if the segments at the beginning or
#'   end of the route are connected to the rest of the route at the top of the
#'   coordinate matrix, bottom of the coordinate matrix, or if the route is
#'   contained to just one segment, respectively. }
#'   }
#' @name rivernetwork 
#' @rdname rivernetwork
#' @aliases rivernetwork-class
#' @exportClass rivernetwork
#' @author Matt Tyers
NULL

#' The "riverdensity" Class
#'
#' A class that holds density information computed from point data along a river network.
#'
#' Created by \link{makeriverdensity} from point data and a river network.  Contains all information for plotting in \link{plot.riverdensity}.
#'
#'@section Elements:
#'  \describe{
#'    \item{\code{densities}:}{Object of class \code{"list"}. Each list element corresponds to a unique value of survey.  Each element is itself of class \code{"list"}, with each element corresponding to a segment from the associated river network.  Each element is a vector of class \code{"numeric"}, with values equal to the scaled densities calculated at the river network vertices stored in \code{$densverts} of the associated river network segment.}
#'    \item{\code{endptverts}:}{List of vectors of class \code{"numeric"}.  Each list element is a vector of the vertices of the endpoints of the subsegments considered for density calculation.  Each list element corresponds to a river segment from the associated river network.}
#'    \item{\code{densverts}:}{List of vectors of class \code{"numeric"}.  Each element is a vector of the vertices of the points of the subsegments considered for density calculation, that were used for density calculation.  Each list element corresponds to a river segment from the associated river network.}
#'    \item{\code{pointsegs}:}{Vector of class \code{"numeric"}.  Defined as the segment numbers of the point data used for density calculation.}
#'    \item{\code{pointverts}:}{Vector of class \code{"numeric"}.  Defined as the vertex numbers of the point data used for density calculation.}
#'    \item{\code{survey}:}{Vector of class \code{"numeric"} or class \code{"character"}.  Defined as the survey identifiers associated with the point data used for density calculation.}
#'    \item{\code{rivers}:}{Object of class \code{"rivernetwork"} ; see \link{rivernetwork-class}}.
#'   }
#' @name riverdensity
#' @rdname riverdensity
#' @aliases riverdensity-class
#' @exportClass riverdensity
#' @author Matt Tyers
NULL

#' The "homerange" Class
#'
#' A class that holds information computed from the \link{homerange} function.  Contains all information for plotting in \link{plot.homerange}.
#'
#'@section Elements:
#'  \describe{
#'    \item{\code{ranges}:}{Object of class \code{"data.frame"}. Contains a column of the identifiers for each individual, and a column of the associated home ranges.}
#'    \item{\code{subseg_n}:}{List of the number of times each subsegment was traveled.  The first level of the list corresponds to individual, the second level to river segment.}
#'    \item{\code{subseg_length}:}{List of lengths of each subsegment.}
#'    \item{\code{seg, vert, unique, rivers}:}{All inputs from the original \link{homerange} call.}
#'   }
#' @name homerange-class
#' @rdname homerange-class
#' @exportClass homerange
#' @author Matt Tyers
NULL
