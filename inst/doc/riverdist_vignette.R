## ----eval=FALSE----------------------------------------------------------
#  library(riverdist)
#  MyRivernetwork <- line2network(path=".", layer="MyShapefile")
#  
#  # Re-projecting in Alaska Albers Equal Area projection:
#  AKalbers <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154
#      +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
#  
#  MyRivernetwork <- line2network(path=".", layer="MyShapefile", reproject=AKalbers)

## ----fig.width=5,fig.height=5--------------------------------------------
library(riverdist)
data(Gulk)
plot(x=Gulk)

## ----eval=FALSE----------------------------------------------------------
#  abstreams_fixed <- cleanup(abstreams0)

## ----fig.width=5,fig.height=5--------------------------------------------
topologydots(rivers=Gulk)

## ----fig.width=5,fig.height=3--------------------------------------------
data(fakefish)
fakefish_riv <- xy2segvert(x=fakefish$x, y=fakefish$y, rivers=Gulk)
head(fakefish_riv)  # a look at the first few rows of the output
hist(fakefish_riv$snapdist, main="snapping distance (m)")

## ----eval=FALSE----------------------------------------------------------
#  segvert_from_shp <- pointshp2segvert(path=".", layer="MyPointShapefile", rivers=MyRivernetwork)

## ----fig.width=5,fig.height=5--------------------------------------------
zoomtoseg(seg=c(11, 14), rivers=Gulk)
points(fakefish$x, fakefish$y, pch=16, col="red")
riverpoints(seg=fakefish_riv$seg, vert=fakefish_riv$vert, rivers=Gulk, pch=15, col="blue")

## ----fig.width=5,fig.height=5--------------------------------------------
# starting location: segment 7, vertex 49
# ending location: segment 14, vertex 121
detectroute(start=7, end=14, rivers=Gulk)
riverdistance(startseg=7, startvert=49, endseg=14, endvert=121, rivers=Gulk, map=TRUE)

## ------------------------------------------------------------------------
data(smallset)
smallset
riverdistanceseq(unique=smallset$id, survey=smallset$flight, seg=smallset$seg, 
                   vert=smallset$vert, rivers=Gulk)

## ------------------------------------------------------------------------
riverdistancematbysurvey(indiv=1, unique=smallset$id, survey=smallset$flight,
      seg=smallset$seg, vert=smallset$vert, rivers=Gulk, full=FALSE)

## ----fig.width=7,fig.height=3.5------------------------------------------
# calculating observed minimum home range for all individuals
par(mfrow=c(1,3))
homerange(unique=smallset$id, seg=smallset$seg, vert=smallset$vert, rivers=Gulk, map=TRUE)

## ------------------------------------------------------------------------
dmat <- riverdistancemat(smallset$seg,smallset$vert,Gulk)
round(dmat)[1:7,1:7]  # only showing the first 7 rows & columns for clarity

## ------------------------------------------------------------------------
logi1 <- (smallset$seg==2)
obsID <- paste0("id",smallset$id,"-flight",smallset$flight)  # constructing observation labels
riverdistancemat(seg=smallset$seg, vert=smallset$vert, rivers=Gulk, logical=logi1, ID=obsID)

## ------------------------------------------------------------------------
streamlocs.seg <- c(2,2,2)
streamlocs.vert <- c(10,20,30)
streamlocs.ID <- c("loc A","loc B","loc C")

logi2 <- (smallset$seg==2)
obsID <- paste0("id",smallset$id,"-flight",smallset$flight)

riverdistancetofrom(seg1=streamlocs.seg, vert1=streamlocs.vert, seg2=smallset$seg, 
                    vert2=smallset$vert, ID1=streamlocs.ID, ID2=obsID, logical2=logi2, 
                    rivers=Gulk)

## ----fig.width=5,fig.height=5--------------------------------------------
showends(seg=1,rivers=Gulk)
Gulk1 <- setmouth(seg=1, vert=1, rivers=Gulk)

## ----fig.width=5,fig.height=4--------------------------------------------
zoomtoseg(seg=c(6,3), rivers=Gulk)
riverpoints(seg=c(6,4), vert=c(250,250), col=4, pch=15, rivers=Gulk1)
#riverdistance(startseg=6, endseg=4, startvert=250, endvert=250, rivers=Gulk1, map=TRUE)
text(c(859122.4, 872104.1), c(6964127.4,6969741.0), pos=c(3, 4), labels=c("beginning", "end"))
riverdirection(startseg=6, endseg=4, startvert=250, endvert=250, rivers=Gulk1)
upstream(startseg=6, endseg=4, startvert=250, endvert=250, rivers=Gulk1, net=FALSE)
upstream(startseg=6, endseg=4, startvert=250, endvert=250, rivers=Gulk1, net=TRUE)
upstream(startseg=6, endseg=4, startvert=250, endvert=250, rivers=Gulk1, flowconnected=TRUE)

## ------------------------------------------------------------------------
data(abstreams)

tstart <- Sys.time()
riverdistance(startseg=120, startvert=10, endseg=131, endvert=10, rivers=abstreams, 
              algorithm="sequential")
Sys.time()- tstart

tstart <- Sys.time()
riverdistance(startseg=120, startvert=10, endseg=131, endvert=10, rivers=abstreams, 
              algorithm="Dijkstra")
Sys.time()- tstart

tstart <- Sys.time()
riverdistance(startseg=120, startvert=10, endseg=131, endvert=10, rivers=abstreams)

# Note: it is not necessary to specify the algorithm here: the distance function
# will automatically select the fastest algorithm unless otherwise specified.
Sys.time()- tstart

## ----eval=FALSE----------------------------------------------------------
#  data(Gulk, fakefish)
#  fakefish_density <- makeriverdensity(seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk,
#    survey=fakefish$flight.date, resolution=2000, bw=10000)
#  par(mfrow=c(3,3))
#  plot(x=fakefish_density, ramp="blue", dark=0.85, maxlwd=15,
#                   whichplots=c(1:8,10))   # showing only nine plots for clarity

## ----fig.width=7.5,fig.height=8.5,echo=FALSE, eval=FALSE-----------------
#  data(fakefish_density)
#  par(mfrow=c(3,4))
#  plot(x=fakefish_density, ramp="blue", dark=0.85, maxlwd=15)

## ----fig.width=7.5,fig.height=10,echo=FALSE------------------------------
data(fakefish_density)
par(mfrow=c(3,3))
plot(x=fakefish_density, ramp="blue", dark=0.85, maxlwd=15,
                 whichplots=c(1:8,10))

## ----fig.width=5,fig.height=4--------------------------------------------
x <- homerange(unique=fakefish$flight.date, seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk)
x
plot(x$range, type='b', xaxt='n', xlab="", ylab="range (m)", cex.axis=.6)
axis(1,at=1:10,labels=sort(unique(fakefish$flight.date)), cex.axis=.6, las=3)

## ----fig.width=7.5,fig.height=7.5----------------------------------------
par(mfrow=c(3,3))
kfunc(seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk, survey=fakefish$flight.date,
     maxdist=200000, whichplots=c(1:8,10))   # showing only nine plots for clarity

## ------------------------------------------------------------------------
mouthdist(seg=fakefish$seg[1], vert=fakefish$vert[1], rivers=Gulk)
x <- mouthdistbysurvey(unique=fakefish$fish.id, survey=fakefish$flight,
    seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk)

round(x)

## ----fig.width=7.5,fig.height=5------------------------------------------
x <- mouthdistbysurvey(unique=fakefish$fish.id, survey=fakefish$flight.date,
    seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk)

par(mfrow=c(1,2))
plotseq(seqbysurvey=x)
plotseq(seqbysurvey=x, type="dotline")
par(mfrow=c(1,1))

## ----fig.width=7.5,fig.height=5------------------------------------------
x <- upstreamseq(unique=fakefish$fish.id, survey=fakefish$flight, seg=fakefish$seg,
                 vert=fakefish$vert, rivers=Gulk)
round(x)
par(mfrow=c(1,2))
plotseq(seqbysurvey=x)

plotseq(seqbysurvey=x, type="dotplot")
abline(h=0, lty=3)

## ----fig.width=5,fig.height=10-------------------------------------------
matbysurveylist <- matbysurveylist(unique=fakefish$fish.id, survey=fakefish$flight, seg=fakefish$seg, 
                         vert=fakefish$vert, rivers=Gulk)
plotmatbysurveylist(matbysurveylist)

## ----eval=FALSE----------------------------------------------------------
#  data(abstreams0)  # a messy river network
#  abstreams_fixed <- cleanup(abstreams0)  # fixing many problems

## ----fig.width=7.5,fig.height=4, warning=FALSE, message=FALSE------------
Gulk_trim1 <- trimriver(trim=c(10,11,12,13,14), rivers=Gulk)
Gulk_trim2 <- trimriver(trimto=c(10,11,12,13,14), rivers=Gulk)

par(mfrow=c(1,3))
plot(x=Gulk, main="original")
plot(x=Gulk_trim1, main="trim=c(10,11,12,13,14)")
plot(x=Gulk_trim2, main="trimto=c(10,11,12,13,14)")

## ----fig.width=5,fig.height=5, message=FALSE-----------------------------
data(Kenai3)
x <- c(174185, 172304, 173803, 176013)
y <- c(1173471, 1173345, 1163638, 1164801)
Kenai3.buf1 <- trimtopoints(x=x, y=y, rivers=Kenai3, method="snap")
Kenai3.buf2 <- trimtopoints(x=x, y=y, rivers=Kenai3, method="snaproute")
Kenai3.buf3 <- trimtopoints(x=x, y=y, rivers=Kenai3, method="buffer", dist=5000)

plot(x=Kenai3, main="original")
points(x, y, pch=15, col=4)
legend(par("usr")[1], par("usr")[4], legend="points to buffer around", pch=15, col=4, cex=.6)

## ----fig.width=7.5,fig.height=3.5----------------------------------------

par(mfrow=c(1,3))
plot(x=Kenai3.buf1, main="snap")
points(x, y, pch=15, col=4)
plot(x=Kenai3.buf2, main="snaproute")
points(x, y, pch=15, col=4)
plot(x=Kenai3.buf3, main="buffer, dist=5000")
points(x, y, pch=15, col=4)

## ----fig.width=7.5,fig.height=5, message=FALSE---------------------------
data(Koyukuk2)
Koy_subset <- trimriver(trimto=c(30,28,29,3,19,27,4),rivers=Koyukuk2)
Koy_subset <- setmouth(seg=1,vert=427,rivers=Koy_subset)

Koy_subset_trim <- removeunconnected(Koy_subset)

par(mfrow=c(1,2))
plot(x=Koy_subset, main="original")
plot(x=Koy_subset_trim, main="unconnected segments removed")

## ----fig.width=7.5,fig.height=5, warning=FALSE, message=FALSE------------
data(Kenai2)
Kenai2_sub <- trimriver(trimto=c(26,157,141,69,3,160,2,35,102,18,64,86,49,103,61
                                 ,43,183,72,47,176), rivers=Kenai2)

Kenai2_sub_dissolve <- dissolve(rivers=Kenai2_sub)

par(mfrow=c(1,2))
plot(x=Kenai2_sub, main="original")
plot(x=Kenai2_sub_dissolve, main="dissolved")

## ----fig.width=7.5,fig.height=4.5, warning=FALSE, message=FALSE----------
data(Koyukuk1)

Koyukuk1.split <- splitsegments(rivers=Koyukuk1)

par(mfrow=c(1,2))
topologydots(rivers=Koyukuk1, main="original")
topologydots(rivers=Koyukuk1.split, main="split")

## ----fig.width=7.5,fig.height=4, warning=FALSE, message=FALSE------------
data(Koyukuk0)

Koyukuk0.1 <- connectsegs(connect=21, connectto=20, rivers=Koyukuk0)

par(mfrow=c(1,2))
plot(Koyukuk0, ylim=c(1930500,1931500), xlim=c(194900,195100), main="original")
topologydots(Koyukuk0, add=TRUE)

plot(Koyukuk0.1,ylim=c(1930500,1931500), xlim=c(194900,195100), main="connected")
topologydots(Koyukuk0.1, add=TRUE)

## ----eval=FALSE----------------------------------------------------------
#  data(abstreams0)
#  abstreams2 <- removemicrosegs(abstreams0)

## ----fig.width=7.5,fig.height=4, warning=FALSE, message=FALSE------------
data(Kenai3)
Kenai3split <- addverts(Kenai3, mindist=200)

par(mfrow=c(1,2))
zoomtoseg(seg=c(47,74,78), rivers=Kenai3, main="segment 74 vertices")
points(Kenai3$lines[[74]]) 

zoomtoseg(seg=c(47,74,78), rivers=Kenai3split, main="adding points every 200m")
points(Kenai3split$lines[[74]])  

## ----fig.width=7.5,fig.height=5------------------------------------------
data(Gulk, KilleyW)
par(mfrow=c(1,2))
plot(x=Gulk, main="Gulkana River")
plot(x=KilleyW, ylim=c(1164500, 1168500), main="Killey River West")

checkbraided(rivers=Gulk, progress=FALSE)
checkbraided(rivers=KilleyW, progress=FALSE)

## ----fig.width=5,fig.height=5--------------------------------------------
Kenai3.subset <- trimriver(trimto=c(22,2,70,30,15,98,96,89,52,3), rivers=Kenai3)
plot(x=Kenai3.subset)

checkbraided(startseg=1, endseg=7, rivers=Kenai3.subset)
checkbraided(startseg=1, endseg=5, rivers=Kenai3.subset)

## ------------------------------------------------------------------------
Killey.dists <- riverdistancelist(startseg=1,endseg=16,startvert=25,endvert=25,
   rivers=KilleyW,reps=1000)
Killey.dists  # 18 routes are detected.

## ----fig.width=7.5,fig.height=5------------------------------------------
par(mfrow=c(1,2))
plot(x=KilleyW, ylim=c(1164500, 1168500), main="shortest route")
riverdistance(startvert=25, endvert=25, path=Killey.dists$routes[[1]], 
              rivers=KilleyW, map=TRUE, add=TRUE)
plot(KilleyW, ylim=c(1164500, 1168500), main="longest route")
riverdistance(startvert=25, endvert=25, path=Killey.dists$routes[[18]], 
              rivers=KilleyW, map=TRUE, add=TRUE)

## ------------------------------------------------------------------------
detectroute(start=1, end=16, rivers=KilleyW)
Killey.dists$routes[[1]]  #calculated above

