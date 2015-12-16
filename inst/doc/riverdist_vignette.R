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

## ----fig.width=5,fig.height=5--------------------------------------------
data(fakefish)
fakefish_riv <- xy2segvert(x=fakefish$x, y=fakefish$y, rivers=Gulk)
head(fakefish_riv)  # a look at what it returned

## ----eval=FALSE----------------------------------------------------------
#  segvert_from_shp <- pointshp2segvert(path=".", layer="MyPointShapefile", rivers=MyRivernetwork)

## ----fig.width=5,fig.height=5--------------------------------------------
zoomtoseg(seg=c(11, 14), rivers=Gulk)
points(fakefish$x, fakefish$y, pch=16, col="red")
riverpoints(seg=fakefish_riv$seg, vert=fakefish_riv$vert, rivers=Gulk, pch=15, col="blue")

## ----fig.width=5,fig.height=5--------------------------------------------
# starting location: segment 7, vertex 49
# ending location: segment 14, vertex 121
riverdistance(startseg=7, startvert=49, endseg=14, endvert=121, rivers=Gulk, map=TRUE)
detectroute(start=7, end=14, rivers=Gulk)

## ----fig.width=7,fig.height=7--------------------------------------------
head(fakefish)
riverdistanceseq(unique=fakefish$fish.id, survey=fakefish$flight, seg=fakefish$seg, 
                   vert=fakefish$vert, rivers=Gulk)

## ------------------------------------------------------------------------
riverdistancematobs(indiv=2, ID=fakefish$fish.id, survey=fakefish$flight,
      seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk, full=FALSE)

## ------------------------------------------------------------------------
riverdistancematobs(indiv=1, ID=fakefish$fish.id, survey=fakefish$flight,
      seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk, full=TRUE)

## ----fig.width=5,fig.height=5--------------------------------------------
# calculating observed minimum home range for all individuals
homerange(unique=fakefish$fish.id, seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk)

# calculating and mapping observed minimum home range for individual 15
homerange(unique=fakefish$fish.id[fakefish$fish.id==15], seg=fakefish$seg[fakefish$fish.id==15],
           vert=fakefish$vert[fakefish$fish.id==15], rivers=Gulk, map=TRUE)

## ------------------------------------------------------------------------
logi1 <- fakefish$flight.date==as.Date("2015-11-25")
riverdistancemat(seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk, logical=logi1)

## ------------------------------------------------------------------------
streamlocs.seg <- c(1,8,11)
streamlocs.vert <- c(50,70,90)
streamlocs.ID <- c("loc A","loc B","loc C")

logi2 <- fakefish$flight.date==as.Date("2015-11-25")

riverdistancetofrom(seg1=streamlocs.seg, vert1=streamlocs.vert, seg2=fakefish$seg, 
                    vert2=fakefish$vert, ID1=streamlocs.ID, logical2=logi2, rivers=Gulk)

## ----fig.width=5,fig.height=5--------------------------------------------
showends(seg=1,rivers=Gulk)
Gulk1 <- setmouth(seg=1, vert=1, rivers=Gulk)

## ----fig.width=7,fig.height=7--------------------------------------------
riverdistance(startseg=6, endseg=4, startvert=250, endvert=250, rivers=Gulk1, map=TRUE)
text(c(859122.4, 872104.1), c(6964127.4,6969741.0), pos=c(3, 4), labels=c("beginning", "end"))
riverdirection(startseg=6, endseg=4, startvert=250, endvert=250, rivers=Gulk1)
upstream(startseg=6, endseg=4, startvert=250, endvert=250, rivers=Gulk1, net=FALSE)
upstream(startseg=6, endseg=4, startvert=250, endvert=250, rivers=Gulk1, net=TRUE)
upstream(startseg=6, endseg=4, startvert=250, endvert=250, rivers=Gulk1, flowconnected=TRUE)

## ------------------------------------------------------------------------
mouthdist(seg=fakefish$seg[1], vert=fakefish$vert[1], rivers=Gulk)
x <- mouthdistobs(unique=fakefish$fish.id, survey=fakefish$flight.date,
    seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk)

head(x)

## ----fig.width=7.5,fig.height=4------------------------------------------
x <- mouthdistobs(unique=fakefish$fish.id, survey=fakefish$flight.date,
    seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk)

par(mfrow=c(1,2))
plotseq(seqobs=x)
plotseq(seqobs=x, type="dotline")
par(mfrow=c(1,1))

## ----fig.width=7,fig.height=7--------------------------------------------
data(abstreams)
plot(x=abstreams)
riverpoints(seg=c(120,131), vert=c(10, 10), rivers=abstreams, pch=15, col=2)

abstreams_nosegroutes <- abstreams
abstreams_nosegroutes$segroutes <- NULL # removing the segment routes for demonstration
abstreams_nosegroutes0 <- abstreams_nosegroutes
abstreams_nosegroutes0$braided <- NA   # removing the braiding information for demonstration

# distance calculation before defining braiding or building routes
tstart <- Sys.time()
riverdistance(startseg=120, startvert=10, endseg=131, endvert=10, rivers=abstreams_nosegroutes0)
tend <- Sys.time()
tend - tstart

# distance calculation before building routes
tstart <- Sys.time()
riverdistance(startseg=120, startvert=10, endseg=131, endvert=10, rivers=abstreams_nosegroutes)
tend <- Sys.time()
tend - tstart

# distance calculation after building routes
abstreams_segroutes <- buildsegroutes(abstreams_nosegroutes)
tstart <- Sys.time()
riverdistance(startseg=120, startvert=10, endseg=131, endvert=10, rivers=abstreams_segroutes)
tend <- Sys.time()
tend - tstart

## ----eval=FALSE----------------------------------------------------------
#  data(abstreams0)  # a messy river network
#  abstreams_fixed <- cleanup(abstreams0)  # fixing many problems

## ----fig.width=7,fig.height=7--------------------------------------------
data(Kenai1)
plot(x=Kenai1)

Kenai1.trim <- trimriver(trim=c(46,32,115,174,169,114,124,142,80), rivers=Kenai1)
plot(x=Kenai1.trim)

## ----fig.width=7,fig.height=7--------------------------------------------
plot(x=Kenai1)

## ----fig.width=5,fig.height=5--------------------------------------------
Kenai1.trim.2 <- trimriver(trimto=c(20,57,118,183,45,162,39,98,19), rivers=Kenai1)
plot(x=Kenai1.trim.2)

## ----fig.width=7,fig.height=7--------------------------------------------
data(Kenai3)
x <- c(174185, 172304, 173803, 176013)
y <- c(1173471, 1173345, 1163638, 1164801)

plot(Kenai3)
points(x, y, pch=15, col=4)
legend(par("usr")[1], par("usr")[4], legend="points to buffer around", pch=15, col=4, cex=.6)

## ----fig.width=5,fig.height=5--------------------------------------------
Kenai3.buf1 <- trimtopoints(x=x, y=y, rivers=Kenai3, method="snap")
plot(x=Kenai3.buf1)
points(x, y, pch=15, col=4)

## ----fig.width=5,fig.height=5--------------------------------------------
Kenai3.buf2 <- trimtopoints(x=x, y=y, rivers=Kenai3, method="snaproute")
plot(x=Kenai3.buf2)
points(x, y, pch=15, col=4)

## ----fig.width=5,fig.height=5--------------------------------------------
Kenai3.buf3 <- trimtopoints(x=x, y=y, rivers=Kenai3, method="buffer", dist=5000)
plot(x=Kenai3.buf3)
points(x, y, pch=15, col=4)

## ----fig.width=5,fig.height=5--------------------------------------------
data(Kenai1)
Kenai1.1 <- dissolve(Kenai1)
Kenai1.1 <- setmouth(seg=63,vert=40,rivers=Kenai1.1)

Kenai1.2 <- removeunconnected(Kenai1.1)
plot(Kenai1.2)

## ----fig.width=7,fig.height=7--------------------------------------------
data(Kenai2)
plot(x=Kenai2)

Kenai2.dissolve <- dissolve(rivers=Kenai2)
plot(x=Kenai2.dissolve)

## ----fig.width=7,fig.height=7--------------------------------------------
data(Koyukuk1)
topologydots(rivers=Koyukuk1)

Koyukuk1.split <- splitsegments(rivers=Koyukuk1)
topologydots(rivers=Koyukuk1.split)

## ----fig.width=5,fig.height=5--------------------------------------------
data(Koyukuk0)
plot(Koyukuk0, ylim=c(1930500,1931500), xlim=c(194900,195100))
topologydots(Koyukuk0, add=TRUE)

Koyukuk0.1 <- connectsegs(connect=21, connectto=20, rivers=Koyukuk0)
plot(Koyukuk0.1,ylim=c(1930500,1931500), xlim=c(194900,195100))
topologydots(Koyukuk0.1, add=TRUE)

## ----eval=FALSE----------------------------------------------------------
#  data(abstreams1)
#  abstreams2 <- removemicrosegs(abstreams1)

## ----fig.width=5,fig.height=5--------------------------------------------
data(Gulk)
plot(x=Gulk)
checkbraided(rivers=Gulk)

data(KilleyW)
plot(x=KilleyW, ylim=c(1164500, 1168500))
checkbraided(rivers=KilleyW)

## ----fig.width=5,fig.height=5--------------------------------------------
Kenai3.subset <- trimriver(trimto=c(18,1,64,27,104,93,91,83,45,2), rivers=Kenai3)
plot(x=Kenai3.subset)

checkbraided(startseg=1, endseg=7, rivers=Kenai3.subset)
checkbraided(startseg=1, endseg=5, rivers=Kenai3.subset)

## ------------------------------------------------------------------------
Killey.dists <- riverdistancelist(startseg=1,endseg=16,startvert=25,endvert=25,
   rivers=KilleyW,reps=1000)
Killey.dists  # 18 routes are detected.

## ----fig.width=5,fig.height=5--------------------------------------------
plot(x=KilleyW, ylim=c(1164500, 1168500))
riverdistance(startvert=25, endvert=25, path=Killey.dists$routes[[1]], 
              rivers=KilleyW, map=TRUE, add=TRUE)
plot(KilleyW, ylim=c(1164500, 1168500))
riverdistance(startvert=25, endvert=25, path=Killey.dists$routes[[18]], 
              rivers=KilleyW, map=TRUE, add=TRUE)

