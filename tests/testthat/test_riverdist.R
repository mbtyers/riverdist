flipsegs <- function(rivers,whichflip="all") {
  length <- length(rivers$lines)
  if(whichflip=="sample") flip <- sample(c(T,F),size=length,replace=T)
  if(whichflip=="half") flip <- c(rep(T,(floor(length/2))),rep(F,(length-floor(length/2))))
  if(whichflip=="all") flip <- rep(T,length)
  for(i in 1:length(rivers$lines)) {
    if(flip[i]) {
      rivers$lines[[i]] <- rivers$lines[[i]][nrow(rivers$lines[[i]]):1,]
    }
  }
  lines <- rivers$lines
  tolerance <- rivers$tolerance
  connections <- rivers$connections
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
  rivers$connections <- connections
  if(!is.na(rivers$mouth$mouth.seg)) {
    if(flip[rivers$mouth$mouth.seg]) rivers$mouth$mouth.vert <- ifelse(rivers$mouth$mouth.vert==1,nrow(rivers$lines[[rivers$mouth$mouth.seg]]),1)
  }
  rivers <- addcumuldist(rivers)
  if(!is.null(rivers$segroutes)) rivers <- buildsegroutes(rivers)
  if(!is.null(rivers$distlookup)) rivers <- buildlookup(rivers)
  
  rivers <- update_sf(rivers)
  
  return(rivers)
}

flipdataverts <- function(seg,vert,rivers) {
  seglengths <- sapply(rivers$lines,nrow)
  newverts <- seglengths[seg]-vert+1
  return(newverts)
}

Gulk_flip <- flipsegs(Gulk)
fakefish_flip_verts <- flipdataverts(fakefish$seg,fakefish$vert,Gulk)

data(Gulk)
test_that("distance",{
  expect_equal(riverdistance(startseg=7, startvert=49, endseg=14, endvert=121, rivers=Gulk, algorithm="Dijkstra"), 155435.2, tolerance=0.001)
  expect_equal(riverdistance(startseg=7, startvert=49, endseg=14, endvert=121, rivers=Gulk, algorithm="sequential"), 155435.2, tolerance=0.001)
  expect_equal(riverdistance(startseg=7, startvert=49, endseg=14, endvert=121, rivers=Gulk, algorithm="segroutes"), 155435.2, tolerance=0.001)
  expect_equal(riverdistance(startseg=7, startvert=49, endseg=14, endvert=121, rivers=Gulk), 155435.2, tolerance=0.001)
  expect_equal(riverdistance(startseg=1, startvert=49, endseg=14, endvert=27, rivers=Gulk), 155105.9, tolerance=0.001)
  expect_error(riverdistance(startseg=77, startvert=49, endseg=14, endvert=121, rivers=Gulk))
  expect_error(riverdistance(startseg=7, startvert=149, endseg=14, endvert=121, rivers=Gulk))
  expect_equal(riverdistance(startseg=1, endseg=3, startvert=20, endvert=20, rivers=Gulk, algorithm="Dijkstra"), 76375.44, tolerance=0.01)  #end-beginning
  expect_equal(riverdistance(startseg=3, endseg=1, startvert=20, endvert=20, rivers=Gulk, algorithm="Dijkstra"), 76375.44, tolerance=0.01)  #beginning-end
  expect_equal(riverdistance(startseg=3, endseg=4, startvert=20, endvert=20, rivers=Gulk, algorithm="Dijkstra"), 2996.144, tolerance=0.01)  #beginning-beginning
  expect_equal(riverdistance(startseg=1, endseg=3, startvert=20, endvert=20, rivers=Gulk, algorithm="sequential"), 76375.44, tolerance=0.01)  #end-beginning
  expect_equal(riverdistance(startseg=3, endseg=1, startvert=20, endvert=20, rivers=Gulk, algorithm="sequential"), 76375.44, tolerance=0.01)  #beginning-end
  expect_equal(riverdistance(startseg=3, endseg=4, startvert=20, endvert=20, rivers=Gulk, algorithm="sequential"), 2996.144, tolerance=0.01)  #beginning-beginning
  expect_equal(riverdistance(startseg=1, endseg=3, startvert=20, endvert=20, rivers=Gulk, algorithm="segroutes"), 76375.44, tolerance=0.01)  #end-beginning
  expect_equal(riverdistance(startseg=3, endseg=1, startvert=20, endvert=20, rivers=Gulk, algorithm="segroutes"), 76375.44, tolerance=0.01)  #beginning-end
  expect_equal(riverdistance(startseg=3, endseg=4, startvert=20, endvert=20, rivers=Gulk, algorithm="segroutes"), 2996.144, tolerance=0.01)  #beginning-beginning
  expect_equal(riverdistance(startseg=1, endseg=3, startvert=20, endvert=20, rivers=Gulk), 76375.44, tolerance=0.01)  #end-beginning
  expect_equal(riverdistance(startseg=3, endseg=1, startvert=20, endvert=20, rivers=Gulk), 76375.44, tolerance=0.01)  #beginning-end
  expect_equal(riverdistance(startseg=3, endseg=4, startvert=20, endvert=20, rivers=Gulk), 2996.144, tolerance=0.01)  #beginning-beginning
  expect_equal(sum(riverdistancemat(seg=fakefish$seg,vert=fakefish$vert,rivers=Gulk)),638495319,tolerance=0.1)
  expect_equal(sum(riverdistancemat(seg=fakefish$seg,vert=fakefish$vert,rivers=Gulk,algorithm="segroutes")),638495319,tolerance=0.1)
  expect_equal(sum(riverdistancemat(seg=fakefish$seg,vert=fakefish$vert,rivers=Gulk,algorithm="Dijkstra")),638495319,tolerance=0.1)
  expect_equal(sum(riverdistancemat(seg=fakefish$seg,vert=fakefish_flip_verts,rivers=Gulk_flip)),638495319,tolerance=0.1)
  expect_equal(sum(riverdistancemat(seg=fakefish$seg,vert=fakefish_flip_verts,rivers=Gulk_flip,algorithm="segroutes")),638495319,tolerance=0.1)
  expect_equal(sum(riverdistancemat(seg=fakefish$seg,vert=fakefish_flip_verts,rivers=Gulk_flip,algorithm="Dijkstra")),638495319,tolerance=0.1)
  expect_equal(riverdistance(startseg=1,startvert=nrow(Gulk$lines[[1]]),endseg=3,endvert=99,rivers=Gulk),Gulk$cumuldist[[3]][99],tolerance=0.0001)
  expect_equal(riverdistance(startseg=3,startvert=1,endseg=3,endvert=99,rivers=Gulk),Gulk$cumuldist[[3]][99],tolerance=0.0001)
  expect_equal(riverdistance(endseg=1,endvert=nrow(Gulk$lines[[1]]),startseg=3,startvert=99,rivers=Gulk),Gulk$cumuldist[[3]][99],tolerance=0.0001)
  expect_equal(riverdistance(endseg=3,endvert=1,startseg=3,startvert=99,rivers=Gulk),Gulk$cumuldist[[3]][99],tolerance=0.0001)
  expect_equal(riverdistance(startseg=1,startvert=nrow(Gulk$lines[[1]]),endseg=3,endvert=99,rivers=Gulk,algorithm="segroutes"),Gulk$cumuldist[[3]][99],tolerance=0.0001)
  expect_equal(riverdistance(startseg=3,startvert=1,endseg=3,endvert=99,rivers=Gulk,algorithm="segroutes"),Gulk$cumuldist[[3]][99],tolerance=0.0001)
  expect_equal(riverdistance(endseg=1,endvert=nrow(Gulk$lines[[1]]),startseg=3,startvert=99,rivers=Gulk,algorithm="segroutes"),Gulk$cumuldist[[3]][99],tolerance=0.0001)
  expect_equal(riverdistance(endseg=3,endvert=1,startseg=3,startvert=99,rivers=Gulk,algorithm="segroutes"),Gulk$cumuldist[[3]][99],tolerance=0.0001)
  expect_equal(riverdistance(startseg=1,startvert=nrow(Gulk$lines[[1]]),endseg=3,endvert=99,rivers=Gulk,algorithm="Dijkstra"),Gulk$cumuldist[[3]][99],tolerance=0.0001)
  expect_equal(riverdistance(startseg=3,startvert=1,endseg=3,endvert=99,rivers=Gulk,algorithm="Dijkstra"),Gulk$cumuldist[[3]][99],tolerance=0.0001)
  expect_equal(riverdistance(endseg=1,endvert=nrow(Gulk$lines[[1]]),startseg=3,startvert=99,rivers=Gulk,algorithm="Dijkstra"),Gulk$cumuldist[[3]][99],tolerance=0.0001)
  expect_equal(riverdistance(endseg=3,endvert=1,startseg=3,startvert=99,rivers=Gulk,algorithm="Dijkstra"),Gulk$cumuldist[[3]][99],tolerance=0.0001)
})

data(fakefish)
fakefish.riv <- xy2segvert(fakefish$x, fakefish$y, rivers=Gulk)
test_that("xy2segvert",{
  expect_equal(fakefish.riv$seg,fakefish$seg)       
  expect_equal(fakefish.riv$vert,fakefish$vert)
})

Gulk1 <- buildsegroutes(Gulk)
Gulk_flip1 <- buildsegroutes(Gulk_flip)
data(abstreams)
abstreams_nosegroutes <- abstreams
abstreams_nosegroutes$segroutes <- NULL
abstreams_nosegroutes$distlookup <- NULL
abstreams_nosegroutes1 <- buildsegroutes(abstreams_nosegroutes)
test_that("buildsegroutes",{
  expect_equal(unlist(Gulk1$segroutes),c(1,1,3,5,2,1,3,1,4,1,3,5,1,3,6,1,3,6,7,1,3,6,8,1,4,9,1,4,10,1,4,10,11,1,4,10,12,1,4,10,11,13,1,4,10,11,14))
  expect_equal(unlist(Gulk_flip1$segroutes),c(1,1,3,5,2,1,3,1,4,1,3,5,1,3,6,1,3,6,7,1,3,6,8,1,4,9,1,4,10,1,4,10,11,1,4,10,12,1,4,10,11,13,1,4,10,11,14))
  expect_equal(abstreams_nosegroutes1$segroutes,abstreams$segroutes)
  expect_equal(abstreams_nosegroutes1$distlookup,abstreams$distlookup)
  expect_equal(Gulk$cumuldist,addcumuldist(Gulk)$cumuldist)
})

data(Kenai3)
Kenai3.1 <- setmouth(seg=68,vert=40,rivers=Kenai3)
Kenai3.subset <- suppressWarnings(trimriver(trimto=c(22,2,70,30,15,98,96,89,52,3), rivers=Kenai3))
test_that("checkbraided",{
  expect_false(checkbraidedTF(rivers=Gulk, toreturn="logical"))
  expect_true(checkbraidedTF(rivers=Kenai3.1, toreturn="logical"))
  })

data(abstreams)
test_that("detectroute",{
  expect_equal(detectroute(start=1,end=9,rivers=Gulk),c(1,4,9))
  expect_error(detectroute(start=1,end=99,rivers=Gulk))
  expect_equal(detectroute(start=120,end=111,rivers=abstreams),c(120,103,106,109,112,116,124,132,134,133,135,142,153,152,144,136,127,115,114,107,108,111))
  expect_equal(detectroute(start=120,end=111,rivers=abstreams,algorithm="Dijkstra"),c(120,103,106,109,112,116,124,132,134,133,135,142,153,152,144,136,127,115,114,107,108,111))
  expect_equal(detectroute(start=120,end=111,rivers=abstreams,algorithm="sequential"),c(120,103,106,109,112,116,124,132,134,133,135,142,153,152,144,136,127,115,114,107,108,111))
  expect_equal(detectroute(start=116,end=14,rivers=abstreams),detectroute(start=116,end=14,rivers=abstreams_nosegroutes))
})

data(Kenai2)
Kenai3flip <- flipsegs(Kenai3)
Kenai2flipdis <- dissolve(flipsegs(Kenai2))
Kenai3flip$sp <- NULL
Kenai2flipdis$sp <- NULL
test_that("dissolve",{
  expect_equal(dissolve(Kenai2),Kenai3)   
  expect_equal(length(dissolve(Gulk)$segroutes),13,tolerance=0.001)
  expect_equal(sum(dissolve(Gulk)$distlookup$middist),8360513,tolerance=1)
  expect_equal(sum(dissolve(Gulk)$distlookup$endtop,na.rm=T),126,tolerance=0.001)
  expect_equal(sum(dissolve(Gulk)$distlookup$starttop,na.rm=T),126,tolerance=0.001)
  expect_equal(Kenai3flip,Kenai2flipdis)
})

hr <- homerange(unique=fakefish$fish.id,seg=fakefish$seg,vert=fakefish$vert,survey=fakefish$flight,rivers=Gulk)
hr_flip <- homerange(unique=fakefish$fish.id,seg=fakefish$seg,vert=fakefish_flip_verts,survey=fakefish$flight,rivers=Gulk_flip)
hr_flipflip <- hr_flip
for(i in 1:length(hr$subseg_n)) {
  for(j in 1:length(hr$subseg_length)) {
    hr_flipflip$subseg_n[[i]][[j]] <- hr_flipflip$subseg_n[[i]][[j]][length(hr_flipflip$subseg_n[[i]][[j]]):1]
  }
}
for(j in 1:length(hr$subseg_length)) {
  hr_flipflip$subseg_length[[j]] <- hr_flipflip$subseg_length[[j]][length(hr_flipflip$subseg_length[[j]]):1]
}
hr_overlap <- homerangeoverlap(hr)
test_that("homerange",{
  expect_equal(hr$ranges[,1], c(1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))
  expect_equal(hr$ranges[,2], c(165698.89,94833.30,232892.91,141143.68,138765.21,145436.14,113141.15,113860.33,101682.66,156097.77,97081.89,177000.52,179146.30,149433.33,139167.26,179123.34,70523.57,151396.99,174099.34),tolerance=0.01)
  expect_equal(names(hr$ranges),c("ID","range"))
  expect_error(homerange(unique=1:10,seg=fakefish$seg,vert=fakefish$vert,rivers=Gulk))
  expect_equal(hr$ranges,hr_flip$ranges)
  expect_equal(hr$subseg_n,hr_flipflip$subseg_n)
  expect_equal(hr$subseg_length,hr_flipflip$subseg_length)
  expect_equal(sum(unlist(hr$subseg_n)),71171)
  expect_equal(sum(unlist(hr$subseg_length)),371439.8,tolerance=0.1)
  expect_equal(sum(hr_overlap$either),69995271,tolerance=0.1)
  expect_equal(sum(hr_overlap$both),33384663,tolerance=0.1)
  expect_equal(sum(hr_overlap$prop_both),177.729,tolerance=0.1)
})

test_that("isflowconnected",{
  expect_true(isflowconnected(seg1=130,seg2=158,rivers=abstreams))
  expect_true(isflowconnected(seg1=130,seg2=158,rivers=abstreams_nosegroutes))
  expect_false(isflowconnected(seg1=130,seg2=104,rivers=abstreams))
  expect_false(isflowconnected(seg1=130,seg2=104,rivers=abstreams_nosegroutes))
})

test_that("mouthdist",{
  expect_equal(mouthdist(4,19,abstreams),92745.93,tolerance=0.001)
  expect_equal(mouthdist(4,19,abstreams_nosegroutes),92745.93,tolerance=0.001)
  expect_error(mouthdist(4,19,Kenai3))
})

test_that("direction",{
  expect_equal(riverdirection(startseg=7, startvert=49, endseg=14, endvert=121, rivers=Gulk), "up")
  expect_equal(riverdirection(startseg=12, startvert=49, endseg=3, endvert=27, rivers=Gulk), "down")
  expect_equal(riverdirection(startseg=12, startvert=49, endseg=12, endvert=49, rivers=Gulk), "0")
  expect_true(is.na(riverdirection(startseg=7, startvert=49, endseg=14, endvert=121, flowconnected=T, rivers=Gulk)))
  expect_false(is.na(riverdirection(startseg=7, startvert=49, endseg=1, endvert=121, flowconnected=T, rivers=Gulk)))
  expect_error(riverdirection(startseg=77, startvert=49, endseg=14, endvert=121, rivers=Gulk))
})

test_that("upstream",{
  expect_equal(upstream(startseg=7, startvert=49, endseg=14, endvert=121, rivers=Gulk), 155435.2, tolerance=0.001)
  expect_equal(upstream(startseg=12, startvert=49, endseg=3, endvert=27, rivers=Gulk), -61647.23, tolerance=0.001)
  expect_equal(upstream(startseg=7, startvert=49, endseg=14, endvert=121, rivers=Gulk, net=T), 18764.18, tolerance=0.001)
  expect_equal(upstream(startseg=12, startvert=49, endseg=3, endvert=27, rivers=Gulk, net=T), -57735.08, tolerance=0.001)
  expect_equal(upstream(startseg=12, startvert=49, endseg=12, endvert=49, rivers=Gulk), 0, tolerance=0.001)
  expect_true(is.na(upstream(startseg=7, startvert=49, endseg=14, endvert=121, flowconnected=T, rivers=Gulk)))
  expect_false(is.na(upstream(startseg=7, startvert=49, endseg=1, endvert=121, flowconnected=T, rivers=Gulk)))
  expect_error(upstream(startseg=77, startvert=49, endseg=14, endvert=121, rivers=Gulk))
})

dm <- riverdistancemat(fakefish$seg,fakefish$vert, logical=(fakefish$flight.date==as.Date("2015-11-25")), rivers=Gulk)
um <- upstreammat(fakefish$seg,fakefish$vert, logical=(fakefish$flight.date==as.Date("2015-11-25")), rivers=Gulk)
dirm <- riverdirectionmat(fakefish$seg,fakefish$vert, logical=(fakefish$flight.date==as.Date("2015-11-25")), rivers=Gulk)
test_that("mats",{
  expect_equal(sum(dm),5027666,tolerance=0.001)
  expect_equal(sum(um[,1]),-583799.3,tolerance=0.001)
  expect_equal(dirm[2,1],"down")
  expect_equal(dirm[1,2],"up")
  expect_equal(dirm[1,1],"0")
  expect_equal(row.names(dm),c("91",  "92",  "93",  "94",  "95",  "96",  "97",  "98",  "99",  "100"))
  expect_equal(row.names(um),c("91",  "92",  "93",  "94",  "95",  "96",  "97",  "98",  "99",  "100"))
  expect_equal(row.names(dirm),c("91",  "92",  "93",  "94",  "95",  "96",  "97",  "98",  "99",  "100"))
})

ds <- riverdistanceseq(unique=fakefish$fish.id, survey=fakefish$flight, seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk)
us <- upstreamseq(unique=fakefish$fish.id, survey=fakefish$flight, seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk)
dirs <- riverdirectionseq(unique=fakefish$fish.id, survey=fakefish$flight, seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk)
test_that("seqs",{
  expect_equal(ds[1,8],54220.046,tolerance=0.001)
  expect_equal(sum(unlist(ds)[!is.na(unlist(ds))]),3145402,tolerance=.1)
  expect_equal(sum(unlist(us)[!is.na(unlist(us))]),49838.52,tolerance=.01)
  expect_equal(us[1,8],-54220.046,tolerance=0.001)
  expect_equal(as.character(dirs[1,8]),"down")
  expect_true(is.na(ds[1,2]))
  expect_true(is.na(us[1,2]))
  expect_true(is.na(dirs[1,2]))
  expect_equal(names(ds),c("1 to 2","2 to 3","3 to 4","4 to 5","5 to 6","6 to 7","7 to 8","8 to 9","9 to 10"))
  expect_equal(names(us),c("1 to 2","2 to 3","3 to 4","4 to 5","5 to 6","6 to 7","7 to 8","8 to 9","9 to 10"))
  expect_equal(names(dirs),c("1 to 2","2 to 3","3 to 4","4 to 5","5 to 6","6 to 7","7 to 8","8 to 9","9 to 10"))
  expect_equal(row.names(ds),c("1",  "3",  "4",  "6",  "7",  "8",  "9",  "10", "11", "13", "14", "15", "16", "17", "18", "19", "20"))
  expect_equal(row.names(us),c("1",  "3",  "4",  "6",  "7",  "8",  "9",  "10", "11", "13", "14", "15", "16", "17", "18", "19", "20"))
  expect_equal(row.names(dirs),c("1",  "3",  "4",  "6",  "7",  "8",  "9",  "10", "11", "13", "14", "15", "16", "17", "18", "19", "20"))
})

streamlocs.seg <- c(1,8,11)
streamlocs.vert <- c(50,70,90)
streamlocs.ID <- c("loc A","loc B","loc C")
logi2 <- fakefish$flight.date==as.Date("2015-11-25")
dt <- riverdistancetofrom(seg1=streamlocs.seg, vert1=streamlocs.vert, ID1=streamlocs.ID, seg2=fakefish$seg, vert2=fakefish$vert, logical2=logi2, rivers=Gulk1)
ut <- upstreamtofrom(seg1=streamlocs.seg, vert1=streamlocs.vert, ID1=streamlocs.ID, seg2=fakefish$seg, vert2=fakefish$vert, logical2=logi2, rivers=Gulk1)
dirt <- riverdirectiontofrom(seg1=streamlocs.seg, vert1=streamlocs.vert, ID1=streamlocs.ID, seg2=fakefish$seg, vert2=fakefish$vert, logical2=logi2, rivers=Gulk1)
test_that("tofrom",{
  expect_equal(sum(dt),2446906,tolerance=0.001)
  expect_equal(sum(ut),-1205628,tolerance=0.001)
  expect_equal(dirt[1,4],"up")
  expect_equal(dirt[2,4],"down")
  expect_equal(dimnames(dt)[[2]],c("91",  "92",  "93",  "94",  "95",  "96",  "97",  "98",  "99",  "100"))
  expect_equal(dimnames(ut)[[2]],c("91",  "92",  "93",  "94",  "95",  "96",  "97",  "98",  "99",  "100"))
  expect_equal(dimnames(dirt)[[2]],c("91",  "92",  "93",  "94",  "95",  "96",  "97",  "98",  "99",  "100"))
  expect_equal(row.names(dt),c("loc A","loc B","loc C"))
  expect_equal(row.names(ut),c("loc A","loc B","loc C"))
  expect_equal(row.names(dirt),c("loc A","loc B","loc C"))
})

data(KilleyW)
Killey.dists <- riverdistancelist(startseg=1,endseg=16,startvert=25,endvert=25,rivers=KilleyW,reps=1000)
test_that("routedistlist",{
  expect_equal(Killey.dists$routes[[1]],c(1,2,4,15,16))
  expect_equal(Killey.dists$routes[[3]],c(1,3,4,15,16))
  expect_equal(sum(Killey.dists$distances),115582.5,tolerance=0.001)
})

KilleyW <- setmouth(seg=1,vert=288,rivers=KilleyW)
test_that("set",{
  expect_equal(KilleyW$mouth$mouth.seg,1)
  expect_equal(KilleyW$mouth$mouth.vert,288)
})

Killey2 <- sequenceverts(rivers=KilleyW)
test_that("sequence",{
  expect_equal(Killey2$lines[[1]][1,1],184621.8,tolerance=0.001)
  expect_equal(KilleyW$lines[[1]][1,1],183649.9,tolerance=0.001)
})

data(Koyukuk1,Koyukuk2)
Koyukuk1a <- splitsegments(rivers=Koyukuk1)
Koyukuk1b <- splitsegments(rivers=Koyukuk1, splitthese=c(7,7,7), splitthemat=c(14,5,12))
Koyukuk1c <- splitsegments(rivers=Koyukuk1, splitthese=c(7,7,7), splitthemat=c(14,5,12), append=T)
Koyukuk0.2 <- connectsegs(connect=c(20,21,22), connectto=c(21,22,23), 
                          nearestvert=c(FALSE,FALSE,TRUE), rivers=Koyukuk0)
test_that("splitsegments",{
  expect_equal(Koyukuk1a,Koyukuk2)
  expect_equal(length(Koyukuk1b$lines), 20)
  expect_equal(sum(Koyukuk1b$lines[[20]]), 505784711)
  expect_equal(length(Koyukuk1c$lines), 20)
  expect_equal(sum(Koyukuk1c$lines[[20]]), 673547226)
  expect_equal(Koyukuk0.2$connections[20,21],2)
  expect_equal(Koyukuk0.2$connections[21,22],2)
  expect_equal(Koyukuk0.2$connections[22,23],2)
  expect_equal(Koyukuk0.2$connections[23,27],3)
  expect_equal(Koyukuk0.2$connections[22,27],1)
})

Gulk3 <- Gulk
Gulk3$segroutes <- NULL             
Gulk3$distlookup <- NULL
Gulk.trim <- trimriver(trim=1:4,rivers=Gulk3)
Gulk.trimto <- trimriver(trimto=1:4,rivers=Gulk3)
data(Koyukuk0)
Koyukuk0a <- trimriver(trimto=c(1,2,9,10,17:23),rivers=Koyukuk0)
test_that("trimriver",{
  expect_equal(Gulk.trimto$lines,Gulk$lines[1:4])
  expect_equal(Gulk.trimto$connections,Gulk$connections[1:4,1:4])
  expect_equal(Gulk.trimto$lengths,Gulk$lengths[1:4])
  expect_equal(Gulk.trimto$names,Gulk$names[1:4])
  # expect_equal(Gulk.trimto$lineID,Gulk$lineID[1:4,])
  expect_equal(Gulk.trimto$sequenced,Gulk$sequenced)
  expect_equal(Gulk.trimto$tolerance,Gulk$tolerance)
  expect_equal(Gulk.trimto$units,Gulk$units)
  expect_equal(Gulk.trimto$mouth$mouth.seg,1)
  expect_equal(Gulk.trim$lines,Gulk$lines[5:14])
  expect_equal(Gulk.trim$connections,Gulk$connections[5:14,5:14])
  expect_equal(Gulk.trim$lengths,Gulk$lengths[5:14])
  expect_equal(Gulk.trim$names,Gulk$names[5:14])
  expect_equal(Gulk.trim$sequenced,Gulk$sequenced)
  expect_equal(Gulk.trim$tolerance,Gulk$tolerance)
  expect_equal(Gulk.trim$units,Gulk$units)
  expect_true(is.na(Gulk.trim$mouth$mouth.seg))
  expect_equal(length(Koyukuk0a$lines),11)
  expect_equal(dim(Koyukuk0a$connections),c(11,11))
})

x <- c(174185, 172304, 173803, 176013)
y <- c(1173471, 1173345, 1163638, 1164801)
Kenai3 <- setmouth(seg=68,vert=40,rivers=Kenai3)
Kenai3.buf1 <- trimtopoints(x=x, y=y, method="snap", rivers=Kenai3)
Kenai3.buf2 <- trimtopoints(x=x, y=y, method="snaproute", rivers=Kenai3)
Kenai3.buf3 <- trimtopoints(x=x, y=y, method="buffer", dist=5000, rivers=Kenai3)
test_that("trimtopoints",{
  expect_equal(length(Kenai3.buf1$lines),2)
  expect_equal(length(Kenai3.buf1$lengths),2)
  expect_equal(length(Kenai3.buf1$names),2)
  expect_equal(dim(Kenai3.buf1$connections),c(2,2))
  expect_equal(length(Kenai3.buf2$lines),6)
  expect_equal(length(Kenai3.buf2$lengths),6)
  expect_equal(length(Kenai3.buf2$names),6)
  expect_equal(dim(Kenai3.buf2$connections),c(6,6))
  expect_equal(length(Kenai3.buf3$lines),26)
  expect_equal(length(Kenai3.buf3$lengths),26)
  expect_equal(length(Kenai3.buf3$names),26)
  expect_equal(dim(Kenai3.buf3$connections),c(26,26))
  expect_true(is.na(Kenai3.buf1$mouth$mouth.seg))
  expect_true(is.na(Kenai3.buf2$mouth$mouth.seg))
  expect_equal(Kenai3.buf3$mouth$mouth.seg,20)
  expect_equal(Kenai3.buf3$mouth$mouth.vert,40)
})

data(Kenai1)
Kenai1a <- dissolve(Kenai1)
Kenai1a$mouth$mouth.seg <- 71
Kenai1a$mouth$mouth.vert <- 40

segs <- c(38,71,89,12)
verts <- c(1,1,1,1)

test_that("stopiferror, flowconnected",{  
  expect_error(riverdistance(startseg=segs[1],endseg=segs[2],startvert=verts[1],endvert=verts[2],rivers=Kenai1a))
  expect_true(is.na(riverdistance(startseg=segs[1],endseg=segs[2],startvert=verts[1],endvert=verts[2],rivers=Kenai1a,stopiferror=F)))
  expect_equal(riverdistance(startseg=segs[3],endseg=segs[2],startvert=verts[3],endvert=verts[2],rivers=Kenai1a,stopiferror=F),2648.679,tolerance=0.001)
  expect_error(riverdirection(startseg=segs[1],endseg=segs[2],startvert=verts[1],endvert=verts[2],rivers=Kenai1a))
  expect_true(is.na(riverdirection(startseg=segs[3],endseg=segs[4],startvert=verts[3],endvert=verts[4],rivers=Kenai1a,flowconnected=T)))
  expect_true(is.na(riverdirection(startseg=segs[1],endseg=segs[2],startvert=verts[1],endvert=verts[2],rivers=Kenai1a,stopiferror=F)))
  expect_equal(riverdirection(startseg=segs[2],endseg=segs[3],startvert=verts[2],endvert=verts[3],rivers=Kenai1a,stopiferror=F,flowconnected=T),"up") 
  expect_error(upstream(startseg=segs[1],endseg=segs[2],startvert=verts[1],endvert=verts[2],rivers=Kenai1a))
  expect_true(is.na(upstream(startseg=segs[3],endseg=segs[4],startvert=verts[3],endvert=verts[4],rivers=Kenai1a,flowconnected=T)))
  expect_true(is.na(upstream(startseg=segs[1],endseg=segs[2],startvert=verts[1],endvert=verts[2],rivers=Kenai1a,stopiferror=F)))
  expect_equal(upstream(startseg=segs[2],endseg=segs[3],startvert=verts[2],endvert=verts[3],rivers=Kenai1a,stopiferror=F,flowconnected=T),2648.679,tolerance=0.001) 
})

data(abstreams0)
Gulk <- setmouth(seg=1,vert=1,rivers=Gulk)
Gulk1 <- trimriver(trim=10,rivers=Gulk3)
Gulk2 <- removeunconnected(Gulk1)
test_that("cleanup funcs",{
  expect_equal(length(removeduplicates(abstreams0)$lines),202)
  expect_equal(length(removemicrosegs(abstreams0)$lines),179)
  expect_equal(Gulk2,trimriver(trimto=1:9,rivers=Gulk3))
})

filepath <- system.file("extdata", package="riverdist")
sf <- sf::read_sf(dsn = filepath, layer = "Gulk_UTM5")
ptshp <- pointshp2segvert(path=filepath, layer="fakefish_UTM5", rivers=Gulk)
Gulktest <- line2network(path=filepath, layer="Gulk_UTM5")
Gulktest_sf <- line2network(path=filepath, layer="Gulk_UTM5")
test_that("line2network and pointshp2segvert works", {
  expect_equal(length(line2network(path=filepath, layer="Gulk_UTM5")$lines),14)
  expect_equal(length(line2network(sf=sf)$lines), 14)
  expect_equal(dim(ptshp),c(100,8))
  expect_equal(sum(ptshp[,1:2]),27095)
  expect_equal(sum(unlist(Gulktest$lines)), sum(unlist(Gulk$lines)))
  expect_true(isTRUE(all.equal(Gulktest$connections, Gulk$connections)))
  expect_true(isTRUE(all.equal(Gulktest$lengths, Gulk$lengths)))
  expect_equal(sum(unlist(Gulktest_sf$lines)), sum(unlist(Gulk$lines)))
  expect_true(isTRUE(all.equal(Gulktest_sf$connections, Gulk$connections)))
  expect_true(isTRUE(all.equal(Gulktest_sf$lengths, Gulk$lengths)))
  expect_silent(plot(Gulktest$sf_current))
  expect_true(inherits(Gulktest$sf_current, "sf"))
  expect_true(inherits(Gulktest$sf_current, "data.frame"))
  expect_equal(length(Gulktest$sf_current$geometry[[1]]), 14)
  expect_equal(as.character(sf::st_geometry_type(Gulktest$sf_current$geometry[[1]])), "MULTILINESTRING")
  expect_equal(dim(Gulktest$sf_current$geometry[[1]][[1]]), c(812, 2))
  expect_equal(Gulktest, Gulktest_sf)
  expect_equal(Gulktest$lines, Gulk$lines)
  expect_equal(Gulktest$lengths, Gulk$lengths)
  expect_equal(sum(unlist(Gulktest$sf_current$geometry)), sum(unlist(Gulk$lines)))
  expect_equal(sum(unlist(Gulktest_sf$sf_current$geometry)), sum(unlist(Gulk$lines)))
}) 

test_that("matbysurvey", {
  expect_equal(dim(riverdistancematbysurvey(indiv=1, unique=fakefish$fish.id, survey=fakefish$flight,
                                       seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk, full=FALSE)),c(7,7))
  expect_equal(dim(riverdistancematbysurvey(indiv=1, unique=fakefish$fish.id, survey=fakefish$flight,
                                       seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk, full=TRUE)),c(10,10))
  expect_equal(sum(riverdistancematbysurvey(indiv=1, unique=fakefish$fish.id, survey=fakefish$flight,
                                       seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk, full=FALSE),na.rm=TRUE),2694810,tolerance=0.001)
}) 

fakefish_sub <- subset(fakefish,vert<40)
fakesubdens <- makeriverdensity(seg=fakefish_sub$seg,vert=fakefish_sub$vert,survey=fakefish_sub$flight.date,rivers=Gulk)
test_that("riverdensity", {
  expect_equal(length(fakesubdens$densities),7)
  expect_equal(length(fakesubdens$densities[[1]]),14)
  expect_equal(sum(unlist(fakesubdens$densities)),0.04737915,tolerance=0.000001)
})

test_that("addverts", {
  expect_equal(dim(Gulk$lines[[1]]),c(812,2))
  expect_equal(dim(addverts(rivers=Gulk,mindist=100)$lines[[1]]),c(1130,2))
})

asdf<-kfunc(seg=fakefish$seg, vert=fakefish$vert, rivers=Gulk, survey=fakefish$flight,envreps=100, maxdist=200000, returnoutput=TRUE)   
test_that("kfunc", {
  expect_equal(length(asdf),4,tolerance=0.001)
  expect_equal(length(asdf$lines),10,tolerance=0.001)
  expect_equal(length(asdf$env_low),10,tolerance=0.001)
  expect_equal(length(asdf$env_high),10,tolerance=0.001)
  expect_equal(length(asdf$dists),100,tolerance=0.001)
  expect_equal(sum(unlist(asdf$lines)),73357.78,tolerance=0.001)
  expect_equal(sum(asdf$dists),10000000,tolerance=0.001)
})

K2 <- trimriver(trimto=c(2,30,70,15),rivers=Kenai3)
K2f <- flipsegs(K2,whichflip="half")
K2l <- buildlookup(K2)
K2fl <- buildlookup(K2f)
test_that("connections 5 and 6", {
  expect_equal(K2$connections[2,3],5,tolerance=0.001)
  expect_equal(K2l$connections[3,2],5,tolerance=0.001)
  expect_equal(riverdistance(startseg=2, endseg=3, startvert=5, endvert=5, rivers=K2), 486.0265, tolerance=0.001)
  expect_equal(riverdistance(startseg=2, endseg=3, startvert=5, endvert=5, rivers=K2l), 486.0265, tolerance=0.001)
  expect_equal(riverdistance(startseg=2, endseg=3, startvert=1, endvert=1, rivers=K2), 0, tolerance=0.001)
  expect_equal(riverdistance(startseg=2, endseg=3, startvert=1, endvert=1, rivers=K2l), 0, tolerance=0.001)
  expect_equal(riverdistance(startseg=2, endseg=3, startvert=24, endvert=24, rivers=K2), 113.4867, tolerance=0.001)
  expect_equal(riverdistance(startseg=2, endseg=3, startvert=24, endvert=24, rivers=K2l), 113.4867, tolerance=0.001)
  expect_equal(riverdistance(startseg=2, endseg=3, startvert=25, endvert=25, rivers=K2), 0, tolerance=0.001)
  expect_equal(riverdistance(startseg=2, endseg=3, startvert=25, endvert=25, rivers=K2l), 0, tolerance=0.001)
  expect_equal(K2f$connections[2,3],6,tolerance=0.001)
  expect_equal(K2fl$connections[3,2],6,tolerance=0.001)
  expect_equal(riverdistance(startseg=2, endseg=3, startvert=21, endvert=5, rivers=K2f), 486.0265, tolerance=0.001)
  expect_equal(riverdistance(startseg=2, endseg=3, startvert=21, endvert=5, rivers=K2fl), 486.0265, tolerance=0.001)
  expect_equal(riverdistance(startseg=2, endseg=3, startvert=25, endvert=1, rivers=K2f), 0, tolerance=0.001)
  expect_equal(riverdistance(startseg=2, endseg=3, startvert=25, endvert=1, rivers=K2fl), 0, tolerance=0.001)
  expect_equal(riverdistance(startseg=2, endseg=3, startvert=2, endvert=24, rivers=K2f), 113.4867, tolerance=0.001)
  expect_equal(riverdistance(startseg=2, endseg=3, startvert=2, endvert=24, rivers=K2fl), 113.4867, tolerance=0.001)
  expect_equal(riverdistance(startseg=2, endseg=3, startvert=1, endvert=25, rivers=K2f), 0, tolerance=0.001)
  expect_equal(riverdistance(startseg=2, endseg=3, startvert=1, endvert=25, rivers=K2fl), 0, tolerance=0.001)
})

filepath <- system.file("extdata", package="riverdist")
WFT2 <- line2network(path=filepath, layer="West_Fork_Trib2")
test_that("more line2network", {
  expect_equal(length(line2network(sf=abstreams0$sf)$lines), 179)
  expect_equal(length(line2network(sf=Koyukuk0$sf)$lines), 26)
  expect_equal(length(line2network(sf=Koyukuk1$sf)$lines), 17)
  expect_equal(length(line2network(sf=Koyukuk2$sf_current)$lines), 31)
  expect_equal(length(line2network(sf=Kenai1$sf)$lines), 152)
  expect_error(line2network(path=filepath, layer="fakefish_UTM5"), "Invalid input.  Either specified shapefile is not a linear feature, 
         or not all geometry types are LINESTRING or MULTILINESTRING.")
  expect_equal(length(WFT2$lines), 1)
  expect_equal(nrow(WFT2$lines[[1]]), 61)
  expect_equal(riverdistance(startseg=1, endseg=1, startvert=10, endvert=20, rivers=WFT2), 673.6803)
})
