## timing change methods

## definition of move depends on Jprob
move <- function(n, Jprob){
  #print(n)
  #print(Jprob)
  if (length(Jprob) == 1){
    #JJ <- Jprob * 2 + 1
    return(floor(runif(n) * (Jprob * 2 + 1) - Jprob))
  } else {
    Jprob <- Jprob/sum(Jprob) ## already normalized
    JJ <- length(Jprob)
    Jdist <- cumsum(Jprob)
    ## drop intervals with zero prob
    nonzero <- which(Jprob != 0)
    Jprob <- Jprob[nonzero]
    Jdist <- Jdist[nonzero]
    JL <- length(Jdist)
    p <- runif(n)
    #print(p)
    res <- sapply(p, function(x) {which(Jdist >= x)[1]})
    #print(res)
    #print(nonzero)
    res[res<=0] <- 1
    res[res > JL] <- JL ## should never happen
    res <- nonzero[res] - (JJ + 1)/2
    #print(res)
    return(res)
  }
}

setklist <- function(k, nt){
  if (is.null(k))
    k <- nt
  else if (k <0) {
    k <- floor(runif(1,min=0,max=nt)+.5)
  }
  if (k < nt)
    klist <- sample(nt,k)
  else
    klist <- 1:nt
  return(klist)
}

## method="keep.dss"
ch.dur <- function(sduri, Jprob, jfixed=TRUE, k=NULL){
  #J <- length(Jprob)
  #JJ <- 2*J+1
  nt <- sum(!is.na(sduri))-1

  ##change <- floor(runif(nt) * JJ - J)
  if (nt > 0) {
    klist <- setklist(k, nt)
    change <- integer(nt)
    if (jfixed)
      change[klist] <- move(1, Jprob=Jprob)
    else
      change[klist] <- move(length(klist), Jprob=Jprob)
    #print(change)
    for (i in klist) {
      if (change[i] < 0)
        change[i] <- -min(-change[i], sduri[i] - 1)
      else if (change[i] >0)
        change[i] <- min(change[i], sduri[i+1] - 1)

      sduri[i] <- sduri[i] + change[i]
      sduri[i+1] <- sduri[i+1] - change[i]
    }
  }
  return(sduri)
}


## method="indep"
ch.dur.indep <- function(sduri, Jprob, jfixed=FALSE, k=NULL){
  #J <- length(Jprob)
  #JJ <- 2*J+1
  nt <- sum(!is.na(sduri))-1

  ##change <- floor(runif(nt) * JJ - J)
  if (nt > 0) {
    klist <- setklist(k, nt)
    change <- integer(nt)
    if (jfixed)
      change[klist] <- move(1, Jprob=Jprob)
    else
      change[klist] <- move(length(klist), Jprob=Jprob)
    #print(change)
    for (i in klist) {
      if (change[i] < 0){
        chg <- -change[i]
        sduri[i+1] <- sduri[i+1] + chg
        if (sduri[i] >= chg){
          sduri[i] <- sduri[i] - chg
        } else {
          rchg <- chg - sduri[i]
          sduri[i] <-0
          ii <- i-1
          while (ii > 0 & rchg > 0){
            dur <- sduri[ii]
            sduri[ii] <- max(0, dur - rchg)
            rchg <- rchg - (dur - sduri[ii])
            ii <- ii-1
          }
          if (rchg > 0) { ## must adjust change at i+1
            sduri[i+1] <- sduri[i+1] - rchg
          }
        }
      }
      else {
        chg <- change[i]
        sduri[i] <- sduri[i] + chg
        if (sduri[i+1] >= chg){
          sduri[i+1] <- sduri[i+1] - chg
        } else {
          rchg <- chg - sduri[i+1]
          sduri[i+1] <-0
          ii <- i+2
          while (ii < (nt + 1) & rchg > 0){
            dur <- sduri[ii]
            sduri[ii] <- max(0, dur - rchg)
            rchg <- rchg - (dur - sduri[ii])
            ii <- ii+1
          }
          if (rchg > 0) { ## must adjust change at i
            sduri[i] <- sduri[i] - rchg
          }
        }
      }
    }
  }
  return(sduri)
}

## method="relative"
ch.dur.relat <- function(sduri, Jprob, jfixed=FALSE, k=NULL){
  #J <- length(Jprob)
  #JJ <- 2*J+1
  nt <- sum(!is.na(sduri))-1

  ##change <- floor(runif(nt) * JJ - J)
  if (nt > 0) {
    klist <- setklist(k, nt)
    change <- integer(nt)
    if (jfixed)
      change[klist] <- move(1, Jprob=Jprob)
    else
      change[klist] <- move(length(klist), Jprob=Jprob)
    #print(change)
    for (i in klist) {
      if (change[i] < 0){
        chg <- -change[i]
        sduri[nt+1] <- sduri[nt+1] + chg
        if (sduri[i] >= chg){
          sduri[i] <- sduri[i] - chg
        } else {
          rchg <- chg - sduri[i]
          sduri[i] <-0
          ii <- i-1
          while (ii > 0 & rchg > 0){
            dur <- sduri[ii]
            sduri[ii] <- max(0, dur - rchg)
            rchg <- rchg - (dur - sduri[ii])
            ii <- ii-1
          }
          if (rchg > 0) { ## must adjust change at nt+1 (last spell)
            sduri[nt+1] <- sduri[nt+1] - rchg
          }
        }
      }
      else {
        chg <- change[i]
        sduri[i] <- sduri[i] + chg
        if (sduri[nt+1] >= chg){
          sduri[nt+1] <- sduri[nt+1] - chg
        } else {
          rchg <- chg - sduri[nt+1]
          sduri[nt+1] <-0
          ii <- nt
          while (ii > i & rchg > 0){
            dur <- sduri[ii]
            sduri[ii] <- max(0, dur - rchg)
            rchg <- rchg - (dur - sduri[ii])
            ii <- ii-1
          }
          if (rchg > 0) { ## must adjust change at i
            sduri[i] <- sduri[i] - rchg
          }
        }
      }
    }
  }
  return(sduri)
}

