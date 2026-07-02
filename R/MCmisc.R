## build table of i and j corresponding to position in dist vector

getijtable <- function(k){
  x <- NULL
  y <- NULL
  for (i in 1:(k-1)){
    x <- c(x,rep(i,k-i))
    y <- c(y,(i+1):k)
  }
  cbind(x,y)
}


## position in dist object given i and j
didx <- function(i,j,k){
  k*(i-1) - i*(i-1)/2 + j-i
}

MCdistrib <- function(dss,dur,id2=1:2, Jprob, R=100, model="keep.dss",
                      jfixed=FALSE, kchanges=NULL,
                      method="LCS", ...){
  Jprob <- Jprob/sum(Jprob)
  if (model=="keep.dss")
    ch.meth <- 1
  else if (model=="indep")
    ch.meth <- 2
  else if (model=="relative")
    ch.meth <- 3
  else
    stop("Unkown model value")

  #seq2sets <- seqdistMCSE:::seqMCset(dss[id2,],dur[id2,], Jprob=Jprob, R=R)
  #MC.dist <- seqdistMCSE:::seqdistsple(seq2sets, method=method, what="MC distances", ...)
  seq2sets <- seqMCset(dss[id2,],dur[id2,], Jprob=Jprob, R=R, ch.meth=ch.meth,
                       jfixed=jfixed, kchanges=kchanges)
  MC.dist <- seqdistsple(seq2sets, method=method, what="MC distances", ...)
  return(MC.dist)
}

# wrank adapted from wCorr package
w.rank <- function(x, w=rep(1,length(x))) {
  # sort by x so we can just traverse once
  ord <- order(x)
  ##rord <- (1:length(x))[order(ord)] # reverse order
  rord <- order(ord) # reverse order
  xp <- x[ord] # x, permuted
  wp <- w[ord] # weights, permuted
  rnk <- rep(NA, length(x)) # blank ranks vector
  # setup first itteration
  t1 <- 0 # total weight of lower ranked elements
  i <- 1 # index
  t2 <- 0 # total weight of tied elements (including self)
  n <- 0 # number of tied elements
  while(i < length(x)) {
    t2 <- t2 + wp[i] # tied weight increases by this unit
    n <- n + 1
    if(xp[i+1] != xp[i]) { # the next one is not a tie
      # find the rank of all tied elements
      rnki <- t1 + (1 + (t2-1)/2)
      # push that rank to all tied units
      for(ii in 1:n) {
        rnk[i-ii+1] <- rnki
      }
      # reset for next iteration
      t1 <- t1 + t2 # new total weight for lower values
      t2 <- 0 # new tied weight starts at 0
      n <- 0
    }
    i <- i + 1
  }
  # final row
  t2 <- t2 + wp[i] # add final weight to tied weight
  rnki <- t1 + (1 + (t2-1)/2) # final rank
  # push that rank to all final tied units
  for(ii in 1:n) {
    rnk[i-ii+1] <- rnki
  }
  # order by incoming index, so put in the original order
  rnk[rord]
}

###############

set_ncores <- function(core){
  if (is.null(core) || (is.character(core) & core != "auto"))
    stop('core must be numeric or "auto"')
  if (core=="auto" || core > 2){
    available.cores <- parallelly::availableCores(logical = TRUE)
    if (core == "auto") {
      core <- available.cores - 1
    } else {
      core <- round(min(core, available.cores-1))
    }
  }
  core <- max(core,1)
}
