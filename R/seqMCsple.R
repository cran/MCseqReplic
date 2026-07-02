## seqMCset joint sets of sequences replicated from the pair of sequences
## defined by the rows of dss and sdur
##  dss and sdur two-row matrices
seqMCset <- function(dss,sdur,Jprob,R=10,ch.meth,
                     jfixed, kchanges){
  mdss <- as.matrix(dss)
  ncs <- ncol(mdss)
  set1 <- set2 <- matrix(nrow=R,ncol=ncs)
  #dur1NA <- is.na(sdur[1,])
  #dur2NA <- is.na(sdur[2,])
  for (i in 1:R){
    if (ch.meth==1)
      dur <- t(apply(sdur,1,ch.dur,Jprob=Jprob, jfixed = jfixed, k=kchanges))
    else if (ch.meth==2)
      dur <- t(apply(sdur,1,ch.dur.indep,Jprob=Jprob, jfixed = jfixed, k=kchanges))
    else if (ch.meth==3)
      dur <- t(apply(sdur,1,ch.dur.relat,Jprob=Jprob, jfixed = jfixed, k=kchanges))

    if (nrow(dur)==1) dur <- t(dur)
    wdur1 <- which(dur[1,]>0)
    wdur2 <- which(dur[2,]>0)
    set1[i,wdur1] <- paste0(mdss[1,wdur1],"/",dur[1,wdur1])
    set2[i,wdur2] <- paste0(mdss[2,wdur2],"/",dur[2,wdur2])
    #set1[i,dur1NA] <- NA
    #set2[i,dur2NA] <- NA
    #set1[i,dur[1,]<=0] <- NA
    #set2[i,dur[2,]<=0] <- NA
  }
  #print(set1[1:3,])
  #print(set2[1:3,])
  sts1 <- seqformat(set1, from="SPS", to="STS", SPS.in=list(xfix="", sdsep="/"), stsep="-")
  sts2 <- seqformat(set2, from="SPS", to="STS", SPS.in=list(xfix="", sdsep="/"), stsep="-")
  rownames(sts1) <- paste0(1:R,"B1")
  rownames(sts2) <- paste0(1:R,"B2")

  #print(length(alphabet(dss)))
  seq2sets <- suppressMessages(
    seqdef(rbind(sts1,sts2), alphabet=alphabet(dss))
  )
  #print(length(alphabet(seq2sets)))
  #print("seqdef done!")
  return(seq2sets)
}

## seqdistsple list of MC distances between a pair of sequences or summary of the distances
## importFrom stats mean sd

seqdistsple <- function(seq2sets, method, what="distMC", ...){

  #indel <- list(...)[["indel"]]
  #sm <- list(...)[["sm"]]
  #print(length(indel))
  #print(dim(sm))
  #print(length(alphabet(seq2sets)))

  MCdi <- suppressMessages(seqdist(seq2sets,
                                   refseq=list(1:(nrow(seq2sets)/2),(nrow(seq2sets)/2+1):nrow(seq2sets)),
                                   method=method,
                                   ...))
  distMC <- list(mean = mean(MCdi), std = sd(MCdi), distMC = MCdi)
  class(distMC) <- c(class(distMC),"distMC")
  if (what=="distMC") return(distMC)
  else return(MCdi)
}
