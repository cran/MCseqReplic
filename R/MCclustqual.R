#' @title Cluster quality measures by MC-sets
#'
#' @description
#' Cluster quality measures for a range of number of groups by MC-replicated set.
#'
#'
#' @param disslist List of MC-dissimilarity matrices (or \code{dist} objects).
## #' @param diss.o
#' @param ncluster integer vector. Range of number of groups. Default is \code{2:10}.
#' @param clustmeth character. Clustering method. Either \code{"PAM"} (default) or a \code{stats::\link{hclust}} method.
#' @param weights vector of doubles. Case weights. If \code{NULL} (default), equal weights are used.
#' @param core Integer. Number of cores for parallel computing.
#' @param snow Logical. If \code{TRUE}, \code{doSNOW} is used for parallel computing, otherwise \code{doParallel} is used.
#' @param silent Logical. Should waiting and timing messages be hidden?
#' @param ... Further arguments passed to clustering functions.
#'
## #' @param what String. One of the indexes computed by \code{aricode} package such as {\code{"ARI"} (adjusted Rand index, default), \code{"RI"} (Rand index), \code{"VI"} (variation of information), \code{"NVI"}. Can also be "all \code{"all"}.
#'
#' @details
#' When \code{attr(MCdisslist,"obs")} is \code{TRUE}, the last element of \code{disslist} is treated as the dissimilarity matrix of the observed sequences.
#'
#' @return A list with two lists: \code{qual.tab}, list of tables of cluster quality statistics per MC-dissimilarity matrix, and \code{qual.max} list of cluster number $k$ for which the statistics reach their maximum (minimum for HC), \code{max.freq}, the frequency table of maximum over the MC-replicated sets, and \code{qual.obs}, cluster quality indexes for the observed sequences.
#'
#' @seealso \code{\link[WeightedCluster]{as.clustrange}},  \code{\link[WeightedCluster]{wcKMedRange}}
#'
#'
#' @importFrom WeightedCluster wcKMedRange as.clustrange
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makePSOCKcluster stopCluster makeCluster detectCores
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach %dopar% %do% %:%
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats hclust
#'
#' @export
#'
#' @examples
#' ## mini test data, 6 sequences of length 4, 4 unique sequences
#' exdata <- read.table(text="
#'                 a a b b
#'                 a a b b
#'                 b b a a
#'                 a c c b
#'                 b b a c
#'                 b b a c
#'                 ")
#' weights=rep(1, nrow(exdata))
#' s.exdata <- seqdef(exdata, weights = weights, id=paste("id",1:nrow(exdata), sep=""))
#'
#' ## 3 altered sequence datasets
#' set.seed(25)
#' altseq.list <- MCseqReplicate(s.exdata, J=1, R=3)
#' ## list of dissimilarity matrices
#' disslist <- MCdisslist(altseq.list, method="LCS")
#' diss.o <- seqdist(s.exdata, method="LCS")
#' ## cluster per MC-dissimilarity matrices
#' res <- MCclustqual(disslist,ncluster=3)
#' res
#'
#'
MCclustqual <- function(disslist, ncluster=10, clustmeth="PAM",  weights=NULL,
                        core=1, snow=TRUE, silent=FALSE, ...){

  N <- length(disslist)
  isobs <- attr(disslist,"obs")
  if (is.null(isobs)) isobs <- FALSE

  if (inherits(disslist[[1]],"dist")){
    ncases <- attr(disslist[[1]],"Size")
  } else {
    ncases <- nrow(disslist[[1]])
  }
  if (is.null(weights)) {
    weights <- rep(1, ncases)
    weighted=FALSE
  } else {
    if (length(weights) != ncases)
      stop("length of weights not equal to number of cases!")
    if (all(weights==1))
      weighted <- FALSE
    else weighted <- TRUE
  }

  if (!silent & (core==1 | snow)) {
    #cat("\n")
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    pb <- txtProgressBar(min=1, max=N, initial=1, style=3)
  }
  else opts <- list()

  KMedRangeStats <- function(diss, kvals, weights, ...){
    wcKMedRange(diss, kvals=kvals, weights=weights, ...)$stats
  }
  hclustRange <- function(diss, clustmeth, ncluster, weights, ... ){
    hcl <- hclust(as.dist(diss), method=clustmeth, members=weights)
    as.clustrange(hcl, diss=diss, ncluster=ncluster, weights=weights, ...)$stats
  }
  if (core==1) {
    if (clustmeth %in% c("PAM","KM")){
      qlist <- lapply(disslist, KMedRangeStats, kvals=2:ncluster, weights=weights, ...)
    } else {
      qlist <- lapply(disslist, hclustRange, clustmeth=clustmeth, ncluster=ncluster, weights=weights, ...)
    }
    if (!silent) close(pb)
  } else { ## parallel computing
    cl <- makeCluster(core, type="SOCK")
    if (snow)
      registerDoSNOW(cl)
    else
      registerDoParallel(cl)

    qlist <- foreach (i = 1:N, .combine=c, .packages=c("WeightedCluster"), .options.snow=opts) %dopar% {
      if (clustmeth %in% c("PAM","KM")){
        list(KMedRangeStats(disslist[[i]], kvals=2:ncluster, weights=weights, ...))
      } else {
        list(hclustRange(disslist[[i]], clustmeth=clustmeth, ncluster=ncluster, weights=weights, ...))
      }

    }
    stopCluster(cl)
    if(!silent & snow) {
      close(pb)
    }
  }

  if (isobs){
    qobs <- qlist[[N]]
    qlist <- qlist[-N]
    N <- N-1
  } else {
    qobs <- NULL
  }


  which.max.vect <- function(m) {
    stats <- m
    stats[,"HC"] <- -stats[,"HC"] ## For HC we seek the min
    apply(stats,2,which.max) + 1
    }
  max.list <- lapply(qlist, which.max.vect)

  ## distribution of max
  tabmax <- do.call(rbind,max.list)
  tabmaxf <- data.frame(factor(tabmax[,1],levels=2:ncluster))
  #print(tabmaxf)
  maxdist <- data.frame(table(tabmaxf[,1]))
  #print(maxdist)
  rownamesd <- maxdist[,1]
  maxdist <- maxdist[,2]
  for (i in 2:ncol(tabmax)){
    tabmaxf <- cbind(tabmaxf,factor(tabmax[,i],levels=2:ncluster))
    maxdist <- cbind(maxdist,data.frame(table(tabmaxf[,i]))[,2])
  }
  rownames(maxdist) <- rownamesd
  colnames(tabmaxf) <- colnames(maxdist) <- colnames(tabmax)

  ret <- list(qual.tab = qlist, qual.max = tabmaxf, max.freq=maxdist, qual.obs=qobs)

  class(ret) <- c(class(ret),"MCclustQ")
  attr(ret,"obs") <- attr(disslist,"obs")
  return(ret)
}

