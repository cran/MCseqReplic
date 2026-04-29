#' @title Correlation between 1st MDS factor of observed and MC-simulated distances
#'
#' @param disslist List of matrices or dist objects: the MC-replicated dissimilarities
#' @param diss.o Matrix or dist object: Observed dissimilarities
#' @param method String. One of \code{"Spearman"} (default) and \code{"Pearson"}.
#' @param weights vector of doubles. Case weights. If \code{NULL} (default), equal weights are used.
#' @param what String. One of \code{"corr"} (correlations, default), \code{"mds"} (list of mds scores), and \code{"both"}.
#' @param core Integer. Number of cores for parallel computing.
#' @param snow Logical. If \code{TRUE}, \code{doSNOW} is used for parallel computing, otherwise \code{doParallel} is used.
#' @param silent Logical. Should waiting and timing messages be hidden?
#'
#' @details
#' When \code{diss.o=NULL}, the last element of \code{disslist} is taken as \code{diss.o} and the other elements as sets of MC-replicated dissimilarities.
#'
#' @return when \code{what="corr"}, vector of correlation between mds of dissimilarities in MC-replicated sets, when \code{what="mds"}, of first mds scores, and when \code{what="both"}, list with \code{corr} as first element and \code{mdslist}, the list of mds scores as second element.
#'
#'
#' @importFrom vegan  wcmdscale
#' @importFrom wCorr  weightedCorr
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makePSOCKcluster stopCluster makeCluster detectCores
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach %dopar% %do% %:%
#' @importFrom utils txtProgressBar setTxtProgressBar
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
#' disslist <- MCdisslist(altseq.list)
#' MCmdscorr(disslist)
#'
#'
MCmdscorr <- function(disslist, diss.o=NULL, method="Spearman", weights=NULL, what="corr",
                      core=1, snow=TRUE, silent=FALSE){
  message(method," correlations between 1st MDS factor of observed and replicated distances")
  N <- length(disslist)
  if (is.null(diss.o)){
    message("Extracting diss.o from disslist")
    diss.o <- disslist[[N]]
    disslist <- disslist[-N]
    N <- N-1
  }

  ## setting weights
  if (inherits(diss.o,'dist')){
    ncases <- attr(diss.o,'Size')
  } else {
    ncases <- nrow(diss.o)
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

  message("Computing first MDS scores of observed dissimilarities")
  mds.o <- vegan::wcmdscale(diss.o, k=1, w=weights)
  message("Computing first MDS scores of replicated dissimilarities")
  if (!silent & (core==1 | snow)) {
    #cat("\n")
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    pb <- txtProgressBar(min=1, max=N, initial=1, style=3)
  }
  else opts <- list()

  if (core==1){
    #mdslist <- lapply(disslist,wcmdscale,k=1, w=weights)
    mdslist <- list()
    for (i in 1:N){
      if (!silent) setTxtProgressBar(pb, i)
      mdslist[[i]] <- vegan::wcmdscale(disslist[[i]],k=1, w=weights)
    }
    if (!silent) close(pb)
  }
  else { ## parallel computing
    #pb <- txtProgressBar(min=1, max=k-1, initial=1, style=3)
    #cl <- makePSOCKcluster(core)
    cl <- makeCluster(core, type="SOCK")
    if (snow)
      registerDoSNOW(cl)
    else
      registerDoParallel(cl)

    mdslist <- foreach (i = 1:N, .combine=c, .options.snow=opts) %dopar% {
      list(vegan::wcmdscale(disslist[[i]],k=1, w=weights))
    }

    stopCluster(cl)
    if(!silent & snow) {
      close(pb)
    }
  }


  if (what != "mds"){
    message("Computing correlations")
    corr <- sapply(mdslist, function(x){wCorr::weightedCorr(x=x, y=mds.o, method=method, weights=weights)})
    for (i in which(corr < 0)) {
      corr[i] <- -corr[i]
      mdslist[[i]] <- -mdslist[[i]]
    }
  }
  mdslist <- c(mdslist, list(mds.o=mds.o))

  if (what == "mds")
    return(mdslist)
  else if (what %in% c("cor","corr"))
    return(corr)
  else if (what == "both"){
    return(list(corr=corr, mdslist = mdslist))
  }


}


