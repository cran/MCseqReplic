#' @title Correlation between observed and MC-simulated distances
#'
#' @param disslist List of matrices or dist objects: the MC-replicated dissimilarities
#' @param diss.o Matrix or dist object: Observed dissimilarities
#' @param method String. One of \code{"Spearman"} (default) and \code{"Pearson"}.
#' @param weights vector of doubles. Case weights. If \code{NULL} (default), equal weights are used.
## #' @param what String. One of \code{"corr"} (correlations, default).
## #' @param core Integer. Number of cores for parallel computing.
## #' @param snow Logical. If \code{TRUE}, \code{doSNOW} is used for parallel computing, otherwise \code{doParallel} is used.
## #' @param silent Logical. Should waiting and timing messages be hidden?
#'
#' @details
#' When \code{diss.o=NULL}, the last element of \code{disslist} is taken as \code{diss.o} and the other elements as sets of MC-replicated dissimilarities.
#'
#' @return vector of correlation between observed and MC-dissimilarities.
#'
#'
#' @importFrom wCorr  weightedCorr
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats as.dist
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
#' altseq.list <- MCseqReplicate(s.exdata, J=1, R=3, include.obs=TRUE)
#' ## list of dissimilarity matrices
#' disslist <- MCdisslist(altseq.list)
#' MCdisscorr(disslist)
#'
#'
MCdisscorr <- function(disslist, diss.o=NULL, method="Spearman", weights=NULL){

  toref <- attr(disslist,"toref")
  if (is.null(toref)) toref <- FALSE

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

  if (toref){
    nc <- ncol(diss.o)
    corr <- list()
    for (i in 1:nc){
      corr[[i]] <- sapply(disslist, function(x){wCorr::weightedCorr(x=x[,i], y=diss.o[,i], method=method, weights=weights)})
    }
    names(corr) <- colnames(diss.o)

  } else {
    w <- as.dist(tcrossprod(weights))
    corr <- sapply(disslist, function(x){wCorr::weightedCorr(x=as.dist(x), y=as.dist(diss.o), method=method, weights=w)})
  }

  return(corr)
}


