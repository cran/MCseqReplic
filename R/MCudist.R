#' @title Dissimilarities between unique replicated sequences
#'
#' @description
#' Returns the dissimilarity matrix (or \code{dist} object) between merged replicated sequences with the disaggregation indexes as attribute.
#'
#' @param MCrseqdata list of replicated \code{stslist} state sequence datasets (all of same size and with same alphabet)
#' @param method string. Name of distance method (see \code{seqdist}).
#' @param seqref state sequence object of class \code{stslist}. Fixed reference sequences.
#' @param ... Further arguments passed to \code{seqdist}
#'
#' @export
#' @importFrom WeightedCluster wcAggregateCases
#'
#' @returns object of class \code{u.diss} (pairwise dissimilarities between unique sequences) with two attributes: \code{sdx}, inverted aggregation indexes, \code{N}, number of datasets, and \code{obs}, logical indicating whether \code{k=N} corresponds to observed sequences.
#' @seealso \code{\link{MCExtractDist}}
#' @export
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
#' (altseq.list <- MCseqReplicate(s.exdata, J=1, R=3))
#'
#' MCnunique(altseq.list, check=TRUE)
#'
#'
#' u.diss <- MCudist(altseq.list, method="LCS", full.matrix=FALSE)
#' ## Dissimilarities within first MC-set
#' MCExtractDist(u.diss, 1)
#'
#' ## list of dissimilarity matrices
#' disslist <- MCdisslist(altseq.list, use.udiss=TRUE)
#'
MCudist <- function(MCrseqdata, method="LCS", seqref=NULL, ...){

  N <- length(MCrseqdata)
  # collecting list of replicated sets into a single object
  MCrseqdata <- do.call(rbind, MCrseqdata)

  aggCases <- wcAggregateCases(MCrseqdata, weights=attr(MCrseqdata,"weights"))
  ## u.seqdata is set of unique sequences
  u.seqdata <- MCrseqdata[aggCases[["aggIndex"]],]
  sdx <- aggCases[["disaggIndex"]]
  names(sdx) <- rownames(MCrseqdata)

  ## matrices of dissimilarities between unique sequences
  if (is.null(seqref)){
    u.diss <- #suppressMessages(
      seqdist(u.seqdata, method=method, ...)
    #)
    attr(u.diss,"toref") <- FALSE
  } else { ## dist to ref
    u.diss <- seqdistToRef(u.seqdata, ref=seqref, method=method, ... )
    attr(u.diss,"toref") <- TRUE
  }

  attr(u.diss,"sdx") <- sdx
  attr(u.diss,"N") <- N
  attr(u.diss,"obs") <- attr(MCrseqdata,"obs")
  class(u.diss) <- c(class(u.diss),"u.diss")

  return(u.diss)
}

#' @title Extract k-th dissimilarity matrix from u.diss
#'
#' @param u.diss \code{u.diss} object returned by \code{MCudist}: dissimilarities between unique replicated sequences.
#' @param k integer. Subset index number for which the dissimilarity matrix must be extracted
#' @param full.matrix logical. If \code{FALSE}, the distance matrix is returned as a \code{dist} object. Ignored for distances to reference sequences.
#' @returns a dissimilarity matrix or distance object.
#' @seealso \code{\link{MCudist}}
#'
#' @export

MCExtractDist <- function(u.diss, k, full.matrix=FALSE){
  toref <- attr(u.diss,"toref")
  if (is.null(toref)) toref <- FALSE
  # length of seqdata
  sdx <- attr(u.diss,"sdx")
  N <- attr(u.diss,"N")
  n <- length(sdx)/N
  if (floor(n)*N != length(sdx))
    stop("length(sdx) must be a multiple of N")
  idk <-  (1:n)+((k-1)*n)
  if (toref){
    diss <- u.diss[sdx[idk],,drop=FALSE]
    rownames(diss) <- names(sdx[idk])
  } else {
    u.diss <- as.matrix(u.diss)
    diss <- u.diss[sdx[idk],sdx[idk],drop=FALSE]
    dimnames(diss) <- list(names(sdx[idk]),names(sdx[idk]))
  }
  if (!full.matrix & !toref){
    diss <- as.dist(diss)
  }
  return(diss)
}

#' @title Number of unique replicated sequences
#'
#' @param MCrseqdata list of replicated \code{stslist} state sequence datasets (all of same size and with same alphabet.
#' @param check logical. When \code{TRUE}, check if the number of unique replicated sequences is less than  \code{sqrt(number of replicated datasets) * (sample size)}?
#'
#'
#' @return \code{nu} number of unique replicated sequences and, when \code{check=TRUE}, \code{u.ok} the check result.
#' @importFrom WeightedCluster wcAggregateCases
#' @export
#' @seealso \code{\link{MCudist}}, \code{\link{MCseqdistSE}}
#'
MCnunique <- function(MCrseqdata, check=FALSE){
  N <- length(MCrseqdata)
  n <- nrow(MCrseqdata[[1]])

  # collecting list of replicated sets into a single object
  MCrseqdata <- do.call(rbind, MCrseqdata)
  aggCases <- wcAggregateCases(MCrseqdata, weights=attr(MCrseqdata,"weights"))
  ## u.seqdata is set of unique sequences
  u.seqdata <- MCrseqdata[aggCases[["aggIndex"]],]
  nu <- nrow(u.seqdata)
  if (check){
    u.ok <- (nu < n*sqrt(N))
    return(list(nu=nu, u.ok=u.ok))
  } else return(nu)
}
