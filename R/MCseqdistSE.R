#' @title
#' Distance standard errors derived from sets of MC-replicated sequences
#'
#' @description
#' Computes the mean and standard deviation of
#' each element of the pairwise distance matrix
#' across sets of MC-replicated sequences.
#'
#' @param MCrseqdata list of MC-replicated sequence datasets of class \code{stslist}. The last element is supposed to be the observed dataset.
#' @param dissrepl list, string, or object of class \code{u.diss}. If a list, list of same length as \code{MCrseqdata}. List of dissimilarity matrices or \code{dist} objects. If a character string, a method name for computing the dissimilarities with \code{\link{MCudist}}. Can also be an object of class \code{u.diss} previously computed with \code{MCudist}.
#' @param udiss logical. When \code{dissrepl} is a distance method, should distance be computed with \code{\link{MCudist}}. See details.
#' @param full.matrix logical. Should dissimilarities be organized in matrix form? Default is \code{FALSE} in which case dissimilarity matrices are converted into \code{dist} objects. If \code{TRUE}, dissimilarity \code{dist} objects are converted into matrices.
#' @param ... additional arguments passed to \code{\link{MCudist}} or \code{\link{MCdisslist}} when \code{dissrepl} is a method name.
#'
#' @details
#' Providing \code{u.diss} distances or computing distances with \code{MCudist} may be faster and can save space when the number of unique replicated sequences is smaller than the sample size times the squared root of R, which can be checked with \code{\link{MCnunique}}. When the number of unique replicated sequences largely exceeds the threshold, it is more efficient to compute distance matrices separately for each updated set of sequences with \code{\link{MCdisslist}} or by setting \code{udiss=FALSE}.
#'
#'
#' @export
#' @importFrom TraMineR seqdist
#'
#' @return Five objects:\cr
#'  \code{MCmean} Mean of distance objects over replicated sets of sequences.\cr
#'  \code{MCsd} Standard deviation of distances over replicated sets of sequences.\cr
#'  In addition, when the observed distances are provided as last element of the \code{dissrepl} list:\cr
#'  \code{MCbias} Difference between observed distance and \code{MCmean}\cr
#'  \code{MCse} Standard error of individual distances.\cr
#'  \code{MCmse} Mean square error of individual distances.\cr
#'  The five objects are of class \code{dist} when \code{attr(MCrseqdata,"toref")==FALSE} and matrices otherwise.
#'
#' @seealso \code{\link{MCseqReplicate}}, \code{\link{MCdisslist}}, \code{\link{MCudist}}, \code{\link{print.distMC}}, \code{\link{summary.distMC}}
#' @examples
#' # example code
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
#' ## 3 MC-replicated sequence datasets
#' altseq.list <- MCseqReplicate(s.exdata, J=1, R=3, include.obs=TRUE)
#' ## list of dissimilarity matrices
#' disslist <- MCdisslist(altseq.list, method="HAM")
#'
#' MCdselist <- MCseqdistSE(disslist)
#' print(MCdselist)
#'
#' MCratioslist <- MCratios(MCdselist)
#' print(MCratioslist)
#'

MCseqdistSE <- function(dissrepl="LCS", MCrseqdata=NULL, udiss = FALSE, full.matrix=FALSE, ...){

  if (is.null(MCrseqdata) & !inherits(dissrepl,"udist") & !is.list(dissrepl)){
    stop("MCrseqdata cannot be NULL when dissrepl is a distance method")
  }
  if (!is.list(MCrseqdata) & !is.null(MCrseqdata)){
    stop("MCrseqdata must be NULL or a list of stslist objects")
  }

  toref <- attr(dissrepl,"toref")
  if (is.null(toref)) toref <- FALSE

  if (is.character(dissrepl)){
    if (udiss){
      dissrepl <- suppressMessages(MCudist(MCrseqdata, method=dissrepl, ...))
    }
    else {
      #dissrepl <- suppressMessages(lapply(MCrseqdata, seqdist,
      #                                    method=dissrepl, full.matrix=full.matrix, ...))
      dissrepl <- MCdisslist(MCrseqdata, method=dissrepl, full.matrix=full.matrix, ...)
    }
  }
  else if (!is.list(dissrepl) & !inherits(dissrepl,"u.diss"))
    stop("Bad dissrepl type!")
  else if (!is.null(MCrseqdata)) message("MCrseqdata ignored because distances are provided!")

  attrs <- attributes(dissrepl)
  if (is.list(dissrepl) & !toref) {
    if (full.matrix)
      dissrepl <- lapply(dissrepl, as.matrix)
    else if (!full.matrix)
      dissrepl <- lapply(dissrepl, as.dist)
    attributes(dissrepl) <- attrs
  }

  ## Removing obs if present.

  diss.o <- NULL

  #print(attr(dissrepl,"obs"))

  if (is.null(attr(dissrepl,"obs"))) attr(dissrepl,"obs") <- FALSE
  if (inherits(dissrepl,"u.diss")){
    sdx <- attr(dissrepl,"sdx")
    N <- attr(dissrepl,"N")
    if (attr(dissrepl,"obs")){
      diss.o <- MCExtractDist(dissrepl, k=N)
      N <- N-1
    }
    sumDiss <- MCExtractDist(dissrepl, k=1)
    idnames <- sub("R1-","",attr(sumDiss,"Labels"))
    sumDiss2 <- sumDiss^2
    for (k in 2:N){
      sumDiss <- sumDiss + MCExtractDist(dissrepl, k=k)
      sumDiss2 <- sumDiss2 + MCExtractDist(dissrepl, k=k)^2
    }
    MCmean <- sumDiss/N
    MCsd <- sumDiss2/N - MCmean^2
  } else if (is.list(dissrepl)){
    if (attr(dissrepl,"obs")){
      diss.o <- dissrepl[[length(dissrepl)]]
      dissrepl <- dissrepl[-length(dissrepl)]
    }
    if (toref){
      idnames <- sub("R1-","",rownames(dissrepl[[1]]))
      #refnames <- paste0("ref.",colnames(dissrepl[[1]]))
    } else {
      idnames <- sub("R1-","",attr(dissrepl[[1]],"Labels"))
    }
    N <- length(dissrepl)
    MCmean <- Reduce("+", dissrepl)/N
    dissSqlist <- lapply(dissrepl, function(x){x^2})
    meanSquaredDiss <- Reduce("+", dissSqlist)/N
    MCsd <- meanSquaredDiss - MCmean^2
  }

  MCsd[MCsd<0] <- 0 ## negative values can result from numerical precision
  MCsd <- sqrt(MCsd*N/(N-1))

  attributes(dissrepl) <- attrs
  attr(dissrepl,"toref") <- toref

  if (toref){
    rownames(MCmean) <- rownames(MCsd) <-idnames
    #colnames(MCmean) <- colnames(MCsd) <- refnames
  }
  else {
    MCmean <- as.dist(MCmean)
    MCsd <- as.dist(MCsd)
    attr(MCmean,"Labels") <- attr(MCsd,"Labels") <- idnames
  }
  ret <- list(MC.mean=MCmean, MC.sd=MCsd, N=N)

  if(!is.null(diss.o)){ ## compute bias
    MCbias <- diss.o - MCmean
    MCmse <- MCsd^2 + MCbias^2
    MCse <- sqrt(MCmse)
    if (toref) {
      rownames(MCbias) <- rownames(MCmse) <- rownames(diss.o) <- idnames
      #colnames(MCbias) <- colnames(MCmse) <- colnames(diss.o) <- refnames
    }
    else {
      MCbias <- as.dist(MCbias)
      MCmse <- as.dist(MCmse)
      diss.o <- as.dist(diss.o)
      attr(MCbias,"Labels") <- attr(MCmse,"Labels") <- attr(diss.o,"Labels") <- idnames
    }
    ret <- c(ret, list(diss.o = diss.o, MC.bias=MCbias, MC.se = MCse, MC.mse=MCmse))
  }

  attr(ret,"obs") <- attr(dissrepl,"obs")
  attr(ret,"toref") <- toref
  class(ret) <- c(class(ret),list("distMC"))
  return(ret)

}



