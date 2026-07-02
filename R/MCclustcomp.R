#' @title Comparing MC-clusters with cluster of observed data
#'
#' @description
#' Cluster comparison indices (CCI) between clusters based on observed data and each of the MC-replicated partitions.
#'
#'
#' @param clustlist List of MC-replicated vectors of cluster memberships.
#' @param clust.o Cluster memberships based on observed dissimilarities.
#' @param weights vector of doubles. Case weights. If \code{NULL} (default), equal weights are used.
#' @param AMI logical. AMI is more costly! Should AMI also be computed. Deafult is FALSE.
## #' @param what String. One of the indexes computed by \code{aricode} package such as {\code{"ARI"} (adjusted Rand index, default), \code{"RI"} (Rand index), \code{"VI"} (variation of information), \code{"NVI"}. Can also be "all \code{"all"}.
#'
#' @details
#' When \code{diss.o=NULL}, the last element of \code{disslist} is taken as \code{diss.o} and the other elements as sets of MC-replicated dissimilarities.
#'
#' For a description of the CCIs, see \insertCite{SundqvistChiquetRigalli2022CS}{MCseqReplic}.
#'
#' @return A table with in columns the list of comparison scores provided by \code{aricode::\link[aricode]{compare_clustering}} for each replicated set, except Chi2, which is replaced by Cramer's V.
#'
#' @seealso \code{\link[aricode]{compare_clustering}}
#'
#' @references
#'   \insertNoCite{ChiquetRigaillSundqvistDervieuxEtAl2023aricodeEfficientComputations}{MCseqReplic}
#'   \insertNoCite{SundqvistChiquetRigalli2022CS}{MCseqReplic}
#'   \insertAllCited{}
#'
## #' @importFrom mclust adjustedRandIndex
#' @importFrom aricode compare_clustering
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
#' library(WeightedCluster)
#' clust.o <- wcKMedoids(diss.o, k=2, cluster.only=TRUE)
#' clustlist <- lapply(disslist, wcKMedoids, k=2, cluster.only=TRUE)
#' res <- MCclustcomp(clustlist, clust.o=clust.o)
#' res
#'
MCclustcomp <- function(clustlist, clust.o=NULL, weights=NULL, AMI=FALSE){
  message(" Cluster comparison scores between MC-Clusters and clusters of observed data")

  N <- length(clustlist)
  if (is.null(clust.o)){
    message("Extracting clust.o from clustlist")
    clust.o <- clustlist[[N]]
    clustlist <- clustlist[-N]
    N <- N-1
  }
  n <- length(clust.o)
  k <- length(unique(clust.o))

  if (!is.null(weights))
    message("Weights currently not supported! Weights are ignored.")
    ## setting weights
  # ncases <- length(clust.o)
  # if (is.null(weights)) {
  #   weights <- rep(1, ncases)
  #   weighted=FALSE
  # } else {
  #   if (length(weights) != ncases)
  #     stop("length of weights not equal to number of cases!")
  #   if (all(weights==1))
  #     weighted <- FALSE
  #   else weighted <- TRUE
  # }

  ## should write own function for weighted data
  ## aricode::compare_clustering does not support weights

  ret <- lapply(clustlist, function(x) {as.matrix(
    unlist(aricode::compare_clustering(x,clust.o,AMI=AMI)))})
  ret <- do.call(cbind,ret)
  ret["Chi2",] <- sqrt(ret["Chi2",]/(n*(k-1)))
  rownames(ret)[rownames(ret)=="Chi2"] <- "V"
  ret
}

