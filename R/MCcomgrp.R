#' @title Group comparison by MC-sets
#'
#' @description
#' Collects statistics and p-values for the comparison of groups of sequences.
#'
#' @details
#' The function collects the values of R2 and its p-value \insertCite{StuderRitschardGabadinhoMuller2011SMR}{MCseqReplic} returned by \code{TraMineR::\link[TraMineR]{dissassoc}} and the values of LRT, its p-value, and delta BIC \insertCite{LiaoFasang2020SM}{MCseqReplic} returned by \code{TraMineRextras::\link[TraMineRextras]{dissCompare}}. Since \code{dissCompare} works only with two groups, only R2 and its p-value are returned when there are more than two groups.
#'
#' Except for \code{group} and \code{weights}, \code{dissassoc} and \code{dissCompare} are called by default with the default values of their arguments. This can be changed by passing the wanted arguments as a list to \code{dissassoc.args} and \code{dissCompare.args}.
#'
#' The R2 and its p-value are computed by \code{dissassoc}, which computes th p-value using permutation tests. The default number of permutation is \code{R=1000} but this can be changed by means of the \code{dissassoc.args} argument, for example, by passing \code{dissassoc.args = list(R=500)}.
#'
#' The LRT and delta BIC are computed by \code{dissCompare}, which computes the LRT for samples of \code{s} data, with \code{s} possibly greater than the number of observed data. When \code{s=0} (default in \code{MCcompgrp}), no sampling is applied. \code{dissCompare} computes the p-value of LRT using the appropriate Chi-square distribution. In case of multiple samples, i.e. when \code{s} is smaller than the greatest group size, \code{BFopt=1} is used by default. \code{BFopt=NULL} could generate unpredictable results in that case.
#'
#'
#'
#' @param disslist list of dissimilarity matrices or \code{dist} objects.
#' @param group vector of group memberships of length equal to number of rows of the dissimilarity matrices.
#' @param weights vector of case weights
#' @param dissassoc.args list of additional arguments passed to \code{TraMineR::dissassoc}.
#' @param dissCompare.args list of additional arguments passed to \code{TraMineRextras::dissCompare}.
#' @param verbose logical. Should messages be printed?
#'
#' @importFrom TraMineRextras dissCompare
#'
#' @export
#'
#' @references
#'   \insertNoCite{RitschardLiao2026IJoSRM}{MCseqReplic}
#'   \insertAllCited{}
#'
#' @examples
#' ## mini test data, 6 sequences of length 4, 4 unique sequences
#' exdata <- read.table(text="t1 t2 t3 t4 sex
#'                 a a b b f
#'                 a a b b f
#'                 b b a a f
#'                 a c c b m
#'                 b b a c m
#'                 b b a c m
#'                 ", header=TRUE)
#' weights=rep(1, nrow(exdata))
#' s.exdata <- seqdef(exdata[,1:4], weights = weights, id=paste("id",1:nrow(exdata), sep=""))
#'
#' ## 3 altered sequence datasets
#' set.seed(25)
#' altseq.list <- MCseqReplicate(s.exdata, J=1, R=3)
#' ## list of dissimilarity matrices
#' disslist <- MCdisslist(altseq.list, method="LCS")
#' ## Group comparison per MC-dissimilarity matrices
#' res <- MCcompgrp(disslist,group=exdata$sex)
#' res

MCcompgrp <- function(disslist, group, weights=NULL,
                      dissassoc.args=list(), dissCompare.args=list(),
                      verbose=TRUE){

  statBIC <- length(unique(group)) == 2

  dissassoc.args <- c(list(group=group, weights=weights),dissassoc.args)
  lapply.args <- c(list(X=disslist, FUN=dissassoc),dissassoc.args)
  assoc.list <- do.call(lapply, lapply.args)

  r2pval <- function(x){x$stat[3,,drop=FALSE]}
  mcomp <- matrix(unlist(t(sapply(assoc.list,r2pval))),ncol=2,byrow=FALSE)
  colnames(mcomp) <- c("R2","p(R2)")

  if (statBIC) {
    dissCompare.args <- c(list(group=group, weights=weights),dissCompare.args)
    if (is.null(dissCompare.args[["BFopt"]]))
      dissCompare.args[["BFopt"]]<-1
    if (is.null(dissCompare.args[["s"]]))
      dissCompare.args[["s"]] <- 0
    lapply.args <- c(list(X=disslist, FUN=TraMineRextras::dissCompare), dissCompare.args)
    #print(lapply.args)
    suppressMessages(comp.list <- do.call(lapply, lapply.args))
    mcomp2 <- matrix(unlist(comp.list),nrow=length(comp.list),byrow=TRUE)
    #colnames(mcomp) <- colnames(bf.comp.list[[1]])
    colnames(mcomp2) <- c("LRT","p(LRT)","Delta BIC","Bayes Factor")
    mcomp <- cbind(mcomp,mcomp2)
  } else {
    if (verbose) message("More than 2 groups, LRT and BIC not computed! ")
  }
  rm(dissassoc.args,dissCompare.args,lapply.args)

  return(mcomp)

}
