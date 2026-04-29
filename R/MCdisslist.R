#' @title List of dissimilarity matrices
#'
#' @description
#'  Compute the dissimilarity matrix for each of the provided sets of sequences.
#'
#' @param MCrseqdata List of state sequence objects of class \code{stslist}.
#' @param method string. Name of a distance method (see \code{\link[TraMineR]{seqdist}}).
#' @param seqref state sequence object of class \code{stslist}. Fixed reference sequences.
#' @param use.udiss logical. Should computation be based on unique sequences?
#' @param full.matrix logical. Should pairwise distances be returned in matrix form? If \code{FALSE} (default),
#' a list of \code{dist} objects is returned. Applies only when \code{seqref=NULL}.
#' @param ... further arguments passed to \code{seqdist}.
#'
#' @details
#' When \code{use.udiss=TRUE}, the function first computes dissimilarities between unique merged replicated sequences through a single call to \code{seqdist()} and the set of dissimilarity matrices are then extracted from the resulting distance matrix. This is generally faster when the number of unique merged replicated sequences is less than \code{sqrt(number of replicated datasets) * (sample size)}, which can be checked with \code{\link{MCnunique}}.
#'
#' @return list of dissimilarity matrices or \code{dist} objects with logical attribute \code{"obs"}, which is \code{TRUE} when the list includes the dissimilarities between observed sequences as last element.
#' @seealso \code{\link{MCseqReplicate}}, \code{\link{MCudist}} and examples in their help pages.
#' @export

MCdisslist <- function(MCrseqdata, method="LCS", seqref=NULL, full.matrix=FALSE, use.udiss=FALSE, ...){
  if (is.null(seqref)){
    if (use.udiss){
      udiss <- suppressMessages(MCudist(MCrseqdata, method=method, full.matrix=full.matrix, ...))
      disslist <- list()
      for (k in 1:length(MCrseqdata)){
        disslist[[k]] <- MCExtractDist(udiss,k,full.matrix=full.matrix)
      }
    } else {
      disslist <- suppressMessages(lapply(MCrseqdata, seqdist, method=method,
                                        full.matrix=full.matrix, ...))
    }
    attr(disslist,"toref") <- FALSE
  } else { ## distances to fixed reference sequences
    if (use.udiss){
      udiss <- suppressMessages(MCudist(MCrseqdata, method=method, seqref=seqref, ...))
      disslist <- list()
      for (k in 1:length(MCrseqdata)){
        disslist[[k]] <- MCExtractDist(udiss,k)
      }

    } else {
    disslist <- suppressMessages(lapply(MCrseqdata, seqdistToRef,
                                        ref = seqref,
                                        method=method,
                                        ...))
    }
    attr(disslist,"toref") <- TRUE
  }
  attr(disslist,"obs") <- attr(MCrseqdata,"obs")
  return(disslist)
}

seqdistToRef <- function(s, ref, method, ...){
  nseq <- nrow(s)
  nref <- nrow(ref)
  rownames(ref) <- paste0("ref.",rownames(ref))
  seqr <- rbind(s,ref)
  dist <- seqdist(seqr, refseq = list(1:nseq,nseq+(1:nref)), method=method, ...)
}
