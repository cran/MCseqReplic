#############
printdistmat <- function(distmat,toref,n,nn,...){
  if (n<=0 | nn < n+1)
    print(distmat,...)
  else {
    cat("For first ",n," sequences\n")
    if(toref)
      print(distmat[1:n,], ...)
    else
      print(as.dist(as.matrix(distmat)[1:n,1:n]), ...)
    cat("...\n")
  }
}


#' Print method for distMC objects
#'
#' Prints, for each pair of the first \code{n} sequences, the mean and/or the standard deviation of the MC-replicated distances between sequences. When available, ratios are also printed by default.
#'
#' @param x \code{distMC} object as returned by \code{MCseqdistSE}.
#' @param n Integer. Number of first sequences. Default is 6. If \code{n==0} or there are less than \code{n} sequences, results are printed for all pairs of sequences.
#' @param what character string. One of \code{"mean"}, \code{"sd"}, \code{"bias"}, \code{"both"}, and \code{"all"} (default). When \code{"all"}, ratios, when present are printed together with the mean and standard deviation. When \code{"both"}, means and standard deviations are printed.
#' @param ... further arguments passed to or from other methods.
#'
#' @author Gilbert Ritschard
#' @seealso \code{\link{MCseqdistSE}}, \code{\link{summary.distMC}}.
#'
#' @return Last printed table, a \code{matrix} when \code{toref} attribute is \code{TRUE} and a \code{dist} object otherwise.
#'
#' @export
#'
#' @exportS3Method
#' print distMC

print.distMC <- function(x, n=6, what="all", ...){
  toref <- if (is.null(attr(x,"toref"))) FALSE else attr(x,"toref")
  toreftext <- " and refs,"
  if (toref){
    ns <- nrow(x[["MC.mean"]])
    toreftext <- "and refs,"
  } else {
    toreftext <- ","
    ns <- attr(x[["MC.mean"]],"Size")
  }
  N <- x[["N"]]
  R <- x[["R"]]
  RR <- if (is.null(R)) N else R*R
    cat(ns, " sequences,",RR," dissimilarity MC-replications\n")

  if (!what %in% c("all","both","mean","sd","obs","se","bias"))
    stop(' what should be one of "all","both","mean","sd","bias"')
  cat("Dissimilarity between MC-simulated sequences",toreftext,RR," dissimilarity replications")
  if (!is.null(x[["MC.sd"]])){
    if (attr(x,"obs") & what %in% c("obs","all") ){
      cat("\n diss.o: Observed dissimilarities\n")
      printdistmat(x[["diss.o"]],toref,n,ns,...)
    }
  }
  if (what %in% c("mean","both","all") ){
    cat("\n MC.mean: Mean of dissimilarities \n")
    printdistmat(x[["MC.mean"]],toref,n,ns,...)
  }
  if (!is.null(x[["MC.sd"]]) & what %in% c("sd","both","all") ){
    cat("\n MC.sd: Standard deviation of dissimilarities\n")
    printdistmat(x[["MC.sd"]],toref,n,ns,...)
  }
  if (!is.null(x[["MC.se"]]) & what %in% c("se","all") ){
    cat("\n MC.se: Standard eror of dissimilarities\n")
    printdistmat(x[["MC.se"]],toref,n,ns,...)
  }
  if (!is.null(x[["diss.z"]]) & what == "all"){
    cat("\n diss.z: Ratios diss/MC.se \n")
    printdistmat(x[["diss.z"]],toref,n,ns,...)
  }

  if (length(x[["mean.se"]]) > 0 & what == "all"){
    cat("\n MC.mean.z: Ratios MC.mean/mean.se \n")
    printdistmat(x[["MC.mean.z"]],toref,n,ns,...)

    cat("\n mean.se: Standard error of mean simulated dissimilarities \n")
    printdistmat(x[["mean.se"]],toref,n,ns,...)
  }
  if (!is.null(x[["MC.bias"]]) & what %in% c("all","bias")){
    cat("\n MC.bias: Bias (observed minus MC.mean) \n")
    printdistmat(x[["MC.bias"]],toref,n,ns,...)

#    cat("\n MC.mse: Mean squared error \n")
#    printdistmat(x[["MC.mse"]],toref,n,ns,...)
  }
}

#' Summary method for distMC objects
#'
#' Prints summary statistics of the observed dissimilarity \code{diss}, the mean \code{MC.mean}, standard deviation \code{MC.sd}, and standard error of dissimilarities between MC-replicated sequences, and the ratios  \code{diss/MC.se} and \code{MC.mean/MC.se}. Reported statistics concern all distances between original sequences.
#'
#' @param object \code{distMC} object as returned by \code{MCseqdistSE}.
#' @param ... further arguments passed to or from other methods.
#' @param silent logical: Should additional info be displayed?
#'
#' @author Gilbert Ritschard
#' @seealso \code{\link{MCseqdistSE}}, \code{\link{print.distMC}}
#'
#' @export
#'
#' @exportS3Method
#' summary distMC
#'
#' @return  \code{fivenumb} table with the statistics (min, Q1, med, Q3, max) of the observed dissimilarities, the mean, standard deviation, and standard error of the MC-simulated dissimilarities, standardized ratios, MC-bias and mean squared errors when available.



summary.distMC <- function(object, ..., silent=FALSE){
  ## need weights for stats
  toref <- attr(object,"toref")
  if (is.null(toref)) toref <- FALSE
  if (toref){
    n <- nrow(object[["MC.mean"]])
    nref <- ncol(object[["MC.mean"]])
  } else {
    n <- attr(object[["MC.mean"]],"Size")
  }
  w.u <- object[["weights"]]
  if (is.null(w.u) | length(w.u) == 1) w.u <- rep(1, n)

  #sd.w <- as.dist(w.u %*% t(w.u))
  sd.w <- as.dist(tcrossprod(w.u)) # same as above but faster

  diss.o <- object[["diss.o"]]
  R <- object[["R"]]
  N <- object[["N"]]
  #diss.z <- object[["diss.z"]]
  #MC.mean.z <- object[["MC.mean.z"]]

  if (toref){
    fivenumb <- list()
    for (i in 1:nref){
      if (any(object[["MC.mean"]] > 0) & n > 1 & !is.na(sum(w.u)) ) {
        fivenumb[[i]] <- matrix(TraMineR:::wtd.fivenum.tmr(diss.o[,i], weights=w.u), ncol=1)
        fivenumb[[i]] <- cbind(fivenumb[[i]],matrix(TraMineR:::wtd.fivenum.tmr(object[["MC.mean"]][,i], weights=w.u), ncol=1))
        cnames <- c("diss", "MC.mean")
        if (!is.null(object[["MC.sd"]])){
          fivenumb[[i]] <- cbind(fivenumb[[i]],matrix(TraMineR:::wtd.fivenum.tmr(object[["MC.sd"]][,i], weights=w.u), ncol=1))
          cnames <- c(cnames, "MC.sd")
        }
        if (!is.null(object[["mean.se"]])){
          fivenumb[[i]] <- cbind(fivenumb[[i]],matrix(TraMineR:::wtd.fivenum.tmr(object[["mean.se"]][,i], weights=w.u), ncol=1))
          fivenumb[[i]] <- cbind(fivenumb[[i]],matrix(TraMineR:::wtd.fivenum.tmr(object[["MC.mean.z"]][,i], weights=w.u), ncol=1))
          cnames <- c(cnames, "mean.se", "MC.mean/mean.se")
        }
        if (!is.null(object[["MC.bias"]])){
          fivenumb[[i]] <- cbind(fivenumb[[i]],matrix(TraMineR:::wtd.fivenum.tmr(object[["MC.se"]][,i], weights=w.u), ncol=1))
          fivenumb[[i]] <- cbind(fivenumb[[i]],matrix(TraMineR:::wtd.fivenum.tmr(object[["MC.bias"]][,i], weights=w.u), ncol=1))
          # fivenumb[[i]] <- cbind(fivenumb[[i]],matrix(TraMineR:::wtd.fivenum.tmr(object[["MC.mse"]][,i], weights=w.u), ncol=1))
          cnames <- c(cnames, "MC.se", "MC.bias")
        }
        colnames(fivenumb[[i]]) = cnames

        #rownames(fivenumb) = c("Lower notch","Lower hinge","Median","Upper hinge","Upper notch")
        rownames(fivenumb[[i]]) = c("Min","Q1","Median","Q3","Max")
      } else {
        fivenumb[[i]] <- NULL
      }
    }
    names(fivenumb) <- colnames(object[["MC.mean"]])
  } else {
    if (any(object[["MC.mean"]] > 0) & n > 1 & !is.na(sum(sd.w)) ) {
      fivenumb <- matrix(TraMineR:::wtd.fivenum.tmr(diss.o, weights=sd.w), ncol=1)
      fivenumb <- cbind(fivenumb,matrix(TraMineR:::wtd.fivenum.tmr(object[["MC.mean"]], weights=sd.w), ncol=1))
      cnames <- c("diss", "MC.mean")
      if (!is.null(object[["MC.sd"]])){
        fivenumb <- cbind(fivenumb,matrix(TraMineR:::wtd.fivenum.tmr(object[["MC.sd"]], weights=sd.w), ncol=1))
        cnames <- c(cnames, "MC.sd")
      }
      if (length(object[["mean.se"]])>0){
        fivenumb <- cbind(fivenumb,matrix(TraMineR:::wtd.fivenum.tmr(object[["mean.se"]], weights=sd.w), ncol=1))
        fivenumb <- cbind(fivenumb,matrix(TraMineR:::wtd.fivenum.tmr(object[["MC.mean.z"]], weights=sd.w), ncol=1))
        cnames <- c(cnames, "mean.se", "MC.mean/mean.se")
      }
      if (!is.null(object[["MC.bias"]])){
        fivenumb <- cbind(fivenumb,matrix(TraMineR:::wtd.fivenum.tmr(object[["MC.se"]], weights=sd.w), ncol=1))
        fivenumb <- cbind(fivenumb,matrix(TraMineR:::wtd.fivenum.tmr(object[["MC.bias"]], weights=sd.w), ncol=1))
        # fivenumb <- cbind(fivenumb,matrix(TraMineR:::wtd.fivenum.tmr(object[["MC.mse"]], weights=sd.w), ncol=1))
        cnames <- c(cnames, "MC.se", "MC.bias")
      }
      colnames(fivenumb) = cnames

      #rownames(fivenumb) = c("Lower notch","Lower hinge","Median","Upper hinge","Upper notch")
      rownames(fivenumb) = c("Min","Q1","Median","Q3","Max")
    } else {
      fivenumb <- NULL
    }
  }
  if (!silent) {
    if (is.null(N))
      cat(n," sequences, R =",R,", ", R^2, " MC-simulated dissimilarities per observed dissimilarity \n")
    else
      cat(n," sequences, N =",N," MC-simulated dissimilarities per observed dissimilarity \n")
  }

  return(fivenumb)
}


#' Print method for MCratios objects
#'
#' Prints ratios for each pair of the first \code{n} sequences.
#'
#' @param x \code{MCratios} object as returned by \code{MCratios}.
#' @param n Integer. Number of first sequences. Default is 6. If \code{n==0} or there are less than \code{n} sequences, results are printed for all pairs of sequences.
#' @param what character string. One of \code{"all"} (default), \code{"diss"}, \code{"mean"}, and \code{"se"} .
#' @param ... further arguments passed to or from other methods.
#'
#' @author Gilbert Ritschard
#' @seealso \code{\link{MCratios}}. %\code{\link{seqdistMCSE}},
#' @return Last printed table, a \code{matrix} when \code{toref} attribute is \code{TRUE} and a \code{dist} object otherwise.
#'
#' @export
#'
#' @exportS3Method
#' print MCratios

print.MCratios <- function(x, n=6, what="all", ...){

  toref <- attr(x,"toref")
  if(is.null(toref)) toref <- FALSE

  if (toref)
    nn <- nrow(x[["diss.z"]])
  else
    nn <- attr(x[["diss.z"]],"Size")
  if (!what %in% c("all","diss","mean","se"))
    stop(' what should be one of "all","diss","mean","se"')
  cat(nn, " sequences\n")
  cat("Dissimilarity between MC-simulated sequences")
  if (what %in% c("diss","all") ){
    cat("\n diss.z: Ratios diss/MC.se \n")
    printdistmat(x[["diss.z"]],toref,n,nn,...)
  }
  if (what %in% c("mean","all") ){
    cat("\n MC.mean.z: Ratios MC.mean/mean.se \n")
    printdistmat(x[["MC.mean.z"]],toref,n,nn,...)
  }
  if (what %in% c("se","all") ){
    cat("\n mean.se: Standard error of mean simulated dissimilarities \n")
    printdistmat(x[["mean.se"]],toref,n,nn)
  }
}


#' Summary method for MCratios objects
#'
#' Prints summary statistics of the ratios  \code{diss/MC.se} and \code{MC.mean/MC.se}. Reported statistics concern all distances between original sequences.
#'
#' @param object \code{MCratios} object as returned by \code{MCratios}.
#' @param weights vector of doubles. Case weights.
#' @param ... further arguments passed to or from other methods.
#' @param silent logical: Should additional info be displayed?
#' @param thresh real: threshold for counting ratios less than \code{thresh}
#'
#' @author Gilbert Ritschard
#' @seealso \code{\link{MCseqdistSE}}, \code{\link{print.distMC}}
#'
#' @export
#'
#' @exportS3Method
#' summary MCratios
#'
#' @return  \code{fivenumb} table with the statistics (min, Q1, med, Q3, max) of \code{mean.se} and the standardized ratios \code{diss.z} and \code{MC.mean.z}.

summary.MCratios <- function(object, ..., weights=NULL, silent=FALSE, thresh=2 ){
  toref <- attr(object,"toref")
  if(is.null(toref)) toref <- FALSE

  ## need weights for stats
  if (toref){
    n <-  nrow(object[["MC.mean.z"]])
    nc <- ncol(object[["MC.mean.z"]])
  } else {
    n <- attr(object[["MC.mean.z"]],"Size")
  }
  w.u <- weights
  #print(n)
  #print(w.u)
  if (is.null(w.u) | length(w.u) == 1) w.u <- rep(1, n)

  #sd.w <- as.dist(w.u %*% t(w.u))
  sd.w <- as.dist(tcrossprod(w.u)) # same as above but faster

  diss.o <- object[["diss.o"]]
  R <- object[["R"]]
  N <- object[["N"]]
  if (toref) {
    if (any(object[["mean.se"]] > 0) & n > 1 & !is.na(sum(w.u)) ) {
      fivenumb <- list()
      pLTt <- nLTt <- NULL
      ndiss <- n
      for (i in 1:nc) {
        fivenumb[[i]] <- matrix(TraMineR:::wtd.fivenum.tmr(object[["mean.se"]][,i], weights=w.u), ncol=1)
        fivenumb[[i]] <- cbind(fivenumb[[i]],matrix(TraMineR:::wtd.fivenum.tmr(object[["diss.z"]][,i], weights=w.u), ncol=1))
        fivenumb[[i]] <- cbind(fivenumb[[i]],matrix(TraMineR:::wtd.fivenum.tmr(object[["MC.mean.z"]][,i], weights=w.u), ncol=1))
        colnames(fivenumb[[i]]) = c("mean.se", "diss/MC.se", "MC.mean/mean.se")
        nnLTt <- sum(object[["diss.z"]][,i] < thresh,na.rm=TRUE)
        pLTt <- c(pLTt,nnLTt/n)
        nLTt <- rbind(nLTt,nnLTt)

        #rownames(fivenumb) = c("Lower notch","Lower hinge","Median","Upper hinge","Upper notch")
        rownames(fivenumb[[i]]) = c("Min","Q1","Median","Q3","Max")
      }
    } else {
      fivenumb <- NULL
      warning("summary not applicable!")
    }
    names(fivenumb) <- colnames(object[["mean.se"]])
  } else {
    if (any(object[["mean.se"]] > 0) & n > 1 & !is.na(sum(sd.w)) ) {
      fivenumb <- matrix(TraMineR:::wtd.fivenum.tmr(object[["mean.se"]], weights=sd.w), ncol=1)
      fivenumb <- cbind(fivenumb,matrix(TraMineR:::wtd.fivenum.tmr(object[["diss.z"]], weights=sd.w), ncol=1))
      fivenumb <- cbind(fivenumb,matrix(TraMineR:::wtd.fivenum.tmr(object[["MC.mean.z"]], weights=sd.w), ncol=1))
      colnames(fivenumb) = c("mean.se", "diss/MC.se", "MC.mean/mean.se")
      nLTt <- sum(object[["diss.z"]] < thresh,na.rm=TRUE)
      pLTt <- nLTt/(n*(n-1)/2)
      ndiss <- n*(n-1)/2
      #rownames(fivenumb) = c("Lower notch","Lower hinge","Median","Upper hinge","Upper notch")
      rownames(fivenumb) = c("Min","Q1","Median","Q3","Max")
    } else {
      fivenumb <- NULL
      warning("summary not applicable!")
    }
  }
  if (!silent) {
    if (is.null(N))
      message(n," sequences, R = ",R,", ", R^2, " MC-simulated dissimilarities per observed dissimilarity")
    else
      message(n," sequences, N = ",N," MC-simulated dissimilarities per observed dissimilarity")
  }

  return(list(fivenumb=fivenumb,prop.LT.thresh=pLTt,n.LT.thresh=nLTt,ndiss=ndiss,thresh=thresh))
}

