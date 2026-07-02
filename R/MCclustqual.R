#' @title Cluster quality measures by MC-sets
#'
#' @description
#' Cluster quality measures for a range of number of groups by MC-replicated set.
#'
#'
#' @param disslist List of MC-dissimilarity matrices (or \code{dist} objects).
## #' @param diss.o
#' @param ncluster integer vector. Maximum number of groups. Default is \code{10}. CQIs are computed for the range \code{2:ncluster}
#' @param clustmeth character. Clustering method. Either \code{"PAM"} (default) or a \code{stats::\link{hclust}} method.
#' @param weights vector of doubles. Case weights. If \code{NULL} (default), equal weights are used.
#' @param core Integer or \code{"auto"}. Number of cores for parallel computing If \code{"auto"}, the maximum available cores are used.
#' @param snow Logical. If \code{TRUE}, \code{doSNOW} is used for parallel computing, otherwise \code{doParallel} is used.
#' @param silent Logical. Deprecated, use \code{!verbose} instead!
#' @param verbose Logical. Should waiting and timing messages be printed?
#' @param ... Further arguments passed to clustering or plot functions.
#'
## #' @param what String. One of the indexes computed by \code{aricode} package such as {\code{"ARI"} (adjusted Rand index, default), \code{"RI"} (Rand index), \code{"VI"} (variation of information), \code{"NVI"}. Can also be "all \code{"all"}.
#'
#' @details
#' When \code{attr(MCdisslist,"obs")} is \code{TRUE}, the last element of \code{disslist} is treated as the dissimilarity matrix of the observed sequences.
#'
#' \code{MCclustqual} computes the range of CQI values for all the CQIs included in the \code{stats} element returned by \code{WeightedCluster::\link[WeightedCluster]{wcClusterQuality}}.
#'
#' @return List of length 3: \cr
#' - \code{qual.tab}: list of tables of cluster quality statistics per MC-dissimilarity matrix,\cr
#' - \code{qual.max}: table of cluster number $k$ for which the statistics reach their maximum (minimum for HC) by MC-sets and observed sequence set (rows),\cr
#' - \code{max.freq}: the frequency table of optimal $k$ over the MC-replicated sets, and\cr
#'
#' @seealso \code{\link[WeightedCluster]{as.clustrange}},  \code{\link[WeightedCluster]{wcKMedRange}}
#'
#' @references
#'   \insertNoCite{Studer2013}{MCseqReplic}
#'   \insertAllCited{}
#'
#'
#'
#' @importFrom WeightedCluster wcKMedRange as.clustrange
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makePSOCKcluster stopCluster makeCluster detectCores
#' @importFrom parallelly availableCores
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach %dopar% %do% %:%
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats hclust as.dist
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
#' res <- MCclustqual(disslist,ncluster=3, verbose=FALSE)
#' res
#' ggplotMCcqi(res,"PBC")
#'
#'
MCclustqual <- function(disslist, ncluster=10, clustmeth="PAM",  weights=NULL,
                        core=1, snow=TRUE, verbose=!silent, silent=FALSE, ...){

  N <- length(disslist)
  isobs <- attr(disslist,"obs")
  if (is.null(isobs)) isobs <- FALSE

  if (inherits(disslist[[1]],"dist")){
    ncases <- attr(disslist[[1]],"Size")
  } else {
    ncases <- nrow(disslist[[1]])
  }
  ncluster <- min(ncases,ncluster)
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

  core <- set_ncores(core)

  if (verbose & (core==1 | snow)) {
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
    if (verbose) close(pb)
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
    if(verbose & snow) {
      close(pb)
    }
  }

  if (isobs){
    irange <- 1:(N-1)
    #qobs <- qlist[[N]]
    #qlist <- qlist[-N]
    #N <- N-1
  } else {
    #qobs <- NULL
    irange <- 1:N
  }

  which.max.vect <- function(m) {
    stats <- m
    stats[,"HC"] <- -stats[,"HC"] ## For HC we seek the min
    apply(stats,2,which.max) + 1
    }
  max.list <- lapply(qlist, which.max.vect)

  ## distribution of max over MC-sets
  tabmax <- do.call(rbind,max.list)
  tabmaxf <- data.frame(factor(tabmax[,1],levels=2:ncluster))
  #print(tabmaxf)
  maxdist <- data.frame(table(tabmaxf[irange,1]))
  #print(maxdist)
  rownamesd <- maxdist[,1]
  maxdist <- maxdist[,2]
  for (i in 2:ncol(tabmax)){
    tabmaxf <- cbind(tabmaxf,factor(tabmax[,i],levels=2:ncluster))
    maxdist <- cbind(maxdist,data.frame(table(tabmaxf[irange,i]))[,2])
  }
  rownames(maxdist) <- rownamesd
  colnames(tabmaxf) <- colnames(maxdist) <- colnames(tabmax)
  rownames(tabmaxf) <- paste0("MC",1:N)
  if (isobs) rownames(tabmaxf)[N] <- "Obs"

  ret <- list(qual.tab = qlist, qual.max = tabmaxf, max.freq=maxdist)

  class(ret) <- c(class(ret),"MCclustQ")
  attr(ret,"obs") <- attr(disslist,"obs")
  return(ret)
}


##################

#' @rdname MCclustqual
#'
## Plot range of values of a selected CQI by MC-sets
#'
#' @description \code{ggplotMCcqi} makes a ggplot of the range of values of the selected CQI by MC-sets and for the observed sequences. When \code{attr(data,"obs")} is \code{TRUE}, the range of CQI values for the observed sequences is also plotted.
#'
#' @param data an \code{MCclustQ} object as returned by \code{MCclustqual}
#' @param cqi string. The name of the selected CQI.
#' @param meancqi logical. Should the range of mean values of the selected CQI be plotted?
#' @param scalelwd double. Line width scale value.
#' @param linecolor vector of three line colors in the order Mean, Obs, MCset. If \code{NULL}, default colors are used and if of length less than 3, default colors are used for the first elements.
## @param ... Further \code{ggplot} arguments.
#'
#' @import ggplot2
#' @importFrom tidyr gather
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#' @export
#'
#' @returns \code{ggplotMCcqi} returns the ggplot object.
#'
#'
ggplotMCcqi <- function(data, cqi="PBC", meancqi = TRUE, scalelwd = 1, linecolor = NULL,  ...){
  obs <- attr(data,"obs")
  extract.cqi <- function(m, cqi) {
    m[[cqi]]
  }
  cqi.list <- lapply(data$qual.tab, extract.cqi, cqi=cqi)
  cqi.data <- t(do.call(rbind,cqi.list))
  nr <- nrow(cqi.data)
  nc <- ncol(cqi.data)
  colnames(cqi.data) <- paste0("MC",1:nc)

  ccol <- c("black","red","steelblue")
  if (length(linecolor) == 3)
    ccol <- linecolor
  else if (length(linecolor) == 2)
    ccol[2:3] <- linecolor
  else if (length(linecolor) == 1)
    ccol[1] <- linecolor
  else if (!is.null(linecolor))
    stop ("length(linecolor) > 3")

  clab <- c("Mean","Obs","MCsets")
  colsel <- c(meancqi,obs,TRUE)

  if (meancqi){
    cqi.data <- cbind(cqi.data, Mean = apply(cqi.data[,1:(nc-obs)], 1, mean))
  }

  k=rep(1:nrow(cqi.data) + 1,nc + meancqi)
  mycol <- c(rep(ccol[3],length(k)))
  size <- rep(scalelwd,length(k))
  if (obs){
    colnames(cqi.data)[nc] <- "Obs"
    mycol[(length(k)-(1+meancqi)*nr+1):(length(k)-meancqi*nr)] <- rep(ccol[2],nr)
    size[(length(k)-(1+meancqi)*nr+1):(length(k)-meancqi*nr)] <- rep(scalelwd*1.3,nr)
  }
  if (meancqi){
    mycol[(length(k)-nr+1):(length(k))] <- rep(ccol[1],nr)
    size[(length(k)-nr+1):(length(k))] <- rep(scalelwd*1.3,nr)
  }

  cqi.long <- cqi.data %>%
    as.data.frame %>%
    gather(key="key", value="value") ##%>%
  ##  mutate(k = k, mycol = mycol, size = size)

  ggp <-ggplot(cqi.long, aes(x=k, y=.data$value, key=.data$key, colour=mycol, ...)) +
    geom_line(linewidth=size) +
    scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(min(x), (max(x) + 1) * 1.1))))) +
    scale_color_manual(name=cqi, labels = clab[colsel], values = ccol[colsel]) +
    ggtitle(paste(cqi, "by number k of clusters"))

  ggp
}


##### Print method for object MCclustQ
#' @rdname MCclustqual
#'
## Print method for MCratios objects
#'
#' @description The  print method only prints by default the \code{qual.max} and \code{max.freq} tables of MCclustQ objects.
#'
#' @param x \code{MCclustQ} object as returned by \code{MCclustqual}.
#' @param all logical. Should tables by MC-sets also be printed? Default is \code{FALSE}.
#' @param nMC numeric. Maximal number of MC-sets for which optimal size by CQIs are printed. Default is 5.
#' @param ... further arguments passed to or from other methods.
#'
#' @author Gilbert Ritschard
#' @seealso \code{\link{MCclustqual}}.
#' @return the print method returns the last printed tables.
#'
#' @export
#'
#' @exportS3Method
#' print MCclustQ

print.MCclustQ <- function(x, all=FALSE, nMC=5, ...){
  if (all) print(x[1])
  if (nrow(x$qual.max) < (nMC + attr(x,"obs"))) {
    nMC <- nrow(x$qual.max) - attr(x,"obs")
    cat("$qual.max\n")
  } else {
    cat(paste("$qual.max (first",nMC,"MC-sets)\n"))
  }
  if (attr(x,"obs"))
    print(x[[2]][c(1:nMC,nrow(x$qual.max)),],...)
  else
    print(x[[2]][1:nMC,],...)
  cat(" \n")
  print(x[3], ...)
}



