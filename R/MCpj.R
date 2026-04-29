#' @title Generate distribution of timing errors
#'
#' @description
#' Generates a distribution of timing errors that complies with the provided expected size of non-zero timing errors and the expected probability of no error.
#'
#' @details
#' Currently \code{\link{MCseqReplicate}} expects a vector Pj with same number of backward and forward error values. To comply with this, the shorter side of Pj is by default filled with zeros.
#'
#'
#'
#' @param Emean scalar or vector of size two. Expected size of non-zero timing errors. If a vector, the first value is used for negative errors and the second value for positive errors. If a scalar, the value is used for both negative and positive errors. Values must be strictly greater than 1.
#' @param pzero number in range [0,1]. Probability of no-error. If \code{NULL} (default), \code{pzero} is set to the the greatest probability of zero between the right and left side Poisson distributions.
#' @param maxterr integer. Maximal error size to consider. Default is 10.
#' @param pinterv control value used for solving numerically an implicit function. Default is .99 and should be increased in case the zero of the implicit function cannot be found because of ending values of same sign.
#' @param fill.short.side logical. Should the shortest side be filled with zeros to equal length of the other side. Default is \code{TRUE}.
#'
#' @returns The vector of probabilities Pj with the computed \code{lambda} values as attribute.
#'
#' @seealso \code{\link{MCseqReplicate}}
#'
#' @export
#'
#' @importFrom stats dpois uniroot
#'
#' @examples
#' # expected timing error of 1.2 on each side
#' MCpj(Emean=1.2, pzero=.4)
#'
#' # expected backward timing error higher than for forward errors
#' MCpj(Emean=c(3.5,1.2), pzero=.4)
#'
#'

MCpj <- function(Emean, pzero=NULL, maxterr=10, pinterv=.99, fill.short.side=TRUE){

  if (length(Emean) > 2){
    warning("length of Emean greater than 2, only first 2 values are used!")
    vEmean <- Emean[1:2]
  }

  if (length(Emean)==2){
    nEmean <- 2
  } else {
    nEmean <- 1
  }

  lpois <- list()
  vEmean <- Emean
  pz <- 0
  lbda <- NULL
  for (i in 1:nEmean) {
    Emean <- vEmean[i]

    if (Emean <= 1){
      stop("Emean must be strictly greater than 1")
    }

    # Define the implicit function F(x, y) = 0
    implicit_fun <- function(lambda, expected) {
      return(expected - lambda/(1-exp(-lambda)))
    }

    # Find the root for lambda in the interval around the expected value
    interv <- c(Emean - Emean * pinterv, Emean + Emean * pinterv)
    # uniroot() requires the function to be defined such that the interval endpoints
    # produce values with opposite signs.

    tryCatch({
      solution_lambda <- uniroot(f = implicit_fun, interval = interv, expected = Emean)
    }, error = function(e) {
      print(paste("Error finding root:", e$message))
    })

    lambda <- max(solution_lambda$root,.001)
    pois <- dpois(1:maxterr, lambda)
    pz <- max(pz, dpois(0, lambda))
    lbda <- c(lbda, lambda)

    ## assuming round(pois,2)[1:floor(lambda)] > 0

    #pois.c <- pois[(floor(lambda) + 1):length(pois)]
    #rpois <- c(pois[1:max(floor(lambda),1)],pois.c[round(pois.c,2) > 0])

    lpois[[i]] <- pois
    if (length(pois) > ceiling(lambda)) {
      pois.hend <- pois[(ceiling(lambda) + 1):length(pois)]
      # Keep only ending values with prob > .005
      lpois[[i]] <- c(pois[1:ceiling(lambda)],  pois.hend[pois.hend >.005])
    }
  }

  if (nEmean == 1){
    ## lpois[[2]] id prob on rigth side
    lpois[[2]] <- lpois[[1]]
  }
  if (fill.short.side & length(lpois[[1]]) != length(lpois[[2]])){
    # make both vectors the same length
    lgth <- c(length(lpois[[1]]),length(lpois[[2]]))
    idmin <- which.min(lgth)
    zeros <- rep(0, abs(lgth[1]-lgth[2]))
    lpois[[idmin]] <- c(lpois[[idmin]],zeros)
  }

  ## for left side we reverse order of prob
  lpois[[1]] <- lpois[[1]][length(lpois[[1]]):1]

  ## if pzero is NULL, we set it to the computed value pz
  if (is.null(pzero)) pzero <- pz

  ## normalize to sum up to 1
  sumpois <- sum(c(lpois[[1]],lpois[[2]]))
  lpois[[1]] <- lpois[[1]] * (1-pzero)/sumpois
  lpois[[2]] <- lpois[[2]] * (1-pzero)/sumpois


  Pj <- c(lpois[[1]],pzero,lpois[[2]])
  attr(Pj, "lambda") <- lbda

  return(Pj)
}



