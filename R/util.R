#' bsxfun {pracma} with single expansion (real Matlab style)
#' @param func the function used by bsxfun
#' @param x a matrix
#' @param y a vector need to be expanded
#' @param  expandByRow applies only when x is a square matrix
#' @return value of func
#' @export
bsxfun.se <- function(func, x, y, expandByRow=TRUE) {

  if(length(y) == 1) return(pracma::arrayfun(func, x, y)) else
    stopifnot(nrow(x) == length(y) || ncol(x) == length(y))

  expandCol <- nrow(x) == length(y)
  expandRow <- ncol(x) == length(y)

  if(expandCol & expandRow & expandByRow) expandCol <- FALSE
  if(expandCol & expandRow & !expandByRow) expandRow <- FALSE

  # repeat row (if dim2expand = 1, then length(y) = ncol(x))
  if(expandRow) y.repmat <- matrix(rep(as.numeric(y), each=nrow(x)), nrow=nrow(x))

  # repeat col (if dim2expand = 2, then length(y) = nrow(x))
  if(expandCol) y.repmat <- matrix(rep(as.numeric(y), ncol(x)), ncol=ncol(x))

  pracma::bsxfun(func, x, y.repmat)
}

#' Compute log(sum(exp(x),dim)) while avoiding numerical underflow
#' @param x a matrix
#' @param margin used for apply
#' @return log(sum(exp(x),dim)) a matrix of the sample size as x
logsumexp <- function(x, margin=1) {

  stopifnot(is.matrix(x))

  # subtract the largest in each column
  y <- apply(x, margin, max)

  if (nrow(x) == ncol(x)) {
    x <- bsxfun.se("-", x, y, expandByRow = F)
  } else {
    x <- bsxfun.se("-", x, y)
  }

  s <- y + log(apply(exp(x), margin, sum))

  i <- which(!is.finite(s))

  if(length(i) > 0) s[i] <- y[i]

  s
}


#' Convert ccf to vaf
#' @param x ccf
#' @param t purity
#' @param cn normal total copy number
#' @param cr total copy number in reference population
#' @param cv total copy number in variant population
#' @param bv mutation multplicity
#' @param epi sequencing error
#' @return vaf
#' @export
cp2ap <- function(x, t, cn, cr, cv, bv, epi = 1e-3){
  p1 <- (1-t)*cn
  p2 <- t*(1-x)*cr
  p3 <- t*x*cv
  p <- p1+p2+p3
  u1 <- epi
  u2 <- epi
  u3 <- (bv/cv)*(1-epi)
  y <- p1*u1+p2*u2+p3*u3
  return(y/p)
}

#' Convert vaf to ccf using pyclone style mapping
#' @param x vaf
#' @param t purity
#' @param cn normal total copy number
#' @param cr total copy number in reference population
#' @param cv total copy number in variant population
#' @param bv mutation multplicity
#' @param epi sequencing error
#' @param constraint whether constriant the transformation
#' @param lower lower bound of the constriant
#' @param upper upper bound of the constriant
#' @return ccf
#' @export
MapVaf2CcfPyClone <- function(x, t, cn, cr, cv, bv,
                              epi = 1e-3, constraint = T,
                              lower = -Inf, upper = Inf) {
  un <- epi
  ur <- epi
  if (bv == cv) {
    uv <-  1-epi
  } else if ( bv == 0) {
    uv <- epi
  } else {
    uv <- (bv/cv)*(1-epi)
  }

  ccf <- ( (1-t) * cn * (un - x) + t * cr *(ur - x) ) /
    ( t * cr * (ur - x) - t * cv * (uv - x) )

  if (constraint) {
    if (is.na(ccf)) {
      return(as.numeric(NA))
    } else if (ccf < 0.9 && bv > 1) {
      return(as.numeric(NA))
    } else if (ccf < upper && ccf > lower) {
      return(ccf)
    } else {
      return(as.numeric(NA))
    }
  } else {
    return(ccf)
  }
}

############ Vector dot product ############
# handle single row matrix by multiplying each value
# but not sum them up
dot.ext <- function(x,y,mydim) {

  if(missing(mydim)) pracma::dot(x,y) else {

    if(1 %in% pracma::size(x) & mydim == 1) x * y else pracma::dot(x,y)
  }
}

############ Logarithmic Multivariate Gamma function ############
# Compute logarithm multivariate Gamma function.
# Gamma_p(x) = pi^(p(p-1)/4) prod_(j=1)^p Gamma(x+(1-j)/2)
# log Gamma_p(x) = p(p-1)/4 log pi + sum_(j=1)^p log Gamma(x+(1-j)/2)
logmvgamma <- function(x, d) {

  s <- pracma::size(x)

  x <- matrix(as.numeric(x), nrow=1)

  x <- bsxfun.se("+", kronecker(matrix(1,d,1), x), (1 - matrix(1:d))/2)

  y <- d*(d-1)/4*log(pi) + colSums(lgamma(x))

  y <- matrix(as.numeric(y), nrow=s[1], ncol=s[2])

  y
}
