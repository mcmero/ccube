############ Sort the model  ############
# Sort model paramters in increasing order of averaged means
# of d variables
sort_components_smm <- function(model) {

  idx <- order(apply(model$normalMean, 2, mean))

  model$normalMean <- model$normalMean[, idx]

  model$responsibility <- model$responsibility[, idx]

  model$logResponsibility <- model$logResponsibility[, idx]

  model$dirichletConcentration <- model$dirichletConcentration[idx]

  model$normalRelativePrecision <- model$normalRelativePrecision[idx]

  model$whishartDof <- model$whishartDof[idx]

  model$invWhishartScale <- model$invWhishartScale[,,idx]

  model
}


#' Variational Bayesian Student-t mixture model (VB-SMM)
#' @param data N x D data matrix
#' @param init k (1 x 1) or label (1 x n, 1<=label(i)<=k) or center (d x k)
#' @param prior prior parameters
#' @param tol VBEM convergence threshold
#' @param maxiter VBEM maximum iteration
#' @param verbose show progress
#' @return a list containing model parameters
#' @export
vbsmm <- function(data, init=2, prior, tol=1e-20, maxiter=1e3, verbose=FALSE) {

	data <- as.matrix(data)

	n <- nrow(data)
	d <- ncol(data)

	X <- t(data) # Work with D by N for convenience

	message(sprintf("Running VB-SMM on a %d-by-%d data with %d clusters ...\n", n, d, init))

	if(missing(prior)) {

		if( length(init) == 1 ) {
			# more general prior with equal dirichletConcentration
		  k <-  init
		  ranged <- range(data)
		  prior <- list(
		    dirichletConcentration = rep(1e-2,k),
		    normalRelativePrecision = 1,
		    normalMean = t(as.matrix(rep(0.5, k))),
		    #normalMean = t(as.matrix(seq(ranged[1],ranged[2],length.out = k))),
		    whishartDof = d,
		    invWhishartScale = var(data)/k, # M = inv(W)
		    gammaDof = rep(1,k)
			)
		}
	} else {

		stopifnot(
			all(names(prior) %in% c("dirichletConcentration","normalRelativePrecision","normalMean","whishartDof","invWhishartScale")) &
			all(sapply(prior, is.numeric)) & nrow(prior$normalMean) == d &
			ncol(prior$normalMean) == 1 &
			nrow(prior$normalMean) == d & ncol(prior$invWhishartScale) == d)
	}

	# lower variational bound (objective function)
	L <- rep(-Inf, maxiter)
	converged <- FALSE
	t <- 1

	model <- list()

	initParams <- initialization_smm(X, init, prior) # initialize responsibility and hidden scale
	model$responsibility <- initParams$R
	model$gammaAlpha <- initParams$gammaAlpha
	model$gammaBeta <- initParams$gammaBeta
	model$hiddenScale <- initParams$U
	model$gammaDof <- prior$gammaDof

	while(!converged & t < maxiter) {
	  t <- t + 1
		model <- VariationalMaximimizationStep_smm(X, model, prior)
		model <- VarationalExpectationStep_smm(X, model)
		L[t] <- variationalLowerBound_smm(X, model, prior)/n
		converged <- abs(L[t] - L[t-1]) < tol * abs(L[t])
		if(verbose) message(sprintf("VB-EM-%d: L = %.6f \r", t, L[t]))
	}

	L <- L[2:t]

	model <- sort_components_smm(model)

	if (init > 1) {
	  label <- rep(0, n)
	  label <- apply(model$responsibility, 1, which.max)
	  nk <- colSums(model$responsibility)
	  Epi <- (model$dirichletConcentration + nk) / (k*prior$dirichletConcentration + n)
	  model$Epi <- Epi/sum(Epi)
	} else {
	  label <- rep(1, n)
	  model$Epi <- 1
	}


	# unique to have consecutive label eg 2,3,6 changed to 1,2,3
	# label <- match(label, sort(unique(label)))

	if(converged) message(sprintf("Converged in %d steps.\n", t-1)) else
	  warnings(sprintf("Not converged in %d steps.\n", maxiter))

	list(label=label, R=model$responsibility, mu=model$normalMean, full.model=model, L=L)
}


############ Initialization of responsibility (intialization) ############

initialization_smm <- function(X, init, prior) {

  d <- nrow(X)

  n <- ncol(X)

  stopifnot(length(init) %in% c(1, n) ||
              (nrow(init) == d  & ncol(init) == k))

  if(length(init) == 1) { # init = k gaussian components

    k <- init

    idx <- sample(1:n, k)

    m <- X[,idx,drop=F]


    res <- kmeans(t(X), k)
    label <- res$cluster
    normalMean <- t(res$centers)

    R <- as.matrix(Matrix::sparseMatrix(1:n, label, x=1))

  } else {

    if(length(init) == n) { # initialize with labels

      label <- init
      k <- max(label)
      R <- as.matrix(sparseMatrix(1:n, label, x=1))

    } else {

      if(!is.null(dim(init))) {

        if(nrow(init) == d  & ncol(init) == k) { # initialize with centers

          k <- ncol(init)
          m <- init
          label <- apply(bsxfun.se("-", pracma::crossprod(m, X),
                                   as.matrix(dot.ext(m,m,1)/2)), 2, which.max)
          R <- as.matrix(Matrix::sparseMatrix(1:n, label, x=1))

        } else stop(message("Invalid init."))
      }
    }
  }

  gammaDof <- prior$gammaDof
  dgammaDof2 <- (d+gammaDof)/2
  gammaBeta <- array(0, dim=c(n,k))
  invWhishartScale <- prior$invWhishartScale
  whishartDof <- prior$whishartDof
  normalRelativePrecision <- prior$normalRelativePrecision

  for(i in 1:k) {

    U <- chol(invWhishartScale)

    Q <- solve(t(U), bsxfun.se("-", X, normalMean[,i]))

    QQ <- dot.ext(Q, Q, 1)

    gammaBeta[,i] <- (d/(normalRelativePrecision) + whishartDof *QQ + gammaDof[i])/2
  }

  gammaAlpha <- dgammaDof2

  U = bsxfun.se("*", 1/gammaBeta, gammaAlpha)

  return(list(R = R, gammaAlpha = gammaAlpha, gammaBeta = gammaBeta, U = U))
}


############ Variational-Maximimization ############
VariationalMaximimizationStep_smm <- function(X, model, prior) {

  dirichletConcentration0 <- prior$dirichletConcentration
  normalRelativePrecision0 <- prior$normalRelativePrecision
  normalMean0 <- prior$normalMean
  whishartDof0 <- prior$whishartDof
  invWhishartScale0 <- prior$invWhishartScale
  responsibility <- model$responsibility
  hiddenScale <- model$hiddenScale

  scaledResponsibility <- responsibility * hiddenScale

  dataSumResponsibility <- colSums(responsibility) # 10.51
  dataSumScaledResponsibility <- colSums(scaledResponsibility)
  dirichletConcentration <- dirichletConcentration0 + dataSumResponsibility # 10.58
  scaledResponsibilityWeightedX <- X %*% scaledResponsibility
  normalRelativePrecision <- normalRelativePrecision0 + dataSumScaledResponsibility # 10.60

  normalMean <- bsxfun.se("*", bsxfun.se("+", scaledResponsibilityWeightedX, normalRelativePrecision0 * normalMean0), 1/normalRelativePrecision) # 10.61
  whishartDof <- whishartDof0 + dataSumResponsibility # 10.63 (NB: no 1 in the matlab code)

  dimensionOfData <- nrow(normalMean)
  numberOfComponents <- ncol(normalMean)

  invWhishartScale <- array(0, c(dimensionOfData, dimensionOfData, numberOfComponents))
  sqrtScaledResponsibility <- sqrt(scaledResponsibility)

  averageScaledResponsibilityWeightedX <- bsxfun.se("*", scaledResponsibilityWeightedX, 1/dataSumScaledResponsibility) # 10.52
  averageScaledResponsibilityWeightedXNormalMean0 <- bsxfun.se("-", averageScaledResponsibilityWeightedX, normalMean0)
  rescaledRelativePrecision <- (normalRelativePrecision0 * dataSumScaledResponsibility) / normalRelativePrecision

  for(i in 1:numberOfComponents) {

    Xs <- bsxfun.se("*", bsxfun.se("-", X, averageScaledResponsibilityWeightedX[,i]), t(sqrtScaledResponsibility[,i]))

    averageScaledResponsibilityWeightedXNormalMean0i <- averageScaledResponsibilityWeightedXNormalMean0[,i] # 10.62

    invWhishartScale[,,i] <- invWhishartScale0 +  Matrix::tcrossprod(Xs, Xs) +
      rescaledRelativePrecision[i] *
      Matrix::tcrossprod(averageScaledResponsibilityWeightedXNormalMean0i, averageScaledResponsibilityWeightedXNormalMean0i)
  }

  gammaAlpha <- model$gammaAlpha
  gammaBeta <- model$gammaBeta
  expectedlogHiddenScale <- bsxfun.se("+", -log(gammaBeta), digamma(gammaAlpha))
  tmp <- colSums(responsibility*(expectedlogHiddenScale - hiddenScale)) / dataSumResponsibility

  gammaDof <- model$gammaDof

  for (i in 1:numberOfComponents) {
    gammaDof[i] <- suppressWarnings( uniroot(
      function(x) { log(x/2) + 1 - digamma(x/2) + tmp[i] },
      c(2,200), extendInt = "yes")$root )
  }

  model$dirichletConcentration <- dirichletConcentration
  model$normalRelativePrecision <- normalRelativePrecision
  model$normalMean <- normalMean
  model$whishartDof <- whishartDof
  model$invWhishartScale <- invWhishartScale # Whishart: M = inv(W)
  model$gammaDof <- gammaDof

  model
}


############ Variational-Expectation ############

VarationalExpectationStep_smm <- function(X, model, no.weights = FALSE) {

  dirichletConcentration <- model$dirichletConcentration	# Dirichlet
  normalRelativePrecision <- model$normalRelativePrecision	# Gaussian
  normalMean <- model$normalMean			# Gasusian
  whishartDof <- model$whishartDof			# Whishart
  invWhishartScale <- model$invWhishartScale	# Whishart: inv(W) = V'*V

  gammaDof <- model$gammaDof # Gamma hyper

  if (!is.array(normalMean) ){
    normalMean <- array(normalMean, c(1, length(normalMean)))
  }

  if (!is.array(invWhishartScale) ){
    invWhishartScale <- array(invWhishartScale, c(1, 1, length(invWhishartScale)))
  }

  n <- ncol(X)
  d <- nrow(normalMean)
  k <- ncol(normalMean)




  logW <- array(0, dim=c(1,k))
  EQ <- array(0, dim=c(n,k))
  gammaBeta <- array(0, dim=c(n,k))

  for(i in 1:k) {

    U <- chol(invWhishartScale[,,i])

    logW[i] <- -2 * sum(log(diag(U)))

    Q <- solve(t(U), bsxfun.se("-", X, normalMean[,i]))

    QQ <- dot.ext(Q, Q, 1)

    EQ[,i] <- d/(normalRelativePrecision[i]* gammaDof[i]) + whishartDof[i]/gammaDof[i] * QQ	# eq (19)  bishop 10.64


    gammaBeta[,i] <- (d/(normalRelativePrecision[i]) + whishartDof[i] *QQ + gammaDof[i])/2
  }

  vd <- bsxfun.se("-", matrix(rep(whishartDof+1, d),nrow=d,byrow=T), as.matrix(1:d))/2

  ElogLambda <- colSums(digamma(vd)) + d*log(2) + logW	# 10.65
  Elogpi <- digamma(dirichletConcentration) - digamma(sum(dirichletConcentration))			# 10.66

  dgammaDof <- d+gammaDof

  logRho1 <- bsxfun.se("*", log(EQ+1), -dgammaDof)


  if (no.weights) {
    logRho <- 0.5 * bsxfun.se("+", logRho1,
                              + ElogLambda - d*log(gammaDof * pi) +
                                2 *lgamma((d+gammaDof)/2) - 2*lgamma(gammaDof/2)) #eq (19)  10.46
  } else {

    logRho <- 0.5 * bsxfun.se("+", logRho1,
                              2*Elogpi + ElogLambda - d*log(gammaDof * pi) +
                                2 *lgamma((d+gammaDof)/2) - 2*lgamma(gammaDof/2)) #eq (19)  10.46


  }

  # ke: add bound to avoid numerical issue
  #mapper <- logRho < -700
  #logRho <- mapper * -700 + (!mapper) * logRho
  if (n==k) {
    logR <- bsxfun.se("-", logRho, logsumexp(logRho, 1), expandByRow = F)	# 10.49
  } else {
    logR <- bsxfun.se("-", logRho, logsumexp(logRho, 1))	# 10.49
  }

  R <- exp(logR)

  model$logResponsibility <- logR
  model$responsibility<- R

  gammaAlpha <- dgammaDof/2
  hiddenScale <- bsxfun.se("*", 1/gammaBeta, gammaAlpha)

  model$gammaAlpha <- gammaAlpha
  model$gammaBeta <- gammaBeta
  model$hiddenScale <- hiddenScale

  model
}



############ Variational-(lower)-Bound Evaluation ############

variationalLowerBound_smm <- function(X, model, prior) {

  dirichletConcentration0 <- prior$dirichletConcentration
  normalRelativePrecision0 <- prior$normalRelativePrecision
  m0 <- prior$normalMean
  v0 <- prior$whishartDof
  M0 <- prior$invWhishartScale

  dirichletConcentration <- model$dirichletConcentration	# Dirichlet
  normalRelativePrecision <- model$normalRelativePrecision	# Gaussian
  m <- model$normalMean			# Gasusian
  v <- model$whishartDof			# Whishart
  M <- model$invWhishartScale			# Whishart: inv(W) = V'*V
  R <- model$responsibility
  logR <- model$logResponsibility
  hiddenScale <- model$hiddenScale

  d <- nrow(m)
  k <- ncol(m)
  nk <- colSums(R)									# 10.51
  RU <- R*hiddenScale
  omegak <- colSums(RU)

  Elogpi <- digamma(dirichletConcentration) - digamma(sum(dirichletConcentration))		# 10.66

  Epz = pracma::dot(nk, Elogpi)								# 10.72
  Eqz = pracma::dot(as.numeric(R), as.numeric(logR))			# 10.75
  # logCalpha0 = lgamma(k * alpha0) - k * lgamma(alpha0) # for scalar alpha0
  logCdirichletConcentration0 = lgamma(sum(dirichletConcentration0)) - sum(lgamma(dirichletConcentration0))
  # Eppi <- logCalpha0+(alpha0-1)*sum(Elogpi) # for scalar alpha0
  Eppi <- logCdirichletConcentration0+pracma::dot(dirichletConcentration0-1, Elogpi)			# 10.73
  logCdirichletConcentration <- lgamma(sum(dirichletConcentration))-sum(lgamma(dirichletConcentration))
  Eqpi = pracma::dot(dirichletConcentration-1, Elogpi) + logCdirichletConcentration				# 10.76

  # part of 10.70
  L <- Epz - Eqz + Eppi - Eqpi

  U0 <- chol(M0)
  sqrtR <- sqrt(R)
  sqrtRU <- sqrt(RU)
  xbar <- bsxfun.se("*", X %*% (R*hiddenScale), 1/omegak)				# 10.52

  logW <- array(0, dim = c(1, k))
  trSW <- array(0, dim = c(1, k))
  trM0W <- array(0, dim = c(1, k))
  xbarmWxbarm <- array(0, dim = c(1, k))
  mm0Wmm0 <- array(0, dim = c(1, k))

  for(i in 1:k) {

    U <- chol(M[,,i])
    logW[i] <- -2 * sum(log(diag(U)))

    Xs <- bsxfun.se("*", bsxfun.se("-", X, as.matrix(xbar[,i,drop=F])), t(sqrtRU[,i,drop=F]))
    V <- chol(tcrossprod(Xs, Xs)/omegak[i])
    Q <- solve(U, V)
    # equivalent to tr(SW)=trace(S/M)
    trSW[i] <- pracma::dot(as.numeric(Q), as.numeric(Q))
    Q <- solve(U, U0)
    trM0W[i] <- pracma::dot(as.numeric(Q), as.numeric(Q))

    q <- solve(t(U), xbar[,i,drop=F]-m[,i,drop=F])
    xbarmWxbarm[i] = pracma::dot(q, q)
    #q <- solve(t(U), m[,i,drop=F]-m0)
    q <- solve(t(U), m[,i,drop=F]-m0[,i,drop=F]) #ke: allow multivariate prior
    mm0Wmm0[i] <- pracma::dot(q, q)
  }

  vd <- bsxfun.se("-", matrix(rep(v+1, d),nrow=d,byrow=T), as.matrix(1:d))/2
  ElogLambda <- colSums(digamma(vd)) + d*log(2) + logW	# 10.65

  # first half of 10.74
  Epmu <- sum(d*log(normalRelativePrecision0/(2*pi))+ElogLambda-d*normalRelativePrecision0/normalRelativePrecision-normalRelativePrecision0*(v*mm0Wmm0))/2
  logB0 <- v0*sum(log(diag(U0)))-0.5*v0*d*log(2)-logmvgamma(0.5*v0,d) # B.79
  # second half of 10.74
  EpLambda <- k*logB0+0.5*(v0-d-1)*sum(ElogLambda)-0.5*pracma::dot(v,trM0W)


  Eqmu <- 0.5*sum(ElogLambda+d*log(normalRelativePrecision/(2*pi)))-0.5*d*k	# 10.77 (1/2)
  logB <- (-v) * (logW+d*log(2))/2 - logmvgamma(0.5*v, d)	# B.79
  EqLambda <- 0.5*sum((v-d-1)*ElogLambda-v*d)+sum(logB) 	# 10.77 (2/2)

  gammaAlpha <- model$gammaAlpha
  gammaBeta <- model$gammaBeta
  ElogU <- bsxfun.se("+", -log(gammaBeta), digamma(gammaAlpha))
  Uterm1 <- colSums(R * ElogU) / nk
  Uterm2 <- omegak/nk

  EpX <- 0.5*pracma::dot(nk, ElogLambda +
                   d*Uterm1 - d*Uterm2/normalRelativePrecision - v*Uterm2*trSW -v*xbarmWxbarm/nk-d*log(2*pi))	# eq (41) 10.71

  gammaDof <- model$gammaDof

  gammaDof2 <- gammaDof/2

  term1 <- bsxfun.se("*", ElogU, gammaDof2 - 1) - bsxfun.se("*", hiddenScale, gammaDof2)

  term2 <- bsxfun.se("+", term1, gammaDof2 * log(gammaDof2) - lgamma(gammaDof2) )

  Epu <- pracma::dot(as.numeric(R), as.numeric(term2))

  term1 <- -lgamma(gammaAlpha) + (gammaAlpha - 1) * digamma(gammaAlpha) - gammaAlpha
  term2 <- bsxfun.se("+", log(gammaBeta), term1)
  Equ <- pracma::dot(as.numeric(R), as.numeric(term2))

  L <- L+Epmu-Eqmu+EpLambda-EqLambda + EpX + Epu - Equ	# 10.70

  L
}









