############ Sort the model  ############
# Sort model paramters in increasing order of averaged means
# of d variables
sort_components <- function(model) {

  idx <- order(apply(model$ccfMean, 2, mean))

  model$ccfMean <- model$ccfMean[, idx]

  model$ccfCov <- model$ccfCov[,idx]

  model$responsibility <- model$responsibility[, idx]

  model$logResponsibility <- model$logResponsibility[, idx]

  model$dirichletConcentration <- model$dirichletConcentration[idx]

  model
}


#' Run Ccube with model 6: Normal-Binomial
#' @param mydata mutation data frame
#' @param epi sequencing error
#' @param init scalar number of clusters, or vector all possible cluster numbers
#' @param prior prior parameters
#' @param tol VBEM convergence threshold
#' @param maxiter VBEM maximum iteration
#' @param fit_mult whether or not to estimate multiplicities
#' @param fit_hyper whether or not to estimat hyperparameters
#' @param use methods to get rough estimates of ccf
#' @param verbose show progress
#' @export
#' @return a list containing model parameters
ccube_m6 <- function(mydata, epi=1e-3, init=2, prior, tol=1e-20, maxiter=1e3, fit_mult = F, fit_hyper = T, use = c("use_base", "use_one"),verbose=FALSE) {

  stopifnot(
    all(c("var_counts","ref_counts","normal_cn",
          "major_cn","minor_cn","purity") %in% names(mydata)))
  mydata <- GetCcf(mydata, use=use)

  dn <- mydata$ref_counts + mydata$var_counts
  bn <- mydata$var_counts
  cn <- unique(mydata$normal_cn)
  cr <- mydata$major_cn + mydata$minor_cn
  major_cn <- mydata$major_cn
  purity <- unique(mydata$purity)
  bv <- mydata$mult
  rawCcf <- mydata$ccf
  rawCcf <- as.matrix(rawCcf)
	n <- nrow(rawCcf)
	d <- ncol(rawCcf)

	X <- t(rawCcf) # Work with D by N for convenience

	message(sprintf("Running VB-SMM on a %d-by-%d data with %d clusters ...\n", n, d, init))

	if(missing(prior)) {

		if( length(init) == 1 ) {
			# General prior with equal dirichletConcentration
		  k <-  init
		  prior <- list(
		    dirichletConcentration = 1e-2,
		    normalMean = 1,
		    invWhishartScale = var(rawCcf)*(d+1)# M = inv(W)
			)
		}
	} else {

		stopifnot(
			all(names(prior) %in% c("dirichletConcentration","normalMean","invWhishartScale")) &
			all(sapply(prior, is.numeric)) & nrow(prior$normalMean) == d &
			ncol(prior$normalMean) == 1 &
			nrow(prior$normalMean) == d & ncol(prior$invWhishartScale) == d)
	}

	# variational lower bound (objective function)
	L <- rep(-Inf, maxiter)
	converged <- FALSE
	degenerated <- FALSE
	vbiter <- 1

	model <- list()

	initParams <- initialization(X, init, prior) # initialize responsibility and hidden scale
	model$responsibility <- initParams$R
	model$ccfMean <- initParams$ccfMean
	model$ccfCov <- initParams$ccfCov
	model$bv <- bv
	model$dirichletConcentration0 <- prior$dirichletConcentration
	model$normalMean <- prior$normalMean
	model$invWhishartScale <- prior$invWhishartScale

	while(!converged & vbiter < maxiter & !degenerated) {
	  vbiter <- vbiter + 1
		model <- VariationalMaximimizationStep(bn, dn, cn, cr, major_cn, epi, purity, model,
		                                       fit_mult = fit_mult, fit_hyper = fit_hyper)
		model <- VarationalExpectationStep(bn, dn, cn, cr, epi, purity, model)
		L[vbiter] <- variationalLowerBound(bn, dn, cn, cr, epi, purity, model)/n
		converged <- abs(L[vbiter] - L[vbiter-1]) < (tol) * abs(L[vbiter])
		degenerated <- (L[vbiter] - L[vbiter-1]) < 0
		#degenerated = F
		if(verbose) cat(sprintf("\rVB-EM-%d: L = %.8f \r", vbiter, L[vbiter]))
	}

	L <- L[2:vbiter]

	model <- sort_components(model)

	if (init > 1) {
	  label <- rep(0, n)
	  label <- apply(model$responsibility, 1, which.max)
	  nk <- colSums(model$responsibility)
	  Epi <- (model$dirichletConcentration + nk) / (k*model$dirichletConcentration0 + n)
	  model$Epi <- Epi/sum(Epi)
	} else {
	  label <- rep(1, n)
	  model$Epi <- 1
	}


	if(converged) cat(sprintf("\nConverged in %d steps.\n", vbiter-1)) else
	  cat(sprintf("Not converged in %d steps.\n", maxiter))
	if(degenerated)
	  cat(sprintf("Degenerated at %d steps.\n", vbiter-1))

	list(label=label, R=model$responsibility, mu=model$ccfMean, full.model=model, L=L)
}

#' Get a rough estimate of ccf and multiplicities
#' @param mydata mutation data frame
#' @param use methods for get rough estimates of ccf
#' @return mydata mutation data frame
#' @export
GetCcf <- function(mydata, use = c("use_base", "use_one")) {
  GetMult <- function(x, y, z) {
    index <- which(x==z)
    mean(y[index])
  }
  if (!"total_counts" %in% names(mydata)) {
    mydata <- dplyr::mutate(mydata, total_counts = ref_counts + var_counts)
  }

  if (!"total_cn" %in% names(mydata)) {
    mydata <- dplyr::mutate(mydata, total_cn = major_cn + minor_cn)
  }

  if (use=="use_base") {
    mydata <- dplyr::mutate(dplyr::rowwise(mydata),
                            ccf1 = MapVaf2CcfPyClone(var_counts/total_counts,
                                                     purity,
                                                     normal_cn,
                                                     total_cn,
                                                     total_cn,
                                                     major_cn, lower = 0, upper = 2),
                            ccf2 = MapVaf2CcfPyClone(var_counts/total_counts,
                                                     purity,
                                                     normal_cn,
                                                     total_cn,
                                                     total_cn,
                                                     minor_cn, lower = 0, upper = 2),
                            ccf3 = MapVaf2CcfPyClone(var_counts/total_counts,
                                                     purity,
                                                     normal_cn,
                                                     total_cn,
                                                     total_cn,
                                                     1, lower = 0, upper = 2))

    mydata <- dplyr::mutate(dplyr::rowwise(mydata), ccf = mean(c(ccf1, ccf2, ccf3), na.rm = T))
    dd <- dplyr::filter(mydata, ccf1 != ccf2 | ccf1 != ccf3 | ccf2 != ccf3)
    dd <- dplyr::mutate(dplyr::rowwise(dd), ccf = min(c(ccf1, ccf2, ccf3), na.rm = T ))
    mydata[mydata$mutation_id %in% dd$mutation_id,]$ccf = dd$ccf

    mydata <- dplyr::mutate(dplyr::rowwise(mydata),
                            mult = GetMult(c(ccf1, ccf2, ccf3),
                                           c(major_cn, minor_cn, 1), ccf))

    dd1 <- dplyr::filter(mydata, is.na(ccf))
    if (nrow(dd1) > 0) {
      dd1 <- dplyr::mutate(dplyr::rowwise(dd1),
                           ccf = MapVaf2CcfPyClone(var_counts/total_counts,
                                                   purity,
                                                   normal_cn,
                                                   total_cn,
                                                   total_cn,
                                                   1, constraint = F))


      mydata[mydata$mutation_id %in% dd1$mutation_id,]$ccf = dd1$ccf
      mydata[mydata$mutation_id %in% dd1$mutation_id,]$mult = 1
    }
  }

  if (use == "use_one") {
    mydata$mult = 1
    mydata <- dplyr::mutate(dplyr::rowwise(mydata),
                         ccf = MapVaf2CcfPyClone(var_counts/total_counts,
                                                 purity,
                                                 normal_cn,
                                                 total_cn,
                                                 total_cn,
                                                 1, constraint = F))
  }


 return(mydata)
}


############ helper function to initialize responsibility and other parameters ############
initialization <- function(X, init, prior) {

  d <- nrow(X)

  n <- ncol(X)

  stopifnot(length(init) %in% c(1, n) ||
              (nrow(init) == d  & ncol(init) == k))

  k <- init
  res <- kmeans(t(X), k)
  label <- res$cluster
  normalMean <- t(res$centers)
  if (length(which(normalMean > 1) ) >0 ){
    normalMean[which(normalMean>1)] = 1
  }
  R <- as.matrix(Matrix::sparseMatrix(1:n, label, x=1))
  invWhishartScale <- prior$invWhishartScale
  ccfMean = normalMean
  ccfCov = array(invWhishartScale, dim = c(d,k))

  return(list(R = R,
              ccfMean = ccfMean, ccfCov = ccfCov))
}


############ Variational-Maximimization ############
VariationalMaximimizationStep <- function(bn, dn, cn, cr, major_cn, epi, purity, model, fit_mult = F, fit_hyper = T) {

  bv <- model$bv
  dirichletConcentration0 <- model$dirichletConcentration0
  responsibility <- model$responsibility
  normalMean <- model$normalMean
  invWhishartScale <- model$invWhishartScale



  # Gaussian approximation for each clone
  ccfMean = model$ccfMean
  ccfCov = model$ccfCov

  k <- length(ccfMean)
    Bn = (1-purity)*cn + purity*cr
    Cn = purity*(bv*(1-epi) - cr*epi)
    for (i in 1:k){
      term1 = 1/invWhishartScale*normalMean
      term2 = 1/invWhishartScale
      upper <- 1
      lower <- 0
      tmp <- NULL
      jj <- 0
      while (!is.numeric(tmp)) {
        if (jj >=  1000) {
          break
        }
        jj <- jj + 1
        if (jj > 1) {upper <- upper + 0.1}
        tmp <- try(suppressWarnings( uniroot(
          function(x) {
            An = (1-purity)*cn*epi + purity*(1-x)*cr*epi + purity*x*bv*(1-epi)
            term4 = (bn/An - (dn-bn)/(Bn-An)) * Cn
            term5 = term2*x
            return(sum(responsibility[,i]*term4) - term5 + term1)
          },
          c(lower, upper), extendInt = "yes")$root), T)
      }
      ccfMean[i] <- tmp

      An = (1-purity)*cn*epi + purity*(1-ccfMean[i])*cr*epi + purity*ccfMean[i]*bv*(1-epi)
      ccfCov[i] = solve(term2 + sum(responsibility[,i]*(bn/(An^2) + (dn-bn)/((Bn-An)^2) * Cn^2)))
    }


  dataSumResponsibility <- colSums(responsibility) # 10.51
  dirichletConcentration <- dirichletConcentration0 + dataSumResponsibility # 10.58


  if (length(which(ccfMean > 1) ) >0 ){
    ccfMean[which(ccfMean>1)] = 1
  }
  if (length(which(ccfMean < 1e-20)) > 0){
    ccfMean[which(ccfMean<1e-20)] = 1e-20
  }


  numberOfDataPoints <- nrow(responsibility)
  if (fit_mult) {
    for (ii in 1:numberOfDataPoints) {
      qq <- rep(0, major_cn[ii])
      bvPool <- 0:major_cn[ii]
      for (jj in seq_along(bvPool) ) {
        aa <- purity * (bvPool[jj] *(1-epi) -cr[ii]*epi) / ((1-purity)*cn + purity * cr[ii])
        aa2 <- aa^2
        bb <- epi
        term1 <- sum(responsibility[ii, ] * bn[ii] * (log (aa * ccfMean +bb) - aa2*ccfCov/(2 * (aa * ccfMean +bb)^2 ) ))
        term2 <- sum(responsibility[ii, ] * (dn[ii] - bn[ii]) * (log (1 - aa * ccfMean - bb) - aa2*ccfCov/(2 * (1 - aa * ccfMean -bb)^2)  ))
        term3 <- sum( responsibility[ii, ]*log(choose(dn[ii], bn[ii]) ))
        qq[jj] <- term1 + term2 + term3
      }
      bv[ii] <- bvPool[which.max(qq)]

    }
    model$bv <- bv
  }

  # estimate hyper-parameters
  if (fit_hyper) {
#     Elogpi <- sum(digamma(dirichletConcentration) - digamma(sum(dirichletConcentration)))
#     tmp <- NULL
#     jj <- 0
#     upper <- 1e-2
#     lower <- 1e-99
#     while (!is.numeric(tmp)) {
#       stopifnot(jj < 1000)
#       jj <- jj + 1
#       if (jj > 1) {upper <- upper + 1e-3}
#       tmp <- try(suppressWarnings( uniroot(
#         function(x) {
#           term1 <- digamma(k*x) - digamma(x)
#           return(term1+Elogpi/k)
#         },
#         c(lower, upper), extendInt = "yes")$root), T)
#     }
#     model$dirichletConcentration0 <- tmp
    model$normalMean <- mean(ccfMean)
    model$invWhishartScale <-mean((ccfMean - normalMean)^2 + ccfCov)
  }

  model$dirichletConcentration <- dirichletConcentration
  model$ccfMean <- ccfMean
  model$ccfCov <- ccfCov

  model
}


############ Variational-Expectation ############
VarationalExpectationStep <- function(bn, dn, cn, cr, epi, purity, model, no.weights = FALSE) {

  bv <- model$bv
  dirichletConcentration <- model$dirichletConcentration	# Dirichlet
  ccfMean <- model$ccfMean
  ccfCov <- model$ccfCov
  Epi <- model$Epi

  n <- length(bn)
  d <- nrow(ccfMean)
  k <- ncol(ccfMean)


  Epbk <- array(0, dim=c(n,k))

  aa <- purity * (bv *(1-epi) -cr*epi) / ((1-purity)*cn + purity * cr)
  aa2 <- aa^2
  bb <- epi

  for(i in 1:k) {
    Epbk[,i] <- bn * (log (aa * ccfMean[i] +bb) - aa2*ccfCov[,i]/(2 * (aa * ccfMean[i] +bb)^2 ) ) +
      (dn - bn) * (log (1 - aa * ccfMean[i] - bb) - aa2*ccfCov[,i]/(2 * (1 - aa * ccfMean[i] -bb)^2)  )
  }
  Elogpi <- digamma(dirichletConcentration) - digamma(sum(dirichletConcentration))
  if (no.weights) {
    logRho <-  Epbk  #eq (19)  10.46
  } else {

    logRho <- bsxfun.se("+", Epbk, Elogpi)
  }

  if (n==k) {
    logR <- bsxfun.se("-", logRho, logsumexp(logRho, 1), expandByRow = F)	# 10.49
  } else {
    logR <- bsxfun.se("-", logRho, logsumexp(logRho, 1))	# 10.49
  }

  R <- exp(logR)

  model$logResponsibility <- logR
  model$responsibility<- R

  model
}



############ Variational-(lower)-Bound Evaluation ############
variationalLowerBound <- function(bn, dn, cn, cr, epi, purity, model) {

  bv <- model$bv
  dirichletConcentration0 <- model$dirichletConcentration0
  m0 <- model$normalMean
  M0 <- model$invWhishartScale

  dirichletConcentration <- model$dirichletConcentration	# Dirichlet
  m <- model$ccfMean			# Gasusian
  M <- model$ccfCov
  R <- model$responsibility
  logR <- model$logResponsibility
  ccfMean <- model$ccfMean
  ccfCov <- model$ccfCov
  nk <- colSums(R)									# 10.51

  n <- length(bn)
  d <- nrow(m)
  k <- ncol(m)

  Elogpi <- digamma(dirichletConcentration) - digamma(sum(dirichletConcentration))		# 10.66

  Epz = pracma::dot(nk, Elogpi)								# 10.72
  Eqz = pracma::dot(as.numeric(R), as.numeric(logR))			# 10.75
  logCdirichletConcentration0 = lgamma(k*dirichletConcentration0) - k*lgamma(dirichletConcentration0)
  Eppi <- logCdirichletConcentration0+(dirichletConcentration0-1)*sum(Elogpi)
  logCdirichletConcentration <- lgamma(sum(dirichletConcentration))-sum(lgamma(dirichletConcentration))
  Eqpi = pracma::dot(dirichletConcentration-1, Elogpi) + logCdirichletConcentration				# 10.76

  # part of 10.70
  L <- Epz - Eqz + Eppi - Eqpi

  U0 <- chol(M0)
  trM0W <- array(0, dim = c(1, k))
  mm0Wmm0 <- array(0, dim = c(1, k))

  Epbk <- array(0, dim=c(1,k))

  aa <- purity * (bv *(1-epi) -cr*epi) / ((1-purity)*cn + purity * cr)
  aa2 <- aa^2
  bb <- epi

  for(i in 1:k) {

    Epbk[,i] <- sum(R[,i]*(bn * (log (aa * ccfMean[i] +bb) - aa2*ccfCov[,i]/(2 * (aa * ccfMean[i] +bb)^2 ) ) +
      (dn - bn) * (log (1 - aa * ccfMean[i] - bb) - aa2*ccfCov[,i]/(2 * (1 - aa * ccfMean[i] -bb)^2)))) +
      sum(R[,i]*log(choose(dn,bn)))
    q <- solve(t(U0), m[,i,drop=F]-m0)
    mm0Wmm0[i] <- pracma::dot(q, q)
    U <- chol(M[,i])
    Q <- solve(U0, U)
    trM0W[i] <- pracma::dot(as.numeric(Q), as.numeric(Q))
  }

  # first half of 10.74
  Epmu <- sum(-log(2*pi)-log(M)-mm0Wmm0-trM0W)/2

  Eqmu <- -0.5*sum(1+log(2*pi)+log(M))	# 10.77 (1/2)

  L <- L + Epmu - Eqmu + sum(Epbk)	# 10.70

  L
}

cull <- function(model) {
  ccfMean <- model$ccfMean
  ccfCov <- model$ccfCov
  Epi <- model$Epi
  k <- ncol(model$ccfMean)
  kl <- rep(0, k)
  for (i in 1:k-1) {
    kl[i] <- monomvn::kl.norm(mu1=ccfMean[,i], S1 = ccfCov[,i], mu2 = ccfMean[,i+1], S2 = ccfCov[,i+1])
  }
  tol <- 10
  cullCluster <- kl < tol
}


#' Estimate purity
#' @param data
#' @return purity
#' @export
GetPurity <- function(mydata) {
  vtox<-function(v,nA,nB,tA,tB)
  {
    (nB - nB*v - nA*v) / (tA*v + tB*v + nB - tB -nA*v - nB*v)
  }
  tmpdata <- dplyr::filter(mydata, major_cn == minor_cn & major_cn != 0 )
  if (nrow(tmpdata) == 0) {
    return(NA)
  }
  tmpdata <- dplyr::mutate(dplyr::rowwise(tmpdata), ploidy = major_cn + minor_cn)
  tmpdata <- dplyr::mutate(tmpdata, vaf = var_counts/(var_counts+ref_counts))
  tmpdata <- dplyr::mutate(dplyr::rowwise(tmpdata), cp = vtox(vaf, 2,0,ploidy/2,ploidy/2))
  tmpdata <- dplyr::filter(tmpdata, !is.infinite(cp) & !is.na(cp))
  K = 6
  if (K >= nrow(tmpdata)) {
    K = nrow(tmpdata) - 1
    if (K == 0) {K = 1}
  }
  if (nrow(tmpdata) == 1) {
    purity <- tmpdata$cp
  } else {
    res <- vbsmm(tmpdata$cp, init = K, tol = 1e-5,  verbose = F)
    pool <- res$mu[unique(res$label)]
    ww <- (res$full.model$Epi[unique(res$label)]>1.5e-2)
    pool1 <- pool[ww]
    maxCp <- max(pool1[pool1 <=1])
    purity  <- if (maxCp > 1) 1 else maxCp
  }
  
  return(purity)
}





