#' Run Ccube Core function for SV
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
CcubeSVCore <- function(mydata, epi=1e-3, init=2, prior, tol=1e-20, maxiter=1e3, fit_mult = F, fit_hyper = T, use = c("use_base", "use_one"),verbose=T) {

  stopifnot(
    all(c("var_counts1","ref_counts1","normal_cn",
          "purity", "var_counts2","ref_counts2","subclonal_cn1", "subclonal_cn2") %in% names(mydata)))

  stopifnot(
    all(c("major_cn1_sub1","major_cn1_sub2","minor_cn1_sub1", "minor_cn1_sub2",
          "frac_cn1_sub1", "frac_cn1_sub2",
          "major_cn2_sub1","major_cn2_sub2","minor_cn2_sub1", "minor_cn2_sub2",
          "frac_cn2_sub1", "frac_cn2_sub2")
        %in% names(mydata)))

  mydata <- GetCcf_sv(mydata, use=use)

  dn1 <- mydata$ref_counts1 + mydata$var_counts1
  bn1 <- mydata$var_counts1
  cn <- unique(mydata$normal_cn)
  cr1 <- mydata$total_cn1
  max_mult_cn1_sub1 <- mydata$total_cn1_sub1
  max_mult_cn1_sub2 <- mydata$total_cn1_sub2
  frac_cn1_sub1 <- mydata$frac_cn1_sub1
  frac_cn1_sub2 <- mydata$frac_cn1_sub2
  bv1 <- mydata$mult1
  bv1_sub1 <- rep(-100, length(max_mult_cn1_sub1))
  bv1_sub2 <- rep(-100, length(max_mult_cn1_sub2))

  dn2 <- mydata$ref_counts2 + mydata$var_counts2
  bn2 <- mydata$var_counts2
  cr2 <- mydata$total_cn2
  max_mult_cn2_sub1 <- mydata$total_cn2_sub1
  max_mult_cn2_sub2 <- mydata$total_cn2_sub2
  frac_cn2_sub1 <- mydata$frac_cn2_sub1
  frac_cn2_sub2 <- mydata$frac_cn2_sub2
  bv2 <- mydata$mult2
  bv2_sub1 <- rep(-100, length(max_mult_cn2_sub1))
  bv2_sub2 <- rep(-100, length(max_mult_cn2_sub2))

  subclonal_cn1 <- mydata$subclonal_cn1
  subclonal_cn2 <- mydata$subclonal_cn2

  purity <- unique(mydata$purity)
  rawCcf <- mydata$ccf
  rawCcf <- as.matrix(rawCcf)
  n <- nrow(rawCcf)
  d <- ncol(rawCcf)

  X <- t(rawCcf) # Work with D by N for convenience

  message(sprintf("Running VB-VB-Normal-Binomial on a %d-by-%d data with %d clusters ...\n", n, d, init))

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
  model$bv1 <- bv1
  model$bv2 <- bv2
  model$bv1_sub1 <- bv1_sub1
  model$bv1_sub2 <- bv1_sub2
  model$bv2_sub1 <- bv2_sub1
  model$bv2_sub2 <- bv2_sub2
  model$dirichletConcentration0 <- prior$dirichletConcentration
  model$normalMean <- prior$normalMean
  model$invWhishartScale <- prior$invWhishartScale

  while(!converged & vbiter < maxiter & !degenerated) {
    vbiter <- vbiter + 1
    model <- VariationalMaximimizationStep_sv(bn1, dn1, cn, cr1, max_mult_cn1_sub1, max_mult_cn1_sub2,
                                              frac_cn1_sub1, frac_cn1_sub2, subclonal_cn1,
                                              bn2, dn2, cr2, max_mult_cn2_sub1, max_mult_cn2_sub2,
                                              frac_cn2_sub1, frac_cn2_sub2, subclonal_cn2,
                                              epi, purity, model,
                                              fit_mult = fit_mult, fit_hyper = fit_hyper)
    model <- VarationalExpectationStep_sv(bn1, dn1, cn, cr1,
                                                   bn2, dn2, cr2,
                                                   epi, purity, model)
    L[vbiter] <- VariationalLowerBound_sv(bn1, dn1, cn, cr1,
                                          bn2, dn2, cr2,
                                          epi, purity, model)/n
    converged <- abs(L[vbiter] - L[vbiter-1]) < (tol) * abs(L[vbiter])
    degenerated <- (L[vbiter] - L[vbiter-1]) < 0
    #degenerated = F
    if(verbose) cat(sprintf("\rVB-EM-%d: L = %.8f \r", vbiter, L[vbiter]))
  }

  L <- L[2:vbiter]

  model <- SortClusters(model)

  if (init > 1) {
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

  list(label=label, full.model=model, L=L)
}



#' Get a rough estimate of ccf and multiplicities
#' @param mydata mutation data frame
#' @param use methods for get rough estimates of ccf
#' @return mydata mutation data frame
#' @export
GetCcf_sv <- function(mydata, use = c("use_base", "use_one")) {
  GetMult <- function(x, y, z) {
    index <- which(x==z)
    mean(y[index])
  }
  if (!"total_counts1" %in% names(mydata)) {
    mydata <- dplyr::mutate(mydata, total_counts1 = ref_counts1 + var_counts1)
  }

  if (!"total_cn1" %in% names(mydata)) {
    mydata <- dplyr::mutate(mydata, total_cn1 = frac_cn1_sub1 * (major_cn1_sub1 + minor_cn1_sub1) +
                              frac_cn1_sub2 * (major_cn1_sub2 + minor_cn1_sub2) )
  }

  if (!"major_cn1" %in% names(mydata)) {
    mydata <- dplyr::mutate(mydata, major_cn1 = frac_cn1_sub1 * major_cn1_sub1 +
                              frac_cn1_sub2 * major_cn1_sub2)
  }

  if (!"minor_cn1" %in% names(mydata)) {
    mydata <- dplyr::mutate(mydata, minor_cn1 = frac_cn1_sub1 * minor_cn1_sub1 +
                              frac_cn1_sub2 * minor_cn1_sub2 )
  }

  if (!"total_counts2" %in% names(mydata)) {
    mydata <- dplyr::mutate(mydata, total_counts1 = ref_counts2 + var_counts2)
  }

  if (!"total_cn2" %in% names(mydata)) {
    mydata <- dplyr::mutate(mydata, total_cn2 = frac_cn2_sub1 * (major_cn2_sub1 + minor_cn2_sub1) +
                              frac_cn2_sub2 * (major_cn2_sub2 + minor_cn2_sub2))
  }

  if (!"major_cn2" %in% names(mydata)) {
    mydata <- dplyr::mutate(mydata, major_cn2 = frac_cn2_sub1 * major_cn2_sub1 +
                              frac_cn2_sub2 * major_cn2_sub2)
  }

  if (!"minor_cn2" %in% names(mydata)) {
    mydata <- dplyr::mutate(mydata, minor_cn2 = frac_cn2_sub1 * minor_cn2_sub1 +
                              frac_cn2_sub2 * minor_cn2_sub2 )
  }
  if (use=="use_base") {
    mydata <- dplyr::mutate(dplyr::rowwise(mydata),
                            ccf1_1 = MapVaf2CcfPyClone(var_counts1/total_counts1,
                                                     purity,
                                                     normal_cn,
                                                     total_cn1,
                                                     total_cn1,
                                                     major_cn1, lower = 0, upper = 2),
                            ccf1_2 = MapVaf2CcfPyClone(var_counts1/total_counts1,
                                                     purity,
                                                     normal_cn,
                                                     total_cn1,
                                                     total_cn1,
                                                     minor_cn1, lower = 0, upper = 2),
                            ccf1_3 = MapVaf2CcfPyClone(var_counts1/total_counts1,
                                                     purity,
                                                     normal_cn,
                                                     total_cn1,
                                                     total_cn1,
                                                     1, lower = 0, upper = 2),
                            ccf2_1 = MapVaf2CcfPyClone(var_counts2/total_counts2,
                                                       purity,
                                                       normal_cn,
                                                       total_cn2,
                                                       total_cn2,
                                                       major_cn2, lower = 0, upper = 2),
                            ccf2_2 = MapVaf2CcfPyClone(var_counts2/total_counts2,
                                                       purity,
                                                       normal_cn,
                                                       total_cn2,
                                                       total_cn2,
                                                       minor_cn2, lower = 0, upper = 2),
                            ccf2_3 = MapVaf2CcfPyClone(var_counts2/total_counts2,
                                                       purity,
                                                       normal_cn,
                                                       total_cn2,
                                                       total_cn2,
                                                       1, lower = 0, upper = 2)
                            )

    mydata <- dplyr::mutate(dplyr::rowwise(mydata), ccf1 = mean(c(ccf1_1, ccf1_2, ccf1_3), na.rm = T),
                            ccf2 = mean(c(ccf2_1, ccf2_2, ccf2_3), na.rm = T))

    dd <- dplyr::filter(mydata, ccf1_1 != ccf1_2 | ccf1_1 != ccf1_3 | ccf1_2 != ccf1_3)

    if (nrow(dd) > 0) {
      dd <- dplyr::mutate(dplyr::rowwise(dd), ccf1 = min(c(ccf1_1, ccf1_2, ccf1_3), na.rm = T ))
      mydata[mydata$mutation_id %in% dd$mutation_id,]$ccf1 = dd$ccf1
    }

    dd <- dplyr::filter(mydata, ccf2_1 != ccf2_2 | ccf2_1 != ccf2_3 | ccf2_2 != ccf2_3)

    if (nrow(dd) > 0) {
      dd <- dplyr::mutate(dplyr::rowwise(dd), ccf2 = min(c(ccf2_1, ccf2_2, ccf2_3), na.rm = T ))
      mydata[mydata$mutation_id %in% dd$mutation_id,]$ccf2 = dd$ccf2
    }


    mydata <- dplyr::mutate(dplyr::rowwise(mydata),
                            mult1 = GetMult(c(ccf1_1, ccf1_2, ccf1_3),
                                           c(major_cn1, minor_cn1, 1), ccf1),
                            mult2 = GetMult(c(ccf2_1, ccf2_2, ccf2_3),
                                            c(major_cn2, minor_cn2, 1), ccf2) )

    dd1 <- dplyr::filter(mydata, is.na(ccf1))
    if (nrow(dd1) > 0) {
      dd1 <- dplyr::mutate(dplyr::rowwise(dd1),
                           ccf1 = MapVaf2CcfPyClone(var_counts1/total_counts1,
                                                   purity,
                                                   normal_cn,
                                                   total_cn1,
                                                   total_cn1,
                                                   1, constraint = F))


      mydata[mydata$mutation_id %in% dd1$mutation_id,]$ccf1 = dd1$ccf1
      mydata[mydata$mutation_id %in% dd1$mutation_id,]$mult1 = 1
    }

    dd1 <- dplyr::filter(mydata, is.na(ccf2))
    if (nrow(dd1) > 0) {
      dd1 <- dplyr::mutate(dplyr::rowwise(dd1),
                           ccf2 = MapVaf2CcfPyClone(var_counts2/total_counts2,
                                                    purity,
                                                    normal_cn,
                                                    total_cn2,
                                                    total_cn2,
                                                    1, constraint = F))


      mydata[mydata$mutation_id %in% dd1$mutation_id,]$ccf2 = dd1$ccf2
      mydata[mydata$mutation_id %in% dd1$mutation_id,]$mult2 = 1
    }

  }

  if (use == "use_one") {
    mydata$mult1 = 1
    mydata$mult2 = 1
    mydata <- dplyr::mutate(dplyr::rowwise(mydata),
                         ccf1 = MapVaf2CcfPyClone(var_counts1/total_counts1,
                                                 purity,
                                                 normal_cn,
                                                 total_cn1,
                                                 total_cn1,
                                                 1, constraint = F),
                         ccf2 = MapVaf2CcfPyClone(var_counts2/total_counts2,
                                                  purity,
                                                  normal_cn,
                                                  total_cn2,
                                                  total_cn2,
                                                  1, constraint = F))
  }

 mydata <- dplyr::mutate(dplyr::rowwise(mydata), ccf = mean( ccf1, ccf2 ) )
 return(mydata)
}


neg_ELBO_bv <- function(x, responsibility, bn, dn, cr, cn, purity, ccfMean, ccfCov, epi=1e-3) {

  w <- purity * (x *(1-epi) -cr*epi) / ((1-purity)*cn + purity * cr)
  ww <- w^2
  ef <- w * ccfMean +epi

  term1 <- bn* (log (ef) - ww*ccfCov/(2 * ef^2 ) )
  term2 <- (dn- bn) * (log (1 - ef) - ww*ccfCov/(2 * (1 - ef)^2))

  L =  sum( responsibility * (term1 + term2 ) )

  return( - L)
}


neg_ELBO_bv_g <- function(x, responsibility, bn, dn, cr, cn, purity, ccfMean, ccfCov, epi=1e-3) {

  w <- purity * (x *(1-epi) -cr*epi) / ((1-purity)*cn + purity * cr)
  ww <- w^2
  ef <- w * ccfMean +epi
  term1 = bn
  term2 = dn-bn

  dLdEf = term1 * (1/ef + ww*ccfCov/ ef^3 ) + term2 * ( -1/(1-ef) - ww*ccfCov / (1-ef)^3 )
  dLdW = term1 * ( - w*ccfCov / ef^2  ) + term2 * ( - w*ccfCov/(1-ef)^2 )
  dEfdW = ccfMean

  dWdBv =  purity * (1-epi)  / ((1-purity)*cn + purity * cr)

  g = sum( responsibility* ( dLdEf*dEfdW*dWdBv + dLdW*dWdBv ) )

  g = -g
}

ELBO_bv_g <- function(x, responsibility, bn, dn, cr, cn, purity, ccfMean, ccfCov, epi=1e-3) {

  w <- purity * (x *(1-epi) -cr*epi) / ((1-purity)*cn + purity * cr)
  ww <- w^2
  ef <- w * ccfMean +epi
  term1 = bn
  term2 = dn-bn

  dLdEf = term1 * (1/ef + ww*ccfCov/ ef^3 ) + term2 * ( -1/(1-ef) - ww*ccfCov / (1-ef)^3 )
  dLdW = term1 * ( - w*ccfCov / ef^2  ) + term2 * ( - w*ccfCov/(1-ef)^2 )
  dEfdW = ccfMean

  dWdBv =  purity * (1-epi)  / ((1-purity)*cn + purity * cr)

  g = sum( responsibility* ( dLdEf*dEfdW*dWdBv + dLdW*dWdBv ) )
}


neg_ELBO_ccf <- function(x, bn1, dn1, cr1, bv1, bn2, dn2, cr2, bv2,
                         responsibility, normalMean, invWhishartScale, cn, purity, epi=1e-3) {

  w1 = ( purity*(bv1*(1-epi) - cr1*epi) ) / ( (1-purity)*cn + purity*cr1 )
  w2 = ( purity*(bv2*(1-epi) - cr2*epi) ) / ( (1-purity)*cn + purity*cr2 )

  ef1 = w1 * x + epi
  ef2 = w2 * x + epi


  term1 <- bn1* log (ef1) + (dn1- bn1) * (log (1 - ef1) )
  term2 <- bn2* log (ef2) + (dn2- bn2) * (log (1 - ef2) )

  L =  sum( responsibility * (term1 + term2 ) ) -
    (x-normalMean)^2/(2 * invWhishartScale)

  return( - L)
}

neg_ELBO_ccf_g <- function(x, bn1, dn1, cr1, bv1, bn2, dn2, cr2, bv2,
                           responsibility, normalMean, invWhishartScale, cn, purity, epi=1e-3) {

  w1 = ( purity*(bv1*(1-epi) - cr1*epi) ) / ( (1-purity)*cn + purity*cr1 )
  w2 = ( purity*(bv2*(1-epi) - cr2*epi) ) / ( (1-purity)*cn + purity*cr2 )

  ef1 = w1 * x + epi
  ef2 = w2 * x + epi

  dLdEf1 = bn1/ef1  + (dn1- bn1) / (1 - ef1)
  dLdEf2 = bn2/ef2  + (dn2- bn2) / (1 - ef2)

  dLdCcf = (normalMean-x)  /invWhishartScale
  dEf1dCcf = w1
  dEf2dCcf = w2

  g = sum( responsibility *(  dLdEf1* dEf1dCcf + dLdEf2* dEf2dCcf) ) + dLdCcf

  g = -g
}


my_repmat <- function(x, n) {
  # this is much faster than pracma::repmat
  xMat = matrix(rep(as.numeric(x), each=n), nrow=n)
  return( xMat )
}


############ Variational-Maximimization ############
VariationalMaximimizationStep_sv <- function(bn1, dn1, cn, cr1, max_mult_cn1_sub1, max_mult_cn1_sub2,
                                             frac_cn1_sub1, frac_cn1_sub2, subclonal_cn1,
                                             bn2, dn2, cr2, max_mult_cn2_sub1, max_mult_cn2_sub2,
                                             frac_cn2_sub1, frac_cn2_sub2, subclonal_cn2,
                                             epi, purity, model,
                                             fit_mult = T, fit_hyper = T) {

  bv1 = model$bv1
  bv2 = model$bv2
  bv1_sub1 = model$bv1_sub1
  bv1_sub2 = model$bv1_sub2
  bv2_sub1 = model$bv2_sub1
  bv2_sub2 = model$bv2_sub2

  dirichletConcentration0 <- model$dirichletConcentration0
  responsibility <- model$responsibility
  normalMean <- model$normalMean
  invWhishartScale <- model$invWhishartScale

  # Gaussian approximation for each clone
  ccfMean = model$ccfMean
  ccfCov = model$ccfCov

  k <- length(ccfMean)

  w1 = ( purity*(bv1*(1-epi) - cr1*epi) ) / ( (1-purity)*cn + purity*cr1 )
  w2 = ( purity*(bv2*(1-epi) - cr2*epi) ) / ( (1-purity)*cn + purity*cr2 )

  ccfMeanOld <- ccfMean

    for (i in 1:k){


      term1 = 1/invWhishartScale*normalMean
      term2 = 1/invWhishartScale
      upper <- 1
      lower <- 0
      tmp <- NULL
      jj <- 0
      while (!is.numeric(tmp)) {
        if (jj >=  1000) {
          tmp <- ccfMeanOld[i]
          break
        }
        jj <- jj + 1

        if (jj > 1) {upper <- upper + 0.1}
        tmp <- try(suppressWarnings( uniroot(
          function(x) {

            ef1 = w1 * x + epi
            ef2 = w2 * x + epi

            term4 = (bn1/ef1 - (dn1-bn1)/(1 - ef1)) * w1 + (bn2/ef2 - (dn2-bn2)/(1- ef2)) * w2

            term5 = term2*x
            return(sum(responsibility[,i]*term4) - term5 + term1)
          },
          c(lower, upper), extendInt = "yes")$root), T)
      }

      ccfMean[i] <- tmp

      ef1 = w1 * ccfMean[i] + epi
      ef2 = w2 * ccfMean[i] + epi

      ccfCov[i] = solve(1/invWhishartScale + sum(responsibility[,i]*( w1^2 * (bn1/(ef1)^2 + (dn1-bn1)/( (1-ef1)^2 ) ) +
                                                           w2^2 * (bn2/(ef2)^2 + (dn2-bn2)/( (1-ef2)^2 ) ) )
                                    )
                        )
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
      # TODO: General implementation of more than two subclonal copy number populations

      # break point 1
      if (subclonal_cn1[ii]) {


        sub_cn1_mults = pracma::meshgrid(0:max_mult_cn1_sub1[ii], 0:max_mult_cn1_sub2[ii])
        bvPool1 <- frac_cn1_sub1[ii] * sub_cn1_mults$X + frac_cn1_sub2[ii] * sub_cn1_mults$Y
        #bvPool1[ bvPool1<1 ] = NA
        bvPool1Mat <- t(my_repmat(as.vector(bvPool1), length(ccfMean)))
        ccfMeanMat <- my_repmat(ccfMean, length(bvPool1))
        ccfCovMat <- my_repmat(ccfCov, length(bvPool1))
        respMat <- my_repmat(responsibility[ii, ], length(bvPool1))
        w1 <- purity * (bvPool1Mat*(1-epi) -cr1[ii]*epi) / ((1-purity)*cn + purity * cr1[ii])
        w1w1 <- w1^2
        ef1 <- w1 * ccfMeanMat +epi
        term1_breakpoint1 <- bn1[ii] * (log (ef1) - w1w1*ccfCovMat/(2 * ef1^2 ) )
        term2_breakpoint1 <- (dn1[ii] - bn1[ii]) * (log (1 - ef1) - w1w1*ccfCovMat/(2 * (1 - ef1)^2)  )
        term3_breakpoint1 <- logChoose(dn1[ii], bn1[ii])
        qq1 <- rowSums ( respMat *  (term1_breakpoint1 + term2_breakpoint1 + term3_breakpoint1))
        qq1[1] = NA # remove both multiplicities being zero
        maxQq1 <- which.max(qq1)
        bv1[ii] <- bvPool1[maxQq1]
        bv1_sub1[ii] <- sub_cn1_mults$X[maxQq1]
        bv1_sub2[ii] <- sub_cn1_mults$Y[maxQq1]


        # m_upper1 = (1-epi)*( (1-purity)*cn + purity*cr1[ii] ) / (ccfMean*purity)*(1-epi) +
        #   epi*cr1[ii]/(1-epi)
        # m_lower1 = (0-epi)*( (1-purity)*cn + purity*cr1[ii] ) / (ccfMean*purity)*(1-epi) +
        #   epi*cr1[ii]/(1-epi)
        # bv1_old = bv1[i]
        # upper <- min( c(m_upper1, major_cn1[ii]))
        # lower <- max(c(m_lower1,1))
        # if (lower >= upper & upper == major_cn1[ii]) {
        #   bv1[ii] = major_cn1[ii]
        # } else {
        #   tmp <- try(suppressWarnings( uniroot(
        #     ELBO_bv_g, c(lower, upper),
        #     bn = bn1[ii], dn=dn1[ii], cr=cr1[ii], cn=cn,
        #     purity=purity, ccfMean=ccfMean, ccfCov=ccfCov,
        #     epi=epi, responsibility = responsibility[ii,])$root), T)
        #
        #   if (!is.numeric(tmp)) {
        #     bv1[ii] <- min( bv1_old, upper)
        #   } else {
        #     bv1[ii] <- tmp
        #   }
        # }
        #

      } else {
        # clonal cn
        bvPool1 <- 1:max_mult_cn1_sub1[ii]

        qq1 <- rep(NA, length(bvPool1))
        for (jj in seq_along(bvPool1) ) {

          w1 <- purity * (bvPool1[jj] *(1-epi) -cr1[ii]*epi) / ((1-purity)*cn + purity * cr1[ii])
          w1w1 <- w1^2
          ef1 <- w1 * ccfMean +epi

          term1_breakpoint1 <- bn1[ii] * (log (ef1) - w1w1*ccfCov/(2 * ef1^2 ) )
          term2_breakpoint1 <- (dn1[ii] - bn1[ii]) * (log (1 - ef1) - w1w1*ccfCov/(2 * (1 - ef1)^2)  )
          term3_breakpoint1 <- logChoose(dn1[ii], bn1[ii])

          qq1[jj] <- sum ( responsibility[ii, ] *  (term1_breakpoint1 + term2_breakpoint1 + term3_breakpoint1)  )
        }
        bv1[ii] <- bvPool1[which.max(qq1)]

      }


      # break point 2
      if (subclonal_cn2[ii]) {
        sub_cn2_mults = pracma::meshgrid(0:max_mult_cn2_sub1[ii], 0:max_mult_cn2_sub2[ii])
        bvPool2 <- frac_cn2_sub1[ii] * sub_cn2_mults$X + frac_cn2_sub2[ii] * sub_cn2_mults$Y
        #bvPool2[ bvPool2<1 ] = NA
        bvPool2Mat <- t(my_repmat(as.vector(bvPool2), length(ccfMean)) )
        ccfMeanMat <- my_repmat(ccfMean, length(bvPool2))
        ccfCovMat <- my_repmat(ccfCov, length(bvPool2))
        respMat <- my_repmat(responsibility[ii, ], length(bvPool2))
        w2 <- purity * (bvPool2Mat *(1-epi) -cr2[ii]*epi) / ((1-purity)*cn + purity * cr2[ii])
        w2w2 <- w2^2
        ef2 <- w2 * ccfMeanMat +epi
        term1_breakpoint2 <- bn2[ii] * (log (ef2) - w2w2*ccfCovMat/(2 * ef2^2 ) )
        term2_breakpoint2 <- (dn2[ii] - bn2[ii]) * (log (1 - ef2) - w2w2*ccfCovMat/(2 * (1 - ef2)^2))
        term3_breakpoint2 <- logChoose(dn2[ii], bn2[ii])
        qq2 <- rowSums( respMat *  (term1_breakpoint2 + term2_breakpoint2 + term3_breakpoint2))
        qq2[1] = NA
        maxQq2 <- which.max(qq2)
        bv2[ii] <- bvPool2[maxQq2]
        bv2_sub1[ii] <- sub_cn2_mults$X[maxQq2]
        bv2_sub2[ii] <- sub_cn2_mults$Y[maxQq2]


        # m_upper2 = (1-epi)*( (1-purity)*cn + purity*cr2[ii] ) / (ccfMean*purity)*(1-epi) +
        #   epi*cr2[ii]/(1-epi)
        # m_lower2 = (0-epi)*( (1-purity)*cn + purity*cr2[ii] ) / (ccfMean*purity)*(1-epi) +
        #   epi*cr2[ii]/(1-epi)
        # bv2_old <- bv2[ii]
        # upper <- min(c(m_upper2, major_cn2[ii]))
        # lower <- max( c(m_lower2, 1) )
        # if (lower >= upper & upper == major_cn2[ii] ) {
        #   bv2[ii] = major_cn2[ii]
        # } else {
        #   tmp <- try(suppressWarnings( uniroot(
        #     ELBO_bv_g, c(lower, upper),
        #     bn = bn2[ii], dn=dn2[ii], cr=cr2[ii], cn=cn,
        #     purity=purity, ccfMean=ccfMean, ccfCov=ccfCov,
        #     epi=epi, responsibility = responsibility[ii,])$root), T)
        #
        #
        #   if (!is.numeric(tmp)) {
        #     bv2[ii] <- min( bv2_old, upper)
        #   } else {
        #     bv2[ii] <- tmp
        #   }
        # }

      } else {
        # clonal cn
        bvPool2 <- 1:max_mult_cn2_sub1[ii]

        qq2 <- rep(NA, length(bvPool2))
        for (jj in seq_along(bvPool2) ) {

          w2 <- purity * (bvPool2[jj] *(1-epi) -cr2[ii]*epi) / ((1-purity)*cn + purity * cr2[ii])
          w2w2 <- w2^2
          ef2 <- w2 * ccfMean +epi

          term1_breakpoint2 <- bn2[ii] * (log (ef2) - w2w2*ccfCov/(2 * ef2^2 ) )
          term2_breakpoint2 <- (dn2[ii] - bn2[ii]) * (log (1 - ef2) - w2w2*ccfCov/(2 * (1 - ef2)^2)  )
          term3_breakpoint2 <- logChoose(dn2[ii], bn2[ii])

          qq2[jj] <- sum ( responsibility[ii, ] *  (term1_breakpoint2 + term2_breakpoint2 + term3_breakpoint2)  )
        }
        bv2[ii] <- bvPool2[which.max(qq2)]
      }

    }
    model$bv1 <- bv1
    model$bv2 <- bv2
    model$bv1_sub1 <- bv1_sub1
    model$bv1_sub2 <- bv1_sub2
    model$bv2_sub1 <- bv2_sub1
    model$bv2_sub2 <- bv2_sub2
  }

  # estimate hyper-parameters
  if (fit_hyper) {
    # Elogpi <- sum(digamma(dirichletConcentration) - digamma(sum(dirichletConcentration)))
    # tmpOld <- model$dirichletConcentration0
    # jj <- 0
    # upper <- 1e-2
    # lower <- 1e-99
    # while (!is.numeric(tmp)) {
    #   if (jj > 1000) {
    #     break
    #   }
    #   jj <- jj + 1
    #   if (jj > 1) {upper <- upper + 1e-3}
    #   tmp <- try(suppressWarnings( uniroot(
    #     function(x) {
    #       term1 <- digamma(k*x) - digamma(x)
    #       return(term1+Elogpi/k)
    #     },
    #     c(lower, upper), extendInt = "yes")$root), T)
    # }
    # if (jj > 1000) {
    #   model$dirichletConcentration0 <- tmpOld
    #   cat(jj)
    # } else {
    #   model$dirichletConcentration0 <- tmp
    # }
    model$normalMean <- mean(ccfMean)
    model$invWhishartScale <-mean((ccfMean - model$normalMean)^2 + ccfCov)
  }

  model$dirichletConcentration <- dirichletConcentration
  model$ccfMean <- ccfMean
  model$ccfCov <- ccfCov

  model
}


############ Variational-Expectation ############
VarationalExpectationStep_sv <- function(bn1, dn1, cn, cr1,
                                         bn2, dn2, cr2,
                                         epi, purity, model, no.weights = FALSE) {

  bv1 <- model$bv1
  bv2 <- model$bv2
  dirichletConcentration <- model$dirichletConcentration	# Dirichlet
  ccfMean <- model$ccfMean
  ccfCov <- model$ccfCov
  Epi <- model$Epi

  n <- length(bn1)
  if (is.matrix(ccfMean)) {
    d <- nrow(ccfMean)
    k <- ncol(ccfMean)
  } else {
    d <- 1
    k <- length(ccfMean)
  }

  Epbk <- array(0, dim=c(n,k))

  w1 <- purity * (bv1 *(1-epi) -cr1*epi) / ((1-purity)*cn + purity * cr1)
  w1w1 <- w1^2

  w2 <- purity * (bv2 *(1-epi) -cr2*epi) / ((1-purity)*cn + purity * cr2)
  w2w2 <- w2^2

  # todo: only works for d = 1
  for(i in 1:k) {
    ef1 <- w1 * ccfMean[i] +epi
    ef2 <- w2 * ccfMean[i] +epi

    Epbk[,i] <- bn1 * (log (ef1) - w1w1*ccfCov[i]/(2 * ef1^2 ) ) +
      (dn1 - bn1) * (log (1 - ef1) - w1w1*ccfCov[i]/(2 * (1 - ef1)^2) ) +
      bn2 * (log (ef2) - w2w2*ccfCov[i]/(2 * ef2^2 ) ) +
      (dn2 - bn2) * (log (1 - ef2) - w2w2*ccfCov[i]/(2 * (1 - ef2)^2) )
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
VariationalLowerBound_sv <- function(bn1, dn1, cn, cr1,
                                     bn2, dn2, cr2,
                                     epi, purity, model) {

  bv1 <- model$bv1
  bv2 <- model$bv2
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

  n <- length(bn1)

  if (is.matrix(m)) {
    d <- nrow(m)
    k <- ncol(m)
  } else {
    d <- 1
    k <- length(m)
  }


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

  w1 <- purity * (bv1 *(1-epi) -cr1*epi) / ((1-purity)*cn + purity * cr1)
  w1w1 <- w1^2

  w2 <- purity * (bv2 *(1-epi) -cr2*epi) / ((1-purity)*cn + purity * cr2)
  w2w2 <- w2^2

  for(i in 1:k) {
    ef1 <- w1 * ccfMean[i] +epi
    ef2 <- w2 * ccfMean[i] +epi

    Epbk[,i] <- sum(R[,i]*( bn1 * (log (ef1) - w1w1*ccfCov[i]/(2 * ef1^2 ) ) +
                            (dn1 - bn1) * (log (1 - ef1) - w1w1*ccfCov[i]/(2 * (1 - ef1)^2) ) +
                            bn2 * (log (ef2) - w2w2*ccfCov[i]/(2 * ef2^2 ) ) +
                            (dn2 - bn2) * (log (1 - ef2) - w2w2*ccfCov[i]/(2 * (1 - ef2)^2) ) ) ) +
      sum( R[,i]* (  logChoose(dn1, bn1) + logChoose(dn2, bn2) ) )
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

#' Estimate purity
#' @param mydata ccube data
#' @param wgd whole genome duplication status
#' @param K total number of cluster in student't mixture
#' @param th threshold on weights to be eligible for clonal cluster
#' @return purity
#' @export
GetPurity_sv <- function(mydata, wgd=F, K = 6, th = 1.5e-2) {
  vtox<-function(v,nA,nB,tA,tB)
  {
    (nB - nB*v - nA*v) / (tA*v + tB*v + nB - tB -nA*v - nB*v)
  }

  tmpdata <- dplyr::filter(mydata, major_cn == minor_cn & major_cn != 0)

  if (!wgd) {
    tmpdata <- dplyr::filter(mydata, major_cn == minor_cn & major_cn == 1)
  }

  if (nrow(tmpdata) == 0) {
    return(NA)
  }
  tmpdata <- dplyr::mutate(dplyr::rowwise(tmpdata), ploidy = major_cn + minor_cn)
  tmpdata <- dplyr::mutate(tmpdata, vaf = var_counts/(var_counts+ref_counts))
  tmpdata <- dplyr::mutate(dplyr::rowwise(tmpdata), cp = vtox(vaf, 2, 0, ploidy/2, ploidy/2))
  tmpdata <- dplyr::filter(tmpdata, !is.infinite(cp) & !is.na(cp))

  if (K >= nrow(tmpdata)) {
    K = nrow(tmpdata) - 1
    if (K == 0) {K = 1}
  }
  if (nrow(tmpdata) == 1) {
    purity <- tmpdata$cp
  } else {
    res <- vbsmm(tmpdata$cp, init = K, tol = 1e-5,  verbose = F)
    uniqLabels <- sort(unique(res$label))
    ww <- uniqLabels[which( (table(res$label)/length(res$label) ) > th)]
    pool <- res$mu[ww]
    maxCp <- max(pool[pool<=1])
    purity  <- if (maxCp > 1) 1 else maxCp
  }

  return(purity)
}

#' Make Ccube results plot
#' @param ssm data
#' @param res Ccube result list
#' @param myColors colors
#' @param printPlot output flag
#' @param fn output file name
#' @return NULL
#' @export
MakeCcubeStdPlot_sv <- function(ssm, res, myColors=gg_color_hue(10), printPlot = F, fn = NULL) {

  if (printPlot) {
    pdf(fn, width=8, height=8)
  }

  par(mfrow=c(2,2))
  plot(ssm$ccube_ccf1, ssm$vaf1, col = myColors[res$label],
       xlab = "cancer cell fraction", ylab = "variant allele frequecy",
       main = "break point 1: ccf vs vaf \n(colored by cluster memebership)")
  cellularity <- unique(ssm$purity)
  ssm$total_cn1 =ssm$frac_cn1_sub1 * (ssm$major_cn1_sub1 + ssm$minor_cn1_sub1) +
    ssm$frac_cn1_sub2 *(ssm$major_cn1_sub2 + ssm$minor_cn1_sub2)
  uniqueTotCn = unique(ssm$total_cn1)
  xx = seq(0,2, length.out = 100)
  for (cn in uniqueTotCn) {
    for (i in 1:cn) {
      lines(MapVaf2CcfPyClone(xx, cellularity, 2, cn, cn, i, constraint = F), xx, lty = 6, col = 80)
    }
  }

  plot(ssm$ccube_ccf2, ssm$vaf2, col = myColors[res$label],
       xlab = "cancer cell fraction", ylab = "variant allele frequecy",
       main = "break point 2: ccf vs vaf \n(colored by cluster memebership)")
  cellularity <- unique(ssm$purity)
  ssm$total_cn2 =ssm$frac_cn2_sub1 * (ssm$major_cn2_sub1 + ssm$minor_cn2_sub1) +
    ssm$frac_cn2_sub2 *(ssm$major_cn2_sub2 + ssm$minor_cn2_sub2)
  uniqueTotCn = unique(ssm$total_cn2)
  xx = seq(0,2, length.out = 100)
  for (cn in uniqueTotCn) {
    for (i in 1:cn) {
      lines(MapVaf2CcfPyClone(xx, cellularity, 2, cn, cn, i, constraint = F), xx, lty = 6, col = 80)
    }
  }

  if (is.matrix(res$full.model$ccfMean)) {
    Emu <- res$full.model$ccfMean[,]
  } else {
    Emu <- res$full.model$ccfMean
  }

  if (is.matrix(res$full.model$ccfCov)) {
    Esigma <- res$full.model$ccfCov[,]
  } else {
    Esigma <- res$full.model$ccfCov
  }

  Epi <- res$full.model$Epi

  params <- data.frame(Emu, Esigma, Epi)
  xx <- seq(0,2,  length.out = 1000)
  ll <- 0

  for (j in seq_len(nrow(params))) {
    ll <- ll + params[j,]$Epi * dnorm(xx, mean = params[j,]$Emu, sd = sqrt(params[j,]$Esigma))
  }

  hist(ssm$ccube_ccf1, density=20, breaks=20, prob=TRUE,
       main = "ccf histograms \n(red: break point 1, blue: break point 2) + \ncluster uncertainties",
       xlab = "cancer cell fraction", col =  "red")
  hist(ssm$ccube_ccf2, density=20, breaks=20, prob=TRUE, col = "blue", add = T)
  lines(xx,ll, lwd=3, col = "darkred")

  numSnv <- table(res$label)
  uniqLabels = unique(res$label)
  names(numSnv) <- as.character(format(round(Emu[sort(uniqLabels)], 2), nsmall = 2))
  barplot(numSnv, las = 2, col = myColors[sort(uniqLabels)],
          xlab = "cluster mean", ylab="number of variants",
          main = "cluster prevalence")
  if (printPlot) {
    dev.off()
  }

}

#' Remove empty clusters
#' @param res Ccube result list
#' @param ssm Ccub input data
#' @param epi sequencing error
#' @return Ccube result list
#' @export
CullEmptyClusters_sv <- function(res, ssm, useEstep = T, epi = 1e-3) {

  idx <- which(! seq_along(res$full.model$ccfMean) %in% unique(res$label) )

  if (useEstep) {
    return(RemoveClusterAndReassignVariantsWithEstep_sv(res = res, removeIdx = idx, ssm = ssm, epi = epi))
  } else {
    return(RemoveClusterAndReassignVariants_sv(res = res, removeIdx = idx, ssm = ssm))
  }
}

#' Remove small clusters
#' @param res Ccube result list
#' @param ssm Ccub input data
#' @param th threshold for small clusters
#' @param epi sequencing error
#' @return Ccube result list
#' @export
CullSmallClusters_sv <- function(res, ssm, th = 1e-2, epi = 1e-3,  useEstep = T) {

  tt <- table(res$label)
  idx <- which( tt/sum(tt) < th )

  if (useEstep) {
    return(RemoveClusterAndReassignVariantsWithEstep_sv(res = res, removeIdx = idx, ssm = ssm, epi = epi))
  } else {
    return(RemoveClusterAndReassignVariants_sv(res = res, removeIdx = idx, ssm = ssm))
  }
}

#' Remove a (or more) cluster and reassign its data if the cluster is nonempty
#' @param res Ccube result list
#' @param  removeIdx clusters to remove
#' @param ssm data
#' @param label assigned labels if res doesn't have label variable
#' @return Ccube result list
#' @export
RemoveClusterAndReassignVariants_sv <- function(res, removeIdx, ssm = NULL, label = NULL) {

  if (length(removeIdx) == 0) {
    return(res)
  }

  if (!is.null(res$label)) {
    uniqLabels <- sort(unique(res$label))
  } else {
    uniqLabels <- sort(unique(label))
  }

  remainedLabels <- uniqLabels[which(!uniqLabels %in% removeIdx)]
  uniqLabels[ which(uniqLabels %in% removeIdx) ] <- NA
  newLabels <- match(res$label, uniqLabels)
  reassignIdx <- which(is.na(newLabels))
  reassignSsm <- ssm[reassignIdx, ]

  if (! is.null(res$full.model) ) {
    res$full.model$ccfMean <- res$full.model$ccfMean[-removeIdx]
    res$full.model$ccfCov <- res$full.model$ccfCov[-removeIdx]

    logRho <- res$full.model$logResponsibility[, -removeIdx]

    if (!is.matrix(logRho)) {
      logRho <- as.matrix(logRho)
    }

    res$full.model$logResponsibility <-
      if (length(res$label) ==length(res$full.model$ccfMean)) {
        bsxfun.se("-", logRho, logsumexp(logRho, 1), expandByRow = F)	# 10.49
      } else {
        bsxfun.se("-", logRho, logsumexp(logRho, 1))	# 10.49
      }

    if(!is.null(ssm) & length(reassignIdx) >0) {
      reassignList <- Assign(reassignSsm$ccube_ccf, res$full.model$ccfMean, res$full.model$ccfCov)
      res$full.model$logResponsibility[reassignIdx, ] <- reassignList$logR
    }

    res$full.model$responsibility <- exp(res$full.model$logResponsibility)

    res$label <- apply(res$full.model$responsibility, 1, which.max)

    res$full.model$dirichletConcentration <- res$full.model$dirichletConcentration0 + colSums(res$full.model$responsibility)

    if (!is.null(res$R)) {
      res$R=res$full.model$responsibility
    }

    if (!is.null(res$mu)) {
      res$mu=res$full.model$ccfMean
    }
    if (! is.null(res$full.model$Epi)) {
      Epi <- (res$full.model$dirichletConcentration + colSums(res$full.model$responsibility)) /
        (length(res$full.model$ccfMean) * res$full.model$dirichletConcentration0 + length(res$label))
      res$full.model$Epi <- Epi/sum(Epi)
    }

  } else {

    res$ccfMean <- res$ccfMean[-removeIdx]
    res$ccfCov <- res$ccfCov[-removeIdx]

    logRho <- res$logResponsibility[, -removeIdx]

    if (!is.matrix(logRho)) {
      logRho <- as.matrix(logRho)
    }

    res$logResponsibility <-
      if (length(label) ==length(res$ccfMean)) {
        bsxfun.se("-", logRho, logsumexp(logRho, 1), expandByRow = F)	# 10.49
      } else {
        bsxfun.se("-", logRho, logsumexp(logRho, 1))	# 10.49
      }

    if(!is.null(ssm) & length(reassignIdx) >0) {
      reassignList <- Assign(reassignSsm$ccube_ccf, res$ccfMean, res$ccfCov)
      res$logResponsibility[reassignIdx, ] <- reassignList$logR
    }

    res$responsibility <- exp(res$logResponsibility)

    res$dirichletConcentration <- res$dirichletConcentration0 + colSums(res$responsibility)
  }

  res
}

#' Remove a (or more) cluster and reassign its data if the cluster is nonempty
#' @param res Ccube result list
#' @param  removeIdx clusters to remove
#' @param ssm data
#' @param label assigned labels if res doesn't have label variable
#' @param epi sequencing error
#' @return Ccube result list
#' @export
RemoveClusterAndReassignVariantsWithEstep_sv <- function(res, removeIdx, ssm = NULL, label = NULL, epi = 1e-3) {

  if (length(removeIdx) == 0) {
    return(res)
  }

  if (!is.null(res$label)) {
    uniqLabels <- sort(unique(res$label))
  } else {
    uniqLabels <- sort(unique(label))
  }

  remainedLabels <- uniqLabels[which(!uniqLabels %in% removeIdx)]
  uniqLabels[ which(uniqLabels %in% removeIdx) ] <- NA
  newLabels <- match(res$label, uniqLabels)
  reassignIdx <- which(is.na(newLabels))
  reassignSsm <- ssm[reassignIdx, ]

  if (! is.null(res$full.model) ) {
    res$full.model$ccfMean <- t(as.matrix(res$full.model$ccfMean[-removeIdx]))
    res$full.model$ccfCov <- t(as.matrix(res$full.model$ccfCov[-removeIdx]))

    logRho <- res$full.model$logResponsibility[, -removeIdx]

    if (!is.matrix(logRho)) {
      logRho <- as.matrix(logRho)
    }

    res$full.model$logResponsibility <-
      if (length(res$label) ==length(res$full.model$ccfMean)) {
        bsxfun.se("-", logRho, logsumexp(logRho, 1), expandByRow = F)	# 10.49
      } else {
        bsxfun.se("-", logRho, logsumexp(logRho, 1))	# 10.49
      }

    res$full.model$responsibility <- exp(res$full.model$logResponsibility)
    if(!is.null(ssm)) {
      res$full.model$dirichletConcentration <- res$full.model$dirichletConcentration0 + colSums(res$full.model$responsibility)
      res$full.model <- VarationalExpectationStep_sv(bn1 = ssm$var_counts1,
                                                     dn1 = ssm$ref_counts1 + ssm$var_counts1,
                                                     cn = unique(ssm$normal_cn),
                                                     cr1 = ssm$frac_cn1_sub1 * (ssm$major_cn1_sub1 + ssm$minor_cn1_sub1) +
                                                       ssm$frac_cn1_sub2 *(ssm$major_cn1_sub2 + ssm$minor_cn1_sub2),
                                                     bn2 = ssm$var_counts2,
                                                     dn2 = ssm$ref_counts2 + ssm$var_counts2,
                                                     cr2 = ssm$frac_cn2_sub1 * (ssm$major_cn2_sub1 + ssm$minor_cn2_sub1) +
                                                       ssm$frac_cn2_sub2 *(ssm$major_cn2_sub2 + ssm$minor_cn2_sub2),
                                                     epi = epi,
                                                     purity = unique(ssm$purity),
                                                     model = res$full.model)

      res$full.model <- SortClusters(res$full.model)

      res$label <- apply(res$full.model$responsibility, 1, which.max)
      res$full.model$dirichletConcentration <- res$full.model$dirichletConcentration0 + colSums(res$full.model$responsibility)
    }



    if (!is.null(res$R)) {
      res$R=res$full.model$responsibility
    }

    if (!is.null(res$mu)) {
      res$mu=res$full.model$ccfMean
    }
    if (! is.null(res$full.model$Epi)) {
      Epi <- (res$full.model$dirichletConcentration + colSums(res$full.model$responsibility)) /
        (length(res$full.model$ccfMean) * res$full.model$dirichletConcentration0 + length(res$label))
      res$full.model$Epi <- Epi/sum(Epi)
    }

  }

  res
}

#' Remove a (or more) cluster and reassign its data if the cluster is nonempty
#' @param res Ccube result list
#' @param removeIdx clusters to remove
#' @param ssm data
#' @param label assigned labels if res doesn't have label variable
#' @param tol stopping condition
#' @param maxiter maximum iteration
#' @param epi sequencing error
#' @param verbose show progress
#' @param fit_mult flag to estimate multiplicities
#' @param fit_hyper flag to estimate hyperparameters
#' @return Ccube result list
#' @export
RemoveClusterAndReassignVariantsWithEMsteps_sv <- function(res, removeIdx, ssm = NULL, label = NULL, tol = 1e-8, maxiter = 100, epi = 1e-3, verbose = F,
                                                        fit_mult = T, fit_hyper = T) {

  if (length(removeIdx) == 0) {
    return(res)
  }

  if (!is.null(res$label)) {
    uniqLabels <- sort(unique(res$label))
  } else {
    uniqLabels <- sort(unique(label))
  }

  remainedLabels <- uniqLabels[which(!uniqLabels %in% removeIdx)]
  uniqLabels[ which(uniqLabels %in% removeIdx) ] <- NA
  newLabels <- match(res$label, uniqLabels)
  reassignIdx <- which(is.na(newLabels))
  reassignSsm <- ssm[reassignIdx, ]

  if (! is.null(res$full.model) ) {
    res$full.model$ccfMean <- t(as.matrix(res$full.model$ccfMean[-removeIdx]))
    res$full.model$ccfCov <-  t(as.matrix(res$full.model$ccfCov[-removeIdx]))

    logRho <- res$full.model$logResponsibility[, -removeIdx]

    if (!is.matrix(logRho)) {
      logRho <- as.matrix(logRho)
    }

    res$full.model$logResponsibility <-
      if (length(res$label) ==length(res$full.model$ccfMean)) {
        bsxfun.se("-", logRho, logsumexp(logRho, 1), expandByRow = F)	# 10.49
      } else {
        bsxfun.se("-", logRho, logsumexp(logRho, 1))	# 10.49
      }

    res$full.model$responsibility <- exp(res$full.model$logResponsibility)
    res$full.model$dirichletConcentration <- res$full.model$dirichletConcentration0 + colSums(res$full.model$responsibility)

    if(!is.null(ssm)) {
      ll = rep(-Inf, maxiter)
      vbiter = 1
      converged = F
      degenerated = F
      while (!converged & vbiter < maxiter & !degenerated) {
        vbiter = vbiter + 1

        res$full.model <- VarationalExpectationStep_sv(bn1 = ssm$var_counts1,
                                                       dn1 = ssm$ref_counts1 + ssm$var_counts1,
                                                       cn = unique(ssm$normal_cn),
                                                       cr1 = ssm$frac_cn1_sub1 * (ssm$major_cn1_sub1 + ssm$minor_cn1_sub1) +
                                                         ssm$frac_cn1_sub2 *(ssm$major_cn1_sub2 + ssm$minor_cn1_sub2),
                                                       bn2 = ssm$var_counts2,
                                                       dn2 = ssm$ref_counts2 + ssm$var_counts2,
                                                       cr2 = ssm$frac_cn2_sub1 * (ssm$major_cn2_sub1 + ssm$minor_cn2_sub1) +
                                                         ssm$frac_cn2_sub2 *(ssm$major_cn2_sub2 + ssm$minor_cn2_sub2),
                                                       epi = epi,
                                                       purity = unique(ssm$purity),
                                                       model = res$full.model)

        res$full.model <- VariationalMaximimizationStep_sv(bn1 = ssm$var_counts1,
                                                        dn1 = ssm$ref_counts1 + ssm$var_counts1,
                                                        cn = unique(ssm$normal_cn),
                                                        cr1 = ssm$frac_cn1_sub1 * (ssm$major_cn1_sub1 + ssm$minor_cn1_sub1) +
                                                          ssm$frac_cn1_sub2 *(ssm$major_cn1_sub2 + ssm$minor_cn1_sub2),
                                                        max_mult_cn1_sub1 = ssm$major_cn1_sub1 + ssm$minor_cn1_sub1,
                                                        max_mult_cn1_sub2 = ssm$major_cn1_sub2 + ssm$minor_cn1_sub2,
                                                        frac_cn1_sub1 = ssm$frac_cn1_sub1,
                                                        frac_cn1_sub2 = ssm$frac_cn1_sub2,
                                                        subclonal_cn1 = ssm$subclonal_cn1,
                                                        bn2 = ssm$var_counts2,
                                                        dn2 = ssm$ref_counts2 + ssm$var_counts2,
                                                        cr2 = ssm$frac_cn2_sub1 * (ssm$major_cn2_sub1 + ssm$minor_cn2_sub1) +
                                                          ssm$frac_cn2_sub2 *(ssm$major_cn2_sub2 + ssm$minor_cn2_sub2),
                                                        max_mult_cn2_sub1 = ssm$major_cn2_sub1 + ssm$minor_cn2_sub1,
                                                        max_mult_cn2_sub2 = ssm$major_cn2_sub2 + ssm$minor_cn2_sub2,
                                                        frac_cn2_sub1 = ssm$frac_cn2_sub1,
                                                        frac_cn2_sub2 = ssm$frac_cn2_sub2,
                                                        subclonal_cn2 = ssm$subclonal_cn2,
                                                        epi = epi,
                                                        purity = unique(ssm$purity),
                                                        model = res$full.model,fit_mult = fit_mult,
                                                        fit_hyper = fit_hyper)

        ll[vbiter] = VariationalLowerBound_sv(bn1 = ssm$var_counts1,
                                              dn1 = ssm$ref_counts1 + ssm$var_counts1,
                                              cn = unique(ssm$normal_cn),
                                              cr1 = ssm$frac_cn1_sub1 * (ssm$major_cn1_sub1 + ssm$minor_cn1_sub1) +
                                                ssm$frac_cn1_sub2 *(ssm$major_cn1_sub2 + ssm$minor_cn1_sub2),
                                              bn2 = ssm$var_counts2,
                                              dn2 = ssm$ref_counts2 + ssm$var_counts2,
                                              cr2 = ssm$frac_cn2_sub1 * (ssm$major_cn2_sub1 + ssm$minor_cn2_sub1) +
                                                ssm$frac_cn2_sub2 *(ssm$major_cn2_sub2 + ssm$minor_cn2_sub2),
                                              epi = epi,
                                              purity = unique(ssm$purity),
                                              model = res$full.model)/n
        converged <- abs(ll[vbiter] - ll[vbiter-1]) < (tol * abs(ll[vbiter]))
        degenerated <- (ll[vbiter] - ll[vbiter-1]) < 0
        if(verbose) cat(sprintf("\rVB-EM-%d: L = %.8f \r", vbiter, ll[vbiter]))
      }


      res$full.model <- SortClusters(res$full.model)

      res$label <- apply(res$full.model$responsibility, 1, which.max)
    }



    if (!is.null(res$R)) {
      res$R=res$full.model$responsibility
    }

    if (!is.null(res$mu)) {
      res$mu=res$full.model$ccfMean
    }
    if (! is.null(res$full.model$Epi)) {
      Epi <- (res$full.model$dirichletConcentration + colSums(res$full.model$responsibility)) /
        (length(res$full.model$ccfMean) * res$full.model$dirichletConcentration0 + length(res$label))
      res$full.model$Epi <- Epi/sum(Epi)
    }

  }

  res
}



#' Merge clusters
#' @param res ccube results list
#' @param ssm ccube data frame
#' @param tol stopping condition in VBEM
#' @param maxiter maximum iteration in VBEM
#' @param epi sequencing error
#' @return res ccube results list
#' @export
MergeClusters_sv <- function(res = res, ssm = ssm, tol = 1e-8, maxiter = 100, epi = 1e-3) {

  res <- CullEmptyClusters(res = res, ssm = ssm)

  if (is.matrix(res$full.model$ccfMean)) {
    ccfCentersMap = res$full.model$ccfMean[,]
  } else {
    ccfCentersMap = res$full.model$ccfMean
  }

  clusterWeights <- as.data.frame(table(res$label), stringsAsFactors = F)
  clusterWeights <- dplyr::mutate(clusterWeights, Var1 = as.integer(Var1))
  ccfDistMat <- dist(ccfCentersMap)
  idx <- as.data.frame(which(as.matrix( ccfDistMat )< 0.1, arr.ind = TRUE), stringsAsFactors =F)
  idx <- dplyr::filter(idx, row != col)
  idx <- dplyr::left_join(idx, clusterWeights, by = c("col"="Var1"))
  idx <- dplyr::rename(idx, col_weights = Freq)
  idx <- dplyr::left_join(idx, clusterWeights, by = c("row"="Var1"))
  idx <- dplyr::rename(idx, row_weights = Freq)
  idx <- dplyr::mutate(rowwise(idx), remove_idx = c(row, col)[which.min(c(row_weights, col_weights))] )
  removeIdx <- unique(idx$remove_idx)
  res$mergeCluster <- length(removeIdx)>0
  return(RemoveClusterAndReassignVariantsWithEMsteps_sv(res = res,
                                                     removeIdx = removeIdx, ssm = ssm,
                                                     tol = tol, maxiter = maxiter, epi = epi))
}

#' Write files in PCAWG-11 formats, works for both CcubeCore and ccube_m6 output
#' @param ssm data
#' @param res Ccube result list
#' @param resultFolder path to file
#' @param sampleName sample name
#' @return NULL
#' @export
WritePcawgFormats_sv <- function(ssm, res, resultFolder, sampleName,
                              allFormats = F, basicFormats = T,
                              outputMult = F, outputAssignProb = F, outputAssign = F,
                              outputSubStruct = F, outputCcm = F, outputCcmIdx = F) {
  if (basicFormats) {
    outputMult <- T
    outputAssign <- T
    outputSubStruct <- T
  }

  if (allFormats) {
    outputMult <- T
    outputAssignProb <- T
    outputAssign <- T
    outputSubStruct <- T
    outputCcm <- T
    outputCcmIdx <- T
  }

  ## output calibration format
  uniqLabels <- unique(res$label)
  dir.create(resultFolder, recursive = T)

  id <- do.call(rbind, strsplit(as.character(ssm$mutation_id), "_", fixed = T))

  # Multiplicity
  if (outputMult) {
    mult <- data.frame(chr = id[,1], pos = id[,2])
    mult$tumour_copynumber <- ssm$major_cn+ssm$minor_cn
    mult$multiplicity <- ssm$ccube_mult
    fn <- paste0(resultFolder, "/",
                 sampleName, "_multiplicity.txt")
    write.table(mult, file = fn, sep = "\t", row.names = F, quote = F)
    shellCommand <- paste0("gzip -f ", fn)
    system(shellCommand, intern = TRUE)
    rm(mult)
  }


  # Assignment Prob
  if(outputAssignProb){
    mutAssign <- data.frame(chr = id[,1], pos = id[,2])

    if (length(uniqLabels) == 1) {
      mutR = data.frame(res$full.model$responsibility)
      colnames(mutR) <- "cluster_1"
    } else {
      mutR <- data.frame(res$full.model$responsibility[, sort(uniqLabels)])
      colnames(mutR) <- paste0("cluster_", seq_along(uniqLabels))
    }

    mutAssign <- data.frame(mutAssign, mutR)
    fn <- paste0(resultFolder, "/",
                 sampleName, "_assignment_probability_table.txt")
    write.table(mutAssign, file = fn, sep = "\t", row.names = F, quote = F)
    shellCommand <- paste0("gzip -f ", fn)
    system(shellCommand, intern = TRUE)
    rm(mutR, mutAssign)
  }

  # Assignment
  if (outputAssign) {
    mutAssign <- data.frame(chr = id[,1], pos = id[,2])
    if (length(uniqLabels) == 1) {
      mutAssign$cluster = 1
    } else {
      mutAssign$cluster <- res$label
    }
    fn <- paste0(resultFolder, "/",
                 sampleName, "_mutation_assignments.txt")
    write.table(mutAssign, file = fn, sep = "\t", row.names = F, quote = F)
    shellCommand <- paste0("gzip -f ", fn)
    system(shellCommand, intern = TRUE)
    rm(mutAssign)
  }

  # subclonal structure
  if (outputSubStruct) {
    cellularity <- unique(ssm$purity)
    clusterCertainty <- as.data.frame(table(res$label), stringsAsFactors = F)
    clusterCertainty <- dplyr::rename(clusterCertainty, cluster = Var1, n_ssms = Freq)
    clusterCertainty$proportion <- res$full.model$ccfMean[as.integer(clusterCertainty$cluster)] * cellularity
    clusterCertainty$cluster <- seq_along(uniqLabels)
    fn <- paste0(resultFolder, "/",
                 sampleName, "_subclonal_structure.txt")
    write.table(clusterCertainty, file = fn, sep = "\t", row.names = F, quote = F)
    shellCommand <- paste0("gzip -f ", fn)
    system(shellCommand, intern = TRUE)
  }


  # ccm
  if (outputCcm) {
    coClustMat <- Matrix::tcrossprod(res$full.model$responsibility)
    diag(coClustMat) <- 1
    fn = paste0(resultFolder, "/", sampleName, "_coassignment_probabilities.txt")
    write.table(coClustMat, file = fn, sep = "\t", col.names = F, row.names = F, quote = F )
    shellCommand <- paste0("gzip -f ", fn)
    system(shellCommand, intern = TRUE)
    rm(coClustMat)
  }


  # ccm index
  if (outputCcmIdx) {
    indexFile <- cbind(id[,1], id[, 2], seq_along(id[,1]))
    colnames(indexFile) <- c("chr", "pos", "col")
    fn = paste0(resultFolder, "/", sampleName, "_index.txt")
    write.table(indexFile, file = fn, sep = "\t", row.names = F, quote = F )
    shellCommand <- paste0("gzip -f ", fn)
    system(shellCommand, intern = TRUE)
  }

}
