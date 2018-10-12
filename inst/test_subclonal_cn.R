rm(list = ls())
library(dplyr)
library(ccube)
library(doParallel)
source("R/ccube_sv.R")
source("R/util.R")
registerDoParallel(cores=3)

set.seed(1234)

numSv <- 500
maxiter <- 100
init =5  # number of clusters
numOfClusterPool = 6
numOfRepeat = 1
baseDepth = 40
ccfSet <- c(1, 0.3, 0.7) # true ccf pool
ccfTrue <- sample(ccfSet, numSv, c(0.5,0.3,0.2), replace = T)

purity <- 0.8


cnPoolMaj <- c(1,2,3,4)
cnPoolMin <- c(0,1,2)


# 1st break point
cnProfile1 <- sample(cnPoolMaj, numSv, c(0.25, 0.25, 0.25,0.25), replace =T)
cnProfile2 <- sample(cnPoolMin, numSv, c(1/3, 1/3, 1/3), replace =T)

cnProfileTot <- cnProfile1 + cnProfile2

cnProfile <- cbind(cnProfile1, cnProfile2, cnProfileTot )
cnProfile <- t( apply(cnProfile, 1, sort) )


mydata <- data.frame(mutation_id = paste0("ss", seq_len(numSv)) ,
                     ccf_true = ccfTrue,
                     minor_cn1 = cnProfile[,1],
                     major_cn1 = cnProfile[,2],
                     total_cn1 = cnProfile[,3],
                     stringsAsFactors = F)


mydata$purity <- purity
mydata$normal_cn <- 2

mydata <- mutate(rowwise(mydata),
                 mult1 = sample(c(1,if (major_cn1 ==1) { 1 } else {major_cn1}), 1),
                 vaf1 = cp2ap(ccf_true, purity, normal_cn, total_cn1, total_cn1, mult1),
                 total_counts1 = rpois(1, total_cn1/2 * baseDepth),
                 var_counts1 = rbinom(1, total_counts1, vaf1),
                 ref_counts1 = total_counts1 - var_counts1)

# 2nd break point

cnProfile1 <- sample(cnPoolMaj, numSv, c(0.5, 0.2, 0.2,0.1), replace =T)
cnProfile2 <- sample(cnPoolMin, numSv, c(0.3, 0.5, 0.2), replace =T)

cnProfileTot <- cnProfile1 + cnProfile2

cnProfile <- cbind(cnProfile1, cnProfile2, cnProfileTot )
cnProfile <- t( apply(cnProfile, 1, sort) )

mydata$minor_cn2 = cnProfile[,1]
mydata$major_cn2 = cnProfile[,2]
mydata$total_cn2 = cnProfile[,3]

mydata <- mutate(rowwise(mydata),
                 mult2 = sample(c(1,if (major_cn2 ==1) { 1 } else {major_cn2}), 1),
                 vaf2 = cp2ap(ccf_true, purity, normal_cn, total_cn2, total_cn2, mult2),
                 total_counts2 = rpois(1, total_cn1/2 * baseDepth ),
                 var_counts2 = rbinom(1, total_counts2, vaf2),
                 ref_counts2 = total_counts2 - var_counts2
                 )
mydata$subclonal_cn = c(rep(T,numSv/2), rep(T, numSv/2))


mydata <- GetCcf_sv(mydata, use="use_base")

dn1 <- mydata$ref_counts1 + mydata$var_counts1
bn1 <- mydata$var_counts1
cn <- unique(mydata$normal_cn)
cr1 <- mydata$major_cn1 + mydata$minor_cn1
major_cn1 <- mydata$major_cn1
bv1 <- mydata$mult1

dn2 <- mydata$ref_counts2 + mydata$var_counts2
bn2 <- mydata$var_counts2
cr2 <- mydata$major_cn2 + mydata$minor_cn2
major_cn2 <- mydata$major_cn2
bv2 <- mydata$mult2

subclonal_cn <- mydata$subclonal_cn




purity <- unique(mydata$purity)
rawCcf <- mydata$ccf
rawCcf <- as.matrix(rawCcf)
n <- nrow(rawCcf)
d <- ncol(rawCcf)

X <- t(rawCcf)

k <-  numOfClusterPool
prior <- list(
  dirichletConcentration = 1e2,
  normalMean = 0.5,
  invWhishartScale = var(rawCcf)*(d+1)# M = inv(W)
)


L <- rep(-Inf, maxiter)
converged <- FALSE
degenerated <- FALSE
vbiter <- 1

model <- list()

initParams <- ccube:::initialization(X, init, prior) # initialize responsibility and hidden scale
model$responsibility <- initParams$R
model$ccfMean <- initParams$ccfMean
model$ccfCov <- initParams$ccfCov
model$bv1 <- bv1
model$bv2 <- bv2
model$dirichletConcentration0 <- prior$dirichletConcentration
model$normalMean <- prior$normalMean
model$invWhishartScale <- prior$invWhishartScale
for (vbiter in 1:maxiter) {
  cat(vbiter ,'\n')
  model <- VariationalMaximimizationStep_sv(bn1, dn1, cn, cr1, major_cn1,
                                            bn2, dn2, cr2, major_cn2,
                                            epi=1e-3, purity, subclonal_cn, model,
                                            fit_mult = T, fit_hyper = T)

  model <- VarationalExpectationStep_sv(bn1, dn1, cn, cr1,
                                        bn2, dn2, cr2,
                                        epi=1e-3, purity, model)

  L[vbiter] <- VariationalLowerBound_sv(bn1, dn1, cn, cr1,
                                        bn2, dn2, cr2,
                                        epi=1e-3, purity, model)/n

}

model <- ccube:::SortClusters(model)

if (init > 1) {
  label <- apply(model$responsibility, 1, which.max)
  nk <- colSums(model$responsibility)
  Epi <- (model$dirichletConcentration + nk) / (k*model$dirichletConcentration0 + n)
  model$Epi <- Epi/sum(Epi)
} else {
  label <- rep(1, n)
  model$Epi <- 1
}

res = list(label=label, full.model=model, L=L)

mydata$ccube_ccf_mean <- res$full.model$ccfMean[res$label]
mydata$ccube_double_mult1 <- res$full.model$bv1
mydata$ccube_double_mult2 <- res$full.model$bv2
mydata <- mutate(rowwise(mydata),
                 vaf1 = var_counts1/(var_counts1+ref_counts1),
                 ccube_double_ccf1 = MapVaf2CcfPyClone(vaf1,
                                                       purity,
                                                       normal_cn,
                                                       major_cn1 + minor_cn1,
                                                       major_cn1 + minor_cn1,
                                                       ccube_double_mult1,
                                                       constraint=F),
                 vaf2 = var_counts2/(var_counts2+ref_counts2),
                 ccube_double_ccf2 = MapVaf2CcfPyClone(vaf2,
                                                       purity,
                                                       normal_cn,
                                                       major_cn2 + minor_cn2,
                                                       major_cn2 + minor_cn2,
                                                       ccube_double_mult2,
                                                       constraint=F)
)


mydata <- mutate(rowwise(mydata),
                 vaf1 = var_counts1/(var_counts1+ref_counts1),
                 ccube_double_ccf1 = MapVaf2CcfPyClone(vaf1,
                                                       purity,
                                                       normal_cn,
                                                       major_cn1 + minor_cn1,
                                                       major_cn1 + minor_cn1,
                                                       ccube_double_mult1,
                                                       constraint=F),
                 true_obs_ccf1 = MapVaf2CcfPyClone(vaf1,
                                                   purity,
                                                   normal_cn,
                                                   major_cn1 + minor_cn1,
                                                   major_cn1 + minor_cn1,
                                                   mult1,
                                                   constraint=F),
                 vaf2 = var_counts2/(var_counts2+ref_counts2),
                 ccube_double_ccf2 = MapVaf2CcfPyClone(vaf2,
                                                       purity,
                                                       normal_cn,
                                                       major_cn2 + minor_cn2,
                                                       major_cn2 + minor_cn2,
                                                       ccube_double_mult2,
                                                       constraint=F),
                 true_obs_ccf2 = MapVaf2CcfPyClone(vaf2,
                                                   purity,
                                                   normal_cn,
                                                   major_cn2 + minor_cn2,
                                                   major_cn2 + minor_cn2,
                                                   mult2,
                                                   constraint=F)
)


mydata$ccube_ccf_mean <- res$full.model$ccfMean[res$label]
mydata$ccube_double_mult1 <- res$full.model$bv1
mydata$ccube_double_mult2 <- res$full.model$bv2
mydata <- mutate(rowwise(mydata),
              vaf1 = var_counts1/(var_counts1+ref_counts1),
              ccube_double_ccf1 = MapVaf2CcfPyClone(vaf1,
                                             purity,
                                             normal_cn,
                                             major_cn1 + minor_cn1,
                                             major_cn1 + minor_cn1,
                                             ccube_double_mult1,
                                             constraint=F),
              vaf2 = var_counts2/(var_counts2+ref_counts2),
              ccube_double_ccf2 = MapVaf2CcfPyClone(vaf2,
                                             purity,
                                             normal_cn,
                                             major_cn2 + minor_cn2,
                                             major_cn2 + minor_cn2,
                                             ccube_double_mult2,
                                             constraint=F),
              ccube_ccf1 = ccube_double_ccf1,
              ccube_ccf2 = ccube_double_ccf2
)

mydata$ccube_mult1 <- res$full.model$bv1
mydata$ccube_mult2 <- res$full.model$bv2

MakeCcubeStdPlot_sv(mydata, res)

label1 = res$label

myColors=gg_color_hue(10)

par(mfrow=c(1,2))
plot(mydata$true_obs_ccf1, mydata$ccube_double_ccf1, col = myColors[label1],
     xlim = c(0, max( c(mydata$true_obs_ccf1, mydata$ccube_double_ccf1) ) ),
     ylim = c(0, max( c(mydata$true_obs_ccf1, mydata$ccube_double_ccf1) ) ),
     xlab = "true ccf", ylab = "estimated ccf", main = "double model: 1st break point")
points( seq(0, max( c(mydata$true_obs_ccf1, mydata$ccube_double_ccf1) ), length.out = 100 ),
        seq(0, max( c(mydata$true_obs_ccf1, mydata$ccube_double_ccf1) ), length.out = 100 ),
        type = "l" )

plot(mydata$true_obs_ccf2, mydata$ccube_double_ccf2, col = myColors[label1],
     xlim = c(0, max( c(mydata$true_obs_ccf2, mydata$ccube_double_ccf2) ) ),
     ylim = c(0, max( c(mydata$true_obs_ccf2, mydata$ccube_double_ccf2) ) ),
     xlab = "true ccf", ylab = "estimated ccf", main = "double model: 2nd break point"
)

points( seq(0, max( c(mydata$true_obs_ccf2, mydata$ccube_double_ccf2) ), length.out = 100 ),
        seq(0, max( c(mydata$true_obs_ccf2, mydata$ccube_double_ccf2) ), length.out = 100 ),
        type = "l" )
