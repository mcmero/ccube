rm(list = ls())
library(dplyr)
library(gtools)
library(doParallel)
library(foreach)
registerDoParallel(cores=4)


ComputeBic = function(mllh, nobs, nparam)
{
  bic = mllh - ((nparam)/2 * log(nobs))
  return(bic)
}

numSnv <- 500
ccfSet <- c(0.96, 0.6, 0.4)
ccfTrue <- sample(ccfSet, numSnv, c(0.6,0.2,0.2), replace = T)

purity <- 0.47
cnPool <- as.integer(c(1,2,3))
totalCN <- sample(cnPool, numSnv, c(0.2, 0.5, 0.3), replace =T)

mydata <- data.frame(mutation_id = seq_along(ccfTrue),
                     true_ccf_label=ccfTrue,
                     total_cn = totalCN, normal_cn = as.integer(2), purity = purity,
                     stringsAsFactors = F)

mydata <- mutate(rowwise(mydata),
                 true_mult = sample(1:(if (total_cn ==1) {1} else {total_cn-1}), 1) )

mydata <- mutate(rowwise(mydata),
                 minor_cn =if (total_cn==1) {as.integer(0)}
                 else{sample(1:floor(total_cn/2),1)} )

mydata <- mydata <- mutate(rowwise(mydata),
                           major_cn = total_cn - minor_cn)

mydata <- mutate(rowwise(mydata),
                 vaf = cp2ap(true_ccf_label, purity, normal_cn, total_cn, total_cn, true_mult))

totReads <- rpois(numSnv, 50)
mydata$total_counts <- totReads
mydata <- mutate(rowwise(mydata),
                 var_counts = rbinom(1, total_counts, vaf),
                 ref_counts = total_counts - var_counts)
mydata <- mutate(rowwise(mydata), true_ccf = MapVaf2CcfPyClone(var_counts/total_counts,
                                                                     purity,
                                                                     normal_cn,
                                                                     total_cn,
                                                                     total_cn,
                                                                     true_mult, constraint = F))


# source("./ccube_m4.R")
# res0 <-  ccube_m4(bn, dn, bv, cn, cr, rawCcf, purity, epi, init=K, tol = 1e-10, maxiter = 1e3,  verbose = F)

K =3
# results1 <- foreach(n = 1:K, .combine = c) %dopar% {
#   source("./BinaryDataMixture.R")
#   list(FitBinomialMixtureCcf(mydata$varReads,
#                              mydata$totReads, purity, 2 , mydata$cnProfile,
#                              mydata$cnProfile,
#                              mydata$mut, n, 1e-10, 100, estmate_mut = T))
# }
#
# maxLbIndex <- which.max(Map( function(x) max( ComputeBic(x$llh, numSnv, length(x$Eparam))), results1))
# lb1 <- unlist(Map( function(x) max(ComputeBic(x$llh, numSnv, length(x$Eparam))), results1))
# res1 <- results1[[maxLbIndex]]
#
# #
# vmeasure1 <- bitphylogenyR::vmeasureR(mydata$ccfTrue, res1$pred_label)
# cat("v-measure:", vmeasure1[,3])
# cat("\nEparam:", res1$Eparam)
# cat("\nEpi: ", res1$Epi)
# mydata$cluster1_ccf <- res1$Eparam[res1$pred_label]
# mydata$cluster1_mult <- res1$Emult
#
#



# results2 <- foreach(n = 1:K, .combine = c) %dopar% {
#   list(ccube_m6(mydata, epi=1e-3,
#                 init=n, tol = 1e-10, maxiter = 1e3,
#                 fit_mult = T, fit_hyper = T, use = "use_one", verbose = T))
#
# }
#
#
#
# maxLbIndex <- which.max(Map( function(x) max(x$L), results2))
# lb2 <- unlist(Map( function(x) max(x$L), results2))
# # lb2 <- unlist(Map( function(x) max(ComputeBic(x$L, numSnv, 2*length(x$mu))), results2))
# res2 <- results2[[maxLbIndex]]
#
res2 = ccube_m6(mydata, epi=1e-3,
         init=5, tol = 1e-10, maxiter = 1e3,
         fit_mult = T, fit_hyper = T, use = "use_base", verbose = T)

vmeasure2 <- bitphylogenyR::vmeasureR(mydata$true_ccf, res2$label)
cat("\nv-measure:", vmeasure2[,3])
cat("\nEparam:", res2$mu)
cat("\nEpi: ", res2$full.model$Epi)
mydata$mean_ccf <- res2$mu[res2$label]
mydata$mult <- res2$full.model$bv
#
mydata <- mutate(rowwise(mydata),
                 ccube_ccf = MapVaf2CcfPyClone(var_counts/(var_counts+ref_counts),
                                         purity,
                                         normal_cn,
                                         major_cn + minor_cn,
                                         major_cn + minor_cn,
                                         mult,
                                         constraint=F))
# k <-5
# source("./ccube_m6.R")
# res <- ccube_m6(bn, dn, bv, cn, cr, rawCcf, purity, epi, init=k, tol = 1e-10, maxiter = 1e3, estimate_mult =F,  verbose = T)
# vmeasure <- bitphylogenyR::vmeasureR(mydata$ccfTrue, res$label)
#
# source("./vbsmm.R")
# res2 = vbsmm(mydata$ccf2, init = 3, tol = 1e-5,  verbose = T)
# vmeasure2 <- bitphylogenyR::vmeasureR(mydata$ccfTrue, res2$label)
# # #
# res3 = FitBinomialMixtureCcf(mydata$varReads, mydata$totReads, purity, 2 , mydata$cnProfile,
#                              mydata$cnProfile, mydata$mut, k, 1e-5, 100, estmate_mut = T)
# vmeasure3 <- bitphylogenyR::vmeasureR(mydata$ccfTrue, res3$pred_label)

# ComputeCCfEstRisk <- function(d, x, t, cnNormal, cnRef, cnVar, bv) {
#   r1 = numeric(length = d)
#   r2 = r1
#   vaf <- cp2ap(x, t, cnNormal, cnRef, cnVar, bv)
#   for (b in 1:d) {
#     ccfEst <- MapVaf2CcfPyClone(b/d, t, cnNormal, cnRef, cnVar, bv, constraint = F)
#     r1[b] = (x - ccfEst)^2 * vaf ^ b * (1 - vaf)^(d-b)
#     r2[b] = choose(d, b)
#   }
#   rr = sum(r1 * r2, na.rm = T)
# }
#
#
# gg_color_hue <- function(n) {
#   hues = seq(15, 375, length=n+1)
#   hcl(h=hues, l=65, c=100)[1:n]
# }
#
# myColors <- gg_color_hue(10)
#
# ii = 1
# xx = seq(0,1, length.out = 100)
# dd = sapply(1:mydata$cnProfile[ii],
#             function(i)
#               sapply(xx, function(x)
#                 ComputeCCfEstRisk(mydata$totReads[ii], x, purity, 2,
#                                   mydata$cnProfile[ii], mydata$cnProfile[ii], i)))
# plot(xx, dd[,1], col=myColors[1])
# for (i in 2 : ncol(dd)) {
#   points(xx, dd[,i], col= myColors[i])
# }

#
#
# res3 <- vbsmm(mydata$ccf, K, tol = 1e-5, verbose = T)
# vmeasure3 <- bitphylogenyR::vmeasureR(mydata$ccfTrue, res3$label)
# cat("v-measure:", vmeasure3[,3])
# cat("\nEparam:", res3$mu)
# cat("\nEpi: ", res3$full.model$Epi)
# mydata$cluster3_ccf <- res3$mu[res3$label]

# n <- 3000
# d <- rpois(n, 50)
# f <- sort(sample(c(1.0, 0.6, 0.2), n, replace = T))
# b <- sapply(1:n, function(i) rbinom(1, d[i], f[i]))
#
# res1 <- FitBinomialMixtureVaf(b, d, 3, 1e-3, 100)
# vmeasure <- vmeasureR(f, res1$pred_label)
# cat("v-measure:", vmeasure[,3])
# cat("\nEparam:", res1$Eparam)
# cat("\nEpi: ", res1$Epi)
