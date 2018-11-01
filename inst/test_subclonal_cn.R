rm(list = ls())
library(dplyr)
library(ccube)
library(doParallel)

GenerateCopyNumberProfile <- function(cnPoolMaj, cnPoolMin,cnPoolMajFractions, cnPoolMinFractions, numVariants){
  tmp1 <- sample(cnPoolMaj, numVariants, cnPoolMajFractions, replace =T)
  tmp2 <- sample(cnPoolMin, numVariants, cnPoolMinFractions, replace =T)
  cnProfileTot <- tmp1 + tmp2
  cnProfile <- cbind(tmp1, tmp2, cnProfileTot )
  cnProfile <- t( apply(cnProfile, 1, sort) )
  return(cnProfile)
}

GenerateSubClonalCNProfile <- function(cnPoolMaj, cnPoolMin,
                                       cnPoolMajFractions, cnPoolMinFractions,
                                       numVariants, subclonal, ccfCN) {
  cnProfile <- matrix(-100, nrow = numVariants, ncol = 8)
  for (ii in 1:numVariants) {
    if (subclonal[ii]) {
      tmp1 = tmp2 = NA
      while (  identical(tmp1, tmp2) ) {
        tmp1 = GenerateCopyNumberProfile(cnPoolMaj, cnPoolMin,cnPoolMajFractions, cnPoolMinFractions, 1)
        tmp2 = GenerateCopyNumberProfile(cnPoolMaj, cnPoolMin,cnPoolMajFractions, cnPoolMinFractions, 1)
      }
      cnProfile[ii, ] = c(tmp1, ccfCN[1], tmp2, ccfCN[2])
    } else {
      cnProfile[ii, 1:4] = c( GenerateCopyNumberProfile(cnPoolMaj, cnPoolMin,
                                                        cnPoolMajFractions, cnPoolMinFractions,
                                                        1), 1)
      cnProfile[ii, 8] = 0
    }
  }

  return(cnProfile)
}

registerDoParallel(cores=3)
set.seed(1234)

numSv <- 500
numOfClusterPool = 1:6
numOfRepeat = 1
baseDepth = 50

ccfCN <- c(0.7, 0.3)

ccfSet <- c(1, 0.3, 0.7) # true ccf pool
ccfTrue <- sample(ccfSet, numSv, c(0.5,0.3,0.2), replace = T)

purity <- 0.8
cnPoolMaj <- c(1,2,3,4)
cnPoolMin <- c(0,1,2)
cnPoolMajFractions <- c(0.25, 0.25, 0.25,0.25)
cnPoolMinFractions <- c(1/3, 1/3, 1/3)



# 1st break point

subclonal_cn1 <- sample(c(T, F), numSv, c(1,3), replace = T)
cnProfile = GenerateSubClonalCNProfile(cnPoolMaj, cnPoolMin,
                           cnPoolMajFractions, cnPoolMinFractions,
                           numSv, subclonal_cn1, ccfCN)

mydata <- data.frame(mutation_id = paste0("ss", seq_len(numSv)) ,
                     ccf_true = ccfTrue,
                     minor_cn1_sub1 = cnProfile[,1],
                     major_cn1_sub1 = cnProfile[,2],
                     total_cn1_sub1 = cnProfile[,3],
                     frac_cn1_sub1 = cnProfile[,4],
                     minor_cn1_sub2 = cnProfile[,5],
                     major_cn1_sub2 = cnProfile[,6],
                     total_cn1_sub2 = cnProfile[,7],
                     frac_cn1_sub2 = cnProfile[,8],
                     stringsAsFactors = F)

mydata$purity <- purity
mydata$normal_cn <- 2
mydata <- mutate(rowwise(mydata),
                 true_mult1_sub1 = sample(c(1,if (major_cn1_sub1 ==1) { 1 } else {major_cn1_sub1}), 1),
                 true_mult1_sub2 = if ( major_cn1_sub2 == -100 ) {
                   -100
                   } else {
                     sample(c(1,if (major_cn1_sub2 ==1) { 1 } else {major_cn1_sub2}), 1)
                   },
                 true_mult1 = frac_cn1_sub1 * true_mult1_sub1 + frac_cn1_sub2 * true_mult1_sub2,
                 total_cn1 = frac_cn1_sub1 * total_cn1_sub1 + frac_cn1_sub2 * total_cn1_sub2,
                 vaf1 = cp2ap(ccf_true, purity, normal_cn,
                              total_cn1,
                              total_cn1,
                              true_mult1),
                 total_counts1 = rpois(1, total_cn1/2 * baseDepth),
                 var_counts1 = rbinom(1, total_counts1, vaf1),
                 ref_counts1 = total_counts1 - var_counts1)

# 2nd break point
subclonal_cn2 <- sample(c(T, F), numSv, c(1,3), replace = T)
cnProfile = GenerateSubClonalCNProfile(cnPoolMaj, cnPoolMin,
                                       cnPoolMajFractions, cnPoolMinFractions,
                                       numSv, subclonal_cn2, ccfCN)

mydata$minor_cn2_sub1 = cnProfile[,1]
mydata$major_cn2_sub1 = cnProfile[,2]
mydata$total_cn2_sub1 = cnProfile[,3]
mydata$frac_cn2_sub1 = cnProfile[,4]
mydata$minor_cn2_sub2 = cnProfile[,5]
mydata$major_cn2_sub2 = cnProfile[,6]
mydata$total_cn2_sub2 = cnProfile[,7]
mydata$frac_cn2_sub2 = cnProfile[,8]

mydata <- mutate(rowwise(mydata),
                 true_mult2_sub1 = sample(c(1,if (major_cn2_sub1 ==1) { 1 } else {major_cn2_sub1}), 1),
                 true_mult2_sub2 = if ( major_cn2_sub2 == -100 ) {
                   -100
                 } else {
                   sample(c(1,if (major_cn2_sub2 ==1) { 1 } else {major_cn2_sub2}), 1)
                 },
                 true_mult2 = frac_cn2_sub1 * true_mult2_sub1 + frac_cn2_sub2 * true_mult2_sub2,
                 total_cn2 = frac_cn2_sub1 * total_cn2_sub1 + frac_cn2_sub2 * total_cn2_sub2,
                 vaf2 = cp2ap(ccf_true, purity, normal_cn,
                              total_cn2,
                              total_cn2,
                              true_mult2),
                 total_counts2 = rpois(1, total_cn2/2 * baseDepth),
                 var_counts2 = rbinom(1, total_counts2, vaf2),
                 ref_counts2 = total_counts2 - var_counts2
                 )

mydata$subclonal_cn1 = subclonal_cn1
mydata$subclonal_cn2 = subclonal_cn2


doubleBreakPtsRes <- RunCcubePipeline(dataFolder = "~/Dropbox/for_marek/", sampleName = "sv-test-sample",
                                      ssm = mydata, modelSV = T,
                                      numOfClusterPool = numOfClusterPool, numOfRepeat = numOfRepeat,
                                      runAnalysis = T, runQC = T,
                                      ccubeResultRDataFile = "~/Dropbox/for_marek/ccube_sv_subclonal_results.RData", multiCore = T,
                                      basicFormats = F, allFormats = F, returnAll = T)


fn1 = "~/Desktop/double_break_points_subclonal_results.pdf"
MakeCcubeStdPlot_sv(res = doubleBreakPtsRes$res, ssm = doubleBreakPtsRes$ssm, printPlot = T, fn = fn1)



mydata <- mutate(rowwise(mydata),
                 vaf1 = var_counts1/(var_counts1+ref_counts1),
                 true_obs_ccf1 = MapVaf2CcfPyClone(vaf1,
                                                   purity,
                                                   normal_cn,
                                                   total_cn1,
                                                   total_cn1,
                                                   true_mult1,
                                                   constraint=F),
                 vaf2 = var_counts2/(var_counts2+ref_counts2),
                 true_obs_ccf2 = MapVaf2CcfPyClone(vaf2,
                                                   purity,
                                                   normal_cn,
                                                   total_cn2,
                                                   total_cn2,
                                                   true_mult2,
                                                   constraint=F)
)

mydata$ccube_double_mult1 = doubleBreakPtsRes$res$full.model$bv1
mydata$ccube_double_mult2 = doubleBreakPtsRes$res$full.model$bv2

mydata <- mutate(rowwise(mydata),
              vaf1 = var_counts1/(var_counts1+ref_counts1),
              ccube_double_ccf1 = MapVaf2CcfPyClone(vaf1,
                                             purity,
                                             normal_cn,
                                             total_cn1,
                                             total_cn1,
                                             ccube_double_mult1,
                                             constraint=F),
              vaf2 = var_counts2/(var_counts2+ref_counts2),
              ccube_double_ccf2 = MapVaf2CcfPyClone(vaf2,
                                             purity,
                                             normal_cn,
                                             total_cn2,
                                             total_cn2,
                                             ccube_double_mult2,
                                             constraint=F),
              ccube_ccf1 = ccube_double_ccf1,
              ccube_ccf2 = ccube_double_ccf2
)


label1 = doubleBreakPtsRes$res$label

myColors=gg_color_hue(10)

fn = "~/Desktop/event_ccf_comparsions_subclonal_30_70.pdf"
pdf(fn, width=8, height=4)
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
dev.off()


fn = "~/Desktop/cluster_ccf_times_multiplicity_comparsions_subclonal_30_70.pdf"
pdf(fn, width=8, height=4)
par(mfrow=c(1,2))

mydata$ccube_ccf_mean <- doubleBreakPtsRes$res$full.model$ccfMean[doubleBreakPtsRes$res$label]
mydata$true_cluster_ccf1_mult1 = mydata$ccf_true*mydata$true_mult1
mydata$true_cluster_ccf2_mult2 = mydata$ccf_true*mydata$true_mult2
mydata$ccube_cluster_ccf1_mult1 = mydata$ccube_ccf_mean*mydata$ccube_double_mult1
mydata$ccube_cluster_ccf2_mult2 = mydata$ccube_ccf_mean*mydata$ccube_double_mult2


plot(mydata$true_cluster_ccf1_mult1, mydata$ccube_cluster_ccf1_mult1, col = myColors[label1],
     xlim = c(0, max( c(mydata$true_cluster_ccf1_mult1, mydata$ccube_cluster_ccf1_mult1) ) ),
     ylim = c(0, max( c(mydata$true_cluster_ccf1_mult1, mydata$ccube_cluster_ccf1_mult1) ) ),
     xlab = "true ccf cluster mean* true multiplicity", ylab = "estimated ccf cluster mean * estimated multiplicity",
     main = "double model: 1st break point")
points( seq(0, max( c(mydata$true_cluster_ccf1_mult1, mydata$ccube_cluster_ccf1_mult1) ), length.out = 100 ),
        seq(0, max( c(mydata$true_cluster_ccf1_mult1, mydata$ccube_cluster_ccf1_mult1) ), length.out = 100 ),
        type = "l" )

plot(mydata$true_cluster_ccf2_mult2, mydata$ccube_cluster_ccf2_mult2, col = myColors[label1],
     xlim = c(0, max( c(mydata$true_cluster_ccf2_mult2, mydata$ccube_cluster_ccf2_mult2) ) ),
     ylim = c(0, max( c(mydata$true_cluster_ccf2_mult2, mydata$ccube_cluster_ccf2_mult2) ) ),
     xlab = "true ccf cluster mean* true multiplicity", ylab = "estimated ccf cluster mean * estimated multiplicity", main = "double model: 2nd break point"
)

points( seq(0, max( c(mydata$true_cluster_ccf2_mult2, mydata$ccube_cluster_ccf2_mult2) ), length.out = 100 ),
        seq(0, max( c(mydata$true_cluster_ccf2_mult2, mydata$ccube_cluster_ccf2_mult2) ), length.out = 100 ),
        type = "l" )
dev.off()
