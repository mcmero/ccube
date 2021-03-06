rm(list = ls())
library(dplyr)
library(ccube)
library(doParallel)
library(ggplot2)
library(tidyr)
library(gridExtra)

registerDoParallel(cores=3)
set.seed(1234)

numSv <- 150
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
                 true_vaf1 = cp2ap(ccf_true, purity, normal_cn,
                              total_cn1,
                              total_cn1,
                              true_mult1),
                 total_counts1 = rpois(1, total_cn1/2 * baseDepth),
                 var_counts1 = rbinom(1, total_counts1, true_vaf1),
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
                 true_vaf2 = cp2ap(ccf_true, purity, normal_cn,
                              total_cn2,
                              total_cn2,
                              true_mult2),
                 total_counts2 = rpois(1, total_cn2/2 * baseDepth),
                 var_counts2 = rbinom(1, total_counts2, true_vaf2),
                 ref_counts2 = total_counts2 - var_counts2
                 )

mydata$subclonal_cn1 = subclonal_cn1
mydata$subclonal_cn2 = subclonal_cn2


doubleBreakPtsRes <- RunCcubePipeline(ssm = mydata, modelSV = T,
                                      numOfClusterPool = numOfClusterPool, numOfRepeat = numOfRepeat,
                                      runAnalysisSnap = T, runQC = T, maxiter = 100)


fn1 = "~/Desktop/double_break_points_subclonal_results_03_07_snap.pdf"
MakeCcubeStdPlot_sv(res = doubleBreakPtsRes$res, ssm = doubleBreakPtsRes$ssm, printPlot = T, fn = fn1)

mydata <- doubleBreakPtsRes$ssm

mydata <- mutate(rowwise(mydata),
                 true_obs_ccf1 = MapVaf2CcfPyClone(vaf1,
                                                   purity,
                                                   normal_cn,
                                                   total_cn1,
                                                   total_cn1,
                                                   true_mult1,
                                                   constraint=F),
                 true_obs_ccf2 = MapVaf2CcfPyClone(vaf2,
                                                   purity,
                                                   normal_cn,
                                                   total_cn2,
                                                   total_cn2,
                                                   true_mult2,
                                                   constraint=F)
)


label1 = doubleBreakPtsRes$res$label

myColors=gg_color_hue(10)

fn = "~/Desktop/event_ccf_comparsions_subclonal_30_70.pdf"
pdf(fn, width=8, height=4)
par(mfrow=c(1,2))
plot(mydata$true_obs_ccf1, mydata$ccube_ccf1, col = myColors[label1],
     xlim = c(0, max( c(mydata$true_obs_ccf1, mydata$ccube_ccf1) ) ),
     ylim = c(0, max( c(mydata$true_obs_ccf1, mydata$ccube_ccf1) ) ),
     xlab = "true ccf", ylab = "estimated ccf", main = "double model: 1st break point")
points( seq(0, max( c(mydata$true_obs_ccf1, mydata$ccube_ccf1) ), length.out = 100 ),
        seq(0, max( c(mydata$true_obs_ccf1, mydata$ccube_ccf1) ), length.out = 100 ),
        type = "l" )

plot(mydata$true_obs_ccf2, mydata$ccube_ccf2, col = myColors[label1],
     xlim = c(0, max( c(mydata$true_obs_ccf2, mydata$ccube_ccf2) ) ),
     ylim = c(0, max( c(mydata$true_obs_ccf2, mydata$ccube_ccf2) ) ),
     xlab = "true ccf", ylab = "estimated ccf", main = "double model: 2nd break point"
)

points( seq(0, max( c(mydata$true_obs_ccf2, mydata$ccube_ccf2) ), length.out = 100 ),
        seq(0, max( c(mydata$true_obs_ccf2, mydata$ccube_ccf2) ), length.out = 100 ),
        type = "l" )
dev.off()


fn = "~/Desktop/cluster_ccf_times_multiplicity_comparsions_subclonal_30_70.pdf"
pdf(fn, width=8, height=4)
par(mfrow=c(1,2))

mydata$true_cluster_ccf1_mult1 = mydata$ccf_true*mydata$true_mult1
mydata$true_cluster_ccf2_mult2 = mydata$ccf_true*mydata$true_mult2
mydata$ccube_cluster_ccf1_mult1 = mydata$ccube_ccf_mean*mydata$ccube_mult1
mydata$ccube_cluster_ccf2_mult2 = mydata$ccube_ccf_mean*mydata$ccube_mult2


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


mydata$error_mult1 =  mydata$ccube_mult1 - mydata$true_mult1
mydata$error_mult2 =  mydata$ccube_mult2 - mydata$true_mult2


selectedData <- mydata[, c("mutation_id","ccube_ccf_mean", "error_mult1", "error_mult2", "true_mult1", "true_mult2", "total_cn1", "total_cn2")]

fn = "~/Desktop/mults.pdf"
pdf(fn, width=8, height=8)

selectedData1 = gather(selectedData, key, value, -mutation_id, -total_cn1, -true_mult1, -ccube_ccf_mean)
tt1 = filter( selectedData1, key %in% c("error_mult1") )
g1 = ggplot(tt1, aes(y = value, x = as.factor(true_mult1), fill = as.factor(ccube_ccf_mean))) + geom_boxplot() +
  xlab("true_mult1") +  ylab("error") + theme(legend.position="none")

selectedData2 = gather(selectedData, key, value, -mutation_id, -total_cn2, -true_mult2, -ccube_ccf_mean)
tt2 = filter( selectedData2, key %in% c("error_mult2") )
g2 = ggplot(tt2, aes(y = value, x = as.factor(true_mult2), fill = as.factor(ccube_ccf_mean))) + geom_boxplot() + xlab("true_mult2") +
  ylab("error") + scale_fill_discrete(name="ccube cluster mean", labels= round(doubleBreakPtsRes$res$full.model$ccfMean,2) ) +
  theme(legend.position="bottom")

grid.arrange(g1, g2, nrow = 2)


dev.off()

