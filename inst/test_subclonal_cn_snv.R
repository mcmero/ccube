rm(list = ls())
library(dplyr)
library(ccube)
library(doParallel)
library(ggplot2)
library(tidyr)
library(gridExtra)

registerDoParallel(cores=3)
set.seed(1234)

numSv <- 500
numOfClusterPool = 1:6
numOfRepeat = 1
baseDepth = 50

ccfCN <- c(0.7, 0.3)

ccfSet <- c(1, 0.3, 0.7) # true ccf pool
ccfTrue <- sample(ccfSet, numSv, c(0.5,0.2,0.3), replace = T)

purity <- 0.8
cnPoolMaj <- c(1,2,3,4)
cnPoolMin <- c(0,1,2)
cnPoolMajFractions <- c(0.25, 0.25, 0.25,0.25)
cnPoolMinFractions <- c(1/3, 1/3, 1/3)

subclonal_cn <- sample(c(T, F), numSv, c(1,5), replace = T)
cnProfile = GenerateSubClonalCNProfile(cnPoolMaj, cnPoolMin,
                           cnPoolMajFractions, cnPoolMinFractions,
                           numSv, subclonal_cn, ccfCN)

mydata <- data.frame(mutation_id = paste0("ss", seq_len(numSv)) ,
                     ccf_true = ccfTrue,
                     minor_cn_sub1 = cnProfile[,1],
                     major_cn_sub1 = cnProfile[,2],
                     total_cn_sub1 = cnProfile[,3],
                     frac_cn_sub1 = cnProfile[,4],
                     minor_cn_sub2 = cnProfile[,5],
                     major_cn_sub2 = cnProfile[,6],
                     total_cn_sub2 = cnProfile[,7],
                     frac_cn_sub2 = cnProfile[,8],
                     stringsAsFactors = F)

mydata$purity <- purity
mydata$normal_cn <- 2
mydata$subclonal_cn <- subclonal_cn
mydata <- mutate(rowwise(mydata),
                 true_mult_sub1 = sample(c(1,if (major_cn_sub1 ==1) { 1 } else {major_cn_sub1}), 1),
                 true_mult_sub2 = if ( major_cn_sub2 == -100 ) {
                   -100
                   } else {
                     sample(c(1,if (major_cn_sub2 ==1) { 1 } else {major_cn_sub2}), 1)
                   },
                 true_mult = frac_cn_sub1 * true_mult_sub1 + frac_cn_sub2 * true_mult_sub2,
                 total_cn = frac_cn_sub1 * total_cn_sub1 + frac_cn_sub2 * total_cn_sub2,
                 vaf = cp2ap(ccf_true, purity, normal_cn,
                              total_cn,
                              total_cn,
                              true_mult),
                 total_counts = rpois(1, total_cn/2 * baseDepth),
                 var_counts = rbinom(1, total_counts, vaf),
                 ref_counts = total_counts - var_counts)


ccubeRes <- RunCcubePipeline(ssm = mydata, numOfClusterPool = numOfClusterPool, numOfRepeat = numOfRepeat,
                             runAnalysis = T, runQC = T, multiCore = T,
                             basicFormats = F, allFormats = F, returnAll = T)


fn1 = "~/Desktop/snv_subclonal_results.pdf"
MakeCcubeStdPlot(res = ccubeRes$res, ssm = ccubeRes$ssm, printPlot = T, fn = fn1)

mydata = ccubeRes$ssm
mydata <- dplyr::mutate(rowwise(mydata),
                 true_obs_ccf = MapVaf2CcfPyClone(vaf,
                                                   purity,
                                                   normal_cn,
                                                   total_cn,
                                                   total_cn,
                                                   true_mult,
                                                   constraint=F)
)

label = ccubeRes$res$label
myColors=gg_color_hue(10)

fn = "~/Desktop/snv_event_ccf_comparsions_subclonal_30_70.pdf"
pdf(fn, width=8, height=4)
par(mfrow=c(1,2))
plot(mydata$true_obs_ccf, mydata$ccube_ccf, col = myColors[label],
     xlim = c(0, max( c(mydata$true_obs_ccf, mydata$ccube_ccf) ) ),
     ylim = c(0, max( c(mydata$true_obs_ccf, mydata$ccube_ccf) ) ),
     xlab = "true ccf", ylab = "estimated ccf", main = "SNV model")
points( seq(0, max( c(mydata$true_obs_ccf, mydata$ccube_ccf) ), length.out = 100 ),
        seq(0, max( c(mydata$true_obs_ccf, mydata$ccube_ccf) ), length.out = 100 ),
        type = "l" )

mydata$true_cluster_ccf_mult = mydata$ccf_true*mydata$true_mult
mydata$ccube_cluster_ccf_mult = mydata$ccube_ccf_mean*mydata$ccube_mult

plot(mydata$true_cluster_ccf_mult, mydata$ccube_cluster_ccf_mult, col = myColors[label],
     xlim = c(0, max( c(mydata$true_cluster_ccf_mult, mydata$ccube_cluster_ccf_mult) ) ),
     ylim = c(0, max( c(mydata$true_cluster_ccf_mult, mydata$ccube_cluster_ccf_mult) ) ),
     xlab = "true ccf cluster mean * \n true multiplicity", ylab = "estimated ccf cluster mean *  \n estimated multiplicity",
     main = "SNV model")
points( seq(0, max( c(mydata$true_cluster_ccf_mult, mydata$ccube_cluster_ccf_mult) ), length.out = 100 ),
        seq(0, max( c(mydata$true_cluster_ccf_mult, mydata$ccube_cluster_ccf_mult) ), length.out = 100 ),
        type = "l" )

dev.off()

mydata$error_mult =  mydata$ccube_mult - mydata$true_mult
selectedData <- mydata[, c("mutation_id","ccube_ccf_mean", "true_mult", "total_cn", "error_mult")]
fn = "~/Desktop/snv_mults.pdf"
pdf(fn, width=8, height=4)

selectedData1 = gather(selectedData, key, value, -mutation_id, -total_cn, -true_mult, -ccube_ccf_mean)
tt1 = filter( selectedData1, key %in% c("error_mult") )
g1 = ggplot(tt1, aes(y = value, x = as.factor(true_mult), fill = as.factor(ccube_ccf_mean))) + geom_boxplot() +
  xlab("true_mult") +  ylab("error") + theme(legend.position="none")
print(g1)
dev.off()

