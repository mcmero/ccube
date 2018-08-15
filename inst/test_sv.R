rm(list = ls())
library(dplyr)
library(ccube)


set.seed(1234)

numSv <- 100
maxiter <- 200
init =5  # number of clusters
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
                 ref_counts2 = total_counts2 - var_counts2)



# SV prototype
res <- CcubeSVCore(mydata = mydata, init = init, fit_mult = T, use = "use_base", verbose = T)

label = res$label

mydata$ccube_mult1 <- res$full.model$bv1
mydata$ccube_mult2 <- res$full.model$bv2
mydata <- mutate(rowwise(mydata),
                  vaf1 = var_counts1/(var_counts1+ref_counts1),
                  ccube_ccf1 = MapVaf2CcfPyClone(vaf1,
                                            purity,
                                            normal_cn,
                                            major_cn1 + minor_cn1,
                                            major_cn1 + minor_cn1,
                                            ccube_mult1,
                                            constraint=F),
                  true_obs_ccf1 = MapVaf2CcfPyClone(vaf1,
                                                  purity,
                                                  normal_cn,
                                                  major_cn1 + minor_cn1,
                                                  major_cn1 + minor_cn1,
                                                  mult1,
                                                  constraint=F),
                 vaf2 = var_counts2/(var_counts2+ref_counts2),
                 ccube_ccf2 = MapVaf2CcfPyClone(vaf2,
                                                purity,
                                                normal_cn,
                                                major_cn2 + minor_cn2,
                                                major_cn2 + minor_cn2,
                                                ccube_mult2,
                                                constraint=F),
                 true_obs_ccf2 = MapVaf2CcfPyClone(vaf2,
                                                   purity,
                                                   normal_cn,
                                                   major_cn2 + minor_cn2,
                                                   major_cn2 + minor_cn2,
                                                   mult2,
                                                   constraint=F)
                 )
myColors=gg_color_hue(10)

fn1 = "~/Desktop/double_break_points_results.pdf"
pdf(fn1, width=8, height=8)

par(mfrow=c(2,2))
plot(mydata$true_obs_ccf1, mydata$ccube_ccf1, col = myColors[label],
     xlim = c(0, max( c(mydata$true_obs_ccf1, mydata$ccube_ccf1) ) ),
     ylim = c(0, max( c(mydata$true_obs_ccf1, mydata$ccube_ccf1) ) ),
     xlab = "true ccf", ylab = "estimated ccf", main = "1st break point")
points( seq(0, max( c(mydata$true_obs_ccf1, mydata$ccube_ccf1) ), length.out = 100 ),
          seq(0, max( c(mydata$true_obs_ccf1, mydata$ccube_ccf1) ), length.out = 100 ),
        type = "l" )

plot(mydata$true_obs_ccf2, mydata$ccube_ccf2, col = myColors[label],
     xlim = c(0, max( c(mydata$true_obs_ccf2, mydata$ccube_ccf2) ) ),
     ylim = c(0, max( c(mydata$true_obs_ccf2, mydata$ccube_ccf2) ) ),
     xlab = "true ccf", ylab = "estimated ccf", main = "2nd break point"
)

points( seq(0, max( c(mydata$true_obs_ccf2, mydata$ccube_ccf2) ), length.out = 100 ),
        seq(0, max( c(mydata$true_obs_ccf2, mydata$ccube_ccf2) ), length.out = 100 ),
        type = "l" )

tableSv <- table(label)
uniqLabels = unique(label)
names(tableSv) <- as.character(format(round(res$full.model$ccfMean[sort(uniqLabels)], 2), nsmall = 2))
barplot(tableSv, las = 2, col = myColors[sort(uniqLabels)],
        xlab = "cluster mean", ylab="number of variants",
        main = "cluster prevalence")

plot(res$L, col = myColors[5], type = "p", ylab = "ELBO", xlab = "iteration")

dev.off()

# Run ccube seperately

mydata1 = mydata[, c("mutation_id", "var_counts1", "ref_counts1", "major_cn1", "minor_cn1", "purity", "normal_cn")]
mydata1 = dplyr::rename(mydata1, var_counts = var_counts1, ref_counts = ref_counts1, major_cn = major_cn1, minor_cn = minor_cn1)
res1 <- CcubeCore(mydata = mydata1, init = 5, fit_mult = T, use = "use_one", verbose = T)
mydata1$ccube_ccf_mean <- res1$full.model$ccfMean[res1$label]
mydata1$ccube_mult <- res1$full.model$bv
mydata1 <- mutate(rowwise(mydata1),
              vaf = var_counts/(var_counts+ref_counts),
              ccube_ccf = MapVaf2CcfPyClone(vaf,
                                            purity,
                                            normal_cn,
                                            major_cn + minor_cn,
                                            major_cn + minor_cn,
                                            ccube_mult,
                                            constraint=F) )

fn2 = "~/Desktop/ccube_1st_breakpoint.pdf"
MakeCcubeStdPlot(mydata1, res1, printPlot = T, fn = fn2)

mydata2= mydata[, c("mutation_id", "var_counts2", "ref_counts2", "major_cn2", "minor_cn2", "purity", "normal_cn")]
mydata2 = dplyr::rename(mydata2, var_counts = var_counts2, ref_counts = ref_counts2, major_cn = major_cn2, minor_cn = minor_cn2)
res2 <- CcubeCore(mydata = mydata2, init = 5, fit_mult = T, use = "use_one", verbose = T)
mydata2$ccube_ccf_mean <- res2$full.model$ccfMean[res2$label]
mydata2$ccube_mult <- res2$full.model$bv
mydata2 <- mutate(rowwise(mydata2),
                  vaf = var_counts/(var_counts+ref_counts),
                  ccube_ccf = MapVaf2CcfPyClone(vaf,
                                                purity,
                                                normal_cn,
                                                major_cn + minor_cn,
                                                major_cn + minor_cn,
                                                ccube_mult,
                                                constraint=F) )
fn3 = "~/Desktop/ccube_2nd_breakpoint.pdf"
MakeCcubeStdPlot(mydata2, res2, printPlot = T, fn = fn3)
