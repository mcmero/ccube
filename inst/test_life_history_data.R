rm(list = ls())
library(dplyr)
library(gtools)
library(doParallel)
library(foreach)
library(ccube)
registerDoParallel(cores=4)

# Test life history data

load("~/Dropbox/DirichletPrimer/data/life_history_snv.rdata")
library(colorspace)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
myColors <- gg_color_hue(10)

mydata <- filter(snv, chr!="X")



tt <- GetPurity(mydata)



res = ccube_m6(mydata, epi=1e-3,
         init=5, tol = 1e-10, maxiter = 1e3,
         fit_mult = T, fit_hyper = T, use = "use_base", verbose = T)


cat("\nEparam:", res$mu)
cat("\nEpi: ", res$full.model$Epi)
mydata$mean_ccf <- res$mu[res$label]
mydata$mult <- res$full.model$bv
#
mydata <- mutate(rowwise(mydata),
                 ccube_ccf = MapVaf2CcfPyClone(var_counts/(var_counts+ref_counts),
                                         purity,
                                         normal_cn,
                                         major_cn + minor_cn,
                                         major_cn + minor_cn,
                                         mult,
                                         constraint=F))

par(mfrow=c(2,2))
uniqLabels <- unique(res$label)
plot(mydata$ccube_ccf, mydata$vaf, col = myColors[res$label],
     xlab = "cancer cell fraction", ylab = "variant allele frequecy",
     main = "ccf vs vaf (colored by cluster memebership)")
mydata$total_cn =mydata$major_cn+mydata$minor_cn
uniqueTotCn = unique(mydata$total_cn)
xx = seq(0,2, length.out = 100)
for (cn in uniqueTotCn) {
  for (i in 1:cn) {
    points(MapVaf2CcfPyClone(xx, unique(mydata$purity), 2, cn, cn, i, constraint = F), xx, type = 'l')
  }
}


Emu <- res$full.model$ccfMean
Esigma <- res$full.model$ccfCov
Epi <- res$full.model$Epi
params <- data.frame(Emu, Esigma, Epi)
xx <- seq(0, 2,  length.out = 1000)
ll <- 0
ll1 <- 0

for (j in seq_len(nrow(params))) {
  ll <- ll + params[j,]$Epi * dnorm(xx, mean = params[j,]$Emu, sd = sqrt(params[j,]$Esigma))
}

hist(mydata$ccube_ccf, density=20, breaks=20, prob=TRUE,
     main = "ccf histogram +
     fitted marginal denstiy (red)",
     xlab = "cancer cell fraction")
lines(xx,ll, lwd=2, col = "darkred")

names(Epi) <- as.character(format(round(Emu, 2), nsmall = 2))
barplot(Epi[sort(uniqLabels)], las = 2, col = myColors[sort(uniqLabels)],
        xlab = "cluster mean", ylab="expected weight",
        main = "cluster weights")
