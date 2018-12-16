rm(list = ls())
library(dplyr)
library(ccube)
library(doParallel)

registerDoParallel(cores=3)

set.seed(1234)
purity <-0.98
numOfClusterPool = 1:6
numOfRepeat = 1
mydata <- read.delim("~/Dropbox/for_adam/G-R1-P_pyclone_input_350reads.tsv", stringsAsFactors = F)

mydata$purity <- purity

results <- RunCcubePipeline(ssm = mydata, numOfClusterPool = numOfClusterPool,
                            numOfRepeat = numOfRepeat,runAnalysis = T, runQC = T,
                            multiCore = T)



fn = "~/Dropbox/for_adam/G-R1-P_ccube.pdf"
MakeCcubeStdPlot(ssm = results$ssm, res = results$res, printPlot = T, fn = fn)
