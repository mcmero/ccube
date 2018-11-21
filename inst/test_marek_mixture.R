rm(list = ls())
library(dplyr)
library(ccube)
library(doParallel)
library(ggplot2)
library(tidyr)
library(gridExtra)

registerDoParallel(cores=3)


all_mixture_results=read.delim("~/Downloads/all_mixture_results.tsv", stringsAsFactors=FALSE)

mydata <- filter (all_mixture_results, mix %in% c('0.1-0.9'))

numOfClusterPool = 1:6
numOfRepeat = 1


#res = CcubeSVCore(mydata = mydata, init=1, fit_mult = T, fit_hyper = T, use = "use_base", verbose = F)

doubleBreakPtsRes <- RunCcubePipeline(dataFolder = "~/Dropbox/for_marek/", sampleName = "mixture_0.8_0.2",
                                      ssm = mydata, modelSV = T,
                                      numOfClusterPool = numOfClusterPool, numOfRepeat = numOfRepeat,
                                      runAnalysis = T, runQC = T,
                                      ccubeResultRDataFile = "~/Dropbox/for_marek/ccube_mixture_0.8_0.2_results.RData",
                                      multiCore = T,
                                      basicFormats = F, allFormats = F, returnAll = T)


fn1 = "~/mixture_0.8_0.2_results.pdf"
MakeCcubeStdPlot_sv(res = doubleBreakPtsRes$res, ssm = doubleBreakPtsRes$ssm, printPlot = T, fn = fn1)


