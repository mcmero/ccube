rm(list = ls())
library(dplyr)
library(ccube)
library(doParallel)
library(ggplot2)
library(tidyr)
library(gridExtra)
options(stringsAsFactors = F)

dirPath <- "~/OneDrive - University of Glasgow/Geoff_CopyNumber/test_subclonal_consensus//"
testCasesFolders <- dir(dirPath, full.names = T)
testCasesNames <- dir(dirPath)


numOfClusterPool = 1:7
numOfRepeat = 1

#for (ii in seq_along(testCasesFolders) ) {

  ii = which(testCasesNames =="804ffa2e-158b-447d-945c-707684134c87")

  cat (ii, "\n")
  testCasesFolder <- testCasesFolders[ii]
  testCasesName <- testCasesNames[ii]

  fnSNV <- dir(testCasesFolder, full.names = T, pattern = "snv")
  mydataSNV <- read.delim(fnSNV)
  mydataSNV <- mydataSNV[complete.cases(mydataSNV), ]

  mydataSNV$frac_cn_sub1 <- 1
  mydataSNV$frac_cn_sub2 <- 1 - mydataSNV$frac_cn_sub1

  mydataSNV1 <- ccube:::CheckAndPrepareCcubeInupts(mydataSNV)
  mydataSNV2 <- GetCcf(mydataSNV1, use = "use_base")


  res <- CcubeCore(mydata = mydataSNV2, init = max(numOfClusterPool), fit_mult = T,
                   use = "use_base", maxiter = 10, verbose = T)

  # resSNV <- RunCcubePipeline(ssm = mydataSNV, numOfClusterPool = numOfClusterPool, numOfRepeat = numOfRepeat,
  #                            runAnalysisSnap = T, runQC = T, maxiter = 100)
  #
  # fn1 <- paste0("~/Desktop/consensus_subclonal_cn_debug/", testCasesName, "snap_snv.pdf"  )
  # MakeCcubeStdPlot(res = resSNV$res, ssm = resSNV$ssm, printPlot = T, fn = fn1)
  #
  # fnSV <- dir(testCasesFolder, full.names = T, pattern = "sv")
  # mydataSV <- read.delim(fnSV)
  # mydataSV <- mydataSV[complete.cases(mydataSV), ]
  #
  #
  # resSV <- RunCcubePipeline(ssm = mydataSV, numOfClusterPool = numOfClusterPool,
  #                           numOfRepeat = numOfRepeat, modelSV = T,
  #                           runAnalysisSnap = T, runQC = T, maxiter = 100)
  #
  # fn1 <- paste0("~/Desktop/consensus_subclonal_cn_debug/", testCasesName, "snap_sv.pdf"  )
  #
  # MakeCcubeStdPlot_sv(res = resSV$res, ssm = resSV$ssm, printPlot = T, fn = fn1)
#}
