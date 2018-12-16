rm(list = ls())
library(dplyr)
library(ccube)
library(doParallel)
library(ggplot2)
library(tidyr)
library(gridExtra)

options(stringsAsFactors = F)

registerDoParallel(cores=3)
set.seed(1234)


numOfClusterPool = 1:7
numOfRepeat = 10


icgc <- "~/debug_pcawg_samples_sv/svclone_paper_20181127/"
dataFolders <- dir(icgc, full.names = T)
sampleNames <- dir(icgc)


bugSamples <- c()
bugdataFolders <- c()
for (ii in seq_along(sampleNames)) {

  sampleName <- sampleNames[ii]
  dataFolder <- dataFolders[ii]

  if (! "ccube_out/ccube_sv_results.RData" %in% dir(dataFolder, recursive = T) ) {
    bugSamples <- c(bugSamples, sampleName)
    bugdataFolders <- c(bugdataFolders, dataFolder)
  }
}


allDebugFolder <- "~/debug_pcawg_samples_sv/debuged_samples"
problemSamples <- c("36e1d9cc-32ec-4a0a-8fb1-c46f058a6fb8", "60413de1-6cd2-4f74-8180-3bdd394d6d16",
                   "7dc5f8ba-0080-43d3-8426-bd527a970761", "7eac4710-c622-11e3-bf01-24c6515278c0",
                    "ea43434b-197e-48ac-ae2e-46bc7f3776de", "fc249113-83d4-4abe-8c80-a4f7305dcd91")

  ii = which(sampleNames == problemSamples[6])
  cat(ii, "\n")

  sampleName <- sampleNames[ii]
  dataFolder <- dataFolders[ii]
  resultsFolder <- paste0(allDebugFolder,"/", sampleName)

  if (! dir.exists(resultsFolder) ) {
    dir.create(resultsFolder, recursive = T)
  }

  inputFn <- paste0(dataFolder, "/", sampleName, "_ccube_sv_input.txt")

  mydata <- try( read.delim(inputFn), T)



  ccubeRes <-
    RunCcubePipeline(ssm = mydata, modelSV = T,
                     numOfClusterPool = numOfClusterPool,
                     numOfRepeat = numOfRepeat, multiCore =T,
                     runAnalysisSnap = T, runQC = T, maxiter = 10)


    fn = paste0(resultsFolder, "/ccube_sv_results.RDdata")
    save(ccubeRes, file = fn)
    fn = paste0(resultsFolder, "/ccube_sv_results.pdf")
    MakeCcubeStdPlot_sv(res = ccubeRes$res, ssm = ccubeRes$ssm, printPlot = T, fn = fn)




