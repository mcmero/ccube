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

allDebugFolder <- "~/debug_pcawg_samples_sv/debuged_samples"

sampleName <- "abd2d959-d5ed-4eb3-9759-67eb1aa23325"

resultsFolder <- paste0(allDebugFolder,"/", sampleName)

if (! dir.exists(resultsFolder) ) {
  dir.create(resultsFolder, recursive = T)
}

inputFn <- "~/OneDrive - University of Glasgow/Geoff_CopyNumber/abd2d959-d5ed-4eb3-9759-67eb1aa23325_ccube_sv_input.txt"

mydata <- try( read.delim(inputFn), T)

if ( is.data.frame(mydata) ) {

  if (nrow(mydata) == 1) {
    singleEventSamples <- c(singleEventSamples, sampleName)
    cat("single event", file = paste0(resultsFolder, "/single_event_sample")  )
  }

  ccubeRes <-
    RunCcubePipeline(ssm = mydata, modelSV = T,
                     numOfClusterPool = numOfClusterPool,
                     numOfRepeat = numOfRepeat, multiCore =T,
                     runAnalysis = T, runQC = T)

  if (is.list (ccubeRes) ) {
    fn = paste0(resultsFolder, "/ccube_sv_results.RDdata")
    save(ccubeRes, file = fn)
    fn = paste0(resultsFolder, "/ccube_sv_results.pdf")
    MakeCcubeStdPlot_sv(res = ccubeRes$res, ssm = ccubeRes$ssm, printPlot = T, fn = fn)
  } else {
    cat(ccubeRes, file = paste0(resultsFolder, "/bug_info_ccube.txt"))

    if (nrow(mydata) > 1) {
      problemSamples <- c(problemSamples, sampleName)
      cat("problematic sample", file = paste0(resultsFolder, "/problematic_sample")  )
    }

  }

} else {
  cat(mydata, file = paste0(resultsFolder, "/bug_info_mydata.txt"))

}



