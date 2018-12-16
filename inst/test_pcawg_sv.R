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
numOfRepeat = 1


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


allDebugFolder <- "~/debug_pcawg_samples_sv/debuged_samples_all_full_run_low_iter_10_test"

if (! dir.exists(allDebugFolder) ) {
  dir.create(allDebugFolder)
}


singleEventSamples <- c()
problemSamples <- c()
inputProblemSamples <- c()
for (ii in seq_along(sampleNames)) {

  cat(ii, "\n")

  sampleName <- sampleNames[ii]
  dataFolder <- dataFolders[ii]
  resultsFolder <- paste0(allDebugFolder,"/", sampleName)

  if (! dir.exists(resultsFolder) ) {
    dir.create(resultsFolder)
  }

  inputFn <- paste0(dataFolder, "/", sampleName, "_ccube_sv_input.txt")

  mydata <- try( read.delim(inputFn), T)





  if ( is.data.frame(mydata) ) {

    if (nrow(mydata) == 1) {
      singleEventSamples <- c(singleEventSamples, sampleName)
      cat("single event", file = paste0(resultsFolder, "/single_event_sample")  )
    }

    ccubeRes <- try(
      RunCcubePipeline(ssm = mydata, modelSV = T,
                       numOfClusterPool = numOfClusterPool,
                       numOfRepeat = numOfRepeat, multiCore =T,
                       runAnalysis = T, runQC = T, maxiter = 10),
      T
    )

    if (is.list (ccubeRes) ) {
      fn = paste0(resultsFolder, "/ccube_sv_results.RDdata")
      save(ccubeRes, file = fn)

      if (nrow(mydata) > 1) {
        fn = paste0(resultsFolder, "/ccube_sv_results.pdf")
        MakeCcubeStdPlot_sv(res = ccubeRes$res, ssm = ccubeRes$ssm, printPlot = T, fn = fn)
      }

    } else {
      cat(ccubeRes, file = paste0(resultsFolder, "/bug_info_ccube.txt"))

      if (nrow(mydata) > 1) {
        problemSamples <- c(problemSamples, sampleName)
        cat("problematic sample", file = paste0(resultsFolder, "/problematic_sample")  )
      }

    }

  } else {
    inputProblemSamples <- c(inputProblemSamples, sampleName)
    cat(mydata, file = paste0(resultsFolder, "/bug_info_mydata.txt"))

  }

}



