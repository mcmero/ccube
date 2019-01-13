rm(list = ls())
library(dplyr)
library(ccube)
library(doParallel)
library(ggplot2)
library(tidyr)
library(gridExtra)


# Test bad sample

samplePath <- "~/OneDrive - University of Glasgow/Geoff_CopyNumber/post_assign/bad_post_assign/"
sampleSNVRdata <- dir(samplePath, pattern = "ccube_snv_results.RData", full.names = T)
sampleSVRdata <- dir(samplePath, pattern = "ccube_sv_results.RData", full.names = T)

load(sampleSNVRdata)
load(sampleSVRdata)

svRes1 <- doubleBreakPtsRes
snvRes1 <- snvRes

mydata <- svRes1$ssm
referenceRes <- snvRes1$res
postAssignRes <- AssignWithCcube_sv(referenceRes, mydata, verbose = T)
annotatedSsm <- AnnotateCcubeResults_sv(mydata, postAssignRes)

fn = "~/Desktop/bad_post_assign_sample_model_based.pdf"
MakeCcubeStdPlot_sv(ssm = annotatedSsm, res = postAssignRes, printPlot = T,fn=fn)



# Test good sample

samplePath <- "~/OneDrive - University of Glasgow/Geoff_CopyNumber/post_assign/good_post_assign/"
sampleSNVRdata <- dir(samplePath, pattern = "ccube_snv_results.RData", full.names = T)
sampleSVRdata <- dir(samplePath, pattern = "ccube_sv_results.RData", full.names = T)

load(sampleSNVRdata)
load(sampleSVRdata)

svRes2 <- doubleBreakPtsRes
snvRes2 <- snvRes

mydata <- svRes2$ssm
referenceRes <- snvRes2$res
postAssignRes <- AssignWithCcube_sv(referenceRes, mydata, verbose = T)
annotatedSsm <- AnnotateCcubeResults_sv(mydata, postAssignRes)

fn = "~/Desktop/good_post_assign_sample_model_based.pdf"
MakeCcubeStdPlot_sv(ssm = annotatedSsm, res = postAssignRes, printPlot = T,fn=fn)
