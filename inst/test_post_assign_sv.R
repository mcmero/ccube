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

# Post-assign with combined results
postAssignRes <- RunPostAssignPipeline(snvRes = snvRes1$res, svRes = svRes1$res, mydata = mydata)
fn = "~/Desktop/bad_sample_post_assign_sample_combined_model_based.pdf"
MakeCcubeStdPlot_sv(ssm = postAssignRes$ssm, res = postAssignRes$res, printPlot = T,fn=fn)

# Post-assign with SNV results only
postAssignRes <- RunPostAssignPipeline(snvRes = snvRes1$res, mydata = mydata)
fn = "~/Desktop/bad_sample_post_assign_sample_snv_model_based.pdf"
MakeCcubeStdPlot_sv(ssm = postAssignRes$ssm, res = postAssignRes$res, printPlot = T,fn=fn)


# Test good sample
samplePath <- "~/OneDrive - University of Glasgow/Geoff_CopyNumber/post_assign/good_post_assign/"
sampleSNVRdata <- dir(samplePath, pattern = "ccube_snv_results.RData", full.names = T)
sampleSVRdata <- dir(samplePath, pattern = "ccube_sv_results.RData", full.names = T)

load(sampleSNVRdata)
load(sampleSVRdata)

svRes2 <- doubleBreakPtsRes
snvRes2 <- snvRes
mydata <- svRes2$ssm

# Post-assign with combined results
postAssignRes <- RunPostAssignPipeline(snvRes = snvRes2$res, svRes = svRes2$res, mydata = mydata)
fn = "~/Desktop/good_sample_post_assign_combined_model_based.pdf"
MakeCcubeStdPlot_sv(ssm = postAssignRes$ssm, res = postAssignRes$res, printPlot = T,fn=fn)

# Post-assign with snv results only
postAssignRes <- RunPostAssignPipeline(snvRes = snvRes2$res, mydata = mydata)
fn = "~/Desktop/good_sample_post_assign_snv_model_based.pdf"
MakeCcubeStdPlot_sv(ssm = postAssignRes$ssm, res = postAssignRes$res, printPlot = T,fn=fn)
