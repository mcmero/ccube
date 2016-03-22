#!/home/yuan03/R-3.2.0/bin/Rscript
library(doMPI)

args <- commandArgs(trailingOnly = TRUE)
numofCores <- as.integer(args[1]) - 1
print(numofCores)
cl <- startMPIcluster(count = numofCores)
registerDoMPI(cl)
library(dplyr)
library(foreach)

# data path
icgc <- "/lustre/fmlab/yuan03/PCAWG/db_analysis/"
#icgc <- '~/Desktop/LustreCI/PCAWG/bb_batch08/'   

i <- as.integer(args[2])
caller <- as.character(args[3])
jobList <- as.character(args[4])
muse_tier <- as.integer(args[5])

ascatFiles <- dir(paste0(icgc, "ascat_calls/"))
ascatFiles1 <- dir(paste0(icgc, "ascat_calls_1/"))
sampleNameAscat <- Reduce(c, Map(function(x) x[1], strsplit(ascatFiles, '_', fixed=TRUE)))
sampleNameAscat1 <- Reduce(c, Map(function(x) x[1], strsplit(ascatFiles1, '_', fixed=TRUE))) 

ascatFiles <- dir(paste0(icgc, "ascat_calls/"), full.names = T)
ascatFiles1 <- dir(paste0(icgc, "ascat_calls_1/"), full.names = T)

if (caller == "dkfz") {
  vcfFiles <- dir(paste0(icgc, "dkfz_snv_vcf/"))
  sampleNameVcf <- Reduce(c, Map(function(x) x[1], strsplit(vcfFiles, '.', fixed=TRUE)), c())
  dkfzSampleIndexWithBB <- which(sampleNameVcf %in% sampleNameAscat)
  sampleNameVcf <- sampleNameVcf[dkfzSampleIndexWithBB]
} else if (caller == "mutect"){
  vcfFiles <- dir(paste0(icgc, "mutect_snv_vcf/"))
  sampleNameVcf <- Reduce(c, Map(function(x) x[1], strsplit(vcfFiles, '.', fixed=TRUE)), c())
  mutectSampleIndexWithBB <- which(sampleNameVcf %in% sampleNameAscat)
  sampleNameVcf <- sampleNameVcf[mutectSampleIndexWithBB]
} else if (caller == "mutect_smchet"){
  vcfFiles <- dir(paste0(icgc, "mutect_smchet_snv_vcf/"))
  sampleNameVcf <- Reduce(c, Map(function(x) x[1], strsplit(vcfFiles, '.', fixed=TRUE)), c())
  mutectSampleIndexWithBB <- which(sampleNameVcf %in% sampleNameAscat)
  sampleNameVcf <- sampleNameVcf[mutectSampleIndexWithBB]
} else if (caller == "muse"){
  vcfFiles <- dir(paste0(icgc, "muse_snv_vcf/"))
  sampleNameVcf <- Reduce(c, Map(function(x) x[1], strsplit(vcfFiles, '.', fixed=TRUE)), c())
  museSampleIndexWithBB <- which(sampleNameVcf %in% sampleNameAscat)
  sampleNameVcf <- sampleNameVcf[museSampleIndexWithBB]
} else {
  vcfFiles <- dir(paste0(icgc, "sanger_snv_vcf/"))
  sampleNameVcf <- Reduce(c, Map(function(x) x[1], strsplit(vcfFiles, '.', fixed=TRUE)), c())
  sangerSampleIndexWithBB <- which(sampleNameVcf %in% sampleNameAscat)
  sampleNameVcf <- sampleNameVcf[sangerSampleIndexWithBB]
}

if (is.na(jobList)) {
  sampleList <- sampleNameVcf
} else {
  sampleList <- read.csv(jobList, sep="", stringsAsFactors=FALSE)$x
}
sampleName <- sampleList[i]

if (sampleName %in% sampleNameAscat) {
  ascatFile <- ascatFiles[which(sampleNameAscat == sampleName)]
  cellularity <- read.table(ascatFile, header = T)$cellularity
  if (is.null(cellularity)) {
    cellularity <- read.table(ascatFile, header = T)$purity
  }
} else {
  ascatFile <- ascatFiles1[which(sampleNameAscat1 == sampleName)]
  cellularity <- read.table(ascatFile, header = T)$rho[2]
}

if (is.null(cellularity)) {
  cellularity <- read.table(ascatFile, stringsAsFactors=FALSE)$V1
}

dataFolders <- dir(paste0(icgc, "test_snv/"), full.names = T)
dataFolderNames <- dir(paste0(icgc, "test_snv/"))
dataFolder <- dataFolders[which(dataFolderNames == sampleName)]

if (caller == "muse") {
  if (is.na(muse_tier)) {
    muse_tier <- 5
  }
  load( paste0(dataFolder, "/ssm_no_chrxy_no_subclonal_cn_", 
               caller, "_tier_", muse_tier, ".RData") )
} else if (caller == "dkfz") {
  load(paste0(dataFolder, "/ssm_no_chrxy_no_subclonal_cn_dkfz.RData"))
} else if (caller == "mutect"){
  load(paste0(dataFolder, "/ssm_no_chrxy_no_subclonal_cn_mutect.RData"))
} else if (caller == "mutect_smchet"){
  load(paste0(dataFolder, "/ssm_no_chrxy_no_subclonal_cn_mutect_smchet.RData"))
} else {
  load(paste0(dataFolder, "/ssm_no_chrxy_no_subclonal_cn.RData"))
}


library(ccube)

ssm$purity <- GetPurity(ssm)
kk <- 6
if (kk > nrow(ssm)){
  kk <- nrow(ssm) - 1
}
rr <- 10
iterSetting <- data.frame( sort(rep(seq(1, kk, length.out = kk), rr)) )
  
results <- foreach(n = 1:nrow(iterSetting), .combine = c) %dopar% {
  library(ccube)
  k <- iterSetting[n, ]
  cat(k)
  list(ccube_m6(ssm, epi=1e-3,
         init=k, tol = 1e-10, maxiter = 1e3,
         fit_mult = T, fit_hyper = T, use = "use_base", verbose = F))
}

maxLbIndex <- which.max(Map( function(x) max(x$L), results))
lb <- unlist(Map( function(x) max(x$L), results))
res <- results[[maxLbIndex]]
  
if (caller == "muse") {
  if (is.na(muse_tier)) {
    muse_tier <- 0
  }
  save(ssm, ccf, results, res, lb, 
       file = paste0(dataFolder,
                     "/smm_map_cluster_res_no_chrxy_no_subclonal_cn_1_", 
                     caller, "_tier_", muse_tier, ".RData"))
} else if (caller == "dkfz") {
  save(ssm, results, res, lb, file = paste0(dataFolder,"/smm_map_cluster_res_no_chrxy_no_subclonal_cn_1_dkfz.RData"))
} else if (caller == "mutect"){
  save(ssm, results, res, lb, file = paste0(dataFolder,"/smm_map_cluster_res_no_chrxy_no_subclonal_cn_1_mutect.RData"))
} else if (caller == "mutect_smchet"){
  save(ssm, results, res, lb, file = paste0(dataFolder,"/smm_map_cluster_res_no_chrxy_no_subclonal_cn_1_mutect_smchet.RData"))
} else {
  save(ssm, results, res, lb, file = paste0(dataFolder,"/smm_map_cluster_res_no_chrxy_no_subclonal_cn_1.RData"))
}
cat("\n", "finished number", i, "in", length(sampleList), "\n", 
    file = paste0("progress_run_smm_", caller), append = T)

closeCluster(cl)
mpi.quit()
