#!/home/yuan03/R-3.2.0/bin/Rscript
library(dplyr)
icgc <- '/lustre/fmlab/yuan03/PCAWG/db_analysis/'
#icgc <- "~/Desktop/LustreCI//dream_challenge/"

args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])
caller <- as.character(args[2])
jobList <- as.character(args[3])
muse_tier <- as.integer(args[4])

# Get vcf with battenberg calls
ascatFiles <- dir(paste0(icgc, "ascat_calls/"))
ascatFiles1 <- dir(paste0(icgc, "ascat_calls_1/"))
batternbergFiles <- dir(paste0(icgc, "battenberg_calls/"))

sampleNameBB <- Reduce(c, Map(function(x) x[1], strsplit(batternbergFiles, '_', fixed=TRUE)), c())
sampleNameAscat <- Reduce(c, Map(function(x) x[1], strsplit(ascatFiles, '_', fixed=TRUE)), c())
sampleNameAscat1 <- Reduce(c, Map(function(x) x[1], strsplit(ascatFiles1, '_', fixed=TRUE)), c()) 
allAscat <- union(sampleNameAscat, sampleNameAscat1)
sampleBBAscatUnion <- intersect(sampleNameBB, allAscat)

if (caller == "dkfz") {
  vcfFiles <- dir(paste0(icgc, "dkfz_snv_vcf/"))
  sampleNameVcf <- Reduce(c, Map(function(x) x[1], strsplit(vcfFiles, '.', fixed=TRUE)), c())
  dkfzSampleIndexWithBB <- which(sampleNameVcf %in% sampleBBAscatUnion)
  sampleNameVcf <- sampleNameVcf[dkfzSampleIndexWithBB]
  vcfFiles <- dir(paste0(icgc, "dkfz_snv_vcf/"), full.names = T)
  vcfFiles <- vcfFiles[dkfzSampleIndexWithBB]
} else if (caller == "mutect") {
  vcfFiles <- dir(paste0(icgc, "mutect_snv_vcf/"))
  sampleNameVcf <- Reduce(c, Map(function(x) x[1], strsplit(vcfFiles, '.', fixed=TRUE)), c())
  mutectSampleIndexWithBB <- which(sampleNameVcf %in% sampleBBAscatUnion)
  sampleNameVcf <- sampleNameVcf[mutectSampleIndexWithBB]
  vcfFiles <- dir(paste0(icgc, "mutect_snv_vcf/"), full.names = T)
  vcfFiles <- vcfFiles[mutectSampleIndexWithBB]
} else if (caller == "mutect_smchet") {
  vcfFiles <- dir(paste0(icgc, "mutect_smchet_snv_vcf/"))
  sampleNameVcf <- Reduce(c, Map(function(x) x[1], strsplit(vcfFiles, '.', fixed=TRUE)), c())
  mutectSampleIndexWithBB <- which(sampleNameVcf %in% sampleBBAscatUnion)
  sampleNameVcf <- sampleNameVcf[mutectSampleIndexWithBB]
  vcfFiles <- dir(paste0(icgc, "mutect_smchet_snv_vcf/"), full.names = T)
  vcfFiles <- vcfFiles[mutectSampleIndexWithBB]
} else if (caller == "muse") {
  vcfFiles <- dir(paste0(icgc, "muse_snv_vcf/"))
  sampleNameVcf <- Reduce(c, Map(function(x) x[1], strsplit(vcfFiles, '.', fixed=TRUE)), c())
  museSampleIndexWithBB <- which(sampleNameVcf %in% sampleBBAscatUnion)
  sampleNameVcf <- sampleNameVcf[museSampleIndexWithBB]
  vcfFiles <- dir(paste0(icgc, "broad_muse_snv_vcf/"), full.names = T)
  vcfFiles <- vcfFiles[museSampleIndexWithBB]
} else {
  vcfFiles <- dir(paste0(icgc, "sanger_snv_vcf/"))
  sampleNameVcf <- Reduce(c, Map(function(x) x[1], strsplit(vcfFiles, '.', fixed=TRUE)), c())
  sangerSampleIndexWithBB <- which(sampleNameVcf %in% sampleBBAscatUnion)
  sampleNameVcf <- sampleNameVcf[sangerSampleIndexWithBB]
  vcfFiles <- dir(paste0(icgc, "sanger_snv_vcf/"), full.names = T)
  vcfFiles <- vcfFiles[sangerSampleIndexWithBB]
}

if (is.na(jobList)) {
  sampleList <- sampleNameVcf
} else {
  sampleList <- read.csv(jobList, sep="", stringsAsFactors=FALSE)$x
}
sampleName <- sampleList[i]

dataFolders <- dir(paste0(icgc, "test_snv/"), full.names = T)
dataFolderNames <- dir(paste0(icgc, "test_snv/"))
dataFolder <- dataFolders[which(dataFolderNames == sampleName)]

if (caller == "muse") {
  if (is.na(muse_tier)) {
    muse_tier = 5
  }
  ssm <- read.delim(paste0(dataFolder, "/ssm_data_no_chrxy_with_frac_", 
                           caller, "_tier_", muse_tier, ".txt"), 
                    stringsAsFactors = F)
  cnv <- read.delim(paste0(dataFolder, "/cnv_data_no_chrxy_with_frac_", 
                           caller, "_tier_", muse_tier, ".txt"), 
                    stringsAsFactors = F) 
} else if (caller == "dkfz") {
  ssm <- read.delim(paste0(dataFolder, "/ssm_data_no_chrxy_with_frac_dkfz.txt"), 
                    stringsAsFactors = F)
  cnv <- read.delim(paste0(dataFolder, "/cnv_data_no_chrxy_with_frac_dkfz.txt"), 
                    stringsAsFactors = F) 
} else if (caller == "mutect"){
  ssm <- read.delim(paste0(dataFolder, "/ssm_data_no_chrxy_with_frac_mutect.txt"), 
                    stringsAsFactors = F)
  cnv <- read.delim(paste0(dataFolder, "/cnv_data_no_chrxy_with_frac_mutect.txt"), 
                    stringsAsFactors = F) 
} else if (caller == "mutect_smchet"){
  ssm <- read.delim(paste0(dataFolder, "/ssm_data_no_chrxy_with_frac_mutect_smchet.txt"), 
                    stringsAsFactors = F)
  cnv <- read.delim(paste0(dataFolder, "/cnv_data_no_chrxy_with_frac_mutect_smchet.txt"), 
                    stringsAsFactors = F) 
} else {
  ssm <- read.delim(paste0(dataFolder, "/ssm_data_no_chrxy_with_frac.txt"), 
                    stringsAsFactors = F)
  cnv <- read.delim(paste0(dataFolder, "/cnv_data_no_chrxy_with_frac.txt"), 
                    stringsAsFactors = F) 
}

ssm$major_cn = 1
ssm$minor_cn = 1
ssm$cn_frac = 1
ssm$mu_r <- NULL
ssm$mu_v <- NULL
ssm$cn_ref <- NULL
ssm$cn_tot <- NULL

cnv <- na.omit(cnv)
cnv <- filter(cnv, ssms!="" )

if (nrow(cnv)>0) {
  cnvtmp1 <- strsplit(as.character(cnv$ssms), ";")
  for (j in seq_len(nrow(cnv))) {
    if (length(cnvtmp1[[j]])==0) { next }
    cnvtmp1[[j]] = paste(cnvtmp1[[j]], cnv[j,]$frac, sep="," )
  }
  cnvtmp1 <- unlist(cnvtmp1)
  cnvtmp2 <- Reduce(
    rbind, strsplit(cnvtmp1, ",")
  )
  
  if (is.null(dim(cnvtmp2) )) {
    cnvtmp2 = as.data.frame(t(cnvtmp2), stringsAsFactors=F)
  } else {
    cnvtmp2 = as.data.frame(cnvtmp2, stringsAsFactors=F)
  }
  
  for (j in 2:ncol(cnvtmp2)) {
    cnvtmp2[,j] = as.numeric(cnvtmp2[,j])
  }
  
  ssm <- left_join(ssm, cnvtmp2, by=c("id"="V1"))
  ssm$major_cn <- ssm$V3
  ssm$minor_cn <- ssm$V2
  ssm$cn_frac <- ssm$V4
  
  ssm$V2 <- NULL
  ssm$V3 <- NULL
  ssm$V4 <- NULL
  
  ssm[is.na(ssm[,5]), 5] = 1
  ssm[is.na(ssm[,6]), 6] = 1
  ssm[is.na(ssm[,7]), 7] = 1
} 

clonalCnFrac <- sum(ssm$cn_frac==1)/nrow(ssm)
ssm <- filter(ssm, cn_frac==1)


maxSnv <- 30000
if (nrow(ssm) > maxSnv) {
  ssm <- sample_n(ssm, maxSnv)
}
ssm$normal_cn = 2
ssm <- rename(ssm, ref_counts=a, total_counts=d, var_counts=d-a)

if (caller == "muse") {
  save(ssm, clonalCnFrac, file = paste0(dataFolder, "/ssm_no_chrxy_no_subclonal_cn_", 
                                        caller, "_tier_", muse_tier, ".RData"))
} else if (caller == "dkfz") {
  save(ssm, clonalCnFrac, file = paste0(dataFolder, "/ssm_no_chrxy_no_subclonal_cn_dkfz.RData"))
} else if (caller == "mutect"){
  save(ssm, clonalCnFrac, file = paste0(dataFolder, "/ssm_no_chrxy_no_subclonal_cn_mutect.RData"))
} else if (caller == "mutect_smchet"){
  save(ssm, clonalCnFrac, file = paste0(dataFolder, "/ssm_no_chrxy_no_subclonal_cn_mutect_smchet.RData"))
} else {
  save(ssm, clonalCnFrac, file = paste0(dataFolder, "/ssm_no_chrxy_no_subclonal_cn.RData"))
}

cat("\n finished", i, "in", length(dataFolders), "\n", file = "progress_prepare_data_ssm",append = T)
