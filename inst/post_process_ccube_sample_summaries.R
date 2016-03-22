#!/home/yuan03/R-3.2.0/bin/Rscript
library(ggplot2)
library(colorspace)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])
caller <- as.character(args[2])
jobList <- as.character(args[3])
muse_tier <- as.integer(args[4])


if (caller == "muse") {
  if (is.na(muse_tier)) {
    muse_tier <- 5
  }
  codeName <- paste0("results/", caller, "_tier_", muse_tier, "_battenberg_ccube_v0.4/")
} else {
  codeName <- paste0("results/", caller, "_battenberg_ccube_v0.4/")
}

icgc <- "/lustre/fmlab/yuan03/PCAWG/db_analysis/"
#icgc <- '~/Desktop/LustreCI/PCAWG/broad_simulated_samples/batch_2015-09-24/'
dir.create(paste0(icgc, codeName), recursive = T)


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
  sangerSampleIndexWithBB <- which(sampleNameVcf %in% sampleNameAscat)
  sampleNameVcf <- sampleNameVcf[sangerSampleIndexWithBB]
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

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

myColors <- gg_color_hue(10)


library(ccube)

dataFolders <- dir(paste0(icgc, "test_snv/"), full.names = T)
dataFolderNames <- dir(paste0(icgc, "test_snv/"))
dataFolder <- dataFolders[which(dataFolderNames == sampleName)]

if (caller == "muse") {
 
  resultsFile <- paste0(dataFolder, 
                        "/smm_map_cluster_res_no_chrxy_no_subclonal_cn_1_", 
                        caller, "_tier_", muse_tier, ".RData")
  
} else if (caller == "dkfz") {
  resultsFile <- paste0(dataFolder, 
                        "/smm_map_cluster_res_no_chrxy_no_subclonal_cn_1_dkfz.RData")
} else if (caller == "mutect"){
  resultsFile <- paste0(dataFolder, 
                        "/smm_map_cluster_res_no_chrxy_no_subclonal_cn_1_mutect.RData")
} else if (caller == "mutect_smchet"){
  resultsFile <- paste0(dataFolder, 
                        "/smm_map_cluster_res_no_chrxy_no_subclonal_cn_1_mutect_smchet.RData")
} else {
  resultsFile <- paste0(dataFolder, 
                        "/smm_map_cluster_res_no_chrxy_no_subclonal_cn_1.RData")
}

if (file.exists(resultsFile)) {
  
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
  
  load(resultsFile)  
  uniqLabels <- unique(res$label)
  
  resultsFolder <- paste0(icgc, codeName, sampleName)
  dir.create(resultsFolder, recursive = T)
  
  id <- Reduce(rbind, strsplit(as.character(ssm$gene), "_", fixed = T), c())
  
  # Multiplicity
  
  ssm$mult <- res$bv
  
  mult <- data.frame(chr = id[,1], pos = id[,2])
  mult$tumour_copynumber <- ssm$major_cn+ssm$minor_cn
  mult$multiplicity <- ssm$mult
  fn <- paste0(resultsFolder, "/", 
               sampleName, "_multiplicity.txt")
  write.table(mult, file = fn, sep = "\t", row.names = F, quote = F)
  shellCommand <- paste0("gzip -f ", fn)
  system(shellCommand, intern = TRUE)
  
  mutAssign <- data.frame(chr = id[,1], pos = id[,2])
  
  if (length(uniqLabels) == 1) {
    mutR = data.frame(res$R)
    colnames(mutR) <- "cluster_1"
  } else {
    mutR <- data.frame(res$R[, sort(uniqLabels)]) 
    colnames(mutR) <- paste0("cluster_", seq_along(uniqLabels))
  }
  
  mutAssign <- data.frame(mutAssign, mutR)
  fn <- paste0(resultsFolder, "/", 
               sampleName, "_assignment_probability_table.txt")
  write.table(mutAssign, file = fn, sep = "\t", row.names = F, quote = F)
  shellCommand <- paste0("gzip -f ", fn)
  system(shellCommand, intern = TRUE)

  clusterCertainty <- as.data.frame(table(res$label), stringsAsFactors = F)
  clusterCertainty <- rename(clusterCertainty, cluster = Var1, n_ssms = Freq)
  clusterCertainty$proportion <- res$mu[as.integer(clusterCertainty$cluster)] * cellularity
  clusterCertainty$cluster <- seq_along(uniqLabels) 
#    clusterCertainty <- data.frame(chr = id[,1], pos = id[,2])
#    clusterCertainty$most_likely_assignment <- match(res$label, sort(unique(res$label)))	
#    clusterCertainty$average_ccf <- res$mu[res$label]
#    ci_95 <- 2*sqrt(res$full.model$invWhishartScale*(res$full.model$normalRelativePrecision+1) 
#                    / (res$full.model$whishartDof * res$full.model$normalRelativePrecision^2 ))
#    clusterCertainty$lower_95_ci <- res$mu[res$label] - ci_95[res$label]
#    clusterCertainty$upper_95_ci <- res$mu[res$label] + ci_95[res$label]
  fn <- paste0(resultsFolder, "/", 
               sampleName, "_subclonal_structure.txt")
  write.table(clusterCertainty, file = fn, sep = "\t", row.names = F, quote = F)
  shellCommand <- paste0("gzip -f ", fn)
  system(shellCommand, intern = TRUE)

  # summary graph
  fn = paste0(resultsFolder, "/", 
              sampleName, "_results_summary.pdf")
  ppi <- 500
  pdf(fn, width=8, height=8)
  par(mfrow=c(2,2))

  vaf <- ssm[ssm$id %in% ccf$ssm.id, ]$vaf
  plot(ccf$ccf, vaf, col = myColors[res$label], 
       xlab = "cancer cell fraction", ylab = "variant allele frequecy", 
       main = "ccf vs vaf (colored by cluster memebership)")
  ssm$tot_cn =ssm$cn_maj+ssm$cn_min
  uniqueTotCn = unique(ssm$tot_cn)
  xx = seq(0,2, length.out = 100)
  for (cn in uniqueTotCn) {
    for (i in 1:cn) {
      points(MapVaf2CcfPyClone(xx, cellularity, 2, cn, cn, i, constraint = F), xx, type = 'l')
    }
  }
  
  
  Emu <- res$full.model$normalMean
  Esigma <- res$full.model$invWhishartScale/res$full.model$whishartDof
  Epi <- res$full.model$Epi
  EgammaDof <- res$full.model$gammaDof
  
  params <- data.frame(Emu, Esigma, Epi)
  xx <- seq(range(ccf$ccf)[1],range(ccf$ccf)[2],  length.out = 1000)
  ll <- 0
  ll1 <- 0 
  
  for (j in seq_len(nrow(params))) {
    ll <- ll + params[j,]$Epi * dnorm(xx, mean = params[j,]$Emu, sd = sqrt(params[j,]$Esigma))
  }
  
  hist(ccf$ccf, density=20, breaks=20, prob=TRUE, 
             main = "ccf histogram +
       fitted marginal denstiy (red)",
             xlab = "cancer cell fraction")
  lines(xx,ll, lwd=2, col = "darkred")
  
  names(Epi) <- as.character(format(round(Emu, 2), nsmall = 2))
  barplot(Epi[sort(uniqLabels)], las = 2, col = myColors[sort(uniqLabels)], 
          xlab = "cluster mean", ylab="expected weight", 
          main = "cluster weights")
  
  dev.off()
  
} 

