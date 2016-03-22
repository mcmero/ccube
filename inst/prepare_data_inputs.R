library(dplyr)
icgc <- '/lustre/fmlab/yuan03/PCAWG/db_analysis/'
#icgc <- "~/Desktop/LustreCI/dream_challenge/"

args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])
caller <- as.character(args[2])
jobList <- as.character(args[3])
muse_tier <- as.integer(args[4])

if (caller == "muse") {
  if (is.na(muse_tier)) {
    muse_tier = 5
  }
}

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

# Generate ssm_data and cnv_data with create_phylowgs_inputs.py
ascatFiles <- dir(paste0(icgc, "ascat_calls/"), full.names = T)
ascatFiles1 <- dir(paste0(icgc, "ascat_calls_1/"), full.names = T)
batternbergFiles <- dir(paste0(icgc, "battenberg_calls/"), full.names = T)

dir.create(paste0(icgc, "test_snv/"))

if (is.na(jobList)) {
  sampleList <- sampleNameVcf
} else {
  sampleList <- read.csv(jobList, sep="", stringsAsFactors=FALSE)$x
}
sampleName <- sampleList[i]

batternbergFile <- batternbergFiles[which(sampleNameBB == sampleName)]
if (sampleName %in% sampleNameAscat) {
  ascatFile <- ascatFiles[which(sampleNameAscat == sampleName)]
} else {
  ascatFile <- ascatFiles1[which(sampleNameAscat1 == sampleName)]
}

dir.create(paste(icgc, "test_snv/", sampleName, sep =""))
if ("cellularity" %in% unlist(strsplit(ascatFile, split = "_" )) ) {
  cellularity <- read.table(ascatFile, header = T)$cellularity
} else if ("PP" %in% unlist(strsplit(ascatFile, split = "_" ))) {
  cellularity <- read.table(ascatFile, header = T)$purity
} else {
  cellularity <- read.table(ascatFile, header = T)$rho[2]
}

if (is.null(cellularity)) {
  cellularity <- read.table(ascatFile,stringsAsFactors=FALSE)$V1
}

if (cellularity > 1) {
  cellularity = 1
}

vcfFile <- vcfFiles[which(sampleNameVcf == sampleName)]

shellCommandSanger <- paste(
  "python create_ccfclust_inputs.py -v sanger", 
  " -b ", batternbergFile, 
  " -c ", cellularity, 
  " --output-cnvs ", paste0(icgc, "test_snv/", sampleName, "/cnv_data_no_chrxy_with_frac.txt"), 
  " --output-variants ", paste0(icgc, "test_snv/", sampleName, "/ssm_data_no_chrxy_with_frac.txt"),
  " ", vcfFile, sep = ""
)


shellCommandMuse <- paste(
  "python create_ccfclust_inputs.py -v muse", 
  " -b ", batternbergFile, 
  " -c ", cellularity, 
  " --output-cnvs ", paste0(icgc, "test_snv/", sampleName, 
                            "/cnv_data_no_chrxy_with_frac_muse_tier_", 
                            muse_tier, ".txt"), 
  " --output-variants ", paste0(icgc, "test_snv/", sampleName, 
                                "/ssm_data_no_chrxy_with_frac_muse_tier_",
                                muse_tier, ".txt"),
  " --muse-tier ", muse_tier,
  " --tumor-sample ", sampleName,
  " ", vcfFiles[i], sep = ""
)

shellCommandDkfz <- paste(
  "python create_ccfclust_inputs.py -v dkfz", 
  " -b ", batternbergFile, 
  " -c ", cellularity, 
  " --output-cnvs ", paste0(icgc, "test_snv/", sampleName, 
                            "/cnv_data_no_chrxy_with_frac_dkfz.txt"), 
  " --output-variants ", paste0(icgc, "test_snv/", sampleName, 
                                "/ssm_data_no_chrxy_with_frac_dkfz.txt"),
  " --tumor-sample ", sampleName,
  " ", vcfFiles[i], sep = ""
)

shellCommandMutect <- paste(
  "python create_ccfclust_inputs.py -v mutect_pcawg", 
  " -b ", batternbergFile, 
  " -c ", cellularity, 
  " --output-cnvs ", paste0(icgc, "test_snv/", sampleName, 
                            "/cnv_data_no_chrxy_with_frac_mutect.txt"), 
  " --output-variants ", paste0(icgc, "test_snv/", sampleName, 
                                "/ssm_data_no_chrxy_with_frac_mutect.txt"),
  " --tumor-sample ", sampleName,
  " ", vcfFiles[i], sep = ""
)

shellCommandMutectSmcHet <- paste(
  "python create_ccfclust_inputs.py -v mutect_smchet", 
  " -b ", batternbergFile, 
  " -c ", cellularity, 
  " --output-cnvs ", paste0(icgc, "test_snv/", sampleName, 
                            "/cnv_data_no_chrxy_with_frac_mutect_smchet.txt"), 
  " --output-variants ", paste0(icgc, "test_snv/", sampleName, 
                                "/ssm_data_no_chrxy_with_frac_mutect_smchet.txt"),
  " ", vcfFiles[i], sep = ""
)

if (caller == "muse") {
  cat(shellCommandMuse, "\n")
  system(shellCommandMuse, intern = TRUE)
} else if (caller == "dkfz") {
  cat(shellCommandDkfz, "\n")
  system(shellCommandDkfz, intern = TRUE)
} else if (caller == "mutect"){
  cat(shellCommandMutect, "\n")
  system(shellCommandMutect, intern = TRUE)
} else if (caller == "mutect_smchet") {
  cat(shellCommandMutectSmcHet, "\n")
  system(shellCommandMutectSmcHet, intern = TRUE)
} else {
  cat(shellCommandSanger, "\n")
  system(shellCommandSanger, intern = TRUE)
}

cat(i, "in", length(vcfFiles), "\n", file = "progress_prepare_data_inputs", append = T)
