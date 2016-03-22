#!/home/yuan03/R-3.2.0/bin/Rscript

library(dplyr)
library(ccube)

args <- commandArgs(trailingOnly = TRUE)
ssm_file <- as.character(args[1])
cnv_file <- as.character(args[2])


ssm <- read.delim(ssm_file,
                  stringsAsFactors = F)
cnv <- read.delim(cnv_file,
                  stringsAsFactors = F)


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
ssm <- rename(ssm, ref_counts=a, total_counts=d)
ssm <- mutate(ssm, var_counts=d-a)
ssm$purity <- GetPurity(ssm)

write.table(ssm, file = "ssm_ccube.txt", sep = "\t", row.names = F, quote = F)
write.table(unique(ssm$purity), file = "1A.txt", sep = "\t", row.names = F, quote = F)


