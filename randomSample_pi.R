#!/usr/bin/Rscript --vanilla

# This script produces a plot from Bayesian Skyline Data output by BEAST

# Load required libraries
library(reshape)

randomSample <- function(x, bin1, bin2, bin3, bin4) {
    rand1 <- bin1[sample(nrow(bin1), 14),]
    bin1 <- bin1[-as.numeric(rownames(rand1)),]
    rownames(bin1) <- seq(length=nrow(bin1))
    assign('bin1', bin1, envir=.GlobalEnv)
    rand2 <- bin2[sample(nrow(bin2), 3),]
    bin2 <- bin2[-as.numeric(rownames(rand2)),]
    rownames(bin2) <- seq(length=nrow(bin2))
    assign('bin2', bin2, envir=.GlobalEnv)
    rand3 <- bin3[sample(nrow(bin3), 2),]
    bin3 <- bin3[-as.numeric(rownames(rand3)),]
    rownames(bin3) <- seq(length=nrow(bin3))
    assign('bin3', bin3, envir=.GlobalEnv)
    rand4 <- bin4[sample(nrow(bin4), 1),]
    bin4 <- bin4[-as.numeric(rownames(rand4)),]
    rownames(bin4) <- seq(length=nrow(bin4))
    assign('bin4', bin4, envir=.GlobalEnv)
    filename <- paste("sample", x, ".txt", sep="")
    write.table(rand1, file = filename, row.names=F, col.names=F, 
                quote=F, sep= "\t")
    write.table(rand2, file = filename, row.names=F, col.names=F, 
                quote=F, sep= "\t", append=T)
    write.table(rand3, file = filename, row.names=F, col.names=F, 
                quote=F, sep= "\t", append=T)
    write.table(rand4, file = filename, row.names=F, col.names=F, 
                quote=F, sep= "\t", append=T)
}

# get command line arguments for input file name
args <- commandArgs(trailingOnly = TRUE)

genes <- read.table(args[1], skip = 1, stringsAsFactors=F)
LGTgenes <- scan(file = args[2], what = character())
cat(str(LGTgenes))
# remove genes that have LGT breakpoints
cat(nrow(genes), sep = "\n")
genes <- genes[!genes$V1 %in% LGTgenes,]
cat(nrow(genes))
# remove genes with pi too low to be informative for BEAST

genes <- genes[genes$V3 > quantile(genes$V3,0.25),] 

# remove genes with very high values of pi (> 5 standard deviations above mean)

genes <- genes[genes$V3 < (5*sd(genes$V3) + mean(genes$V3)),]

# put data into bins to sample from
bin1 <- genes[genes$V3 < quantile(genes$V3,0.7),] 
bin2 <- genes[genes$V3 > quantile(genes$V3,0.7) & 
    genes$V3 < quantile(genes$V3,0.85),]
bin3 <- genes[genes$V3 > quantile(genes$V3,0.85) &
    genes$V3 < quantile(genes$V3,0.95),]
bin4 <- genes[genes$V3 > quantile(genes$V3,0.95),]

for (i in 1:20) {
    randomSample(i, bin1, bin2, bin3, bin4)
}
