#!/usr/bin/Rscript --vanilla

# This script produces plots of the rate for genes with HPD & a distribution
# of corrected rates

# Load required libraries
library(ggplot2)

exportPlot <- function(gplot, filename, width=2, height=1.5) {
    ggsave(paste(filename,'.pdf',sep=""), gplot, width=width, height=height)
    postscript(file=paste(filename,'.eps',sep=""), width=width, height=height)
    print(gplot)
    dev.off()
    png(file=paste(filename,'.png', sep=""), width=width*100, height=height*100)
    print(gplot)
    dev.off()
}

args <- commandArgs(trailingOnly = TRUE)
r <- read.table(args[1])
correction <- as.numeric(args[2])
r <- r[r$V4 < 1 & r$V2 < 1 & r$V7 > 200,]
r$correctedRate <- r$V4*correction
r$V1 <- factor(r$V1, levels=r$V1[order(r$V4)])

cat(summary(r$correctedRate))
cat('\n')
cat(quantile(r$correctedRate, c(0.025,0.975)))
cat('\n')
cat(sd(r$correctedRate))
cat('\n')

p1 <- ggplot(r, aes(x=correctedRate)) + geom_histogram() + theme_bw() +
    xlab("Rate") + scale_x_log10()
exportPlot(p1, "correctedRateHistogram", width=3, height=3)
p2 <- ggplot(r, aes(x=V1, y=V4)) + geom_point() + 
    geom_errorbar(aes(ymin=V5, ymax=V6)) + theme_bw() + 
    xlab("Gene") + ylab("Rate") + scale_y_log10() + 
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
exportPlot(p2, "medianHPDRatesPlot", width=4, height=4)
