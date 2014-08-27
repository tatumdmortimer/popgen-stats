#!/usr/bin/Rscript --vanilla

# This script produces a plot from Bayesian Skyline Data output by BEAST
# BSP data must be output from Tracer

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

# get command line arguments for input file name
args <- commandArgs(trailingOnly = TRUE)

BSP <- read.table(args[1], header = TRUE, skip = 1, stringsAsFactors=F)

BSP$Time <- as.numeric(BSP$Time)
BSP$Mean <- as.numeric(BSP$Mean)
BSP$Median <- as.numeric(BSP$Median)
BSP$Upper <- as.numeric(BSP$Upper)
BSP$Lower <- as.numeric(BSP$Lower)

plot <- ggplot(BSP, aes(x=Time)) +         
        geom_ribbon(aes(ymin=Lower, ymax=Upper)) + 
        geom_line(aes(y=Median)) +
        theme_bw() + 
        scale_y_log10()

exportPlot(plot, "BSP", width = 3.5, height = 2)
