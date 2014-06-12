#!/usr/bin/Rscript --vanilla

library(PopGenome)
library(ggplot2)

# This script will compute and plot an unfolded site frequency sepectrum given
# an alignment and an outgroup

# Export plot function
ExportPlot <- function(gplot, filename, width=2, height=1.5) {
    ggsave(paste(filename, '.pdf', sep=""), gplot, width=width, height=height)
    postscript(file=paste(filename, '.eps', sep=""), width=width, height=height)
    print(gplot)
    dev.off()
    png(file = paste(filename, '.png', sep=""), width =
        width * 100, height = height * 100)
    print(gplot)
    dev.off()
}

#get command line arguments
args <- commandArgs(trailingOnly = TRUE)

#import alignment
alignment <- readData(args[1])

#set outgroup
alignment <- set.outgroup(alignment, args[2])

#calculate sfs
alignment <- detail.stats(alignment)

alleleFreqs <- alignment@region.stats@minor.allele.freqs[[1]]
freq.table <- list()
freq.table[[1]] <- table(alleleFreqs)

sfs <- data.frame(freq.table)
print(sfs)

p <- qplot(x=sfs$alleleFreqs, y=sfs$Freq) + geom_bar()
ExportPlot(p, "sfs", width = 4, height=4)
