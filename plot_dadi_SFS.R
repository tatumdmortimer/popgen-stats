#!/usr/bin/Rscript --vanilla

# This script produces plots comparing the SFS of various models from dadi

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

observed_SFS <- read.table("observedSFS.txt")
neutral_SFS <- read.table("neutralModelSFS.txt")

SFS <- rbind(data.frame(Model="Observed", 
        Count=observed_SFS$V1,
        Frequency=as.numeric(row.names(observed_SFS))),
    data.frame(Model="Neutral", 
        Count= neutral_SFS$V1,
        Frequency=as.numeric(row.names(observed_SFS))))

if (file.exists("expansionModelSFS.txt")) {
    expansion_SFS <- read.table("expansionModelSFS.txt")
    SFS <- rbind(data.frame(Model = SFS$Model, 
        Count=SFS$Count,
        Frequency=SFS$Frequency),
    data.frame(Model = "Expansion", 
        Count = expansion_SFS$V1, 
        Frequency=as.numeric(row.names(expansion_SFS))))

}

if (file.exists("growthModelSFS.txt")) { 
    growth_SFS <- read.table("growthModelSFS.txt")
    SFS <- rbind(data.frame(Model = SFS$Model, 
            Count=SFS$Count,
            Frequency=SFS$Frequency),
    data.frame(Model = "Exponential Growth", 
            Count = growth_SFS$V1, 
            Frequency=as.numeric(row.names(growth_SFS))))

}

sfs_plot <- ggplot(SFS, aes(x=Frequency, y=Count, fill=Model)) +
    geom_bar(stat="identity", position="dodge") + 
    theme_bw() + 
    scale_fill_brewer(palette="Paired") 

exportPlot(sfs_plot, "dadiSFS", width = 4, height = 3)
