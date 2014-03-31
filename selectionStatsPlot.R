#!/usr/bin/Rscript --vanilla

# This script produces plots for various measures output by selectionStats.py
# Usage: selectionStatsPlot.R <filename>

# Load required libraries
library(ggplot2)
library(tools)
library(reshape)

# Read data from file provided by user
args <- commandArgs(trailingOnly = TRUE)
data <- read.delim(args[1], header = TRUE)
fileBase <- file_path_sans_ext(args[1])

# Plot theta
tiff(filename = paste(fileBase, "_thetaHistogram.tiff", sep = ""))
ggplot(data, aes(x=Theta)) + 
        geom_histogram() + 
        ggtitle("Theta Distribution")
dev.off()

# Plot pi
tiff(filename = paste(fileBase, "_piHistogram.tiff", sep = "")) 
ggplot(data, aes(x=Pi)) + 
        geom_histogram() + 
        ggtitle("Nucleotide Diversity Distribution")
dev.off()

# Plot piN/piS
data$piNpiS <- data$PiN/data$PiS
tiff(filename = paste(fileBase, "_piNpiSHistogram.tiff", sep = ""))
ggplot(data, aes(x=piNpiS)) + 
        geom_histogram() + 
        ggtitle("PiN/PiS Distribution")
dev.off()

# Plot Tajima's D
data$TajimasD <- as.numeric(as.character(data$TajimasD))
head(data$TajimasD)
tiff(filename = paste(fileBase, "_TajimasD.tiff", sep = "")) 
ggplot(data, aes(x=TajimasD)) + 
        geom_histogram() + 
        ggtitle("Tajima's D Distribution")
dev.off()

# Plot Neutrality Index
statNumber <- ncol(data)
NI <- data[,c(1, seq(8, ncol(data), 2))]
NI.m <- melt(NI)
NI.m <- NI.m[ which(NI.m$value > 0 & NI.m$variable != 'NI_SRS000032'),]
summary(NI.m$value)
tiff(filename = paste(fileBase, "_NI.tiff", sep = ""), 
        width = 480, height = 5000)
ggplot(NI.m, aes(variable, Alignment)) +
        geom_tile(aes(fill = value), colour = "white") + 
        scale_fill_gradient2(trans = 'log')
dev.off()
tiff(filename = paste(fileBase, "_NIHistogram.tiff", sep = "")) 
ggplot(NI.m, aes(x=value)) + 
        geom_histogram() + 
        scale_x_log10() +
        ggtitle("Neutrality Index Distribution")
dev.off()

