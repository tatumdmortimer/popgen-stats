#!/usr/bin/Rscript --vanilla

# This script produces plots for various measures output by selectionStats.py
# Usage: selectionStatsPlot.R <filename>

# Load required libraries
library(ggplot2)
library(tools)

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

# Plot MK tests

