rm(list=ls())
# Script to plot LL grid from prfreq output
### Cut and paste the following for all analyses
setwd("C:\\Documents and Settings\\Caitlin\\My Documents\\prfWork\\prfreq_neutral")
setwd("C:\\Documents and Settings\\Caitlin\\My Documents\\prfWork\\Input_Output_files\\5_6_11")
#setwd("C:\\Documents and Settings\\Caitlin\\My Documents\\prfWork\\prfreq_selection")
#setwd("/Users/jgranka/Documents/TB/selection/data/prfreq_selection")
#set my own working directory here

# neutral
#param_grid = read.table("exp_rv_s_folded_grid_2.txt", skip = 16) # has to have the
                                        # last lines removed  # skip = 13
param_grid = read.table("grid_v4.txt", skip = 21)
param_grid = read.table("grid.txt")
#param_grid3 = read.table("grid3.txt", skip = 20)
#delete everything after the grid of values then specify how many lines to skip before getting to the grid of ML values 
#specify the filename here
                                        # selection
#param_grid = read.table("ptmass_lethal_exp_allSmu_grid_2.txt", skip = 29) # has to have the
                                        # last lines removed
###second set of instructions for all analyses
# Make matrix for heat maps
# columns for params
col1 = 1  # these would be changed if other params are examined and output
col2 = 2
LLcol = 5 #5 for neutral, 4 for selection this is where you specify the column containing the likelihoods

param1 = unique(param_grid[,col1])
param2 = unique(param_grid[,col2])
LL_matrix = matrix(ncol = length(param2), nrow = length(param1))

# Make matrix:
for (i in seq(1, length(param1))){
  # these will be the columns
  for (j in seq(1, length(param2))){
    LL_matrix[i, j] = param_grid[length(param2)*(i-1) + j, LLcol]
  }
}

### Run the above for both demographic inference and selection inference

# Demography grid:
# Make contour plot
#pdf("exp_allSmu_s_plot_2.pdf")
filled.contour(-LL_matrix, x = seq(1,nrow(LL_matrix)), y =
               seq(1,ncol(LL_matrix)), key.title = title("log LL",
               cex.main = 0.8), key.axes = axis(4, cex.axis = 0.7),
               color.palette = topo.colors, nlevels = 50, plot.axes =
               {axis(1, at = seq(1, nrow(LL_matrix)), labels =
                     as.character(round(param1, digits = 4)), cex.axis = 0.7)
axis(2, at = seq(1, ncol(LL_matrix)), labels =
     as.character(round(param2, digits = 4)),
               cex.axis = 0.7)}, xlab = expression(omega), ylab =
               expression(tau))
#dev.off()
#note that filled.contour does not support multiple panels
               
#use the following for very fine grids where you don't want to have each tick mark and label shown
#eps("group2-grid3.eps")
filled.contour(-LL_matrix, x = seq(1,nrow(LL_matrix)), y =
               seq(1,ncol(LL_matrix)), key.title = title("log LL",
               cex.main = 0.8), key.axes = axis(4, cex.axis = 0.7),
               color.palette = topo.colors, nlevels = 50, plot.axes =
               {axis(1, at = seq(1, nrow(LL_matrix), length.out = 10), labels =
                     as.character(round(seq(min(param1), max(param1),length.out = 10), digits=4)), cex.axis = 0.7)
axis(2, at = seq(1, ncol(LL_matrix), length.out=10), labels =
      as.character(round(seq(min(param2), max(param2),length.out = 10), digits=4)),
               cex.axis = 0.7)}, xlab = expression(omega), ylab =
               expression(tau))

#Another alternative:
require(lattice)
contourplot(-LL_matrix, aspect="fill", region=TRUE, col.regions=topo.colors, contour=TRUE, labels=FALSE, cuts=50)
#need to do some work to figure out how to fix the axis labels and [maybe] draw only the top contour
#colorkey(col=topo.colors, at=50,tick.number=50)
               
# Selection:
# Make contour plot for selection parameters
pdf("ptmass_neutral_exp_allSmu_2.pdf")
filled.contour(-LL_matrix, x = seq(1,nrow(LL_matrix)), y =
               seq(1,ncol(LL_matrix)), key.title = title("log LL",
               cex.main = 0.8), key.axes = axis(4, cex.axis = 0.7),
               color.palette = topo.colors, nlevels = 50, plot.axes =
               {axis(1, at = seq(1, nrow(LL_matrix)), labels =
                     as.character(round(param1, digits = 4)), cex.axis = 0.7)
axis(2, at = seq(1, ncol(LL_matrix)), labels =
     as.character(round(param2, digits = 4)),
               cex.axis = 0.7)}, xlab = expression(alpha), ylab =
               expression(beta))
dev.off()
## xlab = "Proportion Neutral (p)", ylab =
##                expression(gamma))

# IF not a grid
setwd("C:/Documents and Settings/Caitlin/My Documents/prfWork/Input_Output_files/5_10_11")
param_grid = read.table("grid.txt") # , skip = 31; has to have the
                                        # last lines removed
col1 = 1  # these would be changed if other params are examined and output
LLcol = 4

#pdf("PtMass_exp_allSmu_plot_2.pdf")
plot(param_grid[,col1], -param_grid[,LLcol], type = 'l', xlab = expression(gamma), ylab = "Log Likelihood", main =
     "Point Mass Inference")
#dev.off()

# Plot a gamma
alpha = 0.8
beta = 0.1
x = seq(-10, 10, by = 0.1)
y = dgamma(x, alpha, beta)
plot(-x, y, type = 'l', xlab = expression(gamma), ylab = "PDF")

#######################################################################################
#Find 95% CIs from a grid file
#This is the LL of the parameter MLE:
MLE_LL = -61156.210232
#For a demography file you need to change the column of the param_grid
#For a 1 parameter model, you need estimates within 1 LL unit
CI <- which(param_grid[,4] <= (MLE_LL + 2), arr.ind= TRUE); param_grid[CI,]

