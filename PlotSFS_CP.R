# Script to plot SFS
rm(list=ls())
# Script to fold SFS (given the fixed diffs also)
fold_sfs = function(sfs){
  folded_sfs = array(0, floor(length(sfs)/2))
  for (i in seq(1, length(sfs)-1)){
    if (i <= floor(length(sfs)/2)){
      folded_sfs[i] = sfs[i]
    } else{
      folded_sfs[length(sfs)-i] = folded_sfs[length(sfs)-i] + sfs[i]
    }
  }
  return(folded_sfs)
}
#Here we read in all the files:

#setwd("C:\\Documents and Settings\\Caitlin\\My Documents\\prfWork\\Input_Output_files\\5_5_11v2")
observed = read.table("nSNP")[,1]
neutral = read.table("output_constant", skip = 3, sep = "=")[-47,2]
neutral = read.table("output_constant", skip = 3, sep = "=")[,2]
# NB If there are an odd number of sequences, the last frequency
# will need to be moved to the top of the file and skipped or else the function won't work
# same goes for the files below
exp = read.table("output_fixed_expansion", skip = 4, sep = "=")[-47,2]
exp = read.table("output_fixed_expansion", skip = 4, sep = "=")[,2]
sel = read.table("output_fixed_expansion_ptmass", skip = 4, sep = "=")[-47,2]
sel = read.table("output_fixed_expansion_ptmass", skip = 4, sep = "=")[,2]
sel2 = read.table("output_fixed_expansion_ptmass_neut", skip = 4, sep = "=")[-47,2]
sel2 = read.table("output_fixed_expansion_ptmass_neut", skip = 4, sep = "=")[,2]

#You need to fold all the SFS first, because you will disregard the last column of fixed differences:

folded_neutral<- fold_sfs(neutral)
folded_exp<- fold_sfs(exp)
folded_sel<- fold_sfs(sel)
folded_sel2<- fold_sfs(sel2)
neutral_scaled = folded_neutral*(sum(observed)/sum(folded_neutral))
exp_scaled = folded_exp*(sum(observed)/sum(folded_exp))
sel_scaled = folded_sel*(sum(observed)/sum(folded_sel))
sel2_scaled = folded_sel2*(sum(observed)/sum(folded_sel2))

#Make a matrix of observed + 3 models:
sfs_matrix = matrix(c(observed, neutral_scaled,
  exp_scaled, sel_scaled), byrow =
  T, nrow = 4)

#Make a matrix of observed + 4 models:
sfs_matrix = matrix(c(observed, neutral_scaled, exp_scaled, sel_scaled, sel2_scaled), byrow =
  T, nrow = 5)

#Drop the neutral, constant sized model:
sfs_matrix = matrix(c(observed, exp_scaled, sel_scaled, sel2_scaled), byrow =
  T, nrow = 4)

#Plot with 3 models:

par(pty="s")
barplot(sfs_matrix, beside = T, names.arg =
        as.character(seq(1, ncol(sfs_matrix))), col = c(1,2,3,4), main =
        "Site Frequency Spectrum: nSNPs, group 2", xlab = "Frequency", ylab =
        "Number of SNPs")
legend("topright", legend = c("Observed", "Neutral Constant", "Neutral Expansion",
                     "Expansion + Selection"), fill = c(1,2,3,4))

#Plot with 4 models:
palette(rainbow(20))#to set rainbow palette
palette("default")#to set default colors

plot(1:20, pch=CIRCLE<-16, cex=1:20, col=rainbow(20))#you can look at palette
par(mfrow=c(1,1))
barplot(sfs_matrix, beside = T, names.arg =
        as.character(seq(1, ncol(sfs_matrix))), col = c(1,2,3,4,5), main =
        "Observed and Expected nSNPs", xlab = "Frequency", ylab =
        "Number of SNPs")
legend("topright", legend = c("Observed", "Neutral Constant", "Neutral Expansion",
                     "Expansion + Selection", "Expansion + 2 param"), fill = c(1,2,3,4,5))

#Drop the neutral, constant sized model:
barplot(sfs_matrix, beside = T, names.arg =
        as.character(seq(1, ncol(sfs_matrix))), col = c(1,2,3,4), main =
        "Observed & Expected SFS (nSNPs)", xlab = "Frequency", ylab =
        "Number of SNPs")
legend("topright", legend = c("Observed", "Neutral Expansion",
                     "Expansion + Selection", "Expansion + 2 param"), fill = c(1,2,3,4))

#Save the figure:

dev.print(device=postscript, "Observed_VS_Expected2.eps", onefile=FALSE, horizontal=FALSE)

dev.off()
#this resets all the par settings to default
################################################
#TO MAKE A SIMPLE FIGURE OF OBSERVED NSNP AND SSNP SFS

#par(mfrow=c(1,2))
layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths = c(1.4,1))
nSNP_one = read.table("Gp1_non_syn_SFS_10_20_11")[,1]
sSNP_one = read.table("Gp1_syn_SFS_10_20_11")[,1]

sfs_matrix1 = matrix(c(nSNP_one, sSNP_one), byrow =
  T, nrow = 2)

barplot(sfs_matrix1, beside = T, names.arg =
        as.character(seq(1, ncol(sfs_matrix1))), col = c("blue", "red"), main =
        "Group 1", xlab = "Frequency", ylab =
        "Number of SNPs", ylim = c(0,4000))
legend("topright", legend = c("non-synonymous", "synonymous"), fill = c("blue", "red"))

nSNP_two = read.table("Gp2_non_syn_SFS_10_20_11")[,1]
sSNP_two = read.table("Gp2_syn_SFS_10_20_11")[,1]

sfs_matrix2 = matrix(c(nSNP_two, sSNP_two), byrow =
  T, nrow = 2)

barplot(sfs_matrix2, beside = T, names.arg =
        as.character(seq(1, ncol(sfs_matrix2))), col = c("blue", "red"), main =
        "Group 2", xlab = "Frequency", ylim = c(0,2000))
#legend("topright", legend = c("non-synonymous", "synonymous"), fill = c("blue", "red"))

