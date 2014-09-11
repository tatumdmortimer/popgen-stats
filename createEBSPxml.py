#!/usr/bin/env python

import sys
import os
import getopt
import glob
from Bio import SeqIO

# This script reads in a directory of fasta alignments 
# and creates the input xml for an EBSP analysis using BEAST v 1.8

def usage():
    print "createEBSPxml.py <directory of fasta alignments> <file name prefix>"

if len(sys.argv) != 3:
    usage()
    sys.exit()

directory = sys.argv[1]
prefix = sys.argv[2]

alignments = glob.glob(directory + '*.fasta')
out_xml = open(prefix + '.xml', 'w')

out_xml.write("<beast>\n")

if len(alignments) == 0:
    print ("Directory must contain files ending with .fasta")

alignmentNames = []
taxa = []

# get taxa names from first alignment
handle = open(alignments[0], "r")
for record in SeqIO.parse(handle, "fasta"):
    taxa.append(record.id)
handle.close()

# write taxa section of xml file

out_xml.write("\t<taxa id=\"taxa\">\n")
for t in taxa:
    out_xml.write("\t\t<taxon id=\"%s\"/>\n" % t)
out_xml.write("\t</taxa>\n")

# write alignments to xml file

out_xml.write("\n")
for i,a in enumerate(alignments):
    alignmentNames.append(os.path.basename(a).split('.')[0])
    out_xml.write("\t<alignment id=\"alignment%s\" dataType=\"nucleotide\">\n" 
                    % str(i+1))
    handle = open(a, "r")
    for record in SeqIO.parse(handle, "fasta"):
        out_xml.write("\t\t<sequence>\n")
        out_xml.write("\t\t\t<taxon idref=\"%s\"/>\n" % record.id)
        out_xml.write("\t\t\t%s\n" % record.seq)
        out_xml.write("\t\t</sequence>\n")
    out_xml.write("\t</alignment>\n")

#write patterns to xml file

out_xml.write("\n")
for j,n in enumerate(alignmentNames):
    out_xml.write("\t<patterns id=\"%s.patterns\" from=\"1\" strip=\"false\">\n"
                    % n)
    out_xml.write("\t\t<alignment idref=\"alignment%s\"/>\n" % str(j+1))
    out_xml.write("\t</patterns>\n")

# write starting trees to xml file

out_xml.write("\n")
out_xml.write("\t<constantSize id=\"initialDemo\" units=\"substitutions\">\n")
out_xml.write("\t\t<populationSize>\n")
out_xml.write("\t\t\t<parameter id=\"initialDemo.popSize\" value=\"100.0\"/>\n")
out_xml.write("\t\t</populationSize>\n")
out_xml.write("\t</constantSize>\n")

for name in alignmentNames:
    out_xml.write("\t<coalescentSimulator id=\"%s.startingTree\">\n" % name)
    out_xml.write("\t\t<taxa idref=\"taxa\"/>\n")
    out_xml.write("\t\t<constantSize idref=\"initialDemo\"/>\n")
    out_xml.write("\t</coalescentSimulator>\n")

# write tree model to xml file

for aN in alignmentNames:
    out_xml.write("\t<treeModel id=\"%s.treeModel\">\n" % aN)
    out_xml.write("\t\t<coalescentTree idref=\"%s.startingTree\"/>\n" % aN)
    out_xml.write("\t\t<rootHeight>\n")
    out_xml.write("\t\t\t<parameter id=\"%s.treeModel.rootHeight\"/>\n" % aN)
    out_xml.write("\t\t</rootHeight>\n")
    out_xml.write("\t\t<nodeHeights internalNodes=\"true\">\n")
    out_xml.write("\t\t\t<parameter id=\"%s.treeModel.internalNodeHeights\"/>\n"
                    % aN)
    out_xml.write("\t\t</nodeHeights>\n")
    out_xml.write("\t\t<nodeHeights internalNodes=\"true\" rootNode=\"true\">\n")
    out_xml.write("\t\t\t<parameter id=\"%s.treeModel.allInternalNodeHeights\"/>\n"
                    % aN)
    out_xml.write("\t\t</nodeHeights>\n")
    out_xml.write("\t</treeModel>\n")
   
# write EBSP to xml file

out_xml.write("\n")
out_xml.write("\t<variableDemographic id=\"demographic\" type=\"linear\" useMidpoints=\"true\">\n")
out_xml.write("\t\t<populationSizes>\n")
out_xml.write("\t\t\t<parameter id=\"demographic.popSize\" value=\"0.023\"/>\n")
out_xml.write("\t\t</populationSizes>\n")
out_xml.write("\t\t<indicators>\n")
out_xml.write("\t\t\t<parameter id=\"demographic.indicators\" value=\"0.0\"/>\n"
                )
out_xml.write("\t\t</indicators>\n")
out_xml.write("\t\t<trees>\n")
for align in alignmentNames:
    out_xml.write("\t\t\t<ptree ploidy=\"1.0\">\n")
    out_xml.write("\t\t\t\t<treeModel idref=\"%s.treeModel\"/>\n" % align)
    out_xml.write("\t\t\t</ptree>\n")
out_xml.write("\t\t</trees>\n")
out_xml.write("\t</variableDemographic>\n")

out_xml.write("\t<coalescentLikelihood id=\"coalescent\">\n")
out_xml.write("\t\t<model>\n")
out_xml.write("\t\t\t<variableDemographic idref=\"demographic\"/>\n")
out_xml.write("\t\t</model>\n")
out_xml.write("\t</coalescentLikelihood>\n")
out_xml.write("\t<sumStatistic id=\"demographic.populationSizeChanges\" \
elementwise=\"true\">\n")
out_xml.write("\t\t<parameter idref=\"demographic.indicators\"/>\n")
out_xml.write("\t</sumStatistic>\n")
out_xml.write("\t<exponentialDistributionModel \
id=\"demographic.populationMeanDist\">\n")
out_xml.write("\t\t<mean>\n")
out_xml.write("\t\t\t<parameter id=\"demographic.populationMean\" \
value=\"0.023\"/>\n")
out_xml.write("\t\t</mean>\n")
out_xml.write("\t</exponentialDistributionModel>\n")

# write clock model to xml file
for aName in alignmentNames:
    out_xml.write("\t<strictClockBranchRates id=\"%s.branchRates\">\n" % aName)
    out_xml.write("\t\t<rate>\n")
    out_xml.write("\t\t\t<parameter id=\"%s.clock.rate\" value=\"1.0\" \
lower=\"0.0\"/>\n" % aName)
    out_xml.write("\t\t</rate>\n")
    out_xml.write("\t</strictClockBranchRates>\n")

# write substitution and site models to xml file
out_xml.write("\n")
# for alignName in alignmentNames:
#     out_xml.write("\t<gtrModel id=\"%s.gtr\">\n" % alignName)
#     out_xml.write("\t\t<frequencies>\n")
#     out_xml.write("\t\t\t<frequencyModel dataType=\"nucleotide\">\n")
#     out_xml.write("\t\t\t\t<frequencies>\n")
#     out_xml.write("\t\t\t\t\t<parameter id=\"%s.frequencies\" value=\"0.25 \
# 0.25 0.25 0.25\"/>\n" % alignName)
#     out_xml.write("\t\t\t\t</frequencies>\n")
#     out_xml.write("\t\t\t</frequencyModel>\n")
#     out_xml.write("\t\t</frequencies>\n")
#     out_xml.write("\t\t<rateAC>\n")
#     out_xml.write("\t\t\t<parameter id=\"%s.ac\" value=\"1.0\" \
# lower=\"0.0\"/>\n" % alignName)
#     out_xml.write("\t\t</rateAC>\n")
#     out_xml.write("\t\t<rateAG>\n")
#     out_xml.write("\t\t\t<parameter id=\"%s.ag\" value=\"1.0\" \
# lower=\"0.0\"/>\n" % alignName)
#     out_xml.write("\t\t</rateAG>\n")
#     out_xml.write("\t\t<rateAT>\n")
#     out_xml.write("\t\t\t<parameter id=\"%s.at\" value=\"1.0\" \
# lower=\"0.0\"/>\n" % alignName)
#     out_xml.write("\t\t</rateAT>\n")
#     out_xml.write("\t\t<rateCG>\n")
#     out_xml.write("\t\t\t<parameter id=\"%s.cg\" value=\"1.0\" \
# lower=\"0.0\"/>\n" % alignName)
#     out_xml.write("\t\t</rateCG>\n")
#     out_xml.write("\t\t<rateGT>\n")
#     out_xml.write("\t\t\t<parameter id=\"%s.gt\" value=\"1.0\" \
# lower=\"0.0\"/>\n" % alignName)
#     out_xml.write("\t\t</rateGT>\n")
#     out_xml.write("\t</gtrModel>\n")
#     out_xml.write("\t<siteModel  id=\"%s.siteModel\">\n" % alignName)
#     out_xml.write("\t\t<substitutionModel>\n")
#     out_xml.write("\t\t\t<gtrModel idref=\"%s.gtr\"/>\n" % alignName)
#     out_xml.write("\t\t</substitutionModel>\n")
#     out_xml.write("\t\t<gammaShape gammaCategories=\"4\">\n")
#     out_xml.write("\t\t\t<parameter id=\"%s.alpha\" value=\"0.5\" \
# lower=\"0.0\"/>\n" % alignName)
#     out_xml.write("\t\t</gammaShape>\n")
#     out_xml.write("\t\t<proportionInvariant>\n")
#     out_xml.write("\t\t\t<parameter id=\"%s.pInv\" value=\"0.5\" \
# lower=\"0.0\" upper=\"1.0\"/>\n" % alignName)
#     out_xml.write("\t\t</proportionInvariant>\n")
#     out_xml.write("\t</siteModel>\n")

for alignName in alignmentNames:
    out_xml.write("\t<HKYModel id=\"%s.hky\">\n" % alignName)
    out_xml.write("\t\t<frequencies>\n")
    out_xml.write("\t\t\t<frequencyModel dataType=\"nucleotide\">\n")
    out_xml.write("\t\t\t\t<frequencies>\n")
    out_xml.write("\t\t\t\t\t<parameter id=\"%s.frequencies\" value=\"0.25 \
0.25 0.25 0.25\"/>\n" % alignName)
    out_xml.write("\t\t\t\t</frequencies>\n")
    out_xml.write("\t\t\t</frequencyModel>\n")
    out_xml.write("\t\t</frequencies>\n")
    out_xml.write("\t\t<kappa>\n")
    out_xml.write("\t\t\t<parameter id=\"%s.kappa\" value=\"2.0\" \
lower=\"0.0\"/>\n" % alignName)
    out_xml.write("\t\t</kappa>\n")
    out_xml.write("\t</HKYModel>\n")
    out_xml.write("\t<siteModel  id=\"%s.siteModel\">\n" % alignName)
    out_xml.write("\t\t<substitutionModel>\n")
    out_xml.write("\t\t\t<HKYModel idref=\"%s.hky\"/>\n" % alignName)
    out_xml.write("\t\t</substitutionModel>\n")
    out_xml.write("\t\t<gammaShape gammaCategories=\"4\">\n")
    out_xml.write("\t\t\t<parameter id=\"%s.alpha\" value=\"0.5\" \
lower=\"0.0\"/>\n" % alignName)
    out_xml.write("\t\t</gammaShape>\n")
    out_xml.write("\t\t<proportionInvariant>\n")
    out_xml.write("\t\t\t<parameter id=\"%s.pInv\" value=\"0.5\" \
lower=\"0.0\" upper=\"1.0\"/>\n" % alignName)
    out_xml.write("\t\t</proportionInvariant>\n")
    out_xml.write("\t</siteModel>\n")

# write tree likelihood to xml file

out_xml.write("\n")
for alName in alignmentNames:
    out_xml.write("\t<treeLikelihood id=\"%s.treeLikelihood\" \
useAmbiguities=\"false\">\n" % alName)
    out_xml.write("\t\t<patterns idref=\"%s.patterns\"/>\n" % alName)
    out_xml.write("\t\t<treeModel idref=\"%s.treeModel\"/>\n" % alName)
    out_xml.write("\t\t<siteModel idref=\"%s.siteModel\"/>\n" % alName)
    out_xml.write("\t\t<strictClockBranchRates idref=\"%s.branchRates\"/>\n" 
                    % alName)
    out_xml.write("\t</treeLikelihood>\n")

# write operators to xml file

out_xml.write("\n")
out_xml.write("\t<operators id=\"operators\" \
optimizationSchedule=\"default\">\n")
# for b in alignmentNames:
#     for c in ["ac", "ag", "at", "cg", "gt", "alpha", "pInv"]:
#         out_xml.write("\t\t<scaleOperator scaleFactor=\"0.75\" \
# weight=\"0.1\">\n")
#         out_xml.write("\t\t\t<parameter idref=\"%s.%s\"/>\n" % (b,c))
#         out_xml.write("\t\t</scaleOperator>\n")
#     out_xml.write("\t\t<deltaExchange delta=\"0.01\" weight=\"0.1\">\n")
#     out_xml.write("\t\t\t<parameter idref=\"%s.frequencies\"/>\n" % b)
#     out_xml.write("\t\t</deltaExchange>\n")

for b in alignmentNames:
    for c in ["alpha", "pInv", "kappa"]:
        out_xml.write("\t\t<scaleOperator scaleFactor=\"0.75\" \
weight=\"0.1\">\n")
        out_xml.write("\t\t\t<parameter idref=\"%s.%s\"/>\n" % (b,c))
        out_xml.write("\t\t</scaleOperator>\n")
    out_xml.write("\t\t<deltaExchange delta=\"0.01\" weight=\"0.1\">\n")
    out_xml.write("\t\t\t<parameter idref=\"%s.frequencies\"/>\n" % b)
    out_xml.write("\t\t</deltaExchange>\n")

for d in alignmentNames:
    out_xml.write("\t\t<scaleOperator scaleFactor=\"0.75\" weight=\"3\">\n")
    out_xml.write("\t\t\t<parameter idref=\"%s.clock.rate\"/>\n" % d)
    out_xml.write("\t\t</scaleOperator>\n")
out_xml.write("\t\t<upDownOperator scaleFactor=\"0.75\" weight=\"30\">\n")
out_xml.write("\t\t\t<up>\n")
for index,e in enumerate(alignmentNames):
    if index > 0:
        out_xml.write("\t\t\t\t<parameter idref=\"%s.clock.rate\"/>\n" % e)
out_xml.write("\t\t\t</up>\n")
out_xml.write("\t\t\t<down>\n")
out_xml.write("\t\t\t\t<parameter idref=\"demographic.populationMean\"/>\n")
out_xml.write("\t\t\t\t<parameter idref=\"demographic.popSize\"/>\n")
for f in alignmentNames:
    out_xml.write("\t\t\t\t<parameter \
idref=\"%s.treeModel.allInternalNodeHeights\"/>\n" % f)
out_xml.write("\t\t\t</down>\n")
out_xml.write("\t\t</upDownOperator>\n")
for g in alignmentNames:
    out_xml.write("\t\t<subtreeSlide size=\"0.0017\" gaussian=\"true\" \
weight=\"15\">\n")
    out_xml.write("\t\t\t<treeModel idref=\"%s.treeModel\"/>\n" % g)
    out_xml.write("\t\t</subtreeSlide>\n")
    out_xml.write("\t\t<narrowExchange weight=\"15\">\n")
    out_xml.write("\t\t\t<treeModel idref=\"%s.treeModel\"/>\n" % g)
    out_xml.write("\t\t</narrowExchange>\n")
    out_xml.write("\t\t<wideExchange weight=\"3\">\n")
    out_xml.write("\t\t\t<treeModel idref=\"%s.treeModel\"/>\n" % g)
    out_xml.write("\t\t</wideExchange>\n")
    out_xml.write("\t\t<wilsonBalding weight=\"3\">\n")
    out_xml.write("\t\t\t<treeModel idref=\"%s.treeModel\"/>\n" % g)
    out_xml.write("\t\t</wilsonBalding>\n")
    out_xml.write("\t\t<scaleOperator scaleFactor=\"0.75\" \
weight=\"3\">\n")
    out_xml.write("\t\t\t<parameter idref=\"%s.treeModel.rootHeight\"/>\n" % g)
    out_xml.write("\t\t</scaleOperator>\n")
    out_xml.write("\t\t<uniformOperator weight=\"30\">\n")
    out_xml.write("\t\t\t<parameter \
idref=\"%s.treeModel.internalNodeHeights\"/>\n" % g)
    out_xml.write("\t\t</uniformOperator>\n")
out_xml.write("\t\t<scaleOperator scaleFactor=\"0.9\" weight=\"3\">\n\
\t\t\t<parameter idref=\"demographic.populationMean\"/>\n\
\t\t</scaleOperator>\n\
\t\t<sampleNonActiveOperator weight=\"40\">\n\
\t\t\t<distribution>\n\
\t\t\t\t<parameter idref=\"demographic.populationMeanDist\"/>\n\
\t\t\t</distribution>\n\
\t\t\t<data>\n\
\t\t\t\t<parameter idref=\"demographic.popSize\"/>\n\
\t\t\t</data>\n\
\t\t\t<indicators>\n\
\t\t\t\t<parameter idref=\"demographic.indicators\"/>\n\
\t\t\t</indicators>\n\
\t\t</sampleNonActiveOperator>\n\
\t\t<bitFlipOperator weight=\"100\">\n\
\t\t\t<parameter idref=\"demographic.indicators\"/>\n\
\t\t</bitFlipOperator>\n\
\t\t<scaleOperator scaleFactor=\"0.5\" weight=\"60\">\n\
\t\t\t<parameter idref=\"demographic.popSize\"/>\n\
\t\t\t<indicators pickoneprob=\"1.0\">\n\
\t\t\t\t<parameter idref=\"demographic.indicators\"/>\n\
\t\t\t</indicators>\n\
\t\t</scaleOperator>\n")
for h in alignmentNames:
    out_xml.write("\t\t<upDownOperator scaleFactor=\"0.75\" weight=\"3\">\n")
    out_xml.write("\t\t\t<up>\n")
    out_xml.write("\t\t\t\t<parameter idref=\"%s.clock.rate\"/>\n" % h)
    out_xml.write("\t\t\t</up>\n")
    out_xml.write("\t\t\t<down>\n")
    out_xml.write("\t\t\t\t<parameter \
idref=\"%s.treeModel.allInternalNodeHeights\"/>\n" % h)
    out_xml.write("\t\t\t</down>\n")
    out_xml.write("\t\t</upDownOperator>\n")
out_xml.write("\t</operators>\n")

# write MCMC to xml

out_xml.write("\t<mcmc id=\"mcmc\" chainLength=\"100000000\" \
autoOptimize=\"true\" operatorAnalysis=\"%s.ops\">\n" % prefix)
out_xml.write("\t\t<posterior id=\"posterior\">\n")
out_xml.write("\t\t\t<prior id=\"prior\">\n")
# for k in alignmentNames:
#     out_xml.write("\t\t\t\t<gammaPrior shape=\"0.05\" scale=\"10.0\" \
# offset=\"0.0\">\n")
#     out_xml.write("\t\t\t\t\t<parameter idref=\"%s.ac\"/>\n" % k)
#     out_xml.write("\t\t\t\t</gammaPrior>\n")
#     out_xml.write("\t\t\t\t<gammaPrior shape=\"0.05\" scale=\"20.0\" \
# offset=\"0.0\">\n")
#     out_xml.write("\t\t\t\t\t<parameter idref=\"%s.ag\"/>\n" % k)
#     out_xml.write("\t\t\t\t</gammaPrior>\n")
#     out_xml.write("\t\t\t\t<gammaPrior shape=\"0.05\" scale=\"10.0\" \
# offset=\"0.0\">\n")
#     out_xml.write("\t\t\t\t\t<parameter idref=\"%s.at\"/>\n" % k)
#     out_xml.write("\t\t\t\t</gammaPrior>\n")
#     out_xml.write("\t\t\t\t<gammaPrior shape=\"0.05\" scale=\"10.0\" \
# offset=\"0.0\">\n")
#     out_xml.write("\t\t\t\t\t<parameter idref=\"%s.cg\"/>\n" % k)
#     out_xml.write("\t\t\t\t</gammaPrior>\n")
#     out_xml.write("\t\t\t\t<gammaPrior shape=\"0.05\" scale=\"10.0\" \
# offset=\"0.0\">\n")
#     out_xml.write("\t\t\t\t\t<parameter idref=\"%s.gt\"/>\n" % k)
#     out_xml.write("\t\t\t\t</gammaPrior>\n")
#     out_xml.write("\t\t\t\t<uniformPrior lower=\"0.0\" upper=\"1.0\">\n")
#     out_xml.write("\t\t\t\t\t<parameter idref=\"%s.frequencies\"/>\n" % k)
#     out_xml.write("\t\t\t\t</uniformPrior>\n")
#     out_xml.write("\t\t\t\t<exponentialPrior mean=\"0.5\" offset=\"0.0\">\n")
#     out_xml.write("\t\t\t\t\t<parameter idref=\"%s.alpha\"/>\n" % k)
#     out_xml.write("\t\t\t\t</exponentialPrior>\n")
#     out_xml.write("\t\t\t\t<uniformPrior lower=\"0.0\" upper=\"1.0\">\n")
#     out_xml.write("\t\t\t\t\t<parameter idref=\"%s.pInv\"/>\n" % k)
#     out_xml.write("\t\t\t\t</uniformPrior>\n")

for k in alignmentNames:
    out_xml.write("\t\t\t\t<logNormalPrior mean=\"1.0\" stdev=\"1.25\" \
offset=\"0.0\" meanInRealSpace=\"false\">\n")
    out_xml.write("\t\t\t\t\t<parameter idref=\"%s.kappa\"/>\n" % k)
    out_xml.write("\t\t\t\t</logNormalPrior>\n")
    out_xml.write("\t\t\t\t<uniformPrior lower=\"0.0\" upper=\"1.0\">\n")
    out_xml.write("\t\t\t\t\t<parameter idref=\"%s.frequencies\"/>\n" % k)
    out_xml.write("\t\t\t\t</uniformPrior>\n")
    out_xml.write("\t\t\t\t<exponentialPrior mean=\"0.5\" offset=\"0.0\">\n")
    out_xml.write("\t\t\t\t\t<parameter idref=\"%s.alpha\"/>\n" % k)
    out_xml.write("\t\t\t\t</exponentialPrior>\n")
    out_xml.write("\t\t\t\t<uniformPrior lower=\"0.0\" upper=\"1.0\">\n")
    out_xml.write("\t\t\t\t\t<parameter idref=\"%s.pInv\"/>\n" % k)
    out_xml.write("\t\t\t\t</uniformPrior>\n")

for l in alignmentNames:
    out_xml.write("\t\t\t\t<uniformPrior lower=\"0.0\" upper=\"5.0E18\">\n")
    out_xml.write("\t\t\t\t\t<parameter idref=\"%s.clock.rate\"/>\n" % l)
    out_xml.write("\t\t\t\t</uniformPrior>\n")

out_xml.write("\t\t\t\t<poissonPrior mean=\"0.693147\" offset=\"0.0\">\n\
\t\t\t\t\t<statistic idref=\"demographic.populationSizeChanges\"/>\n\
\t\t\t\t</poissonPrior>\n\
\t\t\t\t<oneOnXPrior>\n\
\t\t\t\t\t<parameter idref=\"demographic.populationMean\"/>\n\
\t\t\t\t</oneOnXPrior>\n\
\t\t\t\t<coalescentLikelihood idref=\"coalescent\"/>\n\
\t\t\t\t<mixedDistributionLikelihood>\n\
\t\t\t\t\t<distribution0>\n\
\t\t\t\t\t\t<exponentialDistributionModel \
idref=\"demographic.populationMeanDist\"/>\n\
\t\t\t\t\t</distribution0>\n\
\t\t\t\t\t<distribution1>\n\
\t\t\t\t\t\t<exponentialDistributionModel \
idref=\"demographic.populationMeanDist\"/>\n\
\t\t\t\t\t</distribution1>\n\
\t\t\t\t\t<data>\n\
\t\t\t\t\t\t<parameter idref=\"demographic.popSize\"/>\n\
\t\t\t\t\t</data>\n\
\t\t\t\t\t<indicators>\n\
\t\t\t\t\t\t<parameter idref=\"demographic.indicators\"/>\n\
\t\t\t\t\t</indicators>\n\
\t\t\t\t</mixedDistributionLikelihood>\n\
\t\t\t</prior>\n")
out_xml.write("\t\t\t<likelihood id=\"likelihood\">\n")
for m in alignmentNames:
    out_xml.write("\t\t\t\t<treeLikelihood idref=\"%s.treeLikelihood\"/>\n" % m)

out_xml.write("\t\t\t</likelihood>\n")
out_xml.write("\t\t</posterior>\n")
out_xml.write("\t\t<operators idref=\"operators\"/>\n")


# write logs to xml

out_xml.write("\n")
out_xml.write("\t\t<log id=\"screenLog\" logEvery=\"10000\">\n\
\t\t\t<column label=\"Posterior\" dp=\"4\" width=\"12\">\n\
\t\t\t\t<posterior idref=\"posterior\"/>\n\
\t\t\t</column>\n\
\t\t\t<column label=\"Prior\" dp=\"4\" width=\"12\">\n\
\t\t\t\t<prior idref=\"prior\"/>\n\
\t\t\t</column>\n\
\t\t\t<column label=\"Likelihood\" dp=\"4\" width=\"12\">\n\
\t\t\t\t<likelihood idref=\"likelihood\"/>\n\
\t\t\t</column>\n\
\t\t</log>\n")

out_xml.write("\t\t<log id=\"fileLog\" logEvery=\"10000\" fileName=\"%s.log\" \
overwrite=\"false\">\n" % prefix)
out_xml.write("\t\t\t<posterior idref=\"posterior\"/>\n")
out_xml.write("\t\t\t<prior idref=\"prior\"/>\n")
out_xml.write("\t\t\t<likelihood idref=\"likelihood\"/>")
for o in alignmentNames:
    out_xml.write("\t\t\t<parameter idref=\"%s.treeModel.rootHeight\"/>\n" % o)
out_xml.write("\t\t\t<sumStatistic \
idref=\"demographic.populationSizeChanges\"/>\n")
out_xml.write("\t\t\t<parameter idref=\"demographic.populationMean\"/>\n")
out_xml.write("\t\t\t<parameter idref=\"demographic.popSize\"/>\n")
out_xml.write("\t\t\t<parameter idref=\"demographic.indicators\"/>\n")
for p in alignmentNames:
#    out_xml.write("\t\t\t<parameter idref=\"%s.ac\"/>\n" % p)
#    out_xml.write("\t\t\t<parameter idref=\"%s.ag\"/>\n" % p)
#    out_xml.write("\t\t\t<parameter idref=\"%s.at\"/>\n" % p)
#    out_xml.write("\t\t\t<parameter idref=\"%s.cg\"/>\n" % p)
    out_xml.write("\t\t\t<parameter idref=\"%s.kappa\"/>\n" % p)
    out_xml.write("\t\t\t<parameter idref=\"%s.frequencies\"/>\n" % p)
    out_xml.write("\t\t\t<parameter idref=\"%s.alpha\"/>\n" % p)
    out_xml.write("\t\t\t<parameter idref=\"%s.pInv\"/>\n" % p)
for q in alignmentNames:
    out_xml.write("\t\t\t<parameter idref=\"%s.clock.rate\"/>\n" % q)
for r in alignmentNames:
    out_xml.write("\t\t\t<treeLikelihood idref=\"%s.treeLikelihood\"/>\n" % r)
out_xml.write("\t\t\t<coalescentLikelihood idref=\"coalescent\"/>\n")
out_xml.write("\t\t</log>\n")
for s in alignmentNames:
    out_xml.write("\t\t<logTree id=\"%s.treeFileLog\" logEvery=\"10000\" \
nexusFormat=\"true\" fileName=\"%s.%s.trees\" sortTranslationTable=\"true\">\n"
% (s,prefix,s))
    out_xml.write("\t\t\t<treeModel idref=\"%s.treeModel\"/>\n" % s)
    out_xml.write("\t\t\t<trait name=\"rate\" tag=\"%s.rate\">\n" % s)
    out_xml.write("\t\t\t\t<strictClockBranchRates idref=\"%s.branchRates\"/>\n"
                    % s)
    out_xml.write("\t\t\t</trait>\n")
    out_xml.write("\t\t\t<posterior idref=\"posterior\"/>\n")
    out_xml.write("\t\t</logTree>\n")
out_xml.write("\t</mcmc>\n")

# write .csv info to xml file
out_xml.write("\n")
out_xml.write("\t<report>\n\
\t\t<property name=\"timer\">\n\
\t\t\t<mcmc idref=\"mcmc\"/>\n\
\t\t</property>\n\
\t</report>\n")
out_xml.write("\t<VDAnalysis id=\"demographic.analysis\" burnIn=\"0.1\" \
useMidpoints=\"true\">\n")
out_xml.write("\t\t<logFileName>\n\t\t\t%s.log\n\t\t</logFileName>\n" % prefix)
out_xml.write("\t\t<treeFileNames>\n")
for t in alignmentNames:
    out_xml.write("\t\t\t<treeOfLoci>\n")
    out_xml.write("\t\t\t\t%s.%s.trees\n" % (prefix,t))
    out_xml.write("\t\t\t</treeOfLoci>\n")
out_xml.write("\t\t</treeFileNames>\n")
out_xml.write("\t\t<populationModelType>\n")
out_xml.write("\t\t\tlinear\n")
out_xml.write("\t\t</populationModelType>\n")
out_xml.write("\t\t<populationFirstColumn>\n")
out_xml.write("\t\t\tdemographic.popSize1\n")
out_xml.write("\t\t</populationFirstColumn>\n")
out_xml.write("\t\t<indicatorsFirstColumn>\n")
out_xml.write("\t\t\tdemographic.indicators1\n")
out_xml.write("\t\t</indicatorsFirstColumn>\n")
out_xml.write("\t</VDAnalysis>\n")
out_xml.write("\t<CSVexport fileName=\"%s.csv\" separator=\",\">\n" % prefix)
out_xml.write("\t\t<columns>\n")
out_xml.write("\t\t\t<VDAnalysis idref=\"demographic.analysis\"/>\n")
out_xml.write("\t\t</columns>\n")
out_xml.write("\t</CSVexport>\n")
out_xml.write("</beast>\n")

out_xml.close()
