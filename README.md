# popgen-stats

Scripts to calculate population genetics statistics

### selectionStats.py
This script calculates diversity and selection statistics for fasta alignments. Alignments must represent in frame coding sequences. The script will calculate theta, pi, piN, piS, and Tajima's D. If an outgroup(s) is provided, results of the McDonald Kreitman test will also be output. A single fasta alignment or a directory of fasta alignments can be provided to the program.

Requirements: EggLib with Bio++ (http://egglib.sourceforge.net/)

Current Versions: Python v 2.7.3, EggLib v 2.1.7, Bio++ v 2.1.0

Usage:
selectionStats.py -a [alignment] -d [directory] -o ["list, of, outgroups"]

### selectionStatsPlot.R
This script plots the results of selectionStats.py

Requirements: ggplot2 (http://ggplot2.org/)

Current Versions: R v 3.0.2, ggplot2 v 0.9.3.1

Usage:
selectionStatsPlot.R [input filename]

### slidingWindowStats.py
This calculates theta, pi, and Tajima's D in a sliding window across an alignment. A default window size of 1000 bp and window step of 300 bp is used unless otherwise specified by the user.

Requirements: EggLib 

Current Versions: Python 2.7.6, EggLib 2.1.7

Usage:
slidingWindowStats.py -a [alignment] -w [window width (default 1000)] -s [window step (default 300)]

### genewisePAML.py
This script runs codeml analysis from PAML on an alignment or directory of alignments. It compares the likelihood of the M1a (nearly neutral) and M2a (positive selection) models. Genes with a significant Likelihood Ratio Test are output to a file called pamlResults.txt.

Requirements: EggLib, PAML, Biopython. Phylip

Current Versions: Python 2.7.6, EggLib 2.1.7, Biopython 1.63, PAML 4.8, Phylip 3.6

Usage:
genewisePAML.py -a [alignment] -d [directory]

### BSP.R
This script plots the Bayesian Skyline Plot data that can be exported from
Tracer

Requirements: ggplot2

Current Versions: R 3.1.0, ggplot2 1.0.0

Usage:
BSP.R [BSP data file]

### createEBSPxml.py
This script creates the input xml file for BEAST v 1.8 from a directory of 
fasta alignments for an Extended Bayesian Skyline Plot analysis. Currently uses
a strict clock and GTR+Gamma+I subsitution model.

Requirements: Biopython

Current Versions: Python 2.7.6, Biopython 1.63

Usage:
createEBSPxml.py [directory] [prefix for files]

### LRT.py
This script performs a likelihood ratio test given two likelihoods and degrees
of freedom.

Requirements: Biopython

Current Versions: Python 2.7.6, Biopython 1.63

Usage:
LRT.py [likelihood 1] [likelihood 2] [degrees of freedom]

### dadiAnalysis.py
This script uses dadi to perform demographic analysis on a site frequency 
spectrum. Test standard neutral model, expansion model, exponential growth
model, bottleneck model, and bottlegrowth model. 

Requirements: Biopython, NumPy, dadi (code.google.com/p/dadi)

Current Versions: Python 2.7.6, Biopython 1.63, NumPy 1.6.1, dadi 1.6.3

Usage:
dadiAnalysis.py [input SFS]

### plot_dadi_SFS.R
Plots the SFS results from dadiAnalysis.py

Requirements: ggplot2

Current Versions: R 3.1.0, ggplot2 1.0.0

Usage:
plot_dadi_SFS.R

### siteFrequencySpectra.py
Reads vcf that has been annotated with snpEff. The outgroup should be the
reference sequence in the VCF. Outputs synonymous, nonsynonymous, and combined
site frequency spectra.

Current Versions: Python 2.7.6

Usage:
siteFrequencySpectra.py -v [VCF] -n [number of strains without outgroup]

### siteFrequencySpectrum.R
Outputs unfolded site frequency spectrum given an alignment and outgroup.

Requirements: PopGenome, ggplot2

Current Versions: R 3.1.0, ggplot2 1.0.0, PopGenome 2.0.8

Usage:
siteFrequencySpectrum.R [directory containing alignment]

### codonUsage.py
Calculates the CAI (codon adaptation index) of a directory of fasta files and puts all the results in one file.
All the fasta files should contain sequences from the same organisms.The idea is to calculate the CAI for the core genome
of a species, so each fasta file is one gene with multiple strains within it.
Requires a codon usage table created with EMBOSS cusp (http://emboss.sourceforge.net/apps/cvs/emboss/apps/cusp.html).

Requirements: EMBOSS cai (http://emboss.sourceforge.net/apps/cvs/emboss/apps/cai.html)

Current Versions: EMBOSS v 6.5.7

Usage:
codonUsage.py [directory containing .fasta files] [codon usage table]
