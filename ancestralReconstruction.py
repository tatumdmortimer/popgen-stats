#!/usr/bin/python

import sys, os, argparse, subprocess, shlex, glob
from datetime import datetime
from multiprocessing import cpu_count
from multiprocessing.dummy import Pool as ThreadPool
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

################################################################################
# This script finds core genes from a a list of genomes and a groups file in
# the format output by OrthoMCL. It aligns the core genes and makes a tree from
# the alignment. It uses Lazarus to run PAML ancestral reconstruction, and it   
# concatenates the ancestral gene sequences.
# 
# Program Requirements: translatorX, mafft, lazarus, raxml, biopython
# Input: OrthoMCL groups file, list of genomes (outgroup last), nucleotide
#       sequences for genes in the genomes
################################################################################


TRANSLATOR_X_PATH = "/opt/PepPrograms/translatorx_vLocal.pl"
LAZARUS_PATH = "/opt/PepPrograms/project-lazarus/lazarus.py"
RAXML_PATH = "/opt/PepPrograms/standard-RAxML/raxmlHPC-PTHREADS-AVX"

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
            os.path.abspath(os.path.expanduser(values)))

def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

def is_dir(dirname):
    """Checks if a path is a directory"""
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname

def is_file(filename):
    """Checks if a file exists"""
    if not os.path.isfile(filename):
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename

def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Ancestral reconstruction of\
 core genome')
    parser.add_argument("groups", help="OrthoMCL groups file", action=FullPaths,
        type=is_file)
    parser.add_argument("genomes", help="File listing genomes to be included in\
 the analysis- outgroup last", action=FullPaths, type=is_file)
    parser.add_argument("genes", 
        help="Directory with .fasta files of nucleotide sequences for genomes",
        action=FullPaths, type=is_dir)
    parser.add_argument("-t", "--threads", 
        help="Number of threads to use (default: 2)",
        type=int, default=2, choices=range(2, cpu_count()))
    return parser.parse_args()

def check_paths():
    for i in [TRANSLATOR_X_PATH, LAZARUS_PATH, RAXML_PATH]:
        if not os.path.isfile(i):
            msg = "{0} does not exist".format(i)
            print msg
            sys.exit()

def call_with_log(cmd):
    """Calls a system command with the subprocess module. Redirects both stdout
    and stderr to a log file"""
    cmd = cmd.format(**(kvmap))
    logfile = open(wd + current_datetime+".log", "a+")
    logfile.write("Executing command: " + cmd + "\n")
    logfile.flush()
    ret = subprocess.call(shlex.split(cmd), stdout=logfile, stderr=logfile)
    if(ret != 0):
        print("Pipeline did not complete successfully. \n Command : \n\n" + 
            cmd + "\n\n returned with non-zero code: " + str(ret))
    logfile.close()

def read_groups_file(inFileName):
    """ Read in groups file and create dictionary of group name and proteins in
    group"""
    print "Reading groups file"
    inFile = open(inFileName, 'r')
    groupsDict = {}
    for line in inFile:
        line = line.strip()
        entries = line.split(':')
        groupName = entries[0]
        groupProteins = entries[1][1:].split(' ')
        groupsDict[groupName] = groupProteins
    inFile.close()
    return groupsDict

def get_core_genes(groupsDict, genomes):
    """ Gets core genes for genomes in list """
    coreGenes = set()
    for group in groupsDict:
        genomeList = []
        proteinList = groupsDict[group]
        for protein in proteinList:
            ids = protein.split('|')
            genomeID = ids[0]
            genomeList.append(genomeID)
        genomeSet = set(genomeList)
        if set(genomes.keys()).issubset(genomeSet):
            if len(genomeList) == len(genomeSet):
                coreGenes.add(group)
    return coreGenes

def make_unaligned_fasta(dnaDirectory, groupsDict, coreGenes, genomes, og):
    """ Reads through files in provided directory to find gene sequences that
    match the proteins in the groups dictionary"""
    print "Collecting core genes"
    def make_fasta(group):
        proteins = groupsDict[group]
        out = open(group + '/' + group + '.fasta', 'w')
        records = []
        outgroup_gene = og
        ingroup_genes = []
        for protein in proteins:
            seqID = protein.split('|')[0]
            if seqID in genomes:
                protein = protein.split('|')[1]
                newRec = seqRecordDict[protein]
                newRec.description = ""
                records.append(newRec)
                if og in newRec.id:
                    outgroup_gene = newRec.id
                else:
                    ingroup_genes.append(newRec.id)
        SeqIO.write(records, out, 'fasta')
        return (group, ingroup_genes, outgroup_gene)
    files = listdir_fullpath(dnaDirectory)
    seqRecordDict = {}
    seqIDs = []
    for f in files:
        handle = open(f, 'r')
        for record in SeqIO.parse(handle, 'fasta'):
            seqRecordDict[record.id] = record
    pool = ThreadPool(args.threads)
    seqIDs = pool.map(make_fasta, coreGenes)
    pool.close()
    pool.join()
    return seqIDs

def align_gene_sequences(coreGenes):
    """ Use MAFFT to align gene sequences"""
    print "Aligning core genes"
    def run_translatorX(infile):
        outfile = "%s/alignment/%s" % (os.path.dirname(infile), 
            os.path.splitext(os.path.split(infile)[1])[0])
        call_with_log(TRANSLATOR_X_PATH + " -i %s -o %s -p F"
            % (infile, outfile))
        return outfile + ".nt_ali.fasta"
    files = [x + "/" + x + ".fasta" for x in coreGenes]
    pool = ThreadPool(args.threads)
    outfileList = pool.map(run_translatorX, files)
    pool.close()
    pool.join()

def make_trees(coreGenes):
    """ Use RAxML to calculate maximum likelihood phylogeny """
    print "Running RAxML"
    def run_raxml(coreGene):
        alignFile = "%s/alignment/%s.nt_ali.fasta" % (coreGene, coreGene)
        outdir = os.path.abspath("%s/tree/" % coreGene)
        name = "ml_" + os.path.splitext(os.path.basename(alignFile))[0]
        call_with_log(RAXML_PATH + " -T 2 -m GTRGAMMA -# 20 -p 123 -s %s -w %s \
-n %s" % (alignFile, outdir, name))
        return name
    pool = ThreadPool(args.threads/2)
    outNames = pool.map(run_raxml, coreGenes)
    pool.close()
    pool.join()

def ancestral_reconstruction(outgroup_genes):
    """ Use Lazarus wrapper to run paml ancestral reconstruction """
    print "Ancestral Reconstruction"
    def run_lazarus(outgroup_gene):
        coreGene, ingroup, outgroup = outgroup_gene
        align = "%s/alignment/%s.nt_ali.fasta" % (coreGene, coreGene)
        tree = "%s/tree/RAxML_bestTree.ml_%s.nt_ali" % (coreGene, coreGene)
        model = "/opt/PepPrograms/paml4.8/dat/wag.dat"
        outdir = os.path.abspath("%s/ancestral/" % coreGene)
        call_with_log(LAZARUS_PATH + " --codeml --outputdir %s --verbose 9 \
--alignment %s --tree %s --model %s --asrv 4 --getanc --ingroup %s --outgroup \
%s" % (outdir, align, tree, model, "[%s]" % (",".join(ingroup)), 
"[%s]" % (outgroup)))
    pool = ThreadPool(args.threads)
    pool.map(run_lazarus, outgroup_genes) 
    pool.close()
    pool.join()


def concatenate_genes(coreGenes):
    """ Using Biopython, concatenate ancestral genes into one file"""
    print "Concatenating ancestral reconstructions"
    def parse_lazarus_output(coreGene):
        ancRecFile = open("%s/ancestral/ancestor.out.txt" % coreGene, "r")
        for i,line in enumerate(ancRecFile):
            if i == 13:
                record = SeqRecord(Seq(line.strip(), IUPAC.ambiguous_dna), 
                    id=coreGene, description="")
        return record
    pool = ThreadPool(args.threads)
    recs = pool.map(parse_lazarus_output, coreGenes)
    SeqIO.write(recs, "ancestralGenes.fa", "fasta")

                

check_paths()
args = get_args()
current_datetime = datetime.today().strftime("%d-%m-%Y-%H%M")
wd = os.getcwd() + "/"
kvmap = {'projectname':'coreAlignment'}
genomes = {}
orderedGenomes = []
with open(args.genomes, "r") as inFile:
    for line in inFile:
        genomes[line.strip()[-4:]] = line.strip()
        orderedGenomes.append(line.strip())
og = orderedGenomes[-1]
groupsDict = read_groups_file(args.groups)
coreGenes = get_core_genes(groupsDict, genomes)
for n in coreGenes:
    try:
        os.mkdir(n)
        os.mkdir(n + "/alignment")
        os.mkdir(n + "/tree")
        os.mkdir(n + "/ancestral")
    except OSError as e:
        if e.errno != os.errno.EEXIST:
            raise
outgroup_genes = make_unaligned_fasta(args.genes, groupsDict, coreGenes, 
    genomes, og)
align_gene_sequences(coreGenes)
make_trees(coreGenes)
ancestral_reconstruction(outgroup_genes)
concatenate_genes(coreGenes)
