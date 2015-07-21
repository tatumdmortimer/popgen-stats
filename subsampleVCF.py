#!/usr/bin/env python

import sys
import os
import argparse
import random
from Bio import AlignIO

# This script randomly samples a VCF & alignment to remove gapped regions and
# create SFS

def get_args():
    # parse command line arguments
    parser = argparse.ArgumentParser(description='Subsamples a VCF to create \
SFS')
    parser.add_argument("vcf", help="VCF to be sampled",
                    type=argparse.FileType('r'))
    parser.add_argument("n", help="Length of output SFS", type=int)
    parser.add_argument("align", help="Alignment to find gaps",
                    type=argparse.FileType('r'))
    parser.add_argument("-p", "--prfreq", help="Create prfreq input files",
                        action="store_true")
    parser.add_argument("-d", "--dadi", help="Create dadi input files",
                        action="store_true")
    return parser.parse_args()
    
def find(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]

def create_vcf_dictionary(vcf):
    vcfDict = {}
    for line in vcf:
        if line[0] == "#":
            continue
        line = line.strip().split()
        pos = line[1]
        ALT = line[4]
        if len(ALT) > 1:    #skip positions that aren't biallelic
            continue
        INFO = line[7]
        outgroup = line[9]
        if outgroup != ".":
            print "VCF in incorrect format."
            print "Outgroup should be reference & first strain in alleles"
            sys.exit()
        alleles = line[10:]
        if "SILENT" in INFO:
            vcfDict[pos] = ["S", alleles]
        if "MISSENSE" in INFO:
            vcfDict[pos] = ["NS", alleles]
    vcf.close()
    return vcfDict

def add_gaps(align, vcfDict):
    alignment = AlignIO.read(align, "fasta")
    for i,seqRecord in enumerate(alignment):
        gapIndex = find(seqRecord.seq, '-')
        for snp in vcfDict:
            if int(snp) - 1 in gapIndex:
                vcfDict[snp][1][i-1] = "-"
    for snp in vcfDict:
        vcfDict[snp][1] = [s for s in vcfDict[snp][1] if s != "-"]
    return vcfDict

def create_sfs(vcfDict, n):
    all_sfs = [0]*(n+1)
    syn_sfs = [0]*(n+1)
    nonSyn_sfs = [0]*(n+1)
    for snp in vcfDict:
        if len(vcfDict[snp][1]) < n:
            continue
        sample = random.sample(vcfDict[snp][1], n)
        freq = len(sample) - sample.count('.') 
        all_sfs[freq] += 1
        if vcfDict[snp][0] == "S":
            syn_sfs[freq] += 1
        else:
            nonSyn_sfs[freq] += 1

    return (all_sfs, syn_sfs, nonSyn_sfs)

def write_dadi(sfs, sfsName):
    dadi_out = open("{0}_dadi.txt".format(sfsName), "w")
    dadi_out.write("# SFS subsampled\n")
    dadi_out.write("# Outgroup: Troy\n")
    dadi_out.write("%i unfolded\n" % (len(sfs)))
    dadi_out.write(" ".join(str(k) for k in sfs) + "\n")
    dadi_out.close()

def write_prfreq(sfs, sfsName):
    prfreq_out = open("{0}_prfreq.txt".format(sfsName), "w")
    for m in sfs[1:-1]:
        prfreq_out.write(str(m) + "\n")

args = get_args()

vcfDict = create_vcf_dictionary(args.vcf)
vcfDict = add_gaps(args.align, vcfDict)
all_sfs, syn_sfs, nonSyn_sfs = create_sfs(vcfDict, args.n)

if args.dadi:
    write_dadi(all_sfs, "sfs")
    write_dadi(syn_sfs, "synonymous_sfs")
    write_dadi(nonSyn_sfs, "nonSynonymous_sfs")
if args.prfreq:
    write_prfreq(all_sfs, "sfs")
    write_prfreq(syn_sfs, "synonymous_sfs")
    write_prfreq(nonSyn_sfs, "nonSynonymous_sfs")
