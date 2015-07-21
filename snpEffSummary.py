#!/usr/bin/python

import sys, os, argparse
import pandas
import vcf
from datetime import datetime

################################################################################
# This script parse output of snpEff (annotated vcf and genes file) and produces
# a summary of the synonymous and nonsynonymous variants
################################################################################

ORTHOMCL_PATH = "/opt/PepPrograms/orthomclSoftware-v2.0.9/bin/"

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
    parser = argparse.ArgumentParser(description='Summarize snpEff output')
    parser.add_argument("vcf", help="annotated vcf file", action=FullPaths,
        type=is_file)
    parser.add_argument("gene", help="genes output file", action=FullPaths,
        type=is_file)
    return parser.parse_args()

args = get_args()  

def summarize_genes(genesFile):
    genes = pandas.read_csv(genesFile, sep = "\t", skiprows=1)
    genes = genes[["#GeneName", "GeneId", "variants_effect_missense_variant",
        "variants_effect_synonymous_variant", "variants_effect_stop_gained", 
        "variants_effect_stop_lost+splice_region_variant"]]
    outGenesFile = os.path.splitext(genesFile)[0] + "_reduced.txt"
    genes.to_csv(path_or_buf=outGenesFile, sep = "\t")

def summarize_vcf(vcf_file):
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    outVcfFile = os.path.splitext(vcf_file)[0] + "_vcf_summary.txt"
    out = open(outVcfFile, "w")
    out.write("Position\tRef\tAlt\tGene\tType\n")
    for record in vcf_reader:
        if record.ALT[0] == None:
            continue
        out.write("{0}\t{1}\t{2}\t".format(str(record.POS), record.REF,
            record.ALT[0]))
        for i in record.INFO['ANN']:
            if "missense_variant" in i:
                gene = i.split(",")[0].split("|")[4]
                out.write(gene + "\tnon-synonymous\n")
                break
            elif "synonymous_variant" in i:
                gene = i.split(",")[0].split("|")[4]
                out.write(gene + "\tsynonymous\n")
                break
            elif "intergenic" in i:
                out.write("-\tintergenic\n")
                break
            elif "stop_gained" in i:
                gene = i.split(",")[0].split("|")[4]
                out.write(gene + "\tstop-gained\n")
                break
            elif "stop_lost" in i:
                gene = i.split(",")[0].split("|")[4]
                out.write(gene + "\tstop-lost\n")
                break

    out.close()



summarize_genes(args.gene)
summarize_vcf(args.vcf)
