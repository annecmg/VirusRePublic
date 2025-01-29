#!/usr/bin/env python3
"""Author: Devin van Valkengoed, Anne Kupczok

Date: 29-Jan-2025
Description: This script is used to obtain genomes from the original
             rnaviralSPAdes assembly files.

Usage: python3 obtain_complete_genomes.py -input_accessions <input_tsv>
       -nt_output <nuc_out>

Key arguments:
    input_tsv -- str, absolute pathway to a tsv file containing information
                  about the desired sequences. The columns should be as
                  followed: accession_id \t NODE_name \t optional rest.
    input_folder -- str, assembly folder for the pipeline, typically:
                  output/assembly/rnaviralspades/
    nuc_output -- str, absolute pathway to the folder where all output
                  nucleotide FASTA files should be written to.

"""
# Import statements
import argparse
import os
import sys
import subprocess
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO


# Commandline argument parsing
def parsing_cmd():
    """Uses argparse to determine the given commandline options

    Returns:
        parser -- argparse.Namespace object
    """
    parser = argparse.ArgumentParser(
        description="This script is used to obtain nucleotide sequence of "
                    "interest from rnaviralSPAdes assembly files.")
    parser.add_argument("-input_accessions",
                        help="Str -- Absolute pathway to a text file "
                             "containing the file directories, "
                             "fasta header and additional description for "
                             "every accession of interest. The file should "
                             "be tab separated and every new line should "
                             "look like this: "
                             "accession_id \t NODE_name \t rest (optional)",
                        dest="input",
                        required=True)
    parser.add_argument("-input_folder",
                        help="str, assembly folder for the pipeline, typically: "
                                      "output/assembly/rnaviralspades/",
                        dest="input_folder", required=True)
    parser.add_argument("-nt_output",
                        help="Str -- Absolute pathway to the directory where "
                             "all FASTA files containing the complete "
                             "nucleotide sequence should be stored.",
                        dest="nt_output",
                        required=True)

    return parser

def get_record(contigf,node):

    for rec in SeqIO.parse(contigf,"fasta"):
        if rec.id==node:
            return rec

def main():
    """Main function that wraps all the code in this python script"""
    # Step 0: parse the cmd line arguments
    parser = parsing_cmd()
    args = parser.parse_args()

    # Step 1: obtain the location of all assembly using the input file
    input_file=args.input
    input_folder=args.input_folder
    output=args.nt_output

    outf=open(output,'w')
    for line in open(input_file):
        spl=line.split()
        acc=spl[0]
        node=spl[1]
        rec=get_record("{}/{}/contigs.fasta".format(input_folder,acc),node)
        rec.id="{}_{}".format(acc,rec.id)
        SeqIO.write(rec,outf,"fasta")

if __name__ == "__main__":
    main()
