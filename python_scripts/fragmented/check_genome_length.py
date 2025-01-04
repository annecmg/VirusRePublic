#!/usr/bin/env python3
"""Author: Devin van Valkengoed

Date: 14-Aug-2023
Description: This script is used to remove gaps from FASTA files and
             determine whether sequences should be included in further
             analyses based on a given sequence length.

Usage: python3 check_genome_length.py -indir <./input_dir/> -min_len <min_len>
       -outdir_included <./output_included/> -outdir_excluded
       <./output_excluded>

Key arguments:
    indir -- str, absolute pathway to the directory containing one or more
             input FASTA files. Only files ending on ".fasta" will be used
             in the analysis.
    min_len -- int, the minimum number of non-X characters a sequence
               should contain to be included.
    outdir_included -- str, absolute pathway to the directory where the
                       FASTA files of included, cleaned genomes should be
                       stored.
    outdir_exluded -- str, absolute pathway to the directory where the
                      FASTA files of exluded, cleaned genomes should be
                      stored.
"""
# Import statements
import argparse
import sys
import glob
import Bio
from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Seq import Seq


# Commandline argument parsing
def parsing_cmd():
    """Uses argparse to determine the given commandline options

    Returns:
        parser -- argparse.Namespace object
    """
    parser = argparse.ArgumentParser(
        description="This script is used to remove gaps from given "
                    "FASTA sequences and determine whether sequences should "
                    "be included in further analyses based on their "
                    "sequence length.\n"
                    "The output of this script is: 1) the cleaned sequences "
                    "with original file name, 2) included_sequences.txt -- "
                    "all file names of sequences that past the length filter, "
                    "3) excluded_sequences.txt -- all file names of "
                    "sequences that did not pass the length filter.")
    parser.add_argument("-indir", help="Str -- Absolute pathway to the "
                                       "directory that contains "
                                       "all input files. All files in this "
                                       "directory with '.fasta' extension "
                                       "will be used as input by the "
                                       "script.",
                        metavar="./input_dir/", dest="input", required=True)
    parser.add_argument("-min_len", help="Int -- The minimum number of non "
                                         "'X' and non gap characters for a "
                                         "sequence to be included in the "
                                         "output.",
                        metavar="min_len", dest="min_len", required=True)
    parser.add_argument("-outdir_included",
                        help="Str -- Directory to which all the output "
                             "sequences and statistic files will be writen "
                             "to that are included based on sequence length.",
                        metavar="./output_dir_inc/", dest="output_inc",
                        required=True)
    parser.add_argument("-outdir_excluded",
                        help="Str -- Directory to which all the output "
                             "sequences and statistic files will be writen "
                             "to that are excluded based on sequence length.",
                        metavar="./output_dir_exc/", dest="output_exc",
                        required=True)

    return parser


def get_file_names(input_dir):
    """Obtains all file names in the given directory that end with '*.fasta'

    Key arguments:
        input_dir -- str, absolute pathway to an input directory.

    Returns:
        file_names -- lst, containing all the absolute pathways to all files
                      the input directory.
    """
    return glob.glob(input_dir + "*.fasta")


def parse_fasta(input_file):
    """Parse a fasta file into a Biopython Seq object

    Key arguments:
        input_file -- str, absolute pathway to the file that contains a
                      fasta sequence.
    Returns:
        fasta_seqs -- lst, containing BioSeq record objects for every FASTA
                      sequences that is found in the input file.
    """
    fasta_seqs = []
    with open(input_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            fasta_seqs.append(record)

    return fasta_seqs


def clean_fasta_sequences(fasta_seq, output_included, output_excluded,
                          file_ext, min_len):
    """Removes gaps from given FASTA sequence and counts non X characters

    Key arguments:
        fasta_seq -- lst, containing Biopython Seq records for every fasta
                     sequence that should be checked.
        output_included -- str, pathway to the directory where the FASTA
                           files of included accessions should be written to.
        output_excluded -- str, pathway to the directory where the FASTA
                           files of excluded accessions should be written to.
        file_ext -- str, general file extension that should be appended to
                    the file names of the cleaned FASTA files.
        min_len -- int, the minimum number of characters in the sequences
                   that should be non-X characters for the sequence to be
                   included.

    Returns:
        excluded_genome_ids -- lst, the record ids of the excluded sequences.
        included_genome_ids -- lst, the record ids of the included sequences.
    """
    excluded_genome_ids = []
    included_genome_ids = []
    for accession in fasta_seq:
        accession_id = accession[0].id.split("_")[0]
        out_exc = output_excluded + accession_id + "_exc_" + \
                  file_ext + ".fasta"
        out_inc = output_included + accession_id + "_inc_" + \
                  file_ext + ".fasta"
        for record in accession:
            # Remove the gaps from the sequences
            record.seq = str(record.seq)
            record.seq = Seq(record.seq.replace("-", ""))

            # Count the number of X characters in the sequences
            x_count = record.seq.count("X")
            seq_len = len(record.seq)
            absolute_len = seq_len - x_count

            if absolute_len < min_len:
                # The sequence is excluded based on its length
                excluded_genome_ids.append(record.id)
                with open(out_exc, "a+") as excluded_out:
                    SeqIO.write(record, excluded_out, "fasta")
            else:
                # The sequence is longer than the minimum non-X length
                included_genome_ids.append(record.id)
                with open(out_inc, "a+") as included_out:
                    SeqIO.write(record, included_out, "fasta")

    return excluded_genome_ids, included_genome_ids


def main():
    """Main function if this script that wraps all the code """
    # Step 0: parse cmd line arguments
    parser = parsing_cmd()
    args = parser.parse_args()

    # Step 1: obtain the pathways to all files in the given input directory
    file_paths = get_file_names(args.input)

    # Step 2: open the files and parse the sequences into a Biopython object
    fasta_sequences = []
    for file in file_paths:
        fasta_sequences.append(parse_fasta(file))

    # Step 3: clean the fasta sequences from gaps and count the non-X char
    excluded_ids, included_ids = clean_fasta_sequences(fasta_sequences,
                                                       args.output_inc,
                                                       args.output_exc,
                                                       "patched_genomes",
                                                       int(args.min_len))

    # Step 4: write the included and exclude genome ids to the output files
    output_excluded_file = args.output_exc + "excluded_genome_ids.txt"
    output_included_file = args.output_inc + "included_genome_ids.txt"

    for inc_id in included_ids:
        with open(output_included_file, "a+") as outf_inc:
            outf_inc.write(inc_id + "\n")

    for exc_id in excluded_ids:
        with open(output_excluded_file, "a+") as outf_exc:
            outf_exc.write(exc_id + "\n")


if __name__ == "__main__":
    main()
