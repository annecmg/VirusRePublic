#!/usr/bin/env python3
""" Author: Devin van Valkengoed

Date: 10-Jan-2023
Description: This script is a helper module for DIAMOND_filter.py and
             includes the following functions:
             *
Dependencies:
    * EMBOSS getorf v6.6.0.0
"""
# Import statements
import sys
import subprocess
import os


# Functions
def assembly_parser(spades_assembly, node):
    """Obtains the desired sequence of given node from a SPAdes assembly file

    Key arguments:
        spades_assembly -- str, pathway to a SPAdes assembly file, containing
                           sequences of assembled contigs in either DNA or
                           RNA format.
        node -- str, the name of the desired node of which the sequence should
                be filtered from the file.

    Returns:
        header -- str,
        clean_fasta -- str,
    """
    cmd_parse_sequence = "sed -n '/{}/,/>/p' {} | " \
                         "head -n -1 | tail -n +2".format(node,
                                                          spades_assembly)
    cmd_parse_header = "sed -n '/{}/,/>/p' {} | " \
                       "head -n 1".format(node, spades_assembly)

    sequence = subprocess.check_output(cmd_parse_sequence,
                                       shell=True).decode("UTF-8").strip()
    header = subprocess.check_output(cmd_parse_header,
                                     shell=True).decode("UTF-8").strip()

    # Remove the newlines in the fasta sequence, for translation with Biopython
    clean_fasta = ""
    for line in sequence.split("\n"):
        clean_fasta += line.strip()

    return header, clean_fasta


def write_line(line_list, output_file):
    """Function that writes lines given in a list to the output file

    Key Arguments:
        line -- lst, containing one or multiple strings.
        output_file -- str, pathway to the file to which the line should be
                       appended.

    Returns:
        None
    """
    with open(output_file, 'a+') as output:
        for line in line_list:
            output.write(str(line) + "\n")

    return


def obtain_orf(input_file, output_file, min_orf):
    """Obtains the Open Reading Frame(s) of a given nucleotide sequence

    Key arguments:
        input_file -- str, pathway to the file containing a nucleotide sequence
                      in fasta format.
        output_file -- str, pathway to the file where the protein sequence
                       should be written to.
        min_orf -- int, minimum length in NT for the orfs to be included.
                   Note: this limit should be set in such a way that only 1
                   orf is returned by EMBOSS getorf.

    Returns:
        None
    """
    if not os.path.exists(output_file):
        print(input_file)
        cmd_get_orf = "getorf -sequence {} -minsize {} " \
                      "-outseq {}".format(input_file, str(min_orf),
                                          output_file)
        subprocess.check_output(cmd_get_orf, shell=True)

    return


def getorf_fasta_parser(fasta_file, nt_file, prefix, suffix):
    """Reads a fasta sequence into memory and removes the input file and header

    Key arguments:
        prefix -- str, something to prefix to the header of the fasta
        fasta_file -- str, pathway to the fasta file that should be read.
        suffix -- str, something to append to the header line.
        nt_file -- str, pathway to the file containing the original nucleotide
                   sequence that should also be deleted.

    Returns:
        fasta_sequence -- str, either a NT or AA sequences, parsed from the
                          given file in fasta format.
    """
    nr_orfs = 0
    sequence = ""
    with open(fasta_file, 'r') as input_file:
        for line in input_file:
            line = line.strip()
            if not line:
                continue
            elif line.startswith(">"):
                nr_orfs += 1
                header = "_".join(line.split(" ")[0].split("_")[:-1])
                if prefix:
                    header = ">{}_{}".format(prefix, header[1:])
                if suffix:
                    header = "{}_{}".format(header, suffix)
            else:
                sequence += line

    fasta_sequence = "{}\n{}".format(header, sequence)

    # Remove the original files
    os.remove(fasta_file)
    os.remove(nt_file)

    # Check if only 1 ORF is found
    if nr_orfs != 1:
        raise ValueError("WARNING: please note that the given minimum ORF "
                         "length gives multiple or no ORFs as output!")

    return fasta_sequence


if __name__ == "__main__":
    """This script is a helper module for DIAMOND_filter.py """
    pass

    # Old code to make the script run by itself
    # fasta_header, fasta_seq = assembly_parser(sys.argv[1], sys.argv[2])
    #
    # write_line([fasta_header, fasta_seq], "test_nt_out.del")
    #
    # obtain_orf("test_nt_out.del", "test_protein_out.del", 100)
    #
    # full_fasta = getorf_fasta_parser("test_protein_out.del", "test_nt_out.del", "DRR086508", sys.argv[3])
    # print(full_fasta)
