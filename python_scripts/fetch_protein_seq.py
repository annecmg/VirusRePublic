#!/usr/bin/env python3
"""Author: Devin van Valkengoed

Date: 16-Feb-2023 
Description: This script can be used to obtain protein FASTA sequences when
             provided with a file containing NCBI protein identifiers on
             every new line.

Usage: python3 fetch_protein_seq.py <sequences_in.txt>

    Key arguments:
        sequences_in -- str, pathway to a file with NCBI protein IDs on
                        every new line.

Dependencies:
    * NCBI commandline tools v16.2 (esearch and efetch)
"""
# Import statements
import sys
import subprocess


def read_lines(in_file):
    """Reads the lines of a file into memory

    Key arguments:
        in_file -- str, pathway to a file with sequence accession on every
                   new line.

    Returns:
        sequences -- lst, all new lines of the file in a list.
    """
    sequences = []
    with open(in_file, 'r') as seq_file:
        for line in seq_file:
            sequences.append(line.strip())

    return sequences


def fetch_protein_sequence(identifiers):
    """Obtains the protein sequences in fasta of the given identifiers

    Key arguments:
        identifiers -- lst, protein identifiers to obtain the fasta sequence
                       for.

    Returns:
        None (prints the protein sequences to stdout)
    """
    for accession in identifiers:
        # Check for UniProtKB/SwissProt identifiers
        if "sp" in accession:
            accession = accession.split("|")[1]
        cmd_protein_fetch = "esearch -db protein -query {} | " \
                            "efetch -format fasta".format(accession)

        print(subprocess.check_output(cmd_protein_fetch,
                                      shell=True).decode("UTF-8").strip())

    return


def main():
    """Main function that connects all the code in this script """
    # Step 1: parse the accessions from the given file
    seq = read_lines(sys.argv[1])

    # Step 2: use ncbi commandline tools to find the desired sequences
    fetch_protein_sequence(seq)


if __name__ == "__main__":
    main()
