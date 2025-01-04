#!/usr/bin/env python3
"""Author: Devin van Valkengoed

Date: 23-Dec-2022 
Description: This script is imported as a helper module for main module
'meta_accesion_parser.py' and includes the following functions:

    * determine_kmer -- determines the 'optimal' K-mer values for a SPAdes
                        assembly based on 75% of the read length of a given
                        FastQ file.

Notes:
    * To run this module, package seqkit has to be installed!
"""
# Import statements
import sys
import subprocess


# Functions
def determine_kmer(fastq, original_kmers):
    """Determines the 'optimal' K-mer value for a SPAdes assembly

    Key arguments:
        fastq -- str, pathway to a fastq file that contains one or multiple
                 reads. Can be a .gz file.
        original_kmers -- str, of odd numbers separated by ',' with a maximum
                          of 128.

    Returns:
        kmer_str -- str, a string of integers separated by ',' based on the
                    original integers given as input. Highest K-mer in the str
                    is the closest to 75% of the read length.
                    E.g. original_kmers = "33, 55, 77, 99, 127" and the avg
                    read length = 99.9. kmer_str = "33, 55, 77", as 77 is the
                    highest number closest to 75% of the avg read length.

    TODO:
        * Check that the original_kmers < 128 and odd numbers.
    """
    # Determine the avg read length of the fastq file
    cmd_seqkit = "seqkit stats -a -T {}".format(fastq)

    file_stats = subprocess.check_output(cmd_seqkit, stderr=subprocess.STDOUT,
                                         shell=True).decode("UTF-8")\
                                                    .strip()\
                                                    .split("\n")
    file_stats_cleaned = []
    for line in file_stats:
        file_stats_cleaned.append(line.split("\t"))

    avg_leng = file_stats_cleaned[1][6]

    # Check if the given K-mers are odd and below 128
    original_kmers = [int(x) for x in original_kmers.split(",")]
    if any(kmer > 128 for kmer in original_kmers):
        raise ValueError("One of the given K-mers is greater than 128. SPAdes "
                         "can only take K-mers below 128 as input. Please "
                         "update the provided string of K-mers!")
    elif any(kmer % 2 == 0 for kmer in original_kmers):
        raise ValueError("The string of given K-mers contains a number that "
                         "is even. SPAdes can only take K-mers that are odd "
                         "as input. Please update the provided string of "
                         "K-mers!")

    # Determine the highest kmer that is 75% of read length
    optimal = .75 * float(avg_leng)
    highest_kmer = min(original_kmers, key=lambda x: abs(x - optimal))

    # Create a list of kmers with the optimal kmer as highest one
    kmer_str = ""
    for kmer in original_kmers:
        if kmer <= highest_kmer:
            kmer_str += "{},".format(kmer)
    # Delete the last comma from the string
    kmer_str = kmer_str[:-1]

    return kmer_str


if __name__ == "__main__":
    """This is a helper module for meta_accession_parser.py"""
    pass
