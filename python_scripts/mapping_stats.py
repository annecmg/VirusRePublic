#!/usr/bin/env python3
"""Author: Devin van Valkengoed

Date: 25-Apr-2023
Description: script that can be used to determine the mapping statistics used
             when raw sequencing reads are mapped against a reference genome.

Usage: python3 mapping_stats.py [-h] -a ACCESSIONS -b BASE -u UNAV -o OUTPUT

Required arguments:
    -a/--accessions -- str, pathway to a file containing NCBI SRA accessions
                       on every new line.
    -b/--base -- str, pathway to the base location where the log files are
                 stored.
    -m/--metadata -- str, pathway to the base folder where all metadata
                     files of the provided accessions are stored.
    -u/--unavailable -- str, pathway to a file where the accessions are
                        writen to for which no mapping stats could be
                        determined.
    -o/--output -- str, pathway (and name) of a file to which the output
                   statistics are writen.

Optional arguments:
    -h/--help -- str, show help message and exit.

TODO:
    * Create write output function per TaxID
    *
"""
# Import statements
import subprocess
import argparse
import sys


# Commandline parsing
def parsing_cmd_line():
    """Uses argparse to obtain the commandline arguments """
    parser = argparse.ArgumentParser(description="This script determines the "
                                                 "average mapping statistics "
                                                 "per taxonomic host based "
                                                 "on the log files that are "
                                                 "produced when running the "
                                                 "main 'iflavirus paper' "
                                                 "Snakemake pipeline. ")
    parser.add_argument("-a", "--accessions", type=str, required=True,
                        help="Pathway to a file containing NCBI-SRA "
                             "accessions on every new line, of which the "
                             "mapping percentage should be determined.",
                        dest="accessions")
    parser.add_argument("-b", "--base", type=str, required=True,
                        help="Pathway to a folder containing all the "
                             "logfiles where the mapping percentage is "
                             "stored in.",
                        dest="base")
    parser.add_argument("-m", "--metadata", type=str, required=True,
                        help="str, pathway to the base folder where all "
                             "metadata files of the given accessions are "
                             "stored. ",
                        dest="metadata")
    parser.add_argument("-u", "--unavailable", type=str, required=True,
                        help="Pathway to a file where the accessions of "
                             "which no mapping statistic could be determined "
                             "are writen to.",
                        dest="unav")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Pathway to the output file where the "
                             "statistics should be written to.",
                        dest="output")

    return parser


# Functions
def determine_taxid(metadata):
    """Uses the (local) metadata of an accession to determine the host TaxID

    Key arguments:
        metadata -- str, pathway to the metadata file of the accession of
                    which the TaxID should be determined.

    Returns:
        taxid -- int, taxonomic identifier as determined by NCBI-SRA obtained
                 from the provided metadata.
    """
    cmd_taxid = "cut -d',' -f28 {} | tail +2".format(metadata)

    try:
        taxid = int(subprocess.check_output(cmd_taxid, shell=True).decode(
            "UTF-8").strip())
    except ValueError:
        raise ValueError("TaxID of one of the accessions could not be "
                         "determined, exiting the script!") from None

    return taxid


def determine_mapping_percentage(logfile):
    """Obtains the percentage of reads that map against the host genome

    Key arguments:
        logfile -- str, pathway to the logfile that contains the percentage
                   of mapped reads.

    Returns:
        mapping_perc -- float, percentage of mapped reads.
    """
    cmd_mapping = "grep -F 'overall alignment rate' {}".format(logfile)

    map_perc = subprocess.check_output(cmd_mapping, shell=True).decode(
        "UTF-8").strip()

    # Remove the trailing '%' sign and convert to type float.
    mapping_perc = float(map_perc.split(" ")[0].rstrip("%"))

    return mapping_perc


def line_parser(input_file):
    """Parses every new line of a file as element in a list

    Key arguments:
        input_file -- str, pathway to a .txt file containing information on
                      every new line.

    Returns:
        line_list -- lst, list containing all the lines of the file as elements
    """
    line_list = []
    with open(input_file, 'r') as in_file:
        for line in in_file:
            if not line:
                continue
            line = line.strip()
            if line.startswith("#"):
                continue
            else:
                line_list.append(line)

    return line_list


def write_line(line, output_file):
    """Function that writes the given string to new line in the give file

    Key Arguments:
        line -- str, line that has to be appended to the specified output file.
        output_file -- str, pathway to the file to which the line should be
                       appended.

    Returns:
        None
    """
    with open(output_file, 'a+') as output:
        output.write(str(line) + "\n")

    return


def write_output():
    """Writes the mapping statistics to the output file

     Key arguments:


     Returns:
         None
    """



def main():
    """Main function that wraps all the code of this script """
    # Step 0.0: parse the commandline arguments
    parser = parsing_cmd_line()
    args = parser.parse_args()
    print("Determining the mapping statistics for the given accessions in "
          "file: '{}'. Please wait!".format(args.accessions))

    # Step 1: read all the given accessions
    accessions = line_parser(args.accessions)

    # Step 2: loop through the individual accessions, determine mapping stats
    taxid_stats = {'total': [0, 0], 'unavailable': [0, 0]}
    indivi_stats = {}

    for acc in accessions:
        metadata_loc = "{}{}_metadata.txt".format(args.metadata, acc)
        logfile_location = "{}main_{}.log".format(args.base, acc)

        # Try to determine the TaxID, if not possible, write to unavailable
        try:
            taxid = determine_taxid(metadata_loc)
        except ValueError:
            taxid_stats['unavailable'][0] += 1
            write_line(acc, args.unav)
            continue

        # Try to determine the mapping percentage, if not write to unavailable
        try:
            mapping_perc = determine_mapping_percentage(logfile_location)
        except subprocess.CalledProcessError:
            taxid_stats['unavailable'][0] += 1
            write_line(acc, args.unav)
            continue

        # Create a total count and percentage (per TaxID)
        if taxid not in taxid_stats.keys():
            taxid_stats[taxid] = [1, mapping_perc]
            taxid_stats['total'][0] += 1
            taxid_stats['total'][1] += mapping_perc
        else:
            taxid_stats[taxid][0] += 1
            taxid_stats[taxid][1] += mapping_perc
            taxid_stats['total'][0] += 1
            taxid_stats['total'][1] += mapping_perc

    # Calculate the average mapping percentage per TaxID
    for identifier in taxid_stats:
        taxid_stats[identifier][1] = round(taxid_stats[identifier][1] /
                                           taxid_stats[identifier][0], 2)


    for tax_identifier in taxid_stats:
        print("{}\t{}\t{}".format(tax_identifier,
                                  str(taxid_stats[tax_identifier][0]),
                                  taxid_stats[tax_identifier][1]))





    # Write the output




if __name__ == '__main__':
    main()
