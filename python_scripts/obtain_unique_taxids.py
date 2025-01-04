#!/usr/bin/env python3
"""Author: Devin van Valkengoed

Date: 28-Apr-2023
Description: this script can be used to determine all the unique taxonomic
             identifiers (TaxIDs) within a given list of accessions. The
             metadata of these accessions should be stored locally.

Usage: python3 <obtain_unique_taxids.py> -m <accessions.txt> -u <undefined.txt>

Required arguments:
    accessions.txt -- str, pathway to the file containing the pathway to the
                      metadata of all accessions on every new line of the
                      file.
    undefined.txt -- str, pathway and name of the file where all the
                     filenames should be stored of which the TaxID could not
                     be determined based on the given metadata.

TODO:
"""
# Import statements
import argparse
import subprocess
import sys


# Create a class to be able to print coloured messages to the stdout
class Bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


# Commandline parsing
def parsing_cmd():
    """Uses argparse to determine the given commandline options """
    parser = argparse.ArgumentParser(description="This script can be used to "
                                                 "determine all the unique "
                                                 "taxonomic identifiers ("
                                                 "TaxIDs) within a given list "
                                                 "of accessions. The metadata "
                                                 "of these accessions should "
                                                 "be stored locally. The "
                                                 "unique TaxIDs are printed "
                                                 "to stdout!")
    parser.add_argument("-m", "--metadata", help="A file containing the "
                                                 "relative path to the "
                                                 "metadata of all the "
                                                 "provided accessions on "
                                                 "every new line.",
                        metavar="<~/metadata_pathways.txt>", dest="meta",
                        required=True)
    parser.add_argument("-u", "--undefined", help="Relative path to where "
                                                  "the accessions should be "
                                                  "written to for which no "
                                                  "TaxID could be "
                                                  "determined.",
                        metavar="<~/undefined.txt>", dest="undefined",
                        required=True)

    return parser


# Functions
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


def obtain_taxid(metadata_file, undefined_out):
    """Determines the NCBI Taxonomic identifier of an accession using metadata

    Key arguments:
        metadata_file -- str, pathway to a NCBI-SRA metadatafile format.
        undefined_out -- str, pathway to a file where the file names of the
                         metadata files should be stored of which the TaxID
                         could not be determined.

    Returns:
        taxid -- int, the TaxID belonging to the given metadata file.
    """
    cmd_obtain_taxid = "cut -d',' -f28 {} | tail -1".format(metadata_file)
    taxid = subprocess.check_output(cmd_obtain_taxid, shell=True).decode(
        "UTF-8").strip()

    if taxid.isnumeric():
        taxid = int(taxid)
    else:
        # print("TaxID of: '{}' could not be determined. Obtained: '{}' which "
        #       "is not possible as TaxID.".format(metadata_file, taxid))
        with open(undefined_out, "a+") as outfile:
            outfile.write("{}\t{}\n".format(metadata_file, taxid))
        raise ValueError("TaxID could not be determined!")

    return taxid


def main():
    """Main function that wraps all the code in this script """
    # Step 0.0: parse the cmd line arguments
    parser = parsing_cmd()
    args = parser.parse_args()

    # Step 2: read all the given metadata locations into memory
    meta_locations = line_parser(args.meta)

    # Step 3: determine the TaxIDs for the given accessions
    taxids = []
    for metadata in meta_locations:
        try:
            taxid = obtain_taxid(metadata, args.undefined)
            taxids.append(taxid)
        except ValueError:
            continue

    # Step 4: print all unique TaxIDs to stdout
    for identifier in set(taxids):
        print(identifier)


if __name__ == "__main__":
    main()
