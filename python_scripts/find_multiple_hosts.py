#!/usr/bin/env python3
"""Author: Devin van Valkengoed

Date: 06-Dec-2022 V1
      27-Jan-2023 V2.0
      06-Apr-2023 V2.1 (not excluding species defined at 'subspecies' level
                        anymore)
Description: Script used to determine if given NCBI TaxIDs are
             likely to contain multiple organisms, genus instead of species and
             subspecies instead of species.

Usage: python3 <find_multiple_hosts.py> -t <input_taxids> [--print]

Key argument:
    input_taxids -- str, pathway to a file that contains NCBI TaxID's on every
                    new line.

Optional argument:
    --print -- str, if the commandline option is given, the output will be
               printed to the stdout.
"""
# Import statements
import sys
import re
import subprocess
import argparse
from multiprocessing import Pool


# Commandline parsing
def parsing_cmds():
    """Uses argparse to obtain the commandline arguments """
    parser = argparse.ArgumentParser(description="This script is used to "
                                                 "determine if the hosts, "
                                                 "given by a file that "
                                                 "contains NCBI TaxIDs, are "
                                                 "of valid format. I.e. "
                                                 "single species at genus "
                                                 "level.")
    parser.add_argument("--print", help="If the print flag is given, "
                                        "output will be printed to "
                                        "the screen.",
                        action="store_true", dest="print")
    parser.add_argument("-t", "--taxids", help="str, pathway to a file that "
                                               "contains NCBI Taxonomic "
                                               "Identifiers on every new line",
                        metavar="<TaxID_file>", type=str, dest="taxid_file",
                        required=True)
    parser.add_argument("-o", "--output", help="str, pathway to a file where "
                                               "the abnormal TaxIDs should be "
                                               "writen to.",
                        metavar="<Output_file", type=str, dest="out_file")

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


def taxid_metadata(taxid):
    """Uses NCBI prefetch to obtain the species information based on TaxID

    Key arguments:
        taxid_list -- lst, containing one or multiple TaxID's as elements.

    Returns:
        taxid2name -- lst, with the organism information and the original
                      TaxID. E.g. ['1. Bombyx mori (domestic silkworm),
                      species, moths', '7091']
    """
    # Obtain the species information using efetch
    cmd = 'efetch -db taxonomy -id {}'.format(taxid)
    species = subprocess.check_output(cmd, shell=True).decode("UTF-8")
    species = species.split("\n")

    clean_str = ""
    for element in species:
        if element:
            if clean_str:
                clean_str += str(" " + element.strip())
            else:
                clean_str += str("" + element.strip())
        else:
            continue

    taxid2name = [clean_str, taxid]

    return taxid2name


def fetch_taxid2name(taxid_list):
    """Function that uses Pool to fetch the organism info of multiple TaxID's

    Key arguments:
        taxid_list -- lst, containing NCBI TaxID's as elements.

    Returns:
         abnormal_ids -- lst, list of Taxonomic Identifiers that either
                         contain multiple species, are sub-species or are not
                         IDs at species level.
    """
    with Pool(2) as pool:
        results = pool.map(taxid_metadata, taxid_list)

    pool.close()

    abnormal_ids = []
    for organism in results:
        if " x " in organism[0]:
            abnormal_ids.append(organism[1])
        # elif "sub" in organism[0]:
        #     abnormal_ids.append(organism[1])
        elif "species" not in organism[0]:
            abnormal_ids.append(organism[1])

    return abnormal_ids


def write_output(abnormal_taxids, output_file, original_input):
    """Writes the output containing abnormal TaxIDs to the given output file

    Key arguments:
        abnormal_taxids -- lst, list of Taxonomic Identifiers that either
                           contain multiple species, are sub-species or are not
                           IDs at species level.
        output_file -- str, pathway to the file where the abnormal TaxIDs
                       should be writen to.

    Returns:
          None
    """
    with open(output_file, "w") as out_file:
        out_file.write("# Original file containing NCBI Taxonomic "
                       "Identifiers: '{}' \n".format(original_input))
        out_file.write("# TaxIDs that are not usable (containing multiple "
                       "species, or not defined at (sub-)species level): \n")
        for accession in abnormal_taxids:
            out_file.write(accession + "\n")

    return


def main():
    """Main function that wraps all the code in this script """
    # Step 0.0: parse the commandline arguments
    parser = parsing_cmds()
    args = parser.parse_args()

    # Step 1: parse all the taxids from every new line of the given file
    taxid_lst = line_parser(args.taxid_file)

    # Step 2: obtain the organism information for every TaxID
    abnormal_taxids = fetch_taxid2name(taxid_lst)

    # Step 3: either print or return the abnormal taxids, based on cmd line
    if args.print:
        print("Taxonomic Identifiers that contain multiple species or are not "
              "defined at species level: ")
        for taxid in abnormal_taxids:
            print(taxid)
    else:
        write_output(abnormal_taxids, args.out_file, args.taxid_file)


if __name__ == "__main__":
    main()
