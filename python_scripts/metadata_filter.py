#!/usr/bin/env python3
""" Author: Devin van Valkengoed

Description: python script that can filter accessions based on library layout
             being paired- or single-end. The two different accessions are
             written to separate files.

Usage: python3 -a <accessions> -s <single_file> -p <paired_file> [-q]

Key arguments:
    accessions -- str, pathway to a file containing NCBI-SRA accessions on
                  every new line.
    single_file -- str, pathway to a file where the single-end accessions
                   should be stored.
    paired_file -- str, pathway to a file where the paired-end accessions
                   should be stored.

Optional argument:
    -q -- When this option is give, the script will not print any messages to
          the stdout.
"""
# Import statements
import sys
import subprocess
import argparse
import os.path
import warnings
from multiprocessing import Pool
from itertools import repeat


# Commandline parsing
def parsing_cmd():
    """Uses argparse to obtain the commandline arguments """
    parser = argparse.ArgumentParser(description="This script can filter NCBI-"
                                                 "SRA accessions based on the "
                                                 "sequencing library lay-out. "
                                                 "Paired- and single-end "
                                                 "accessions are written to "
                                                 "different files. ")
    parser.add_argument("-a", help="Pathway to file containing NCBI-SRA "
                                   "accession numbers on every new line. ",
                        dest="accessions", required=True, metavar="<str>")
    parser.add_argument("-s", help="Absolute pathway to the file where the "
                                   "single-end accessions should be stored.",
                        dest="single", required=True, metavar="<str>")
    parser.add_argument("-p", help="Absolute pathway to the file where the "
                                   "paired-end accessions should be stored.",
                        dest="paired", required=True, metavar="<str>")
    parser.add_argument("-u", help="Absolute pathway to a file where the "
                                   "accession should be stored in, for which "
                                   "the library layout can not be determined.",
                        dest="undetermined", required=True, metavar="<str>")
    parser.add_argument("-q", help="Quiet option, no output will be printed "
                                   "to the stdout",
                        dest="quiet", action="store_true")

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
                # Remove possible comments after line
                if "#" in line:
                    line = line.split("#")[0]
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


def check_if_paired(accession, undetermined_accessions):
    """Checks the library layout of an SRA accession based on number of files

    Key arguments:
        accessions -- str, NCBI SRA accession.
        undetermined_accessions -- str, pathway to a file where accession
                                   should be written to, of which the library
                                   layout can not be determined.

    Returns:
        PAIRED -- str, if the data contains two files and thus seems to be
                  paired-end.
        SINGLE -- str, if the data only contains one file and thus seems to
                  be single-end.

    Raises:
        ValueError if there are no files or, if there are more than two files.
        ValueError if the command is not found and 'ffq' is likely not
        installed.
    """
    cmd_liblayout = "ffq --ftp {}".format(accession)

    try:
        metadata = subprocess.check_output(
            cmd_liblayout, shell=True,
            stderr=subprocess.STDOUT).decode("UTF-8")

        if len(metadata.split("}")) == 3:
            return accession, "PAIRED"
        elif len(metadata.split("}")) == 2:
            return accession, "SINGLE"
        else:
            print("Something is wrong here, ffq call returned an unknown "
                  "number of files for accession: '{}'.".format(accession))
            write_line(accession, undetermined_accessions)
    except subprocess.CalledProcessError as error:
        # Check if conda env is active/FFQ is installed
        if error.returncode == 127:
            raise ValueError("Command not found, most likely because the "
                             "conda environment is not active.") from None
        else:
            return error, accession


def check_if_parsed(line, file_names):
    """Determines if the current line already parsed and present in a file

    Key arguments:
        line -- str, accession number of the current provided metadata.
        file_names -- lst, list of strings containing pathways to the files
                      that should be checked for the presence of the provided
                      line.

    Returns:
        parsed -- bool, True if the line is found in one of the files, False if
                  not.
    """

    if type(file_names) == list:
        files = " ".join(file_names)
    else:
        files = file_names
    cmd = "grep '{}' {}".format(line, files)

    try:
        subprocess.check_output(cmd, shell=True)
        parsed = True
    except subprocess.CalledProcessError:
        parsed = False

    return line, parsed


def check_existence(file_list, quiet):
    """ Function that checks if files already exist or not and creates them

    Key Arguments:
        file_list -- lst, containing strings specifying file pathways to check.

    Returns:
        None
    """
    for file in file_list:
        if not os.path.isfile(file):
            cmd = "touch {}".format(file)
            subprocess.check_output(cmd, shell=True)
            if not quiet:
                print("{} does not yet exist, creating it now.".format(file))

    return


def check_liblayout(input_accessions, single_end, paired_end, undetermined,
                    quiet):
    """Writes accessions to the single- or paired-end file

    Key arguments:
        accession_file -- str, pathway to a file containing NCBI-SRA
                          accessions on every new line.
        single_end -- str, pathway to a file to store the single-end
                      accessions in.
        paired_end -- str, pathway to a file to store the paired-end
                      accessions in.
        undetermined -- str, pathway to a file where the accessions should be
                        stored of which the library layout could not be
                        determined.
        quiet -- bool, True if no messages should be printed to the stdout.

    Returns:
        None
    """
    # Create the output files
    check_existence([single_end, paired_end], quiet)

    # Parse the accessions to the right file based on library layout
    not_parsed = []
    for accession in input_accessions:
        # Check if the accession is in either of the files already
        # if not: write to the correct file
        if not check_if_parsed(accession, [single_end, paired_end])[1]:
            not_parsed.append(accession)
        else:
            if not quiet:
                print("Accession: '{}' was already parsed! Skipping"
                      .format(accession))
            # Accession is already parsed
            continue

    # Determine the library layout of the accession that are not parsed yet
    # NOTE! do not increase the bool above '3' as this will give 429 errors
    with Pool(3) as pool:
        layout = pool.starmap(check_if_paired,
                              zip(not_parsed, repeat(undetermined)))
    pool.close()

    for combo in layout:
        max_retries = 0
        while combo[1] not in ["PAIRED", "SINGLE"]:
            print(f"combo {combo} not paired or single")

            if max_retries < 10:
                # Check if the pool did not have any 429 errors
                if type(combo[0]) == subprocess.CalledProcessError:
                    combo = check_if_paired(combo[1], undetermined)
                max_retries += 1
                print(f"Retries used of max 10: '{max_retries}'")
            else:
                warnings.warn("Maximum of library retries exceeded!\n"
                              "Writing accession: '{}' into the "
                              "undermentioned file.".format(combo[0]))
                write_line(combo[0], undetermined)
                break

        # Write the accessions to the correct file
        if combo[1] == "PAIRED":
            write_line(combo[0], paired_end)
            if not quiet:
                print("Writing accession: '{}' into the paired-end "
                      "file.".format(combo[0]))
        elif combo[1] == "SINGLE":
            write_line(combo[0], single_end)
            if not quiet:
                print("Writing accession: '{}' into the single-end "
                      "file.".format(combo[0]))

    return


def main():
    """Main function that wraps all the code in this script """
    # Step 0.0: check the given commandline arguments
    parser = parsing_cmd()
    args = parser.parse_args()

    # Step 1: parse the accession in the input file into a list
    accessions = line_parser(args.accessions)

    # Step 2: check the library layout and write to the correct file
    check_liblayout(accessions, args.single, args.paired, args.undetermined,
                    args.quiet)


if __name__ == "__main__":
    main()
