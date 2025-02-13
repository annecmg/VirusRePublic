#!/usr/bin/env python3
"""Author: Devin van Valkengoed

Date: 25-Apr-2023
      24-Jan-2025: v2.0
      10-Feb-2025: v2.1 (minor changes/improvements mainly on the file paths.)
Description: script that can be used to determine the mapping statistics used
             when raw sequencing reads are mapped against a reference genome.

Usage: python3 mapping_stats.py [-h] -a ACCESSIONS -b BASE -u UNAV -o OUTPUT

Required arguments:
    -a/--accessions -- str, pathway to a file containing NCBI SRA accessions
                       on every new line.
    -u/--unavailable -- str, pathway to a file where the accessions are
                        writen to for which no mapping stats could be
                        determined.
    -o/--output -- str, pathway (and name) of a file to which the output
                   statistics are writen.

Optional arguments:
    -b/--base -- str, pathway to the base location where the log files are
                 stored.
    -h/--help -- str, show help message and exit.

"""
# Import statements
import subprocess
import argparse
from pathlib import Path

# Set the default relative log file directory based on the Snakemake workflow output
SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_LOG_DIR = SCRIPT_DIR / "../output/logs/main_logs"


# Commandline parsing
def parsing_cmd_line():
    """Uses argparse to obtain the commandline arguments """
    parser = argparse.ArgumentParser(description="This script determines the "
                                                 "mapping statistics per accession based "
                                                 "on the log files that are "
                                                 "produced when running the "
                                                 "main 'VirusRePublic' "
                                                 "Snakemake pipeline. ")
    parser.add_argument("-a", "--accessions", type=str, required=True,
                        help="Pathway to a file containing NCBI-SRA "
                             "accessions on every new line, of which the "
                             "mapping percentage should be determined.",
                        dest="accessions")
    parser.add_argument("-u", "--unavailable", type=str, required=True,
                        help="Pathway to a file where the accessions of "
                             "which no mapping statistic could be determined "
                             "are writen to.",
                        dest="unav")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Pathway to the output file where the "
                             "statistics should be written to.",
                        dest="output")
    parser.add_argument("-b", "--base", type=str, default=DEFAULT_LOG_DIR,
                        help="Base directory to the log files that are created by the Snakemake workflow (named "
                             "main_<accession_id>.log)",
                        dest="base")

    return parser


# Functions
def determine_mapping_percentage(logfile):
    """Obtains the percentage of reads that map against the host genome

    Key arguments:
        logfile -- str, pathway to the logfile that contains the percentage
                   of mapped reads.

    Returns:
        mapping_perc -- float, percentage of mapped reads.
    """
    cmd_mapping = "grep -F 'overall alignment rate' {} 2>/dev/null".format(logfile)

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


def write_output(mapping_stats, output_file):
    """Writes the mapping statistics to the output file

    Key arguments:
        mapping_stats -- list of strings, each containing accession and mapping percentage.
        output_file -- str, pathway to the file where the statistics should be written.

    Returns:
        None
    """
    with open(output_file, 'w') as output:
        output.write("accession\tmapping_percentage\n")
        for stat in mapping_stats:
            output.write(stat + "\n")


def main():
    """Main function that wraps all the code of this script """
    # Step 0.0: parse the commandline arguments
    parser = parsing_cmd_line()
    args = parser.parse_args()
    args.base = Path(args.base)
    print("Determining the mapping statistics for the given accessions in "
          "file: '{}'. Please wait!".format(args.accessions))

    # Step 1: read all the given accessions
    accessions = line_parser(args.accessions)

    # Step 2: check if the log files are present
    log_files = list(args.base.glob("*main*.log"))

    if not log_files:
        raise FileNotFoundError("No log files with the name: 'main_<accession_id>.log' were found in the given "
                                "log file directory! Please check if the files are there or give another directory "
                                "using argument: '-b /my_logfile/dir/'")

    # Step 3: loop through the individual accessions, determine mapping stats
    mapping_stats = []

    for acc in accessions:
        logfile_location = args.base / "main_{}.log".format(acc)

        # Try to determine the mapping percentage, if not write to unavailable
        try:
            mapping_perc = determine_mapping_percentage(logfile_location)
            mapping_stats.append(f"{acc}\t{mapping_perc}")
        except subprocess.CalledProcessError:
            with open(args.unav, 'w+') as unavailable:
                unavailable.write(str(acc) + "\n")

        # Step 3: Write the mapping statistics to the output file
        write_output(mapping_stats, args.output)

    print("Mapping statistics have been written to '{}'.".format(args.output))


if __name__ == '__main__':
    main()
