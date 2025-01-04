#!/usr/bin/env python3
"""Author: Devin van Valkengoed

Date: 02-Dec-2022
      19-Apr-2023 V2.0 Update, renaming *.sralite downloaded files to *.sra
Description: python script that can be used in Snakefile_downloading.smk to
    download sra datasets using prefetch. The script will limit the download
    time to 20 min. and restart the download if it takes longer. This is done,
    because the NCBI servers some time get stuck on the download of a single
    accession. When the download is then terminated and re-started, it often
    finishes way faster compared waiting for the initial download to finish!

Usage: python3 <sra_download.py> -a <accessions>

Required arguments:
    sra_download.py -- str, absolute pathway to this python script
    accession_list -- str, one or multiple NCBI-SRA accessions in string format
                      separated by spaces. (E.g. SRR1234 SRR5678 SRR91011)

Optional arguments:
    none
"""
# Import statements
import sys
import argparse
import subprocess


# Create a class to be able to print coloured messages
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
    """ Uses argparse to obtain the commandline arguments """
    parser = argparse.ArgumentParser(description="This script takes a list of "
                                                 "NCBI-SRA accession numbers "
                                                 "and downloads the "
                                                 "accessions. It can be limit "
                                                 "the max download time per "
                                                 "accession, by specifying a "
                                                 "time in seconds. If 'inf' "
                                                 "is given for time, the "
                                                 "script will not limit the "
                                                 "download time. "
                                                 "Default = 1200 (20 min.)",
                                     usage="sra_download.py [-h] -a "
                                           "<SRRXXX [...]> "
                                           "[-t time]")
    parser.add_argument("-a", help="NCBI-SRA accession numbers "
                                   "separated by spaces. E.g. SRR123 SRR456 "
                                   "SRR789",
                        dest="accessions", required=True, metavar="str",
                        nargs="*")
    parser.add_argument("-t", help="Number of seconds that the "
                                   "download per accession is limited to. "
                                   "When 'inf' is given, the download time is "
                                   "not limited. Default = 1200 (20 min.)",
                        dest="time", metavar="int")
    parser.add_argument("-o", help="Output directory for the download.",
                        dest="out", required=True, metavar="str")
    parser.add_argument("-m", "--max_size", help="The maximum file size in Kb "
                                                 "that prefetch is allowed to "
                                                 "download at once. Default = "
                                                 "100000000 Kb (100 Gb).",
                        dest="max_size", metavar="int", type=int)

    return parser


# Functions
def download_sra_data(accessions, time, output_dir, max_size=100000000):
    """ Function that uses NCBI prefetch shell command to download SRA data

    Key arguments:
        accessions -- lst, list of NCBI accessions that have to be downloaded
        time -- int, number of seconds the function waits for the download of
                each accession to finish. If 'inf' is given as time, the
                download will not be restricted to a certain time.
                Default = 1200 (20 minutes).
        output_dir -- str, pathway to the director where the downloaded file
                      should be stored.
        max_size -- int, maximum file size in Kb that prefetch is allowed to
                    download at once. Default = 100000000 Kb (100 Gb).

    Returns:
        remaining -- lst, list of NCBI accessions of which the download did not
                     finish in the set time.
    """
    # Set the time limit to 1200 seconds if no custom limit is given
    if not time:
        time = 1200

    # Download the accessions using NCBI sra-tools 'prefetch' in shell
    remaining = []
    for acc in accessions:
        if time == "inf":
            sys.stderr.write(Bcolors.OKCYAN +
                             "No time limit per accession \n" +
                             Bcolors.ENDC)
            cmd = "prefetch -f all -X {} {} " \
                  "--output-directory {}".format(max_size, acc, output_dir)
            subprocess.check_output(cmd, shell=True)
        else:
            try:
                sys.stderr.write(Bcolors.OKCYAN +
                                 "Given time limit per accession: {} "
                                 "seconds \n".format(time) +
                                 Bcolors.ENDC)
                cmd = ("perl -e 'alarm shift; exec @ARGV' {} prefetch -f all "
                       "-X {} {} --output-directory {}").format(time,
                                                                max_size, acc,
                                                                output_dir)
                subprocess.check_output(cmd, shell=True)

                # Check the extension of the file that was downloaded
                cmd_filetype = "find {} -name '*.sralite' -type f " \
                               "| wc -l".format(output_dir + acc + "/")
                sra_lite = subprocess.check_output(cmd_filetype, shell=True) \
                    .decode("UTF-8").strip()

                if int(sra_lite) > 0:
                    cmd_rename = "mv {}{}/{}.sralite " \
                                 "{}{}/{}.sra".format(output_dir, acc, acc,
                                                      output_dir, acc, acc)

                    print("Renaming the '{}.sralite' file into "
                          "{}.sra using the following "
                          "command: '{}'".format(acc, acc, cmd_rename))
                    subprocess.check_output(cmd_rename, shell=True)
                else:
                    continue
            # Restart when the set time has passed.
            except subprocess.CalledProcessError as errcode:
                if errcode.returncode == 124:
                    remaining.append(acc)
                    sys.stderr.write(Bcolors.WARNING +
                                     "Reached the time limit for this "
                                     "accession! Retry at the end of "
                                     "the cycle. \n" +
                                     Bcolors.ENDC)
                else:
                    sys.stderr.write(errcode.output.decode("UTF-8") + "\n")

    return remaining


def main():
    """ Function that wraps all the code in this script"""
    # Step 0.0: check if the command line input is given
    parser = parsing_cmd()
    args = parser.parse_args()

    # Step 1: Run a loop to download all the provided accessions
    leftover = download_sra_data(args.accessions, args.time, args.out)
    while leftover:
        sys.stderr.write(Bcolors.WARNING +
                         "The download of some accessions was not "
                         "yet finished! \n"
                         "Rerunning the script for the "
                         "following accessions: {} \n".format(leftover) +
                         Bcolors.ENDC)
        leftover = download_sra_data(args.accessions, args.time,
                                     args.out, args.max_size)
    else:
        sys.stderr.write(Bcolors.OKGREEN + "Finished the download of all "
                                           "given accessions \n" +
                         Bcolors.ENDC)


if __name__ == '__main__':
    main()
