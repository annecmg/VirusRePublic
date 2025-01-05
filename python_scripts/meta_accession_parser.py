#!/usr/bin/env python3
"""Author: Devin van Valkengoed

Date: 01-Nov-2022
Description: script that takes a NCBI-SRA metadata format .txt file as input
    and parses this file to obtain an accession number and the LibraryLayout
    (either paired or single based on read type). It then writes the accession
    number to the corresponding output .txt file to sort the experiments based
    on single or paired end data.

Usage: python3 <accession_parser.py> <metadata_file> [options]

Required arguments:
    -m/--metadata -- str, pathway to a .txt file that contains metadata of a
                     sequencing experiment in NCBI-SRA format.

Optional arguments:
    -p/--paired -- str, pathway to a .txt file to store the accession numbers
                   of the samples that have a paired LibraryLayout.
    -s/--single -- str, pathway to a .txt file to store the accession numbers
                   of the samples that have a single LibraryLayout.

TODO:
    * Update this docstring with required and optional arguments.
    * Implement RefSeq sorting helper module!
    * Make sure the pipeline is fully automatic, also being able to use lists
      as input. Make the pipeline automatic until the end of BLASTx search.
        Interpretation of these results will be hard to automate.
    * Re-write function parse_layout to only output an accession and create new
      function to obtain the bool, based on paired- or single-end data.
    * Make the write function optional, so the accessions can also be printed?
    * Write functions to obtain the unique species in the given metadata_dir/
      and write or print them.
    * Automated download of the host genomes? In combination with snakemake.
    * Write different wrapper functions to combine functions for different
      options?
    * Extract the writing function and evoke it with a option in the cmd,
      somehow re-using the function for obtaining accessions and taxids into
      files.
    * Update main docstring based on new input.
    * taxid subcommand can now have both -c and -m flag, find a way to only
      allow -c flag.


OPTIONAL:
    * Give option to provide a list of accessions as input?
"""
# Import statements
import subprocess
import sys
import os.path
import argparse
import os
from datetime import datetime
from pathlib import Path
from check_genome_assemblies import *
from fetch_assembly_path import *
from determine_kmer import *


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
def parsing_cmd_line():
    """ Uses argparse to obtain the commandline arguments """
    parser = argparse.ArgumentParser(description="This script takes a "
                                                 "NCBI-SRA metadata file as "
                                                 "input and can extract "
                                                 "multiple parameters, such "
                                                 "as accession number, "
                                                 "LibraryLayout, etc.")
    parser.add_argument("--print", help="If the print flag is given, "
                                        "the normal behaviour of every "
                                        "subcommand will be overruled "
                                        "and output will be printed to "
                                        "the screen.",
                        action="store_true", dest="print")
    subparsers = parser.add_subparsers(title="Subcommands",
                                       description="Valid subcommands",
                                       help="Please choose one of these sub-"
                                            "commands. To show the available "
                                            "options specific for each "
                                            "subcommand please type: '{} "
                                            "[subcommand] -h'"
                                       .format(os.path.basename(__file__)),
                                       dest='command')
    # Subcommands
    # Library layout parsing
    liblayout_parser = subparsers.add_parser("library",
                                             help="This option parses the "
                                                  "library layout of the "
                                                  "metadata into two files "
                                                  "(Single- and paired-end).")
    liblayout_parser.add_argument("-m", "--metadata", type=str, required=True,
                                  help="str -- pathway to a .txt file that "
                                       "contains metadata of a sequencing "
                                       "experiment in NCBI-SRA format.",
                                  metavar="<metadata>", dest="metadata")
    liblayout_parser.add_argument("-p", "--paired", type=str,
                                  help="str -- pathway to a .txt file to "
                                       "store the accession numbers of the "
                                       "samples that have a paired "
                                       "LibraryLayout.",
                                  metavar="[paired_outfile]", dest="paired")
    liblayout_parser.add_argument("-s", "--single", type=str,
                                  help="str -- pathway to a .txt file to "
                                       "store the accession numbers of the "
                                       "samples that have a single "
                                       "LibraryLayout.",
                                  metavar="[single_outfile]", dest="single")

    # Obtaining taxid
    taxid_parser = subparsers.add_parser("taxid", help="This option parses "
                                                       "the NCBI Taxonomic "
                                                       "identifier from the "
                                                       "provided metadata. "
                                                       "Default is a file "
                                                       "with unique taxid's.")
    taxid_parser.add_argument("-m", "--metadata", type=str,
                              help="str -- pathway to a .txt file that "
                                   "contains metadata of a sequencing "
                                   "experiment in NCBI-SRA format.",
                              metavar="[metadata]", dest="metadata",
                              required=True)
    taxid_parser.add_argument("-o", "--output", type=str,
                              help="Specify the path and name of the "
                                   "output file.", metavar="[output]",
                              dest="output")
    taxid_parser.add_argument("-u", "--undefined",
                              help="Taxonomic Identifiers that, for some "
                                    "reason could not be determined, are "
                                    "written to this file.",
                              dest="undefined", metavar="<undefined>",
                              type=str)
    taxid_parser.add_argument("-a", "--all",
                              help="If this option is given, the output file "
                                   "is not checked for the presence of the "
                                   "current accession. This will create "
                                   "duplicates taxid's in the output file.",
                              action="store_true")

    # Producing TaxID stats from multiple files
    taxid_count_parser = subparsers.add_parser("taxid_count",
                                               help="This option can be used "
                                                    "to obtain statistics "
                                                    "about the taxids "
                                                    "present in a file "
                                                    "containing NCBI-SRA "
                                                    "accessions. Note that "
                                                    "the metadata of these "
                                                    "accessions should be "
                                                    "stored locally. ")
    taxid_count_parser.add_argument("-a", "--accessions", type=str,
                                    required=True,
                                    help="str -- pathway to a .txt file that "
                                         "contains NCBI-SRA accessions on "
                                         "every new line.",
                                    metavar="<accession_file>",
                                    dest="accessions")
    taxid_count_parser.add_argument("-o", "--output", type=str,
                                    help="Specify the path and name of the "
                                         "output file. Only required when the "
                                         "global flag '--print' is not used.",
                                    metavar="[output]", dest="count_outfile")
    taxid_count_parser.add_argument("-b", "--metadata_base",
                                    help="str -- Base pathway to the "
                                         "directory where all the metadata "
                                         "files are stored.",
                                    metavar="[./metadata_directory/]",
                                    dest="metadata_base")
    taxid_count_parser.add_argument("-x", "--extension",
                                    help="str -- default extension that "
                                         "should be used for the all metadata "
                                         "files. Default = _metadata.txt",
                                    metavar="[metadata.txt]", dest="extension")

    # Obtaining a string of K-mers for SPAdes to use during genome assembly
    provide_kmers = subparsers.add_parser("kmers", help="This subcommand can "
                                                        "be used to print a "
                                                        "string of K-mers "
                                                        "that can be used for "
                                                        "a genome assembly "
                                                        "with SPAdes. The "
                                                        "highest K-mer in the "
                                                        "string will be 75%% "
                                                        "of the avg read "
                                                        "length.")
    provide_kmers.add_argument("-r", type=str, help="Specify the path to a "
                                                    "FastQ file containing "
                                                    "reads.",
                               required=True, dest="fastq")
    provide_kmers.add_argument("-k", type=str, help="A string with odd "
                                                    "integers, below 128 "
                                                    "which will be used as "
                                                    "desired K-mers for SPAdes"
                                                    ".",
                               required=True, dest="kmer_string")

    # Parsing genome assemblies
    genome_parser = subparsers.add_parser("genome", help="This option will "
                                                         "check if assembled "
                                                         "genomes are "
                                                         "available for a "
                                                         "given file with "
                                                         "TaxIDs on every new"
                                                         "line.",
                                          description="Change this based on "
                                                      "new options")
    genome_parser.add_argument("-t", "--taxids", type=str,
                               help="Specify the path to a file that contains "
                                    "NCBI TaxID's on every new line.",
                               metavar="<TaxID_file>", dest="taxids")
    genome_parser.add_argument("-q", "--quiet",
                               help="This option will supress the message "
                                    "that is printed to stdout if an "
                                    "accession is already in one of the "
                                    "files.",
                               action="store_true", dest="quiet")
    genome_parser.add_argument("-o1", "--ref_output", type=str,
                               help="Specify the path and name of the "
                                    "output file to write the TaxIDs with "
                                    "RefSeq assembly to. "
                                    "Only required when the "
                                    "global flag '--print' is not used.",
                               metavar="[output]", dest="ref_output")
    genome_parser.add_argument("-o2", "--no_ref_output", type=str,
                               help="Specify the path and name of the "
                                    "output file to write the TaxIDs with "
                                    "no RefSeq assembly (but a general "
                                    "assembly) to. "
                                    "Only required when the "
                                    "global flag '--print' is not used.",
                               metavar="[output]", dest="no_ref_output")
    genome_parser.add_argument("-o3", "--no_assembly_output", type=str,
                               help="Specify the path and name of the "
                                    "output file to write the TaxIDs with "
                                    "no assembly whatsoever to. "
                                    "Only required when the "
                                    "global flag '--print' is not used.",
                               metavar="[output]", dest="no_assem_output")

    return parser


# Functions
def read_metadata(meta_data):
    """Reads the provided metadata file into a nested list

    Key arguments:
        meta_data -- str, pathway to a .txt file providing the metadata of
                     SRA experiments in the standard NCBI format.
    Returns:
        file_lines -- lst, nested list containing two list. First list is the
                      headers of the NCBI-SRA metadata file and the second list
                      contains the values.
    """
    in_file = open(meta_data, 'r')
    file_lines = []

    # Save the two lines of a NCBI-SRA metadata file into a list
    for line in in_file:
        file_lines.append(line.strip().split(","))

    in_file.close()

    return file_lines


def check_format(meta_list, file_name):
    """Check whether the given metadata file is in the format of NCBI-SRA

    Key arguments:
        meta_list -- lst, nested list containing two list. First list is the
                     headers of the NCBI-SRA metadata file and the second list
                     contains the values.
        file_name -- str, name of the input metadata file.

    Raises:
        ValueError -- whenever the provided metadata file is not of NCBI-SRA
                      format.
    """
    # Perform multiple checks to determine whether the metadata is correct
    if len(meta_list) >= 2 and \
            meta_list[0][0] == "Run" and \
            meta_list[1][0] != '':
        return
    else:
        raise ValueError("The provided metadata in {} is not of the correct "
                         "format!".format(file_name))


def obtain_librarylayout(file_name, metadata_list):
    """ Finds the library layout of the metadata

    Key arguments:
        file_name -- str, pathway to the metadata file that has to be parsed.
        metadata_list -- lst, nested list containing two lists.
                         First list is the headers of the NCBI-SRA metadata
                         file and the second list contains the values.

    Returns:
        paired -- str, PAIRED or SINGLE based on librarylayout found in the
                  metadata list.
    """
    # Double check if LibraryLayout is at the position and not empty
    if metadata_list[0][15] == "LibraryLayout":
        if metadata_list[1][15] != '':
            if metadata_list[1][15] == "PAIRED":
                paired = "PAIRED"
            else:
                paired = "SINGLE"
        else:
            raise ValueError("Exception while parsing metadata of {}. "
                             "LibraryLayout position is empty and can thus "
                             "not be determined. Please double check if the "
                             "provided metadata file is in NCBI-SRA format!"
                             .format(file_name))
    else:
        raise ValueError("Exception while parsing metadata of {}. "
                         "LibraryLayout is not at the expected position and "
                         "cannot be determined.".format(file_name))

    return paired


def check_existence(file_list, datestamp=False):
    """ Function that checks if files already exist or not and creates them

    Key Arguments:
        file_list -- lst, containing strings specifying file pathways to check.

    Optional Arguments:
        datestamp -- bool, True if the current date should be appended to
                     the files. Default = False.

    Returns:
        None
    """
    for file in file_list:
        if not os.path.isfile(file):
            print("{} does not yet exist, creating it now.".format(file))
            Path(file).touch()

            if datestamp:
                # Obtain current date and time
                now = datetime.now()
                current_time = now.strftime("%d-%b-%Y %H:%M:%S")
                write_line("# {}".format(current_time), file)

    return


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
    files = " ".join(file_names)
    cmd = "grep '{}' {}".format(line, files)

    try:
        subprocess.check_output(cmd, shell=True)
        parsed = True
    except subprocess.CalledProcessError:
        parsed = False

    return parsed


def write_line(line, output_file):
    """ Function that writes the given string to new line in the give file

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


def obtain_taxid(metadata_file):
    """Determines the NCBI Taxonomic identifier of an accession using metadata

    Key arguments:
        metadata_file -- str, pathway to a NCBI-SRA metadata file format.

    Returns:
        taxid -- int, the TaxID belonging to the given metadata file.
    """
    if os.path.isfile(metadata_file):
        cmd_obtain_taxid = "cut -d',' -f28 {} | tail -1".format(metadata_file)
        taxid = subprocess.check_output(cmd_obtain_taxid,
                                        shell=True).decode("UTF-8").strip()

        if taxid.isnumeric():
            taxid = int(taxid)
        else:
            raise ValueError("TaxID could not be determined!")
    else:
        raise FileNotFoundError("Metadata file does not exist.\n"
                                "Please first download the metadata of all "
                                "included accessions!")

    return taxid


def count_taxids(accession_file, extension=None, base_loc=None):
    """Counts the occurrences of TaxID's based on multiple metadata files

    Key Arguments:
        accession_file -- str, pathway to a file containing NCBI-SRA
                          accessions on every new line.

    Returns
        taxids -- dict, containing the TaxID's as keys and a count as value
                  based on the number of times the TaxID occurs in
                  the metadata.
    """
    metadata_base = str(base_loc)

    # Set the default extension I use on the server, if not specified
    if not extension:
        extension = "_metadata.txt"
    else:
        extension = extension

    # Count the TaxIDs of all the given accessions
    taxids = {}
    counter = 0
    with open(accession_file, 'r') as input_file:
        for accession in input_file:
            meta_data = metadata_base + accession.strip() + extension
            identifier = obtain_taxid(meta_data)
            counter += 1
            if identifier in taxids.keys():
                taxids[identifier] += 1
            else:
                taxids[identifier] = 1

    input_file.close()

    # Sort the dictionary, TaxID with most counts on top
    taxids = dict(sorted(taxids.items(), key=lambda item: item[1],
                         reverse=True))

    return taxids, counter


# Functions that combine the other functions above based on the given cmd
def combine_librarylayout(file_name, metadata_list, accession, printing,
                          p_out=None, s_out=None):
    """ Combines all the functions that are needed to obtain the librarylayout

    Key Arguments:
        file_name -- str, pathway to the metadata file that has to be parsed.
        metadata_list -- lst, nested list containing two lists. First list is
                         the headers of the NCBI-SRA metadata file and the
                         second list contains the corresponding values.
        accession -- str, accession number belonging to the metadata.
        printing -- bool, True if the output needs to be printed, False if the
                    output should be writen to output files.
        p_out -- str, pathway to the file where paired-end accessions should be
                 appended to. Default = None for printing.
        s_out -- str, pathway to the file where single-end accessions should be
                 appended to. Default = None for printing.

    Returns:
        None
    """
    # Obtain the library layout and date
    librarylayout = obtain_librarylayout(file_name, metadata_list)
    now = datetime.now()
    current_time = now.strftime("%d-%b-%Y %H:%M:%S")

    # Print the accession if 'write' is false
    if printing:
        header = ["Accession", "Library Layout"]
        print(current_time)
        print("{:<22}{}\n{:<22}{}".format(header[0], header[1],
                                          accession, librarylayout))
    # Write to provided files if 'print' flag is not given by the user
    else:
        # Create the outputfiles, if not already present
        check_existence([p_out, s_out])
        # Check if the accession is already present in one of the files
        parsed = check_if_parsed(accession, [p_out, s_out])
        if not parsed:
            if librarylayout == "PAIRED":
                write_line(accession, p_out)
            elif librarylayout == "SINGLE":
                write_line(accession, s_out)
        else:
            print(Bcolors.WARNING +
                  "Note: current accession '{}' was already parsed before, "
                  "exiting the accession_parser.py script.".format(accession) +
                  Bcolors.ENDC)
            sys.exit(0)

    return


def combine_taxid(metadata_loc, accession, printing, output_file, ignore,
                  undefined_ids):
    """ Combines all the functions that are needed to obtain the taxid

    Key Arguments
        metadata_loc -- str, pathway to the file containing the metadata in
                        NCBI-SRA format.
        accession -- str, accession number belonging to the metadata.
        printing -- bool, True if output needs to be printed to the terminal,
                    False if output should be writen to the output file.
        output_file -- str, pathway to the output file to which the taxid's
                       should be writen.
        ignore -- bool, False if only unique taxid's should be writen to the
                  output file. True if all taxid's should be writen to the
                  output file, not checking if the taxid is already present.
                  Default only writes unique taxid's.
        undefined_ids -- str, pathway to the output file where the undefined
                         taxonomic identifiers should be writen to.

    Returns
        None
    """
    # Obtain the taxid
    try:
        taxid = int(obtain_taxid(metadata_loc))

        # Print to the screen if desired
        if printing:
            taxid_string = '{:<22}{}\n{:<22}{}'.format("# Accession number",
                                                       "Taxid", accession,
                                                       taxid)
            print(taxid_string)

        # Write to the provided output file if the '--print' flag is not given
        else:
            if not output_file:
                raise ValueError("No output file is given!")
            else:
                # Create the output file if not present
                check_existence([output_file])
            # Check if only unique taxid's can be writen
            if not ignore:
                if check_if_parsed(taxid, [output_file]):
                    print(Bcolors.WARNING +
                          'Note TaxID: {} is already present in the '
                          'output file, continuing.'.format(taxid) +
                          Bcolors.ENDC)
                    pass
                else:
                    write_line(taxid, output_file)
            else:
                write_line(taxid, output_file)
                print(Bcolors.OKCYAN +
                      'Forced appending of TaxID: {}, '
                      'not checking for duplicates!'.format(taxid) +
                      Bcolors.ENDC)
    except ValueError:
        taxid = obtain_taxid(metadata_loc)
        if printing:
            print("The TaxID of accession: '{}' could not be determined for "
                  "some reason. The following TaxID was found: '{}'."
                  .format(accession, taxid))
        else:
            write_line("Accession: '{}' with undefined TaxID: '{}'"
                       .format(accession, taxid), undefined_ids)

    return


def main():
    """Function that wraps all functions in the script"""
    # Step 0.0: check if the command line input is given
    parser = parsing_cmd_line()
    args = parser.parse_args()

    if len(sys.argv) < 3 and '-h' and '--h' not in args:
        print(Bcolors.WARNING + "No commandline arguments were given!" +
              Bcolors.ENDC)
        print(parser.print_help(sys.stderr))
        sys.exit()

    # Step 0.1: first check if accession details are on the correct position
    #           and obtain the accession ID
    try:
        if args.metadata:
            # Check format
            check_format(read_metadata(args.metadata), args.metadata)

            # Read metadata and store accession ID
            meta = read_metadata(args.metadata)
            accession = meta[1][0]
    except AttributeError:
        pass

    # TODO: CREATE THIS MAIN LIBRARY FUNCTION INTO A SEPARATE MODULE PLEASE!!
    # library: print or write the library layout of the metadata
    if "library" in args.command:
        if args.print:
            combine_librarylayout(args.metadata, meta, accession,
                                      args.print, args.paired, args.single)
        elif args.paired:
            if not args.single:
                print(Bcolors.FAIL +
                      "WARNING: when -p option is given, "
                      "-s option is required!" +
                      Bcolors.ENDC, parser.print_help())
                sys.exit(0)
            else:
                combine_librarylayout(args.metadata, meta, accession,
                                      args.print, args.paired, args.single)
        else:
            print(Bcolors.FAIL +
                  "WARNING: when library subcommand is chosen, option -p and "
                  "-s are required!" +
                  Bcolors.ENDC, parser.print_help())
            sys.exit(0)

    # taxid: print or write a single TaxID from a single metadata file
    if "taxid" in args.command and not "count" in args.command:
        "Obtaining the specific taxid from one file."
        combine_taxid(args.metadata, accession, args.print, args.output,
                      args.all, args.undefined)

    # Counting Taxonomic Identifiers (TaxID's) based on multiple input files
    if "taxid_count" in args.command:
        taxid_stats, nr_files = count_taxids(args.accessions,
                                             args.extension,
                                             args.metadata_base)
        header = "# Original file with " \
                 "accessions: {}\n".format(args.accessions) + \
                 "# Total number of provided accessions: " \
                 "{}".format(nr_files)
        now = datetime.now()
        current_time = now.strftime("%d-%b-%Y %H:%M:%S")

        # print to the screen if global flag '--print' is given
        if args.print:
            print("# " + current_time)
            print(header)
            for key in taxid_stats:
                print("{}\t{}".format(key, taxid_stats[key]))
        else:
            # Write to a specified output file, given by flag '-o/--output'
            # Create the output file if not present
            if not args.count_outfile:
                raise argparse.ArgumentTypeError(Bcolors.FAIL +
                                                 "No output file name or path "
                                                 "was specified using the "
                                                 "-o/--output flag. If the "
                                                 "global flag '--print' is "
                                                 "not used, -o is required!" +
                                                 Bcolors.ENDC)
            else:
                check_existence([args.count_outfile])
                write_line(header, args.count_outfile)
                for write_key in taxid_stats:
                    line = "{}\t{}".format(write_key, taxid_stats[write_key])
                    write_line(line, args.count_outfile)

    # kmer: create a string of K-mers with max 75% of read length.
    if "kmer" in args.command:
        print(determine_kmer(args.fastq, args.kmer_string))

    # genome: determine if a genome assembly is available for given TaxID('s)
    # or find the name of the genome assembly for the given TaxID('s).
    if "genome" in args.command:
        if not args.print and not args.ref_output:
            print(Bcolors.FAIL +
                  "WARNING: when global flag '--print' is not given, option -o"
                  " is required!" +
                  Bcolors.ENDC, parser.print_help())
            sys.exit(1)
        elif args.print and args.ref_output:
            print(Bcolors.FAIL +
                  "WARNING: when global flag '--print' is given, option -o "
                  "is not allowed! Please use either one." +
                  Bcolors.ENDC, parser.print_help())
            sys.exit(1)
        elif not args.print:
            if all(value is not None for value in [args.ref_output,
                                                   args.no_ref_output,
                                                   args.no_assem_output]):
                genome_main(args.taxids, args.print, args.quiet,
                            args.ref_output, args.no_ref_output,
                            args.no_assem_output)
            else:
                print(Bcolors.FAIL +
                      "WARNING: please specify all three output files with "
                      "options -o1, -o2 and -o3." +
                      Bcolors.ENDC)
        elif args.print:
            genome_main(args.taxids, args.print, args.quiet,
                        args.ref_output, args.no_ref_output,
                        args.no_assem_output)

    return


if __name__ == '__main__':
    main()
