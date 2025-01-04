#!/usr/bin/env python3
"""Author: Devin van Valkengoed 

Date: 15-Nov-2022
Description: This script is imported as a helper module for main module
'meta_accesion_parser.py' and includes the following functions:

     * parse_taxids -- returns a list of taxids read from every new line in a
                       provided file.
     * sort_on_reference_assembly -- checks availability of a genome assembly
                                     for a given taxid.

Notes:
    * To run this module, NCBI Datasets commandline tools has to be installed!


TODO:
    * Check for /bin/sh: 1: datasets: not found error and print a message that
      the conda env likely is not active.
"""
# Import statements
import sys
import subprocess
import os


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


# Functions
def parse_taxids(taxid_file):
    """Stores TaxID's each provided on a new line in a file, in a list

    Key Arguments:
        taxid_file -- str, pathway to a .txt file containing single TaxID's on
                      each line of the file.

    Returns:
          taxids -- lst, list containing all the TaxID's. (not unique)

    """
    taxids = []
    with open(taxid_file, 'r') as input_file:
        for line in input_file:
            if not line:
                continue
            else:
                line = line.strip()
                if line.startswith("#"):
                    continue
                else:
                    taxids.append(line)

    input_file.close()

    return taxids


def sort_on_reference_assembly(taxid_lst, reference=True):
    """Checks if there is a genome available for the given TaxID's

    Key Arguments:
        taxid_lst -- lst, containing one or multiple NCBI TaxID's.
        reference -- bool, True if the function should check the availability
                     of a reference genome assembly. False, if it should check
                     the availability of all assemblies.

    Returns:
        ref_genome -- lst, TaxID's that have a reference genome assembly
                      available.
        noref_genome -- lst, TaxID's that do not have a reference genome.

    TODO:
        * Print coloured cmd argument when error occurs?
    """
    # Determine if the reference genome assembly or all assemblies are checked
    if reference:
        cmd = "datasets download --no-progressbar genome --preview taxon " \
              "--reference "
    else:
        cmd = "datasets download --no-progressbar genome --preview taxon"

    # Store the assembly metadata of each provided TaxID:
    taxid_metadata = {}
    for taxid in taxid_lst:
        try:
            stdout = subprocess.check_output("{} {}".format(cmd, taxid),
                                             shell=True,
                                             stderr=subprocess.DEVNULL) \
                .decode("UTF-8")
            taxid_metadata[taxid] = stdout
        except subprocess.CalledProcessError as err:
            if err.returncode == 127:
                print(Bcolors.FAIL +
                      "NCBI Datasets cli ran into error 127. Did you turn on "
                      "your conda env? Exiting the script." +
                      Bcolors.ENDC)
                sys.exit(127)
            else:
                raise ValueError("Datasets cli exited with {} exit "
                                 "code".format(err.returncode)) \
                    from None

    # Separate TaxID's based on if there is a reference genome available
    ref_genome = []
    noref_genome = []
    for identifier in taxid_metadata:
        if int(taxid_metadata
               [identifier].strip().split(',')[1].split(':')[1]) > 0:
            ref_genome.append(identifier)
        else:
            noref_genome.append(identifier)

    return ref_genome, noref_genome


def check_if_parsed(line, file_names):
    """Determines if the current line already parsed and present in a file

    Key arguments:
        line -- str, accession number of the current provided metadata.
        file_names -- lst, list of strings containing pathways to the files
                      that should be checked for the presence of the provided
                      line. Names of different files should be separated by
                      spaces.

    Returns:
        parsed -- bool, True if the line is found in one of the files, False if
                  not.
    """
    file_names = file_names.split(" ")
    if len(file_names) > 1:
        files = " ".join(file_names)
    else:
        files = str(file_names[0])
    cmd = "grep -w '{}' {}".format(line, files)

    try:
        subprocess.check_output(cmd, shell=True)
        parsed = True
    except subprocess.CalledProcessError:
        parsed = False

    return parsed


def genome_main(taxid_file, printing, quiet, refseq_out, assembly_out,
                no_assembly_out):
    """ Temp function that includes main function of check_genome_assemblies.py
    """
    taxid_list = parse_taxids(taxid_file)

    if not quiet:
        print("Checking for TaxID genome assembly availability, "
              "please wait.")
        print("Given input file with TaxIDs: '{}'".format(taxid_file))
    reference_assembly, no_ref = sort_on_reference_assembly(taxid_list,
                                                            True)
    assembly, no_assembly = sort_on_reference_assembly(no_ref, False)

    # Print the output to the screen if global flag '--print' is given
    if printing:
        print("TaxID's with RefSeq assembly: {}\n"
              "TaxID's with general assembly (no RefSeq): {}\n"
              "TaxID's with no assembly: {}\n"
              "Number of TaxID's with a RefSeq genome assembly: {}\n"
              "Number of TaxID's with no RefSeq, but assembly: {}\n"
              "Number of TaxID's with no assembly: {}"
              .format(reference_assembly,
                      assembly,
                      no_assembly,
                      len(reference_assembly),
                      len(assembly),
                      len(no_assembly)))
    # Write to the provided output file if the '--print' flag is not given
    else:
        if not quiet:
            print("Written the TaxID's to the following files: \n"
                  "TaxIDs with RefSeq assembly: {}\n"
                  "TaxIDs with general assembly (but no RefSeq): {}\n"
                  "TaxIDs without any assembly: {}".format(refseq_out,
                                                           assembly_out,
                                                           no_assembly_out))

        # Loop through every list and write to the correct file.
        # Also check for duplicates.
        assemblies = [reference_assembly, assembly, no_assembly]
        out_files = [refseq_out, assembly_out, no_assembly_out]

        for idx, element in enumerate(assemblies):
            if len(element) > 0:
                for taxid in element:
                    if not os.path.isfile(out_files[idx]):
                        with open(out_files[idx], 'a+') as output_file:
                            output_file.write(taxid + "\n")
                    elif not check_if_parsed(taxid, out_files[idx]):
                        with open(out_files[idx], 'a+') as output_file:
                            output_file.write(taxid + "\n")
                    elif check_if_parsed(taxid, out_files[idx]):
                        if not quiet:
                            print("\tTaxID: '{}' is already present in: '{}'"
                                  .format(taxid, out_files[idx]))
                        else:
                            continue
            else:
                # List with TaxIDs is empty
                continue

    return


if __name__ == '__main__':
    """This is a helper module for meta_accession_parser.py"""
    pass
