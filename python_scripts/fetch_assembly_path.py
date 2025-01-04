#!/usr/bin/env python3
"""Author: Devin van Valkengoed

Date: 12-Dec-2022 
Description: This script is imported as a helper module for main module
'meta_accesion_parser.py' and includes the following functions:

     * obtain_assembly_metadata -- returns the metadata of the RefSeq genome
                                   assembly of a given Taxonomic Identifier.
     * parse_assembly -- returns the assembly accession and assembly name,
                         using the assembly metadata.

Raises:
    NameError: when the total count of the metadata = 0 and thus it is very
               likely that there is no RefSeq assembly for the given TaxID.

Notes:
    * To run this module, NCBI Datasets commandline tools has to be installed!
"""
# Import statements
import subprocess
import sys


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
def obtain_assembly_metadata(taxid):
    """Fetches the metadata of a genome assembly using the given TaxID

    Key arguments:
        taxid -- int, species specific NCBI Taxonomic Identifier (TaxID).

    Returns:
        assembly_metadata -- str, metadata of the genome assembly belonging to
                             the RefSeq assembly of the given TaxID.
    """
    cmd = "datasets summary genome taxon {} --reference".format(taxid)
    try:
        assembly_metadata = subprocess.check_output(cmd,
                                                    stderr=subprocess.STDOUT,
                                                    shell=True)\
            .decode("UTF-8").strip()
        return assembly_metadata
    except subprocess.CalledProcessError as error:
        if error.returncode == 127:
            raise OSError("Bash command not found. Is NCBI datasets cli "
                          "installed and active? Please check the path "
                          "variable and make sure optional Conda env is "
                          "active!") from None
        else:
            print(error.output)
            sys.exit(error.returncode)


def parse_assembly(taxid, assembly_metadata):
    """ Parses the assembly accession and name from assembly metadata

    Key arguments:
        taxid -- int, species specific NCBI Taxonomic Identifier (TaxID).
        assembly_metadata -- str, metadata of the genome assembly belonging to
                             the RefSeq assembly of the given TaxID.

    Returns:
        assembly_name -- str, name of the genome assembly.
        assembly_accession -- str, accession of the genome assembly.
    """
    # Split the metadata and set empty values for name and accession
    clean_assembly = assembly_metadata.split(",")
    assembly_name = None
    assembly_accession = None

    # Throw an error if the metadata is only one element long (no RefSeq)
    if len(clean_assembly) == 1:
        for ele in clean_assembly:
            if int(ele.strip("{}").split(":")[1].strip()) == 0:
                raise NameError("No RefSeq assembly is available for TaxID: "
                                "{}. Please only provide a list of TaxIDs "
                                "with RefSeq genome assembly!")

    # Determine the assembly name and accession
    for element in clean_assembly:
        # Obtain the assembly accession
        if "current_accession" in element:
            assembly_accession = str(element.split(":")[1]).strip("\"")
        # Obtain the assembly name
        elif "assembly_name" in element and "organelle" not in element:
            assembly_name = str(element.split(":")[1]).strip("\"")

    # Throw an error if the assembly name or accession could not be found
    if not assembly_name:
        print(Bcolors.FAIL +
              "Assembly name could not be determined! Please double check the "
              "metadata of TaxID {} manually.".format(taxid) +
              Bcolors.ENDC)
        sys.exit(1)
    if not assembly_accession:
        print(Bcolors.FAIL +
              "Assembly accession could not be determined! Please double "
              "check the metadata of TaxID {} manually.".format(taxid) +
              Bcolors.ENDC)
        sys.exit(1)

    return assembly_accession, assembly_name


if __name__ == "__main__":
    """This is a helper module for meta_accession_parser.py"""
    pass
