#!/usr/bin/env python3
"""Author: Devin van Valkengoed

Date: 06-Jun-2023
Description: This script can be used to determine the viral completeness of
             assembled genomes using DIAMOND files. A config file is
             required to set viral family specific cut-offs.

Usage: python3 viral_completeness.py -c <config_file> {viral_completeness,
       patching}

Key arguments:
    config_file -- str, pathway to a YAML file that contains all the user
                   set parameter information.
    viral_completeness -- str, subcommand to determine the viral
                          completeness of the given accessions.
    patching -- str, subcommand to perform patching of the contigs found
                for the given accessions.

Run: 'python3 viral_completeness.py viral_completeness/patching -h' to show
     the possible options for every subcommand.
"""
# Import statements
import argparse
import sys
import os.path
import yaml
import subprocess
import itertools
import re
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
from datetime import datetime


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


# Commandline, config and log file functions
def parsing_cmd():
    """Uses argparse to determine the given commandline options

    Returns:
        parser -- argparse.Namespace object
    """
    parser = argparse.ArgumentParser(
        description="This script can be used to 'patch' fragmented genome "
                    "assemblies. The patching is done by trimming and joining "
                    "contigs based on DIAMOND2 nt and protein alignments. "
                    "A config file is required.")
    parser.add_argument("-c", "--config",
                        help="The relative pathway to a config file that "
                             "contains all information about file locations "
                             "user settings. File should be in YAML format.",
                        metavar="<./config.yml>", dest="config", required=True)
    parser.add_argument("-s", "--safe", help="This mode will remove "
                                             "accessions from which the "
                                             "metadata or diamond files are "
                                             "not available. It will prevent "
                                             "the script from running into "
                                             "an error.",
                        action="store_true", dest="safe")
    parser.add_argument("-p", "--print", help="This option will print "
                                              "extra information to stdout "
                                              "during the process. ",
                        action="store_true", dest="print")
    subparser = parser.add_subparsers(title="Subcommands",
                                      description="Valid subcommands",
                                      dest="main_parser")

    # Subcommands
    # Viral completeness (vc)
    vc_parser = subparser.add_parser("viral_completeness",
                                     help="This option will determine the "
                                          "viral completeness of every "
                                          "accession given in the config "
                                          "file based on the parameters "
                                          "provided in the config file. ")
    vc_parser.add_argument("-v", "--viral_fam",
                           help="Names of the viral families that should be "
                                "included in patching and/or while "
                                "determining viral completeness. "
                                "Multiple families can be given, separated by "
                                "spaces. E.g. iflaviridae secoviridae. "
                                "By default the six main families within the "
                                "Picornavirales order will be included: "
                                "Iflaviridae, Dicistroviridae, Marnaviridae, "
                                "Picornaviridae, Secoviridae and "
                                "Calciviridae.",
                           metavar="<included_families>", dest="included_fam",
                           nargs="+")

    # Genome patching
    patching_parser = subparser.add_parser("patching",
                                           help="This option will try to "
                                                "patch the genomes of every "
                                                "diamond file that contains "
                                                "fragmented genomes.")
    patching_parser.add_argument("-v", "--viral_fam",
                                 help="Names of the viral families that "
                                      "should be included in patching "
                                      "and/or while determining viral "
                                      "completeness.  Multiple families can "
                                      "be given, separated by spaces. E.g. "
                                      "iflaviridae secoviridae. By default "
                                      "the six main families within the "
                                      "Picornavirales order will be included: "
                                      "Iflaviridae, Dicistroviridae, "
                                      "Marnaviridae, Picornaviridae, "
                                      "Secoviridae and Calciviridae.",
                                 metavar="<included_families>",
                                 dest="included_fam",
                                 nargs="+")
    patching_parser.add_argument("-r", "--reset",
                                 help="If this option is given, "
                                      "any output files will created again.",
                                 action="store_true", dest="reset")

    return parser


def config_parsing(config_yaml):
    """Parses the information from a config file in YAML format

    Key arguments:
        config_yaml -- str, relative pathway to a config file writen in YAML
                       format

    Returns:
        configs -- dict, containing the elements from the YAML file.
    """
    with open(config_yaml, "r") as config_file:
        configs = yaml.safe_load(config_file)

    return configs


def line_parser(input_file):
    """Parses every new line of a file as element in a list

    Key arguments:
        input_file -- str, pathway to a .txt file containing information on
                      every new line.

    Returns:
        line_list -- lst, list containing all the lines of the file as
                     elements.
    """
    line_list = []
    with open(input_file, 'r') as in_file:
        for line in in_file:
            line = line.strip()
            if line == "":
                continue
            if not line:
                continue
            if line.startswith("#"):
                continue
            else:
                line_list.append(line)

    return line_list


# Functions ###################################################################
def obtain_taxid(metadata_file):
    """Determines the NCBI Taxonomic identifier of an accession using metadata

    Key arguments:
        metadata_file -- str, pathway to a NCBI-SRA metadata file format.

    Returns:
        taxid -- int, the TaxID belonging to the given metadata file.
    """
    cmd_obtain_taxid = "cut -d',' -f28 {} | tail -1".format(metadata_file)
    taxid = subprocess.check_output(cmd_obtain_taxid, shell=True).decode(
        "UTF-8").strip()

    if taxid.isnumeric():
        taxid = int(taxid)
    else:
        raise ValueError("TaxID of: '{}' could not be determined. "
                         "Obtained: '{}' which is not possible as "
                         "TaxID.".format(metadata_file, taxid))

    return taxid


def parse_diamond_files(diamond, taxid):
    """Reads DIAMOND files, in tabular format, into memory

    Key arguments:
        diamond -- str, pathway to a DIAMOND file in BLAST tabular
                   format (-f 6).
        taxid -- str, taxonomic identifier to include in the name of the nodes.

    Returns:
        dmnd_file -- dict, containing the query name as key and the other
                     columns of the DIAMOND file in a list, as value.
    """
    dmnd_file = {}
    with open(diamond, 'r') as input_file:
        for line in input_file:
            if not line:
                continue
            line = line.strip().split("\t")
            key_name = "{}_{}&against&{}".format(line[0], taxid, line[1])
            dmnd_file[key_name] = line[1:]

    return dmnd_file


def filter_length(diamond_dict, min_length, max_length=None):
    """Filters given diamond alignments based on the original contig length

    Key arguments:
        diamond_dict -- dict, containing the query name as key and the other
                        columns of the DIAMOND file in a list, as value.
        min_length -- int, minimum length of the original contig (NT) to be
                      included in the output.
        max_length -- int, maximum length of the original contig (nt) to be
                      included in the output. Default = None and thus no max
                      length filter will be applied.

    Returns:
        length_filtered_alignments -- dict, containing the query name as key
                                      and the other columns of the DIAMOND
                                      file in a list, as value. Only alignments
                                      with an original contig length >
                                      min_length are stored.
        below_min_length -- dict, containing the query name as key and the
                            other columns of the DIAMOND file in a list,
                            as value. Only alignments with an original
                            contig length < min_length are stored.
        above_max_length -- dict, containing the query name as key and the
                            other columns of the DIAMOND file in a list,
                            as value. Only alignments with an original
                            contig length > max_length are stored.
    """
    length_filtered_alignments = {}
    below_min_length = {}
    above_max_length = {}

    for header in diamond_dict:
        original_header = header.split("&")
        length = original_header[0].split("_")[3]
        if not max_length:
            if float(length) >= min_length:
                # Add to filtered dictionary
                length_filtered_alignments[header] = diamond_dict[header]
            else:
                # Add to below minimum ignored
                below_min_length[header] = diamond_dict[header]
        else:
            if min_length <= float(length) <= max_length:
                # Add to filtered dictionary
                length_filtered_alignments[header] = diamond_dict[header]
            elif min_length > float(length):
                # Add to below minimum ignored
                below_min_length[header] = diamond_dict[header]
            elif max_length < float(length):
                # Add to above max ignored
                above_max_length[header] = diamond_dict[header]

    return length_filtered_alignments, below_min_length, above_max_length


def filter_alignment_product(diamond_dict, min_product):
    """Filters diamond alignments based on a minimum '%ID * alignment_len'

    Key arguments:
        diamond_dict -- dict, containing the query name as key and the
                        other columns of the DIAMOND file in a list, as value.
        min_product -- int, minimum product of the %AID * alignment length for
                       the alignments to pass the filter.

    Returns:
        product_filtered_alignments -- dict, containing the query name as key
                                       and the other columns of the DIAMOND
                                       file in a list, as value. Only
                                       alignments with a product of '%ID *
                                       alignment_len' > min_product are stored.
        nr_product_ignored -- int, number of alignments that are discarded from
                              the original diamond file, based on the given
                              minimum product.
    """
    product_filtered_alignments = {}
    nr_product_ignored = 0

    for key in diamond_dict:
        aid_product = float(diamond_dict[key][1]) * int(diamond_dict[key][2])
        if aid_product > min_product:
            product_filtered_alignments[key] = diamond_dict[key]
        else:
            nr_product_ignored += 1

    return product_filtered_alignments, nr_product_ignored


def fetch_largest_product(diamond_dict):
    """Obtains alignment with the largest product of sum(length) * %identity

    Key arguments:
        diamond_dict -- dict, containing the query name as key and the other
                        columns of the DIAMOND file in a list, as value.

    Returns:
        top_alignment -- lst, finds the alignment with the largest product
                         of 'contig assembly length' * '% alignment ID' *
                         'alignment length'.
    """
    # Calculate the product of the alignment %id * contig length
    for key in diamond_dict:
        assembly_len = float(key.split("&")[0].split("_")[3])
        identity = float(diamond_dict[key][1])
        align_len = float(diamond_dict[key][2])
        custom_product = float(identity * assembly_len * align_len)

        diamond_dict[key].append(custom_product)

    # Create a list that only contains the alignment key and the product
    list_of_alignments = []
    for alignment in diamond_dict:
        list_of_alignments.append((alignment, diamond_dict[alignment][11]))

    # Sort the list in descending order
    sorted_alignments = sorted(list_of_alignments, key=lambda x: x[1],
                               reverse=True)

    # Find the alignment with the largest product in the original dictionary
    top_alignment = [sorted_alignments[0][0],
                     diamond_dict[sorted_alignments[0][0]]]

    return top_alignment


def fetch_unique_nodes(diamond_dict):
    """Finds the unique nodes in a given diamond alignment dictionary

    Key arguments:
        diamond_dict -- dict, containing the query name as key and the other
                        columns of the DIAMOND file in a list, as value.

    Returns:
        unique_nodes -- dict, containing unique nodes as keys and an integer
                        that gives the number of alignments for that node, as
                        value. E.g. {'NODE_2': 25} shows that the alignment
                        file contains one unique node (NODE_2) with 25
                        different alignments.
    """
    unique_nodes = []

    for key in diamond_dict:
        node = key.split("&")[0]
        if node not in unique_nodes:
            unique_nodes.append(node)

    return unique_nodes


def get_key(input_dict, element):
    """Function wil obtain the key belonging to the given value present in lst

    Key arguments:
        input_dict -- dict, containing one or multiple keys, with values
                      being lists of elements. These lists will be checked
                      for the given value.
        element -- str, object that should be found in the value lists,
                   inside the dictionary.

    Returns:
        key -- str, if the object is found in the value list,
               the corresponding key wil be returned.
        unclassified -- str, when the element cannot be found inside the
                        dictionary the string: 'unclassified' will be returned.
    """
    for key, val in input_dict.items():
        if element in val:
            return key

    return "unclassified"


def count_target_occurrence(diamond_dict, target):
    """Determines the number of times a given target occurs in a dmnd file

    Key arguments:
        diamond_dict -- dict, containing the query name as key and the other
                        columns of the DIAMOND file in a list, as value.
        target -- str, target identifier which should be counted throughout
                  the diamond file.

    Returns:
         target_count -- int, the number of occurrences of the target
                         throughout the complete diamond file.
    """
    target_count = 0
    for alignment in diamond_dict:
        if diamond_dict[alignment][0] == target:
            target_count += 1

    return target_count


def obtain_target_node(diamond_dict, target):
    """Obtains all node names that have an alignment against the given target

    Key arguments:
        diamond_dict -- dict, containing the query name as key and the other
                        columns of the DIAMOND file in a list, as value.
        target -- str, target identifier for which all query names should be
                  found within the given diamond_dict.

    Returns:
        query_list -- lst, containing strings of the query names that have
                      an alignment against the given target.
    """
    query_list = []
    for alignment in diamond_dict:
        if diamond_dict[alignment][0] == target:
            query_list.append(alignment)

    return query_list


def obtain_target_alignments(diamond_dict, target):
    """Obtains the complete alignment stats for all alignments against target

    Key arguments:
        diamond_dict -- dict, containing the query name as key and the other
                        columns of the DIAMOND file in a list, as value.
        target -- str, the identifier of the target for which all parameters
                  of every alignment should be obtained.

    Returns:
        target_alignments -- dict, nested, containing the alignment queryID
                             as key and a dict as value. This dict contains
                             the key 'dmnd_params' with a list of the
                             diamond alignment parameters as value.
    """
    query_alignments = {}
    for alignment in diamond_dict:
        if diamond_dict[alignment][0] == target:
            query_alignments[alignment] = {"dmnd_params":
                                           diamond_dict[alignment][:-1]}

    return query_alignments


def grep_nt_sequence(nodes, fasta_file, output_file):
    """Obtains the original nt sequences from contigs.fasta files using sed

    Key arguments:
        nodes -- lst, contig header names of which the nt sequence should be
                 obtained.
        fasta_file -- str, pathway to a fasta file which contains the
                      contigs that are provided in the nodes list.
        output_file -- str, output file to which the fasta sequences should
                       be writen to.

    Returns:
        None
    """
    for alignment in nodes:
        original_contig_name = "_".join(alignment.split("&")[0].
                                        split("_")[:-1])
        sed_cmd = "sed -n '/{}/,/>/p' {} | " \
                  "head -n -1 | tail +2".format(original_contig_name,
                                                fasta_file)

        # Use sed to obtain the sequence of every contig
        sequence = subprocess.check_output(sed_cmd,
                                           shell=True).decode("UTF-8").strip()
        # Replace newline characters that sed includes
        sequence = sequence.replace("\n", "")

        # Store the sequence information in a BioPython SeqRecord object
        record = SeqRecord(
            Seq(sequence),
            id=original_contig_name,
            description=""
        )

        # Write the sequences to the same fasta file
        with open(output_file, "a+") as fasta_out:
            SeqIO.write(record, fasta_out, "fasta")

    return


def grep_prot_seq(node, fasta_file):
    """Obtains a protein sequence from a fasta file based on the header name

    Key arguments:
        node -- str, (part of the) header name of the sequence that should be
                obtained from the fasta file.
        fasta_file -- str, pathway to the fasta file containing the desired
                      protein sequences.

    Returns:
        protein_lst -- lst, the protein sequence in a list, containing a
                       single aminoacid on every position.
    """
    # Use sed to obtain the desired sequence from the file
    sed_cmd = "sed -n '/{}/,/>/p' {} | " \
              "head -n -1 | tail +2".format(node, fasta_file)

    prot_sequence = subprocess.check_output(sed_cmd,
                                            shell=True).decode("UTF-8").strip()

    # Replace newline characters that sed includes
    prot_sequence = prot_sequence.replace("\n", "")

    # Put the sequence in a list
    prot_lst = []
    for pos in prot_sequence:
        prot_lst.append(pos)

    return prot_lst


def find_multi_orf_contigs(fasta_input):
    """Finds contigs that have multiple open reading frames (ORFs)

    Key arguments:
        fasta_input -- str, pathway to a fasta file containing protein
                       sequences.

    Returns:
         multi_orf -- bool, True if the fasta file contains multiple ORFs
                            for a single contig. False if every contig has
                            a single ORF.
    """
    unique_fasta_headers = []
    with open(fasta_input, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            header = "_".join(record.id.split("_")[:-1])
            if header in unique_fasta_headers:
                return True
            else:
                unique_fasta_headers.append(header)

    return False


def find_longest_unique(fasta_input):
    """Removes shorter duplicate sequences from a fasta file, based on header

    Key arguments:
        fasta_input -- str, pathway to a fasta file containing nt or protein
                       sequences.

    Returns:
        None

    Will write the deduplicated fasta sequences (based on header name) to a
    txt file in fasta format. If output_fasta=None, then the original file
    will be overwritten.
    """
    # Create a list of tuples, with (sequence_names, length)
    length_info = []
    with open(fasta_input, "r") as fasta_in:
        for record in SeqIO.parse(fasta_in, "fasta"):
            length_info.append(("_".join(record.id.split("_")[:2]),
                                len(record.seq),
                                "_".join(record.id.split("_")[2:])))

    # Sort the list descending on length
    sorted_length = sorted(length_info, key=lambda x: x[1], reverse=True)

    # Determine the longest unique sequences
    unique_sequences = {}
    for sequence in sorted_length:
        if sequence[0] in unique_sequences.keys():
            continue
        else:
            unique_sequences[sequence[0]] = (sequence[1],
                                             sequence[2])

    # Create a list of unique original IDs
    original_header = []
    for seq in unique_sequences:
        original_header.append(str(seq + "_" + unique_sequences[seq][1]))

    # Write the original unique sequences to the output fasta file
    records = list(SeqIO.parse(fasta_input, "fasta"))

    # Empty the original input fasta file
    with open(fasta_input, "w") as fasta_file:
        pass

    # Write the unique contigs to the empty input fasta
    with open(fasta_input, "a+") as output_fasta:
        for rec in records:
            if rec.id in original_header:
                SeqIO.write(rec, output_fasta, "fasta")

    return


def print_viral_completeness(accession, genome_nr, genome_completeness,
                             viral_family, best_alignment):
    """Generalised print function for the check_viral_completeness function

    Key arguments:
        accession -- str, NCBI accession ID.
        genome_completeness -- str, the completeness of the obtained genome.
        viral_family -- str, the viral family of the best target.
        best_alignment -- str,best alignment containing all DIAMOND format
                          6 output parameters.

    Returns:
        None
    """
    print("{}\t{}\t{}\t{}\t{}".format(accession, genome_nr,
                                      genome_completeness,
                                      viral_family, best_alignment))

    return


def obtain_orf(input_file, prot_file, min_orf, log_file, check=False):
    """Obtains the Open Reading Frame(s) of (a) given nucleotide sequence(s)

    Key arguments:
        input_file -- str, pathway to the file containing a nucleotide sequence
                      in fasta format.
        prot_file -- str, pathway to the file where the protein sequence
                       should be written to.
        min_orf -- int, minimum length in NT for the orfs to be included.
                   Note: this limit should be set in such a way that only 1
                   orf is returned by EMBOSS getorf.
        log_file -- str, absolute pathway to a log file to write messages to.

    Optional arguments:
        check -- bool, if True, this will check that number of output
                 protein sequences is equal to the number of input nucleotide
                 sequences that are given. Default = False.

    Returns:
        None

    Please note that this function will only find the largest possible ORF
    for every given input nt sequence.
    """
    cmd_get_orf = "getorf -sequence {} -minsize {} -auto " \
                  "-outseq {}".format(input_file, str(min_orf),
                                      prot_file)
    subprocess.check_output(cmd_get_orf, shell=True)

    # Optional: check for duplicates and the nr of input and output sequences
    if check:
        # If there are multiple ORFs for a single nt sequence, deduplicate:
        if find_multi_orf_contigs(prot_file):
            c_print("Some provided nt sequences have multiple ORFs, "
                    "now obtaining the largest ORF for every sequence.",
                    None, log_file)
            find_longest_unique(prot_file)

        # Check if the nr of input nt sequences is equal to the nr of prot seqs
        c_print("Only a single ORF present for every given nt sequence. "
                "Double checking the nr of input nt sequences and nr of "
                "output prot sequences.",
                None, log_file)
        cmd_check_nr_input = "grep '>' {} | wc -l".format(input_file)
        nr_input_seqs = \
            int(subprocess.check_output(cmd_check_nr_input,
                                        shell=True).decode("UTF-8"))

        cmd_check_nr_output = "grep '>' {} | wc -l".format(prot_file)
        nr_prot_seqs = \
            int(subprocess.check_output(cmd_check_nr_output,
                                        shell=True).decode("UTF-8"))

        if nr_input_seqs > nr_prot_seqs:
            c_print("\nWARNING: not all nucleotide sequences could be "
                    "translated to a protein sequences. This is because of "
                    "the 'min_prot_len_nt' cut-off in the config "
                    "file. Hence, some nucleotide sequences were excluded "
                    "for patching!\n",
                    None, log_file)
        if nr_input_seqs < nr_prot_seqs:
            raise ValueError(Bcolors.FAIL +
                             "For some input nucleotide sequences, multiple "
                             "ORFs were found and this could not be "
                             "resolved! Number of input sequences: '{}', "
                             "number of translated protein sequences: '{}'"
                             .format(nr_input_seqs, nr_prot_seqs) +
                             Bcolors.ENDC)

    return


def run_diamond_blastp(query_fasta, diamond_db, output_file, cores):
    """This function uses DIAMOND2 blastp option to align protein sequences

    Key arguments:
        query_fasta -- str, pathway to a fasta file containing the protein
                       sequences that should be aligned against the given
                       database.
        diamond_db -- str, pathway to a diamond database containing protein
                      sequences.
        output_file -- str, pathway to a file where the output of the
                       diamond alignment should be stored.
        cores -- int, the number of cores to use for parallel processing.

    Returns:
        None

    DIAMOND2 output is written to the file specified in the output_file
    path. The output will be stored in the DIAMOND -f 6 format.
    """
    cmd_run_blastp = "diamond blastp -d {} -q {} -p {} -o {} " \
                     "-f 6 --quiet".format(diamond_db,
                                           query_fasta,
                                           cores, output_file)

    subprocess.check_output(cmd_run_blastp, shell=True)

    return


def fetch_prot_len(protein_id):
    """Determines the sequence length (aa) of a protein sequence on NCBI

    Key arguments:
        protein_id -- str, identifier of the protein of interest as stored by
                      NCBI protein db.

    Returns:
         protein_len -- int, length of the amino acid sequence of the
                        provided protein.

    Please note that the obtained length is based on the information given
    on the NCBI protein page and doesn't involve a sequence count.
    """
    prot_info_cmd = "efetch -db protein -id {} -format info | grep LOCUS"\
        .format(protein_id)

    prot_info = subprocess.check_output(prot_info_cmd, shell=True).decode(
        "UTF-8").strip().split()

    # Check if the obtained sequence is the right format
    if "LOCUS" not in prot_info[0] or "aa" not in prot_info[3]:
        raise ValueError("The format of the protein information is not "
                         "correct!")
    else:
        protein_len = int(prot_info[2])

    return protein_len


def initiate_dummy_lst(lst_length, fill="X"):
    """Creates a list with the number of elements provided

    Key arguments:
        lst_length -- int, the number of elements that should be placed in
                      the list. Every element is at a separate position.

    Optional arguments:
        fill -- str, the element that should be used to fill the list.
                Default = 'X'.

    Returns:
        dummy_list -- lst, containing the elements given for fill, and the
                      number of elements provided by lst_length.
    """
    dummy_list = [fill] * lst_length

    return dummy_list


def determine_gap_length(diamond_alignment):
    """Determines the gap length based on DIAMOND blast alignment information

    Key arguments:
        diamond_alignment -- lst, containing the DIAMOND blast tabular
                             format parameters (without the query name).

    Returns:
        gap_len -- int, the relative gap between the length of the target
                   and the query based on start and end positions of the
                   alignment.
    """
    # Determine start and end positions of query and target
    q_start = int(diamond_alignment[5]) - 1
    q_end = int(diamond_alignment[6])
    t_start = int(diamond_alignment[7]) - 1
    t_end = int(diamond_alignment[8])

    # Determine the alignment length of the query and target
    q_len = q_end - q_start
    t_len = t_end - t_start

    # Calculate gap length
    gap_len = t_len - q_len

    return gap_len


def c_print(message, colour=None, print_file=None):
    """Custom print function that can differentiate between stdout and logfile

    Key arguments:
        message -- str, the message that should be printed to stdout or
                   writen to the given file.
        colour -- str, one of the options of the dict headers below,
                  determining the optional colour of the message in stdout.
                  Default = None
        print_file -- str, optional file to write the message to.
                      Default = None

    Returns:
        None
    """
    bcolours_dict = {"HEADER": '\033[95m',
                     "OKBLUE": '\033[94m',
                     "OKCYAN": '\033[96m',
                     "OKGREEN": '\033[92m',
                     "WARNING": '\033[93m',
                     "FAIL": '\033[91m',
                     "ENDC": '\033[0m',
                     "BOLD": '\033[1m',
                     "UNDERLINE": '\033[4m',
                     "SKIP": ""}
    if not colour:
        colour = "SKIP"

    if not print_file:
        print(bcolours_dict[colour] +
              message +
              bcolours_dict["ENDC"])
    else:
        with open(print_file, "a+") as out_file:
            print(message,
                  file=out_file)

    return None


def initiate_patch_log(args, log_file, accession, config_dict,
                       config_file_loc):
    """Creates a plain text file to use as logfile during the patching

    Key arguments:
        args -- Namespace, argparse object containing all the arguments
                given on the commandline.
        log_file -- str, pathway and name of the file that should be created.
        accession -- str, NCBI-SRA identifier.
        config_dict -- dict, containing all parameters found in the
                       yml config file as a python dict object.
        config_file_loc -- str, the location to the config file that was used.

    Returns:
        None
    """
    if os.path.exists(log_file):
        rerun = True
    else:
        rerun = False

    with open(log_file, "a+") as out_file:
        # Set the current date and time
        now = datetime.now()
        current_time = now.strftime("%d-%b-%Y %H:%M:%S")
        out_file.write("\n"*2 + current_time + "\n")

        # Print message if this accession has been patched before
        if rerun:
            out_file.write("NOTE: THIS IS A RE-RUN, "
                           "PATCHING WAS DONE BEFORE!\n")

        # Initial logfile message
        start_message = "This is the patch log of accession: '{}'.\n"\
            .format(accession)
        second_line = "All information about the patching process will be " \
                      "stored in this file.\nPlease not that if the script " \
                      "is used again for this accession, the \nnew log " \
                      "information will be appended below.\n\n"
        out_file.write(start_message)
        out_file.write(second_line)

        # Store the given cmd arguments
        out_file.write("The commandline arguments that were provided:\n")
        out_file.write(str(args))
        out_file.write("\n\n")

        # Store all used file locations
        out_file.write("This config file was used during the "
                       "patching: '{}'.\nContaining the following "
                       "parameters: \n".format(config_file_loc))
        out_file.write("-" * 70 + "\n")
        for header in config_dict:
            params = "{}: {}\n".format(header, config_dict[header])
            out_file.write(params)

        out_file.write("-" * 70 + "\n\n")

    return None


def remove_target_alignments(diamond_dict, target_id):
    """Will remove all nodes that have an alignment against the target ID

    Key arguments:
        diamond_dict -- dict, containing the query name as key and the other
                        columns of the DIAMOND file in a list, as value.
        target_id -- str, the name of the target ID for which alignments
                     should be removed.

    Returns:
        updated_diamond_dict -- dict, containing the query name as key and
                                the other columns of the DIAMOND file in a
                                list, as value. The nodes having an
                                alignment against the given target will be
                                removed.
    """
    # Determine the keys that have an alignment against the given target
    keys2pop = []
    for query_id, alignment in diamond_dict.items():
        if alignment[0] == target_id:
            keys2pop.append(query_id)

    # Delete all nodes that have an alignment against the target
    for query in list(diamond_dict):
        if any("_".join(part.split("_")[:3]) in query for part in keys2pop):
            diamond_dict.pop(query)

    return diamond_dict


def check_viral_completeness(diamond_dict, config_dict, accession,
                             known_proteins):
    """Checks the completeness of a virus that is found based on dmnd file

    Key arguments:
        diamond_dict -- dict, containing the query name as key and the other
                        columns of the DIAMOND file in a list, as value.
        config_dict -- dict, containing all parameters found in the
                       yml config file as a python dict object.
        accession -- str, the original NCBI-SRA accession ID.
        known_proteins -- dict, containing viral family names as keys and a
                          list of NCBI protein identifiers belonging to the
                          family as value.

    Returns:
        genome_completeness -- str, specifies the determined genome
                               completeness based on the given cut-offs in
                               the config file.
        best_alignment -- lst, containing the alignment ID as the first
                          element and a list as second element. This list
                          contains all the standard DIAMOND output format 6
                          elements and an additional column containing the
                          product of the %id * alignment length.
        viral_family -- str, name of the viral family to which the target
                        ID belongs to.

    This function is used to check the viral completeness of the given
    diamond file. It checks the completeness of the best alignment using the
    viral specific parameters given in the config file.
    """
    # Determine the best alignment within the given diamond file
    best_alignment = fetch_largest_product(diamond_dict)
    best_target = best_alignment[1][0]
    viral_family = get_key(known_proteins, str(best_target))

    # Check whether the viral family is unclassified
    if viral_family == "unclassified":
        return "unclassified", best_alignment, viral_family

    # Import the cut-offs for every viral family from the config file
    min_assembly_len = config_dict[str(viral_family)]["min_genome_len"]
    max_assembly_len = config_dict[viral_family]["max_genome_len"]
    assembly_len = float(best_alignment[0].split("&")[0].split("_")[3])

    # Check the completeness of the alignment
    # Complete genome
    if min_assembly_len < assembly_len < max_assembly_len:
        return "complete", best_alignment, viral_family

    # Incomplete genomes
    elif assembly_len < min_assembly_len:
        if count_target_occurrence(diamond_dict, best_alignment[1][0]) > 1:
            # Multiple alignments, maybe eligible for patching
            return "fragmented", best_alignment, viral_family
        elif count_target_occurrence(diamond_dict, best_alignment[1][0]) <= 1:
            # Single alignment, incomplete genome
            return "incomplete", best_alignment, viral_family
    elif assembly_len > max_assembly_len:
        # The genome is longer than the maximum genome length of the family
        return "oversized", best_alignment, viral_family
    else:
        raise ValueError("Something went wrong while determining the viral "
                         "completeness of '{}'".format(accession))


def execute_contig_patch(contig_name, genome_template, contig_sequence,
                         alignment_info, log_file):
    """Patches a protein contig into a template protein genome sequence

    Key arguments:
        contig_name -- str, name of the contig that should be patched into
                       the template sequence.
        genome_template -- lst, consisting of "X" characters at every position.
                           The length of the list should correspond to the
                           length of the template genome.
        contig_sequence -- lst, containing a protein sequence in the form of
                           aminoacid characters at every position.
                           e.g. ['M', 'A', 'A', 'I'].
        alignment_info -- lst, containing the standard output format 6
                          alignment parameters of the diamond alignment
                          (except for the query id and bitscore).
                          I.e. [target_id, identity, length, mismatch,
                          gapopen, qstart, qend, tstart, tend, evalue].
        log_file --

    Returns:
        genome_template -- lst, input list where all 'X' characters are
                           substituted by valid amino acid characters where
                           the sequence position overlaps with the alignment
                           information of the given contig.
        contig_used -- bool, true when the given contig is 'patched' into
                       the genome_template, false if the contig is ignored.
        gap_len -- int, the gap between the given contig and the template
                   sequence, based on the diamond alignment information.
    """
    # Determine the number of gap positions that should be inserted
    gap_len = determine_gap_length(alignment_info)
    if gap_len < 0:
        c_print("Note: gap length is negative, so the target takes up more "
                "space than there is available in the range of the "
                "template!\n", None, log_file)
    elif gap_len > 100:
        c_print("NOTE: there is a huge gap between the target and "
                "contig sequence {}, of '{}' positions!".format(contig_name,
                                                                gap_len),
                None, log_file)

    # Determine the part of the query that aligns against the template
    query_start = int(alignment_info[5]) - 1
    query_end = int(alignment_info[6])
    contig_aligned_seq = contig_sequence[query_start:query_end]

    # Performing the patching
    template_start = int(alignment_info[7]) - 1
    template_end = int(alignment_info[8])

    # Print a message of the used parameters to the log file
    c_print("The following parameters were used to try patch contig '{}':\n"
            .format(contig_name) +
            "Length of aligned contig: '{}'\n"
            .format(len(contig_aligned_seq)) +
            "query_start\tquery_end\ttemplate_start\ttemplate_end: \n"
            "{:<16}{:<16}{:<16}{}".format(query_start, query_end,
                                            template_start, template_end),
            None, log_file)

    # Loop over the template and contig sequence and fill the template if X
    contig_used = False
    previous_pos_template = None
    for idx_template, contig_aa, idx_contig in \
            itertools.zip_longest(range(template_start, template_end),
                                  contig_aligned_seq,
                                  range(0, len(contig_aligned_seq)),
                                  fillvalue="-"):
        # If the template position is not filled, then there is room to patch
        if idx_template != "-":
            if genome_template[idx_template] == "X":
                genome_template[idx_template] = contig_aa
                previous_pos_template = idx_template
                contig_used = True
            else:
                # The template is already filled at this position
                continue
        else:
            # There is a gap in the template, so append to the last position
            if not previous_pos_template:
                # The first aa of the current contig should be appended to
                # the last position of the genome template
                previous_pos_template = len(genome_template) - 1
            c_print("Current aa: '{}' at position: '{}' in the contig, "
                    "is appended to the last 'patched' position in the "
                    "template: '{}'.".format(contig_aa, idx_contig,
                                             previous_pos_template),
                    None, log_file)
            genome_template[previous_pos_template] = \
                genome_template[previous_pos_template] + contig_aa

    # Print the layout of the patched genome to the log file
    if contig_used:
        c_print("The current contig was used during patching.\n"
                "The genome template with the current contig patched into "
                "it: \n {}\n".format(genome_template),
                None, log_file)
    else:
        c_print("The contig above was not used to patch the genome "
                "template!\n", None, log_file)

    return genome_template, contig_used, gap_len


def create_patched_genome(sorted_queries, template_list, input_seq_file,
                          log_file):
    """Used to fill a template genome sequence with given (protein) contigs

    Key arguments:
        sorted_queries -- OrderedDict, keys are alignment query IDs,
                          value is a dictionary containing 'dmnd_params' as
                          key and a list of the dmnd parameters as value.
                          E.g. {'NODE_1': {'dmnd_params': ['AHX00961.1',
                          '98.3', etc..]}}. The dict is ordered based on
                          alignment length (3rd element of the dmnd_params).
        template_list -- lst, filled with 'X' characters. The number of
                         elements should be equal to the length of the
                         template genome.
        input_seq_file -- str, absolute pathway to the file containing the
                          fasta format sequences of the nodes given in
                          'sorted_queries'.
        log_file -- str, absolute location of the file that should be used
                    to write log messages to.

    Returns:
        patched_sequence -- str, complete sequence of the patched contigs.
                            Also containing the remaining 'X' characters
                            and gaps ('-').
        used_contigs -- lst, containing the names of the contigs that were
                        used during the patching process.
    """
    used_contigs = []
    for contig in sorted_queries:
        # Obtain the protein sequences of the contigs into a list
        original_header = "_".join(contig.split("_")[:2]) + "_"
        contig_sequence = grep_prot_seq(original_header, input_seq_file)

        # Obtain the alignment information in a list
        alignment_info = sorted_queries[contig]["dmnd_params"]

        # Perform the patching of the sequence
        template_list, store_contig, curr_gap = \
            execute_contig_patch(contig, template_list, contig_sequence,
                                 alignment_info, log_file)

        # Store the name of the contigs that were used for patching
        if store_contig:
            used_contigs.append(contig)

    # Create a string of the sequence
    patched_sequence = "".join(template_list)

    return patched_sequence, used_contigs


def write_to_fasta(sequence, accession, best_target, contig_names, genome_nr,
                   output_file):
    """Will write a given string to the given output file in FASTA format

    Key arguments:
        sequence -- str, the string that should be converted to fasta format.
        accession -- str, the identifier of the accession that should be
                     used in the FASTA header.
        best_target -- str, the identifier of the target to which the
                       sequence has the best alignment. This is included in
                       the FASTA header.
        contig_names -- lst, full names of the nodes of which the sequence
                        consists. Shortened versions of these names are
                        included in the FASTA header.
        genome_nr -- int, the count of the genome, this will be included in
                     the FASTA header.
        output_file -- str, absolute pathway to the file where the FASTA
                       sequence should be writen to.

    Returns:
        None

    Note: this functions writes the sequences to the provided output file in
    FASTA format. The header will look like this: e.g. '>SRRXXXX_pg_nr1
    AET36829.1 [NODE_1_length_4000,NODE_2_length_1000]'. Where pg stands for
    'patched genome' the number is the count of the genomes found for this
    file, then the accession of the best target and lastly a list containing
    the names of all node IDs that were used during patching.
    """
    # Set and shorten the contig names
    identifier = "{}_pg_nr{}".format(accession, genome_nr)

    abbreviated_contig_names = ""
    for contig in contig_names:
        abbreviated_contig_names += "_".join(contig.split("_")[:4]) + ","
    description = best_target + " [" + abbreviated_contig_names[:-1] + "]"

    # Create a BioSeq record of the sequence and write to a file in FASTA form
    record = SeqRecord(Seq(sequence),
                       id=identifier,
                       description=description)

    with open(output_file, "a+") as fasta_out:
        SeqIO.write(record, fasta_out, "fasta")

    return


# Wrapper functions ###########################################################
def viral_completeness_wrapper(diamond_dict, config_dict,
                               accession, viral_filter):
    """Main function that wraps al functions to determine viral completeness

    Key arguments:
        diamond_dict -- dict, containing the query name as key and the other
                        columns of the DIAMOND file in a list, as value.
        config_dict -- dict, containing all parameters found in the
                       yml config file as a python dict object.
        accession -- str, the original NCBI-SRA accession ID.

    Optional arguments:
        viral_filter -- lst, viral family names that should be included in the
                        analysis. All other targets that are found will
                        be excluded. Default = iflaviridae.

    Returns:
        None -- This function prints the alignment info to the stdout.
     """
    # Determine the known proteins for every viral family based on configfile
    known_proteins = {"iflaviridae": line_parser(
                          config_dict["iflaviridae"]["known_proteins"]),
                      "dicistroviridae": line_parser(
                          config_dict["dicistroviridae"]["known_proteins"]),
                      "marnaviridae": line_parser(
                          config_dict["marnaviridae"]["known_proteins"]),
                      "picornaviridae": line_parser(
                          config_dict["picornaviridae"]["known_proteins"]),
                      "secoviridae": line_parser(
                          config_dict["secoviridae"]["known_proteins"]),
                      "calciviridae": ["QKW94212.1"]
                      }

    # Check whether the alignments should be filtered based on %ID * align_len
    if not viral_filter:
        viral_filter = "iflaviridae"
        diamond_dict = \
            filter_alignment_product(diamond_dict,
                                     config_dict["min_p_id_alnlen"])[0]
    else:
        viral_filter.append("unclassified")

    # Check if there are any alignments
    if len(fetch_unique_nodes(diamond_dict)) < 1:
        print_viral_completeness(accession, 0, "no_filtered_alignments",
                                 "NA", "NA" + "\tNA"*11)
        return

    # Loop over the diamond file to determine the viral completeness
    previous_genome_completeness = []
    genome_counter = 0
    while len(fetch_unique_nodes(diamond_dict)) > 0:
        completeness, best_alignment, viral_family = \
            check_viral_completeness(diamond_dict, config_dict, accession,
                                     known_proteins)

        # Obtain the target that was used and format it for printing
        best_target = best_alignment[1][0]
        best_alignment_print = "{}".format(best_alignment[0].split("&")[0])
        for item in best_alignment[1][:11]:
            best_alignment_print += "\t" + str(item)

        # Remove all alignments against the bests target from the dmnd file
        diamond_dict = remove_target_alignments(diamond_dict,
                                                best_target)

        # Check whether the viral genome should be included
        if viral_family not in viral_filter:
            # print("Best target is against an excluded viral family")
            continue

        if completeness not in previous_genome_completeness:
            previous_genome_completeness.append(completeness)
            genome_counter += 1

            # Print the viral completeness of the alignments to the screen
            print_viral_completeness(accession, genome_counter, completeness
                                     + "_1", viral_family,
                                     best_alignment_print)
        else:
            count_occurrence = \
                str(previous_genome_completeness.count(completeness) + 1)
            genome_counter += 1

            # Print the viral completeness of the alignments to the screen
            print_viral_completeness(accession, genome_counter, completeness
                                     + "_" + count_occurrence, viral_family,
                                     best_alignment_print)
            previous_genome_completeness.append(completeness)

        continue

    return


def patching_wrapper(args, diamond_dict, config_dict, config_file_loc,
                     accession, metadata, families_oi, printing, reset=None):
    """Main function to perform the patching of fragmented genome assemblies

    Key arguments:
        args -- Namespace, argparse object containing all the arguments
                given on the commandline.
        diamond_dict -- dict, containing the query name as key and the other
                        columns of the DIAMOND file in a list, as value.
        config_dict -- dict, containing all parameters found in the
                       yml config file as a python dict object.
        config_file_loc -- str, the location to the config file that was used.
        accession -- str, the original NCBI-SRA accession ID.
        metadata -- dict, accession ID as key and a dict as value, internal
                    dict should have keys for all the required data of this
                    accession and values should be the absolute locations of
                    these files.
        families_oi -- list, additional viral families that should be
                       included. Family names should be in str and separated
                       by spaces and commas. E.g. ['iflaviridae',
                       'dicistroviridae']. Default = None (only iflaviridae
                       is included).
        printing -- bool, if true, additional message will be printed to
                    stdout when the script is running.
        reset -- bool, True when the nucleotide and protein sequences should
                 be re-created even when the output files are already
                 present. False to use the files and sequences that are
                 already present. Default = None (re-creating).

    Returns:
        None
    """
    # Set logfile location, and create it
    log_file = config_dict["log_files"] + "/" + accession + "_patching.log"
    initiate_patch_log(args, log_file, accession, config_dict, config_file_loc)

    # Get the known proteins of the included viral families
    # Determine the known proteins for every viral family based on configfile
    known_proteins = {"iflaviridae": line_parser(
        config_dict["iflaviridae"]["known_proteins"]),
        "dicistroviridae": line_parser(
            config_dict["dicistroviridae"]["known_proteins"]),
        "marnaviridae": line_parser(
            config_dict["marnaviridae"]["known_proteins"]),
        "picornaviridae": line_parser(
            config_dict["picornaviridae"]["known_proteins"]),
        "secoviridae": line_parser(
            config_dict["secoviridae"]["known_proteins"]),
        "calciviridae": ["QKW94212.1"]
    }

    # Check the included viral families and filter accordingly
    if not families_oi:
        families_oi = "['iflaviridae']"

    if str(families_oi) == "['iflaviridae']":
        # Filter the diamond alignments on minimum %ID * alignment_length
        diamond_dict = \
            filter_alignment_product(diamond_dict,
                                     config_dict["min_p_id_alnlen"])[0]

    # Filter based on minimum contig assembly length
    diamond_dict = \
        filter_length(diamond_dict,
                      config_dict["best_target"]["min_contig_len"])[0]

    # Check if there are any alignments at all
    if len(fetch_unique_nodes(diamond_dict)) < 1:
        c_print("Skip this accession: '{}', no diamond alignment after "
                "filtering, so patching cannot be performed!"
                .format(accession),
                "WARNING", log_file)
        return

    # Loop diamond file to determine if there are multiple fragmented genomes
    nr_patched_genomes = 0
    nr_complete_genomes = 0
    nr_incomplete_genomes = 0
    nr_excluded_genomes = 0
    nr_unclassified_genomes = 0
    nr_oversized_genomes = 0
    genome_counter = 1
    while len(fetch_unique_nodes(diamond_dict)) > 0:
        # Check the viral completeness present in the given dmnd alignment file
        genome_completeness, best_alignment, viral_family = \
            check_viral_completeness(diamond_dict, config_dict, accession,
                                     known_proteins)
        best_target = best_alignment[1][0]

        # Skip current best target if it doesn't belong to included families
        if viral_family == "unclassified":
            c_print("The genome against target: '{}' is not being patched as "
                    "this genome belongs to a viral family that cannot be "
                    "classified according to the known protein files given "
                    "as input.\n".format(best_target),
                    "WARNING", log_file)
            nr_unclassified_genomes += 1
            genome_counter += 1
            # Remove all alignments against current best target from the dmnd
            diamond_dict = remove_target_alignments(diamond_dict, best_target)
            continue

        if viral_family not in families_oi:
            c_print("The genome against target: '{}' is not being patched as "
                    "this genome belongs to a viral family that is not "
                    "included in the commandline.\n".format(best_target),
                    "WARNING", log_file)
            nr_excluded_genomes += 1
            genome_counter += 1
            # Remove all alignments against current best target from the dmnd
            diamond_dict = remove_target_alignments(diamond_dict, best_target)
            continue

        # Print message to stdout based on cmd line args
        if printing:
            print("For accession: '{}', completeness of genome number '{}' "
                  "is: '{}'.".format(accession,
                                     genome_counter,
                                     genome_completeness))

        # Skip the current best target if the genome is not fragmented
        if genome_completeness != "fragmented":
            c_print("Skip the genome against this target: '{}', "
                    "as genome completeness is: '{}'"
                    .format(best_target, genome_completeness.upper()),
                    "OKGREEN", log_file)
            if genome_completeness == "complete":
                nr_complete_genomes += 1
            elif genome_completeness == "incomplete":
                nr_incomplete_genomes += 1
            elif genome_completeness == "oversized":
                nr_oversized_genomes += 1
            else:
                raise ValueError("ANOTHER VIRAL GENOME COMPLETENESS OPTION, "
                                 "PLEASE ACCOUNT FOR!")
            # Remove all alignments against current best target from the dmnd
            diamond_dict = remove_target_alignments(diamond_dict, best_target)
            continue

        # If viral genome is fragmented and against included family:
        # Obtain all the node names that have a hit against the best target
        queries = obtain_target_alignments(diamond_dict, best_target)
        c_print("-" * 70 + "\n" +
                "Based on the given diamond file fragmented viral genome nr "
                "'{}' was found. "
                "Now starting the patching!\n".format(genome_counter) +
                "Best target: '{}'\n"
                "Nodes against target: {}\n".format(best_target,
                                                    list(queries.keys())),
                None, log_file)

        # Remove all alignments against the current best target from the dmnd
        diamond_dict = remove_target_alignments(diamond_dict, best_target)

        # Obtain the nt sequences of all the query contigs and write to fasta
        output_nt_fasta = "{}/{}_best_nt_genome{}.fasta" \
            .format(config_dict["temp_files"], accession, genome_counter)

        # Check if the nucleotide sequences are already stored
        if os.path.isfile(output_nt_fasta) and not reset:
            c_print("NUCLEOTIDE SKIP (file with FASTA nucleotide sequences "
                    "already exists. Skipping the extraction from the "
                    "assembly file.)",
                    None, log_file)
        else:
            # Clear/create the nucleotide file and store the sequences
            with open(output_nt_fasta, "w"):
                pass
            grep_nt_sequence(queries,
                             metadata[accession]["assembly"],
                             output_nt_fasta)
            c_print("Creating this nucleotide fasta file: '{}' that "
                    "contains the nucleotide sequences of the contigs that "
                    "have an alignment against the best target."
                    .format(output_nt_fasta),
                    None, log_file)

        # Translate the nucleotide sequences into protein sequences
        output_prot_fasta = "{}/{}_prot_genome{}.fasta" \
            .format(config_dict["temp_files"], accession, genome_counter)
        if os.path.isfile(output_prot_fasta) and not reset:
            c_print("PROTEIN SKIP (file with FASTA protein sequences already "
                    "exists. Skipping the translation of the nt sequences.)",
                    None, log_file)
        else:
            with open(output_prot_fasta, "w"):
                pass
            reset = True
            # Obtain the protein sequence of the best target sequence
            obtain_orf(output_nt_fasta,
                       output_prot_fasta,
                       config_dict["best_target"]["min_prot_len_nt"],
                       log_file,
                       True)
            c_print("Creating this protein fasta file: '{}' that contains the "
                    "protein sequences of the contigs that have an alignment "
                    "against the best target: '{}'".format(output_nt_fasta,
                                                           best_target),
                    None, log_file)

        # Blastp using DIAMOND against the picornavirales db
        diamond_db = config_dict["diamond_prot_db"]
        diamond_cores = config_dict["diamond_cores"]
        diamond_outfile = "{}/{}_prot_genome{}.dmnd" \
            .format(config_dict["diamond_blastp_out"],
                    accession,
                    genome_counter)
        if os.path.isfile(diamond_outfile) and not reset:
            c_print("DIAMOND SKIP (file with protein alignments already "
                    "exists. Skipping the blastp of the protein sequences.)",
                    None, log_file)
        else:
            try:
                run_diamond_blastp(output_prot_fasta, diamond_db,
                                   diamond_outfile, diamond_cores)
                c_print("Creating this diamond file: '{}' that contains the "
                        "blastp output of the protein sequences of the "
                        "contigs that have an alignment against the best "
                        "target.\n".format(output_nt_fasta),
                        None, log_file)
            except subprocess.CalledProcessError:
                # Would be nice to include a feature here so script continues
                raise ValueError("This error shouldn't happen! This is likely "
                                 "due to the fact that the file with protein "
                                 "sequences is empty and thus DIAMOND cannot "
                                 "perform alignments. Please consider "
                                 "lowering the minimum protein length in the "
                                 "config file.") from None

        # Read the output of the DIAMOND blastp alignment
        blastp_dmnd_file = parse_diamond_files(diamond_outfile, "")

        # Double check if the best target is the same
        prot_best_target = fetch_largest_product(blastp_dmnd_file)[1][0]
        if best_target != prot_best_target:
            raise ValueError("The nt best target is not the same as the "
                             "protein best target for this accession: '{}'!"
                             .format(accession))

        # Determine the alignments/queries against the best target
        prot_queries = obtain_target_alignments(blastp_dmnd_file, best_target)

        # Sort the queries based on alignment length
        sorted_queries = \
            OrderedDict(sorted(prot_queries.items(),
                               key=lambda x: int(x[1]["dmnd_params"][2]),
                               reverse=True))
        c_print("The protein queries sorted on descending alignment "
                "length. In this order, they will be patched: ",
                None, log_file)
        for query in sorted_queries:
            c_print("{}\t{}".format(query,
                                    sorted_queries[query]["dmnd_params"]),
                    None, log_file)

        # Determine the template genome length and initiate list with Xs
        target_polyp_len = fetch_prot_len(best_target)
        template_list = initiate_dummy_lst(target_polyp_len)
        c_print("The genome template of '{}' which is '{}' amino acids long "
                "will be used.\n".format(best_target, target_polyp_len),
                None, log_file)

        # Perform the patching, using the template and contig protein sequence
        nr_patched_genomes += 1
        c_print("Starting to fill the template sequence with contig "
                "characters.",
                None, log_file)

        patched_sequence, used_contigs = \
            create_patched_genome(sorted_queries, template_list,
                                  output_prot_fasta, log_file)

        # Write sequence to the output in FASTA format and write log message
        fasta_output_file = config_dict["output_fasta_dir"] + "/" + \
                            accession + "_patched_genomes.fasta"
        write_to_fasta(patched_sequence, accession, best_target, used_contigs,
                       genome_counter, fasta_output_file)
        c_print("The patching of the genome with the best alignment against: "
                "{} is done! Now checking for other eligible genomes.\n"
                .format(best_target) +
                "-" * 70 + "\n\n",
                None, log_file)

        # Restart at the top of the while loop to check for other genomes
        genome_counter += 1
        continue

    # Print message to logfile and stdout
    c_print("\n" + "\u2304" * 70 + "\n"
            "\t\t\tFINISHED\n" +
            "Patching is done, no more alignments left in this "
            "diamond file!\n" +
            "STATISTICS FOR THIS ACCESSION:\n"
            "Number of patched genomes: {}\n"
            .format(nr_patched_genomes) +
            "Number of complete genomes: {}\n"
            .format(nr_complete_genomes) +
            "Number of incomplete genomes: {}\n"
            .format(nr_incomplete_genomes) +
            "Number of excluded viral genomes: {}\n"
            .format(nr_excluded_genomes) +
            "Number of unclassified viral genomes: {}\n"
            .format(nr_unclassified_genomes) +
            "Number of oversize viral genomes: {}\n"
            .format(nr_oversized_genomes) +
            "\u2304" * 70,
            None, log_file)

    return


def main():
    """Main function that wraps all the code in this script """
    # Step 0.0: parsing the cmd line arguments and config file
    parser = parsing_cmd()
    args = parser.parse_args()
    config = config_parsing(args.config)
    print(Bcolors.OKGREEN +
          "The script is started correctly, please wait while it is "
          "running!" +
          Bcolors.ENDC)

    # Step 1: obtain the accessions from the input file
    accessions = line_parser(config["input_accessions"])
    if not accessions:
        print(Bcolors.FAIL +
              "No accessions could be obtained from the input file. Please "
              "double check that the given input file containing accessions "
              "is not empty!" +
              Bcolors.ENDC)
        return

    # Step 2: determine all file locations for every accession
    metadata_dict = {}
    for run in accessions:
        metadata_dict[run] = \
            {"metadata": "{}/{}{}".format(config["metadata_loc"],
                                          run,
                                          config["metadata_ext"]),
             "logfile": "{}/{}_patch.log".format(config["log_files"],
                                                 run),
             "diamond": "{}/{}{}".format(config["diamond_loc"],
                                         run,
                                         config["diamond_ext"]),
             "assembly": "{}/{}/{}".format(config["assembly_loc"],
                                           run,
                                           config["assembly_ext"])}

    # Step 2.1: determine if the directory for log, temp and output files exists
    if not os.path.exists(config["log_files"]):
        print(Bcolors.OKGREEN +
              "Directory to store logfiles: '{}' doesn't exits. Creating it "
              "now!".format(config["log_files"]) +
              Bcolors.ENDC)
        subprocess.check_output("mkdir {}".format(config["log_files"]),
                                shell=True)
    if not os.path.exists(config["temp_files"]):
        print(Bcolors.OKGREEN +
              "Directory to store temp files: '{}' doesn't exits. Creating it "
              "now!".format(config["temp_files"]) +
              Bcolors.ENDC)
        subprocess.check_output("mkdir {}".format(config["temp_files"]),
                                shell=True)
    if not os.path.exists(config["output_fasta_dir"]):
        print(Bcolors.OKGREEN +
              "Directory to store output files: '{}' doesn't exits. Creating it "
              "now!".format(config["log_files"]) +
              Bcolors.ENDC)
        subprocess.check_output("mkdir {}".format(config["log_files"]),
                                shell=True)

    # Step 3: append the taxid to the metadata dict of every accession
    metadata_pop = []
    for acc in metadata_dict:
        if args.safe:
            try:
                metadata_dict[acc]["taxid"] = \
                    obtain_taxid(metadata_dict[acc]["metadata"])
            except ValueError:
                print("No metadata for {}".format(acc))
                metadata_pop.append(acc)
                print(metadata_pop)
        else:
            metadata_dict[acc]["taxid"] = \
                obtain_taxid(metadata_dict[acc]["metadata"])

    # Step 3.2: remove accessions for which no metadata is available
    for key in metadata_pop:
        metadata_dict.pop(key)

    # Print message for stdout, if viral completeness should be determined
    if "viral_completeness" in args.main_parser:
        if not args.included_fam:
            print(Bcolors.OKGREEN +
                  "The alignments are filtered for a minimum product "
                  "of 'ID% * alignment length' > {}!"
                  .format(config["min_p_id_alnlen"]),
                  Bcolors.ENDC)
        else:
            print(Bcolors.FAIL +
                  "Viral families include other families than "
                  "Iflaviridae, so no filter for the product of the "
                  "'%ID * alignment_length' could be applied!" +
                  Bcolors.ENDC)

    # Print message to stdout when patching
    if "patching" in args.main_parser:
        print(Bcolors.OKGREEN +
              "Now starting the patching process, please wait!" +
              Bcolors.ENDC)
        if not args.included_fam or \
                str(args.included_fam) == "['iflaviridae']":
            print("The alignments are filtered for a minimum "
                  "product of 'ID% * alignment length' > {}!"
                  .format(config["min_p_id_alnlen"]))

            # Print message about settings to stdout
            print("Minimum contig length filter: {}"
                  .format(config["best_target"]["min_contig_len"]))
            print("Minimum 'EMBOSS getorf' ORF size: {}\n\n"
                  .format(config["best_target"]["min_prot_len_nt"]))
        else:
            print(Bcolors.FAIL +
                  "NOTE: The included viral families contain other "
                  "families than Iflaviridae, so no filter for the product "
                  "of the '%ID * alignment_length' could be applied!" +
                  Bcolors.ENDC)

    # Step 4: loop through diamond files and check completeness or patch
    counter_no_file = 0
    for accession in metadata_dict:
        try:
            diamond_dict = \
                parse_diamond_files(metadata_dict[accession]['diamond'],
                                    metadata_dict[accession]['taxid'])
        except FileNotFoundError:
            # If the safe option is given, script will ignore accessions
            # without diamond file
            if args.safe:
                counter_no_file += 1
                continue
            else:
                raise FileNotFoundError("Accession: '{}' has no "
                                        "dmnd file!".format(accession))

        # Check subcommands, either viral_completeness or patching
        # Check the viral completeness per accession based on the family
        if "viral_completeness" in args.main_parser:
            viral_completeness_wrapper(diamond_dict, config, accession,
                                       args.included_fam)
        # Patch the fragmented genomes
        elif "patching" in args.main_parser:
            patching_wrapper(args, diamond_dict, config, args.config,
                             accession, metadata_dict, args.included_fam,
                             args.print, args.reset)

    # Print message of patching completion to stdout
    if "patching" in args.main_parser:
        print(Bcolors.OKGREEN +
              "Patching is finished!\n"
              "The following config file was used: '{}'"
              .format(args.config) +
              Bcolors.ENDC)


if __name__ == "__main__":
    main()
