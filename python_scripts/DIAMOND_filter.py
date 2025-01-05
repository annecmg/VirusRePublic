#!/usr/bin/env python3
"""Author: Devin van Valkengoed

Description: This script can be used to filter the top hits of a DIAMOND
             alignment using custom parameters given by the user.
             The parameters that can be changed are:
                * Coverage
                * Alignment %ID
                * E-value
                * SPAdes original contig length

Usage: DIAMOND_filter.py [-h] (-d [[dmnd] ...] | -D [input-dir]) (-o <output> | -p) [-q]
       [-id [%AID]] [-l [<min_len> [[max_len] ...]]] [-al <alignment_len>] [-c [coverage]] [-e [E-value]]
       [--nodes | --top | --targets | --top_aid | --self] [--version]

Dependencies
    * Biopython

**Enhancements**:
    * Split the main function in general filtering and advanced filtering.
    * Make it so the output also shows the host species of the accession?
      Combination of efetching the metadata and then running meta_accession\
      parser.py functions to determine the host?
"""
# Import statements
import argparse
from translate_rna import *

# Version and date
__version__ = "2.0"
__date__ = "December 2024"


# Commandline parsing:
def parsing_cmd():
    """This function uses argparse to obtain the given command line arguments

    Returns:
        parser -- argparse.Namespace object
    """
    # Initiate argument parser
    parser = argparse.ArgumentParser(description="This script takes one or "
                                                 "multiple DIAMOND alignment "
                                                 "files in BLAST tabular "
                                                 "format. It will then filter "
                                                 "the alignments files, based "
                                                 "on the parameter settings "
                                                 "given by the user and "
                                                 "output all the alignments "
                                                 "in one file. \n"
                                                 f"Version: {__version__}, {__date__}",
                                     formatter_class=argparse.RawTextHelpFormatter)

    # Provide diamond files with (multiple) string pathways or one complete directory
    in_files = parser.add_mutually_exclusive_group(required=True)
    in_files.add_argument("-d", "--diamond", help="Absolute pathway to one or "
                                                  "multiple DIAMOND alignment "
                                                  "files in tabular output "
                                                  "format, given as a string separated by space.",
                          metavar="[dmnd]", nargs="*",
                          dest="dmnd")
    in_files.add_argument("-D", "--input-dir",
                          help="Absolute pathway to an input directory that contains diamond output files with a "
                               "'.dmnd' extension that should be used as input.",
                          metavar="[input-dir]", dest="dir")

    # Add general arguments
    output_args = parser.add_mutually_exclusive_group(required=True)
    output_args.add_argument("-o", "--output", help="Pathway and name to use for "
                                                    "the output file, given as a string. NOTE: if the given output "
                                                    "already exists, the data will be appended.",
                             metavar="<output>", dest="out")
    output_args.add_argument("-p", "--print", help="Print the output to the stdout.",
                             action="store_true", dest="print")
    parser.add_argument("-q", "--quiet", help="Minimal print statements to "
                                              "the screen will be shown. This will also exclude the header information "
                                              "from being printed to stdout.",
                        action="store_true", dest="quiet")

    # Filtering options
    parser.add_argument("-id", "--percid", help="Floating point number "
                                                "between 0.0 and 100.0 used "
                                                "as the minimum percent "
                                                "alignment identity (%%AID) "
                                                "between the query and the "
                                                "subject. Every alignment > "
                                                "this given %%AID is saved "
                                                "to the output file.",
                        metavar="[%AID]", dest="AID")
    parser.add_argument("-l", "--length",
                        help="Integer, filter for alignments with an original "
                             "contig length equal to or longer than the "
                             "given number. Second argument is optional and "
                             "can be used to filter for a max contig length. \n"
                             "Only alignments with original contig length "
                             "equal to or below the max limit will be stored. "
                             "Original contig length is given by SPAdes in "
                             "number of nucleotides.",
                        metavar=("<min_len>", "[max_len]"), dest="length",
                        nargs="*")
    parser.add_argument("-al", "--alignment_length",
                        help="Integer, filter for alignments with a minimum "
                             "length equal to or longer than the given "
                             "number.",
                        metavar="<alignment_len>", dest="al")
    parser.add_argument("-c", "--coverage", help="Floating point number used "
                                                 "as the minimum K-mer "
                                                 "coverage. All alignments "
                                                 "with coverage > "
                                                 "'--coverage' will be saved "
                                                 "to the output file.",
                        metavar="[coverage]", dest="cov")
    parser.add_argument("-e", "--evalue", help="Floating point number used as "
                                               "the maximum E-value. Every "
                                               "alignment that has a E-value "
                                               "> 'eval' will be ignored.",
                        metavar="[E-value]", dest="eval")

    # Specific filtering options
    advanced_filtering = parser.add_mutually_exclusive_group()
    advanced_filtering.add_argument("--top_aid",
                                    help="Obtains the best alignment for every unique node, based on alignment "
                                         "identity percentage (%% AID).",
                                    action="store_true", dest="top_aid")
    advanced_filtering.add_argument("--top", help="Obtain the contig with the best "
                                                  "hit, based on %%AID for every "
                                                  "given input file. Other commands "
                                                  "can be used to filter for e.g. "
                                                  "alignment length or contig length. ",
                                    action="store_true", dest="top")
    advanced_filtering.add_argument("--nodes", help="When this option is given, "
                                                    "only the unique nodes that are "
                                                    "found after filtering are "
                                                    "given as output. The output format is as: "
                                                    "<file_name> tab <node_name> tab <file_path> ",
                                    action="store_true", dest="nodes")
    advanced_filtering.add_argument("--targets",
                                    help="This option shows all the accessions of the "
                                         "unique targets to which there is an alignment "
                                         "after filtering.",
                                    action="store_true", dest="target")
    advanced_filtering.add_argument("--self", help="When this option is given, self "
                                                   "alignments are removed from the "
                                                   "output based on query and "
                                                   "target name. I.e. if the query "
                                                   "and target protein are the same "
                                                   "and thus give 100%% AID.",
                                    action="store_true", dest="filter_self")

    # Add version
    parser.add_argument("--version", action="version", version=f"{__version__}")

    return parser


def write_line(line: str | list, output_file: str) -> None:
    """
    Writes the given string or list of strings to a new line in the specified file.

    :param line: The content to be written. Can be a string or a list of strings.
    :type line: str or list
    :param output_file: The path to the file where the content should be appended.
    :type output_file: str
    :returns: None
    :rtype: NoneType
    """
    with open(output_file, 'a+') as output:
        if isinstance(line, str):
            output.write(str(line) + "\n")
        elif isinstance(line, list):
            for single_line in line:
                output.write(str(single_line) + "\n")

    return


def parse_location_file(input_file: str) -> list:
    """
    Reads the lines of the given input file containing diamond locations.

    :param input_file: Path to a file with absolute locations of diamond files,
                       each on a new line.
    :type input_file: str
    :returns: A list of unique locations of the diamond files to be used during filtering.
    :rtype: list
    """
    file_list = []
    with open(input_file, 'r') as infile:
        for line in infile:
            if not line:
                continue
            line = line.strip()
            if line not in file_list:
                file_list.append(line)

    return file_list


def parse_location_dir(input_dir: str, ext: str = ".dmnd") -> list:
    """
    Obtains all files in a given directory that have the desired extension.

    :param input_dir: Path to the directory from which files should be extracted.
    :type input_dir: str
    :param ext: Extension of the files to be obtained. Defaults to ".dmnd".
    :type ext: str
    :returns: A list of absolute file paths for files in the directory with the given extension.
    :rtype: list
    """
    # Make sure the given extension starts with a dot
    if not ext.startswith("."):
        ext = "." + ext

    # Find all files with the given extension
    file_list = [
        os.path.abspath(os.path.join(input_dir, file))
        for file in os.listdir(input_dir)
        if file.endswith(ext)
    ]

    return file_list


def line_parser(input_file: str) -> list:
    """
    Parses each new line of a file as an element in a list.

    :param input_file: Path to a .txt file containing information on each new line.
    :type input_file: str
    :returns: A list containing all lines of the file as elements.
    :rtype: list
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


def read_metadata(meta_data: str) -> list:
    """
    Reads the provided metadata file into a nested list.

    :param meta_data: Path to a .txt file providing metadata of SRA experiments
                      in the standard NCBI format.
    :type meta_data: str
    :returns: A nested list containing two lists. The first list holds the
              headers of the NCBI-SRA metadata file, and the second list contains the values.
    :rtype: list
    """
    file_lines = []

    # Save the two lines of a NCBI-SRA metadata file into a list
    with open(meta_data, 'r') as in_file:
        for line in in_file:
            file_lines.append(line.strip().split(","))

    return file_lines


def get_key(input_dict: dict[str, list], element: str) -> str:
    """
    Obtains the key corresponding to a given value present in the dictionary's value lists.

    :param input_dict: A dictionary with keys mapped to lists of elements.
                       These lists are searched for the given value.
    :type input_dict: dict[str, list]
    :param element: The value to be found within the dictionary's value lists.
    :type element: str
    :returns: The corresponding key if the value is found, or "unclassified" if not.
    :rtype: str
    """
    for key, val in input_dict.items():
        if element in val:
            return key

    return "unclassified"


def obtain_taxid(meta_data_list: list[list[str]]) -> str:
    """
    Obtains the taxonomic identifier (TaxID) from an NCBI-SRA accession metadata file.

    :param meta_data_list: A nested list containing two lists:
                           - The first list contains the headers of the NCBI-SRA metadata file.
                           - The second list contains the values.
    :type meta_data_list: list[list[str]]
    :returns: The taxonomic identifier (TaxID) as determined by NCBI-SRA.
    :rtype: str
    :raises ValueError: If the TaxID is not at the expected position, or if no accession number is provided.
    """
    # Check if the current TaxID is at the correct position and not empty
    if meta_data_list[0][27] == "TaxID":
        if meta_data_list[1][0] != '':
            tax_id = meta_data_list[1][27]
        else:
            raise ValueError(
                "Exception while parsing metadata. No accession number is provided "
                "in the metadata file. Please ensure the file is in NCBI-SRA format."
            )
    else:
        raise ValueError(
            "Exception while parsing metadata. TaxID is not at the expected position "
            "and cannot be determined. Please check the metadata file structure."
        )

    return tax_id


def parse_diamond_files(diamond, filter_self, taxid):
    """Reads DIAMOND files, in tabular format, into memory

    Key arguments:
        diamond -- str, pathway to a DIAMOND file in BLAST tabular
                   format (-f 6).
        filter_self -- bool, true if self alignments should be filtered. I.e.
                       proteins aligned against itself and thus 100% AID.
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
            if filter_self:
                if line[0] == line[1]:
                    continue
                else:
                    if taxid==0:
                        key_name = "{}&against&{}".format(line[0],line[1])
                    else:
                        key_name = "{}_{}&against&{}".format(line[0],
                                                         taxid,
                                                         line[1])
                    dmnd_file[key_name] = line[1:]
            else:
                if taxid==0:
                    key_name = "{}&against&{}".format(line[0],line[1])
                else:
                    key_name = "{}_{}&against&{}".format(line[0], taxid, line[1])
                dmnd_file[key_name] = line[1:]

    return dmnd_file


def filter_length(
        diamond_dict: dict[str, list],
        min_length: int,
        max_length: int = None
) -> tuple[dict[str, list], dict[str, list], dict[str, list]]:
    """
    Filters DIAMOND alignments based on the original contig length.

    :param diamond_dict: Dictionary containing the query name as key and other
                         columns of the DIAMOND file as a list, as value.
    :type diamond_dict: dict[str, list]
    :param min_length: Minimum length of the original contig (nt) to be included
                       in the output.
    :type min_length: int
    :param max_length: Maximum length of the original contig (nt) to be included
                       in the output. Defaults to None.
    :type max_length: int, optional
    :returns:
        - length_filtered_alignments: A dictionary containing query names as keys
          and the other columns of the DIAMOND file as values, for alignments that
          meet the length criteria.
        - below_min_length: A dictionary of alignments ignored for being below
          the minimum length.
        - above_max_length: A dictionary of alignments ignored for exceeding
          the maximum length.
    :rtype: tuple[dict[str, list], dict[str, list], dict[str, list]]
    """
    length_filtered_alignments = {}
    below_min_length = {}
    above_max_length = {}

    for header in diamond_dict:
        original_header = header.split("&")
        length = float(original_header[0].split("_")[3])  # Parse the contig length

        if not max_length:  # Only filtering for minimum length
            if length >= min_length:
                length_filtered_alignments[header] = diamond_dict[header]
            else:
                below_min_length[header] = diamond_dict[header]
        else:  # Both minimum and maximum length filtering
            if min_length <= length <= max_length:
                length_filtered_alignments[header] = diamond_dict[header]
            elif min_length > length:
                below_min_length[header] = diamond_dict[header]
            elif max_length < length:
                above_max_length[header] = diamond_dict[header]

    return length_filtered_alignments, below_min_length, above_max_length


def filter_alignment_length(
        diamond_dict: dict[str, list],
        min_align: int
) -> tuple[dict[str, list], int]:
    """
    Filters DIAMOND alignments based on the alignment length.

    :param diamond_dict: Dictionary containing the query name as key and other
                         columns of the DIAMOND file as values in a list.
    :type diamond_dict: dict[str, list]
    :param min_align: Minimum alignment length to be included in the output.
    :type min_align: int
    :returns: A tuple containing:
        - alignment_len_filtered_alignments: Dictionary with query names as keys
          and other columns of the DIAMOND file as values, for alignments that
          meet the length criteria.
        - alignment_len_ignored: The number of alignments ignored due to not
          meeting the minimum alignment length.
    :rtype: tuple[dict[str, list], int]
    """
    alignment_len_ignored = 0
    alignment_len_filtered_alignments = {}

    for header in diamond_dict:
        if int(diamond_dict[header][2]) >= min_align:
            alignment_len_filtered_alignments[header] = diamond_dict[header]
        else:
            alignment_len_ignored += 1

    return alignment_len_filtered_alignments, alignment_len_ignored


def filter_coverage(
        diamond_dict: dict[str, list],
        min_cov: float
) -> tuple[dict[str, list], int]:
    """
    Filters DIAMOND alignments based on a given minimum coverage.

    :param diamond_dict: Dictionary containing the query name as key and other
                         columns of the DIAMOND file as values in a list.
    :type diamond_dict: dict[str, list]
    :param min_cov: Minimum K-mer coverage for alignments to be included in the
                    filtered output.
    :type min_cov: float
    :returns: A tuple containing:
        - cov_filtered_alignments: A dictionary with query names as keys and
          other columns of the DIAMOND file as values, for alignments that meet
          the coverage criteria.
        - cov_ignored_alignments: The number of alignments ignored due to not
          meeting the minimum coverage.
    :rtype: tuple[dict[str, list], int]
    """
    cov_ignored_alignments = 0
    cov_filtered_alignments = {}

    for key in diamond_dict:
        original_key = key.split("&")
        coverage = float(original_key[0].split("_")[-2])

        if coverage >= min_cov:
            cov_filtered_alignments[key] = diamond_dict[key]
        else:
            cov_ignored_alignments += 1

    return cov_filtered_alignments, cov_ignored_alignments


def filter_alignment_id(
        diamond_dict: dict[str, list],
        min_aid: float
) -> tuple[dict[str, list], int]:
    """
    Filters DIAMOND alignments based on a given minimum percentage of alignment identity.

    :param diamond_dict: Dictionary containing the query name as key and other
                         columns of the DIAMOND file as values in a list.
    :type diamond_dict: dict[str, list]
    :param min_aid: Minimum percentage alignment identity between the query and
                    the subject to be included in the filtered output.
    :type min_aid: float
    :returns: A tuple containing:
        - id_filtered_alignments: A dictionary with query names as keys and other
          columns of the DIAMOND file as values, for alignments that meet the
          alignment identity criteria.
        - id_ignored_alignments: The number of alignments ignored due to not
          meeting the minimum alignment identity.
    :rtype: tuple[dict[str, list], int]
    """
    id_ignored_alignments = 0
    id_filtered_alignments = {}

    for key in diamond_dict:
        if float(diamond_dict[key][1]) >= min_aid:
            id_filtered_alignments[key] = diamond_dict[key]
        else:
            id_ignored_alignments += 1

    return id_filtered_alignments, id_ignored_alignments


def filter_e_value(
        diamond_dict: dict[str, list],
        min_evalue: float
) -> tuple[dict[str, list], int]:
    """
    Filters DIAMOND alignments based on a given minimum E-value.

    :param diamond_dict: Dictionary containing the query name as key and other
                         columns of the DIAMOND file as values in a list.
    :type diamond_dict: dict[str, list]
    :param min_evalue: Maximum E-value that alignments can have to be included
                        in the filtered output. This value can be provided in
                        scientific notation (e.g., 5.32e-56).
    :type min_evalue: float
    :returns: A tuple containing:
        - eval_filtered_alignments: A dictionary with query names as keys and
          other columns of the DIAMOND file as values, for alignments that meet
          the E-value criteria.
        - eval_ignored_alignments: The number of alignments ignored due to
          exceeding the E-value threshold.
    :rtype: tuple[dict[str, list], int]
    """
    eval_ignored_alignments = 0
    eval_filtered_alignments = {}

    for key in diamond_dict:
        if float(diamond_dict[key][9]) <= min_evalue:
            eval_filtered_alignments[key] = diamond_dict[key]
        else:
            eval_ignored_alignments += 1

    return eval_filtered_alignments, eval_ignored_alignments


def filter_alignment_product(
        diamond_dict: dict[str, list],
        min_product: float
) -> tuple[dict[str, list], int]:
    """
    Filters DIAMOND alignments based on the minimum product of alignment ID
    and alignment length.

    :param diamond_dict: Dictionary containing the query name as key and other
                         columns of the DIAMOND file as values in a list.
    :type diamond_dict: dict[str, list]
    :param min_product: Minimum product of the alignment ID (%) and alignment
                        length for alignments to pass the filter.
    :type min_product: float
    :returns: A tuple containing:
        - product_filtered_alignments: A dictionary with query names as keys and
          other columns of the DIAMOND file as values, for alignments that meet
          the product criteria.
        - nr_product_ignored: The number of alignments ignored due to not meeting
          the minimum product threshold.
    :rtype: tuple[dict[str, list], int]
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


def fetch_unique_nodes(
        diamond_dict: dict[str, list]
) -> dict[str, int]:
    """
    Finds the unique nodes in a given DIAMOND alignment dictionary.

    :param diamond_dict: Dictionary containing the query name as key and other
                         columns of the DIAMOND file as values in a list.
    :type diamond_dict: dict[str, list]
    :returns: A dictionary with unique nodes as keys and an integer value
              representing the number of alignments for each node.
              Example: {'NODE_2': 25} shows that the alignment file contains
              one unique node (NODE_2) with 25 different alignments.
    :rtype: dict[str, int]
    """
    unique_nodes = {}

    for key in diamond_dict:
        node = key.split("&")[0]
        if node not in unique_nodes:
            unique_nodes[node] = 1
        else:
            unique_nodes[node] += 1

    return unique_nodes


def obtain_unique_targets(
        diamond_dict: dict[str, list]
) -> list[str]:
    """
    Finds all the unique targets that have an alignment in a DIAMOND file.

    :param diamond_dict: Dictionary containing the query name as key and other
                         columns of the DIAMOND file as values in a list.
    :type diamond_dict: dict[str, list]
    :returns: A list containing the identifiers of the unique targets that
              have an alignment in the given DIAMOND file.
    :rtype: list[str]
    """
    unique_targets = []
    for alignment in diamond_dict:
        target = diamond_dict[alignment][0]
        if target not in unique_targets:
            unique_targets.append(target)
        else:
            continue

    return unique_targets


def find_best_alignment(
        diamond_dict: dict[str, list]
) -> tuple[str, list]:
    """
    Finds the best alignment in a DIAMOND file based on the percentage of identical positions.

    :param diamond_dict: Dictionary where keys represent query names, and values are lists containing
                         columns of the DIAMOND file.
    :type diamond_dict: dict[str, list]
    :returns:
        A tuple containing:
        - **top_aid_key** (*str*): The key of the alignment with the highest percentage
          of identical positions.
        - **top_aid_alignment** (*list*): A list containing all columns of the DIAMOND
          output for the best alignment.
    :rtype: tuple[str, list]
    """
    top_aid = 0.0  # Start with the lowest possible alignment identity percentage
    top_aid_key = ""
    top_aid_alignment = []

    for alignment, details in diamond_dict.items():
        current_aid = float(details[1])  # Percentage of identical positions
        if current_aid > top_aid:
            top_aid = current_aid
            top_aid_key = alignment
            top_aid_alignment = details

    return top_aid_key, top_aid_alignment


def print_general(
        file_name: str,
        diamond_dict: dict[str, list],
        nr_total_align: int,
        nr_total_nodes: int,
        ignored_num: int,
        quiet: bool
) -> None:
    """
    Prints a summary of the filtered alignments and detailed alignments to the screen.

    :param file_name: Name of the DIAMOND alignment file.
    :type file_name: str
    :param diamond_dict: Dictionary containing query names as keys and the other columns of
                         the DIAMOND file as values.
    :type diamond_dict: dict[str, list]
    :param nr_total_align: Total number of alignments in the original DIAMOND file.
    :type nr_total_align: int
    :param nr_total_nodes: Total number of unique query nodes in the original DIAMOND file.
    :type nr_total_nodes: int
    :param ignored_num: Number of alignments that did not pass the filters.
    :type ignored_num: int
    :param quiet: If `True`, suppresses printing the summary information.
    :type quiet: bool
    :returns: None
    :rtype: NoneType
    """
    header = (
        "Query accession \t Target accession \t Sequence identity \t "
        "Length \t Mismatches \t Gap openings \t Query start \t "
        "Query end \t Target start \t Target end \t E-value \t Bit score"
    )

    # Collect unique query names
    unique_queries = [
        "_".join(query.split("_")[:2])
        for query in fetch_unique_nodes(diamond_dict).keys()
    ]

    # Print summary information if quiet mode is disabled
    if not quiet:
        print(f"# Name of the alignment file: {file_name}")
        print(f"# Number of total unique queries: {nr_total_nodes}")
        print(
            f"# Number of unique queries that passed the filters: "
            f"{len(fetch_unique_nodes(diamond_dict))}"
        )
        print(f"# The unique query names: {unique_queries}")
        print(f"# Number of total alignments in original file: {nr_total_align}")
        print(f"# Number of alignments that did not pass the set filters: {ignored_num}")
        print("# " + header)

    # Print filtered alignments
    for key, values in diamond_dict.items():
        print(key.split("&")[0], *values, sep="\t")


def write_general_out(file_name: str, diamond_dict: dict[str, list]) -> None:
    """
    Writes the given DIAMOND alignments to a specified output file.

    :param file_name:Path to the output file where the alignments should be written.
    :type file_name: str
    :param diamond_dict: Dictionary containing query names as keys and the other columns
                         of the DIAMOND file as values.
    :type diamond_dict: dict[str, list]
    :returns:
        None
    :rtype: NoneType
    """
    with open(file_name, "a+") as outfile:
        for key, values in diamond_dict.items():
            # Create a tab-separated line
            line = "{}\t{}".format(key.split("&")[0], "\t".join(map(str, values)))
            outfile.write(line + "\n")


def get_top_alignments(diamond_dict: dict[str, list], accession: str) -> list[tuple]:
    """
    Finds the top alignments for each unique node in a DIAMOND file.

    :param diamond_dict: Dictionary containing the query name as key and the other columns
                         of the DIAMOND file as values.
    :type diamond_dict: dict[str, list]
    :param accession: Name of the given accession in the DIAMOND file.
    :type accession: str
    :returns: A list of tuples, where each tuple contains the accession, node name,
              and the corresponding best alignment values.
    :rtype: list[tuple]

    This function identifies the alignment with the highest percentage identity
    (%AID) and maximum alignment length for each unique node in the DIAMOND file.
    """
    # Fetch unique nodes from the diamond dictionary
    unique_nodes = fetch_unique_nodes(diamond_dict)

    # List to store results
    results = []

    # Iterate over each unique node and determine its best alignment
    for node in unique_nodes:
        max_aid_perc = 0
        max_align_length = 0
        max_aid_target = None

        # Check all alignments for the current node
        for aligned_key, aligned_value in diamond_dict.items():
            if node in aligned_key:
                curr_aid = float(aligned_value[1])  # % identity
                curr_align_length = float(aligned_value[2])  # alignment length

                # Update the best alignment based on %AID and alignment length
                if curr_aid > max_aid_perc and curr_align_length > max_align_length:
                    max_aid_perc = curr_aid
                    max_align_length = curr_align_length
                    max_aid_target = aligned_value

        # Print the best alignment for the current node
        if max_aid_target:
            results.append((accession, node, *max_aid_target))

    return results


def main():
    """Main function that wraps all the code in this script """
    # Step 0.0: check if the command line input is given
    parser = parsing_cmd()
    args = parser.parse_args()

    # Step 0.1: check whether all given subcommands are empty
    specific_options = [args.nodes, args.top, args.top_aid, args.target]

    blacklisted_accessions = []
    unique_targets = []
    host_count = []

    # Step 0.2 parse the location to the provided diamond files
    if args.dir:
        diamond_files = parse_location_dir(args.dir)
    elif args.dmnd:
        diamond_files = args.dmnd
    else:
        # Argparse should handle this case
        raise argparse.ArgumentError(message="No diamond file given!",
                                     argument=args.dir)

    for file in diamond_files:
        # Step 0.3: set the ignored alignments to 0
        ignored = 0

        accession_name = file.split("/")[-1].split(".")[0].split("_")[0]

        if accession_name in blacklisted_accessions:
            continue

        # Obtain the metadata location and taxid of the current accession
        #metadata_location = os.getcwd() + "/output/metadata/paired_accessions/" \
        #                                  "{}_metadata.txt".format(accession_name)
        #taxonomic_id = obtain_taxid(read_metadata(metadata_location))
        taxonomic_id=0

        # Step 1.1: read the diamond file(s) into memory
        # Double check if the file exists at the given location
        try:
            alignments = parse_diamond_files(file, args.filter_self,
                                             taxonomic_id)
        except FileNotFoundError:
            if args.quiet:
                continue
            elif not args.quiet:
                print("File: '{}' was not found, continuing!".format(file))
                continue

        nr_total_alignments = len(alignments.keys())
        nr_total_unique_nodes = len(fetch_unique_nodes(alignments))

        # Step 1.2: check if the given alignment files contain spades
        # assembled contigs. Otherwise, the coverage filter will fail
        for alignment_name in alignments.keys():
            if "cov" not in alignment_name and any([args.cov, args.length]):
                raise UserWarning("The given file does not contain contigs "
                                  "assembled by SPAdes (K-mer coverage is not "
                                  "in the contig name). The coverage or length"
                                  " filter can thus not be used. \n"
                                  "Please run the script again without the "
                                  "coverage/length filter or provide a "
                                  "different DIAMOND alignment file.")
            else:
                continue

        # Step 2.1: filter based on original contig length
        if args.length:
            min_length = int(args.length[0])
            if len(args.length) == 2:
                max_length = float(args.length[1])
            else:
                max_length = None
            length_filtered = filter_length(alignments, min_length, max_length)
            alignments = length_filtered[0]
            ignored += len(length_filtered[1])

        # Step 2.1.2: filter for alignment length
        if args.al:
            al_filtered = filter_alignment_length(alignments, int(args.al))
            alignments = al_filtered[0]
            ignored += al_filtered[1]

        # Step 2.2: filter based on coverage
        if args.cov:
            cov_filtered = filter_coverage(alignments, float(args.cov))
            alignments = cov_filtered[0]
            ignored += cov_filtered[1]

        # Step 2.3: filter based on % alignment identity
        if args.AID:
            aid_filtered = filter_alignment_id(alignments, float(args.AID))
            alignments = aid_filtered[0]
            ignored += aid_filtered[1]

        # Step 2.4: filter based on E-value
        if args.eval:
            eval_filtered = filter_e_value(alignments, float(args.eval))
            alignments = eval_filtered[0]
            ignored += eval_filtered[1]

        # Step 3.1: obtain the unique nodes having an alignment against the db
        unique_nodes_filtered = fetch_unique_nodes(alignments)
        nr_unique_nodes_filtered = len(unique_nodes_filtered.keys())

        # Step 3.2: obtain the unique targets
        for target in obtain_unique_targets(alignments):
            if target not in unique_targets:
                unique_targets.append(target)

        # NOTE: all steps below are optional based on the given cmd arguments
        # If no other filtering options are given, the filtered output is printed or writen to the output file
        # Step 4.1a: print the filtered alignments
        if args.print and all(v is False for v in specific_options):
            print_general(file, alignments, nr_total_alignments,
                          nr_total_unique_nodes, ignored, args.quiet)

        # Step 4.1b: write the filtered alignments
        elif args.out and all(v is False for v in specific_options):
            write_general_out(file_name=args.out, diamond_dict=alignments)

        # Step 4.2: print the best alignment, based on %AID
        elif args.top:
            if len(obtain_unique_targets(alignments)) > 0:
                top_alignment = find_best_alignment(alignments)
                top_stats = [top_alignment[0].split("&")[0]] + top_alignment[1]

                # Either print the output or write to file
                if args.print:
                    print(accession_name, *top_stats, sep="\t")
                else:
                    write_line("\t".join([accession_name] + list(map(str, top_stats))), args.out)
            else:
                host_count.append("no_alignments")

        # Step 4.3: unique nodes for every file that pass the filter
        elif args.nodes:
            for node in unique_nodes_filtered.keys():
                acc = file.split("/")[-1]
                if args.print:
                    print(acc + "\t" + node + "\t" + file + "\t")
                else:
                    if args.out:
                        write_line("{} \t {} \t {}".format(acc, node, file),
                                   args.out)
                    else:
                        raise argparse.ArgumentError(message="No output file "
                                                             "was given!",
                                                     argument=args.out)

        # Step 4.4: filter for the top alignment %AID for every node
        elif args.top_aid:
            top_alignments = get_top_alignments(alignments, accession_name)
            for alignment in top_alignments:
                if args.print:
                    print("\t".join(map(str, alignment)))
                else:
                    write_line("\t".join(map(str, alignment)), args.out)

    # Step 4.5: output all unique targets
    if args.target:
        for target in unique_targets:
            if args.print:
                if not args.quiet:
                    print(f"{accession_name}\t{target}")
                else:
                    print(target)
            else:
                if not args.quiet:
                    write_line(f"{accession_name}\t{target}", args.out)
                else:
                    write_line(target, args.out)


if __name__ == "__main__":
    main()
