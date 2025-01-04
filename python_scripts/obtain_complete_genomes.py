#!/usr/bin/env python3
"""Author: Devin van Valkengoed

Date: 14-Aug-2023
Description: This script is used to obtain a desired contig/genome from an
             original rnaviralSPAdes assembly file containing multiple
             sequences. This genome is than writen to a separat FASTA file
             and also translated to a protein sequences and stored in
             another FASTA file.

Usage: python3 obtain_complete_genome.py -input_accessions <input_tsv>
       -nt_output <nuc_out> -prot_output <prot_out>

Key arguments:
    intput_tsv -- str, absolute pathway to a tsv file containing information
                  about the desired sequences. The columns should be as
                  followed: accession_id \t genome_nr \t NODE_name \t
                  best_target.
    nuc_output -- str, absolute pathway to the folder where all output
                  nucleotide FASTA files should be written to.
    prot_output -- str, absolute pathway to the folder where all output
                   protein FASTA files should be written to.
"""
# Import statements
import argparse
import os
import sys
import subprocess
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO


# Commandline argument parsing
def parsing_cmd():
    """Uses argparse to determine the given commandline options

    Returns:
        parser -- argparse.Namespace object
    """
    parser = argparse.ArgumentParser(
        description="This script is used to obtain nucleotide sequence of "
                    "interest from rnaviralSPAdes assembly files. It will "
                    "also translate the nucleotide sequences into protein "
                    "sequences.")
    parser.add_argument("-input_accessions",
                        help="Str -- Absolute pathway to a text file "
                             "containing the file directories, "
                             "fasta header and additional description for "
                             "every accession of interest. The file should "
                             "be tab separated and every new line should "
                             "look like this: <folder/> \t "
                             "<genome_nr> \t <NODE_X> \t "
                             "<best_target>",
                        metavar="input_accessions.txt", dest="input",
                        required=True)
    parser.add_argument("-nt_output",
                        help="Str -- Absolute pathway to the directory where "
                             "all FASTA files containing the complete "
                             "nucleotide sequence should be stored.",
                        metavar="./nt_folder/", dest="nt_output",
                        required=True)
    parser.add_argument("-prot_output",
                        help="Str, -- Absolute pathway to the directory "
                             "where all FASTA files containing the complete "
                             "protein sequence should be stored.",
                        metavar="./prot_folder/", dest="prot_output",
                        required=True)

    return parser


# Functions
def parse_input_file(input_file):
    """Parses the input file to obtain the input folder, node and description

    Key arguments:
        input_file -- str, absolute pathway to an input file containing all
                      information of the files from which sequences should
                      be obtained. The file format should be tab separated
                      and should look like this:
                      <absolute_path/to_folder/> \t
                      <genome_nr> \t
                      <NODE_of_interest> \t
                      <additional_description>

    Returns:
        sequence_locations -- lst, of tuples containing the following
                              elements for every accession: accession_id,
                              file_location, genome_nr, node, description.
    """
    sequence_locations = []
    with open(input_file, "r") as in_file:
        for line in in_file:
            if not line or line.startswith("#"):
                continue
            else:
                accession, genome_nr, node, description = line.strip().split()
                location = "/lustre/BIF/nobackup/valke024/msc_thesis/" \
                           "main_data/assembly/rnaviralspades/" + accession + \
                           "/contigs.fasta"

                sequence_locations.append((accession, location, genome_nr,
                                           node, description))

    return sequence_locations


def grep_nt_sequence(accession, nodes, fasta_file, output_file, description):
    """Obtains the original nt sequences from contigs.fasta files using sed

    Key arguments:
        accession -- str, NCBI-SRA accession ID to use in the FASTA header.
        nodes -- lst, contig header names of which the nt sequence should be
                 obtained.
        fasta_file -- str, pathway to a fasta file which contains the
                      contigs that are provided in the nodes list.
        output_file -- str, output file to which the fasta sequences should
                       be writen to.
        description -- lst, genome number and the best target. These are
                       included in the output FASTA header.

    Returns:
        None
    """
    for alignment in nodes:
        original_node_name = "_".join(alignment.split("_")[:-1])
        tax_id = alignment.split("_")[-1]
        sed_cmd = "sed -n '/{}/,/>/p' {} | " \
                  "head -n -1 | tail +2"\
            .format(original_node_name,
                    fasta_file)

        # Use sed to obtain the sequence of every contig
        sequence = subprocess.check_output(sed_cmd,
                                           shell=True).decode("UTF-8").strip()

        # Replace newline characters that sed includes
        sequence = sequence.replace("\n", "")

        # Store the sequence information in a BioPython SeqRecord object
        record = SeqRecord(
            Seq(sequence),
            id=accession + "_" + description[0],
            description=description[1] + " " + tax_id
        )

        # Write the sequences to the same fasta file
        with open(output_file, "a+") as fasta_out:
            SeqIO.write(record, fasta_out, "fasta")

    return


def main():
    """Main function that wraps all the code in this python script"""
    # Step 0: parse the cmd line arguments
    parser = parsing_cmd()
    args = parser.parse_args()

    # Step 1: obtain the location of all assembly using the input file
    locations = parse_input_file(args.input)

    # Step 2: obtain the nt sequences from the assembly files
    nucleotide_files = []
    all_targets = {}
    for accession_info in locations:
        # Obtain all the file information
        accession_id = accession_info[0]
        assembly_file = accession_info[1]
        genome_nr = accession_info[2]
        node = accession_info[3]
        best_target = accession_info[4]
        description = [genome_nr, best_target]
        all_targets[accession_id + "_" + genome_nr] = best_target

        # Grep the nucleotide sequence of the contig from the assembly file
        output_nt = args.nt_output + accession_id + \
                    "_nt_complete_genomes.fasta"
        grep_nt_sequence(accession_id, [node], assembly_file, output_nt,
                         description)
        nucleotide_files.append(output_nt)

    # Step 3: translate the nucleotide sequences into protein sequences
    unique_nuc_files = set(nucleotide_files)
    protein_files = []

    for nuc_file in unique_nuc_files:
        accession_identifier = nuc_file.split("/")[-1].split("_")[0]
        print(accession_identifier)
        prot_out_file = args.prot_output + accession_identifier +\
                        "_prot_complete_genomes.fasta.del"

        cmd_get_orf = "getorf -sequence {} -minsize 6000 -auto -outseq {}"\
            .format(nuc_file, prot_out_file)
        subprocess.check_output(cmd_get_orf, shell=True)

        protein_files.append(prot_out_file)

    # Step 4: change the headers of the protein FASTA sequences
    for prot_file in protein_files:
        final_file = prot_file[:-4]
        with open(prot_file, "r") as input_file:
            for record in SeqIO.parse(input_file, "fasta"):
                if int(record.id.split("_")[-1]) != 1:
                    raise ValueError("This nucleotide translation yielded "
                                     "multiple proteins. Accession: {}"\
                                     .format(final_file))
                else:
                    record.id = "_".join(record.id.split("_")[:-1])
                    prot_target = all_targets[record.id]
                    record.description = record.description.split(" ")[-1]
                    record.description = prot_target + " " + record.description
                    with open(final_file, "a+") as fasta_out:
                        SeqIO.write(record, fasta_out, "fasta")

    # Step 5: remove the intermediate protein files with extension .del
    for p_file in protein_files:
        os.remove(p_file)


if __name__ == "__main__":
    main()
