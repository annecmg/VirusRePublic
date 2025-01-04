"""
Author: Devin van Valkengoed

Description: Please note that this is a helper module for Snakefile_main.smk
             which is the main workflow of my thesis. This Snakemake file is
             used to process the reads of the given accessions. It will trim
             the reads, map against the host genome and extract the unmapped
             reads.
Date: 21-Dec-2022
This file includes the following rules:
    * fastp -- quality and adapter trims the reads using FastP.
    * mapping -- maps the trimmed reads against the corresponding host genome
    * sam2bam -- converts the mapped .sam files into .bam files to save space
    * bam_extraction -- extracts the unmapped reads, using samtools
    * bam_sort -- sorts the extracted reads
    * convert_bam_to_fastq -- converts the bam files into a FastQ file that can
                              be used for assembly by SPAdes.

NOTE!! The file that contains the provided NCBI-SRA accession numbers given in
the configfile, can only contain accessions that have a host genome with a
RefSeq genome assembly available. (This can be checked using
meta_accession_parser.py script)
"""
######################### Import statements ###################################
import sys
import os.path
import subprocess
from snakemake.utils import min_version
from contextlib import redirect_stdout
sys.path.insert(1, "/" + config["scripts"]["python_root"])
from meta_accession_parser import Bcolors
from meta_accession_parser import read_metadata
from meta_accession_parser import obtain_taxid


############################ Functions ########################################
def read_lines(input_file):
    """Returns a list of all words on every new line in a given file

    Key arguments:
        input_file -- str, pathway to a file that has accessions on every new
                      line.

    Returns:
        samples -- lst, list containing an element for every new line.
    """
    with open(input_file) as f:
        samples = [sample for sample in f.read().split() if len(sample) > 0]
    f.close()

    return samples


def determine_host(accessions):
    """Determines the host of an accession number based on the metadata

    Key arguments:
        accessions -- lst/str, either a list of accessions or a single
                      accession given as a string.

    Returns:
          taxid -- lstr/str, either a list of TaxID's or a single TaxID
                   corresponding to the given accession.


    Please note that this function only works when it is used in this snakefile
    as it also needs a config file and the metadata of the given accessions in
    the exact location specified by the config file.
    """
    taxid = []
    if type(accessions) == list:
        for acc in accessions:
            taxid.append(obtain_taxid("{}_metadata.txt"
            .format(str(config["root_dir"] +
                        config["metadata"]["meta_root"] +
                        config["metadata"]["paired_accessions"] +
                        acc))))

    elif type(accessions) == str:
        taxid.append(obtain_taxid("{}_metadata.txt"
             .format(str(config["root_dir"] +
                         config["metadata"]["meta_root"] +
                         config["metadata"]["paired_accessions"] +
                         accessions))))

    return taxid


######################## Processing of reads ##################################
# Trimming of reads using fastp
rule fastp:
    """Uses FastP to trim the reads of the provided accessions

    TODO:
        * Make the input .gz files?
    """
    input:
        fwd = config["root_dir"] +
              config["raw_reads"] +
              "{accession}_1.fastq",
        rev = config["root_dir"] +
              config["raw_reads"] +
              "{accession}_2.fastq"
    output:
        fwd_trimmed=temp(config["root_dir"] +
                         config["trimmed_dir"] +
                         "{accession}_trimmed_1.fq.gz"),
        rev_trimmed=temp(config["root_dir"] +
                         config["trimmed_dir"] +
                         "{accession}_trimmed_2.fq.gz")
    priority: 15
    params:
        log = config["logs"]["root"] +
              config["logs"]["main"] +
              "main_{accession}.log",
        html_file = config["logs"]["root"] +
                    config["logs"]["trimming"] +
                    "{accession}_fastp.html",
        json_file = config["logs"]["root"] +
                    config["logs"]["trimming"] +
                    "{accession}_fastp.json"
    resources:
        disk_mb = 2
    threads: 3 # Fastp default
    run:
        # Check for or create the logfile
        os.makedirs(os.path.dirname(params.log), exist_ok=True)

        with open(params.log, 'a') as logfile:
            with redirect_stdout(logfile):
                # Obtain the current time
                now = datetime.now()
                current_time = now.strftime("%d-%b-%Y %H:%M:%S")

                # Setting and printing the command
                cmd_fastp = "fastp --in1 {} --in2 {} --out1 {} --out2 {} " \
                            "--dedup --dup_calc_accuracy 6 -j {} -h {}"\
                            .format(input.fwd, input.rev,
                                    output.fwd_trimmed, output.rev_trimmed,
                                    params.json_file,
                                    params.html_file)

                print(Bcolors.OKGREEN +
                      "Current time: {}\n"
                      "Trimming the forward and reverse read of "
                      "accession: '{}' using the following command: '{}'"
                      .format(current_time, wildcards.accession, cmd_fastp) +
                      Bcolors.ENDC, file=logfile)

                print(subprocess.check_output(cmd_fastp,
                                              stderr=subprocess.STDOUT,
                                              shell=True).decode("UTF-8"),
                      file=logfile)
