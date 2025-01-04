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
from datetime import datetime
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

        with open(params.log, 'a+') as logfile:
            with redirect_stdout(logfile):
                # Obtain the current time
                now = datetime.now()
                current_time = now.strftime("%d-%b-%Y %H:%M:%S")

                # Setting and printing the command
                cmd_fastp = "fastp --in1 {} --in2 {} --out1 {} --out2 {} " \
                            "--dedup --dup_calc_accuracy 3 -j {} -h {}"\
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


# Mapping of trimmed reads
rule mapping:
    """HISAT2 Maps the trimmed reads against the corresponding host genome """
    input:
        fwd=config["root_dir"] + config["trimmed_dir"] +
            "{accession}_trimmed_1.fq.gz",
        rev=config["root_dir"] + config["trimmed_dir"] +
            "{accession}_trimmed_2.fq.gz"
    output:
        temp(config["root_dir"] + config["mapped_dir"] + "{accession}.sam")
    priority: 40
    params:
        log = config["logs"]["root"] + config["logs"]["main"] +
              "main_{accession}.log"
    resources:
        disk_mb = 1
    threads: config["mapping"]["threads"] # default = 4
    run:
        with open(params.log,'a') as logfile:
            with redirect_stdout(logfile):
            # Determine the right host genome to map to
                host_taxid = determine_host(wildcards.accession)[0]

                # Map reads of certain hosts against close relatives
                if int(host_taxid) == 213944:
                    close_relative = 110799
                    print(Bcolors.WARNING +
                          "The host genome of '{}' is not available, mapping "
                          "against the assembly of close "
                          "relative '{}'.".format(host_taxid, close_relative) +
                          Bcolors.ENDC, file=logfile)

                    host_index_path = "{}{}/{}_index"\
                        .format(config["root_dir"] +
                                config["host_genomes"],
                                close_relative,
                                close_relative)
                else:
                    # Set the path to the host genome
                    host_index_path = "{}{}/{}_index".\
                        format(config["root_dir"] +
                               config["host_genomes"],
                               host_taxid,
                               host_taxid)

                    # Print messages to the log file
                    print(Bcolors.OKGREEN +
                          "The current accession is: '{}'\n"
                          .format(wildcards.accession),
                        "The reads will be mapped against the host genome "
                        "with TaxID: '{}'\n".format(host_taxid),
                          "Pathway of this host genome: '{}'\n"
                          .format(host_index_path) +
                          Bcolors.ENDC,
                        file=logfile)

                # Running the shell command
                cmd_mapping = "hisat2 -t -p {} -x {} -1 {} -2 {} -S {} " \
                              "--no-temp-splicesite".format(threads,
                                                            host_index_path,
                                                            input.fwd,
                                                            input.rev,
                                                            output)
                print(Bcolors.OKGREEN +
                      "Shell command to perform the mapping: '{}'"
                      .format(cmd_mapping) +
                      Bcolors.ENDC, file=logfile)
                print(subprocess.check_output(cmd_mapping,
                                              stderr=subprocess.STDOUT,
                                              shell=True).decode("UTF-8"),
                      file=logfile)


# Convert the sam files to .bam files
rule sam2bam:
    """Converts the mapped sam files into bam files to save disk space"""
    input:
        config["root_dir"] + config["mapped_dir"] + "{accession}.sam"
    output:
        temp(config["root_dir"] + config["mapped_dir"] +
        "{accession, [A-Za-z0-9]+}.bam")
    priority: 45
    params:
        log=config["logs"]["root"] + config["logs"]["main"] +
            "main_{accession}.log"
    threads:
        config["mapping"]["threads"] # default = 4
    message:
        "Converting the .sam file of the mapped reads of accession: "
        "'{wildcards.accession}' into a .bam file using samtools view."
    run:
        with open(params.log,'a') as logfile:
            with redirect_stdout(logfile):
                cmd_convert = "(samtools view --threads {} -Sbh {} > " \
                              "{}) >> {} 2>&1".format(threads, input, output,
                                                      params.log)

                print(Bcolors.OKGREEN +
                      "Converting the .sam file of accession: '{}' "
                      "into a .bam file using the following command: '{}'.\n"
                      "For the error log, "
                      "see: '{}'\n".format(wildcards.accession,
                                           cmd_convert,
                                           params.log) +
                      Bcolors.ENDC,
                      file=logfile)

                print(subprocess.check_output(cmd_convert,
                    stderr=subprocess.STDOUT,
                    shell=True).decode("UTF-8"),
                    file=logfile)


########################### Read extraction ###################################
# Extracting unmapped reads
rule bam_extraction:
    """Extracts the reads that are not mapped against the reference genome """
    input:
        config["root_dir"] + config["mapped_dir"] +
        "{accession}.bam"
    output:
        temp(config["root_dir"] + config["mapped_dir"] +
        "{accession, [A-Za-z0-9]+}_unmapped.bam")
    priority: 50
    params:
        log = config["logs"]["root"] + config["logs"]["main"] +
              "main_{accession}.log"
    threads:
        config["mapping"]["threads"] # default = 4
    run:
        with open(params.log, 'a') as logfile:
            with redirect_stdout(logfile):
                cmd_extraction = "samtools view --threads {} -h -b -f 4 {} " \
                                 "-o {}".format(threads, input, output)

                print(Bcolors.OKGREEN +
                      "Extracting the unmapped reads of accession '{}' "
                      "using the following command: "
                      "'{}'\n".format(wildcards.accession,
                                      cmd_extraction) +
                      Bcolors.ENDC, file=logfile)

                print(subprocess.check_output(cmd_extraction,
                                              stderr=subprocess.STDOUT,
                                              shell=True).decode("UTF-8"),
                      file=logfile)


# Convert from sam to fastq file
rule convert_bam_to_fastq:
    """Converts the unmapped reads from .bam to FastQ.gz format """
    input:
        config["root_dir"] + config[
            "mapped_dir"] + "{accession}_unmapped.bam"
    output:
        extract_fwd=config["root_dir"] + config[
            "extracted_dir"] + "{accession}_ext_1P.fq.gz",
        extract_rev=config["root_dir"] + config[
            "extracted_dir"] + "{accession}_ext_2P.fq.gz",
        extract_unpaired=config["root_dir"] + config[
            "extracted_dir"] + "{accession}_ext_unpaired.fq.gz"
    priority: 55
    params:
        log=config["logs"]["root"] + config["logs"]["main"] +
            "main_{accession}.log"
    threads:
        config["mapping"]["threads"] # default = 4
    message:
        "Extracting reads from accession: '{wildcards.accession}'"
    run:
        with open(params.log, 'a') as logfile:
            with redirect_stdout(logfile):
                cmd_conversion = "samtools fastq --threads {} -1 {} -2 {} " \
                                 "-s {} -N {}"\
                                 .format(threads,
                                         output.extract_fwd,
                                         output.extract_rev,
                                         output.extract_unpaired,
                                         input)

                print(Bcolors.OKGREEN +
                      "Converting the .bam file of accession: '{}' "
                      "into a FastQ.gz file using the following command: "
                      "'{}'\n".format(wildcards.accession,
                                      cmd_conversion) +
                      Bcolors.ENDC, file=logfile)

                print(subprocess.check_output(cmd_conversion,
                                              stderr=subprocess.STDOUT,
                                              shell=True).decode("UTF-8"),
                      file=logfile)
