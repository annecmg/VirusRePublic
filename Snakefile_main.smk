"""
Author: Devin van Valkengoed
Description: Snakemake file to extract possible viral contigs from sequencing
             data.
Date: 13-Oct-2022
      13-March-2024 V2.0
      Sep-2024 V3.0
Usage: snakemake -s <Snakefile_main.smk> -j 50 -n
       --scheduler greedy --resources api_calls=9 disk_mb=15

This file includes the following rules:
    * all -- Used to download the genome assemblies and produce the final
             output (viral genome assembly and DIAMOND alignment)
    * assembly -- Creates an assembly of the unmapped paired reads using either
                  rnaviralSPAdes or metaSPAdes.
    * diamond -- Aligns all the assembled contigs against a given protein
                 database using DIAMOND2.

Key arguments:
    Snakefile_main.smk -- str, absolute pathway to this Snakefile.

Please also download the config.yaml file, update the parameters in the
    config file according to your needs and update the pointer to the config
    file in the script below.

"""
######################### Config file #########################################
# Location of provided config file
configfile: "config/main_config.yaml"


######################### Import statements ###################################
import sys
import os.path
import subprocess
from datetime import datetime
from snakemake.utils import min_version
sys.path.insert(1, workflow.basedir + "/" + config["scripts"]["python_root"])
from read_lines import read_lines
from meta_accession_parser import Bcolors
from meta_accession_parser import determine_kmer
from meta_accession_parser import read_metadata
from meta_accession_parser import obtain_taxid


########################### Version check #####################################
min_version("5.2.0")


######################### Include statements ##################################
# make sure the location of these scripts is relative to this main script!
include: config["scripts"]["downloading_snakefile"] # download SRA files
include: config["scripts"]["host_genome_snakefile"] # download host genomes
include: config["scripts"]["read_processing_snakefile"] # processing of reads


############################ Functions ########################################
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
            if acc.startswith("#"):
                continue
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


######################## Global variables #####################################
def print_config_info():
    """Print config info to stdout when starting Snakemake script. """
    print("Number of used total accessions: {}".format(len(
        read_lines(config["root_dir"] +
                   config["metadata"]["meta_root"] +
                   config["metadata"]["stats"] +
                   config["metadata"]["filtered_accessions"]))))
    print(f"Root directory set to:{config['root_dir']}\n")

    # Log the used conda env
    conda_env = subprocess.check_output("echo $CONDA_DEFAULT_ENV",
        shell=True).decode("UTF-8")
    print(f"Conda environment that was used running the workflow: {conda_env}")

    # Log the sys.path
    print(f"Sys PATH locations that were used: {sys.path}")
    print(read_lines(config["root_dir"] +
                     config["metadata"]["meta_root"] +
                     config["metadata"]["stats"] +
                     config["metadata"]["filtered_accessions"]))

onstart:
    print_config_info()


######################### Config file checks  #################################
# Checking the assembly tool specified in the provided config file
if config["assembly"]["tool"].lower() not in ["rnaviralspades", "metaspades"]:
    raise Exception("Incorrect assembly tool given: {}, please check "
                      "the assembly tool specified in the config file. "
                      "Use either rnaviralSPAdes or metaSPAdes. Stopping the "
                      "Snakemake process!")
    sys.exit(1)


#################### Defining output for every rule ###########################
# The input that is commented out is mostly used to create intermediate files
# during testing. Please note that the final output (assembly and
# DIAMOND alignment) and the indexing of the host genomes should never be
# commented out when using this pipeline.
rule all:
    input:
        # # Download of NCBI-SRA files (downloading Snakefile)
        # expand(config["root_dir"] + config["download"]["dir"] +
        #        "{accession}/{accession}.sra",
        #        accession=read_lines(config["root_dir"] +
        #                       config["metadata"]["meta_root"] +
        #                       config["metadata"]["stats"] +
        #                       config["metadata"]["filtered_accessions"])),
        # # Splitting of the .sra files into a forward and reverse read
        # expand(config["root_dir"] + config["raw_reads"] +
        #        "{accession}_1.fastq",
        #        accession=read_lines(config["root_dir"] +
        #                       config["metadata"]["meta_root"] +
        #                       config["metadata"]["stats"] +
        #                       config["metadata"]["filtered_accessions"])),
        # # Download of host genome file
        # expand(config["root_dir"] + config["host_genomes"] +
        #        "{taxid}/{taxid}.zip",
        #        taxid=determine_host(read_lines(config["root_dir"] +
        #                       config["metadata"]["meta_root"] +
        #                       config["metadata"]["stats"] +
        #                       config["metadata"]["filtered_accessions"]))),
        # # Unzipping of the RefSeq genome
        # expand(config["root_dir"] + config["host_genomes"] +
        #        "{taxid}/ncbi_dataset/",
        #        taxid=determine_host(read_lines(config["root_dir"] +
        #                       config["metadata"]["meta_root"] +
        #                       config["metadata"]["stats"] +
        #                       config["metadata"]["filtered_accessions"]))),
        # # Rehydrating the genome assemblies
        # expand(config["root_dir"] + config["host_genomes"] +
        #        "{taxid}/{taxid}_genome.fna",
        #        taxid=determine_host(read_lines(config["root_dir"] +
        #                       config["metadata"]["meta_root"] +
        #                       config["metadata"]["stats"] +
        #                       config["metadata"]["filtered_accessions"]))),
        # MAIN Host genome download and indexing -- HAST TO STAY IN RULE ALL!!!
        expand(config["root_dir"] + config["host_genomes"] +
               "{taxid}/{taxid}_index.1.ht2",
               taxid=determine_host(
                   read_lines(config["root_dir"] +
                              config["metadata"]["meta_root"] +
                              config["metadata"]["stats"] +
                              config["metadata"]["filtered_accessions"])
               )
        ),
        # # FastP trimming (read processing Snakefile)
        # expand(config["root_dir"] + config["trimmed_dir"] +
        #        "{accession}_trimmed_{read}.fq.gz",
        #        accession=read_lines(config["root_dir"] +
        #                       config["metadata"]["meta_root"] +
        #                       config["metadata"]["stats"] +
        #                       config["metadata"]["filtered_accessions"]),
        #        read=[1, 2]),
        # # Hisat2 mapping (read processing Snakefile)
        # expand(config["root_dir"] + config["mapped_dir"] +
        #        "{accession}.sam",
        #        accession=read_lines(config["root_dir"] +
        #                       config["metadata"]["meta_root"] +
        #                       config["metadata"]["stats"] +
        #                       config["metadata"]["filtered_accessions"])),
        # # Convert mapped files from .sam to .bam (read processing Snakefile)
        # expand(config["root_dir"] + config["mapped_dir"] +
        #        "{accession}.bam",
        #        accession=read_lines(config["root_dir"] +
        #                       config["metadata"]["meta_root"] +
        #                       config["metadata"]["stats"] +
        #                       config["metadata"]["filtered_accessions"]))
        # # Extract unmapped reads (read processing Snakefile)
        # expand(config["root_dir"] + config["mapped_dir"] +
        #        "{accession}_unmapped.bam",
        #        accession=read_lines(config["root_dir"] +
        #                       config["metadata"]["meta_root"] +
        #                       config["metadata"]["stats"] +
        #                       config["metadata"]["filtered_accessions"])),
        # # Converting from sam to fastq (read processing Snakefile)
        # expand(config["root_dir"] + config["extracted_dir"] +
        #        "{accession}_ext_1P.fq.gz",
        #        accession=read_lines(config["root_dir"] +
        #                       config["metadata"]["meta_root"] +
        #                       config["metadata"]["stats"] +
        #                       config["metadata"]["filtered_accessions"])),
        # rnaviralSPAdes/metaSPAdes assembly
        expand(config["root_dir"] + config["assembly_dir"] +
               "{tool}/{accession}/contigs.fasta",
               tool=config["assembly"]["tool"].lower(),
               accession=read_lines(config["root_dir"] +
                              config["metadata"]["meta_root"] +
                              config["metadata"]["stats"] +
                              config["metadata"]["filtered_accessions"])),
        # DIAMOND2
        expand(config["root_dir"] + config["diamond_dir"] +
               "{accession}_" + config["diamond"]["name"] + ".dmnd",
               accession=read_lines(config["root_dir"] +
                              config["metadata"]["meta_root"] +
                              config["metadata"]["stats"] +
                              config["metadata"]["filtered_accessions"]))


##################### Assembly using SPAdes/Trinity ###########################
rule assembly:
    """Creates an assembly of the unmapped paired reads using SPAdes """
    input:
        fwd = config["root_dir"] +
              config["extracted_dir"] +
              "{accession}_ext_1P.fq.gz",
        rev = config["root_dir"] +
              config["extracted_dir"] +
              "{accession}_ext_2P.fq.gz"
    output:
        config["root_dir"] +
        config["assembly_dir"] +
        "{tool}/{accession}/contigs.fasta"
    priority: 60
    threads: config["assembly"]["threads"] # default = 16
    params:
        outdir = config["root_dir"] +
                 config["assembly_dir"] +
                 "{tool}/{accession}/",
        log = config["logs"]["root"] +
              config["logs"]["main"] +
              "main_{accession}.log",
        kmer_values = config["assembly"]["kmer_values"],
    run:
        with open(params.log, 'a') as logfile:
            # Obtain the current time
            now = datetime.now()
            current_time = now.strftime("%d-%b-%Y %H:%M:%S")

            # Determine the str of K-mers that should be used for SPAdes
            kmers = determine_kmer(input.fwd, params.kmer_values)

            # Check the assembly tool that should be used
            if config["assembly"]["tool"].lower() == "rnaviralspades":
                tool = "--rnaviral"
            else:
                tool = "--meta"

            # Run the assembly using rnaviral- or metaSPAdes.
            cmd_assembly = "spades.py {} -1 {} -2 {} -t {} -k {} -o {}"\
                           .format(tool, input.fwd, input.rev, threads,
                                   kmers, params.outdir)
            # Print the used command to the log file
            print(Bcolors.OKGREEN +
                  "Performing assembly with {}SPAdes using the "
                  "following command: '{}'\n"
                  "Current time: {}\n".format(tool[2:],
                                              cmd_assembly,
                                              current_time) +
                  Bcolors.ENDC, file=logfile)

            # Run the command and check for exit 21 (read length < K) and 255
            try:
                print(subprocess.check_output(cmd_assembly,
                                              stderr=subprocess.STDOUT,
                                              shell=True).decode("UTF-8"),
                      file=logfile)

                # Remove the K*/ directories when the assembly is complete
                cmd_removal = "rm -r {}K*/".format(params.outdir)
                print(Bcolors.OKGREEN +
                      "Removing the redundant ./K*/ directories now the "
                      "assembly is complete, with the following command: {}"
                      .format(cmd_removal) +
                      Bcolors.ENDC, file=logfile)

                print(subprocess.check_output(cmd_removal,
                                              stderr=subprocess.STDOUT,
                                              shell=True).decode("UTF-8"),
                      file=logfile)
            except subprocess.CalledProcessError as error:
                if error.returncode == 21:
                    # If exit code == 21, read length < K. Thus remove
                    # highest K-mer from the given list and retry the
                    # assembly.
                    if len(kmers.split(",")[:-1]) > 0:
                        updated_kmers = ",".join(kmers.split(",")[:-1])
                        highest_kmer = updated_kmers.split(",")[-1]
                        new_cmd_assembly = "spades.py {} -t {} -k {} " \
                                           "--restart-from k{} -o {}"\
                                           .format(tool,
                                                   threads,
                                                   updated_kmers,
                                                   highest_kmer,
                                                   params.outdir)

                        print(Bcolors.OKGREEN +
                              "SPAdes returned exit code 21, this means "
                              "that average read length < K. Now "
                              "continuing the assembly with the "
                              "following command: '{}'"
                              .format(new_cmd_assembly) +
                              Bcolors.ENDC, file=logfile)
                        print(subprocess.check_output(new_cmd_assembly,
                                                      stderr=subprocess
                                                             .STDOUT,
                                                      shell=True)
                                                      .decode("UTF-8"),
                              file=logfile)
                    else:
                        # This only happens when there are no K-mers left
                        raise ValueError("Average read length is too low "
                                         "to use the lowest specified "
                                         "K-mer. Please update the "
                                         "desired K-mer str")
                elif error.returncode == 255:
                    print(Bcolors.FAIL +
                          "SPAdes exited with code 255, when trying to "
                          "assemble the reads of accession {}"
                          "this is most likely "
                          "due to the fact that something is wrong with the "
                          "original reads of this accession. E.g. the reads "
                          "could in reality be single-end when the metadata "
                          "shows they should be paired-end."
                          .format(wildcards.accession) +
                          Bcolors.ENDC, file=logfile)

                    if not os.path.exists(config["root_dir"] +
                                          "corrupted_accessions.txt"):
                        subprocess.check_output("touch {}"
                                                .format(config["root_dir"] +
                                                "corrupted_accessions.txt"),
                                                shell=True)

                    # Write the accessions to a file with corrupted accessions
                    with open(config["root_dir"] +
                              "corrupted_accessions.txt", 'a+') as cor_acc:
                        cor_acc.write(wildcards.accession + "\n")

                    # Create a dummy file to produce the final output
                    subprocess.check_output("touch {}".format(output),
                                            shell=True)
                else:
                    print(error.output)
                    print(error.returncode)
                    print(subprocess.check_output("echo $?", shell=True))
                    sys.exit(error.returncode)


################## BLASTx against protein db using DIAMOND2 ###################
rule diamond:
    """Aligns all the assembled contigs against a given protein database """
    input:
        config["root_dir"] +
        config["assembly_dir"] +
        config["assembly"]["tool"].lower() + "/" +
        "{accession}/contigs.fasta"
    priority: 60
    threads: config["diamond"]["threads"] # default = 2
    params:
        diamond_db = config["diamond"]["db"],
        log = config["logs"]["root"] +
              config["logs"]["main"] +
              "main_{accession}.log"
    output:
        config["root_dir"] +
        config["diamond_dir"] +
        "{accession}_" +
        config["diamond"]["name"] +
        ".dmnd"
    run:
        with open(params.log, 'a') as logfile:
            cmd_diamond = "diamond blastx -d {} -q {} -p {} -o {} -f 6"\
                          .format(params.diamond_db, input,
                                  threads, output)

            # Print the used command to the logfile
            print(Bcolors.OKGREEN +
                  "Aligning the assembled reads of accession '{}' to the "
                  "database that can be found at: '{}'\n"
                  "Used command: '{}'".format(wildcards.accession,
                                              params.diamond_db,
                                              cmd_diamond) +
                  Bcolors.ENDC,
                  file=logfile)

            # Run the command
            try:
                print(subprocess.check_output(cmd_diamond,
                                              stderr=subprocess.STDOUT,
                                              shell=True).decode("UTF-8"),
                      file=logfile)
            except subprocess.CalledProcessError as err:
                if "First line seems to be blank" \
                        in err.output.decode("UTF-8"):
                    print(Bcolors.FAIL +
                          "The assembly file of {} seems to be empty. Likely, "
                          "this happened because some error occurred during "
                          "assembly. Now touching an empty DIAMOND output "
                          "file!".format(wildcards.accession) +
                          Bcolors.ENDC, file=logfile)

                    # Checking if the accession is in the corrupted file
                    if not subprocess.check_output("grep '{}' {}"
                                                   .format(wildcards.accession,
                                                           config["root_dir"] +
                                                           "corrupted_accessions.txt"),
                                                           shell=True):
                        with open(config["root_dir"] +
                                  "corrupted_accessions.txt", "a+") as cor_f:
                            cor_f.write(wildcards.accession + "\n")

                    # Touch the final output to let Snakemake finish
                    subprocess.check_output("touch {}".format(output),
                                            shell=True)
                else:
                    print(err.output, file=logfile)
                    print(err.returncode, file=logfile)
                    print(subprocess.check_output("echo $?",shell=True),
                          file=logfile)
                    sys.exit(err.returncode)
