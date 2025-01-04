"""
Author: Devin van Valkengoed
MSc Thesis Bioinformatic, Bioinformatics group, WUR

Description: Please note that this is a helper module for Snakefile_main.smk
             which is the main workflow of my thesis. This Snakemake file is
             used to automatically determine which host genomes are needed for
             mapping the reads of the given accession. This snakefile will
             download and index these genomes.
Date: 01-Dec-2022

This file includes the following rules:
    * genome_download -- downloads RefSeq genome assemblies, based on the given
                         accessions.
    * unzip_assemblies -- unzips the downloaded NCBI genome files,
                          renames and moves them.
    * rehydrate_assemblies -- Rehydrates the downloaded genome assemblies,
                              to create genomic files.
    * host_genome_indexing -- Creates a HISAT2 index of the genome files.

NOTE!! The file that contains the provided NCBI-SRA accessions given in
the configfile, can only contain accessions that have a host genome with a
RefSeq genome assembly available. (This can be checked using
meta_accession_parser.py script)

TODO:
    * Can the indices be .gz?
    * Import functions and class Bcolors
"""
######################### Import statements ###################################
import sys
import subprocess
import os.path
from datetime import datetime
from contextlib import redirect_stdout
sys.path.insert(1, "/" + config["scripts"]["python_root"])
from meta_accession_parser import Bcolors


####################### Downloading of genomes ################################
# Using NCBI datasets to download the host genomes
rule genome_download:
    """Downloads RefSeq host genomes, based on input NCBI-SRA accessions """
    output:
        temp(config["root_dir"] +
             config["host_genomes"] + "{taxid}/{taxid}.zip")
    priority: 20
    params:
        log = config["logs"]["root"] +
              config["logs"]["download"] + "{taxid}_genome_download.log"
    resources:
        api_calls = 1
    run:
        # Ensure the logs directory exists
        os.makedirs(os.path.dirname(params.log),exist_ok=True)

        with open(params.log, 'a') as logfile:
            with redirect_stdout(logfile):
                # Write the current date and time to the logfile
                now = datetime.now()
                current_time = now.strftime("%d-%b-%Y %H:%M:%S")
                print(current_time)

                # Create shell command
                cmd = (
                    "datasets download genome taxon {} --reference "
                    "\--filename {}{}/{}.zip --dehydrated"
                    .format(wildcards.taxid, config["root_dir"]
                            + config["host_genomes"], wildcards.taxid,
                            wildcards.taxid)
                )

                # Create shell command with log file
                shell_cmd = f"{cmd} &> {params.log}"
                try:
                    print(subprocess.check_output(shell_cmd,
                        stderr=subprocess.STDOUT,
                        shell=True).decode("UTF-8"))
                    print(Bcolors.OKGREEN +
                          "Genome of {} was downloaded successfully!\n"
                          .format(wildcards.taxid) +
                          Bcolors.ENDC)
                except subprocess.CalledProcessError as errcode:
                    if errcode.returncode != 0:
                        print(Bcolors.FAIL +
                              "NCBI datasets command line tool returned "
                              "non-zero exit code. This might be due to the "
                              "fact that given TaxID {} has no reference "
                              "assembly available. \n Please only provide "
                              "accessions with reference genome assembly!\n"
                              "NCBI datasets return message: \n"
                        .format({wildcards.taxid}) +
                              Bcolors.ENDC)
                        stdout = errcode.output
                        print(stdout.decode("UTF-8"))
                        sys.exit(214)
                    else:
                        print(Bcolors.FAIL +
                              "Some other error occurred during genome "
                              "download '{}'".format(errcode) +
                              Bcolors.ENDC)


rule unzip_assemblies:
    """Used to unzip the downloaded NCBI genome files, rename and move them """
    input:
        config["root_dir"] + config["host_genomes"] + "{taxid}/{taxid}.zip"
    output:
        directory(config["root_dir"] + config["host_genomes"] +
                       "{taxid}/ncbi_dataset/")
    priority: 25
    params:
        wdir = config["root_dir"] + config["host_genomes"] + "{taxid}/",
        log = config["logs"]["root"] +
              config["logs"]["download"] + "{taxid}_genome_download.log"
    run:
        with open(params.log, 'a') as log_file:
            print(Bcolors.OKGREEN +
                  "Unzipping the genome files of: '{}'."
                  .format(wildcards.taxid) +
                  Bcolors.ENDC, file=log_file)
            unzip = "unzip {}{}.zip -d {}".format(params.wdir,
                                                  wildcards.taxid,
                                                  params.wdir)
            print(subprocess.check_output(unzip, stderr=subprocess.STDOUT,
                                          shell=True).decode("UTF-8"),
                  file=log_file)


rule rehydrate_assemblies:
    """Rehydrates the downloaded genome assemblies, to create genomic files
    """
    input:
        config["root_dir"] + config["host_genomes"] + "{taxid}/ncbi_dataset/"
    output:
        temp(config["root_dir"] + config["host_genomes"] +
             "{taxid}/{taxid}_genome.fna")
    priority: 30
    params:
        wdir = config["root_dir"] + config["host_genomes"] + "{taxid}/",
        log=config["logs"]["root"] +
            config["logs"]["download"] + "{taxid}_genome_download.log"
    run:
        with open(params.log, 'a') as logf:
            # Rehydrating the NCBI genome files
            rehydrate_cmd = "datasets rehydrate " \
                            "--directory {}".format(params.wdir)
            print(Bcolors.OKGREEN +
                  "Rehydrating the downloaded NCBI file of {} using the "
                  "following command: '{}'."
                  .format(wildcards.taxid, rehydrate_cmd) +
                  Bcolors.ENDC, file=logf)
            print(subprocess.check_output(rehydrate_cmd,
                                          stderr=subprocess.STDOUT,
                                          shell=True).decode("UTF-8"),
                  file=logf)

            # Renaming the files for ease of use
            target = str(params.wdir + wildcards.taxid + "_genome.fna")
            rename_cmd = "find {}ncbi_dataset/data -type f -name " \
                         "*genomic.fna -exec mv ".format(params.wdir) + \
                         "{} " + \
                         "{} \; -quit".format(target)

            print(Bcolors.OKGREEN +
                  "Renaming the genome assembly file of {} into: "
                  "'{}{}_genome.fna'".format(wildcards.taxid,
                                             params.wdir,
                                             wildcards.taxid) +
                  Bcolors.ENDC, file=logf)
            print(subprocess.check_output(rename_cmd, stderr=subprocess.STDOUT,
                                          shell=True).decode("UTF-8"),
                  file=logf)


rule host_genome_indexing:
    """Creates a HISAT2 index of the downloaded genome files """
    input:
        config["root_dir"] + config["host_genomes"] +
        "{taxid}/{taxid}_genome.fna"
    output:
        config["root_dir"] + config["host_genomes"] +
        "{taxid}/{taxid}_index.1.ht2"
    priority: 35
    params:
        log = config["logs"]["root"] +
              config["logs"]["download"] + "{taxid}_genome_download.log"
    run:
        with open(params.log, 'a') as loging_f:
            hisat_build = "hisat2-build -f {} {}{}/{}_index"\
                          .format(input,
                                  config["root_dir"] + config["host_genomes"],
                                  wildcards.taxid,
                                  wildcards.taxid)

            print(Bcolors.OKGREEN +
                  "Building a HISAT2 index for genome file: '{}'.\n"
                  "Only the the index files will be saved.\n"
                  "Command that was used: '{}'".format(input, hisat_build) +
                  Bcolors.ENDC, file=loging_f)

            print(subprocess.check_output(hisat_build,
                  stderr=subprocess.STDOUT, shell=True).decode("UTF-8"),
                  file=loging_f)

            print(Bcolors.OKGREEN +
                  "Finished the download and indexing of genome '{}'!\n"
                  "Current time: \n"
                  "{}\n".format(wildcards.taxid,
                              datetime.now().strftime("%d-%b-%Y %H:%M:%S")) +
                  Bcolors.ENDC, file=loging_f)
