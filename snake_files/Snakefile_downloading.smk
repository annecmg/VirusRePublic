"""
Author: Devin van Valkengoed

Description: Please note that this is a helper module for Snakefile_main.smk
             which is the main workflow. This Snakemake file is
             used to automatically download SRA data.
Date: 01-Dec-2022

This file includes the following rules:
    * download_sra -- downloads the sra file corresponding to the accessions
    * split_sra -- splits the downloaded sra file into a forward (1) and
                   reverse (2) read fastq file.

TODO:
"""
########################## Downloading SRA data ###############################
rule download_sra:
    """Downloads .sra files via custom python script that uses NCBI datasets
    """
    input:
        filtered_accessions = config["root_dir"] +
        config["metadata"]["meta_root"] +
        config["metadata"]["stats"] +
        config["metadata"]["filtered_accessions"]
    output:
        temp(config["root_dir"] +
             config["download"]["dir"] +
             "{accession}/{accession}.sra")
    priority: 0
    params:
        script = config["scripts"]["python_root"] +
                 config["scripts"]["sra_download_script"],
        output_dir = config["root_dir"] +
                     config["download"]["dir"],
        time = config["download"]["timeout"],
    resources:
        api_calls = 1,
        disk_mb = 10
    log:
        logfile = config["logs"]["root"] +
                  config["logs"]["download"] +
                  "{accession}_download.log"
    shell:
        """
        (python3 {params.script} -a {wildcards.accession} -t {params.time}\
        -o {params.output_dir}) &> {log.logfile}
        """

################ Splitting the SRA files into separate reads ##################
rule split_sra:
    """Split downloaded .sra files using fasterq-dump """
    input:
        config["root_dir"] +
        config["download"]["dir"] +
        "{accession}/{accession}.sra"
    output:
        temp(config["root_dir"] +
             config["raw_reads"] +
             "{accession}_1.fastq"),
        temp(config["root_dir"] +
             config["raw_reads"] +
             "{accession}_2.fastq")
    priority: 10
    threads: config["download"]["threads"] # default = 5
    params:
        outdir = config["root_dir"] +
                 config["raw_reads"],
        temp_dir = config["root_dir"] +
                   config["download"]["temp_files"]
    resources:
        disk_mb = 5
    log:
        logfile = config["logs"]["root"] +
                  config["logs"]["download"] +
                  "{accession}_split.log"
    shell:
        """
        (fasterq-dump {input} --split-files -e {threads} -t\
        {params.temp_dir} -O {params.outdir}) &> {log.logfile}
        """
