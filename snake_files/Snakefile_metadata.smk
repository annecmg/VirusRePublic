"""
Author: Devin van Valkengoed
Description: Snakemake file to automatically download the metadata of provided
             accession belonging to NCBI-SRA datasets and filter the data based
             on single- and paired-end library layout.
Date: 03-Nov-2022 V1
      26-Jan-2023 V2.0
      20-Mar-2024 V3.0 - Implementation in main workflow
Usage: snakemake -j 1 -s snake_files/Snakefile_metadata.smk
       --resources api_calls=2 -p --verbose

Key arguments:
    api_calls -- int, maximum number of API calls the script is allowed to make
               per second. Set to 2, or increase to 9 when using a
               NCBI-API key.

Optional arguments:
    cores -- int, number of cores to use for parallel processing. Note that
             this snakemake file cannot make use of parallel processing, so
             there is no benefit of setting cores > 1. The custom python script
             that is used for the FFQ check does use parallel processing.

Dependencies (scripts):
    * meta_accession_parser.py
    * find_multiple_hosts.py

Please also download the config.yaml file, update the parameters in the config
file according to your needs and update the pointer to the config file in the
script below.

TODO:
    * Include the script in Snakefile_main.smk to automatically run the full
      pipeline if the metadata is not downloaded yet.
    * Restructure rules: determine_abnormal and filter_abnormal as filter is
      dependend on the output of determine, but the output of determine is not
      further used.
    * Reset the config location after implementing in main workflow.
"""
######################### Config file #########################################
# Location of provided config file
import subprocess
import sys
import os
from datetime import datetime
from contextlib import redirect_stdout
from pathlib import Path
configfile: "config/main_config.yaml"

# Add python script directory to path for function import
python_scripts_dir = config["root_dir"] + config["scripts"]["python_root"]
sys.path.insert(1, python_scripts_dir)

# Import read lines function
from read_lines import read_lines

############################ Functions ########################################
def read_metadata(meta_data):
    """Reads the provided metadata file into a nested list

    Key arguments:
        meta_data -- str, pathway to a .txt file providing the metadata of
                     SRA experiments in the standard NCBI format.
    Returns:
        file_lines -- lst, nested list containing two list. First list is the
                      headers of the NCBI-SRA metadata file and the second list
                      contains the values.
    """
    in_file = open(meta_data, 'r')
    file_lines = []

    # Save the two lines of a NCBI-SRA metadata file into a list
    for line in in_file:
        file_lines.append(line.strip().split(","))

    in_file.close()

    return file_lines


def obtain_taxid(meta_data_list):
    """Obtains the taxid from an NCBI-SRA accession metadata file

    Key arguments:
        meta_data_list -- lst, nested list containing two list. First list is
                          the headers of the NCBI-SRA metadata file and the
                          second list contains the values.

    Returns:
        tax_id -- str, taxonomic identifier as determined by NCBI-SRA obtained
                  from the provided metadata.
    """
    # Check if current taxid is at the correct position and not empty
    if meta_data_list[0][27] == "TaxID":
        if meta_data_list[1][0] != '':
            tax_id = meta_data_list[1][27]
        else:
            raise ValueError("Exception while parsing metadata. "
                             "No TaxID is provided in the metadata "
                             "file. Please double check if the metadata file "
                             "provided is in NCBI-SRA format!")
    else:
        raise ValueError("Exception while parsing metadata. "
                         "TaxID is not at the expected position "
                         "and cannot be determined.")

    return tax_id


################## Config file checks and variables ###########################
def print_config_info():
    """Print config info to stdout when starting Snakemake script. """
    print("Number of used total accessions: {}".format(len(
        read_lines(config["accession_file"]))))
    print(f"Root directory set to:{config['root_dir']}\n")

    # Log the used conda env
    conda_env = subprocess.check_output("echo $CONDA_DEFAULT_ENV",
        shell=True).decode("UTF-8")
    print(f"Conda environment that was used running the workflow: {conda_env}")

    # Log the sys.path
    print(f"Sys PATH locations that were used: {sys.path}")


onstart:
    print_config_info()

#################### Defining output for every rule ###########################
rule all:
    input:
        # # Filtering based on library layout
        # config["root_dir"] +
        # config["metadata"]["meta_root"] +
        # config["metadata"]["stats"] +
        # "liblayout_parsed.temp",
        # Downloading all the metadata of the given accessions
        # config["root_dir"] +
        # config["metadata"]["meta_root"] +
        # config["metadata"]["stats"] +
        # "downloaded_metadata.txt",
        # # Obtaining the unique TaxIDs
        # config["root_dir"] +
        # config["metadata"]["meta_root"] +
        # config["metadata"]["stats"] +
        # "paired_unique_taxids.txt",
        # Filtering for abnormal TaxIDs
        config["root_dir"] +
        config["metadata"]["meta_root"] +
        config["metadata"]["stats"] +
        config["metadata"]["filtered_accessions"]

############### Filtering accessions based in library layout ##################
rule filter_librlayout:
    """Uses a custom python script that uses FFQ to filter data on library type
    """
    input:
        config["accession_file"]
    priority: 100
    params:
        script = config["root_dir"] +
                 config["scripts"]["python_root"] +
                 config["scripts"]["librlayout_filter_script"],
        undetermined = config["root_dir"] +
                       config["metadata"]["meta_root"] +
                       config["metadata"]["stats"] +
                       config["metadata"]["undetermined"],
        single= config["root_dir"] +
                config["metadata"]["meta_root"] +
                config["metadata"]["stats"] +
                config["metadata"]["single_file"],
    output:
        paired = config["root_dir"] +
                 config["metadata"]["meta_root"] +
                 config["metadata"]["stats"] +
                 config["metadata"]["paired_file"],

    log: config["root_dir"] +
         config["logs"]["root"] +
         config["logs"]["meta_data"] +
         "main_accession_filter.log"
    shell:
        """
        (python3 {params.script} -a {input} -s {params.single} \
        -p {output.paired} -u {params.undetermined}) >> {log} 2>&1
        """


########################## Downloading metadata ###############################
rule metadata_download:
    """Downloads the metadata of the paired-end accessions """
    input:
        paired_accessions = config["root_dir"] +
                            config["metadata"]["meta_root"] +
                            config["metadata"]["stats"] +
                            config["metadata"]["paired_file"]
    output:
        temp(config["root_dir"] +
        config["metadata"]["meta_root"] +
        config["metadata"]["stats"] +
        "downloaded_metadata.txt")
    priority: 50
    params:
        log = config["root_dir"] +
              config["logs"]["root"] +
              config["logs"]["meta_data"] +
              "main_accession_filter.log"
    resources:
        api_calls = 1
    run:
        # Create log file if it is not present yet
        if not os.path.exists(params.log):
            os.makedirs(os.path.dirname(params.log), exist_ok=True)
            with open(params.log, 'w'):
                pass

        with open(params.log, 'a+') as logfile:
            # Obtain current date and time
            now = datetime.now()
            current_time = now.strftime("%d-%b-%Y %H:%M:%S")
            print(current_time, file=logfile)

            # Create directory to store the metadata of paired-end accessions
            # if it doesn't exist yet
            paired_folder = config["root_dir"] + \
                            config["metadata"]["meta_root"] + \
                            config["metadata"]["paired_accessions"]
            if not os.path.exists(paired_folder):
                os.makedirs(paired_folder, exist_ok=True)

            # Obtain the metadata of every accession
            for accession in read_lines(str(input.paired_accessions)):
                metadata_file = "{}{}_metadata.txt"\
                                .format(paired_folder,
                                        accession)
                # If not, create the files
                if not os.path.exists(metadata_file):
                    cmd_download = "esearch -db sra -query {} | efetch " \
                                   "-format runinfo > {}"\
                                   .format(accession, metadata_file)
                    subprocess.check_output(cmd_download, shell=True)
                    print(cmd_download, file=logfile)
                else:
                    print("The metadata file: '{}' "
                          "already exists. Not overwriting the file!"
                          .format(metadata_file), file=logfile)

                # Check if the downloaded file is empty or contains data
                if os.path.getsize(metadata_file) == 0:
                    print("The downloaded file is empty. "
                          "The metadata of accession {} may be "
                          "invalid.".format(accession), file=logfile)

            # Create the final dummy output file for Snakemake
            subprocess.check_output("touch {}".format(output), shell=True)


########### Create a file with unique TaxIDs for genome download ##############
rule obtain_taxids:
    """Creates a list of unique TaxIDs found in the list of accessions """
    input:
        paired_accessions = config["root_dir"] +
                            config["metadata"]["meta_root"] +
                            config["metadata"]["stats"] +
                            config["metadata"]["paired_file"],
        metadata_download = config["root_dir"] +
                            config["metadata"]["meta_root"] +
                            config["metadata"]["stats"] +
                            "downloaded_metadata.txt"
    output:
        unique_taxids = config["root_dir"] +
                        config["metadata"]["meta_root"] +
                        config["metadata"]["stats"] +
                        "paired_unique_taxids.txt",
        undefined_taxids = config["root_dir"] +
                           config["metadata"]["meta_root"] +
                           config["metadata"]["stats"] +
                           "undefined_taxids.txt"
    params:
        script = config["root_dir"] +
                 config["scripts"]["python_root"] +
                 config["scripts"]["accession_parser"],
        metadata_location = config["root_dir"] +
                            config["metadata"]["meta_root"] +
                            config["metadata"]["paired_accessions"],
        log = config["root_dir"] +
              config["logs"]["root"] +
              config["logs"]["meta_data"] +
              "main_accession_filter.log"
    run:
        with open(params.log, 'a+') as logfile, redirect_stdout(logfile):
            with open(str(input.paired_accessions), "r") as accession_file:
                for accession in accession_file:
                    accession = accession.strip()
                    cmd_unique_taxid = "python3 {} taxid -m {}_metadata.txt " \
                                       "-o {} -u {}"\
                        .format(params.script,
                                params.metadata_location + accession,
                                output.unique_taxids,
                                output.undefined_taxids)

                    # Try to obtain the taxonomic identifier of the host
                    try:
                        subprocess.check_output(cmd_unique_taxid, shell=True)
                    except subprocess.CalledProcessError as e:
                        print("The TaxID of accession {} could not be "
                              "determined. Adding this accession to the "
                              "undefined TaxID "
                              "file: {}.".format(accession,
                                                 output.undefined_taxids))
                        with open(output.undefined_taxids, 'a+') as undefined:
                            undefined.write(accession + "\n")

                # Creating dummy output if files are not yet created
                for file in output:
                    if not os.path.exists(file):
                        open(file,'a').close()
                        print("Created a dummy file: {}".format(file),
                              file=logfile)




############## Determine the abnormal TaxIDs present in the data ##############
rule determine_abnormal:
    """Determines TaxIDs that are wrongly defined in the metadata """
    input:
        unique_taxids = config["root_dir"] +
                        config["metadata"]["meta_root"] +
                        config["metadata"]["stats"] +
                        "paired_unique_taxids.txt"
    output:
        abnormal_taxids = config["root_dir"] +
                          config["metadata"]["meta_root"] +
                          config["metadata"]["stats"] +
                          "abnormal_TaxIDs.txt",
        refseq_taxid = config["root_dir"] +
                       config["metadata"]["meta_root"] +
                       config["metadata"]["stats"] +
                       "reference_taxids.txt",
        assembly_taxid = config["root_dir"] +
                         config["metadata"]["meta_root"] +
                         config["metadata"]["stats"] +
                         "genbank_taxids.txt",
        no_assembly_taxid = config["root_dir"] +
                            config["metadata"]["meta_root"] +
                            config["metadata"]["stats"] +
                            "no_assembly_taxids.txt",
    params:
        log = config["root_dir"] +
              config["logs"]["root"] +
              config["logs"]["meta_data"] +
              "main_accession_filter.log",
        filter_script = config["root_dir"] +
                        config["scripts"]["python_root"] +
                        config["scripts"]["filter_host"],
        main_script = config["root_dir"] +
                      config["scripts"]["python_root"] +
                      config["scripts"]["accession_parser"]
    run:
        with open(params.log, "a+") as logfile:
            # Determine the abnormal TaxIDs found in the list of unique TaxIDs
            subprocess.check_output("python3 {} -t {} -o {}"
                                    .format(params.filter_script,
                                            input.unique_taxids,
                                            output.abnormal_taxids),
                                    shell=True)

            # Determine the TaxIDs that have a RefSeq assembly
            subprocess.check_output("python3 {} genome -t {} -o1 {} -o2 {} "
                                    "-o3 {} -q"
                                    .format(params.main_script,
                                            input.unique_taxids,
                                            output.refseq_taxid,
                                            output.assembly_taxid,
                                            output.no_assembly_taxid),
                                            shell=True)

            # Creating dummy output if no accession belong to a category
            for file in output:
                if not os.path.exists(file):
                    open(file, 'a').close()


#################### Filtering for abnormal accessions ########################
rule filter_abnormal:
    """Filters accessions based on genome availability and abnormal TaxIDs """
    input:
        unique_taxids = config["root_dir"] +
                        config["metadata"]["meta_root"] +
                        config["metadata"]["stats"] +
                        "paired_unique_taxids.txt",
        paired_accessions = config["root_dir"] +
                            config["metadata"]["meta_root"] +
                            config["metadata"]["stats"] +
                            config["metadata"]["paired_file"],
        abnormal_taxids = config["root_dir"] +
                          config["metadata"]["meta_root"] +
                          config["metadata"]["stats"] +
                          "abnormal_TaxIDs.txt",
        refseq_taxid = config["root_dir"] +
                       config["metadata"]["meta_root"] +
                       config["metadata"]["stats"] +
                       "reference_taxids.txt",
        genbank_taxid = config["root_dir"] +
                         config["metadata"]["meta_root"] +
                         config["metadata"]["stats"] +
                         "genbank_taxids.txt",
        no_assembly_taxid = config["root_dir"] +
                            config["metadata"]["meta_root"] +
                            config["metadata"]["stats"] +
                            "no_assembly_taxids.txt",
        undefined_taxids= config["root_dir"] +
                          config["metadata"]["meta_root"] +
                          config["metadata"]["stats"] +
                          "undefined_taxids.txt"
    output:
        filtered_accessions = config["root_dir"] +
                              config["metadata"]["meta_root"] +
                              config["metadata"]["stats"] +
                              config["metadata"]["filtered_accessions"],
        abnormal_host_accessions = config["root_dir"] +
                                   config["metadata"]["meta_root"] +
                                   config["metadata"]["stats"] +
                                   "abnormal_host_accessions.txt",
        refseq_accessions = config["root_dir"] +
                            config["metadata"]["meta_root"] +
                            config["metadata"]["stats"] +
                            "reference_accessions.txt",
        genbank_accession = config["root_dir"] +
                            config["metadata"]["meta_root"] +
                            config["metadata"]["stats"] +
                            "genbank_accessions.txt",
        no_assembly_accession = config["root_dir"] +
                                config["metadata"]["meta_root"] +
                                config["metadata"]["stats"] +
                                "no_assembly_accessions.txt"
    params:
        corrupted_meta = config["root_dir"] +
                         config["metadata"]["meta_root"] +
                         config["metadata"]["stats"] +
                         "corrupted_metadata_accessions.txt",
        metadata_location = config["root_dir"] +
                            config["metadata"]["meta_root"] +
                            config["metadata"]["paired_accessions"],
        log = config["root_dir"] +
              config["logs"]["root"] +
              config["logs"]["meta_data"] +
              "main_accession_filter.log"
    run:
        with open(params.log, "a+") as logfile:
            # Read accessions, blacklisted and TaxIDs with RefSeq into memory
            blacklisted_taxids = read_lines(input.abnormal_taxids)
            refseq_taxids = read_lines(input.refseq_taxid)
            genbank_taxids = read_lines(input.genbank_taxid)
            no_assembly_taxids = read_lines(input.no_assembly_taxid)
            accession_list = read_lines(input.paired_accessions)
            undefined_accession = read_lines(input.undefined_taxids)

            # Filter the accessions based on TaxID
            for acc in accession_list:
                # Skip if the TaxID could not be determined
                if acc in undefined_accession:
                    print("TaxID of accession: {} could previously not be "
                          "determined, now skipping".format(acc),
                          file=logfile)
                    continue
                # Obtain the TaxID of the current accession
                current_metadata = read_metadata("{}{}_metadata.txt"
                                              .format(params.metadata_location,
                                                      acc))
                try:
                    current_taxid = obtain_taxid(current_metadata)
                    # Skipp if the TaxID is blacklisted (multiple species e.g.)
                    if current_taxid in blacklisted_taxids:
                        print("{} has TaxID: {} which is in the TaxID "
                              "blacklist, discarding this accession from the "
                              "filtered accession file!".format(acc,
                                                                current_taxid),
                              file=logfile)
                        with open(output.abnormal_host_accessions,"a+") as outfile:
                            outfile.write(acc + "\n")
                    # # Skipp if no RefSeq genome assembly is available for host
                    # elif current_taxid not in refseq_taxids:
                    #     print("{} has TaxID: {} which doesn't have a RefSeq "
                    #           "genome assembly available! Adding the "
                    #           "accession to the filtered file"
                    #     .format(acc,current_taxid),file=logfile)
                    # Host has a RefSeq assembly
                    elif current_taxid in refseq_taxids:
                        print("{} has TaxID: {} which has a RefSeq assembly, "
                              "appending to the filtered and RefSeq "
                              "file!".format(acc, current_taxid), file=logfile)
                        with open(output.filtered_accessions,"a+") as outfile:
                            outfile.write(acc + "\n")
                        with open(output.refseq_accessions,"a+") as outfile:
                            outfile.write(acc + "\n")

                    # Host has Genbank genome assembly
                    elif current_taxid in genbank_taxids:
                        print("{} has TaxID: {} which has a Genbank assembly, "
                              "appending to the filtered and Genbank "
                              "file!".format(acc,current_taxid),file=logfile)
                        with open(output.filtered_accessions,"a+") as outfile:
                            outfile.write(acc + "\n")
                        with open(output.genbank_accession,"a+") as outfile:
                            outfile.write(acc + "\n")

                    # Host has no assembly
                    elif current_taxid in no_assembly_taxids:
                        print("{} has TaxID: {} which has a no genome "
                              "assembly, appending to the filtered and 'no "
                              "assembly' file!".format(acc,current_taxid),
                            file=logfile)
                        with open(output.filtered_accessions,"a+") as outfile:
                            outfile.write(acc + "\n")
                        with open(output.no_assembly_accession,"a+") as outfile:
                            outfile.write(acc + "\n")

                    ## The lines below should now be fixed at the start of
                    # the accession loop.
                    # # Else, raise error
                    # else:
                    #     raise ValueError("No TaxID could be determined for "
                    #                      "accession: {}".format(acc))
                except ValueError:
                    raise ValueError
                # Skipp if the metadata is corrupted/empty
                # except IndexError:
                #     print("The metadata file of {} is corrupted, discarding "
                #           "this accession from the final filtered accession "
                #           "file. Adding the accession to ")
                #     with open(params.corrupted_meta, "a+") as corrupted:
                #         corrupted.write(acc + "\n")
                #     continue

            # Creating dummy output if no accession belong to a category
            for file in output:
                if not os.path.exists(file):
                    open(file,'a').close()
