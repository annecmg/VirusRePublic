# VirusRePublic
This README contains information about how to use the Snakemake workflow 
**VirusRePublic**: VIRUS genome REconstruction from PUBLIC sequencing data.
The workflow will assemble RNA virus genomes from a list of NCBI SRA accession IDs.

# Installation
## Prerequisites

To install the contents of this repository, you need a valid conda and git installation: 

* Please follow the instructions of the corresponding tutorial for conda: https://docs.anaconda.com/miniconda/

* Please follow the instructions of the corresponding tutorial for git: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git

After installing conda, use the following command to install mamba:

```conda install -c conda-forge mamba```


To use the pipeline, you have to clone this repository using Git.

We suggest that you create a specific directory to store ViralGenomeConstrictor.
After accessing the directory through the command line, use the following command:

```git clone https://github.com/annecmg/ViralGenomeConstrictor.git```

After this step, you can configure snakemake as described in the next section.

## Conda environment
For the workflow to be able to work, it is important that all the tools 
that are used and the right versions are available. A conda environment 
file is provided called: `main_environment.yml` which contains all tools 
and version specifics. This file can be used to [set up and activate a 
suitable Conda environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file).

To create the environment, run the code below: 
```
mamba env create --prefix /desired/location/my_env_name -f main_environment.yml
```

Next to the main environment, the workflow needs a separate environment to 
be able to download the SRA sequencing data (some tools are not able to 
work with the latest version of Entrez Direct). This environment is called 
`sra_download.yml`. This Conda environment is created automatically when 
running the workflow. Please make sure that the pointer in the 
`main_config.yaml` file is set to the right location. 

# Running the pipeline

## Prerequisites

To be able to run the Snakemake workflow of this thesis project, there are 
some prerequisites. An extensive explanation for each of these can be found 
in the appendix. 
- [Configuration files](#configuration-files)
- [Conda environment](#conda-environment)
- [NCBI-SRA metadata (parser)](#ncbi-sra-metadata-and-metadata-parser)
- [API calls](#api-calls)
- [Local DIAMOND2 protein database](#local-diamond2-protein-database)

## Snakemake workflow setup and configuration

The main configuration file that is used by the workflow can be found at
`config/main_config.yaml`. This file is used by the workflow to determine where to 
find all the data (e.g. metadata that is downloaded) and where to store all 
the intermediate and final output data (e.g. assembly and DIAMOND files). 
It is advised to download all the scripts in this repository and keep the 
file structure like it is now, only change the in and output directories 
when running the workflow.  
Please make sure that the full path and name to this config file is 
provided in the first line of code in the `Snakefile_main.smk` like this: 
```configfile pointer
configfile: "/my_directory/for_config_files/main_config.yaml"
```
An explanation about all parameters in the file is [here](#main_configyaml).

Change the following two lines in `main_config.yml`:
* `root_dir:` the absolute path to the main folder of this GitHub repo
* `accession_file:` the absolute path to the input libraries.
This file contains library accessions on every new line. 

## Metadata

Before running the main workflow, either manually download the metadata of the desired libraries or 
(preferably) use the standalone snakemake workflow that is included in this 
repository (`/snake_files/Snakefile_metadata.smk`). 
After chancing the parameters in the config file as specified above, this 
Snakefile can be evoked as followed: 

```shell
conda activate vgc_env
snakemake -j 1 -s snake_files/Snakefile_metadata.smk --resources api_calls=2 -p --verbose
```

## Main usage with host genome
The complete Snakemake workflow that includes mapping to the host genome can be started with the following commandline arguments:

```shell
conda activate vgc_env
snakemake -s Snakefile_main.smk -j 8 -p --scheduler greedy --resources api_calls=2 disk_mb=15
```

`-s Snakefile_main.smk` is used to evoke the main Snakemake file of the 
workflow. 
`-j 8` will give the workflow a maximum of 8 cores to work with. 
`-p` (optional) Snakemake option to print all information to stdout. 
`--scheduler greedy` prevents a bug inside Snakemake, where no new rules 
will be activated after some random running time. 
`--resources api_calls=2 disk_mb=15`: disk_mb=15 prevents the workflow from downloading all 
the sequencing data of the given libraries at once. Setting the disk_mb to 
15 will limit the script to download and process the data of roughly 3 
libraries at once. Lowering the number will increase the number of 
libraries that are downloaded at once (and thus the disk size that is used).
The api_calls are further [explained below](#api-calls)

### Test case
The folder `/testing_data` contains files with test 
accessions that can be run to test the installation.

`test_accessions.txt`:
```
SRR11822934 #fragmented genome, host: Spodoptera exigua
SRR9678051 #complete genome, host: Plutella xylostella
SRR1050532 #2 complete and 1 patched genome, host: Spodoptera exigua
SRR7216194 # Single-end accession, will be removed
SRRTESTING # Non-existing accession ID, will be ignored.
```

`test_accessions_nohost.txt`:
```
DRR023324
DRR140180
```

The file `test_accessions_mixed.txt` contains all accessions in the 2 files above.

## Output
The workflow will create two main output folders: 
* `assembly/rnaviralspades/<accession>/contigs.fasta`: for every provided 
  valid library, the workflow will create a fasta file containing contigs
  assembled by [SPAdes](https://github.com/ablab/spades). 
* `diamond/<accession>_picornavirales.dmnd`: for every provided valid 
  library, the workflow will create a DIAMOND2 alignment file (tabular 
  format) against the provided local database.

In addition, supplementary output is created:
* `logs/`: A folder with all log files
* unmapped reads TODO
* `ref_genomes/`: A folder with HISAT2 indexes of the hosts.
These indexes will be re-used in subsequent runs with the same host.
If that is not needed, the folder can be removed.
* `extracted_reads/`: A folder with the unmapped reads for the libraries with an available host genome.

## Running the pipeline without mapping to the host genome (de novo mode)

Also in the de novo mode, first run the script to download the metadata, then start the snakefile:
```shell
conda activate vgc_env
snakemake -j 1 -s snake_files/Snakefile_metadata.smk --resources api_calls=2 -p --verbose
snakemake -s Snakefile_main_denovo.smk -j 8 -p --scheduler greedy --resources api_calls=2 disk_mb=15
```


## Running the pipeline with and without host mapping

You can determine automatically which accessions have a host genome:

1. download the metadata and move the output directory:
```
snakemake -j 1 -s snake_files/Snakefile_metadata.smk --resources api_calls=2 -p --verbose
mv output output-org
```
2. change the `accession_file` in the config file to `output-org/metadata/stats/reference_accessions.txt`, then run the pipeline for the libraries with assembly and move the output directory:
```shell
snakemake -j 1 -s snake_files/Snakefile_metadata.smk --resources api_calls=2 -p --verbose
snakemake -s Snakefile_main.smk -j 8 -p --scheduler greedy --resources api_calls=2 disk_mb=15
mv output output-host
```
3. change the `accession_file` in the config file to `output-org/metadata/stats/no_assembly_accessions.txt`, then run the pipeline for the libraries without assembly and move the output directory:
```shell
snakemake -j 1 -s snake_files/Snakefile_metadata.smk --resources api_calls=2 -p --verbose
snakemake -s Snakefile_main_denovo.smk -j 8 -p --scheduler greedy --resources api_calls=2 disk_mb=15
mv output output-nohost
```

# Appendices and additional information

## Local DIAMOND2 protein database
A local [DIAMOND2](https://github.com/bbuchfink/diamond) protein sequence 
database is required for the workflow to be able to align the assembled 
contigs against. This way, the workflow is able to create the final dmnd 
output files and determine which contigs are of interest and have similarity 
against proteins within the provided database. 

To create the local database, simply download all protein sequences of 
interest into one FASTA file. Using this file and the command below, 
create a local DIAMOND2 database: 
``` 
./diamond makedb --in my_proteins.fasta -d my_diamond_db 
```
`my_proteins.fasta`: the multi-FASTA file containing protein sequences of 
interest. </br> 
`my_diamond_db`: desired name for the binary DIAMOND2 database. 

Finally, make sure to set the pointer in `main_config.yaml` so the workflow 
will use the correct DIAMOND2 database.

## API calls
To automatically download the SRA data, the workflow use the 
[Entrez Direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/) commandline 
tool from NCBI. This tool uses API calls to retrieve the data from the NCBI 
servers. For safety reasons these calls are limited to 3 calls per second. 
The number of API calls can be set by `--resources api_calls=<n>`.
Setting the `api_calls` to a maximum of 2 will limit the snakefile to run a 
maximum of 2 API calls per second and thus prevent any 429 errors from 
occurring (when setting `api_calls=3` it will occasionally still give an 
error). 

Additionally, to increase the number of allowed API calls per second and thus 
the speed of the pipeline, the user can obtain a NCBI API key. 
Update the `NCBI_API_KEY` variable in ~/.profile. 
Please see the NCBI website for a full explanation on how to 
[obtain a NCBI API key](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/).  
When a NCBI API key is used, the snakefiles can be run with `api_calls` set 
to 9 per second. 

## Additional python scripts
Information about the usage of the additional pythons scripts can be found in the respective [directory](python_scripts/README.md).

## Snakemake scripts

The complete workflow consists of two main Snakemake scripts: 
- `Snakefile_main.smk`: main pipeline file.
- `Snakefile_main_denovo.smk`: main pipeline, without the host mapping step. 
  Viral genome assembly is *de novo*.

Snakemake scripts that are included in the main workflow: 
- `snake_files/Snakefile_downloading.smk`: used to download the NCBI-SRA raw sequencing data.
- `snake_files/Snakefile_host_genomes.smk`: used to download the RefSeq host genomes to map 
  against.
- `snake_files/Snakefile_read_processing.smk`: used to process the raw sequencing reads of each
  SRA library. 
- `snake_files/Snakefile_read_processing_denovo.smk`: used to process the raw sequencing reads
  of each SRA library when running the *Snakefile_main_denovo*.  

Additional Snakemake scripts: 
- (optional) `snake_files/Snakefile_metadata_manual.smk`: can be used to automatically 
  download and filter the metadata of the SRA libraries based on library 
  layout and host genome availability.


## Full explanation of the config files
### main_config.yaml
Below, all parameters that can be found in the `main_config.yaml` file are 
further explained. Please note that, for the Snakemake workflow to run, it is 
only needed to adjust the main parameters and the python scripts 
root directory. 

**Main parameters:**
* `root_dir: "my_working/directory/"` </br> 
  Specify the working directory for the workflow. All other directories 
  are based on this.
* `accession_file: "accession_list.txt"` </br>
  Relative pathway to a plain txt file that contains NCBI SRA accessions 
  on every new line. The data of these libraries is used by the workflow.

Scripts:
* `python_root: "folder_containing/python_scripts/"` </br>
  Directory where all the required python scripts can be found. 
* `librlayout_filter_script: "metadata_filter.py"` </br>
  Name and relative path to the script that filters the libraries based on 
  sequencing library layout. 
* `filter_host: "find_multiple_hosts.py"` </br>
  Name and relative path to the script that can filter for TaxIDs that contain 
  multiple species or are not defined at (sub-) species level. 
* `accession_parser: "meta_accession_parser.py"` </br>
  Name and relative path to the script that is used for metadata filtering. 
* `sra_download_script: "sra_download.py"` </br>
  Name and relative path to the script that is used to download the 
  NCBI SRA sequencing data of every provided accession. 
* `downloading_snakefile: "Snakefile_downloading.smk"` </br>
  Name and relative path to the Snakefile that downloads and splits the SRA 
  data. 
* `host_genome_snakefile: "Snakefile_host_genomes.smk"` </br>
  Name and relative path to the Snakefile that is used to automatically 
  download all the required host genome assemblies. 
* `read_processing_snakefile: "Snakefile_read_processing.smk` </br>
  Name and relative path to the Snakefile that is used to process the 
  downloaded sequencing reads (trimming, optional mapping, conversion to bam, 
  extraction, conversion to fastq)
* `read_processing_denovo_snakefile: "Snakefile_read_processing_denovo.smk` 
  </br>
  Name and relative path to the Snakefile that is used to process the 
  downloaded sequencing reads, when using the *de novo* version of the 
  workflow.

Data download:
* `conda_env: "sra_download.yml"` </br>
  Name and relative path to the environment file that contains the latest 
  version of SRA-tools (v3.0.3). This Conda environment is only activated 
  by Snakemake while downloading or splitting the SRA data.
* `dir: "download/"` </br>
  Name and relative path to the directory where all .sra files should be 
  downloaded.  
* `temp_files: "temp_files/"` </br> 
  Name and relative path to the directory where temporary files should be 
  stored that are created by Entrez while downloading the data.
* `timeout: 1200` </br>
  Integer, the number of seconds after which a download of sequencing data 
  should be stopped and re-started. The download is then continued where it 
  was left off. This prevents NCBI from throttling their servers and can 
  speed up the download of the sequencing files. 
* `threads: 5` </br>
  Integer, number of threads that Entrez is allowed to use while downloading. 

Metadata locations:
* `meta_root: "metadata/"` </br>
  Relative path (from root directory above) where the metadata should be 
  stored. 
* `stats: "stats/"` </br> 
  Relative (from meta_root) path to the directory where the statistics that 
  are produced from the metadata, should be stored. E.g. the files 
  containing filtered TaxIDs. 
* `paired_accessions: "paired_accessions/"` </br>
  Name and relative path to the directory where the metadata of the paired-end 
  libraries can be downloaded to. 
* `filtered_accessions: "filtered_accessions.txt"` </br> 
  Name and relative path to the file where the libraries should be stored 
  in which are filtered based on genome availability. 
* `paired_file: "paired_end_accessions.txt"` </br>
  Name and relative path to the file where all IDs should be stored of the 
  paired-end libraries. 
* `single_file: "single_end_accessions.txt"` </br>
  Name and relative path to the file where all IDs should be stored of the 
  single-end libraries. 
* `undeterminded: "undetermined_accessions.txt"` </br>
  Name and relative path to the file where all IDs should be stored of the 
  libraries of which the library layout could not be determined. 

(Temporary) data storage:
* `raw_reads: "raw_reads/"` </br> 
  Name and relative path (based on root_dir above) to the directory where all 
  raw reads should, temporarily, be stored. 
* `trimmed_reads: "trimmed_reads/"` </br> 
  Name and relative path to the directory where all trimmed reads should, 
  temporarily, be stored.
* `mapped_dir: "mapping/"` </br>
  Name and relative path to the directory where all mapping .sam files should, 
  temporarily, be stored.
* `extracted_dir: "extracted_reads/"` </br>
  Name and relative path to the directory where all extracted reads should, 
  temporarily, be stored.
* `assembly_dir: "assembly/"` </br>
  Name and relative path to the base directory where all the assembly 
  directories should be created. For every accession a subdirectory is 
  created. E.g. `./assembly/SRRXYZ/` where all SPAdes output is stored.
* `diamond_dir: "diamond/"` </br>
  Name and relative path to the directory where all the diamond alignment 
  output files should be stored.

Log files:
* `root: "~/main_data/logfiles/"` </br>
  Root pathway to the directory where all logfiles of the workflow should 
  be stored. A separate log file is created 
* `main: "main_logs/"` </br>
  Name and relative path to the directory where all the main log files should
  be stored. The main log file contains all information that the tools that 
  are used by the workflow normally print to the stdout. 
* `meta_data: "meta_data/"` </br>
  Name and relative path to the directory where the log files should be 
  stored that are produced during the filtering based on the metadata. E.g. 
  taxonomic host IDs and library layout.
* `download: "download/"` </br> 
  Name and relative path to the directory where the log files should be 
  stored that are produced during the download and splitting of the public 
  data.
* `trimming: "trimming/"` </br>
  Name and relative path to the directory the log files should be stored 
  that are produced by FastP during trimming of the reads.

Mapping parameters:
* `threads: 4` </br>
  Integer, number of threads that HISAT2 is allowed to use while performing 
  the mapping. 

Assembly parameters:
* `tool: "rnaviralspades"` </br>
  String, name of the tool that should be used to perform the assemblies. 
  Either 'rnaviralspades' or 'metaspades'. 
* `kmer_values: "33,55,77,99,127"` </br>
  The K-mers that are orignially given to SPAdes while performing the 
  mapping. The maximum K-mer length is 127. If the K-mer length > average 
  read length, the workflow will automatically delete the longest K-mer 
  value from the given list of values. 
* `threads: 16` </br>
  Integer, number of threads that SPAdes is allowed to use for a single 
  example at the time. 

DIAMOND2 parameters:
* `db: "./my_directory/diamond_db.dmnd"` </br> 
  Name and absolute path to the DIAMOND2 protein database that should be 
  used during the alignments. 
* `threads: 2` </br>
  Integer, number of threads that DIAMOND2 is allowed to use for every 
  alignment. 


### viral_completeness_config.yml
Below, all parameters that can be found in the `viral_completeness_config.yml` 
file are further explained. This config file should be used when running 
[viral_completeness.py](#viralcompletenesspy) to specify where to find all 
input and where to store output and log files. 

**Main parameters:**
* `input_accessions: "input_accessions.txt"` </br>
  Name and absolute path to a plain .txt file containing NCBI SRA accession 
  IDs on every new line. The data of these libraries will be used by the 
  script. 
* `temp_files: "./intermediate_files"` </br>
  Name and absolute path to a folder where the intermediate output should 
  be stored. (e.g. the FASTA files containing protein translations of the 
  nucleotide contigs given as input). 
* `log_files: "./log_files"` </br>
  Name and absolute path to a folder where the logfiles of the patching 
  process should be stored. 
* `min_p_id_alnlen: 30000` </br>
  Integer, the minimum product of the "%ID * alignment_length" for contigs 
  to be included in the analysis. The default setting of 30.000 is 
  specifically determined for *Iflaviridae* only. Contigs with a product < 
  30.000 often are alignments against non-viral products (no alignment 
  against complete viral polyproteins).
* `min_contig_len: 750` </br>
  Integer, minimum contig length (number of **nucleotides**) for a contig 
  to be included in the analysis. 
* `min_prot_len_nt: 200` </br>
  Integer, minimum length in **amino acids** for a contig to be included 
  after it's ORF is translated. 
* `metadata_loc: "./metadata` </br>
  Name and absolute path to a folder where the script can find the metadata 
  of every library provided in the input file. 
* `metadata_ext: "_metadata.txt"` </br> 
  The extension of every metadata file. The metadata files should be named 
  in a format like this: `./metadata_folder/SRRXYZ_metadata.txt`. 
  Metadata_loc is used as the base folder, then the accessions 
  provided in the input_accessions file and then the metadata_ext. 
* `diamond_prot_db: "./viral_prot_db.dmnd"` </br>
  Name and absolute path to a DIAMOND2 database containing protein 
  sequences of interest. This db will be used to do the blastp alignment 
  against. 
* `diamond_cores: 5` </br>
  Integer, number of cores DIAMOND2 is allowed to use for every blastp 
  alignment. 
* `diamond_blastp_out: "./blastp_output_folder"` </br>
  Name and absolute path to a folder where the DIAMOND2 blastp output files 
  should be stored. 
* `diamond_loc: "./diamond_files"` </br>
  Name and absolute path to a folder where the original input diamond files 
  can be found. 
* `diamond_ext: "_file_extionsion.dmnd"` </br>
  The extension that should be used to find the original diamond files. The 
  following format is used: diamond_loc + accession + diamond_ext. E.g.
  `./diamond_files/SRRXYZ_picornavirales.dmnd`. The accessions 
  provided in the input_accessions file are used as file prefix. 
* `assembly_loc: "./asssembly_files` </br>
  Name and absolute path to the base folder where the assembly files can be 
  found that containg the SPAdes contig (nt) assemblies. The following path 
  will be used: assembly_loc/ + accession/ + assembly_ext. E.g. `.
  /rnaviralspades/SRRXYZ/contigs.fasta`. This is folder format is the 
  general output of SPAdes. Please note that every library provided in 
  the input_accessions file should contain a separate folder with contigs.
  fasta file inside it. 
* `assembly_ext: "contigs.fasta"` </br> 
  The name of the file that contains the contig assemblies. (This name 
  should be the same for all provided accessions). 

Virus specific parameters per viral family: 
* `known_proteins: "./iflaviridae_known_proteins.txt"` </br>
  Name and absolute path to a plain .txt file containing NCBI protein IDs 
  that belong to this specific viral family. 
* `min_genome_len: 8000` </br>
  Integer, the minimum genome length of a contig (number of **nucleotides**)
  to be included in this viral family. 
* `max_genome_len: 12000` </br>
  Integer, the maximum genome length of a contig (number of **nucleotides**)
  to be included in this viral family.
