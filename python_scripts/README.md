# Usage of main python scripts
## Python scripts index
All the python scripts that can be found in this repository.

### Scripts used within the pipeline

- `metadata_filter.py`: filters NCBI-SRA accessions based
  on sequencing library layout (paired-end/single-end).
- [`meta_accession_parser.py`](#meta_accession_parserpy): a script that can
  parse NCBI-SRA metadata based on commandline arguments.
  - `check_genome_assemblies.py`: helper module
  - `determine_kmer.py`: helper module
- `find_multiple_hosts.py`: filter NCBI TaxIDs
  on (sub-) species level.
- `sra_download.py`: time the downloads of the SRA data and restart if necessary.
- `read_lines.py`

### Pre- and postprocessing scripts

- [`DIAMOND_filter.py`](#diamond_filterpy): script to filter DIAMOND2
  alignment files (tabular format).
  - `translate_rna.py`: helper module, translates nucleotide FASTA
    sequences into protein sequences.  
- [`obtain_complete_genomes.py`](#obtain_complete_genomes.py): script to get all detected nucleotide sequences in a FASTA file
- [`fetch_protein_seq.py`](#fetch_protein_seqpy): script to fetch protein sequences in FASTA
  format from NCBI.
- [`viral_completeness.py`](#viral_completenesspy): used to assess viral genome
  completeness in DIAMOND2 alignment files or to patch fragmented viral
  genome assemblies.
- [`mapping_stats.py`](#mappingstats): used to determine host genome mapping
  stats based on the log files of the Snakemake workflow.


## Detailed documentation

Below you can find further detailed information for some of the scripts.

### meta_accession_parser.py

The pipeline only works with paired-end sequencing data and suitable accessions are determined automatically using by parsing SRA metadata.
`Snakefile_downloading` will expect a list of accessions containing libraries that can be both of paired-end or
single-end layout. It will then use the `meta_accession_parser.py`
to determine whether the data is of
paired- or single-end layout and sort the libraries based on this in two
separate text files. It will then continue the workflow and download only the
data of the libraries that contain paired-end reads.

### DIAMOND_filter.py
This standalone python script can be used to filter tabular DIAMOND output files. DIAMOND files can be given as input
using file paths (separated by spaces) or by providing a complete folder in which all DIAMOND files should be used.
The filtered output can be printed to stdout, or written to a given output file. When the option '-q' or '--quiet'
is provided, the script only provides the original tabular output (after filtering) and does not provide headers or file names.

Example usage:
```
python python_scripts/DIAMOND_filter.py --top_aid -D output/diamond/ -l 8000 12000 -o output/diamond.txt
```

```
usage: DIAMOND_filter.py [-h] (-d [[dmnd] ...] | -D [input-dir]) (-o <output> | -p) [-q] [-id [%AID]] [-l [<min_len> [[max_len] ...]]] [-al <alignment_len>] [-c [coverage]] [-e [E-value]] [--nodes | --top | --targets | --top_aid | --self] [--version]

This script takes one or multiple DIAMOND alignment files in BLAST tabular format. It will then filter the alignments files, based on the parameter settings given by the user and output all the alignments in one file.
Version: 2.0, December 2024

options:
  -h, --help            show this help message and exit
  -d, --diamond [[dmnd] ...]
                        Absolute pathway to one or multiple DIAMOND alignment files in tabular output format, given as a string separated by space.
  -D, --input-dir [input-dir]
                        Absolute pathway to an input directory that contains diamond output files with a '.dmnd' extension that should be used as input.
  -o, --output <output>
                        Pathway and name to use for the output file, given as a string. NOTE: if the given output already exists, the data will be appended.
  -p, --print           Print the output to the stdout.
  -q, --quiet           Minimal print statements to the screen will be shown. This will also exclude the header information from being printed to stdout.
  -id, --percid [%AID]  Floating point number between 0.0 and 100.0 used as the minimum percent alignment identity (%AID) between the query and the subject. Every alignment > this given %AID is saved to the output file.
  -l, --length [<min_len> [[max_len] ...]]
                        Integer, filter for alignments with an original contig length equal to or longer than the given number. Second argument is optional and can be used to filter for a max contig length.
                        Only alignments with original contig length equal to or below the max limit will be stored. Original contig length is given by SPAdes in number of nucleotides.
  -al, --alignment_length <alignment_len>
                        Integer, filter for alignments with a minimum length equal to or longer than the given number.
  -c, --coverage [coverage]
                        Floating point number used as the minimum K-mer coverage. All alignments with coverage > '--coverage' will be saved to the output file.
  -e, --evalue [E-value]
                        Floating point number used as the maximum E-value. Every alignment that has a E-value > 'eval' will be ignored.
  --top_aid             Obtains the best alignment for every unique node, based on alignment identity percentage (% AID).
  --top                 Obtain the contig with the best hit, based on %AID for every given input file. Other commands can be used to filter for e.g. alignment length or contig length.
  --nodes               When this option is given, only the unique nodes that are found after filtering are given as output. The output format is as: <file_name> tab <node_name> tab <file_path>
  --targets             This option shows all the accessions of the unique targets to which there is an alignment after filtering.
  --self                When this option is given, self alignments are removed from the output based on query and target name. I.e. if the query and target protein are the same and thus give 100% AID.
  --version             show program's version number and exit
```

### obtain_complete_genomes.py

This standalone script can be used save all detected nucleotide sequences in a FASTA file. It can be run after [`DIAMOND_filter.py`](#diamond_filterpy).

**Note:** The scripts requires Biopython.

Example usage, where `diamond.txt` was generated as in the example usage of [`DIAMOND_filter.py`](#diamond_filterpy):
```
python python_scripts/obtain_complete_genomes.py -input_accessions output/diamond.txt -input_folder output/assembly/rnaviralspades -nt_output output/nt.fasta
```
If you have getorf installed, you can run afterwards, for example:
```
getorf -sequence nt.fasta -outseq prot.fasta -minsize 6000 -find 1
```

```
usage: obtain_complete_genomes.py [-h] -input_accessions INPUT -input_folder INPUT_FOLDER -nt_output NT_OUTPUT

This script is used to obtain nucleotide sequence of interest from rnaviralSPAdes assembly files.

options:
  -h, --help            show this help message and exit
  -input_accessions INPUT
                        Str -- Absolute pathway to a text file containing the file directories, fasta header and additional description for every accession of interest.
                        The file should be tab separated and every new line should look like this: accession_id NODE_name rest (optional)
  -input_folder INPUT_FOLDER
                        str, assembly folder for the pipeline, typically: output/assembly/rnaviralspades/
  -nt_output NT_OUTPUT  Str -- Absolute pathway to the directory where all FASTA files containing the complete nucleotide sequence should be stored.
```


### fetch_protein_seq.py

This standalone script can be used to fetch protein sequences in FASTA format from NCBI.
It will print a multi-fasta to the screen.

Usage: `python3 fetch_protein_seq.py sequences_in.txt > output.fasta`, where
`sequences_in.txt` is a file with an NCBI protein IDs on every line.

### viral_completeness.py
This standalone python script can be used to determine the completeness of
*Iflaviridae* genomes found in the output DIAMOND2 alignment files. It can
also be used to 'patch' fragmented viral genomes. This script needs
`viral_completeness_config.yml` to run. All parameters in this config file
are [explained below](#viralcompletenessconfigyml). Next to this, make sure
that [Biopython](https://biopython.org/) is installed and that the assembly
files given as input for the patching are nucleotide sequences.

For the patching process, the general workflow of the script is as follows:
1. Filter for contigs with length > `min_contig_len` specified in the
   config file.
2. Determine the best target within the complete DIAMOND2 alignment file,
   based on the product of: </br>
   `contig length * % identity * alignment length`.
3. Obtain all the contigs with an alignment against the protein target ID
   of this best target.
4. Determine the viral completeness of the contigs based on the provided
   cut-offs in the config file. </br>
If the genome is fragmented:
5. Obtain the complete nt sequence of the contigs of interest from the
   SPAdes assembly file.
6. Translate the nt contig sequences into a protein sequence using EMBOSS
   getorf.
7. Use DIAMOND2 blastp to determine the exact alignment positions of the
   contig protein sequences against the best target.
8. Parse the alignment positions of the contigs from the DIAMOND2 blastp
   output files and try to concatenate the contigs using the best target
   polyprotein length as a template.
9. If there are contigs left with length > `min_contig_length`, restart at
   step 2.

**Please note** that this script will only work if it is able to determine the
length of every contig. This is done using the FASTA header in SPAdes
format, found at the query ID in the DIAMOND2 alignment files. e.g.
NODE_1_length_6525_cov_5.871304. It is important that the query ID of every
alignment contains the contig length in number of nucleotides at the 4th
position when it is split at underscores ('_'). For an example of a
DIAMOND2 file in the correct format.

Dependencies of the script:
* Configuration file: [viral_completeness_config.yml](#viralcompletenessconfigyml).
* An input file containing NCBI-SRA accession IDs on every new line.
* DIAMOND2 alignment files for these accessions.
* The NCBI metadata files of these accessions.
* Biopython module.
* SPAdes assembly files, containing contigs for every provided accession
  (for patching only).

Determine *Iflaviridae* genome completeness based on DIAMOND2 alignment files:
```
usage: python3 viral_completeness.py -c <./config.yml> [-s] [-p] viral_completeness -v <included_families>

This command can be used to determine the viral genome completeness and print it to the stdout
(in tabular format).

Key arguments:
  -c config.yml         str, relative pathway and name of the config file that can be used by the script
  -v included_families  str, viral families that should be included in the analysis separated by spaces.
                        Choose from the following five: [iflaviridae dicistroviridae marnaviridae
                        picornaviridae secoviridae calciviridae].  

Optional arugments:
  -s      When this option is given, the script will continue running if for some of the input
          accession, no metadata can be found.
  -p      This option will make the script print the output to the screen.  
```

Determine the viral completeness and patch fragmented *Iflaviridae* genomes:
```
usage: python3 viral_completeness.py -c <./config.yml> [-s] [-p] patching [-r] -v <included_families>

This command can be used to patch fragmented genome assemblies, based on the alignment positions
in the provided DIAMOND2 files.

Key arguments:
  -c config.yml         str, relative pathway and name of the config file that can be used by the script
  -v included_families  str, viral families that should be included in the analysis separated by spaces.
                        Choose from the following five: [iflaviridae dicistroviridae marnaviridae
                        picornaviridae secoviridae calciviridae].  

Optional arguments:
  -s      When this option is given, the script will continue running if for some of the input
          accession, no metadata can be found.
  -p      This option will make the script print the output to the screen.
  -r      This option will make the script always re-run all steps.
          Normally, when the script is re-run, it will check for the presence
          of the translated protein and blastp output files. If the files
          are already present, these steps will be skipped, to safe time
          and computation power. This option will always created new
          translated protein files and redo the DIAMOND2 blastp alignments.
```
