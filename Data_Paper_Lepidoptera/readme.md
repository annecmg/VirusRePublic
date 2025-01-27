# Additional data for paper

In this directory, you can find additional data described in the [paper](https://www.biorxiv.org/content/10.1101/2025.01.22.634245v1).

```
non-iflaviridae.fasta
```
Nucleotide sequences for the 479 contigs which have not been assigned to Iflaviridae but show similarity to another Picornavirales genome.
The sequence name gives the SRA library and the Picornavirales family that the contig was assigned to.
The sequence description gives the GenBank accession of the closest related genome for that contig.
Note that the only condition for including the contig is that it has a length of at least 5000nt and a significant DIAMOND hit (e-value threshold 0.001) to a Picornavirales genome.
The alignment length is not considered and might actually be short. Further analyses are needed to confirm that these are indeed viral contigs.

```
patched.fasta
```
Protein sequences of the 240 patched iflavirus genomes. The sequence name provides the SRA library and the sequence description provides the GenBank accession of the protein which was used as a template for patching.

```
complete.fasta
```
Protein sequences of the 1548 complete iflavirus genomes. The sequence name provides the SRA library and the sequence description provides the GenBank accession of the best DIAMOND hit. The description provides the position of the open reading frame on the genome (`complete_nt.fasta`).

```
complete_nt.fasta
```
Nucleotide sequences for the 1548 complete iflavirus genomes.

```
annotation.txt
```
Table with information on the 1812 sequences included in the Lepidoptera data set. 
