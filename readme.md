
# KleTy (Klebsiella typer for the core genome and plasmids)
Typing engine for core genome and plasmids in Klebsiella


# INSTALLATION:

KleTy was developed and tested in Python >=3.8. It depends on several Python libraries: 
~~~~~~~~~~
click
numba
numpy
pandas
biopython
pyarrow
fastparquet
~~~~~~~~~~

All libraries can be installed using pip: 

~~~~~~~~~~
pip install click numba numpy pandas biopython pyarrow fastparquet
~~~~~~~~~~
KleTy also calls two 3rd party programs:

~~~~~~~~~~
ncbi-blast+
diamond
~~~~~~~~~~

Both can be installed via 'apt' in UBUNTU:
~~~~~~~~~~
sudo apt install -y ncbi-blast+ diamond
~~~~~~~~~~

The whole environment can also be installed in conda:


~~~~~~~~~~
conda create --name dty python==3.11
conda activate dty
conda install -c conda-forge biopython numba numpy pandas click pyarrow fastparquet
conda install -c bio-conda blast diamond
~~~~~~~~~~

The installation process normally finishes in <10 minutes. 

NOTE: Please make sure that "makeblastdb", "blastn", and "diamond" are all in the PATH environment variable (can be run without pointing to their actual location). 


# Quick Start (with examples)
## Get allelic and HierCC callings
~~~~~~~~~~~
$ cd /path/to/KleTy/
$ python KleTy.py -q examples/CP015990.fna
~~~~~~~~~~~

The whole calculation finishes in <3 minutes with 8 CPU threads. 


# USAGE:
## KleTy.py - allelic callings and HierCC clusters & species predictions

~~~~~~~~~~~~~~
$ Usage: KleTy.py [OPTIONS]

Options:
  -q, --query TEXT      query genome in fasta or fastq format. May be gzipped.
                        [required]
  -o, --prefix TEXT     prefix for output. default: query filename
  -n, --n_proc INTEGER  number of process to use. default: 8
  -m, --skip_gene       flag to skip AMR/VF searching. default: False
  -g, --skip_cgmlst     flag to skip cgMLST. default: False
  -p, --skip_plasmid    flag to skip plasmid typing. default: False
  --help                Show this message and exit.
~~~~~~~~~~~~~~~~~

# Outputs:
## KleTy generates two files:

~~~~~~~~~~~~~
<prefix_of_input>.KleTy
<prefix_of_input>.alleles.gz
~~~~~~~~~~~~~

### <prefix_of_input>.KleTy contains the genotyping results
~~~~~~~~~~~~~

~~~~~~~~~~~~~



### KleTy also generates an allelic calling results (examples/GCF_005221305.alleles here) like:

~~~~~~~~~~~~~
$ head -4 examples/GCF_005221305.alleles
>CCD84_RS00010 value_md5=f7d49f149ad872340b86448e3602f565 CIGAR=CCD84_RS00010_1:264M39D48M9I12M12I96M3D18M9I1089M accepted=2 reference=GCF_005221305.fna identity=0.846 coordinates=NZ_CP039887.1:1..1569
ATGACGTTAGCTGAATTTTGGCCGCTGTGCCTTCGCCGCCTTCACGAAATGTTGCCTGCCGGGCAGTTTGCGCAATGGATTGCGCCTTTGACCGTGGGCGAAGAAAACGGCGTATGGGTGGTGTATGGTAAAAACCAATTTGCCTGCAATATGCTCAAAAGCCAGTTTGCCGCCAAAATTGACGCCGTGCGTGCCGAATTAGTGCCTCAGCAGGCTGCTTTTGCGTTTAAGCCGGGCGTAGGTACGCATTATGAAATGGCGGCTCAGACTGTTGCGCCGGTGCAAGTGCAAGAGGTCATTGAAGTTGAAGAGTGTGTAGAGCCTGTTCAAATGCCTTTGCAAACTGCTGCGCCAATGGAAGAAAATAGGCCGTCTGAAACGGTTTCCAAACCTGCAGCTGCCATGACGGCTGCCGAGATTTTGGCGCAACGCATGAAAAACCTGCCGCATGAGCCTCAAGTGCAAACTACTGCTTCGGCTGAATCTAAAGCAGTTGCCAAAGCCAAAACCGACGCGCAACACGATGCGGAAGAAGCGCGCTACGAACAGACCAATCTGTCGCGTGACTATACATTTGAAACTTTGGTGGAAGGTAAGGGCAACCGCCTTGCTGCCGCAGCCGCCCAGGCGATTGCTGAAAATCCGGGGCAGGGCTACAACCCGTTTTTCTTATACGGCAGTACCGGTTTGGGTAAAACCCACTTGGTGCAAGCCATCGGCAACGAATTGCTGAAAAACCGTCCTGATGCCAAAGTGCGCTATATGCACTCGGATGACTATATCCGCAGCTTTATGAAGGCGGTGCGCAACAATACTTACGATGTATTCAAGCAACAATACAAACAATATGACCTGCTGATTATCGACGATATCCAGTTCATCAAAGGCAAAGACCGTACGATGGAAGAATTCTTTTATCTGTACAACCATTTTCACAATGAGAAAAAACAACTTATCCTGACGTGCGACGTATTGCCTGCCAAAATCGAAGGTATGGATGACCGTCTCAAATCCCGTTTCTCATGGGGTTTGACTTTGGAACTCGAACCGCCCGAATTGGAAATGCGCGTGGCGATTTTGCAGAAAAAGGCAGAAGCGGCCGGTATCAGTATCGAAGACGAAGCCGCTCTGTTTATCGCCAATCTGATCCGTTCCAATGTGCGTGAGCTGGAAGGCGCGTTCAACCGCGTCAGCGCCAGCAGCCGTTTTATGAACCGTCCTGTCATTGACATGGATTTGGCGCGTACGGCTTTGCAGGATATTATTGCCGAAAAACACAAAGTCATTACCGCCGACATCATCATCGATGCGACAGCCAAATACTACCGTATTAAAATCAGTGATATATTGGGCAAAAAACGTACGCGCAATATTGCCCGTCCGCGCCAAGTTGCCATGAGCTTGACCAAAGAGCTGACCATGCTCAGCCTGCCTTCTATCGGTGATGCCTTTGGCGGTCGCGATCACACGACTGTGATGCACGGTGTCAAAGCGGTGGCGAAACTGCGCGAAGAAGATCCCGAATTGGCGCAAGACTACGAAAAACTGCTGATTTTGATTCAGAACTGA
>CCD84_RS00015 value_md5=e5e0896a24a85d02324abe1d0e1085e7 CIGAR=CCD84_RS00015_1:1104M accepted=2 reference=GCF_005221305.fna identity=0.894 coordinates=NZ_CP039887.1:1676..2779
ATGCTGATTTTACAAGCCGATCGCGACAGTCTGCTCAAGCCGTTGCAAGCCGTTACCGGTATTGTCGAACGTCGCCATACTCTGCCGATTTTGTCTAATGTGTTGCTGGAAAGCAAAGACGGACAAACCAAACTTTTGGCAACCGACTTGGAAATCCAAATCAATACCGCCGGCCCTGAAAGTCAGGCAGGCGATTTCCGTATTACGACCAACGCTAAAAAATTCCAAGACATCCTGCGTGCTCTGCCTGACAGTGCGCTGGTGTCACTGGATTGGGCGGACAACCGTTTGACTCTGCGCGCGGGAAAATCCCGTTTTGCCCTGCAAACCTTGCCGGCTGAAGACTTTCCGTTGATGAGCGTCGGCAGCGACGTCAGTGCGACTTTCTCACTGACTCAAGAAACCTTCAAAACCATGCTTTCGCAAGTGCAATACAGCATGGCAGTTCAAGATATCCGCTATTACCTCAACGGTTTGCTGATGCAGGTTGAAGGTAATCAGCTGCGCCTTGTTGCAACCGACGGCCACCGCCTTGCTTATGCGGCCAGTCAAATTGAAGCAGAACTGCCGAAAACGGAAGTGATCCTGCCGCGTAAAACGGTATTGGAACTCTTCAAGCTGTTGAATAATCCGTCCGAGTCCATCACCGTTGAGCTTTTGGACAATCAAGTACGCTTCCAATGCAATGGCACAACCATTGTCAGCAAAGTCATCGACGGCAAGTTCCCTGACTTTAACCGCGTGATTCCTTTGGATAATGACAAGATTTTCCTCGTATCCCGTACCCAGCTTTTGGGTGCACTCGAGCGTGCCGCCATTCTTGCCAATGAAAAATTCCGCGGCGCACGACTGTTCCTGCAGCCTGGTTTGCTGAGTGTCGTATGTAGCAACAACGAGCAGGAAGAAGCGCGCGAAGAGCTGGAAATCGCTTACCAAGGCGGAGAACTCGAAGTCGGTTTCAACATCGGCTACCTGATGGATGTGTTGCGTAACATCCACTCCGACGATATGCAGCTTGCTTTCGGCGATGCCAACCGTTCAACGCTGTTTACTATGCCGAACAATCCGAACTTCAAATACATCGTAATGCCGATGCGTATTTAA
~~~~~~~~~~~~~

The header of each allelic sequence includes MD5 values, identity to reference alleles and coordinates in the query assembly. It also includes a CIGAR that describes the alignments. 


# Citation and Reproduction Instructions

### Reproduction Instructions
All data required for reproduction of the analysis were distributed in this repository under
https://github.com/zheminzhou/KleTy/tree/main/db


These includes:
* klebsiella.refsets.fas.gz - reference alleles for all pan genes (for calling new alleles)
* klebsiella.cgmlst - A list of core genes used in the dcgMLST scheme
* profile.parq - Allelic profiles of all ~70,000 genomes in parquet format, and can be read using the Pandas library (https://pandas.pydata.org/docs/reference/api/pandas.read_parquet.html). 
* HierCC.tsv.gz - A tab-delimited table consisting of HierCC results for all ~70,000 genomes
* klebsiella.species - A mapping table that specifies correlations between genomes and Klebsiella species
* virulence.fasta - reference sequences for hypervirulence genes as used in kleborate
* AMR* - reference sequences for AMR genes and mutations as used in AMRfinderPlus
* resfams_db - A mapping file that gives more accurate prediction of ESBLs
* plasmids/* - Reference sequences for INC and MOB types. 

# Naclist's modification

I additionally have set up a database for 370,000+ plasmid sequences predicted from 310,000 high quality genomes from NCBI RefSeq database by PlasT (KleTy), and all these genomes have at least one clear information about isolation date, source and geography. This database can be accessed by URLs below or HTMLs [here](https://github.com/Naclist/KleTy/tree/main/plas_test/web) for some simple visualization or statistic.
