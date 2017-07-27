# Generating synthetic datasets

HPG-Hmapper has been tested using synthetic datasets, with 4 million reads and read lengths ranging from 75 nt to 800 nt. To generate the same datasets, use the following procedure.

## Required software and data

* [HPG-Methyl](https://github.com/grev-uv/hpg-methyl) version 3.1 or greater.
* [DWGSIM](https://github.com/nh13/DWGSIM) whole-genome simulator for NGS.
* A reference genome to use as the seed (for example the [GRCh37](ftp://ftp.ensembl.org/pub/release-75//fasta/homo_sapiens/dna/) reference genome).

## Procedure

First, build all the simulation software. From the repository main directory run the command:

```
./build-all
```

When the build process has completed, edit the following fields from the `generate_datasets` script in this directory:

```
FASTQ_OUT_DIR=""
BAM_OUT_DIR=""

(···)

GENOME_AGT_FASTA=""
GENOME_ACT_FASTA=""
GENOME_INDEX_DIR=""

(···)

HPG_METHYL_PATH=""
DWGSIM_PATH=""
```

The next step involves creating a Burrows-Wheeler index for the reference genome using HPG-Methyl. Refeer to its documentation to complete this step.

When the index has been created, generate the datasets running the script in this directory:

```
./generate_datasets
```

The process will yield as output the FASTQ and BAM aligned files for 5mC and 5hmC samples.