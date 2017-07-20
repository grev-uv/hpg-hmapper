# HPG-Hmapper

A parallel software tool for analyzing DNA hydroximethylation data.

Providing two BAM datasets, one generated through bisulphite sequencing (BS-seq) with 5mC information and other generated through TET1 traslocation BS-sequencing (TAB-seq) with 5hmC data, HPG-Hmapper detects and maps the methylated and hydroximethylated regionsm, creating per-chromosome CSV files with a list of all the detected cytosines featuring relevant information and its methylation, hydroximethylation, mutation and coverage information.



# Usage

Starting with two FASTQ files, one with the 5mC sequences and the other with the 5hmC sequences, the first step is aligning them to the reference genome and performing the methylation call. Our recommended software tool for this step is [HPG-Methyl](https://github.com/grev-uv/hpg-methyl). However, other aligners and methylation callers compatible with the Bismark output format can be used.

When the alignment and calling process is completed, call HPG-Hmapper to obtain the methylation and hydroximethylation maps with the command:

```
hpg-hmapper -mc <path_to_bs_seq.bam> -hmc <path_to_tab_seq.bam> -o <out_path>
```

Where `<path_to_bs_seq.bam>` is the BAM file containing the 5mC sequences, `<path_to_tab_seq.bam>` is the BAM file containing the 5 hmC sequences and `<out_path>` is the output directory where the CSV files containing the methylation maps will be stored.

The software is configured to have the best performance in most scenarios, however it can be hand-tuned using the following optional command line parameters.

| Parameter | Description | Possible values | Default |
|:---------:|:-----------:|:---------------:|:-------:|
| `--num-threads` | The number of CPU threads to use. | 1 (single-threaded mode) or 3 to 24 (worker mode) | The number of CPU cores in the machine |
| `--memory` | The maximum RAM memory in MB or GB to use | `X.XXM` or `X.XXG` | 80% of the system memory |
| `--output-format` | Output format for the execution statistics | `text` or `csv` | `text` |
| `--batch-size` | Aligned read batch size for the application queues | > 10000 | 1000000 |
| `--csv-delimiter` | Delimiter for the output CSV columns | A single character (\t and \n are accepted) | Space |
| `--csv-record-delimiter` | Delimiter for the output CSV rows | A single character (\t and \n are accepted) | Line break (`\n`) |

For example, to limit memory consumption to 20 GB and create the methylation maps using commas and semicolons as separators, the optional command line parameters would be:

```
--memory 20.0G --csv-delimiter , --csv-record-delimiter ;
```

# Output format

The generated methylation maps will be stored as CSV files, one per chromosome and per strand (forward and reverse), with four columns using the following format. All values are zero-indexed:

| Name | Description | Possible values |
|:----:|:-----------:|:---------------:|
| Position | Genome position in the chromosome | 0 ~ chromosome length minus one |
| #C | Number of non-methylated cytosines in the position | ≥ 0 |
| #non-C | Number of overlapping non-cytosines in the position | ≥ 0 |
| #5mC | Number of 5-methylated cytosines in the position | ≥ 0 |
| #5hmC | Number of 5-hydroximethylated cytosines in the position | ≥ 0 |

An example output, using tabulators and line breaks as separators would be:

`methylation_map_forward_1.csv`
```
10723   4     249   0     5
10724   191   9     28    26
10726   207   0     109   137
10727   89    0     59    246
```



# System requirements

HPG-Hmapper has been extensively tested using synthetic and real datasets. The package should work correctly on the following set-up.

#### CPU

A **64 bit Intel CPU compatible with SSE4.2** is required. Tests have been conducted with an Intel Xeon E5-2650V4 CPU, using from 1 to 12 threads. 

#### RAM Memory

The peak memory usage ranged between 1 GB and 8 GB for 40 million sequence datasets, with sequence lengths ranging from 75 to 800nt. However, a minimum of **16 GB** is recommended.

When very large datasets (with several hundred million sequences) are being processed, the memory usage for the optimum performance can grow up to tens of GB (or even further). The performance / memory consumption ratio can be controlled through optional command line parameters.



# Building

In the release section of the repository there is a binary package for Linux x86_64 systems.If you are interested in modifying, extending or debugging the software, the following instructions show how to build HPG-Hmapper. To be able to build the software, these software packages must be installed in the system:

| Package | Ubuntu / Debian | Red Hat / Fedora / Centos |
|:-------:|:---------------:|:-------------------------:|
| GNU C Toolchain     | build-essential | make / automake / gcc |
| SConstruct    | scons      | scons                |
| ZLib    | zlib1g-dev      | zlib-devel                |
| Curl    | libcurl4-gnutls-dev | libcurl-devel         |
| libxml  | libxml2-dev     | libxml2-devel             |
| ncurses | libncurses5-dev | ncurses-devel             |
| GNU GSL | libgsl0-dev     | gsl-devel                 |
| check   | check           | check-devel               |

The Samtools software package is not required but it is recommended to transform and view the input files.

When all the required software packages are installed, to build a release executable of HPG-Hmapper use the command:

```
scons
```

To build a debug version of HPG-Hmapper, with extended debug messages and optimization disabled use the command:

```
scons debug=1
```

Besides compiling from the command line, the project contains Visual Studio code macros to build, debug and test the application from the IDE.

# Issues

If you find any bugs, issues, want a specific feature added or need help, feel free to add an issue or extend an existing one. Pull requests are welcome.



# License

HPG-Hmapper is free software and licensed under the GNU General Public License version 2. Check the [COPYING](/COPYING) file for more information.



# Contact

Contact any of the following developers for any enquiry:

* Juanma Orduña ([juan.orduna@uv.es](mailto:juan.orduna@uv.es)).
* Mariano Pérez ([mariano.perez@uv.es](mailto:mariano.perez@uv.es)).
* César González ([cesar.gonzalez-segura@uv.es](mailto:cesar.gonzalez-segura@uv.es)).