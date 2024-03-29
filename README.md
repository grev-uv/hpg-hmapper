# HPG-Hmapper

A parallel software tool for analyzing DNA hydroximethylation data.

If you want to use this tool just now, there is an executable file for Linux x86_64 systems. This compressed file is available at [releases page](../../releases).

Providing two BAM datasets, one generated through bisulphite sequencing (BS-seq) with 5mC information and other generated through TET1 traslocation BS-sequencing (TAB-seq) with 5hmC data, HPG-Hmapper detects and maps the methylated and hydroxymethylated regions, creating per-chromosome CSV files with a list of all the detected cytosines featuring relevant information and its methylation, hydroxymethylation, mutation and coverage information.

# Usage

Starting with two FASTQ files, one with the 5mC sequences and the other with the 5hmC sequences, the first step is aligning them to the reference genome and performing the methylation call. Our recommended software tool for this step is [HPG-Methyl](https://github.com/grev-uv/hpg-methyl). However, other aligners and methylation callers compatible with the Bismark output format can be used.

When the alignment and calling process is completed, call HPG-Hmapper to obtain the methylation and hydroxymethylation maps with the command:

```
$ hpg-hmapper -mc <path_to_bs_seq.bam> -hmc <path_to_tab_seq.bam> -o <out_path>
```

Where `<path_to_bs_seq.bam>` is the BAM file containing the 5mC sequences, `<path_to_tab_seq.bam>` is the BAM file containing the 5 hmC sequences and `<out_path>` is the output directory where the CSV files containing the methylation maps will be stored.

When only methylation map is required, an empty hydroxymethylated bam file is needed. It is possible to download it [here](hmapper/bin/)

The software is configured to have the best performance in most scenarios, however it can be hand-tuned using the following optional command line parameters.

| Parameter | Description | Possible values | Default |
|:---------:|:-----------:|:---------------:|:-------:|
| `--num-threads` | The number of CPU threads to use. | 1 (single-threaded mode) or 3 to 24 (worker mode) | The number of CPU cores in the machine |
| `--memory` | The maximum RAM memory in MB or GB to use | `X.XXM` or `X.XXG` | 80% of the system memory |
| `--output-format` | Output format for the execution statistics | `text` or `csv` | `text` |
| `--batch-size` | Aligned read batch size for the application queues | > 10000 | 1000000 |
| `--csv-delimiter` | Delimiter for the output CSV columns | A single character (\t and \n are accepted) | Space |
| `--csv-record-delimiter` | Delimiter for the output CSV rows | A single character (\t and \n are accepted) | Line break (`\n`) |
| `--quality` | Minimum quality cutoff for the input reads | 0~254 | 20 |
| `-c` `--coverage` | Minimum coverage to validate the methylated position | 1-1000 | 15|
| `-i` `--bwt-index` | Path to directory where the BWT index is stored | /path/to/bwt-index/ | human GRCh37.68 genome |

For example, to limit memory consumption to 20 GB and create the methylation maps using commas and semicolons as separators, the optional command line parameters would be:

```
--memory 20.0G --csv-delimiter , --csv-record-delimiter ;
```

Methylation and coverage graphs can be created with the following command:

```
$ hpg-hmapper-graph-tool <csv_meth_map> <start_position> <end_position> <samples> <chromosome>
```

Two publishing-ready EPS plots will be generated with the data relative to the selected range.

Example datasets can be created using the scripts available in the `datasets` directory. Check the directory readme for instructions on how to generate the datasets.

# Output format

The generated methylation maps will be stored as CSV files, one per chromosome and per strand (forward, reverse and mix), with seven columns using the following format. All values are zero-indexed:

| Name | Description | Possible values |
|:----:|:-----------:|:---------------:|
| Position | Genome position in the chromosome | 0 ~ chromosome length minus one |
| #C | Number of non-methylated cytosines in the position | ≥ 0 |
| #non-C | Number of overlapping non-cytosines, from 5mC files, in the position | ≥ 0 |
| #5mC | Number of 5-methylated cytosines in the position | ≥ 0 |
| #Ch | Number of non-hydroxymethylated cytosines in the position | ≥ 0 |
| #non-Ch | Number of overlapping non-cytosines, from 5hmC files, in the position | ≥ 0 |
| #5hmC | Number of 5-hydroxymethylated cytosines in the position | ≥ 0 |

An example output, using tabulators and line breaks as separators would be:

`methylation_map_mix_13.csv`
```
19778560 0 13 71 1 14 28
19778561 71 13 0 27 14 2
19778565 14 0 74 15 0 31
19778566 75 13 2 27 14 3
```



# System requirements

HPG-Hmapper has been extensively tested using synthetic and real datasets. The package should work correctly on the following set-up.

#### CPU

A **64 bit Intel CPU compatible with SSE4.2** is required. Tests have been conducted with an Intel Xeon E5-2650V4 CPU, using from 1 to 12 threads. 

#### RAM Memory

The peak memory usage ranged between 1 GB and 8 GB for 40 million sequence datasets, with sequence lengths ranging from 75 to 800nt. However, a minimum of **16 GB** is required.

If a memory size smaller than 16 GB is used, then some executions may be aborted with a SIGKILL signal due to memory overflow. These situations can be minimized using the --write-context option if the alignment has been carried out with [hpg-methyl](https://github.com/grev-uv/hpg-methyl).

When very large datasets (with several hundred million sequences) are being processed, the memory usage for the optimum performance can grow up to tens of GB (or even further). The performance / memory consumption ratio can be controlled through optional command line parameters.



# Building

In the [release section](../../releases) of the repository there is a binary package for Linux x86_64 systems. If you are interested in modifying, extending or debugging the software, the following instructions show how to build HPG-Hmapper. To be able to build the software, these software packages must be installed in the system:

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

If the distribution and release of the OS is Ubuntu 14.04, 18.04 or 20.04, please, take a look at Issues section before compile.

The Samtools software package is not required but it is recommended to transform and view the input files. The GNUplot package is required to use the graph generation tool.

To build HPG-Hmapper, the graph generator tool and the tools to create the synthetic datasets use the command:

```
$ build-all
```

When all the required software packages are installed, to build a release executable of only HPG-Hmapper use the commands:

```
$ cd hmapper/
$ scons
```

To build a debug version of HPG-Hmapper, with extended debug messages and optimization disabled use the commands:

```
$ cd hmapper/
$ scons debug=1
```

Besides compiling from the command line, the project contains Visual Studio code macros to build, debug and test the application from the IDE.

# Issues

Firstly, check the distibution and release of the Operating System:

```
$ lsb_release -a
```

If the OS is Ubuntu 14.04, a change is needed to compile properly:

In the file [hmapper/lib/common-libs/containers/test/SConscript](hmapper/lib/common-libs/containers/test/SConscript#L8), the line 8 must be changed deleting the 'subunit' element

```
8 LIBS = ['check', 'curl', 'm', 'z', 'rt'],
```

If the OS is Ubuntu 18.04 or 20.04, a change is needed to compile properly:

In the file [hmapper/lib/common-libs/commons/config/libconfig.c](hmapper/lib/common-libs/commons/config/libconfig.c#L37), the line 37 must be commented

```
37 //#include <xlocale.h>
```

Finally, it can compile following the steps showed at Building section.

If you find any other bugs, issues, want a specific feature added or need help, feel free to add an issue or extend an existing one. Pull requests are welcome.


# License

HPG-Hmapper is free software and licensed under the GNU General Public License version 2. Check the [COPYING](/COPYING) file for more information.



# Contact

Contact any of the following developers for any enquiry:

* Juanma Orduña ([juan.orduna@uv.es](mailto:juan.orduna@uv.es)).
* Mariano Pérez ([mariano.perez@uv.es](mailto:mariano.perez@uv.es)).
* César González ([cesar.gonzalez-segura@uv.es](mailto:cesar.gonzalez-segura@uv.es)).
