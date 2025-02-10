# `LCPan` (Variation Graph Çonstruction Using Locally Consistent Parsing)  
![GitHub release (latest SemVer)](https://img.shields.io/github/v/release/BilkentCompGen/lcpan)
![GitHub last commit](https://img.shields.io/github/last-commit/BilkentCompGen/lcpan)
![GitHub](https://img.shields.io/github/license/BilkentCompGen/lcpan)

This repository provides an implementation of a tool for constructing variation graphs using Locally Consistent Parsing (LCP) tool. It processes genomic sequences from FASTA files, integrates variations from VCF files, and generates variation graphs in rGFA/GFA format for genome assembly and analysis.

## Features

- **Efficient Variation Graph Construction**: Utilizes LCP cores to represent sequences and variations.
- **Variation Representation**: Creates bubble structures in the graph to represent sequence differences.
- **rGFA/GFA Output**: Generates graphs in rGFA/GFA format, suitable for graph-based genome analysis.

## Installation

Clone this repository and use the provided Makefile to build the project.  

```sh
git clone https://github.com/BilkentCompGen/lcpan.git
cd lcpan

# install lcptools
make install

# compile lcpan
make
```

You can run `make clean` command to remove cleanup binaries and the executable.

## Usage

Run the tool using the following command-line options:

```sh
./lcpan -r <fasta_file> -v <vcf_file> -o <output_rgfa_file> -l <lcp_level> -t <thread_number> -p <prefix> [--rgfa | --gfa]
```

- `-r | --ref`: Path to the input FASTA file.
- `-v | --vcf`: Path to the input VCF file.
- `-o | --output`: Path to the output rGFA file.
- `-p | --prefix`: Prefix for the log file [default lcpan].
- `-s`: Output overlapping gfa.
- `-l | --level`: LCP parsing level (integer) [default 5].
- `-t | --thread`: Thread number (integer) [default 1].
- `--gfa`: Output as graphical fragment assembly.
- `--rgfa`: Output as reference gfa [default].
- `--tload-factor`: Ho much workload is assigned per thread relative to the pool size [default 2].
- `-v`: Verbose [default false].

### Merging Files

The `lcpan` tool runs in parallel, hence, it generates multiple output file. Note that these files are dependent, expect the first file (as it stores the partitioned reference genome). At the end of the program execution, you can run `merge.sh lcpan.log` script that will merge all the files.

### Example

```sh
./lcpan -r genome.fasta -v variations.vcf -o output.rgfa -l 3
bash merge.sh lcpan.log
```

This command constructs a variation graph for the input FASTA and VCF files, applying LCP parsing at level 4 using single thread, and saves the result to `output.rgfa` and `output.rgfa.0` files. Then, you need to append content of `output.rgfa.0` file into `output.rgfa` file.

## Licence

`lcpan` is released under the BSD 3-Clause License, which allows for redistribution and use in source and binary forms, with or without modification, under certain conditions. For more detailed terms, please refer to the [license file](https://github.com/BilkentCompGen/lcpan/blob/main/LICENSE).