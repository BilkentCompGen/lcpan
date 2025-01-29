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
./lcpan -f <fasta_file> -v <vcf_file> -r <output_rgfa_file> -l <lcp_level> -t <thread_number>  [--rgfa | --gfa]
```

- `-f`: Path to the input FASTA file.
- `-v`: Path to the input VCF file.
- `-r`: Path to the output rGFA file.
- `--gfa`: Output as graphical fragment assembly.
- `--rgfa`: Output as reference gfa [default].
- `-N`: Output non-overlapping gfa.
- `-l`: LCP parsing level (integer) [default 4].
- `-t`: Thread number (integer) [default 1].

### Merging Files

The `lcpan` tool runs in parallel, hence, it generates multiple output file. Note that these files are dependent, expect the first file (as it stores the partitioned reference genome). At the end of the program execution, you can run `merge.sh lcpan.t<thread-number>.log` script that will merge all the files.

### Example

```sh
./lcpan -f genome.fasta -v variations.vcf -r output.rgfa -l 4
bash merge.sh lcpan.t4.log
```

This command constructs a variation graph for the input FASTA and VCF files, applying LCP parsing at level 4, and saves the result to `output.rgfa` and `output.rgfa.0` files. Then, you need to append content of `output.rgfa.0` file into `output.rgfa` file.

## Licence

`lcpan` is released under the BSD 3-Clause License, which allows for redistribution and use in source and binary forms, with or without modification, under certain conditions. For more detailed terms, please refer to the [license file](https://github.com/BilkentCompGen/lcpan/blob/main/LICENSE).