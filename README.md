# `VG_LCP` (Variation Graph Construction Using Locally Consistent Parsing)  

This repository provides an implementation of a tool for constructing variation graphs using Locally Consistent Parsing (LCP) tool. It processes genomic sequences from FASTA files, integrates variations from VCF files, and generates variation graphs in rGFA/GFA format for genome assembly and analysis.

## Features

- **Efficient Variation Graph Construction**: Utilizes LCP cores to represent sequences and variations.
- **Variation Representation**: Creates bubble structures in the graph to represent sequence differences.
- **rGFA/GFA Output**: Generates graphs in rGFA/GFA format, suitable for graph-based genome analysis.

## Installation

Clone this repository and use the provided Makefile to build the project.  

```sh
git clone https://github.com/BilkentCompGen/lcp_vg.git
cd lcp_vg

# install lcptools
make install

# compile lcp_vg
make
```

You can run `make clean` command to remove cleanup binaries and the executable.

## Usage

Run the tool using the following command-line options:

```sh
./lcp_vg -f <fasta_file> -v <vcf_file> -r <output_rgfa_file> -l <lcp_level> [--rgfa | --gfa]
```

### Merging Files

The `lcp_vg` tool runs in parallel, hence, it generates multiple output file. Note that these files are dependent, expect the first file (as it stores the partitioned reference genome). At the end of the program execution, you can run `merge.sh` script that will merge all the files.

- `-f`: Path to the input FASTA file.
- `-v`: Path to the input VCF file.
- `-r`: Path to the output rGFA file.
- `--gfa`: Output as graphical fragment assembly.
- `--rgfa`: Output as reference gfa [default].
- `-N`: Output non-overlapping gfa.
- `-l`: LCP parsing level (integer) [default 4].
- `-t`: Thread number (integer) [default 1].

### Example

```sh
./lcp_vg -f genome.fasta -v variations.vcf -r output.rgfa -l 4
```

This command constructs a variation graph for the input FASTA and VCF files, applying LCP parsing at level 4, and saves the result to `output.rgfa` and `output.rgfa.0` files. Then, you need to append content of `output.rgfa.0` file into `output.rgfa` file.

## Licence

`lcp_vg` is released under the BSD 3-Clause License, which allows for redistribution and use in source and binary forms, with or without modification, under certain conditions. For more detailed terms, please refer to the [license file](https://github.com/BilkentCompGen/lcp_vg/blob/main/LICENSE).