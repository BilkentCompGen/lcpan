# `LCPan` (Variation Graph Construction Using Locally Consistent Parsing)  
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
git clone --recursive --depth 1 https://github.com/BilkentCompGen/lcpan.git
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
./lcpan [PROGRAM] [OPTIONS]
```

Program:
- `-vg`: Constructs a variation graph using a reference genome and VCF. In this mode, the initial partitioning is done with LCP, and each segment is further divided into sub-segments if variations are present.
- `-vgx`: Constructs an expanded variation graph using a reference genome and VCF. In this mode, each variation is represented by an alternative arc, which connects the latest non-overlapping LCP core to the first LCP core afterward.

Options:

- `-r | --ref`: Path to the input FASTA file.
- `-v | --vcf`: Path to the input VCF file.
- `-p | --prefix`: Prefix for the log and output file [default lcpan].
- `-s | --no-verlap`: Output overlapping gfa.
- `-l | --level`: LCP parsing level (integer) [default 5].
- `-t | --thread`: Thread number (integer) [default 1].
- `-v | --verbose`: Verbose [default false].
- `--gfa`: Output as graphical fragment assembly.
- `--rgfa`: Output as reference gfa [default].
- `--skip-masked`: Skit masked (N) characters. In this mode, segments will contain only nucleotides.
- `--tload-factor`: How much workload is assigned per thread relative to the pool size [default 2].

### Merging Files

The `lcpan` tool runs in parallel, hence, it generates multiple output file. Note that these files are dependent, expect the first file (as it stores the partitioned reference genome). At the end of the program execution, you can run `lcpan-merge.sh lcpan.log` script that will merge all the files.

### Example 1

```sh
./lcpan -vg -r genome.fasta -v variations.vcf -p output -l 4
bash lcpan-merge.sh output.log
```

This command constructs a variation graph for the input FASTA and VCF files, applying LCP parsing at level 4 using single thread, and saves the result to various files. Then, you need to merge the files (which will be done by `lcpan-merge.sh` script).

## Citation
If you use LCPan in your work, please cite:
- LCPan: efficient variation graph construction using Locally Consistent Parsing. Akmuhammet Ashyralyyev, Zülal Bingöl, Begüm Filiz Öz, Kaiyuan Zhu, Salem Malikic, Uzi Vishkin, S. Cenk Sahinalp, Can Alkan. [arXiv: 2511.12205](https://doi.org/10.48550/arXiv.2511.12205), 2025.
 
## License

`lcpan` is released under the BSD 3-Clause License, which allows for redistribution and use in source and binary forms, with or without modification, under certain conditions. For more detailed terms, please refer to the [license file](https://github.com/BilkentCompGen/lcpan/blob/main/LICENSE).
