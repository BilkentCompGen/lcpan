# `VG_LCP` (Variation Graph Construction Using Locally Consistent Parsing)  

This repository provides an implementation of a tool for constructing variation graphs using Locally Consistent Parsing (LCP) tool. It processes genomic sequences from FASTA files, integrates variations from VCF files, and generates variation graphs in rGFA format for genome assembly and analysis.

## Features

- **Efficient Variation Graph Construction**: Utilizes LCP cores to represent sequences and variations.
- **Variation Representation**: Creates bubble structures in the graph to represent sequence differences.
- **rGFA Output**: Generates graphs in rGFA format, suitable for graph-based genome analysis.

## Installation

Clone this repository and use the provided Makefile to build the project.  

1. Clone the repository:
    ```sh
    git clone https://github.com/EgeSrvn/lcp_vg.git
    cd vg-lcp
    ```

2. Install the LCP Tool Repository
    ```sh
    make install
    ```

3. Build the project:
    ```sh
    make compile
    ```

3. Clean the build (if needed):
    ```sh
    make clean
    ```

4. Uninstall the LCP Tool Repository
    ```sh
    make uninstall
    ```

## Usage

Run the tool using the following command-line options:

```sh
./vg_lcp -f <fasta_file> -v <vcf_file> -r <output_rgfa_file> -l <lcp_level>
```

- `-f`: Path to the input FASTA file.
- `-v`: Path to the input VCF file.
- `-r`: Path to the output rGFA file.
- `-l`: LCP parsing level (integer).

### Example

```sh
./vg_lcp -f genome.fasta -v variations.vcf -r output.rgfa -l 4
```

This command constructs a variation graph for the input FASTA and VCF files, applying LCP parsing at level 4, and saves the result to `output.rgfa`.

