# gencomvariable
**Technical Documentation**

Gencomvariable is a command-line tool designed to assemble viral genomes based on a reference.

### Overview
Gencomvariable accepts fastq.gz files generated by the Illumina platform and those generated by Oxford Nanopore Technologies (ONT). The script cleans and performs quality analyses on the sequences before mapping them to a reference that can be supplied directly as a fasta file or through an accession number. The output of the analysis is a consensus genome in fasta format. Additionally, it provides information on various parameters through tables and images, such as the coverage obtained for the assembled genome.

### Installation and Requirements
The script has no software requirements. It is simply executed using the absolute path, as shown below:

#/your/path/#gencomvariable.sh


It also has a help (-h) option to obtain information on the script's requirements, as follows:

#/your/path/#gencomvariable.sh -h


This option will display the following information:

#/your/path/#gencomvariable.sh no option added to this script, please read the following help or use gencomvariable.sh -h:
    a) # use -a to specify Accession reference number to use for mapping.
    m) # use -m FALSE to avoid mapping to a reference and instead use a virus-tailored analysis (currently PTV, with more to come).
    d) # Specify the directory where R1 and R2 reads are located.
    v) # Virus for special analysis. Currently supports HIV, PTV, IMRA Influenza, Illumina and MinION, MPXV, with custom options available.
    c) # Apply a provided reference in fasta format to reads; this option must be used with options -v custom and -m FALSE.
    t) # Provide a TSV file for ONT Influenza sequencing. The first column must be the sample, and the second must be the barcode.
    h | *) # Display help for this pipeline.


### Input Files
To assemble most genomes, the script requires the following files:

1. **Reference**: The reference sequence can be provided in two ways. The first option requires only the accession number of the desired reference genome. In this case, the script uses the ESearch command to find the accession number and download the sequence from the NCBI database. The second option allows a custom analysis by using a manually supplied reference/genome in fasta format.
2. **fastq.gz or fastq**: The script can use paired-end fastq.gz files (R1 and R2) generated through Illumina platforms, where each sequence has a unique identifier. Data generated by the ONT platform must be separated by barcode in directories (the result of demultiplexing). It also requires a comma-separated (csv) or tab-separated (tsv) table where the first column contains the sequence ID (sample ID), and the second column contains the barcode assigned to the sample during DNA library preparation.

### Tools Used in Data Processing
- **fastp**: Used for preprocessing data and quality control of FASTQ files.
- **minimap2**: Used to align fastp results by mapping sequences to the provided reference, resulting in Sequence Alignment Map (SAM) files, which is a text-based format for storing biological sequences aligned to a reference sequence.
- **samtools**: This tool sorts, converts to BAM (Binary Alignment Map), and indexes BAM files, which are binary-compressed versions of SAM files.
- **bedtools**: bedtools utilities handle common genomic tasks, such as finding overlaps and calculating coverage of the assembled genome. It also allows merging, complementing, or intersecting genomic intervals from different files.
- **bcftools**: Used for variant calling and manipulating VCF (Variant Call Format) and BCF (Binary Variant Call Format) files. VCF is a text format for storing variants relative to a reference genome.
  This tool also enables creating a consensus sequence for each analyzed sample, where the sequence incorporates identified variants for the genome. This is achieved using the consensus command.
- **plotgroupedcoveragev2.py**: A Python script that generates coverage plots of analyzed viral sequences or genomes.

### *Assembly of the Eight Genome Segments of the Influenza Virus*
For the influenza virus, the analysis pipeline uses a method developed by the Centers for Disease Control and Prevention (CDC) called IRMA (Iterative Refinement Meta-Assembler).

- **IRMA**: Designed for robust assembly, variant calling, and haplotype estimation (also known as "phasing") of highly variable RNA viruses. IRMA currently has modules for influenza, ebolavirus, and coronavirus. IRMA is free to use and can parallelize calculations for both cluster computing and multicore single-computer setups.
- **coverage4allgenes.R**: This R script combines all separate coverage plots for each influenza virus segment and saves them as pdf files.

### Output Files
At the end of the analysis, depending on the virus, the script creates a directory with results generated at each stage. The directory may contain individual information per sample ID and grouped documents, either as sequence lists in fasta format or as images.

Among the generated directories and files are the following:
- **Directories for each sample**: The results of each sample's analysis are stored in a directory with the corresponding ID.
- **Consensus sequence**: The analyzed sample's sequence in fasta format.
- **fastqreport.html**: A report generated by fastp on the quality of the sequences in the used fastq files.
- **Coverage files directory**: Coverage obtained for all analyzed samples is stored in tables in "csv" format and images in "png" format.
