# Usage information

## Basic execution

Assuming you followed the [setup instructions](installation.md), executing the pipeline using the iGenomes reference for human works as follows:

`nextflow run main.nf --reads 'path/to/*_R1_001.fastq.gz' --genome GRCh37 -profile ccga --gtf mirnas.gtf`

This would start the pipeline using the CCGA profile (which uses Singularity for software provisioning), and will align all reads matching the specified pattern against the human reference GRCh37 included with the iGenomes data. 

Matching miRNA coordinates are provided as a GTF file. 

## Using a custom genome

If you do not wish to use one of the pre-built iGenomes reference assemblies, you have two options:

- Provide an existing STAR index for alignment (STAR version <= 2.6) using the `--star_index /path/to/index` option

or

- Provide the genome assembly as FASTA file using `--fasta /path/to/genome.fasta`; a STAR index will be built for you (but this takes some time).


