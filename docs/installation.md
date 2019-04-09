# Installation

This pipeline uses the Nextflow framework to translate its workflow onto a compute cluster environment. 

To ease the setup process, most of the software is automatically provisioned by the pipeline. However, you will need to have a few things available on your system:

* The [Nextflow](https://github.com/nextflow-io/nextflow/releases) framework. Note that Nextflow requires Java (>= 1.8).

And either:

* The package manager [Conda](https://anaconda.org/)

or

* The container framework [Singularity](https://github.com/sylabs/singularity) (see [documentation](https://www.sylabs.io/docs/) )

## Cluster config files

Nextflow uses config files to be able to e.g. talk to a specific cluster. Some of these are mandatory, some are site-specific and some are optional. A set of config files can be combined into a "profile", which defines a specific compute environment (e.g. your cluster setup). 

To define a new profile, create a new entry in the profile section of `nextflow.config`.

A tpyical setup will require, at minimum, these files:

`base.config` (defines basic behavior of the pipeline)

`your_cluster.config` (defines settings relevant to your cluster)

Please see the included template file [template.config](../conf/template.config) what such a file should contain. 

Critical information includes:

- Your resource manager

- Whether to use Conda or Singularity for software provisioning.

An example of how to use this pipeline with Singularity on a Slurm cluster is given in the profile `ccga`, where the file [rzcluster.config](../conf/rzcluster.config) include the site-specific settings as well the information required to use Singularity. 

More information on how to configure Nextflow, see [here](https://www.nextflow.io/docs/latest/config.html#)

## Genome references

### Genome index

Since this pipelines uses mirTop as part of its processing chain, users are limited to genome assemblies supported in miRBase. 

To speed up the pipelines, please consider downloading your genome assembly of interest to a locatation that is visible on your cluster,
create an index folder for the STAR aligner and provide that location as command line option to the pipeline using `--star_index`. 

Genomes and the ftp link to download them as follows:

Human: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

Mouse: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.4_GRCm38.p2/GCA_000001635.4_GRCm38.p2_genomic.fna.gz

C. elegans: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/002/985/GCA_000002985.3_WBcel235/GCA_000002985.3_WBcel235_genomic.fna.gz

D. melanogaster: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/215/GCA_000001215.4_Release_6_plus_ISO1_MT/GCA_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz

If you do not want to do that, any maybe only want to use this pipeline once, you can skip this step and simply specify the species name as `--genome human` (or your species of interest from the list above). The pipeline will then download the assembly and turn it into a STAR index. However, this step takes a long time.

