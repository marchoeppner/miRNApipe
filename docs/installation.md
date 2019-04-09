# Installation

This pipeline uses the Nextflow framework to translate its workflow onto a compute cluster environment. 

To ease the setup process, most of the software is automatically provisioned by the pipeline. However, you will need to have a few things available on your system:

1) The [Nextflow](https://github.com/nextflow-io/nextflow/releases). Note that Nextflow requires Java (>= 1.8).

And either:

2a) The package manager [Conda](https://anaconda.org/)

or

2b) The container framework [Singularity](https://github.com/sylabs/singularity) (see [documentation](https://www.sylabs.io/docs/) )

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

More information on how to configure Nextflow, see [here](https://www.nextflow.io/docs/latest/config.html#)

## Genome references

### iGenomes

This pipeline is agnostic with respect to the genome you wish to use as a basis for your analysis. We find it easiest to fall back on the iGenomes reference pack as it includes a wide range of pre-formatted references of typically used model systems. 

To get your own iGenomes data, you can refer to [this](https://github.com/ewels/AWS-iGenomes) instruction. 

iGenomes data is organized in a hierarchical folder structure, if which "references" is the top level. You can point your cluster config to this location and the pipeline will be able to pick up the desired genome automatically (assuming it was indeed downloaded). 

For this, please include the following line in your custom config file (see above):

`igenomes_base=/path/to/igenomes/references`

### Your custom genome

If you are using a species not included with iGenomes (or prefer another assembly version), you can of course use that too. Please see the usage instructions on how to specify custom genomes.

