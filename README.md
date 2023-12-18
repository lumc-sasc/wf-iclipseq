# wf-iclipseq
Nextflow pipeline for analyzing iCLIP data.

# Pipeline
TBA: workflow overview.

# Requirements
Nextflow (>23.04.4) and Singularity are required to run this pipeline. Both can be installed using Conda. To download Conda, follow this [tutorial](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html). Check that Conda is up-to-date:

`conda --version`

```plaintext
conda 23.9.0
```

Set up the proper channels, e.g.:

`conda config --add channels defaults`

`conda config --add channels bioconda`

`conda config --add channels conda-forge`

`conda config --set channel_priority strict`


You can install Nextflow using: `conda install -c bioconda nextflow`

Or a specific version: `conda install nextflow=23.10.0`

Make sure to install Singularity from conda-forge and not bioconda, as conda-forge has a more recent version:

`conda install -c conda-forge singularity=3.8.6`

Check the version of the programs to confirm you have a decently up-to-date program:

`nextflow -version`

```plaintext
      N E X T F L O W
      version 23.10.0 build 5889
      created 15-10-2023 15:07 UTC (17:07 CEST)
      cite doi:10.1038/nbt.3820
      http://nextflow.io
```

`singularity --version`

```plaintext
singularity version 3.8.6
```

# Usage
Running the pipeline can be done by using the `nextflow run main.nf` command *after mandatory parameters have been provided.*
TBA: add mandatory parameters.

# Future work
This project is currently ongoing.

# Authors
This pipeline was originally made by Amarise Silié ([@amarisesilie](https://github.com/amarisesilie)) for her MSc internship in Bioinformatics. The pipeline uses modules and subworkflows from [nf-core](https://github.com/nf-core/modules), [nf-core/rnaseq](https://github.com/nf-core/rnaseq) and [goodwright/clipseq](https://github.com/goodwright/clipseq) in addition to local modules and subworkflows.
