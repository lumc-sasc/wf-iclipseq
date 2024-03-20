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

Make a conda environment to download the packages in and activate it:
```
conda create -n nf_env
conda activate nf_env
```

You can install Nextflow using [Bioconda](https://bioconda.github.io/), after setting up the proper channels: `conda install -c bioconda nextflow`

Or a specific version: `conda install -c bioconda nextflow=23.10.0`

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
## Installation
To download the repository, you can use:

```
git clone https://github.com/lumc-sasc/wf-iclipseq.git
```
## Execution
It is mandatory to provide a **samplesheet.csv** in the `input/` folder that has 5 columns and may look like this:

```
sample,fastq_1,fastq_2,control_bam,control_bai
s_1,/path/to/fastq.gz/,,,
s_2,/path/to/R1_fastq.gz/,/path/to/R2_fastq.gz/,,
s_3,/path/to/fastq.gz/,,/path/to/input.bam,/path/to/input.bai
```
Those in **bold** are mandatory to run the pipeline.
- **sample**: sample name, make sure it does not contain periods as that may cause errors
- **fastq_1**: a demultiplexed fastq.gz file
- fastq_2: a second demultiplexed fastq.gz file, in case of paired-end reads
- control_bam: .bam file of input/control data (used for peak calling with PureCLIP/MACS2)
- control_bai: .bai file of input/control data (used for peak calling with PureCLIP)

If control_bam is provided, it is imperative that control_bai is provided as well!
**Note: This pipeline has been tested for single-end reads, and may not support paired-end reads (remains to be tested).**

By default, the reference genome and annotation are specified in the `conf/params.config`. You can use these and move on to running the pipeline or you can replace these with your own.
```
fasta                      = "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
gtf                        = "https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"
```

Finally, you can run the pipeline using the following command:
```
nextflow run main.nf
```
If you are running from a different (sub)directory, make sure you point to the folder `main.nf` is located in. It is recommended to add the `-resume` command when partially rerunning the pipeline. You may also specify the output directory with the `--outdir` parameter.

# Errors
- exit code 137 occurs when the memory limit of a process has been exceeded. For example, SortMeRNA (but also STAR) may run out of resources, especially if your sample is large (e.g. >2GB) and contains mostly rRNA (e.g. >70%) sequences you wish to filter out. You can increase these resources in the `conf/resources.config` file.
- PureCLIP may cause an error. It is usually solved by simply rerunning using the `-resume` command.

# Future work
This project is currently ongoing. Possible improvements include:

- Benchmarking
- Including a test case
- Compatibility with other container runtimes (e.g. Docker, Appcontainer, conda) to improve accessibility (currently only Singularity is supported)
- Testing and improving performance (runtime, speed, resources)
- Error testing

# Authors
This pipeline was originally made by Amarise Silié ([@amarisesilie](https://github.com/amarisesilie)) for her MSc internship in Bioinformatics. The pipeline uses modules and subworkflows from [nf-core](https://github.com/nf-core/modules), [nf-core/rnaseq](https://github.com/nf-core/rnaseq) and [goodwright/clipseq](https://github.com/goodwright/clipseq) in addition to local modules and subworkflows.
