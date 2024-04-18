# wf-iclipseq
wf-iclipseq is a bioinformatic pipeline developed in Nextflow for the analysis of iCLIP data. It takes demultiplexed single-end reads (.fastq.gz), a gene annotation file, and a reference genome file to perform pre-processing (incl. quality control (QC), adapter trimming, rRNA filtering), post-processing (incl. genome alignment and crosslink extraction) and downstream analysis (gene annotation, motif analysis). The gene annotation results can be further visualized, taking the bash and R scripts in the [scripts/](scripts/) folder as reference. 


![Alt text](figures/pipeline_workflow.png?raw=true "Pipeline design")


1. Removing whitespaces and special characters in read IDs (Bash)
2. QC on raw reads ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. Barcode extraction ([UMI-Tools](https://github.com/CGATOxford/UMI-tools))
4. Adapter and quality trimming ([Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
5. Ribosomal RNA filtering ([SortMeRNA](https://github.com/sortmerna/sortmerna) or [Ribodetector](https://github.com/hzi-bifo/RiboDetector)*)
6. Quality control on non-RNA reads ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
7. Genome alignment ([STAR](https://github.com/alexdobin/STAR))
8. Quality control on unmapped reads ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
9. Deduplication ([UMI-Tools](https://github.com/CGATOxford/UMI-tools))
10. Alignment sorting and indexing ([SAMtools](https://sourceforge.net/projects/samtools/files/samtools/))
11. Quality control on deduplicated reads ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
12. Extract crosslink events ([BEDTools](https://github.com/arq5x/bedtools2/))
13. BED to BigWig coverage tracks ([bedGraphToBigWig](https://hgdownload.soe.ucsc.edu/admin/exe/))
14. Peak calling ([PureCLIP](https://github.com/skrakau/PureCLIP) and/or [MACS2](https://github.com/macs3-project/MACS))
15. Gene annotation (Bash & [BEDTools](https://github.com/arq5x/bedtools2/) and/or [HOMER](http://homer.ucsd.edu/homer/download.html))
16. Custom scripts for visualization of results ([R](https://www.r-project.org/) and [various packages](https://github.com/lumc-sasc/wf-iclipseq/blob/main/CITATION.md#r-packages))
17. Motif detection ([STREME](https://meme-suite.org/meme/doc/download.html))
18. Summarizing results up to the last quality control step ([MultiQC](https://multiqc.info/))

Steps in the pipeline (with the exception of STAR alignment, for the time being) can be skipped using the [conf/params.config](conf/params.config) file. The pipeline contains multiple quality control steps which can be used to closely monitor the data before peak calling. By default, quality control before and after trimming is enabled, while the others are disabled.

*not included in the pipeline yet, but the pipeline does contain the nextflow script.

## Default run
By default, the pipeline performs the following things:
1. QC on raw reads by FastQC
2. Quality and adapter trimming by Trim Galore!
3. Genome indexing by STAR
4. Genome alignment by STAR (allowing multimappers, with `--outFilterMultimapNmax 20` and doing chimeric detection)
5. Alignment sorting and indexing by SAMtools
6. Deduplication by UMI-Tools
7. Crosslink extraction
8. Peak calling by PureCLIP
9. Resizing PureCLIP binding sites for gene annotation by 8bp (change this to match your RBPs binding footprint, can also be deactivated)
10. Gene annotation by Bash & BEDTools, and HOMER
11. Motif detection by STREME
12. Report by MultiQC

### Output created by default run
The output of a default run creates the following tree:
```bash
├── 01_fastqc
├── 02_trimgalore
│   └── fastqc
├── 04_star
│   ├── log
│   ├── samtools_stats
│   └── umitools
│       ├── bai
│       ├── bam
│       └── fastqc
├── 05_extract_crosslinks
│   └── 06_crosslink
│       ├── merged
│       ├── neg
│       └── pos
├── 06_peak_calling
│   └── pureclip
├── 07_downstream_analysis
│   ├── 00_sites_with_peakID
│   ├── 01_bedtools_annotations
│   ├── 01_homer_annotations
│   └── 02_streme_motifs
│       └── bedtools_resized
├── 08_multiqc
│   ├── multiqc_data
│   └── multiqc_plots
│       ├── pdf
│       ├── png
│       └── svg
└── pipeline_info
```

If rRNA filtering is enabled, there will be a folder named `03_sortmerna/`.

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

You can install Nextflow using [Bioconda](https://bioconda.github.io/), after setting up the proper channels: `conda install -c conda-forge -c bioconda nextflow`

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
s_3,/path/to/fastq.gz/,,/path/to/input.bam,/path/to/input.bam.bai
```
Those in **bold** are mandatory to run the pipeline.
- **sample**: sample name, make sure it does not contain periods as that may cause errors
- **fastq_1**: a demultiplexed fastq.gz file
- fastq_2: a second demultiplexed fastq.gz file, in case of paired-end reads
- control_bam: .bam file of input/control data (used for peak calling with PureCLIP/MACS2)
- control_bai: .bai file of input/control data (used for peak calling with PureCLIP)

If control_bam is provided, it is imperative that control_bai is provided as well!
**Note: This pipeline has been tested for single-end reads, and may not support paired-end reads (remains to be tested).**

By default, the reference genome and annotation are specified in the [conf/params.config](conf/params.config) file. You can use these and move on to running the pipeline or you can replace these with your own.
```
fasta                      = "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
gtf                        = "https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"
```

Finally, you can run the pipeline using the following command:
```
nextflow run main.nf -profile singularity
```
If you are running from a different (sub)directory, make sure you point to the folder `main.nf` is located in. It is recommended to add the `-resume` command when partially rerunning the pipeline. You may also specify the output directory with the `--outdir` parameter.


# Errors & Tips
- On SLURM, exit code 137 occurs when the memory limit of a process has been exceeded. For example, SortMeRNA (but also STAR) may run out of resources, especially if your sample is large (e.g. >2GB) and contains mostly rRNA (e.g. >70%) sequences you wish to filter out. You can increase these resources in the [conf/resources.config](conf/resources.config) file.
- PureCLIP may cause an error. It is usually solved by simply rerunning using the `-resume` command.
- Genome indexing is computationally expensive and time-consuming. Therefore, we recommend generating these files once by setting `save_reference` parameter in [conf/params.config](conf/params.config) to `true` and then simply running the pipeline as normal the first time/ This saves the indices in the results directory ('04_star/genome'). Afterwards, you can move/copy the STAR index folder to the `input/` folder and specify `star_index =  '/path/to/STAR/index/' in [conf/params.config](conf/params.config), for example: '$projectDir/input/genome/index/star'.

# Future work
Possible improvements include:

- Benchmarking with public iCLIP datasets
- Including a test case
- Including a list of parameters
- Compatibility with other container runtimes (e.g. Docker, Appcontainer, conda) to improve accessibility (currently only Singularity is supported)
- Testing and improving performance (runtime, speed, resources)
- Error testing
- Testing for different operating systems (has only been tested on Linux Ubuntu 22.04.3)
- Defining computational requirements (has only been run on an HPC)
- Specify some commands in its own process file instead of the config file
- Optimize check_samplesheet.py further by implementing error raising (instead of an extra function
- Rename 'fasta' parameter to 'reference_fasta'

# Authors
This pipeline was originally made by Amarise Silié ([@amarisesilie](https://github.com/amarisesilie)) for her MSc internship in Bioinformatics from September 2023 to April 2024. The pipeline uses modules and subworkflows from [nf-core](https://github.com/nf-core/modules), [nf-core/rnaseq](https://github.com/nf-core/rnaseq) and [goodwright/clipseq](https://github.com/goodwright/clipseq) in addition to local modules and subworkflows. The paper by [Busch et al. 2020](https://www.sciencedirect.com/science/article/pii/S1046202318304948?via%3Dihub) was also used. See [CITATION.md](CITATION.md) file for a full list of references.
