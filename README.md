# Supplementary code for public ChIP data analysis.

## Structure of the repository

- `src`: Contains source code used to perform the analysis in the publication.
- `doc`: Relevant documentation to help reproduce the analysis. Includes a main table `datasets.csv` with metadata required to perform the primary analysis. This will be available upon publication.
- `data`: Files used in the analysis that are not a byproduct of its execution. Mostly bed annotations.

After running some of the analyses, you will end up with extra directories:
- `intermediate`: Files generated in the primary analysis that are required for downstream results (main figures of the publication). This includes: bigWig, bam, binned score files, called peaks. Note: Alternatively to running the primary analysis yourself, pre-generated bigwig files will be provided from an alternative source.
- `figures`: Output directory where the publication results go.

## System requirements

This data analysis has been run under a HPC environment (UPPMAX) that runs 
CentOS Linux 7 (Core) and SLURM scheduling. The provided code here is runnable on a local computer. However,
the primary analysis can be computationally expensive, so you can download the generated bigwig files instead
(these will be available upon publication).

Tools used and versions:

    nextflow (version 19.07.0 build 5106)
    bowtie2 (v 2.3.5.1)
    deepTools (v 3.1.0)
    sratools (v 2.9.6-1)
    MACS2 (v 2.1.2)
    bwtool (v1.0)
    bedtools (v2.27.1)
    matplotlib (v 3.1.1)
    BEDOPS (v 2.4.38)
    intervene (v 0.6.4)
    R (v 3.6.0)
    

## Installation guide

The software provided here is not a package and it's therefore not intended for
installation. It relies heavily on the standard wide-used tools mentioned above and
 can be run as long as the required tools are installed. It is available 
for transparency and reproducibility reasons. 

A Nextflow pipeline is provided to re-run the primary analysis. If you want to run it
yourself, you need to fill out the `nextflow.config` file with paths to bowtie
reference index for mm9 and reference genome. 

Then you should be able to run it using nextflow: 

     nextflow run public_chip.nf -c nextflow.config --acclist <your_acclistfile> 

Running the nextflow pipeline will by default output everything on the `intermediate` directory.
This can be used to run the rest of the analyses available under `src`.

A few more relevant intermediate files (peaks calling and aggregation, annotation) are processed
through the corresponding `sh` scripts.

## Running the nextflow pipeline locally

You'll need the accession file list (as provided in the supplementary data of
this publication) and fill out corresponding `nextflow.config` values (an example 
`nextflow.config` file under `sample`):

    params {
        base_dir = './'
        acclist = 'acclist.csv'
        bw_index = 'path to bowtie index of mm9 genome'
        refdir = 'path to mm9 reference genome'
        genomesize = 2150570000
        max_threads = '12'
        binsize = 5000
    }

    profiles {
        standard {
            process.executor = 'local'
        }
    }

## Running a small test example

This workflow processes lines in the accession list one by one. You can check
the kind of output provided just by by providing only the first couple of lines in the
.csv accession list.

## Output 

Nextflow will generate an `intermediate` folder with corresponding subfolders
per dataset:

- `bw`: Bigwig files.
- `bins`: Bed files with 5kb bins signal values.
- `bam`: Bowtie produced alignments.
- `peaks`: Peaks called for some of the samples.

It will also generate logging information on the pipeline run.

Expected running times rely heavily on the type of machine, sizes of datasets,
internet connection and queuing status in case of running this on a HPC system.


