# Supplementary code for public ChIP data analysis.

A Nextflow pipeline along with extra scripts used to generate the main data outputs
used in this publications is provided here for documentation purposes.

## System requirements

This data analysis has been run under a HPC environment (UPPMAX) that runs 
CentOS Linux 7 (Core) and SLURM scheduling. However, it should be straightforward
to run a local version of the necessary scripts, just running the `sbatch` files directly.

Tools used and versions:


    nextflow (version 19.07.0 build 5106)
    bowtie2 (v 2.3.5.1)
    deepTools (v 3.1.0)
    sratools (v 2.9.6-1)
    MACS2 (v 2.1.2)
    bwtool (v1.0)
    bedtools (v2.27.1)
    matplotlib (v 3.1.1) 


## Installation guide

The software provided here is not a package and it's therefore not intended for
installation. It relies heavily on the standard wide-used tools mentioned above and
 can be run as long as the required tools are installed. It is available 
for transparency and reproducibility reasons. 

A Nextflow pipeline is provided to re-run the analysis. If you want to run it
yourself, you need to fill out the `nextflow.config` file with paths to bowtie
reference index for mm9 and reference genome. Also output directory
and work directory can be provided.

Then you should be able to run it using nextflow: ::

     nextflow run public_chip.nf -c nextflow.config --acclist <your_acclistfile> 

Running the nextflow pipeline will give you the main bigwig files also provided as
main data in our publication.

Rest of intermediate relevant files (peaks calling, annotation) are processed
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
the kind of output provided just by by providing only the first n lines in the
.csv accession list.

## Output 

Nextflow will generate an `intermediate` folder with corresponding subfolders
per dataset:

- `bw`: Bigwig files
- `bins`: Bed files with 5kb bins signal values.
- `bam`: Bowtie produced alignments.
- `tdf`: Alternative visualization files. Can be opened in IGV.

It will also generate logging information on the pipeline run.

Expected running times rely heavily on the type of machine, sizes of datasets,
internet connection and queuing status in case of running this on a HPC system.




