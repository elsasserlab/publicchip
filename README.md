# Supplementary code for public ChIP data analysis.

## Structure of the repository

- `src`: Contains source code used to perform the analysis in the publication.
- `doc`: Relevant documentation to help reproduce the analysis. Includes a `Navarro_2020_datasets_metadata.xlsx` table with metadata required to perform the analysis. Additionally, a set of tables `datasets_ebi.tsv`, `datasets_se.tsv` and `datasets_pe.tsv` is available to run with the Nextflow pipeline. 
- `data`: Files used in the analysis that are not a byproduct of its execution (for instance, bed annotations from other publications), or that may be useful per se, like peak annotations.

After running some of the analyses, you will end up with extra directories:
- `intermediate`: Files generated in the primary analysis that are required for downstream results (main figures of the publication). This includes: bigWig, bam, binned score files, called peaks. Note: Alternatively to running the primary analysis yourself, pre-generated bigwig files and other outputs are available at the data repository attached to this publication.
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

R packages required:

    rtracklayer
    GenomicRanges
    ggplot2
    reshape2
    eulerr
    ggpubr
    purrr
    effsize
    pheatmap


## Installation guide

The software provided here is not a package and it's therefore not intended for
installation. It relies heavily on the standard wide-used tools mentioned above and
 can be run as long as the required tools are installed. It is available 
for transparency and reproducibility reasons. 

A Nextflow pipeline is provided to re-run the primary analysis. If you want to run it
yourself, you need to fill out the `nextflow.config` file with paths to bowtie
reference index for mm9 and reference genome. 

Then you should be able to run it using nextflow: 

     nextflow run public_chip.nf -c nextflow.config --acclist doc/datasets_pe.tsv 

**Note**: There are currently three versions of this nextflow pipeline, for running paired-end data, single-end data and
single-end data coming from another sequencing repository. The other versions of the pipeline run with the table in `doc` that has its matching name as `--acclist`, for example:

    nextflow run public_chip_se.nf -c nextflow.config --acclist doc/datasets_se.tsv

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


## Datasets files

This study analyses data from many genomics studies. Relevant metadata needed for its download
and processing through the pipeline will be available under `doc`. The main pipeline 
is `public_chip.nf`, and there are two variants `public_chip_se.nf` (for the processing of
single-ended libraries, specifies on the `datasets_se.csv` file) and `public_ebi.nf` (for
downloading files that are in EBI database instead).

Intermediate data output from each dataset will be under `intermediate/data_type/dataset_id`,
for instance: `intermediate/bw/Deaton_2016`.

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
- `bed`: IAP annotation generated by the consensus run.

It will also generate logging information on the pipeline run.

Expected running times rely heavily on the type of machine, sizes of datasets,
internet connection and queuing status in case of running this on a HPC system.

## Downstream analysis

Code used to generate the figures on the publication is also available in `src`.

Overall: 
- Peaks were called using MACS2 on broadPeak with qval = 0.05 for each replicate. 
The aggregated peak set consists of loci that appear in at least 2 out of the 3
replicates. These are the `reliable` groups available under `peaks/Deaton_2016/aggregated`.
- Peak density heatmaps were calculated using `deepTools computeMatrix`. The
precomputed matrix was then visualized by `peak_density_annotated.py` script. 
- Average heatmaps were calculated using `rtracklayer` and `GenomicRanges` libraries
in R. Corresponding scripts are found per-figure in `src`.
- Bin-based plots also calculated in R and scripts for this are provided.
- Profile plots were calculated by a customized version of ngsplot that
calculates Y axes as 1x FPGC (Fragments per Genome Content) coverage. Very similar results 
can also be obtained by seqplots providing the corresponding bigwig and bedfiles, where
reads were extended to match fragment size in paired-end datasets. If you want to re-run
the customized version of ngsplot, you need to download their original source code and
replace the relevant files with the ones provided. Our customized version includes several
modes of calculating Y axis. The parameter you need to set is `-M relative`, as opposed to
`standard` which reproduces ngsplot usual behavior.
- IAP consensus dataset was calculated integrating several datasets `src/run_iap_consensus.sh`.
This requires the BAM files to be calculated.
- Venn intersection overlap was calculated by `intervene`, and size proportional 
euler diagrams using R package `eulerr`.

