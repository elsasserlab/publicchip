## Supplementary code for public ChIP data analysis.

A Nextflow pipeline is provided to re-run the analysis. If you want to run it yourself, you need to fill out the `nextflow.config` file with paths to bowtie reference index for mm9 and your reference genome. Also output directory and work directory can be provided.

Then you should be able to run it using nextflow: ::

     nextflow run public_chip.nf -c nextflow.config --acclist <your_acclistfile> 

Requires: bowtie2, picardtools, deeptools, igvtools, fastqc, samtools, bedtools.

Accession list data file is provided as part of the supplementary material of this publication.

Extra code used for peak calling and annotation is provided. This was run on a SLURM environment, so it's SLURM-based.
