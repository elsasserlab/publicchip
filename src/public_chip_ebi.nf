#!/usr/bin/env nextflow

// Variant of the pipeline to process single-end datasets

acclist = params.acclist
base_dir = params.base_dir
bw_index = params.bw_index
reference_dir = params.refdir
binsize = params.binsize
max_threads = params.max_threads

Channel
    .fromPath(acclist)
    .splitCsv(by:1, sep:'\t', skip:1,
              header: ['dataset_id', 'ebi_id', 'description', 'filename', 'download_link', 'input', 'layout'])
    .map { item -> [ item['dataset_id'], item['ebi_id'], item['filename'], item['download_link'], item['input'] ] }
    .set { accession_ch }

process download {
    publishDir "${base_dir}/intermediate/fastq/${dataset_id}", mode:'symlink'
    module 'bioinfo-tools'
    cpus = 8

    input:
    set dataset_id, run_id, file_descriptor, download_link, corr_input from accession_ch

    output:
    set dataset_id, file_descriptor, file('*.fastq.gz') into map_ch

    script:
    
    """
    wget -nc -O "${file_descriptor}.fastq.gz" ${download_link}
    
    """
}

process mapReads {
    // publishDir "${base_dir}/analysis/${dataset_id}/bam", mode:'symlink'
    module 'bioinfo-tools:bowtie2:samtools'
    cpus = max_threads

    input:
    set dataset_id, prefix, fastq_r1 from map_ch

    output:
    set dataset_id, prefix, file("*.bam") into bam_ch

    script:
    """
    bowtie2 -p ${max_threads} -x ${bw_index}/mm9 -q -U ${fastq_r1} --fast 2> ${prefix}_bowtie_stats.log | samtools view -bS -F 4  - | samtools sort -o ${prefix}.bam -
    """
}

process cleanUpBAM {
    publishDir "${base_dir}/intermediate/bam/${dataset_id}", mode:'copy'
    module 'bioinfo-tools:picard:samtools:BEDTools'
    cpus 8

    input:
    set dataset_id, prefix, bamfile from bam_ch

    output:
    set dataset_id, prefix, file('*fltd.bam') into cleanbam_bw_ch

    script:
    """
    java -jar \$PICARD_HOME/picard.jar MarkDuplicates REMOVE_DUPLICATES=TRUE I=${bamfile} O=dedup.bam M=${prefix}_MarkDuplicates.txt
    samtools index dedup.bam -@ 4
    bedtools intersect -v -abam dedup.bam -b "${reference_dir}/Blacklist.bed" > ${prefix}.fltd.bam
    rm dedup.bam
    rm dedup.bam.bai
    """
}

process buildBigWig {
    publishDir "${base_dir}/intermediate/bw/${dataset_id}", mode:'copy'
    module 'bioinfo-tools:samtools:deepTools'
    cpus '10'

    input:
    set dataset_id, prefix, bamfile from cleanbam_bw_ch

    output:
    set dataset_id, prefix, file("${prefix}.bw"), file("${prefix}.ext.uniq.bw") into bw_ch

    script:
    genomesize=params.genomesize
    """
    samtools index ${bamfile} -@ 10
    bamCoverage -p max --effectiveGenomeSize ${genomesize} --normalizeUsing RPGC -b ${bamfile} -o ${prefix}.bw
    bamCoverage -p max --effectiveGenomeSize ${genomesize} --normalizeUsing RPGC --minMappingQuality 10 -b ${bamfile} -o ${prefix}.ext.uniq.bw
    rm -f ${bamfile}.bai
    """
}


process bigWigBins {
    publishDir "${base_dir}/intermediate/bins/${dataset_id}"
    module 'bioinfo-tools:BEDTools:deepTools'
    cpus = 10

    input:
    set dataset_id, prefix, bigwig_file, bigwig_uniq_file from bw_ch

    output:
    file "${prefix}_${binsize}_bins.bed"

    script:
    """
    multiBigwigSummary bins -b ${bigwig_file} --binSize ${binsize} --outRawCounts raw_counts.bed -o raw_counts.npz -p max -v --smartLabels
    bedtools sort -i raw_counts.bed > "${prefix}_${binsize}_bins.bed"

    rm raw_counts.bed
    rm raw_counts.npz
    """
}
