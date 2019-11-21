#!/usr/bin/env nextflow

acclist = params.acclist

base_dir = params.base_dir
bw_index = params.bw_index
reference_dir = params.refdir
binsize = params.binsize

max_threads = params.max_threads

Channel
    .fromPath(acclist)
    .splitCsv(by:1, sep:',', skip:1,
              header: ['dataset_id', 'experiment_id','experiment_title','run_id','orig_file_descriptor','clean_file_descriptor', 'input', 'exp_type','organism_name','instrument','submitter','study_id','study_title','sample_id','sample_title','size_mb','total_runs','total_spots','total_bases','library_name','library_strategy','library_source','library_selection'])
    .map { item -> [ item['dataset_id'], item['run_id'], item['clean_file_descriptor'], item['input'] ] }
    .set { accession_ch }

process fastqDump {
    publishDir "${base_dir}/data/download/${dataset_id}", mode:'symlink'
    module 'bioinfo-tools:sratools'
    cpus = 8

    input:
    set dataset_id, run_id, file_descriptor, corr_input from accession_ch

    output:
    set dataset_id, file_descriptor, file('*1.fastq.gz'), file('*2.fastq.gz') into fastq_ch, map_ch

    script:
    // Replace with this line for testing (only dumps a few reads per file)

    // """
    // fastq-dump --outdir . --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip `srapath ${run_id}` -N 10000 -X 20000 --gzip
    // """

    """
    fasterq-dump --outdir . --skip-technical ${run_id}

    pigz -p 8 *1.fastq
    pigz -p 8 *2.fastq

    mv *1.fastq.gz ${file_descriptor}_1.fastq.gz
    mv *2.fastq.gz ${file_descriptor}_2.fastq.gz

    """
}

process fastQC {
    publishDir "${base_dir}/reports/${dataset_id}/qc", mode:'copy'
    module 'bioinfo-tools:FastQC'

    input:
    set dataset_id, prefix, fastq_r1, fastq_r2 from fastq_ch

    output:
    set dataset_id, file('*.html') into fastqc_ch

    script:
    """
    fastqc ${fastq_r1} -o .
    mv *.html ${prefix}_fastqc.html
    rm *.zip
    """
}

process mapReads {
    // publishDir "${base_dir}/analysis/${dataset_id}/bam", mode:'symlink'
    module 'bioinfo-tools:bowtie2:samtools'
    cpus = max_threads

    input:
    set dataset_id, prefix, fastq_r1, fastq_r2 from map_ch

    output:
    set dataset_id, prefix, file("*.bam") into mapqc_ch, bam_ch
    set dataset_id, prefix, file("${prefix}_bowtie_stats.log") into bowtie_ch

    script:
    """
    bowtie2 -p ${max_threads} -x ${bw_index}/mm9 -q -1 ${fastq_r1} -2 ${fastq_r2} --fast 2> ${prefix}_bowtie_stats.log | samtools view -bS -F 4  - | samtools sort -o ${prefix}.bam -
    """
}

process mapQC {
    publishDir "${base_dir}/reports/${dataset_id}/mapqc", mode:'copy'
    module 'bioinfo-tools:samtools:picard'

    // This step runs out of memory easily
    cpus 6

    input:
    set dataset_id, prefix, bamfile from mapqc_ch
    set dataset_id, prefix, bowtiestats from bowtie_ch

    output:
    set dataset_id, prefix, file('*idxstats.txt'), file('*flagstat.txt'), file('*MarkDuplicates.txt'), file('*insert_sizes.txt'), file('*insert_sizes.pdf') into mapqcstats_ch
    set dataset_id, file('dupsline.txt') into dupstats_ch
    set dataset_id, file('fragsize.txt') into fragstats_ch
    set dataset_id, bowtiestats

    script:
    """
    samtools index ${bamfile} -@ 4
    samtools flagstat ${bamfile} > ${prefix}_flagstat.txt
    samtools idxstats ${bamfile} > ${prefix}_idxstats.txt

    java -jar \$PICARD_HOME/picard.jar MarkDuplicates REMOVE_DUPLICATES=TRUE I=${bamfile} O=dedup.bam M=${prefix}_MarkDuplicates.txt
    java -jar \$PICARD_HOME/picard.jar CollectInsertSizeMetrics I=${bamfile} O=${prefix}_insert_sizes.txt H=${prefix}_insert_sizes.pdf M=0.5 STOP_AFTER=10000000

    echo "${prefix} `grep -i unk ${prefix}_MarkDuplicates.txt | cut -f 2-10`" > dupsline.txt
    echo "${prefix} `grep FR ${prefix}_insert_sizes.txt | cut -f 1`" > fragsize.txt

    rm dedup.bam
    """
}

dupstats_ch
    .collectFile(storeDir: "${base_dir}/reports/", newLine:false) { item ->
    [ "${item[0]}_picard_dedup.txt", item[1].getText() ]
}

fragstats_ch
    .collectFile(storeDir: "${base_dir}/reports/", newLine:false) { item ->
        [ "${item[0]}_picard_fragment.txt", item[1].getText() ]
}


process cleanUpBAM {
    publishDir "${base_dir}/intermediate/${dataset_id}/bam", mode:'symlink'
    module 'bioinfo-tools:picard:samtools:BEDTools'
    cpus '6'

    input:
    set dataset_id, prefix, bamfile from bam_ch

    output:
    set dataset_id, prefix, file('*fltd.bam') into cleanbam_tdf_ch, cleanbam_bw_ch

    script:
    """
    java -jar \$PICARD_HOME/picard.jar MarkDuplicates REMOVE_DUPLICATES=TRUE I=${bamfile} O=dedup.bam M=${prefix}_MarkDuplicates.txt
    samtools index dedup.bam -@ 4
    bedtools intersect -v -abam dedup.bam -b "${reference_dir}/Blacklist.bed" > ${prefix}.fltd.bam
    rm dedup.bam
    rm dedup.bam.bai
    """
}

process buildTDF {
    publishDir "${base_dir}/intermediate/${dataset_id}/tdf", mode:'symlink'
    module 'bioinfo-tools:IGVtools'

    input:
    set dataset_id, prefix, bamfile from cleanbam_tdf_ch

    output:
    file '*.tdf'

    script:
    """
    igvtools count -e 60 ${bamfile} ${prefix}.tdf ${bw_index}/mm9.chrom.sizes
    """
}

process buildBigWig {
    publishDir "${base_dir}/intermediate/${dataset_id}/bw", mode:'copy'
    module 'bioinfo-tools:samtools:deepTools'
    cpus '10'

    input:
    set dataset_id, prefix, bamfile from cleanbam_bw_ch

    output:
    set dataset_id, prefix, file('*.bw') into bw_ch

    script:
    genomesize=params.genomesize
    """
    samtools index ${bamfile} -@ 10
    bamCoverage -p max --effectiveGenomeSize ${genomesize} --normalizeUsing RPGC -b ${bamfile} -o ${prefix}.bw

    rm -f ${bamfile}.bai
    """
}

process bigWigBins {
    publishDir "${base_dir}/intermediate/${dataset_id}/bins"
    module 'bioinfo-tools:BEDTools:deepTools'
    cpus = 10

    input:
    set dataset_id, prefix, bigwig_file from bw_ch

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
