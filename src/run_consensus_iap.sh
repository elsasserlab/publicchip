#!/bin/bash

function present_regions() {
    local bamfile=${1}
    local flankfile=${2}
    local iapfile=${3}
    local outfile=${4}
    local uniq_bam_dir=${5}

    basebam=$(basename ${bamfile})
    quality=20

    # Uniquely mapped reads
    samtools view -bq "${quality}" "${bamfile}" > "${uniq_bam_dir}/${basebam}.${quality}.bam"

    # Find reads in the flanks
    bedtools intersect -abam "${uniq_bam_dir}/${basebam}.${quality}.bam" -b "${flankfile}" -f 0.9 -bed -wb | cut -f 13-18 > "IAPEz_intersect_${basebam}.bed"

    # Recover original intervals that overlap present flanks
    bedtools intersect -a "${iapfile}" -b <(bedtools merge -i <(bedtools sort -i "IAPEz_intersect_${basebam}.bed" )) -u > "IAPEz_present_${basebam}.tmp.bed"

    # Sort 
    bedtools sort -i "IAPEz_present_${basebam}.tmp.bed" > "${outfile}"
    rm "IAPEz_present_${basebam}.tmp.bed"
    rm "IAPEz_intersect_${basebam}.bed"
}


outdir="../intermediate/bed/"
uniq_bam_dir="${outdir}/bam_uniq"
evidence_dir="${outdir}/iap_evidence"
consensus_dir="${outdir}/iap_consensus"

mkdir -p "${consensus_dir}"
mkdir -p "${evidence_dir}"
mkdir -p "${uniq_bam_dir}"

deaton=../intermediate/bam/Deaton_2016/ES_H3*.fltd.bam
gdna=../intermediate/bam/Elsasser_2015/*gDNA*.fltd.bam
h3k9=../intermediate/bam/Shi_2019/*H3K9me3*.fltd.bam

original="../data/RepMasker_IAPEz_clustered.mm9.bed"

flankfile="../intermediate/bed/RepMasker_IAPEz_clustered.flanks.250.mm9.bed"
genomes="../data/mm9.chrom.sizes"
flanklen=250

# Make flank file
bedtools slop -i ${original} -g ${genomes} -b -${flanklen} > tmp1.bed
bedtools flank -i tmp1.bed -b ${flanklen} -g ${genomes} > ${flankfile}
rm tmp1.bed

# Collect evidence
for i in ${deaton[@]} ${gdna[@]} ${h3k9[@]}; do
    basebam=$(basename ${i%.*})
    echo ${basebam}
    present_regions "${i}" "${flankfile}" "${original}" "${evidence_dir}/IAPEz_present_${basebam}.bed" "${uniq_bam_dir}"
done

# Take all the curated IAPEz present files and make the merged final dataset
# 1) Merge elements found in any of the samples from the same study
bedops --merge ${evidence_dir}/IAPEz_present_ES_H33*.bed ${evidence_dir}/IAPEz_present_ES_H31*.bed > "${consensus_dir}/kingston_ES_H3_all_merged.bed"
bedops --merge ${evidence_dir}/IAPEz_present_ESC*KO*H3K9me3* > "${consensus_dir}/H3K9m3_KO_all.bed"
bedops --merge ${evidence_dir}/IAPEz_present_ESC*WT*H3K9me3* > "${consensus_dir}/H3K9m3_WT_all.bed"
bedops --merge ${evidence_dir}/IAPEz_present_ESC_H33WT_gDNA* > "${consensus_dir}/ESC_H33WT_gDNA_all.bed"
bedops --merge ${evidence_dir}/IAPEz_present_ESC_H33KO_gDNA* > "${consensus_dir}/ESC_H33KO_gDNA_all.bed"

# 2) From these, use only the ones that are present in every study
bedops --intersect "${consensus_dir}/kingston_ES_H3_all_merged.bed" "${consensus_dir}/H3K9m3_KO_all.bed" "${consensus_dir}/H3K9m3_WT_all.bed" "${consensus_dir}/ESC_H33WT_gDNA_all.bed" "${consensus_dir}/ESC_H33KO_gDNA_all.bed" > "${outdir}/IAPEz_curated_consensus.bed"

# 3) Intersect original to keep strand info and so on
bedtools intersect -a <(bedtools sort -i "${original}") -b <(bedtools sort -i ${outdir}/IAPEz_curated_consensus.bed) -u > ${outdir}/IAPEz_curated_consensus_wstrand.bed

mv ${outdir}/IAPEz_curated_consensus_wstrand.bed ${outdir}/IAPEz_curated_consensus.bed

# 4) Bonus: Some stats on the overlap of these datasets
intervene pairwise -i ${consensus_dir}/H3K9m3_KO_all.bed ${consensus_dir}/H3K9m3_WT_all.bed ${consensus_dir}/kingston_ES_H3_all_merged.bed ${consensus_dir}/ESC_H33KO_gDNA_all.bed ${consensus_dir}/ESC_H33WT_gDNA_all.bed --names H33KO_H3K9me3,H33WT_H3K9me3,kingESC_H33/1,ESC_H33KO_gDNA,ESC_H33WT_gDNA -o ${outdir}
intervene venn -i ${consensus_dir}/*.bed -o ${outdir}
