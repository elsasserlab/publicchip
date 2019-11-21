#!/bin/bash

module load bioinfo-tools
module load BEDTools

binpath="../intermediate/kingston_2016/bins_5000"
otherbins="../intermediate/simon_gen/bins_5000"
bedpath="../data/simon_gen/191028_repeats_curated"

outdir="../results/full_results/191029_annotated_bins"
mkdir -p "${outdir}"

# Merge beds before annotating
cat "${bedpath}/clustered/"*.bed > tmp_repeats_clustered.bed
cat "${bedpath}/individual/"*.bed > tmp_repeats_individual.bed

bedtools sort -i  tmp_repeats_clustered.bed > tmp_repeats_clustered_sorted.bed
bedtools sort -i  tmp_repeats_individual.bed > tmp_repeats_individual_sorted.bed

binfiles="${binpath}/ES_H33_00h_rep1_5000_bins.bed,\
${binpath}/ES_H33_inp_rep1_5000_bins.bed,\
${otherbins}/2019_Shi_ESC1_WT_H3K9me3_5000_bins.bed,\
${otherbins}/2019_Shi_ESC_INPUT_H3K9me3_5000_bins.bed,\
${otherbins}/Bunch_S2_PolII_WT_rep1_5000_bins.bed,\
${otherbins}/Bunch_Input_WT_rep1_5000_bins.bed,\
${otherbins}/Bunch_Total_PolII_WT_rep1_5000_bins.bed,\
${otherbins}/ESC_H33KO_ATAC_5000_bins.bed,\
${otherbins}/ESC_H33WT_ATAC_5000_bins.bed"

python bin_annotator.py --binfiles "${binfiles}" --annotations tmp_repeats_clustered.bed --minoverlap 0.6 --out "${outdir}/H33_0h_vs_H3k9me3_and_more_clustered_repeats_06_overlap.tab"
python bin_annotator.py --binfiles "${binfiles}" --annotations tmp_repeats_individual.bed --minoverlap 0.6 --out "${outdir}/H33_0h_vs_H3k9me3_and_more_individual_repeats_06_overlap.tab"

rm tmp_repeats_clustered.bed
rm tmp_repeats_clustered_sorted.bed
rm tmp_repeats_individual.bed
rm tmp_repeats_individual_sorted.bed
