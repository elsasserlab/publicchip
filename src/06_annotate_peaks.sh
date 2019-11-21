#!/bin/bash

module load bioinfo-tools
module load BEDTools

ourpeakspath="../data/reference/bed"
peakspath="../intermediate/kingston_2016/peaks_aggregated"
bedpath="../data/reference/bed"
outdir="../intermediate/kingston_2016/peaks_annotated"

mkdir -p "${outdir}"

peakfiles="${ourpeakspath}/MACS_H33_NOT_H3K9me3.bed \
${ourpeakspath}/MACS_H33_AND_H3K9me3.bed"

kingstonpeakfiles="${peakspath}/ES_H33_00h_reliable.broadPeak \
${peakspath}/NS_H33_00h_reliable.broadPeak"

for p in ${kingstonpeakfiles}; do
    prefix=$(basename ${p})
    prefix_noext=${prefix%.*}
    # Put some extra columns on aggregated peaks (need to be tagged for this)
    awk -v peak="${prefix_noext}" -F$'\t' 'BEGIN {OFS=FS} {$4=peak; $5=0; $6="."; print;}' "${p}"  > "${outdir}/${prefix_noext}_with_tags.bed"

done

annotation_files="${bedpath}/ESC_ChromHMM15_6fields_number_only.bed,${bedpath}/IAPEz.bed,${outdir}/NS_H33_00h_reliable_with_tags.bed,${outdir}/ES_H33_00h_reliable_with_tags.bed"

for p in ${peakfiles}; do
    basep=$(basename ${p%.*})
    # Any overlap in the case of the peaks, but keep the feature that overlaps
    # the most (for chromhmm mostly)
    set -x
    python bin_annotator.py --binfiles ${p} --annotations ${annotation_files} --minoverlap 0 --out "${outdir}/${basep}_annotated.tab" --notscored
    set +x
done

rm intersection_*.bed
