#!/bin/bash

# 191029: Make an ad_hoc binning of the files (5k) and call peaks (+ reliable calls)

# local outdir="${1}"
# local inbam="${2}"
# local log="${3}"
# local specdir="${4}"
# local qval="${5}"
# local bgbam="${6}"

bam2bg="../data/other/bamtobg.csv"
bamdir="../intermediate/kingston_2016/bam"

while read bam bg; do
    sbatch ./call_peaks.sbatch ../intermediate/kingston_2016 "${bamdir}/${bam}" ../log/call_peaks.log peaks_broad_qval_05 0.05 "${bamdir}/${bg}"
done < "${bam2bg}"
