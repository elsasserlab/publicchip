#!/bin/bash

module load bioinfo-tools
module load BEDTools

aggregate() {
    ## Aggregates peaks of the three replicates.
    ## Calculates peaks that are unique to each, puts them together and removes
    ## Them from the totality of the peaks.
    ## Basically, a peak is kept if it shows up in at least 2 out of 3 replicates.
    local peaksdir="${1}"
    local cellt="${2}"
    local timep="${3}"
    local outfile="${4}"
    local peaktype="${5}"

    local prefix="${peaksdir}/${cellt}_H33_${timep}_rep"
    local suffix=".fltd_peaks.${peaktype}"

    bedtools subtract -A -a "${prefix}1${suffix}" -b "${prefix}2${suffix}" > "tmp1.bed"
    bedtools subtract -A -a "tmp1.bed" -b "${prefix}3${suffix}" > "unique_1.bed"

    bedtools subtract -A -a "${prefix}2${suffix}" -b "${prefix}1${suffix}" > "tmp1.bed"
    bedtools subtract -A -a "tmp1.bed" -b "${prefix}3${suffix}" > "unique_2.bed"

    bedtools subtract -A -a "${prefix}3${suffix}" -b "${prefix}1${suffix}" > "tmp1.bed"
    bedtools subtract -A -a "tmp1.bed" -b "${prefix}2${suffix}" > "unique_3.bed"

    cat "unique_1.bed" "unique_2.bed" "unique_3.bed" > "all_unique.bed"

    bedtools sort -i <( cat "${prefix}1${suffix}" "${prefix}2${suffix}" "${prefix}3${suffix}" ) > "all.bed"

    # Subtract all regions that are only in one replicate. Merge the overlapping ones.
    bedtools merge -i <( bedtools subtract -A -a "all.bed" -b "all_unique.bed" ) > "${outfile}"

    rm "all.bed"
    rm "all_unique.bed"
    rm "tmp1.bed"
    rm "unique_1.bed" "unique_2.bed" "unique_3.bed"

}

intersect_with_bed() {
    local bed_a="${1}"
    local bed_b="${2}"
    local bed_target="${3}"
    local counts_out_file="${4}"
    local row_label="${5}"

    local bed_a_only="a_only.tmp.bed"
    local bed_b_only="b_only.tmp.bed"
    local bed_both="both.tmp.bed"

    bedtools subtract -A -a "${bed_a}" -b "${bed_b}" > "${bed_a_only}"
    bedtools subtract -A -a "${bed_b}" -b "${bed_a}" > "${bed_b_only}"
    bedtools intersect -a "${bed_a}" -b "${bed_b}" > "${bed_both}"

    a_only_count=$(< "${bed_a_only}" wc -l)
    b_only_count=$(< "${bed_b_only}" wc -l)

    both_count=$(< "${bed_both}" wc -l)

    total=$((a_only_count + b_only_count + both_count))

    total_a=$(< "${bed_a}" wc -l)
    total_b=$(< "${bed_b}" wc -l)

    total_a_intersect=$( < <( bedtools intersect -u -a "${bed_target}" -b "${bed_a}" ) wc -l )
    total_b_intersect=$( < <( bedtools intersect -u -a "${bed_target}" -b "${bed_b}" ) wc -l )

    a_only_intersect=$( <  <( bedtools intersect -u -a "${bed_target}" -b "${bed_a_only}" ) wc -l )
    b_only_intersect=$( <  <( bedtools intersect -u -a "${bed_target}" -b "${bed_b_only}" ) wc -l )

    both_intersect=$( <  <( bedtools intersect -u -a "${bed_target}" -b "${bed_both}" ) wc -l )
    # both_intersect=$(( total_a_intersect + total_b_intersect - a_only_intersect - b_only_intersect ))

    echo -e "${row_label}\t${total}\t${total_a} (${total_a_intersect})\t${total_b} (${total_b_intersect})\t${a_only_count} (${a_only_intersect})\t${b_only_count} (${b_only_intersect})\t${both_count} (${both_intersect})" >> "${counts_out_file}"

    rm "${bed_a_only}"
    rm "${bed_b_only}"
    rm "${bed_both}"
}

main() {
    outdir="${1}"
    peaksdir="${2}"

    mkdir -p "${outdir}"
    peaktype="broadPeak"

    bedpath="../data/reference/bed"
    
    for timepoint in 00h 03h 06h 12h; do
        aggregate "${peaksdir}" "ES" "${timepoint}" "${outdir}/ES_H33_${timepoint}_reliable.${peaktype}" "${peaktype}"
        aggregate "${peaksdir}" "NS" "${timepoint}" "${outdir}/NS_H33_${timepoint}_reliable.${peaktype}" "${peaktype}"
    done

}

main "$@"
