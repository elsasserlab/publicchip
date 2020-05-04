#!/bin/bash -l

#SBATCH -A snic2020-15-9
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 20:00:00
#SBATCH -M snowy
#SBATCH -J callpeaks

module load bioinfo-tools
module load deepTools
module load MACS

timestamp() {
    date +"[%F %H:%M:%S]"
}

loginfo() {
    local msg="${1}"
    local logfile="${2}"
    local defaultlevel="INFO"
    local level="${3:-$defaultlevel}"

    printf "%s %-8s %s\n" "$(timestamp)" "${level}" "${msg}" >> "${logfile}"
}

main() {
    local outdir="${1}"
    local inbam="${2}"
    local log="${3}"
    local specdir="${4}"
    local qval="${5}"
    local bgbam="${6}"

    local verbosity="2"

    mkdir -p "${outdir}/${specdir}"

    inbambase=$( basename ${inbam} )
    # outpeaks="${outdir}/${specdir}/${inbambase%.*}_${binsize}.counts.tsv"

    if [ -e "${bgbam}" ]; then
        macs2 callpeak --treatment "${inbam}" --control "${bgbam}" -f BAMPE -g mm --verbose ${verbosity} --outdir "${outdir}/${specdir}" -n "${inbambase%.*}" --broad --broad-cutoff "${qval}"
        loginfo "Done." "${log}" "INFO"

    else
        loginfo "Calling peaks on file ${inbam}. No control file provided..." "${log}" "INFO"
        macs2 callpeak --treatment "${inbam}" -f BAMPE -g mm --verbose ${verbosity} --outdir "${outdir}/${specdir}" -n "${inbambase%.*}_nobg" --tempdir "/scratch/$SLURM_JOB_ID" --broad --broad-cutoff "${qval}"
        loginfo "Done." "${log}" "INFO"
    fi


}

# File which contains corresponding input bam for each bam file
bam2bg="../doc/bamtobg_sample.csv"
bamdir="../intermediate/bam/Deaton_2016"

while read bam bg; do
    main ../intermediate/peaks/Deaton_2016 "${bamdir}/${bam}" ../log/call_peaks.log . 0.05 "${bamdir}/${bg}"
done < "${bam2bg}"
