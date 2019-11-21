#!/bin/bash -l

#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 6:00:00

module load bioinfo-tools
module load deepTools
module load MACS

ncores=4

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
    set -x
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
        macs2 callpeak --treatment "${inbam}" --control "${bgbam}" -f BAMPE -g mm --verbose ${verbosity} --outdir "${outdir}/${specdir}" -n "${inbambase%.*}" --tempdir "/scratch/$SLURM_JOB_ID" --broad --broad-cutoff "${qval}"
        loginfo "Done." "${log}" "INFO"

    else
        loginfo "Calling peaks on file ${inbam}. No control file provided..." "${log}" "INFO"
 	    macs2 callpeak --treatment "${inbam}" -f BAMPE -g mm --verbose ${verbosity} --outdir "${outdir}/${specdir}" -n "${inbambase%.*}_nobg" --tempdir "/scratch/$SLURM_JOB_ID" --broad --broad-cutoff "${qval}"
        loginfo "Done." "${log}" "INFO"
    fi


}

main "$@"
