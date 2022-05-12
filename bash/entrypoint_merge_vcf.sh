#!/usr/bin/bash

#####################################
## Helper script to merge VCF file ##
#####################################
######################
## Achal Neupane    ##
## Date: 04/22/2022 ##
######################

VERSION="1.0"

DATE="/bin/date +%s"
display_date () {
    /bin/date -d @${1} +%Y%m%d_%H%M%S
}

date_diff () {
    earlier=${1}
    later=${2}
    diff=$((${later}-${earlier}))
    if [ ${diff} -gt 86400 ]; then
        date -u -d @$((${diff}-86400)) +"%jd%-Hh%-Mm%-Ss"
    else
        date -u -d @${diff} +"%-Hh%-Mm%-Ss"
    fi
}

trap "echo \"[$(display_date $(${DATE}))] Terminated by SIGTERM \" && sleep 10s" SIGTERM

# Load module
module load bcftools/1.14

start=$(${DATE}); echo "[$(display_date ${start})] bcftools starting"

if [ ! -f "${VCF1}.tbi" ]; then	
echo "Index file missing, Creating index file for ${VCF1}"
bcftools index -t --threads ${THREADS} ${VCF1}
fi

if [ ! -f "${VCF2}.tbi" ]; then	
echo "Index file missing, Creating index file for ${VCF2}"
bcftools index -t --threads ${THREADS} ${VCF2}
fi

echo "RUNNING bcftools merge for: ${VCF1} and ${VCF2}"
bcftools merge --threads ${THREADS} ${VCF1} ${VCF2} -0 -Oz -o ${OUT_DIR}/${MERGED}.vcf.gz
exitcode=$?
end=$(${DATE}); echo "[$(display_date ${end})] bcftools finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"


# Get stats
COUNTS_VCF1="$(zcat ${VCF1} | grep -v '^#' | wc -l)"
COUNTS_VCF2="$(zcat ${VCF2} | grep -v '^#' | wc -l)"
MERGED_SAMPLES=$(zcat ${OUT_DIR}/${MERGED}.vcf.gz |head -5000| grep "^#CHROM" | tr "\t" "\n " | tail -n +10 | uniq | wc -l)
COUNTS="$(zcat ${OUT_DIR}/${MERGED}.vcf.gz | grep -v '^#' | wc -l)"

echo "${COUNTS_VCF1} variants in total in ${VCF1}"
echo "${COUNTS_VCF2} variants in total in ${VCF2}"
echo "${MERGED_SAMPLES} samples in total in merged VCF file from chr${CHR}"
echo "${COUNTS} variants in total in merged VCF file from chr${CHR}"

exit ${exitcode}