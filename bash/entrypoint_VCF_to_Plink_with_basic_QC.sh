#!/usr/bin/bash

###########################################
## Helper script to convert VCF to plink ##
###########################################
## Achal Neupane ##
## Date: 04/25/2025 ##
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2

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

start=$(${DATE}); echo "[$(display_date ${start})] Starting Plink1.9 to convert VCF ---> PLINK"

CURR_STEP="/research_jude/rgs01_jude/groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/plink"

echo "RUNNING VCF to plink conversion: ${VCF}"
${CURR_STEP} --vcf ${VCF} --double-id --vcf-half-call m --keep-allele-order --threads ${THREADS} --make-bed --out ${WORKDIR}/${PLINK_FILE}

echo "RUNNING basic QC for --geno 0.1 and hwe 1e-10"
${CURR_STEP} --bfile ${WORKDIR}/${PLINK_FILE} --allow-no-sex --keep-allele-order --threads ${THREADS} --geno 0.1 --hwe 1e-10 --make-bed --out ${WORKDIR}/${PLINK_FILE}_geno.0.1_hwe.1e-10
rm "${WORKDIR}/${PLINK_FILE}.bed"; rm "${WORKDIR}/${PLINK_FILE}.bim"; rm "${WORKDIR}/${PLINK_FILE}.fam"
exitcode=$?

end=$(${DATE}); echo "[$(display_date ${end})] ${CURR_STEP} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"

exit ${exitcode}