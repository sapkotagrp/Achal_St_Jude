##########################################################
### filtering the vcf-files with bcftools and vcftools ###
###########################################################

# load tools: samtools, vcftools

##############
# Victoria 1 # 
##############

# remove all individuals from family 1 and 3.3, update tags, bgzip file
bcftools view -s ^21553.F1mum11,21644.F1dad11,21675.F2male11,21676.F2male11,21745.F2male11,21746.F2male11,21747.F2male11,21880.F2male11,21881.F2male11,21900.F2male11,21901.F2male11,21902.F2male11,21903.F2male11,21613.F1dad12,21614.F1mum12,21631.F2male12,21632.F2male12,21751.F2male12,21752.F2male12,21753.F2male12,21754.F2male12,21907.F2male12,21908.F2male12,21909.F2male12,21936.F2male12,21937.F2male12,21643.F1mum13,21570.F1mum33,21713.F2male33,21716.F2male33,21717.F2male33,21720.F2male33,21721.F2male33,21722.F2male33,21714.F2male33r Victoria1_final.vcf.gz | bcftools +fill-tags | bgzip -c > Victoria1_final_Fam3.vcf.gz
tabix -fp vcf Victoria1_final_Fam3.vcf.gz

# filter 1 (f1): only keep sites with max. 50% missing data 
bcftools view -i 'AN/(2*N_SAMPLES)>0.5' Victoria1_final_Fam3.vcf.gz | bcftools +fill-tags | bgzip -c > Victoria1_final_Fam3_f1.vcf.gz
tabix -fp vcf Victoria1_final_Fam3_f1.vcf.gz

# calculate mean depth per individual
bsub "vcftools --gzvcf Victoria1_final_Fam3_f1.vcf.gz --depth" 
# calculate proportion of missing data per individual
bsub "vcftools --gzvcf Victoria1_final_Fam3_f1.vcf.gz --missing-indv" 
# merge the two files
paste out.idepth out.imiss > Victoria1_final_indstats.txt
# download it and plot mean depth and mean depth vs missing data per individual in R 

# remove individuals with >50% missing data and <15 mean depth (r1)
bcftools view -s ^22108.F2female34,22283.F2male32,22319.F2male32 Victoria1_final_Fam3_f1.vcf.gz | bcftools +fill-tags | bgzip -c > Victoria1_final_Fam3_f1_r1.vcf.gz
tabix -fp vcf Victoria1_final_Fam3_f1_r1.vcf.gz

# make submssion script with the following filters (f2):
# a) calculate high mean depth cutoff (mpir)
# b) remove (10) snps around indels and set genotypes with DP<10 to missing 
# c) keep only bi-allelic snps and sites with DP<mpir
# d) output new compressed file 
# e) index it

#!/bin/bash
#BSUB -J "filters2"
#BSUB -R "rusage[mem=1000]"
#BSUB -n 1
#BSUB -W 1:00

bcftools stats -d 0,1000000,1 Victoria1_final_Fam3_f1_r1.vcf.gz > Victoria1_final_Fam3_f1_r1.vcf.gz.stats
mpir=$(cat Victoria1_final_Fam3_f1_r1.vcf.gz.stats | grep "^DP" | awk 'BEGIN{q1=0;q3=0}{a+=$7;d+=$6*$3;t+=$6;if(a>=25 && q1==0){q1=$3};if(a>=75 && q3==0){q3=$3}}END{print d/t+1.5*(q3-q1)}')

bcftools filter -g 10 -i 'FORMAT/DP>9' -S . Victoria1_final_Fam3_f1_r1.vcf.gz | \
bcftools view -m2 -M2 -v snps -i 'INFO/DP<'$mpir'' | \
bgzip -c > Victoria1_final_Fam3_f1_r1_f2.vcf.gz
tabix -fp vcf Victoria1_final_Fam3_f1_r1_f2.vcf.gz

# run it
bsub < filters2.lsf

# do stats before and after (check snps, indels, mpir cutoff)
bcftools stats -d 0,1000000,1 Victoria1_final_Fam3_f1_r1.vcf.gz > filters2_statsbefore.txt
bcftools stats -d 0,1000000,1 Victoria1_final_Fam3_f1_r1_f2.vcf.gz > filters2_statsafter.txt

# check individuals for severity of PCR duplicates 
# unzip vcf file first, as required by script
bgzip -c -d Victoria1_final_Fam3_f1_r1_f2.vcf.gz > Victoria1_final_Fam3_f1_r1_f2_checkhets.vcf
# load python module and run it
bsub -n 1 -R "rusage[mem=10000]" -W 4:00 "checkHetsIndvVCF.sh Victoria1_final_Fam3_f1_r1_f2_checkhets.vcf" 
# (script loads vcftools15 and R automatically)

# remove the very worst individuals
bcftools view -s ^21997.F2male32r,22091.F2male32,22094.F2male32,22095.F2male32 Victoria1_final_Fam3_f1_r1_f2.vcf.gz | bcftools +fill-tags | bgzip -c > Victoria1_final_Fam3_f1_r1_f2_r2.vcf.gz
tabix -fp vcf Victoria1_final_Fam3_f1_r1_f2_r2.vcf.gz

# remove sites max. 50% missing data and apply MAF >0.05 filter (f3), update tags, output new compressed file, index it
bcftools view -i 'AN/(2*N_SAMPLES)>0.5 & INFO/MAF>0.05' Victoria1_final_Fam3_f1_r1_f2_r2.vcf.gz | bcftools +fill-tags | bgzip -c > Victoria1_filtered.gz
tabix -fp vcf Victoria1_filtered.gz
# do stats 
bcftools stats Victoria1_filtered.gz > filters3_stats.txt



# for JoinMap input and QTL mapping, subset F2s to positions that are homozygous fixed in grandparents,
# and heterozygous in F1s

# first apply allelic balance filter
# (do it for F0, get homozygous fixed positions and apply these to F2s)
# subset to F0 grandmother 
bcftools view -s 65145.F0mum3 Victoria1_filtered.gz | bcftools +fill-tags | bgzip -c > onlyF0mum.vcf.gz
tabix -fp vcf onlyF0mum.vcf.gz
# decompress for next step
bgzip -c -d onlyF0mum.vcf.gz > onlyF0mum.vcf
# apply allelic balance filter (load Python module first)
allelicBalance.py -i onlyF0mum.vcf -hom -o onlyF0mum_balanced.vcf
# update tags (important here!) and zip again
bcftools +fill-tags onlyF0mum_balanced.vcf | bgzip -c > onlyF0mum_balanced.vcf.gz
tabix -fp vcf onlyF0mum_balanced.vcf.gz

# get all homozygous positions
bcftools view -H -i 'GT="hom"' onlyF0mum_balanced.vcf.gz | awk '{print $1"\t"$2}' > hombalpositionsF0mum

# apply 'hombalpositons' to F1 and get positions that are heterozygous in all 4 F1s (no missing data)
bcftools view -s 21639.F1dad34,21642.F1mum34,22044.F1dad32,22045.F1mum32 -T hombalpositionsF0mum Victoria1_filtered.gz | bcftools +fill-tags | bcftools view -H -i 'COUNT(GT="hom")=0 & NS=4' | awk '{print $1"\t"$2}' > hombalF0mumF1hetpositions

# apply to F2s
bcftools view -s ^21639.F1dad34,21642.F1mum34,22044.F1dad32,22045.F1mum32,65119.F0dad3,65145.F0mum3 -T hombalF0mumF1hetpositions Victoria1_filtered.gz | bcftools +fill-tags | bgzip -c > Victoria1_filtered_hombalF0mumF1het.vcf.gz
tabix -fp vcf Victoria1_filtered_hombalF0mumF1het.vcf.gz

# also apply to F0mum file (since she's not in file with F2s)
bcftools view -T hombalF0mumF1hetpositions onlyF0mum_balanced.vcf.gz | bgzip -c > onlyF0mum_balanced_hombalF0mumF1het.vcf.gz
tabix -fp vcf onlyF0mum_balanced_hombalF0mumF1het.vcf.gz



#############
# Victoria2 #
#############

# to count number of individuals and number of positions (sites) do:
bcftools query -l Victoria2.vcf.gz | wc -l 
bcftools query -f '%POS\n' Victoria2.vcf.gz | wc -l 

# remove two 'undet' individuals, update tags, bgzip file (keep the 'mix' family 3 individuals for now) (r1)
bcftools view -s ^21333.undet,21334.undet Victoria2.vcf.gz | bcftools +fill-tags | bgzip -c > Victoria2_r1.vcf.gz
tabix -fp vcf Victoria2_r1.vcf.gz

# filter 1 (f1): only keep sites with max. 50% missing data 
bcftools view -i 'AN/(2*N_SAMPLES)>0.5' Victoria2_r1.vcf.gz | bcftools +fill-tags | bgzip -c > Victoria2_r1_f1.vcf.gz
tabix -fp vcf Victoria2_r1_f1.vcf.gz

# calculate mean depth per individual
bsub "vcftools --gzvcf Victoria2_r1_f1.vcf.gz --depth" 
# calculate proportion of missing data per individual
bsub "vcftools --gzvcf Victoria2_r1_f1.vcf.gz --missing-indv" 
# merge the two files
paste out.idepth out.imiss > Victoria2_indstats.txt
# download it and plot mean depth and mean depth vs missing data per individual in R 

# remove individuals with >50% missing data and/or <12 mean depth (r2)
bcftools view -s ^21337.F2female1,21435.F2male1,21509.F2male2r,21512.F2female2r Victoria2_r1_f1.vcf.gz | bcftools +fill-tags | bgzip -c > Victoria2_r1_f1_r2.vcf.gz
tabix -fp vcf Victoria2_r1_f1_r2.vcf.gz

# make submssion script with the following filters (f2):
# a) calculate high mean depth cutoff (mpir)
# b) remove (10) snps around indels and set genotypes with DP<10 to missing 
# c) keep only bi-allelic snps and sites with DP<mpir
# d) output new compressed file 
# e) index it

#!/bin/bash
#BSUB -J "filters2"
#BSUB -R "rusage[mem=1000]"
#BSUB -n 1
#BSUB -W 0:30

bcftools stats -d 0,1000000,1 Victoria2_r1_f1_r2.vcf.gz > Victoria2_r1_f1_r2.vcf.gz.stats
mpir=$(cat Victoria2_r1_f1_r2.vcf.gz.stats | grep "^DP" | awk 'BEGIN{q1=0;q3=0}{a+=$7;d+=$6*$3;t+=$6;if(a>=25 && q1==0){q1=$3};if(a>=75 && q3==0){q3=$3}}END{print d/t+1.5*(q3-q1)}')

bcftools filter -g 10 -i 'FORMAT/DP>9' -S . Victoria2_r1_f1_r2.vcf.gz | \
bcftools view -m2 -M2 -v snps -i 'INFO/DP<'$mpir'' | \
bgzip -c > Victoria2_r1_f1_r2_f2.vcf.gz
tabix -fp vcf Victoria2_r1_f1_r2_f2.vcf.gz

# run it
bsub < filters2.lsf

# do stats before and after (check snps, indels, mpir cutoff)
bcftools stats -d 0,1000000,1 Victoria2_r1_f1_r2.vcf.gz > filters2_statsbefore.txt
bcftools stats -d 0,1000000,1 Victoria2_r1_f1_r2_f2.vcf.gz > filters2_statsafter.txt

# check individuals for severity of PCR duplicates 
# unzip vcf file first, as required by script
bgzip -c -d Victoria2_r1_f1_r2_f2.vcf.gz > Victoria2_r1_f1_r2_f2_checkhets.vcf
# load python module and run it
bsub -n 1 -R "rusage[mem=10000]" -W 4:00 "checkHetsIndvVCF.sh Victoria2_r1_f1_r2_f2_checkhets.vcf" 
# (loads vcftools15 and R automatically)

# remove the worst 13 individuals
bcftools view -s ^21343.F2male1r,21350.F2male1r,21415.F2female1,21418.F2female1,21419.F2female1,21421.F2female1,21422.F2female1,21423.F2female1r,21436.F2female1,21437.F2male1,21441.F2male1,21442.F2female1,21451.F2female1 Victoria2_r1_f1_r2_f2.vcf.gz | bcftools +fill-tags | bgzip -c > Victoria2_r1_f1_r2_f2_r3.vcf.gz
tabix -fp vcf Victoria2_r1_f1_r2_f2_r3.vcf.gz

# remove sites max. 50% missing data and apply MAF >0.05 filter (f3), update tags, output new compressed file, index it
bcftools view -i 'AN/(2*N_SAMPLES)>0.5 & INFO/MAF>0.05' Victoria2_r1_f1_r2_f2_r3.vcf.gz | bcftools +fill-tags | bgzip -c > Victoria2_filtered.vcf.gz
tabix -fp vcf Victoria2_filtered.vcf.gz
# do stats 
bcftools stats Victoria2_filtered.vcf.gz > filters3_stats.txt



# for JoinMap input and QTL mapping, subset F2s to positions that are homozygous fixed in grandparents,
# and heterozygous in F1s

# first apply allelic balance filter
# (do it for F0, get homozygous fixed positions and apply these to F2s)
bcftools view -s V7A10.F0maleD,V8F6.F0femaleD Victoria2_filtered.vcf.gz | bcftools +fill-tags | bgzip -c > onlyF0.vcf.gz
tabix -fp vcf onlyF0.vcf.gz
# decompress for next step
bgzip -c -d onlyF0.vcf.gz > onlyF0.vcf
# apply allelic balance filter (load python module first!)
allelicBalance.py -i onlyF0.vcf -hom -o onlyF0_balanced.vcf
# update tags (important here!) and zip again
bcftools +fill-tags onlyF0_balanced.vcf | bgzip -c > onlyF0_balanced.vcf.gz
tabix -fp vcf onlyF0_balanced.vcf.gz

# get homozygous fixed positions 
bcftools view -H -i 'INFO/AC_Het=0 & INFO/AC_Hom=2 & INFO/AF=0.5' onlyF0_balanced.vcf.gz | awk '{print $1"\t"$2}' > homfixbalpositions
# count sites
cat homfixbalpositions | wc -l # 3'473

# apply homfixpositons to F1 and get positions that are heterozygous in both F1s with confident identity
bcftools view -s 21427.F1female1,21428.F1male1 -T homfixbalpositions VVictoria2_filtered.vcf.gz | bcftools +fill-tags | bcftools view -H -i 'COUNT(GT="hom")=0 & NS=2' | awk '{print $1"\t"$2}' > homfixbalF1hetpositions
# count sites
cat homfixbalF1hetpositions | wc -l # 2'358

# apply to F2 and F0
bcftools view -s ^21427.F1female1,21428.F1male1,21320.F1male2r,21321.F1female2 -T homfixbalF1hetpositions Victoria2_filtered.vcf.gz | bcftools +fill-tags | bgzip -c > Victoria2_filtered_homfixbalF1het.vcf.gz
tabix -fp vcf VVictoria2_filtered_homfixbalF1het.vcf.gz

# also apply homfix and F1het positions to F0 file (needed because in Victoria2_filtered_homfixbalF1het.vcf.gz allelic balance filter not applied on F0s, i.e. they still have the old genotypes)
bcftools view -T homfixbalF1hetpositions onlyF0_balanced.vcf.gz | bgzip -c > onlyF0_homfixbalF1het.vcf.gz
tabix -fp vcf onlyF0_homfixbalF1het.vcf.gz



##########
# Malawi #
###########

# filter 1 (f1): only keep sites with max. 50% missing data 
bcftools view -i 'AN/(2*N_SAMPLES)>0.5' Malawi.vcf.gz | bcftools +fill-tags | bgzip -c > Malawi_f1.vcf.gz
tabix -fp vcf Malawi_f1.vcf.gz

# calculate proportion of missing data per individual
vcftools --gzvcf Malawi_f1.vcf.gz --missing-indv 
# calculate mean depth per individual
vcftools --gzvcf Malawi_f1.vcf.gz --depth
# merge the two corresponding files
paste out.idepth out.imiss > Malawi_indstats.txt
# download it and plot mean depth and mean depth vs missing data per individual in R

# remove individuals with >50% missing data and/or <15 mean depth (r1)
bcftools view -s ^2468.CALTAEF2,2230.CALTAEF2,2554.CALTAEF2,2685.CALTAEF2,2788.CALTAEF2,68215.CALTAEF2,68322.CALTAEF2,68237.CALTAEF2,68255.CALTAEF2,68295.CALTAEF2,68304.PROTAE,68305.PROTAE Malawi_f1.vcf.gz | bcftools +fill-tags | bgzip -c > Malawi_f1_r1.vcf.gz
tabix -fp vcf Malawi_f1_r1.vcf.gz

# make submission script with the following filters (f2):
# a) calculate high mean depth cutoff (mpir)
# b) remove (10) snps around indels and set genotypes with DP<10 to missing 
# c) keep only bi-allelic snps and sites with DP<mpir
# d) output new compressed file 
# e) index it

#!/bin/bash
#BSUB -J "filters2"
#BSUB -R "rusage[mem=1000]"
#BSUB -n 1
#BSUB -W 1:00

bcftools stats -d 0,1000000,1 Malawi_f1_r1.vcf.gz > Malawi_f1_r1.vcf.gz.stats
mpir=$(cat Malawi_f1_r1.vcf.gz.stats | grep "^DP" | awk 'BEGIN{q1=0;q3=0}{a+=$7;d+=$6*$3;t+=$6;if(a>=25 && q1==0){q1=$3};if(a>=75 && q3==0){q3=$3}}END{print d/t+1.5*(q3-q1)}')

bcftools filter -g 10 -i 'FORMAT/DP>9' -S . Malawi_f1_r1.vcf.gz | \
bcftools view -m2 -M2 -v snps -i 'INFO/DP<'$mpir'' | \
bgzip -c > Malawi_f1_r1_f2.vcf.gz
tabix -fp vcf Malawi_f1_r1_f2.vcf.gz

# run it
bsub < filters2.lsf

# do stats before and after (check snps, indels, mpir cutoff)
bcftools stats -d 0,1000000,1 Malawi_f1_r1.vcf.gz > filters2_statsbefore.txt
bcftools stats -d 0,1000000,1 Malawi_f1_r1_f2.vcf.gz > filters2_statsafter.txt

# check individuals for severity of PCR duplicates (as above, see Victoria2)
# exclude individuals with bad profiles
bcftools view -s ^1129.ASTCAL,1522.CALTAEF2,1529.CALTAEF2,1532.CALTAEF2,1533.CALTAEF2,1833.CALTAEF2,1897.CALTAEF2,1899.CALTAEF2,1905.CALTAEF2,2211.CALTAEF2,2547.CALTAEF2,2553.CALTAEF2,2688.CALTAEF2,2784.CALTAEF2,2790.CALTAEF2,2791.CALTAEF2,2792.CALTAEF2,2793.CALTAEF2,2811.CALTAEF2,2812.CALTAEF2,68212.CALTAEF2,68227.CALTAEF2,68228.CALTAEF2,68229.CALTAEF2,68236.CALTAEF2,68241.CALTAEF2,68264.CALTAEF2,68274.PROTAE,68282.CALTAEF2,68300.PROTAE,68335.PROTAE,85.ASTCAL,68342.ASTCAL Malawi_f1_r1_f2.vcf.gz | bcftools +fill-tags | bgzip -c > Malawi_f1_r1_f2_r2.vcf.gz
tabix -fp vcf Malawi_f1_r1_f2_r2.vcf.gz

# remove sites max. 50% missing data and APPLY MAF >0.05 filter (f3), update tags, output new compressed file, index it
bcftools view -i 'AN/(2*N_SAMPLES)>0.5 & INFO/MAF>0.05' Malawi_f1_r1_f2_r2.vcf.gz | bcftools +fill-tags | bgzip -c > Malawi_filtered.vcf.gz
tabix -fp vcf Malawi_filtered.vcf.gz
# do stats 
bcftools stats Malawi_filtered.vcf.gz > filters3_stats.txt



# for JoinMap input and QTL mapping, subset F2s to positions that are homozygous fixed in grandparents,
# and heterozygous in F1s

# apply allelic balance filter
# (do it for 'artificial' F0, get homozygous fix positions and apply these to F2s)
# subset to 'artificial' F0s 
bcftools view -s 177.ASTCAL,68271.ASTCAL,68272.ASTCAL,68275.PROTAE,68313.PROTAE,68314.ASTCAL,68315.ASTCAL,68317.ASTCAL,68318.PROTAE,68327.ASTCAL,68330.PROTAE,68350.PROTAE Malawi_filtered.vcf.gz | bcftools +fill-tags | bgzip -c > MalawionlyF0.vcf.gz
tabix -fp vcf MalawionlyF0.vcf.gz
# decompress for next step
bgzip -c -d MalawionlyF0.vcf.gz > MalawionlyF0.vcf
# apply allelic balance filter (load Python module first!)
allelicBalance.py -i MalawionlyF0.vcf -hom -o MalawionlyF0balanced.vcf
# update tags (important here!) and zip again
bcftools +fill-tags MalawionlyF0balanced.vcf | bgzip -c > MalawionlyF0balanced.vcf.gz
tabix -fp vcf MalawionlyF0balanced.vcf.gz

# for this cross get homozygous fixed positions within R-script "Malawi_vcf_to_JoinMap.R"

