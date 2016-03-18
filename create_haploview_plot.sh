#!/bin/bash

# Create haploview input files
printf "Creating haploview input files..."
plink2 --bfile ../genotype/320PHY_2013WS_SEC.geno01.mind01.maf005.cleaned --extract $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps_list --recode HV --out $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps_list &> /dev/null
# generate haploview batch input files
#find $1 -name '*snps_list.chr*' -printf '%P\n' | sort | paste -s -d ' \n' | awk '{print $2 "\t" $1;}'
ls $1/*snps_list.chr* | paste -s -d ' \n' | awk '{print $2 "\t" $1;}' > $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps_list.haploview_batch_input
printf "DONE.\n"


printf "Generating LD plots..."
awk '{system("java -jar ~/bin/Haploview.jar -n -pedfile " $1  " -info " $2 " -skipcheck -blockoutput ALL -dprime -out " ARGV[1] " -svg"); print $0}' $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps_list.haploview_batch_input &> /dev/null
printf "DONE.\n"
