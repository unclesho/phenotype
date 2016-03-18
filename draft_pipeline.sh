#!/bin/bash

### GWAS ###

printf "Running genome-wide associations..."
# 1. Run EMMAX on the trait
emmax -v -d 10 -t ../genotype/320PHY_2013WS_SEC.geno01.mind01.maf005.cleaned -p am1/am1.txt -k ../genotype/320PHY_2013WS_SEC.geno01.mind01.maf005.cleaned.hIBS.kinf -o am1/am1_emmax &> am1/am1_emmax.log
printf "DONE.\n"


# 2. Reformat the EMMAX default result.
printf "Sanitizing theh results..."
{ printf 'SNP\tCHR\tBP\tBETA\tP\n'; awk '{split($1,a,"_"); print $1 "\t" a[2]*1 "\t" a[3] "\t" $2 "\t" $3}' am1/am1_emmax.ps; } > am1/am1_emmax.ps.qqman
printf "DONE.\n"


# 3. Clump
printf "Extracting index and index-linked SNPs..."
plink2 --bfile ../genotype/320PHY_2013WS_SEC.geno01.mind01.maf005.cleaned --make-founders --clump am1/am1_emmax.ps.qqman --clump-snp-field SNP --clump-field P --clump-p1 0.0000005 --clump-kb 200 --clump-r2 0.50 --clump-p2 0.01 --r2 dprime --clump-range ../../../gene_list_MSU7.txt --clump-range-border 20 --out am1/am1_emmax.ps.qqman.emmax_200kb &> am1/am1_emmax.ps.qqman.emmax_200kb.log
printf "DONE.\n"


# 4. Creating a list of significant SNPs
printf "Creating comma-separated list of significant SNPs..."
sed '/^$/d' am1/am1_emmax.ps.qqman.emmax_200kb.clumped | tail -n+2 | awk '{print  $3 "," $12 }' | sed 's/(1)//g' | xargs | sed -e 's/ /,/g' > am1/am1_emmax.ps.qqman.emmax_200kb.clumped.snps
printf "DONE.\n"

printf "Creating a column list of significant SNPs..."
sed '/^$/d' am1/am1_emmax.ps.qqman.emmax_200kb.clumped | tail -n+2 | awk '{print  $3 "," $12 }' | sed 's/(1)//g' | xargs | sed -e 's/ /,/g' | sed 's/,/\n/g' |sort |uniq > am1/am1_emmax.ps.qqman.emmax_200kb.clumped.snps_list
printf "DONE.\n"


# 5. Create a list of the associated genes
printf "Creating a list of the associated or nearby genes..."
sed '/^$/d' am1/am1_emmax.ps.qqman.emmax_200kb.clumped.ranges | tail -n+2 | awk '{print  $7  }' | sed 's/[][]//g'| xargs | sed 's/ /,/g' | sed 's/,/\n/g' |sort |uniq > am1/am1_emmax.ps.qqman.emmax_200kb.clumped.ranges.genes
printf "DONE.\n"

### END OF GWAS RUN ###



### POST-PROCESSING ###

# 6. Given the SNPs from the previous step, create the target list for extracting from the MSU7 reference allele file.
printf "Creating the list of index and index-linked SNPs (TARGET FILE)..."
grep -F -f am1/am1_emmax.ps.qqman.emmax_200kb.clumped.snps_list ../genotype/320PHY_2013WS_SEC.geno01.mind01.maf005.cleaned.bim | sort -k 1n,4n | awk '{print "chr" sprintf("%02d",$1) "\t" $4 "\t" $2 "\t" $5 "\t" $6;}' > am1/am1_emmax.ps.qqman.emmax_200kb.clumped.snps.target_list && bgzip -c am1/am1_emmax.ps.qqman.emmax_200kb.clumped.snps.target_list > am1/am1_emmax.ps.qqman.emmax_200kb.clumped.snps.target_list.gz && tabix -s1 -b2 -e2 am1/am1_emmax.ps.qqman.emmax_200kb.clumped.snps.target_list.gz
printf "DONE.\n"


# 7. Given the SNP target list from step #1, extract the corresponding regions in the genome to know what are the reference alleles. 
#    NOTE: The plink .BIM file does not guarantee that the A2 allele is the reference allele. What it guarantees is that the A1 is **always**
#    the minor allele. So, to address that in making sure that we are annotating against the alternative allele during the SNP annotation, we
#    need to extract that genome regions identified by GWAS and cross reference the .BIM file against this.
printf "Extracting the corresponding regions in the genome..."
tabix ../../../reference_allele_MSU7.txt.gz -R am1/am1_emmax.ps.qqman.emmax_200kb.clumped.snps.target_list.gz | sort -k 1,2n > am1/am1_emmax.ps.qqman.emmax_200kb.clumped.snps.alleles_from_reference
# Note the use of "-R". It's stands for "rocket". :) The other parameter (-T) stands for "turtle".
printf "DONE.\n"


# 8. Column merge the list from the reference genome and the list derived from the gwas.
# Check if the number of lines matched
printf "Creating the annotation input file..."
LINES1=$(wc -l am1/am1_emmax.ps.qqman.emmax_200kb.clumped.snps.alleles_from_reference | awk '{print $1}')
LINES2=$(wc -l am1/am1_emmax.ps.qqman.emmax_200kb.clumped.snps.target_list | awk '{print $1}')
#LINES1=70
#LINES2=70
if [ "$LINES1" -ne "$LINES2" ]; then
	FILE2=$(date '+%s%N')

	printf "ERROR: Reference and target lists have unequal number of lines.\nTerminating.\n"
	exit 1
fi
# Post-condition: Both files have the same number of lines.

# Proceed with the rest of the script
paste am1/am1_emmax.ps.qqman.emmax_200kb.clumped.snps.alleles_from_reference am1/am1_emmax.ps.qqman.emmax_200kb.clumped.snps.target_list | awk '{if ($2 != $5) { print "Mismatched BP in line number ", NR; exit 1; } else { split($3,a,"="); if (a[2] != $8) print $1 "\t" $2 "\t" $2 "\t" a[2] "\t" $8 "\t" ";A1 in .bim is reference allele."; else print $1 "\t" $2 "\t" $2 "\t" a[2] "\t" $7 "\t" ""; }}' > am1/am1_emmax.ps.qqman.emmax_200kb.clumped.snps.avinput


# IMPORTANT: Possibly reformat the "CHROM" field of the *.avinput file because the rice_msu7 format writes it as 'Chr#'
awk '{FS="\t"; sub(/chr/,"",$1); print "Chr" $1*1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' am1/am1_emmax.ps.qqman.emmax_200kb.clumped.snps.avinput > am1/am1_emmax.ps.qqman.emmax_200kb.clumped.snps.avinput.corrected_chroms
printf "DONE.\n"


printf "Annotating..."
# 4. Annotate!
annotate_variation.pl -out am1/am1_emmax.ps.qqman.emmax_200kb.clumped.snps.avinput.corrected_chroms -build msu7b am1/am1_emmax.ps.qqman.emmax_200kb.clumped.snps.avinput.corrected_chroms ~/bin/annovar/rice_msu7/	&> am1/am1_emmax.ps.qqman.emmax_200kb.clumped.snps.avinput.corrected_chroms.log 
printf "DONE.\n"


#if [ $? -eq 0 ]
#then
#  echo "Successfully created file"
#  exit 0
#else
#  echo "Could not create file" >&2
#  exit 1
#fi





