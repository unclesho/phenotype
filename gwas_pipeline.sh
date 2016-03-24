#!/bin/bash

# This script is used to execute the GWAS pipeline.
# Important:
# 1. This script needs to be the at the root of the directory.
# 2. Traits are organized as immediate subdirectories.
# 3. Inside each trait subdirectory should be the plink-formatted phenotype file.
# 4. In a same-level directory (to where this script is stored), should be directory named "genotype"
#    where the required genotype file are stored.
# 5. This is an example of where this file is located:
#    /home/roslen/scratch2/gwas_hdra_v3/run/320PHY_2013WS_SEC/phenotype
# 6. Below is the directory tree. The location where this script is place is indicated by the "<<<***"
#    /home/roslen/scratch2/gwas_hdra_v3
#    |-- genotype
#    |-- phenotype
#    `-- run
#        `-- 320PHY_2013WS_SEC
#            |-- genotype
#            |   `-- scripts
#            `-- phenotype <<<***
#                |-- am1
#                |-- am2
#                |-- app_am
#                |-- app_ap
#                |-- mc_ap
#                `-- sc_ap
#
# Author: Roslen Anacleto (IRRI Grain Quality and Nutrition Center)
# Date: 20160304
#



# Start...

VERSION=0.5.0
# v0.2:
#	- added generation of manhattan and QQ plots
# 	- added generation of haploview plots
#
# v0.3:
#  - added annotation of genes
#
# v0.4
#  - corrected the gwas component. TPED and hIBS.Inf need to be specified
#
# v0.4.1
#  - corrected the case when there are no linked SNPs to the index SNP. sp2=NONE
#  - took care of the case when there are no clumped SNPs.
#  - fixed the bug that creates only one LD plot.
#
# v0.4.3
#  - fixed the reformat the "CHROM" of aviinput file
#  - added create_excel_summary
#
# v0.4.4
#  - added a helper script to create the emmax input files
#  "prepare_emmax_files.sh"
#
# v0.4.5
#  - Corrected the step 6 where "haploview_batch_input" is created
#
# v0.5.0
#  - Corrected the haploview call to output only GABRIEL blocks
#  - Corrected the create_excel_summary.R file to also include the tag SNPs sheet.
#  - Created the create_tagSNP_lanes.R script to output the tagSNP combinations
#  - Created a haploview replacement header
#



printf "GWAS pipeline version $VERSION.\n"


### Perform diagnostics ###

# i. Check if there is a command-argument
if [ "$1" = "" ]; then
	printf "Usage: gwas_pipeline <root_dir_of_trait>\n"
	exit 1
fi
printf "Trait: $1.\n"


# ii. Check if the directory is present
if [ ! -d "$1" ]; then
	  # Control will enter here if $DIRECTORY doesn't exist.
	printf "$1 directory not found.\n"
	exit 2
fi

# Check for the existence of a required files.
#if [ -f .bash_profile ]; then
#    echo "You have a .bash_profile. Things are fine."
#else
#    echo "Yikes! You have no .bash_profile!"
#fi


### GWAS ###

printf "Running genome-wide associations..."
# 1. Run EMMAX on the trait
emmax -v -d 10 -t ../genotype/$2 -p $1/$1.txt -k ../genotype/$2.hIBS.kinf -o $1/$1_emmax &> $1/$1_emmax.log
printf "DONE.\n"



# 2. Reformat the EMMAX default result.
printf "Sanitizing the results..."
{ printf 'SNP\tCHR\tBP\tBETA\tP\n'; awk '{split($1,a,"_"); print $1 "\t" a[2]*1 "\t" a[3] "\t" $2 "\t" $3}' $1/$1_emmax.ps; } > $1/$1_emmax.ps.qqman
printf "DONE.\n"


# 3. Clump
printf "Extracting index and index-linked SNPs..."
plink2 --bfile ../genotype/$2 --make-founders --clump $1/$1_emmax.ps.qqman --clump-snp-field SNP --clump-field P --clump-p1 0.0000005 --clump-kb 200 --clump-r2 0.50 --clump-p2 0.01 --r2 dprime --clump-range ../../../gene_list_MSU7.txt --clump-range-border 20 --out $1/$1_emmax.ps.qqman.emmax_200kb &> $1/$1_emmax.ps.qqman.emmax_200kb.log
printf "DONE.\n"


if [ -f "$1/$1_emmax.ps.qqman.emmax_200kb.clumped" ]; then
	# run the following only when there are clumped SNPs

	# 4. Creating a list of significant SNPs
	printf "Creating comma-separated list of significant SNPs..."
	sed '/^$/d' $1/$1_emmax.ps.qqman.emmax_200kb.clumped | tail -n+2 | awk '{print  $3 "," $12 }' | sed 's/,NONE//g' | sed 's/(1)//g' | xargs | sed -e 's/ /,/g' > $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps
	printf "DONE.\n"

	printf "Creating a column list of significant SNPs..."
	sed '/^$/d' $1/$1_emmax.ps.qqman.emmax_200kb.clumped | tail -n+2 | awk '{print  $3 "," $12 }' | sed 's/,NONE//g' | sed 's/(1)//g' | xargs | sed -e 's/ /,/g' | sed 's/,/\n/g' |sort |uniq > $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps_list
	printf "DONE.\n"
fi


if [ -f "$1/$1_emmax.ps.qqman.emmax_200kb.clumped.ranges" ]; then
	# 5. Create a list of the associated genes
	printf "Creating a list of the associated or nearby genes..."
	sed '/^$/d' $1/$1_emmax.ps.qqman.emmax_200kb.clumped.ranges | tail -n+2 | awk '{print  $7  }' | sed 's/[][]//g'| xargs | sed 's/ /,/g' | sed 's/,/\n/g' |sort |uniq > $1/$1_emmax.ps.qqman.emmax_200kb.clumped.ranges.genes
	printf "DONE.\n"
fi

### END OF GWAS RUN ###



### POST-PROCESSING ###

if [ -f "$1/$1_emmax.ps.qqman.emmax_200kb.clumped" ]; then
# 6. Given the SNPs from the previous step, create the target list for extracting from the MSU7 reference allele file.
	printf "Creating the list of index and index-linked SNPs (TARGET FILE)..."
	grep -F -f $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps_list ../genotype/$2.bim | sort -k 1n,4n | awk '{print "chr" sprintf("%02d",$1) "\t" $4 "\t" $2 "\t" $5 "\t" $6;}' > $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps.target_list && bgzip -c $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps.target_list > $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps.target_list.gz && tabix -s1 -b2 -e2 $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps.target_list.gz
	printf "DONE.\n"

	# 7. Given the SNP target list from step #1, extract the corresponding regions in the genome to know what are the reference alleles.
	#    NOTE: The plink .BIM file does not guarantee that the A2 allele is the reference allele. What it guarantees is that the A1 is **always**
	#    the minor allele. So, to address that in making sure that we are annotating against the alternative allele during the SNP annotation, we
	#    need to extract that genome regions identified by GWAS and cross reference the .BIM file against this.
	printf "Extracting the corresponding regions in the genome..."
	tabix ../../../reference_allele_MSU7.txt.gz -R $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps.target_list.gz | sort -k 1,2n > $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps.alleles_from_reference
	# Note the use of "-R". It's stands for "rocket". :) The other parameter (-T) stands for "turtle".
	printf "DONE.\n"

	# 8. Column merge the list from the reference genome and the list derived from the gwas.
	# Check if the number of lines matched
	printf "Creating the annotation input file..."
	LINES1=$(wc -l $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps.alleles_from_reference | awk '{print $1}')
	LINES2=$(wc -l $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps.target_list | awk '{print $1}')
	#LINES1=70
	#LINES2=70
	if [ "$LINES1" -ne "$LINES2" ]; then
		FILE2=$(date '+%s%N')

		printf "ERROR: Reference and target lists have unequal number of lines.\nTerminating.\n"
		exit 1
	fi
	# Post-condition: Both files have the same number of lines.

	# Proceed with the rest of the script
	paste $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps.alleles_from_reference $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps.target_list | awk '{if ($2 != $5) { print "Mismatched BP in line number ", NR; exit 1; } else { split($3,a,"="); if (a[2] != $8) print $1 "\t" $2 "\t" $2 "\t" a[2] "\t" $8 "\t" ";A1 in .bim is reference allele."; else print $1 "\t" $2 "\t" $2 "\t" a[2] "\t" $7 "\t" ""; }}' > $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps.avinput


	# IMPORTANT: Possibly reformat the "CHROM" field of the *.avinput file because the rice_msu7 format writes it as 'Chr#'
	#awk '{FS="\t"; sub(/chr/,"",$1); print "Chr" $1*1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps.avinput > $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps.avinput.corrected_chroms
	cat $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps.avinput | sed 's/chr0/Chr/g' | sed 's/chr/Chr/g' > $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps.avinput.corrected_chroms
	printf "DONE.\n"


	# 4. Annotate!
	printf "Annotating..."
	annotate_variation.pl -out $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps.avinput.corrected_chroms -build msu7b $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps.avinput.corrected_chroms ~/bin/annovar/rice_msu7/ -aamatrixfile grantham.matrix &> $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps.avinput.corrected_chroms.log
	printf "DONE.\n"


	# 5. Consolidating the annotations
	Rscript consolidate_annotation_files.R $1

	#if [ $? -eq 0 ]
	#then
	#  echo "Successfully created file"
	#  exit 0
	#else
	#  echo "Could not create file" >&2
	#  exit 1
	#fi


	# 6. Create haploview input files in preparation for generating the LD plots.
	printf "Creating haploview input files..."
	plink2 --bfile ../genotype/$2 --extract $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps_list --recode HV --out $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps_list &> /dev/null
	# generate haploview batch input files
	#ls $1/*snps_list.chr* | paste -s -d ' \n' | awk '{print $2 "\t" $1;}' > $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps_list.haploview_batch_input
	find $1 -name '*.ped' -o -name '*.info' |sort | paste -s -d ' \n' | awk '{print $2 "\t" $1;}' > $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps_list.haploview_batch_input
	printf "DONE.\n"


	# 7. Create the LD plots.
	printf "Generating LD plots..."
	awk '{ped=$1; gsub(".ped","",ped); system("java -jar ~/bin/Haploview.jar -n -pedfile " $1  " -info " $2 " -skipcheck -blockoutput GAB -dprime -out " ped " -svg -png -pairwiseTagging -maxDistance 200"); print $0}' $1/$1_emmax.ps.qqman.emmax_200kb.clumped.snps_list.haploview_batch_input &> /dev/null
	printf "DONE.\n"

	# 8.
	Rscript extract_gene_functions.R $1

   # 9. Create tag SNP diagram
   Rscript create_tagSNP_lanes.R $1 $2

   # 10. Create haploview replacement header
   Rscript create_haploview_replacement_header $1 $2
fi

# 6. Create manhattan and qq plot"
Rscript create_gwas_plots.R $1

# 7. Create Excel summary of results
Rscript create_excel_summary.R $1





