# This script creates the excel summary of the gwas run per trait.
#
# Roslen Anacleto
# 20160309
#


# Sys.getenv("R_ZIPCMD", "zip")
Sys.setenv(R_ZIPCMD="/usr/bin/zip")   # without this openxlsx will cry


# Load the necessary libraries
suppressMessages(library(openxlsx))
suppressMessages(library(dplyr))


# Get commandline arguments
args<-base::commandArgs(trailingOnly=T)

if (length(args) != 1) {
   cat("Usage: Rscript create_excel_summary.R <trait>\n")
   quit(save="no", status=1, runLast=T) # non-zero status means an error has occurred.
} else if (!dir.exists(args[1])) {
   cat(paste0("ERROR: '",args[1],"' does not exist.\n"))
   quit(save="no", status=2, runLast=T)
}


# Load the trait and path from the command line paramter
trait=args[1]
path=trait
#cat(paste0("Trait: ", trait,"\n"))


# Set the working directory
#setwd("/Users/ranacleto/Google Drive/02 - Stat analyses/31 - Genomics/50 - GWAS/15-HDRAv3/run/320PHY_2013WS_SEC/phenotype")


### 1. Create the workbook and set font size to 10
#cat("Initializing...")
wb <- createWorkbook()
modifyBaseFont(wb, fontSize=10)


### 2. Create the worksheets under this workbook

# am1/am1_emmax.ps.qqman.emmax_200kb.clumped
addWorksheet(wb, sheetName="SNPs and linked SNPs")

# am1/am1_emmax.ps.qqman.emmax_200kb.clumped.ranges
addWorksheet(wb, sheetName="SNPs and associated genes")

# am1/am1_emmax.ps.qqman.emmax_200kb.clumped.snps.avinput.corrected_chroms.variant_function.consolidated
# Without header
addWorksheet(wb, sheetName="Annotation of SNPs")

# am1/am1_emmax.ps.qqman.emmax_200kb.clumped.ranges.genes.annotation
addWorksheet(wb, sheetName="annotation of genes")
#cat("DONE.\n")

# create the tag SNPs
addWorksheet(wb, sheetName="tag SNPs")


### 3. Write the data into the worksheets
# sheet 1
cat("Compiling SNPs and linked SNPs...")
data <- read.table(paste0(path,"/",trait,"_emmax.ps.qqman.emmax_200kb.clumped"), header=T, stringsAsFactors=F)
writeData(wb=wb, sheet="SNPs and linked SNPs", x=data, colNames=TRUE, rowNames=FALSE)
setColWidths(wb, sheet="SNPs and linked SNPs", cols=1:(dim(data)[2]), widths="auto")
cat("DONE.\n")

# sheet 2
cat("Compiling SNPs and associated genes...")
data <- read.table(paste0(path,"/",trait,"_emmax.ps.qqman.emmax_200kb.clumped.ranges"), header=T, stringsAsFactors=F)
writeData(wb=wb, sheet="SNPs and associated genes", x=data, colNames=TRUE, rowNames=FALSE)
setColWidths(wb, sheet="SNPs and associated genes", cols=1:(dim(data)[2]), widths="auto")
cat("DONE.\n")


# sheet 3
cat("Compiling SNP annotation...")
data <- read.table(paste0(path,"/",trait,"_emmax.ps.qqman.emmax_200kb.clumped.snps.avinput.corrected_chroms.variant_function.consolidated"),
                   header=F, sep="\t", stringsAsFactors=F)
colnames(data) <- c("SNP_location","LOCUS","CHROM","start","end","Ref","Alt","Notes","Line_no","Type_of_AA_change","Transcript_and_AA_change","Notes")
writeData(wb=wb, sheet="Annotation of SNPs", x=data, colNames=TRUE, rowNames=FALSE)
setColWidths(wb, sheet="Annotation of SNPs", cols=1:(dim(data)[2]), widths="auto")
cat("DONE.\n")


# sheet 4
# NOTE: Putting 'quote=""' ensures that any quotation marks ", or ' in any field is ignored.
cat("Compiling gene annotations...")
data <- read.table(paste0(path,"/",trait,"_emmax.ps.qqman.emmax_200kb.clumped.ranges.genes.annotation"), header=T, sep="\t", stringsAsFactors=F, quote="")
writeData(wb=wb, sheet="annotation of genes", x=data, colNames=TRUE, rowNames=FALSE)
setColWidths(wb, sheet="annotation of genes", cols=1:(dim(data)[2]), widths="auto")
cat("DONE.\n")


# sheet 5
# save the TAG SNPs here

# Save the workbook
cat("Saving the summary...")
saveWorkbook(wb, paste0(path,"/",trait,"_emmax.ps.qqman.emmax_200kb.summary.xlsx"), overwrite=TRUE)
cat("DONE.\n")



