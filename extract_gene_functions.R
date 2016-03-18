# Extract gene function from Nipponbare annotation
# 
# Roslen Anacleto
# 20160307
# 

suppressMessages(library(readxl))
suppressMessages(library(dplyr))


# Get commandline arguments
args<-base::commandArgs(trailingOnly=T)

if (length(args) != 1) {
	   cat("Usage: Rscript extract_gene_function.R <trait>\n")
   quit(save="no", status=1, runLast=T) # non-zero status means an error has occurred.
} else if (!dir.exists(args[1])) {
	   cat(paste0("ERROR: '",args[1],"' does not exist.\n"))
   quit(save="no", status=2, runLast=T)
}

TRAIT <- args[1]
PATH <- TRAIT

cat("Retrieving list of genes...")
genes <- read.table(paste0(PATH,"/",TRAIT,"_emmax.ps.qqman.emmax_200kb.clumped.ranges.genes"), header=F, stringsAsFactors=F)
cat("DONE.\n")

cat("Retrieving functional annotation of genes...")
annot <- readxl::read_excel("../../../japonica_annotation_updated.xlsx")
cat("DONE.\n")


# Get the annotation of the genes from the annotation file
cat("Saving annotated genes to disk...")
res <- dplyr::inner_join(genes, annot, by=c("V1"="GENE_ID"))
colnames(res)[1:2] <- c("LOCUS","CHROM")
write.table(res, file=paste0(PATH,"/",TRAIT,"_emmax.ps.qqman.emmax_200kb.clumped.ranges.genes.annotation"), append=F, quote=F, sep="\t", col.names=T, row.names=F)
cat("DONE.\n")




