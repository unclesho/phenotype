# This script visualizes the tag SNPs lane for the LD block.
#
# Roslen Anacleto
# 20160319
#


### Function Definitions
# Function that returns the SNP bp position
getSnpPosition <- function(snpId, snpLookup) {
   bp<-NULL
   for (i in 1:nrow(snpLookup)) {
      if (snpId==snp_lookup[i,1]) {
         bp <- snp_lookup[i,2]
         break
      }
   }
   return(bp)
}

# Just return the index of the snpId in the snp_lookup list at listitem i
getSnpIndex <- function(snpId, snpLookupDf) {
   bp <- NULL
   for (i in 1:nrow(snpLookupDf)) {
      if (snpLookupDf[i, "V1"]==snpId) {
         bp <- i; break;
      }
   }
   return(bp);
}


# Returns true if SNP is a captured allele
isCapturedAllele <- function(snpId, capturedAlleles) {
   res <- FALSE
   for (i in 1:length(capturedAlleles)) {
      if (snpId == capturedAlleles[i]) {
         res <- TRUE
         break;
      }
   }
   return(res)
}

### End of function definitions


# Load the necessary libraries
#suppressMessages(library(openxlsx))
suppressMessages(library(dplyr))
#suppressMessages(library(plotrix))

# Get commandline arguments
args<-base::commandArgs(trailingOnly=T)

if (length(args) != 2) {
   cat("Usage: Rscript create_tagSNP_lanes.R <trait> <genotype_prefix>\n")
   quit(save="no", status=1, runLast=T) # non-zero status means an error has occurred.
} else if (!dir.exists(args[1])) {
   cat(paste0("ERROR: '",args[1],"' does not exist.\n"))
   quit(save="no", status=2, runLast=T)
} else if (!file.exists(paste0("../genotype/",args[2],".bim"))) {
   cat(paste0("ERROR: '",paste0("../genotype/",args[2],".bim"),"' does not exist.\n"))
   quit(save="no", status=2, runLast=T)
}


# Load the trait and path from the command line paramter
trait=args[1]
path=trait
GENOTYPE <- args[2]
#cat(paste0("Trait: ", trait,"\n"))


# Read into memory the BIM file to gain access to the A1 and A2 alleles.
cat("Reading the genotype data...")
bim <- read.table(paste0("../genotype/",GENOTYPE,".bim"),
                  sep="\t", header=F, stringsAsFactors=F)
cat("DONE.\n")


# Create the SNP lookup table
cat("Sorting the genotype data...")
snp_info_files <- list.files(path=path, pattern="*.info$")
snp_lookup <- vector('list', length(snp_info_files))
for (i in 1:length(snp_info_files)) {
   # Store into the list item the SNPs that are listed under the current chrom file.
   df <- read.table(paste0(path,"/",snp_info_files[i]), sep="\t", stringsAsFactors=F)
   df <- arrange(df, V2) # make sure that they are arranged along the chromosome
   df <- dplyr::inner_join(df, bim[, c(2,5:6)], by=c("V1"="V2"))
   snp_lookup[[i]] <- df
}
cat("DONE.\n")
# Post-condition: At this point, snp_lookup contains list items per chromosome,
# and under each list item a dataframe that contains the SNPs.



cat("Compiling tag SNPs...")
tag_files <- list.files(path=path, pattern="*.TAGS$", include.dirs=T)
for (i in 1:length(tag_files)) {
#for (i in 2:2) {
   # Get the contents of the entire TAGS file into a vector
   all_data <- readLines(paste0(path,"/",tag_files[i]))

   # start PDF device here.

   # No of alleles is in line 1, second word
   snps <- as.numeric(strsplit(all_data[1]," ")[[1]][2]) # numbers are hardcoded
   tag_snps <- as.numeric(strsplit(all_data[3]," ")[[1]][2]) # numbers are also hardcoded
   l_edge<- 5
   BOX_SPACE <- 75
   BOX_LANE <- 10
   r_edge<- l_edge+(snps-1)*BOX_SPACE
   HEIGHT <- tag_snps*BOX_SPACE # allocate a plot space of 75 pixels per lane
   PLOT_MARGIN <- 300

   # Check if the info and tag files match for index i.
   if (nrow(snp_lookup[[i]]) != snps) {
      cat(paste0("ERROR: SNPs in ", snp_info_files[i], " does not match ",
                 tag_files[i], ".\n"))
      break; # don't proceed any further
   }

   pdf(file=paste0(path,"/",tag_files[i],"_tagSNPs.pdf"))
   par(mar=c(1,1,1,1))
   #plot(NULL,axes=F,ann=FALSE,xlim=c(l_edge-200, r_edge+200),ylim=c(1, HEIGHT+100))
   #plot(NULL,axes=F,ann=FALSE,xlim=c(l_edge-PLOT_MARGIN, r_edge+PLOT_MARGIN), ylim=c(1, r_edge+PLOT_MARGIN))
   plot(NULL,axes=F,ann=FALSE,xlim=c(1, r_edge+BOX_SPACE), ylim=c(1, r_edge+BOX_SPACE))

   lane <- 1
   offset <- 10 # from the bottom of the plot
   lane_space <-75
   for (j in 1:tag_snps) {
      res <- strsplit(all_data[snps+7+(j-1)], split="\t")
      # Note: res[[1]][1] - holds the tag SNP
      # res[[1]][2] - holds the comma-separated list of captured SNPs

      # Ignore the tag SNP; Use only the captured SNPs
      res2 <- strsplit(res[[1]][2], ",")
      res2[[1]] <- sort(res2[[1]]) # to prevent complications later

      # Get captured alleles
      captured_alleles <- res2[[1]]

      for (k in 1:snps) {
         xbot <-l_edge+(k-1)*BOX_SPACE
         ybot <-offset+(j-1)*BOX_SPACE
         rect(xbot, ybot, xbot+(BOX_SPACE-BOX_LANE), ybot+(BOX_SPACE-BOX_LANE),
              col=ifelse(snp_lookup[[i]][k,1]==res[[1]][1],"red",
                         ifelse(isCapturedAllele(snp_lookup[[i]][k,1], captured_alleles),
                                "green","white"))
              )
      }
   }

   # Plot the snps as boxes including their A2/A1 inside the box
   for(k in 1:snps) {
      xbot <-l_edge+(k-1)*BOX_SPACE
      ybot <-offset+tag_snps*BOX_SPACE
      text(xbot+(BOX_SPACE-BOX_LANE)/2, ybot+3*(BOX_LANE+BOX_SPACE),
           labels=paste0(snp_lookup[[i]][k,1]," (",
                         snp_lookup[[i]][k,4],"/",snp_lookup[[i]][k,3], ")"),
           srt=90, cex=0.6)
   }

   mtext(paste0(trait," - ", tag_files[i]), side=3)

   res <- dev.off() # close the PDF device
}
cat("DONE.\n")





