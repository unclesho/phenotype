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
suppressMessages(library(svglite))
#suppressMessages(library(plotrix))

# Get commandline arguments
args<-base::commandArgs(trailingOnly=T)

if (length(args) != 2) {
   cat("Usage: Rscript create_haploview_replacement_header.R <trait> <genotype_prefix>\n")
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


cat(paste0("Trait: ",trait,"\nGenotype: ",GENOTYPE,"\n"))

# Read the GWAS results data
cat("Reading GWAS results...")
gwas <- read.table(paste0(path,"/",trait,"_emmax.ps.qqman"),
                   header=T, sep="\t", stringsAsFactors=F)
cat("DONE.\n")

# Read into memory the BIM file to gain access to the A1 and A2 alleles.
cat("Reading the genotype data...")
bim <- read.table(paste0("../genotype/",GENOTYPE,".bim"),
                  sep="\t", header=F, stringsAsFactors=F)
cat("DONE.\n")


# Join the gwas and bim files
cat("Joining the genotype data and GWAS results...")
bim_gwas <- dplyr::inner_join(bim, gwas[,c(1,4:5)], by=c("V2"="SNP"))
#.. then remove individual variables to free memory
rm(list=c("gwas","bim"))
cat("DONE.\n")


# Create the SNP lookup table
cat("Sorting the genotype data")
snp_info_files <- list.files(path=path, pattern="*.info$")
snp_lookup <- vector('list', length(snp_info_files))
for (i in 1:length(snp_info_files)) {
   cat(paste0("...",i))
   # Store into the list item the SNPs that are listed under the current chrom file.
   df <- read.table(paste0(path,"/",snp_info_files[i]), sep="\t", stringsAsFactors=F)
   df <- arrange(df, V2) # make sure that they are arranged along the chromosome

   # These are viable options to read in large files. The problem is that they're not as
   # convenient as the read.table option.
   # df <- system(paste0("sort -k2n,2 ", path,"/",snp_info_files[i]), intern=T)
   # df<-system2(command="sort", args=c("-k2n,2",paste0(path,"/",snp_info_files[i])), stdout=T)

   df <- dplyr::inner_join(df, bim_gwas[, c(2,5:8)], by=c("V1"="V2"))
   snp_lookup[[i]] <- df
}
rm(list=c("df","i"))
cat("...DONE.\n")
# Post-condition: At this point, snp_lookup contains list items per chromosome,
# and under each list item a dataframe that contains the SNPs.


# Read the block SNPs
cat("Reading LD blocks")
block_files <- list.files(path=path, pattern="*.GABRIELblocks$")
block_snp <- vector('list', length(block_files)) # represents the GABRIEL block files
for (i in 1:length(block_files)) {
   cat(paste0("...",i))
   # Read in the blocks
   #df <- readLines(paste0(path,"/",block_files[i]))
   blocks <- system(paste0("grep 'BLOCK' ",path,"/",block_files[i]), intern=T)
   block_snp[[i]] <- vector('list', length(blocks)) # represents the defined blocks within the file
   for (j in 1:length(blocks)) {
      row_vec <- unlist(strsplit(blocks[j],"  ")) # convert to vector
      # Note: row_vec[1] is the "BLOCK #."
      #       row_vec[2] is the "MARKERS: #[ #]*
      block_snp[[i]][[j]] <- unlist(strsplit(substring(row_vec[2],10), split=" "))
   }
}
rm(list=c("blocks","row_vec","i"))
cat("...DONE.\n")


cat("Compiling tag SNPs")
tag_files <- list.files(path=path, pattern="*.TAGS$", include.dirs=T)
tag_snp <- vector('list', length(tag_files)) # represents the GABRIEL block files
for (i in 1:length(tag_files)) {
   cat(paste0("...",i))
   # Get the contents of the entire TAGS file into a vector
   all_data <- readLines(paste0(path,"/",tag_files[i]))

   # No of alleles is in line 1, second word
   snps <- as.numeric(strsplit(all_data[1]," ")[[1]][2]) # numbers are hardcoded
   ts <- as.numeric(strsplit(all_data[3]," ")[[1]][2]) # numbers are also hardcoded

   # Check if the info and tag files match for index i.
   if (nrow(snp_lookup[[i]]) != snps) {
      cat(paste0("ERROR: SNPs in ", snp_info_files[i], " does not match ",
                 tag_files[i], ".\n"))
      break; # don't proceed any further
   }

   for (j in 1:ts) {
      res <- unlist(strsplit(all_data[snps+7+(j-1)], split="\t"))
      # Note: res[[1]][1] - holds the tag SNP
      # res[[1]][2] - holds the comma-separated list of captured SNPs

      # Ignore the tag SNP; Use only the captured SNPs
      res2 <- strsplit(res[[1]][2], ",")
      res2[[1]] <- sort(res2[[1]]) # to prevent complications later

      tag_snp[[i]][j] <- res[1]
   }
   tag_snp[[i]] <- sort(tag_snp[[i]])
}
rm(list=c("all_data","snps","ts","res","res2","i","j"))
cat("...DONE.\n")


# Read the annotation file
cat("Reading the MSU7 annotation...")
gff <- read.table("../../../../msu7/all.dir/all.gff3",
                  comment.char="#", sep="\t", header=F, stringsAsFactors=F)
cat("DONE.\n")


# At this point:
# bim_gwas - holds the genotype data with markers associated with effect sizes and p-values
# snp_lookup - a list whose length is same as # of *.info; each item is a df of snps.
# block_snp - contains the SNPs identified in each Gabriel block
# tag_snp - contains the tag SNPs per chromosome
# gff - holds the complete MSU7 annotation


# Generates PDFs of the figures per list item stored in snp_lookup.
cat("Generating the plot")
for (i in 1:length(snp_lookup)) {
   cat(paste0("...",i))
   # Get chromosome using the first SNP in the list.
   chrom <- as.numeric(unlist(strsplit(snp_lookup[[i]][,"V1"][1], split="_"))[2])
   if (chrom<1) {
      cat("Error: Invalid chromosome number. SNP_ID must be 'snp_xx_xxxxxx'.")
      break;
   }
   chrom <- paste0("Chr",  as.numeric(unlist(strsplit(snp_lookup[[i]][,"V1"][1], split="_"))[2]))

   # Set the left and right edges of the plot
   l_edge <- min(snp_lookup[[i]][,"V2"])
   r_edge <- max(snp_lookup[[i]][,"V2"])
   HEIGHT <- 2000

   # compute for the aspect ratio
   aspect_ratio <- diff(c(l_edge, r_edge))/HEIGHT # how many x dots per y dot
   CHROM_POS <- 300 * aspect_ratio


   suppressWarnings(svglite::svglite(file=paste0(path,"/",snp_info_files[i],"_haploview_lanes.svg"),
                                     width=8, height=8))
   par(mar=c(1,1,1,1))
   par(fig=c(0,1,0,0.8), new=TRUE) # create the main figure at the lower 80% portion of figure
   plot(NULL,axes=F,ann=FALSE, asp=1,
        xlim=c(l_edge, r_edge), # set the width of the plot
        ylim=c(1, max(HEIGHT, (r_edge-l_edge)+1))) # set the height of the plot
   lines(x=c(l_edge, r_edge),
         y=c(CHROM_POS, CHROM_POS),
         lwd=2,
         col="gray")

   # Extract and plot the genes in the chromosome
   genes <- gff %>%
      filter(V1==chrom, V3=="gene", ((V4>=l_edge & V4<=r_edge) | (V5>=l_edge & V5<=r_edge)))
   #..then plot the genes
   for (j in 1:nrow(genes)) {
      # get the locus_id
      locus_id <- unlist(strsplit(unlist(strsplit(genes[j,9], split=";"))[1], split="="))[2]
      #..then display it
      text(x = l_edge+(r_edge-l_edge)/(nrow(genes)-1)*(j-1),
           y = CHROM_POS-(100*aspect_ratio),
           labels=locus_id, srt=90,
           cex=0.3, pos=2, offset=0)
      lines(lwd=0.5, #lend=1,
            x=c(l_edge+(r_edge-l_edge)/(nrow(genes)-1)*(j-1),
                l_edge+(r_edge-l_edge)/(nrow(genes)-1)*(j-1)),
            y=c(CHROM_POS-(90*aspect_ratio),
                CHROM_POS-(95*aspect_ratio)))

      if (genes[j,"V4"] < l_edge) {
         # Draw the top, right, and bottom lines only
         rect(xleft=l_edge, ybottom=CHROM_POS-(10*aspect_ratio),
              xright=genes[j,"V5"], ytop=CHROM_POS+(10*aspect_ratio),
              col="red", lwd=0.5, border=NA) # no border
         # Draw top line
         lines(lwd=0.5, x=c(l_edge, genes[j,"V5"]),
               y=c(CHROM_POS+(10*aspect_ratio), CHROM_POS+(10*aspect_ratio)))
         # Draw right line
         lines(lwd=0.5, x=c(genes[j,"V5"], genes[j,"V5"]),
               y=c(CHROM_POS-(10*aspect_ratio), CHROM_POS+(10*aspect_ratio)))
         # Draw bottom line
         lines(lwd=0.5, x=c(l_edge, genes[j,"V5"]),
               y=c(CHROM_POS-(10*aspect_ratio), CHROM_POS-(10*aspect_ratio)))
         # Draw the label line
         lines(lwd=0.5, x=c(genes[j,"V5"], genes[j,"V5"]),
               y=c(CHROM_POS-(15*aspect_ratio), CHROM_POS-(20*aspect_ratio)))
         lines(lwd=0.5, col="darkgray",
               x=c(genes[j,"V5"], l_edge+(r_edge-l_edge)/(nrow(genes)-1)*(j-1)),
               y=c(CHROM_POS-(20*aspect_ratio), CHROM_POS-(90*aspect_ratio)))

      } else if (genes[j, "V5"]> r_edge) {
         # Draw the top, right, and bottom lines only
         rect(xleft=genes[j,"V4"], ybottom=CHROM_POS-(10*aspect_ratio),
              xright=r_edge, ytop=CHROM_POS+(10*aspect_ratio),
              col="red", lwd=0.5, border=NA) # no border
         # Draw top line
         lines(lwd=0.5, x=c(genes[j,"V4"], r_edge),
               y=c(CHROM_POS+(10*aspect_ratio), CHROM_POS+(10*aspect_ratio)))
         # Draw left line
         lines(lwd=0.5, x=c(genes[j,"V4"], genes[j,"V4"]),
               y=c(CHROM_POS-(10*aspect_ratio), CHROM_POS+(10*aspect_ratio)))
         # Draw bottom line
         lines(lwd=0.5, x=c(genes[j,"V4"], r_edge),
               y=c(CHROM_POS-(10*aspect_ratio), CHROM_POS-(10*aspect_ratio)))
         # Draw the label line
         lines(lwd=0.5, x=c(genes[j,"V4"], genes[j,"V4"]),
               y=c(CHROM_POS-(15*aspect_ratio), CHROM_POS-(20*aspect_ratio)))
         lines(lwd=0.5, col="darkgray",
               x=c(genes[j,"V4"], l_edge+(r_edge-l_edge)/(nrow(genes)-1)*(j-1)),
               y=c(CHROM_POS-(20*aspect_ratio), CHROM_POS-(90*aspect_ratio)))


      } else {
         rect(xleft=genes[j,"V4"], ybottom=CHROM_POS-(10*aspect_ratio),
              xright=genes[j,"V5"], ytop=CHROM_POS+(10*aspect_ratio),
              col="red", lwd=0.5)
         # Draw the label line
         lines(lwd=0.5, x=c(genes[j,"V4"], genes[j,"V4"]),
               y=c(CHROM_POS-(15*aspect_ratio), CHROM_POS-(20*aspect_ratio)))
         lines(lwd=0.5, col="darkgray",
               x=c(genes[j,"V4"], l_edge+(r_edge-l_edge)/(nrow(genes)-1)*(j-1)),
               y=c(CHROM_POS-(20*aspect_ratio), CHROM_POS-(90*aspect_ratio)))
      }
   }

   # It's time to mark the SNPs
   for (j in 1:nrow(snp_lookup[[i]])) {
      text(x = l_edge+(r_edge-l_edge)/(nrow(snp_lookup[[i]])-1)*(j-1),
           y = CHROM_POS+(100*aspect_ratio),
           labels=paste0(snp_lookup[[i]][j,1]," (",
                         snp_lookup[[i]][j,4], "/",
                         snp_lookup[[i]][j,3], ")"),
           srt=90, cex=0.3, pos=4, offset=0,
           col=ifelse(snp_lookup[[i]][j,1] %in% tag_snp[[i]],"red","black"),
           font=ifelse(snp_lookup[[i]][j,1] %in% tag_snp[[i]],2,1))
      lines(lwd=0.5, #lend=1,
            x=c(l_edge+(r_edge-l_edge)/(nrow(snp_lookup[[i]])-1)*(j-1),
                l_edge+(r_edge-l_edge)/(nrow(snp_lookup[[i]])-1)*(j-1)),
            y=c(CHROM_POS+(90*aspect_ratio),
                CHROM_POS+(95*aspect_ratio)))

      lines(lwd=0.5, col="darkgray", #lend=1,
            x=c(snp_lookup[[i]][j,2],
                l_edge+(r_edge-l_edge)/(nrow(snp_lookup[[i]])-1)*(j-1)),
            y=c(CHROM_POS+(20*aspect_ratio),
                CHROM_POS+(90*aspect_ratio)))

      lines(lwd=0.5, #lend=1,
            x=c(snp_lookup[[i]][j,2], snp_lookup[[i]][j,2]),
            y=c(CHROM_POS+(15*aspect_ratio),
                CHROM_POS+(20*aspect_ratio)))
   }
   #axis(1);box();

   # create the plot of p-values and allele effects at the top 20%
   par(mar=c(4,4,1,1)+0.1)
   par(fig=c(0,1,.8,1), new=TRUE)
   #boxplot(mtcars$wt, horizontal=TRUE, axes=FALSE)
   snps <- snp_lookup[[i]]
   snps$xpos <- seq(from=snps[1,2], to=snps[nrow(snps),2], length.out=nrow(snps))
   snps$logp <- -log10(snps$P)
   snps$label <- paste0(snps$V1," (",snps$V6,"/",snps$V5,")")
   snps$absBETA <- abs(snps$BETA)
   snps$l_width <- ceiling(1+(snps$absBETA-min(snps$absBETA))*(9/diff(range(snps$absBETA))))
   plot(snps[, c("xpos","logp")],
        ylim=c(0,ceiling(max(snps$logp))),
        xlim=c(min(snps$xpos),max(snps$xpos)),
        type="n", axes=F, xlab="", ylab="")
   for (j in 1:nrow(snps)) {
      lines(x=c(snps[j,"xpos"], snps[j,"xpos"]),
            y=c(0,snps[j, "logp"]),
            lwd=snps[j, "l_width"], lend=1,
            col=ifelse(snps[j,"BETA"]<0,"red","black"))
   }
   axis(1, las=2, cex.axis=0.3, mgp=c(3, .7, 0), at=snps$xpos, labels=snps$label)
   axis(2, las=2, cex.axis=0.3, mgp=c(3, .7, 0))
   mtext(text=expression(-log[10](p)), side=2, line=1, cex=0.3)
   res <- dev.off()
}
cat("...DONE.\n")
rm(list=c("chrom", "genes", "l_edge", "r_edge","i","snps"))


