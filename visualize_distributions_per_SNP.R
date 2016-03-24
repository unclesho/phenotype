# This script is used to visualize the distribution of the phenotype per SNP for each trait.
#
# Roslen Anacleto
# 20160310
#
# Requirements:
# 1. .raw file must be in the genotype directory. It is created by the "--recode A" in plink2.
# 2. .bim file must be in the genotype directory.
# 3. phenotype data - the original set saved as tab-separated values from the xlsx file.
#


# Get commandline arguments
args<-base::commandArgs(trailingOnly=T)

if (length(args) != 2) {
   cat("Usage: Rscript visualize_distributions_per_SNP.R <genotype_prefix> <phenotype_datasheet>\n")
   quit(save="no", status=1, runLast=T) # non-zero status means an error has occurred.
} else if (!file.exists(paste0("../genotype/",args[1],".bim"))) {
   cat(paste0("ERROR: '",paste0("../genotype/",args[1],".bim"),"' does not exist.\n"))
   quit(save="no", status=2, runLast=T)
} else if (!file.exists(paste0("../genotype/",args[1],".raw"))) {
   cat(paste0("ERROR: '",paste0("../genotype/",args[1],".raw"),"' does not exist.\n"))
   quit(save="no", status=2, runLast=T)
} else if (!file.exists(args[2])) {
   cat(paste0("ERROR: '",paste0("../genotype/",args[2],".bim"),"' does not exist.\n"))
   quit(save="no", status=2, runLast=T)
}


# Load the libraries
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(broom))
suppressMessages(library(openxlsx))
suppressMessages(library(svglite))
#Sys.setenv(R_ZIPCMD="/usr/bin/zip")


##### DEFINE ALL FUNCTIONS HERE

getAllele <- function (snpIndex, bimDf, alleleFormat=c("A1A2","A1","A2", "expression")) {
   if (alleleFormat=="A1A2") {
      return (paste(bimDf[snpIndex,c(5,6)], collapse=""))
   } else if (alleleFormat=="A1") {
      return (bimDf[snpIndex,c(5)])
   } else if (alleleFormat=="A2") {
      return (bimDf[snpIndex,c(6)])
   } else if (alleleFormat=="expression") {
      return (paste0("~frac(",bimDf[snpIndex,c(6)],",",bimDf[snpIndex,c(5)],")~"))
   }
}

getAlleleExpression <- function (snpIndex, alleleDose, bimDf) {
   if (alleleDose=="0") {
      return (bimDf[snpIndex,c(6)])
   } else if(alleleDose=="1") {
      return (paste0("~frac(",bimDf[snpIndex,c(6)],",",bimDf[snpIndex,c(5)],")~"))
   } else {
      # case when alleleDose==2
      return (bimDf[snpIndex,c(5)])
   }
}

getHaplotypeLabels <- function(hapCodes, bimDf){
   haplotypeLabels<-array(dim=length(hapCodes))
   for (i in 1:length(haplotypeLabels)) {
      alleleDosage <- unlist(strsplit(hapCodes[i], split=""))
      hapString <- ""
      for (j in 1:length(alleleDosage)) {
         hapString <- paste0(hapString, getAlleleExpression(j, alleleDose=alleleDosage[j], bimDf))
      }
      haplotypeLabels[i] <- parse(text=hapString)
   }

   return(haplotypeLabels)
}


getGenotypeLabels <- function (snpIndex, bimDf) {
   genotypeLabels<-array(3) # it is always 3 A1A1, A1A2, A2A2
   genotypeLabels[1] <- parse(text=paste0(bimDf[snpIndex,c(5)],bimDf[snpIndex,c(5)]))
   genotypeLabels[2] <- parse(text=paste0(bimDf[snpIndex,c(5)],bimDf[snpIndex,c(6)]))
   genotypeLabels[3] <- parse(text=paste0(bimDf[snpIndex,c(6)],bimDf[snpIndex,c(6)]))
   return(genotypeLabels)
}

##### END OF FUNCTION DEFINITION

# MAIN

# Load the trait and path from the command line paramter
GENOTYPE <- args[1]
PHENOTYPE <- args[2]


# Load the BIM file
cat("Loading the main genotype file...")
bim <- read.table(paste0("../genotype/", GENOTYPE,".bim"), header=F, stringsAsFactors=F)
cat("DONE.\n")
# V2 - SNP_ID
# V5 - A1
# V6 - A2

# Load the phenotype data from the original matrix extracted from the source excel file
cat("Loading the phenotype file...")
pheno <- read.table(paste0(PHENOTYPE), header=T, sep="\t", quote="", stringsAsFactors=F)
cat("DONE.\n")

# Remove the FID column
pheno <- pheno[,-1]

# {POSTCONDITION: At this point "pheno" contains a dataframe with an IID column and
# individual columns for the traits.}


# Plot the SNP distribution per chromosome because that's how the haploview file
# is structured.
# NOTE: This current version will only work if there is one haplotype block shown
# per chromosome.


# List down the trait names
TRAITS <- colnames(pheno[,-c(1)]) # Tricky bit: No FID already. Col#1 is the IID field.


# Go trait by trait in plotting the SNP distribution.
# {PRECONDITION: "pheno" still contains the IID column so remove it first from the iterate list.}
for (t in 1:length(TRAITS)) {
   cat(paste0("Trait: ", TRAITS[t],"\n"))

   # set the trait name as path
   path <- TRAITS[t]

   test_files <- list.files(path=path, pattern="*.TESTS$", include.dirs=T)
   # Go chromosome by chromosome
   for (f in 1:length(test_files)) {
      cat(paste0("   ", test_files[f],"\n"))

      # create the raw files
      system(paste0("plink2",
                    " --bfile ../genotype/", GENOTYPE,
                    " --recode A",
                    " --extract ", path, "/", test_files[f],
                    " --out ../genotype/", test_files[f]))

      # Read the genotype data
      cat("Loading the genotype table...")
      geno <- read.table(paste0("../genotype/", test_files[f],".raw"),
                         header=T, quote="", stringsAsFactors=F)
      cat("DONE.\n")

      # Remove the unnecessary columns in genotype file
      geno2 <- geno[,-c(1,3:6)]

      # Serialize geno2 in preparation for the one-way ANOVA
      cat("Serializing the genotype data...")
      geno3 <- geno2 %>% gather(snp,genotype,2:length(colnames(geno2)))
      cat("DONE.\n")

      # REMOVE THE NAs at this point!
      geno3 <- na.omit(geno3)

      # Bring the phenotype into the serialized data
      cat("Appending the phenotype to serialized genotype data...")
      geno3_pheno <- dplyr::inner_join(geno3, pheno, by=c("IID"="iid"))
      #..then convert the snp and genotype columns as factors
      geno3_pheno$snp <- as.factor(geno3_pheno$snp)
      geno3_pheno$genotype <- as.factor(geno3_pheno$genotype)
      cat("DONE.\n")


      # Now, modify geno2 such that those with NA will be represented as "N"
      geno2[is.na(geno2)]<-"N" # positions with NA


      # Create haplotypes
      hap <- as.data.frame(cbind(geno2[,1],apply(geno2[,-1], 1, paste, collapse="")),
                           stringsAsFactors=F)
      colnames(hap) <- c("IID","hap")


      # Get the frequency of each haplotype
      hap_df <- as.data.frame(table(hap$hap)) %>% arrange(desc(Freq))
      colnames(hap_df)[1] <- "Haplotype"

      # Now create the hap_pheno table --- UNFILTERED!!!
      hap_pheno <- dplyr::inner_join(hap, pheno[,c(1,(t+1))], by=c("IID"="iid"))
      colnames(hap_pheno)[2] <- "haplotype"
      write.table(hap_pheno,
                  file=paste0(TRAITS[t],"/",test_files[f],"_haplotype_table.txt"),
                  quote=F, sep="\t", row.names=F, col.names=T, append=F)

      # Create the boxplots per haplotype for TRAITS[t]
      # pdf(file=paste0(TRAITS[t],"/", TRAITS[t],
      #                 "_", test_files[f],
      #                 "_haplotype_boxplot.pdf"),
      #     width=7, height=5)
      pdf(file=paste0(TRAITS[t],"/", TRAITS[t],
                      "_", test_files[f],
                      "_haplotype_boxplot.pdf"),
          width=7, height=5)

      # IMPORTANT: Decide at this point whether or not to show the haplotypes with missing
      # genotype data.
      hap_pheno_table <- as.data.frame(table(hap_pheno$haplotype)) %>% arrange(desc(Freq))
      hap_pheno_filtered <- hap_pheno[-(grep("N", hap_pheno$haplotype)),]


      # Just assign a friendly name to the filter criteria to simplify the appearance
      good_haps <- as.character(hap_pheno_table[hap_pheno_table[,2] >=
                                                   ceiling(nrow(hap_pheno)*0.01),][,1])
      hap_pheno_filtered <- hap_pheno_filtered[hap_pheno_filtered[,2] %in% good_haps,]
      # POSTCONDITIONS:
      # 1. ALL haplotypes with missing calls are eliminated
      # 2. Only haplotypes with with at least 1% occurrence are included in the boxplot.

      # Ready to plot the boxplot for TRAIT[t]!!!
      res<-boxplot(as.formula(paste0(TRAITS[t],"~haplotype")),
                   hap_pheno_filtered, # Already filtered!!!
                   las=2, xaxt='n', cex.axis=0.7, #axes=F,
                   ylab=TRAITS[t], main=paste0("Haplotypes for ",TRAITS[t]),
                   par(mar = c(10, 4, 4, 1) + 0.1))
      axis(1, at=c(1:length(res$names)), labels=getHaplotypeLabels(res$names, bim),
           las=2, cex.axis=0.7)

      # Show the frequency of each haplotype
      mtext(paste("(n=",res$n,")",sep=""), at=c(1:length(res$names)), side=1, line=8, cex=0.7)

      res<-dev.off() # Just prevent the annoying messsage when a device is closed.

      # POSTCONDITION: the haplotype for the trait has already been generated at this point.


      # Create the individual boxplot per genotype per SNP per trait.
      SNPS <- unique(geno3_pheno$snp)
      for (j in 1:length(SNPS)) {
         pdf(file=paste0(path, "/",test_files[f],"_",
                         paste(unlist(strsplit(as.character(SNPS[j]),
                                               "_"))[1:3],collapse="_"),
                         "_boxplot.pdf"), width=5, height=5)

         # Start of replacement code that shows A1A1, A1A2, A2A2 coding
         res<-boxplot(as.formula(paste0(TRAITS[t]," ~ genotype")),
                      geno3_pheno %>% filter(snp==SNPS[j]),
                      las=1, xaxt='n',
                      # xlab=paste0("Dosage of A1 allele (",
                      #             unlist(strsplit(as.character(SNPS[j]),"_"))[4],")"),
                      xlab="",
                      ylab=paste0(TRAITS[t]),
                      main=paste0(paste(unlist(strsplit(as.character(SNPS[j]),"_"))[1:3],
                                        collapse="_"), " on ",TRAITS[t]))
         axis(1, at=c(1:length(res$names)), labels=getGenotypeLabels(j, bim), las=1)
         # End of replacement code

         # Show the frequency of each haplotype
         mtext(paste("(n=",res$n,")",sep=""), at=c(1:length(res$names)), side=1, line=2, cex=0.7)

         res<-dev.off() # suppress the message from dev.off()
      } # for (j...)
   } # for (f...)
} # for (t...)
