# This script is used to visualize the distribution of the phenotype per SNP for each trait.
#
# Roslen Anacleto
# 20160310
#
# Requirements:
# 1. .raw file must be in the genotype directory. It is created by the "--recode A" in plink2.
# 2. .bim file must be in the genotype directory.
# 3. phenotype data
#


GENOTYPE <- "320PHY_2013WS_SEC_LOC_Os06g04200.3.mind01"
PHENOTYPE <- "320PHY_2013WS_SEC_LOC_Os06g04200.3"



suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(broom))
suppressMessages(library(openxlsx))
suppressMessages(library(svglite))
#Sys.setenv(R_ZIPCMD="/usr/bin/zip")

# Point the stats directory
setwd("/Users/ranacleto/Google Drive/02 - Stat analyses/31 - Genomics/50 - GWAS/15-HDRAv3/run/320PHY_2013WS_SEC_LOC_Os06g04200.3/stats")
#setwd("/home/roslen/scratch2/gwas_hdra_v3/run/320PHY_2013WS_SEC_LOC_Os06g04200.3/stats")

# Read the genotype data
cat("Loading the genotype table...")
geno <- read.table(paste0("../genotype/", GENOTYPE,".raw"), header=T, quote="", stringsAsFactors=F)
cat("DONE.\n")


# Load the phenotype data from the original matrix extracted from the source excel file
cat("Loading the phenotype file...")
pheno <- read.table(paste0("../phenotype/", PHENOTYPE,".txt"), header=T, sep="\t", quote="", stringsAsFactors=F)
cat("DONE.\n")

# Remove the FID column
pheno <- pheno[,-1]


# Join the genotype and phenotype
#geno_pheno <- dplyr::inner_join(data, pheno, by=c("IID"="iid"))

# Remove the unnecessary columns in genotype file
geno2 <- geno[,-c(1,3:6)]


# Load the BIM file
bim <- read.table(paste0("../genotype/", GENOTYPE,".bim"), header=F, stringsAsFactors=F)
# V2 - SNP_ID
# V5 - A1
# V6 - A2



# Create haplotypes
hap <- as.data.frame(cbind(geno2[,1],apply(geno2[,-1], 1, paste, collapse="")), stringsAsFactors=F)
colnames(hap) <- c("IID","hap")


# Get the frequency of each haplotype
hap_df <- as.data.frame(table(hap$hap)) %>% arrange(desc(Freq))
colnames(hap_df)[1] <- "Haplotype"

# Sort the haplotypes according to frequency
#hap_df %>% arrange(desc(Freq))

hap_pheno <- dplyr::inner_join(hap, pheno, by=c("IID"="iid"))
TRAITS <- colnames(pheno[,-c(1)]) # starts at column 4 in geno3_pheno


# OLD BUT WORKING VERSION:
# Create the boxplots per haplotype per trait
# for (i in 1:length(TRAITS)) {
#     cat(paste0("Processing ",TRAITS[i],"..."))
#     pdf(file=paste0("images/",TRAITS[i],"_haplotype_boxplot.pdf"), width=7, height=5)
#     boxplot(as.formula(paste0(TRAITS[i],"~hap")),hap_pheno,
#             las=2, ylab=TRAITS[i], main=paste0("Haplotypes for ",TRAITS[i]),
#             par(mar = c(8, 4, 4, 1) + 0.1))
#     mtext(paste("(n=",table(hap_pheno$hap),")",sep=""),
#           at=c(1:15), side=1, line=5, cex=0.7)
#     mtext("Haplotypes based on minor allele dosage", side=1, line=6.5)
#     res<-dev.off()
#     cat("DONE.\n")
# }



 labs <- c("AT~frac(T,C)~G","GGAA","TTAA","A~over(A,C)~AA")
 vec <- array()
 for (i in 1:15){ vec[i] <- parse(text=labs[(i%%4)+1]) }
# plot(c(1:4), c(1:4), axes=F, xlab="Haplotype",ylab="", par(mar=c(8,3,2,1)))
# axis(1, at=c(1:4), labels=vec, las=2)
# mtext(parse(text=labs[1]), side=3)



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
         cat(paste0(hapString,"\n")) # remove this in the final version. Used for diagnostics only.
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



# Create the boxplots per haplotype per trait
 for (i in 1:length(TRAITS)) {
     cat(paste0("Processing ",TRAITS[i],"..."))
     pdf(file=paste0("images/",TRAITS[i],"_haplotype_boxplot.pdf"), width=7, height=5)
     # res<-boxplot(as.formula(paste0(TRAITS[i],"~hap")), hap_pheno,
     #         las=2, ylab=TRAITS[i], main=paste0("Haplotypes for ",TRAITS[i]),
     #         par(mar = c(8, 4, 4, 1) + 0.1))

     # Start of replacement code to label the x-axis with the nucleotide sequence
     res<-boxplot(as.formula(paste0(TRAITS[i],"~hap")), hap_pheno,
                  las=2, xaxt='n', #axes=F,
                  ylab=TRAITS[i], main=paste0("Haplotypes for ",TRAITS[i]),
                  par(mar = c(8, 4, 4, 1) + 0.1))
     axis(1, at=c(1:length(res$names)), labels=getHaplotypeLabels(res$names, bim), las=2, cex.axis=0.8)
     # End of replacement code for allelic dosage plot

     # Show the frequency of each haplotype
     mtext(paste("(n=",res$n,")",sep=""), at=c(1:length(res$names)), side=1, line=5, cex=0.7)

     # Show a note that these haplotypes are based on A1 allele dosage
     #mtext("Haplotypes based on minor allele dosage", side=1, line=6.5)

     res<-dev.off() # Just prevent the annoying messsage when a device is closed.
     cat("DONE.\n")
 }




# Serialize geno2 in preparation for the one-way ANOVA
cat("Serializing the genotype data...")
#geno3 <- geno2 %>% gather(snp,genotype,snp_06_1766647_C:snp_06_1769727_T)
geno3 <- geno2 %>% gather(snp,genotype,2:length(colnames(geno2)))
cat("DONE.\n")


# Bring the phenotype into the serialized data
cat("Appending the phenotype to serialized genotype data...")
geno3_pheno <- dplyr::inner_join(geno3, pheno, by=c("IID"="iid"))
#..then convert the snp and genotype columns as factors
geno3_pheno$snp <- as.factor(geno3_pheno$snp)
geno3_pheno$genotype <- as.factor(geno3_pheno$genotype)
cat("DONE.\n")




# Create the individual boxplot per genotype per SNP per trait.
SNPS <- unique(geno3_pheno$snp)
#TRAITS <- colnames(geno3_pheno[,-c(1:3)]) # starts at column 4 in geno3_pheno
results <- vector('list', length(TRAITS)) # preallocate memory
names(results) <- TRAITS # assign the trait vector as names of the list items
for (i in 1:length(TRAITS)){
    cat(paste0("Processing ",TRAITS[i],"..."))
    #df <- data.frame(snp=character(), R2=double(), SE=double(), p=double(), log10p=double()) # empty data frame
    for (j in 1:length(SNPS)) {
        pdf(file=paste0("images/",TRAITS[i],"_", paste(unlist(strsplit(as.character(SNPS[j]),"_"))[1:3],collapse="_"),
                        "_boxplot.pdf"), width=5, height=5)

        # PREVIOUS WORKING VERSION BUT SHOWING ONLY THE ALLELE DOSAGE.
        # res<-boxplot(as.formula(paste0(TRAITS[i]," ~ genotype")), geno3_pheno %>% filter(snp==SNPS[j]),
        #         xlab=paste0("Dosage of A1 allele (", unlist(strsplit(as.character(SNPS[j]),"_"))[4],")"),
        #         ylab=paste0(TRAITS[i]),
        #         main=paste0(paste(unlist(strsplit(as.character(SNPS[j]),"_"))[1:3],collapse="_"),
        #                     " on ",TRAITS[i]))
        # #freq <- as.matrix(table((geno3_pheno %>% filter(snp==SNPS[j]))$genotype))
        # mtext(paste("(n=",(table((geno3_pheno %>% filter(snp==SNPS[j]))$genotype)),")",sep=""),
        #       at=c(1,2,3),side=1, line=2, cex=0.7)

        # Start of replacement code that shows A1A1, A1A2, A2A2 coding
        res<-boxplot(as.formula(paste0(TRAITS[i]," ~ genotype")), geno3_pheno %>% filter(snp==SNPS[j]),
                     las=1, xaxt='n',
                     xlab=paste0("Dosage of A1 allele (", unlist(strsplit(as.character(SNPS[j]),"_"))[4],")"),
                     ylab=paste0(TRAITS[i]),
                     main=paste0(paste(unlist(strsplit(as.character(SNPS[j]),"_"))[1:3],collapse="_"),
                                 " on ",TRAITS[i]))
        axis(1, at=c(1:length(res$names)), labels=getGenotypeLabels(j, bim), las=1)
        # End of replacement code

        # Show the frequency of each haplotype
        mtext(paste("(n=",res$n,")",sep=""), at=c(1:length(res$names)), side=1, line=2, cex=0.7)

        res<-dev.off()
    }
    results[[i]] <- df # store the compiled association for SNPj into the ith element of results
    cat("DONE.\n")
}


