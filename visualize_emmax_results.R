# This script is used to visualize the output of emmax using qqman.
#
# Roslen Anacleto
# 20160304
#


library(qqman)
library(svglite)

setwd("/Users/ranacleto/Google Drive/02 - Stat analyses/31 - Genomics/50 - GWAS/15-HDRAv3/run/320PHY_2013WS_SEC/phenotype")


# Count how many SNPs on each chromosome
#as.data.frame(table(data$CHR))


# For reference in creating Manhattan plots, see
# http://www.gettinggeneticsdone.com/2014/05/qqman-r-package-for-qq-and-manhattan-plots-for-gwas-results.html

snps <- 122785
genomewideline <- -log10(0.05/snps)
suggestiveline <- F
# Note: The genomewide line indicates the minimum -log10(p) that is considered as a signal of significant
# association. To recompute it back to the minimum association p, p=10^(-genomewideline).

ylim <-10
xlim <-ylim

traits <- c("am1","am2","app_am","app_ap","mc_ap","sc_ap")
trait_labels <- c("Amylose1","Amylose2","Apparent Amylose","Apparent Amylopectin",
                  "Medium-chain Amylopectin","Short-chain Amylopectin")

for (i in 1:6) {
   cat(paste0("reading ",traits[i],"/",traits[i],"_emmax.ps.qqman..."))
   data <- read.table(paste0(traits[i],"/",traits[i],"_emmax.ps.qqman"), header=T, stringsAsFactors = F)
   cat("DONE.\n")

   cat("generating manhattan plot...")
   #svglite(file="am1_emmax.svg", width=5, height=5)
   png(filename=paste0(traits[i],"/",traits[i],"_emmax_manhattan.png"))
   manhattan(data, ylim=c(0,ylim), main=trait_labels[i],
             suggestiveline = suggestiveline, genomewideline = genomewideline)
   dev.off()
   cat("DONE.\n")

   cat("generating QQ plot...")
   #svglite(file="am1_emmax_qq.svg", width=5, height=5)
   png(filename=paste0(traits[i],"/",traits[i],"_emmax_qq.png"))
   qq(data$P, main = paste0("Q-Q plot (",trait_labels[i],")"), xlim=c(0,xlim), ylim=c(0,ylim))
   abline(0,1)
   dev.off()
   cat("DONE.\n\n")
}

