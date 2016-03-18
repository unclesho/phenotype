# Create the manhattan and qq plot
#
# Roslen Anacleto
# 20160306
#

suppressMessages(library(qqman))


# Get commandline arguments
args<-base::commandArgs(trailingOnly=T)

# Note: **basename** and **dirname** are useful in extracting the traitname and fullpath from the
# commandline parameter. However, since the requirement is to specify only the trait name and
# since the filesystem is organized such that all trait names are individual subdirectory names
# on the same level as this script, then there is no need to extract the basename and dirname.

#cat(is.null(args))
#cat("\n")

if (length(args) != 1) {
   cat("Usage: create_gwas_plots.R <trait>\n")
   quit(save="no", status=1, runLast=T) # non-zero status means an error has occurred.
} else if (!dir.exists(args[1])) {
   cat(paste0("ERROR: '",args[1],"' does not exist.\n"))
   quit(save="no", status=2, runLast=T)
}


# Load the trait and path from the command line paramter
trait=args[1]
path=trait


# Load the data
cat(paste0("Loading ",trait," data..."))
data <- read.table(paste0(path,"/",trait,"_emmax.ps.qqman"), header=T, sep="\t", stringsAsFactors=F)
data <- data[,-c(4)]
cat("DONE.\n")


# Set parameters based on the data
snps <- nrow(data)
genomewideline <- -log10(0.05/snps)
suggestiveline <- F
ylim <-10
xlim <-ylim


cat(paste0("generating ",paste0(path,"/",trait,"_emmax_manhattan.png"),"..."))
#svglite(file="am1_emmax.svg", width=5, height=5)
png(filename=paste0(path,"/",trait,"_emmax_manhattan.png"))
manhattan(data, ylim=c(0,ylim), main=trait,
          suggestiveline = suggestiveline, genomewideline = genomewideline)
res <- dev.off()
cat("DONE.\n")

cat(paste0("generating ",paste0(path,"/",trait,"_emmax_qq.png"),"..."))
#svglite(file="am1_emmax_qq.svg", width=5, height=5)
png(filename=paste0(path,"/",trait,"_emmax_qq.png"))
qq(data$P, main = paste0("Q-Q plot (",trait,")"), xlim=c(0,xlim), ylim=c(0,ylim))
abline(0,1)
res <- dev.off()
cat("DONE.\n")
