# This is the code to show the distribution of the six data points forwarded to me by
# Vito regarding the 320PHY_2013WS_SEC data.
#
# Roslen Anacleto
# 20160304
#


library(dplyr)


# Set the working directory
setwd("/Users/ranacleto/Google Drive/02 - Stat analyses/31 - Genomics/50 - GWAS/15-HDRAv3/run/320PHY_2013WS_SEC/phenotype")

# Load the data
data <- read.table("320PHY_2013WS_SEC.txt", header=T, stringsAsFactors = F, sep="\t")

# Use this as basis to determine who in the "data" are to be excluded in the visualization.
iid <- read.table("am1/am1.txt", header=F, stringsAsFactors = F, sep="\t")

# Get the intersection of the original data source and the filtered data source
res <- dplyr::inner_join(iid, data, by=c("V1"="fid","V2"="iid"))

# visualize the histograms
traits <- c("am1","am2","mc_ap","sc_ap","app_am","app_ap")
for (i in 1:length(traits)) {
  pdf(file=paste0("figures/",traits[i],"_hist.pdf"),width=5, height=5)
  hist(data[,i+2], main=traits[i])
  dev.off()
}


