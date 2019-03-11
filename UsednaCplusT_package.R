#PACKAGE dnaCplusT: Clustering plus hierarchical Testing for DNA copy number data

#Code by Mark van de Wiel <mark.vdwiel@vumc.nl>
# If the input is a CGHregions object; suggested parameter settings for CGHregions: 0.005 (very little information loss); 
# Then convert CGHregions object by running:
library(CGHbase)
library(dnaCplusT)
#load("cghreg.Rdata")
#datainfo <- data.frame(chromosome=chromosomes(cghreg),bpstart=bpstart(cghreg),nclone= nclone(cghreg))
#datacgh <- regions(cghreg)

#If the input is an aCGH regions file of which the first three colums represent base-pair start, chromosome and 
#number of clones in the region:
#setwd("C:\\VUData\\Saskia\\LungData")
datacgh <- read.delim("0.01_RegTxtFile_Mark_predictie.txt",header=TRUE, sep="\t", na.strings = c("NA","#N/A"),comment.char="%")
datainfo <- data.frame(chromosome=datacgh[,2],bpstart = datacgh[,3], bpend = datacgh[,4],nclone=datacgh[,5]) 
datacgh <- datacgh[,-(1:5)]


#data(datacgh)
#data(datainfo)
pvals <- pvalschisq(datacgh,group=c(7,30),groupnames=c("MSI+", "CIN+"),niter=10000) #10.000 iterations recommended for final results
#save(pvals,file="pvals10000.Rdata")
clust <- acghClustering(datacgh,datainfo)
clusters <-clust[,1]
pvalsclust <- pvclust(clusters,pvals)
pvalsreg<-pvals[[1]]
testresults <- acghHierarchTest(clusters,pvalsreg,pvalsclust,alpha=0.1,datainfo)
plot_dnaCplusT(clust,testresults,datacgh,datainfo,group=c(7,30),alpha=0.1)
