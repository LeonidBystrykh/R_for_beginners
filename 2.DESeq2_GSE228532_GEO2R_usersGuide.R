# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Differential expression analysis with DESeq2
#This is modified version initially based on GEO2R, then modified according to the DESeq2 users guide
library(DESeq2)
library(GEOquery) #to fetch sample ata
library("pheatmap")
#BiocManager::install("vsn")
library("vsn")

# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE228532", "file=GSE228532_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")
#OR use path to download file then
tbl<-read.csv("~/Downloads/GSE228532_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="\t", row.names = 1)

# load gene annotations from GEO online
apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
rownames(annot) <- annot$GeneID

#import sample info using getGEO
s<-getGEO("GSE228532")
samples<-s[["GSE228532_series_matrix.txt.gz"]]@phenoData@data
colnames(samples)
sample_info<-samples[,c("title", "cell line:ch1","cell type:ch1","treatment:ch1")]
class(sample_info)
#sample_info<-as.data.frame(samples[,"title"])

# just make sample table based on the info in title column
gs <- factor(c(rep(1,3),rep(2,3),rep(3,3)))
groups <- make.names(c("ctrl","CM_ctl","CM_cc"))
levels(gs) <- groups
sample_info <- data.frame(Group = gs, row.names = colnames(tbl))
sample_info$batch<-factor(c(rep(c(1,2,3),3)))
str(sample_info)

# pre-filter low count genes
# keep genes with at least N counts > 10, where N = size of smallest group
keep <- rowSums( tbl >= 10 ) >= min(table(gs))
tbl <- tbl[keep, ]

#proper place for boxplot
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
boxplot(log2(tbl+1), main="GSE228532", ylab="log2(counts)", las=2, 
        col=as.numeric(sample_info$Group), xlim=c(0,11))
legend("topright", groups, fill=palette(), bty="n")

#major DESeq2 lines
ds <- DESeqDataSetFromMatrix(countData=tbl, colData=sample_info, design= ~Group+batch+0)
#figures
#from DESeq2 users guide 
vsd <- vst(ds, blind=FALSE) #Variance stabilizing transformation
rld <- rlog(ds, blind=FALSE) #Regularized log transformation
plotPCA(vsd, intgroup="Group")

#check similarity between samples using pheatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         show_colnames = F,
         show_rownames = F,
         annotation_col = sample_info)

#another dispersion plot from DESeq2 users guide
meanSdPlot(assay(vsd)) #from vsn library

#run analysis
ds<-DESeq(ds, test="Wald")#, reduced= ~1)#, test="LRT")

#samples after normalization
boxplot(log2(counts(ds, normalize=T)), main="GSE228532", ylab="log2(counts)", las=2, 
        col=as.numeric(sample_info$Group))

ds$Group
counts(ds, normalized=T)

#dispersion plot from DESeq2
plotDispEsts(ds, main="GSE228532 Dispersion Estimates")

#results, select the contrast 
r <- results(ds, alpha=0.05, pAdjustMethod ="fdr")# default results
r
summary(r)

#extracting results more slowly
resultsNames(ds)
contrast=c("Group", "CM_cc","CM_ctl") #use this if you want groups comparison
#contrast=c("batch","2","3") #or this if you interested in reproducibility in batches
r <- results(ds, alpha=0.05, pAdjustMethod ="fdr", contrast=contrast)

tT <- r[order(r$padj),] 
tT <- merge(as.data.frame(tT), annot, by=0, sort=F)
tT <- subset(tT, select=c("GeneID","padj","pvalue","stat","baseMean","Symbol","Description"))

#missing line in GEO2R script
tops<-tT[tT$padj<0.05,]
#finally you can save the results table with DEGs only
#write.table(tops, file=stdout(), row.names=F, sep="\t")

#Figures, which you might find useful

# volcano plot
old.pal <- palette(c("#00BFFF", "#FF3030")) # low-hi colors
par(mar=c(4,4,2,1), cex.main=1.5)
plot(r$log2FoldChange, -log10(r$padj), main=paste(groups[2], "vs", groups[3]),
     xlab="log2FC", ylab="-log10(Padj)", pch=20, cex=0.5)
with(subset(r, padj<0.05 & abs(log2FoldChange) >= 0),
     points(log2FoldChange, -log10(padj), pch=20, col=(sign(log2FoldChange) + 3)/2, cex=1))
legend("bottomleft", title=paste("Padj<", 0.05, sep=""), legend=c("down", "up"), pch=20,col=1:2)

#MA or MD figure from DESeq2
plotMA(r, main=paste(groups[2], "vs", groups[3], cex=2))

#heatmap selected genes
s<-tT$GeneID[1:50]
df<-as.data.frame(counts(ds, normalized=T)[as.character(s),])
rownames(df)<-tT$Symbol[1:50]
pheatmap(df, scale="row", cluster_cols=F, annotation_col = sample_info,
         show_colnames = F, fontsize_row = 5)

#plot individual expression
plotCounts(dds=ds, gene="3429", intgroup="Group", normalized=T, main="IFI27",
           col=c(rep("orange",3),rep("red",3),rep("brown",3)), pch=19)

