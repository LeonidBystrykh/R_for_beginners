#data as in alternative GLM-edgeR. script #5
#another aging in Hs HSC
#groups old vs young, male vs female, absolute ages known
#example using DEGreport figures in DESeq2

library(DESeq2)
library(GEOquery)
library(DEGreport)
library(Glimma)
library(org.Hs.eg.db)
columns(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)

#get data: go to NCBI GEO, download and untar GSE104406_RAW.tar
setwd("~/Documents/NCBI_GEO/Hs_aging/GSE104406_RAW")
files<-list.files(pattern=".gz")
files
data<-edgeR::readDGE(files, columns=c(1,2)) #column 4 is redundant
head(data$counts,5)
table(is.na(data$counts))
data$counts<-na.omit(data$counts)
dim(data$counts)

#simplify ensembl names, #add gene annotation at this stage
genes<-gsub("\\..*","",rownames(data))
table(duplicated(genes))
rownames(data)<-genes
annots <- select(org.Hs.eg.db, keys=genes,
                 columns=c("SYMBOL","GENENAME"), keytype="ENSEMBL")

#sample info table
samples<-getGEO("GSE104406",GSEMatrix =TRUE, AnnotGPL = FALSE)
sample_info<-samples[["GSE104406_series_matrix.txt.gz"]]@phenoData@data
colnames(sample_info)
reduced<-sample_info[c("title","geo_accession",
                       "donor age:ch1","donor sex:ch1")]
colnames(data)
rownames(reduced)<-colnames(data)
reduced$group<-"young"
old<-grep("Old", reduced$title)
reduced$group[old]<-"old"
colnames(reduced)<-c("title","geo", "age", "sex","group")
reduced$age[5]<-31

#check outliers by read counts and gene counts
header="raw GSE104406"
palette(c("#00BFFF", "#FF3035"))
plot(colSums(data$counts), colSums(data$counts>0), 
     main=header, xlab="reads", ylab="genes",
     col=as.numeric(as.factor(reduced$group)), pch=19)
text(colSums(data$counts), colSums(data$counts>0), reduced$title, cex=0.5)

#from DEGreport. It shows ratio of expression to the arerage gene
#as kind of a test that most of the genes do not change
degCheckFactors(log2(data$counts+1)) #from degreport
#major lines
#ds <- DESeqDataSetFromMatrix(countData=data, colData=reduced, design= ~group+sex+0)
#OR
dds <- DESeqDataSetFromMatrix(countData=data, colData=reduced, design= ~group+sex)
dim(counts(dds))
#filter as recommended in DESeq2 guide
smallestGroupSize <- 8
keep <- rowSums(counts(dds) >= 100) >= smallestGroupSize
dds <- dds[keep,]

#generic boxplot
par(mar=c(8,2,2,2))
color<-as.numeric(as.factor(colData(dds)$group))
boxplot(log2(counts(dds)+1), col=color, las=2,
        main=header, names=reduced$title, cex=0.25)

#run DE analysis
dds<-DESeq(dds, test="Wald")
# check FC change hypothesis
degCheckFactors(log2(counts(dds)+1))+ #from degreport
  ggtitle("Test for gene expression rations to the average gene, raw")
degCheckFactors(log2(counts(dds, normalized=T)+1))+ #from degreport
  ggtitle("Test for gene expression rations to the average gene, normed")
#note colData changed, mystery column  "replaceable"
cd<-as.data.frame(colData(dds))

#tried to remove replaceable but dropped
#filter<-cd$replaceable==T
#dds<-ds[,filter]
#cannot use dds, have to rerun DESeqDataSet

#from DESeq2 help
se <- SummarizedExperiment(log2(counts(dds, normalized=T) + 1),
                           colData=colData(dds))
# the call to DESeqTransform() is needed to
# trigger our plotPCA method. or use vsd transform 
plotPCA( DESeqTransform( se ), intgroup="group", ntop=1000) +
  ggtitle("Testing GSE104406 data")

glimmaMDS(dds, labels= dds$title) #to the reader

#select contrast and get the result
resultsNames(dds)
r <- results(dds, alpha=0.05, pAdjustMethod ="fdr", name="group_young_vs_old")
plotMA(r)
#glimmaMA(dds) it does default groups which is not what I want
#from DEGreport
degs <- degComps(dds, combs = "group", contrast = list("group_young_vs_old",c("group", "young", "old")))
names(degs)
head(deg(degs[[1]]),5)
degMA(degs[[1]], diff = 2, limit = 3, raw = T) +
  ggtitle(" MA from DEGreport")
#volcano in style of DEGreport
r<-r[order(r$padj),] 
r[["id"]] <- row.names(r)
show <- as.data.frame(r[1:100, c("log2FoldChange", "padj", "id")])
degVolcano(r[,c("log2FoldChange", "padj")], plot_text = show)+
  ggtitle("Volcano from DEGreport")

#finalize the data
merged=merge(as.data.frame(r), annots, by.x=0, by.y="ENSEMBL")
tT<-na.omit(merged)
tT<-tT[order(tT$padj),] 
colnames(tT)
show <- as.data.frame(tT[1:100, c("log2FoldChange", "padj", "SYMBOL")])
degVolcano(r[,c("log2FoldChange", "padj")], plot_text = show)

#make tops
tops<-tT[tT$padj<0.05,]
hist(log2(tops$baseMean))
#individual genes
r
degPlot(dds = dds, res = r, n = 6, xs = "group")+
  ggtitle("top 6 from results table")
#the same in one plot
degPlotWide(dds, rownames(dds)[1:5], group="group")

#testing genes by gender
data(geneInfo)
degSignature(dds,geneInfo, group = "group")+
  ggtitle("Groups by age")
degSignature(dds,geneInfo, group = "sex")+
  ggtitle("Groups by gender")

#try to find trends
ma = assay(rlog(dds))[row.names(r)[1:100],]
design(dds)
ma[1:5,1:5]
meta<-as.data.frame(colData(dds))
r <- degPatterns(ma,  meta, time = "group", plot=T)
res <- degPatterns(ma, colData(dds), time = "sex", col="group", plot = T) #this line failed

#write for later analysis 
setwd("~/Documents/R_scripts/RNAseq/DESeq2/DEGreport")
write.table(ma, "DEG_expressions_GSE104406.tsv", sep="\t")
write.table(as.data.frame(colData(dds)),"col_data_GSE104406.tsv", sep="\t")

#from DEGreport
#testing for randomness of DEGs across means and vars
counts <- counts(dds, normalized = TRUE)
degQC(counts, meta[["group"]], pvalue = r[["pvalue"]])

#testing covariates across PCA scales
resCov <- degCovariates(log2(counts(dds)+0.5),
                        colData(dds))
#covariates in a sample table
cor <- degCorCov(colData(dds))
