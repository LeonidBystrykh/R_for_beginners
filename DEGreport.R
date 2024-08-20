#testing DEGreport
#BiocManager::install("DEGreport")
library(DEGreport)
data(humanGender)
library(DESeq2)
library(Glimma)
idx <- c(1:10, 75:85)
dds <- DESeqDataSetFromMatrix(assays(humanGender)[[1]][1:1000, idx],
                              colData(humanGender)[idx,], design=~group)
colData(dds)
dds <- DESeq(dds)

res <- results(dds)
res
counts <- counts(dds, normalized = TRUE)
dim(counts)
design <- as.data.frame(colData(dds))
#compare changes with respect to average gene
degCheckFactors(counts) #no option for colors?
#check whether DEGs depend on the average expression or variance (they should not) 
degQC(counts, design[["group"]], pvalue = res[["pvalue"]])
#check for covariates in the data
resCov <- degCovariates(log2(counts(dds)+0.5),
                        colData(dds))
#Also, the correlation among covariates and metrics from the analysis can be tested
cor <- degCorCov(colData(dds))
#it is pretty evident correlation between lib.size and sizeFactor, proof of principle

#analysis at the dds level
degs <- degComps(dds, combs = "group",
                 contrast = list("group_Male_vs_Female",
                                 c("group", "Female", "Male")))
names(degs)
#you can call any contrast available
deg(degs[[1]])
#MA plot
glimmaMA(dds)
degMA(degs[[1]], diff = 2, limit = 3, raw = TRUE)
#note, original plotMA works both with dds and with res data
res<-results(dds)
design(dds)
res
plotMA(res)

#trying volcano plot
res[["id"]] <- row.names(res)
show <- as.data.frame(res[1:10, c("log2FoldChange", "padj", "id")])
degVolcano(res[,c("log2FoldChange", "padj")], plot_text = show)

#plot expression of individual genes
degPlot(dds = dds, res = res, n = 6, xs = "group")
#other version of the same
degPlotWide(dds, rownames(dds)[1:5], group="group")

#Markers can be used to show whether different conditions are enriched in different markers. 
#For instance, in this example, Females and Males show different total expression for chromosome X/Y markers

data(geneInfo)
degSignature(humanGender, geneInfo, group = "group")
ma = assay(rlog(dds))[row.names(res),]
design(dds)
r <- degPatterns(ma,  design, time = "group")
