#alternative GLM-edgeR
#for data another aging in Hs HSC
library(edgeR)
library(GEOquery)

#get data
setwd("~/Documents/NCBI_GEO/Hs_aging/GSE104406_RAW")
files<-list.files(pattern=".gz")
files
data<-readDGE(files, columns=c(1,2)) #column 4 is redundant
head(data$counts,5)
table(is.na(data$counts))
data$counts<-na.omit(data$counts)
#sample info
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
write.table(reduced, "samples.txt", sep="\t")
#major lines
y <- DGEList(counts=data, samples=reduced)
y$samples
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
boxplot(log2(y$counts+1), col=as.numeric(as.factor(y$samples$group))+1, las=2, names=1:20)
y <- normLibSizes(y)
genes<-colSums(y$counts>0)
reads<-colSums(y$counts)
plot(reads, genes, col=as.numeric(as.factor(y$samples$group))+1, pch=19)
plotMDS(cpm(y, log=T), col=as.numeric(as.factor(y$samples$group))+1, labels=1:20)
#make design
design <- model.matrix(~0+group+sex, y$samples)
#dispersion
y <- estimateDisp(y, design) 
plotBCV(y)

#OR
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
#OR
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
fit$design
cont.matrix <- makeContrasts(groupold - groupyoung, 
                             levels=design)
#colnames(fit$coefficients)
#result <- glmLRT(fit, coef=2)
result <- glmLRT(fit, contrast=cont.matrix)
plotMD(result) #alike to plotMA
topTags(result, n=5)
tT<-topTags(result, n=Inf)


tops<-tT[tT$table$FDR<0.05,]
dim(tops$table)

#annotate
library(org.Hs.eg.db)
columns(org.Hs.eg.db)
library(AnnotationDbi)
#cut rownames at dot
e_keys<-gsub('\\..*', '',rownames(tops$table))
head(e_keys,5)
head(tops, 5)
annots <- select(org.Hs.eg.db, keys=e_keys,
                 columns=c("SYMBOL","GENENAME"), keytype="ENSEMBL")
tops$table$en_gene<-e_keys
merged<-merge(annots, tops, by.x="ENSEMBL", by.y="en_gene")
merged<-na.omit(merged)
