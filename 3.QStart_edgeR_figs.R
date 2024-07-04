#edgeR-Qstart

library(edgeR)
library(stringr)
#get data
setwd("~/Documents/NCBI_GEO/GSE115348")
list.files()
x <- read.delim("GSE115348_countTable_HPC.txt.gz")
#filter duplicated
filter<-duplicated(x$gene)
data<-x[!filter,]
rownames(data)<-data$gene
data<-data[,-1]
#make sample table and groups within
colnames(data)
samples<-t(data.frame(str_split(colnames(data), "_",4)))
colnames(samples)<-c("source", "age", "gender","type")
rownames(samples)<-colnames(data)
samples[,2]<- gsub("age","",samples[,2])
samples<-data.frame(samples)
samples$group<-"old"
samples$group[samples$age<30]<-"young"
#makeDGEList
y <- DGEList(counts=data, samples=samples)
y$samples
#figures
hist(log2(y$counts+1))
par(mar=c(9,3,3,3))
boxplot(log2(y$counts+1), las=2, cex.axis=0.7, col=as.numeric(as.factor(y$samples$group))+1)
reads<-y$samples$lib.size
genes<-colSums(y$counts>0)
plot(reads, genes)
text(reads, genes, 1:11)
plotMDS(y$counts, labels=1:11, col=as.numeric(as.factor(y$samples$group))+1, cex=2)
#filter
dim(y$counts)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
min(rowSums(y$counts))
#normalize
y <- normLibSizes(y)
hist(log2(rowSums(cpm(y))))
#figures for normed data
lognormed<-cpm(y, log=T)
boxplot(lognormed, las=2, cex.axis=0.7, col=as.numeric(as.factor(y$samples$group))+1)
reads<-colSums(lognormed)
genes<-colSums(lognormed>0)
plot(reads, genes, pch=19, cex=0.5, col="red")
text(reads, genes, 1:11)
plotMDS(lognormed, labels=1:11, col=as.numeric(as.factor(y$samples$group))+1, cex=1)
#design groups and fit
design <- model.matrix(~0+group, y$samples)
design
fit <- glmQLFit(y,design)
#show dispersion
plotQLDisp(fit) #replacement for plotSA
#contrasts and significance
fit$design
cont.matrix <- makeContrasts(groupold - groupyoung, 
                             levels=design)
qlf <- glmQLFTest(fit,contrast= cont.matrix)
#make results table
tT<-topTags(qlf, n=Inf)
#tops<-tT[tT$table$PValue<0.001,]
#dim(tops)
#tops$table
#alternative results table
result<-decideTests(qlf)
summary(result)
plotMD(qlf, main="Old to Young") #for fold changes
filter<-rownames(result[abs(result)==1,])
tops<-tT[filter,]
tops$table


