#Limma-edgeR
#Use example lines from limma, chapter RNAseq

library(edgeR)
library(GEOquery)
#get the data
setwd("~/Documents/R_scripts/NCBI_GEO/data")
data_table<-read.csv("GSE133281_DataRaw_ReadCounts.txt.gz", sep="\t")
filter<-duplicated(data_table$X)
data<-data_table[!filter,]
rownames(data)<-data$X
data<-data[,-1]

#sample info
samples<-getGEO("GSE133281",GSEMatrix =TRUE, AnnotGPL = FALSE)
sam<-samples[["GSE133281_series_matrix.txt.gz"]]@phenoData@data
colnames(sam)
selected<-c("title","geo_accession","age:ch1","cell type:ch1","disease diagnosis:ch1"  )
samples<-sam[,selected]
rownames(samples)<-samples$title
rownames(samples)
filter<-grep("MDS", rownames(samples))
samples$`age:ch1`[filter]="unk"
colnames(samples)<-c("title","geo_acc","age","cell_type","disease")

#align info for samples and data
colnames(data) %in% samples$title
filter<-colnames(data) %in% samples$title
data<-data[,filter]
samples$title %in% colnames(data)
filter<-samples$title %in% colnames(data)
samples<-samples[filter,]
data<-data[,samples$title]


#the limma_edgeR lines and my supporting lines and illustrations
dge <- DGEList(counts=data, samples=samples)
hist(log2(dge$counts+1))
dge$samples
#check outliers 
readcounts<-colSums(dge$counts)
genecounts<-colSums(dge$counts>0)
plot(genecounts, readcounts)
text(genecounts, readcounts, colnames(dge$counts))
plot(readcounts, dge$samples$lib.size)
dge$samples

#back to major lines, model
design=model.matrix(~age+0, dge$samples)
#filter
keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes=FALSE]
#see the raw samples
boxplot(log2(dge$counts+1), las=2, col=as.numeric(as.factor(dge$samples$age))+1)
#normalize
dge <- calcNormFactors(dge)
#extra filtering
dge$samples
#log-transform
logCPM <- cpm(dge, log=TRUE, prior.count=3)
#min(logCPM)
#hist(logCPM)
#
logCPM[logCPM<0]=0
v<-voom(dge, plot=T)
#visualize again
readcounts<-colSums(logCPM)
genecounts<-colSums(logCPM>0)
plot(genecounts, readcounts, pch=19, col="red", cex=0.5)
text(genecounts, readcounts, colnames(dge$counts))
#and more visuals
boxplot(logCPM, las=2, col=as.numeric(as.factor(dge$samples$age))+1)
plotMDS(logCPM, col=as.numeric(as.factor(dge$samples$age))+1)

#fit and final table
fit <- lmFit(v, design)
fit$design
cont.matrix <- makeContrasts(ageAged - ageYoung, 
                             levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
plotSA(fit2, main="All samples") #that's mine again
fit2 <- eBayes(fit2, 0.01)
result<-decideTests(fit2, p.value=0.05)
summary(result)
names<-rownames(result[abs(result[,1])==1,])

tT<-topTable(fit2, number=Inf)
tops<-tT[names,]
tops
#save the data
getwd()
write.csv(tops, "DE_limma-edgeR-voom.csv")

#compare
list.files()
regular<-read.csv("DE_limma-edgeR-more10.csv", row.names=1)
voom<-read.csv("DE_limma-edgeR-voom.csv", row.names=1)
#install.packages("ggvenn")
library(ggvenn)
x <- list(reg=row.names(regular),
          voo=row.names(voom))
x
length(x)
ggvenn(x, show_elements = F)
#4 is the maximum


