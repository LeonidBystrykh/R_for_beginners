#Trying time in DEReq2
#example from https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html

library(DESeq2)
library(ggplot2)
#BiocManager::install("fission")
library("fission")

data("fission")
dim(fission)
head(fission, 5)
class(fission)
fission$minute
#define Time as continuous variable
library(splines)
X <- ns(as.numeric(fission$minute), df=3)
#OR
X <- poly(as.numeric(fission$minute), degree=3)
#add to the data samples description
fission$X<-X
fission$X

###Here you choose either time as categorical or time as continuous var
#use time as groups. I added replicate instead of interactions
ddsTC <- DESeqDataSet(fission, ~ strain + minute + replicate)

#use time as continuous
ddsTC <- DESeqDataSet(fission, ~ X[,1] +X[,2]+ X[,3]+strain+replicate)

#Check content, remove empty lines
dim(counts(ddsTC))
colData(ddsTC)
#ddsTC$sizeFactor #NULL
#head(ddsTC@assays@data$counts,5)
filter<-rowSums(counts(ddsTC))>0
ddsTC<-ddsTC[filter,]
palette(c('#e41a1c','#377eb8','#4daf4a','#ff7f00','#984ea3','#76ef33','#a65628','#078ea8','#99aa77'))
boxplot(log2(counts(ddsTC, normalize=F)+1), las=2, 
        col=as.numeric(as.factor(colData(ddsTC)$strain )))

samples<-colData(ddsTC)
table(samples$strain )
table(samples$minute )

#this block is about use of transformed time as continuous variable
ddsTC<-DESeq(ddsTC, test="Wald")
resultsNames(ddsTC)
res<-results(ddsTC, list(c("X...1.","X...2.","X...3.")))
res$symbol <- mcols(ddsTC)$symbol
hist(res$padj[res$padj>0] )
#go to tops section

#for categorical factors within defined model
ddsTC <- DESeq(ddsTC, test="Wald")
resultsNames(ddsTC)
#either use contrast
res<-results(ddsTC, contrast=c("strain","mut","wt"))
head(res,5)
#or for multiple comparisons use list
res<-results(ddsTC, list(c("minute_15_vs_0","minute_30_vs_0",    
                           "minute_60_vs_0","minute_120_vs_0","minute_180_vs_0")))


#tops section, can be used for all res data
res$symbol <- mcols(ddsTC)$symbol
res<-na.omit(res)
tops<-res[res$padj<0.05, ]
tops<-tops[order(tops$padj),]
head(tops, 5)
dim(tops) 

#save the tops, give appropriate name
getwd()
setwd("~/Documents/R_scripts/RNAseq/DESeq2/time")
write.table(tops, "tops_Xspline_time.tsv", sep="\t")

#plot data for one gene
gene="SPAP8A3.04c"
main="hsp9"
plotCounts(ddsTC, 
           gene=gene,main=main,
           col=as.numeric(as.factor(samples$strain)), pch=19,
           intgroup = c("minute","strain"), returnData = F)

#or in ggplot stype
fiss <- plotCounts(ddsTC, #which.min(resTC$padj), #the last in the tops
                   gene=gene,main=main,
                   intgroup = c("minute","strain"), returnData = TRUE)
fiss$minute <- as.numeric(as.character(fiss$minute))
fiss
ggplot(fiss,
       aes(x = minute, y = count, color = strain, group = strain)) + 
  geom_point() + 
  stat_summary(fun=mean, geom="line") +
  scale_y_log10()+
  ggtitle(main)

#compare tops for continuous and categorical times
setwd("~/Documents/R_scripts/RNAseq/DESeq2/time")
list.files(pattern=".tsv")
cont<-read.table("tops_Xspline_time.tsv", sep="\t")
cate<-read.table("tops_cat_time.tsv", sep="\t")
library(eulerr)
v=venn(list(Conti=rownames(cont),
             Categ=rownames(cate)))
plot(v, counts=TRUE, main="Compare time dependent tops")
#Calculate significance of the overlap
Conti=rownames(cont)
Categ=rownames(cate)
all<-rownames(ddsTC)
ov<-Conti[Conti %in% Categ]
ov
#used stackoverflow solution
#https://stackoverflow.com/questions/18340123/calculate-venn-diagram-hypergeometric-p-value-using-r
#Thanks to Ferdinand.kraft

install.packages("gmp")
require(gmp)

enrich_pvalue <- function(N, A, B, k)
{
  m <- A + k
  n <- B + k
  i <- k:min(m,n)
  
  as.numeric( sum(chooseZ(m,i)*chooseZ(N-m,n-i))/chooseZ(N,n) )
}
enrich_pvalue(length(all), length(Conti)-length(ov), 
              length(Categ)-length(ov), length(ov))
[1] 9.926154e-168
#there is also phyper solution in Stackoverflow (link is there)
#there is also RVenn library with significance, but it did not work for me
#https://cran.r-project.org/web/packages/RVenn/vignettes/vignette.html
#there is also giant lib to test significance of gene lists but I had no time for it
#https://cran.case.edu/web/packages/GiANT/vignettes/giant_package_vignette.pdf

#for clustering and classification you need counts
library(Rfast, quietly =T) #for correlations and cora
library(scales) #for alpha
colorz<-c('#e41a1c','#377eb8','#4daf4a','#ff7f00','#984ea3','#76ef33','#a65628','#078ea8','#99aa77')

ov_counts<-subset(log2(counts(ddsTC, normalized=T)+1), rownames(counts(ddsTC)) %in% ov)
#then do MDS by correlation
hist(rowMax(ov_counts)-rowMin(ov_counts))
filter<-(rowMax(ov_counts)-rowMin(ov_counts))>3
table(filter)
input<-ov_counts[filter,]

mds.cor<-function(x){1-cora(t(x))} #for MDS
coords<-as.data.frame(cmdscale(mds.cor(input)), x.ret=F) #all 6000 are too many!
colnames(coords) <- c("Dim.1", "Dim.2")
plot(coords$Dim.1,coords$Dim.2, xlab="Dim.1", ylab="Dim.2",
     pch=19)
#try clustering
# Dissimilarity matrix
d <- dist(coords[,1:2], method = "euclidean")
# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete" ) #you can also make a tree if not too long list
clust<-as.data.frame(cutree(hc1, k=9))
colnames(clust)="clusters"
plot(coords$Dim.1,coords$Dim.2, xlab="Dim.1", ylab="Dim.2",
     pch=19, col=clust$clusters)
text(0,0,"Clusters colors", col="grey")
for (i in 1:9){
  text(0,-i/10,i,col=i)
}
#all CLUSTERS together
#par(bg = "white")
par(mfrow=c(3,3))
par(mar=c(2,2,2,2))
high<-round(max(as.matrix(input)),0)+1
low<-round(min(as.matrix(input)),0)-1
minutes<-as.numeric(samples$minute)
s<-grep("wt", samples$strain)
s<-grep("mut", samples$strain)
for (i in 1:9){
  clus1=subset(input, rownames(input) %in% rownames(clust)[clust$clusters==i])
  plot(minutes[s], colMeans(clus1)[s], col="grey", ylim=c(low, high),
       main=paste("Cluster", i))
  for (j in 1:dim(clus1)[1]){ 
    lines(minutes[s], clus1[j,s],   col=alpha(colorz[i], 
                    0.3))}
}
tops
clust
final_results_table<-merge(as.data.frame(tops), clust, by=0)
final_results_table[final_results_table$symbol=="hsp9",]
