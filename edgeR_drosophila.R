#example edgeR for time dependent DGE analysis
#this is example from the major edgeR tutorial, chapter 4.9, page 114 (version 2024)
#BiocManager::install("org.Dm.eg.db")
library(org.Dm.eg.db)
library(AnnotationDbi)
library(edgeR)
CountFile <- "http://bowtie-bio.sourceforge.net/recount/countTables/modencodefly_count_table.txt"
Counts <- read.delim(CountFile, row.names=1)

SampleFile <- "http://bowtie-bio.sourceforge.net/recount/phenotypeTables/modencodefly_phenodata.txt"
Samples <- read.delim(SampleFile, row.names=1, sep=" ", stringsAsFactors=FALSE)
#check data
PooledCounts <- sumTechReps(Counts, ID=Samples$stage)
dim(PooledCounts)
colnames(PooledCounts)
#generate time scale
Hours <- seq(from=2, to=24, by=2)
Time <- paste0(Hours,"hrs")
Time
colnames(PooledCounts)[1:12]
#create DGEList variable. Only Embryos taken
y <- DGEList(counts=PooledCounts[,1:12], group=Time)
y$samples
#filter by expression
keep <- filterByExpr(y)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
#add gene annotation
Anno <- select(org.Dm.eg.db, keys=rownames(y), keytype="FLYBASE",
               column="SYMBOL")
rownames(y$counts)
genes<-subset(Anno, Anno$FLYBASE %in% rownames(y$counts))
table(genes$FLYBASE==rownames(y$counts))
head(genes, 5)
head(y$counts, 5)
y$genes<-genes$SYMBOL
#normalise
y <- calcNormFactors(y)
y$samples
#look at data
plotMDS(y, labels=Hours)
#design computes orthogonal polynomials
#polynomial spline, three coefficients mean linear, quadratic and cubic coefficients
X <- poly(Hours, degree=3)
design <- model.matrix(~ X)
#or using spine function, it is cubic spline fit if df==3
library(splines)
X <- ns(Hours, df=3) #ns generates matrix of cubic spines
X
design <- model.matrix(~ X)
design
#estimate common dispersion
y <- estimateDisp(y, design)
plotBCV(y)
#fit
fit <- glmQLFit(y, design, robust=TRUE)

#In a time course experiment, we are looking for genes that change expression level over time.
#Here, the design matrix uses 3 natural spline basis vectors to model smooth changes over
#time, without assuming any particular pattern to the trend. We test for a trend by conducting
#F-tests on 3 df for each gene:

fit$coefficients
fit <- glmQLFTest(fit, coef=2:4)
tt<-as.data.frame(topTags(fit, n=Inf))
tops <- subset(as.data.frame(topTags(fit, n=Inf)), FDR<0.001)
#polytops<-tops
cubitops<-tops
library(eulerr)
v=euler(list(Poly=rownames(polytops),
               Cubi=rownames(cubitops)))
plot(v, counts=TRUE)
head(tops)
#Finally, we visualize the fitted spline curves for the top four genes. We start by computing
#the observed and fitted log-CPM values for each gene:
logCPM.obs <- cpm(y, log=TRUE, prior.count=fit$prior.count)
logCPM.fit <- cpm(fit$fitted.values, log=TRUE)
head(logCPM.fit)
#We then loop through the first four genes in the topTags table, plotting the observed and
#fitted values for each gene:
par(mfrow=c(3,3))
for(i in 4:12) {FlybaseID <- row.names(tops)[i] 
  Symbol <- tops$Symbol[i] 
  logCPM.obs.i <- logCPM.obs[FlybaseID,] 
  logCPM.fit.i <- logCPM.fit[FlybaseID,] 
  plot(Hours, logCPM.obs.i, ylab="log-CPM", pch=16, main=tops$ID[i]) 
  lines(Hours, logCPM.fit.i, col="red", lwd=2)}

#try to cluster time profiles
library(TMixClust)
input<-subset(logCPM.fit, rownames(logCPM.fit) %in% rownames(tops))[1:500,]
dim(input)
cluster_obj = TMixClust(input, nb_clusters = 9)
#cluster_obj[1:5]
plot_silhouette(cluster_obj)
#visualize
palette(c('#e41a1c','#377eb8','#4daf4a','#ff7f00','#984ea3','#76ef33','#a65628','#078ea8','#99aa77'))
colorz<-c('#e41a1c','#377eb8','#4daf4a','#ff7f00','#984ea3','#76ef33','#a65628','#078ea8','#99aa77')
m=9
clus1<-input[cluster_obj$em_cluster_assignment==m,]
#plot individual cluster
plot_time_series_df(clus1, plot_title =paste("Cluster",m), 
                    time_points = cluster_obj$ts_time_points, data_color=m)
#try MDS approach
#define similarity functions 
install.packages("Rfast")
library(Rfast, quietly =T) #for correlations and cora
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
m=6
#rownames(clust1)[clust1$clusters==m]
clus1=subset(input, rownames(input) %in% rownames(clust)[clust$clusters==m])
plot_time_series_df(clus1, plot_title =paste("Cluster",m), 
                    time_points = cluster_obj$ts_time_points, data_color=m)
#or all together
par(bg = "white")
par(mfrow=c(3,3))
high<-round(max(as.matrix(input)),0)+1
low<-round(min(as.matrix(input)),0)-1
library(scales)
for (i in 1:9){
  clus1=subset(input, rownames(input) %in% rownames(clust)[clust$clusters==i])
  plot(Hours, colMeans(clus1), col="grey", ylim=c(low, high),
       main=paste("Cluster", i))
  for (j in 1:dim(clus1)[1]){ 
    lines(Hours, clus1[j,], col=alpha(colorz[i], 0.3))}
}
tops
clust
final_results_table<-merge(tops, clust, by=0)
#That's more or less the end