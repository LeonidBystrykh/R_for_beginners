#testing absolute ages in combination with categorical variable (sex)
#aging in Hs HSC GSE104406, 20 samples in groups young and aged (old)
# by LV Bystrykh, Sept 2024. 
#There are some glitches left. Final illustrations need further details
library(edgeR)
library(org.Hs.eg.db)
library(AnnotationDbi)

#this is a short cut to the data load. 
#Full detais are in the script 5.alt_glm_edgeR.R
setwd("~/Documents/NCBI_GEO/Hs_aging/GSE104406_RAW")
#get counts table
counts<-read.table("GSE104406_counts.tsv")
#sample info
samples<-read.table("samples.txt")
hist(samples$age)

#major lines
y <- DGEList(counts=counts, samples=samples)
y$samples
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
#palette from colorbrewer
palette(c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628'))
#brief check with boxplot
boxplot(log2(y$counts+1), col=as.numeric(as.factor(y$samples$group)), 
        las=2, names=1:20)
#normalization 
y <- normLibSizes(y)
boxplot(log2(cpm(y)), col=as.numeric(as.factor(y$samples$group)), 
        las=2, names=1:20)
#MDS plot
plotMDS(cpm(y, log=T), col=as.numeric(as.factor(y$samples$group)), labels=1:20)
#simplify rownames and correct colnames in the sample table
e_keys<-gsub('\\..*', '',rownames(cpm(y)))
logCPM.obs <- cpm(y, log=TRUE)#, prior.count=fit$prior.count)
rownames(logCPM.obs)<-e_keys
rownames(y$counts)<-e_keys

#make design
#option 1: use age groups
design <- model.matrix(~0+group+sex, y$samples)
#OR
X <- poly(as.numeric(y$samples$age), degree=3)
#option 2: use real ages
design <- model.matrix(~ X + sex, y$samples)

design
#dispersion
y <- estimateDisp(y, design) 
plotBCV(y)

#fit
fit <- glmQLFit(y, design, robust=TRUE) #note, it was glmFit()
fit$design

#do this for groups (if old and young is used)
cont.matrix <- makeContrasts(groupold - groupyoung, 
                            levels=design)
head(fit$coefficients,5)
fit <- glmQLFTest(fit, contrast=cont.matrix)
plotMD(fit) #alike to plotMA

#do this for ages (no age groups in design!)
fit <- glmQLFTest(fit, coef=2:4)

#converted both designs into this
tt<-as.data.frame(topTags(fit, n=Inf))
tops <- subset(as.data.frame(topTags(fit, n=Inf)), FDR<0.05)

#annotate tops and save
annots <- select(org.Hs.eg.db, keys=rownames(tops),
                 columns=c("SYMBOL","GENENAME"), keytype="ENSEMBL")
merged<-merge(annots, tops, by.x="ENSEMBL", by.y=0)
merged<-na.omit(merged)
#Choose one option out of 2
grouptops<-merged #if age groups were used
agetops<-merged #if absolute ages were used

#compare tops
library(eulerr)
v=euler(list(Ages=rownames(agetops),
             Groups=rownames(grouptops)))
plot(v, quantities = TRUE, counts=TRUE, main="Consistency of tops")
#write.table(merged, "DE_FDR005.tsv", sep="\t")

#From here we work on the illustrations of the found agetops genes
#get observed and fit values. Note: logCPM.obs was already done above
logCPM.fit <- cpm(fit$fitted.values, log=TRUE)
cpm_obs<-subset(logCPM.obs, rownames(logCPM.obs) %in% merged$ENSEMBL)
cpm_fit<-subset(logCPM.fit, rownames(logCPM.obs) %in% merged$ENSEMBL)
#This is my way of getting Ages and Sex into the expression tables
#It is not ideal, but works for me. You kan do it differently
Ages<-as.numeric(y$samples$age)
Sex<-as.numeric(as.factor(y$samples$sex)) #1 for f 2 for m
cpm_obs<-rbind(cpm_obs, Ages, Sex)
cpm_fit<-rbind(cpm_fit, Ages, Sex)
#order
cpm_obs<-cpm_obs[,order(cpm_obs["Ages",])]
cpm_fit<-cpm_fit[,order(cpm_fit["Ages",])]

#plot one gene in a time
i=2
plot(cpm_obs["Ages",],cpm_obs[i,], main=merged$SYMBOL[i], col=cpm_obs["Sex",], pch=19)
M<-cpm_fit[,cpm_fit["Sex",]==2]
lines(M["Ages",],M[i,], col=2, lwd=2)
F<-cpm_fit[,cpm_fit["Sex",]==1]
lines(F["Ages",],F[i,], col=1, lwd=2)

#more genes in a loop, all 9
par(mfrow=c(3,3))
for (i in 1:9){
  plot(cpm_obs["Ages",],cpm_obs[i,], main=merged$SYMBOL[i], col=cpm_obs["Sex",], pch=19)
  M<-cpm_fit[,cpm_fit["Sex",]==2]
  lines(M["Ages",],M[i,], col=2, lwd=2)
  F<-cpm_fit[,cpm_fit["Sex",]==1]
  lines(F["Ages",],F[i,], col=1, lwd=2)
}
dim(t(cpm_obs))

#try all or cluster in groups and show clusters profiles
library(TMixClust)
#truncate expression table to expression rows only
size<-dim(cpm_fit)
input<-as.data.frame(cpm_fit[1:(size[1]-2),])
#plot all genes together (it's a mess)
plot_time_series_df(input, time_points=as.numeric(cpm_fit["Ages",]))

#do clustering by expression profiles, number of clusters is arbitrary                  
cluster_obj = TMixClust(input,time_points = as.numeric(cpm_fit["Ages",]), nb_clusters = 4)
plot_silhouette(cluster_obj)
colr<-cluster_obj$em_cluster_assignment
colr
plotMDS(t(input), 
        cex=0.5,
       # label=1:size[1]-2, 
        col=colr)

#show one cluster in a time
m=3
clus1<-cpm_fit[cluster_obj$em_cluster_assignment==m,]
dim(clus1)
#plot individual cluster
plot_time_series_df(clus1, plot_title =paste("Cluster",m), 
                    time_points = cluster_obj$ts_time_points, data_color=m)

#do all clusters
sex<-cpm_fit["Sex",]
sex
size<-dim(input)
size
input$cluster<-colr #stitch precalculated cluster values to the input table
#this will work for 4 clusters, otherwise change accordingly
par(mfrow=c(2,2))
#this is just a regular plot() function in a loop
#you can make it with lines, but then you use a scheme as above to split data first on M and F groups
for (i in 1:4){
  sub<-input[input$cluster==i,]
  size=dim(sub)
  plot(as.numeric(cpm_fit["Ages",]), colMeans(sub[,1:(size[2]-1)]),
       main=paste("Cluster", i, "( n=", size[1], ")"),
       pch=15+as.numeric(sex), col=i+sex-1, ylab="cpm", xlab="age")
}
#One of the two final tables to save if it is important.
merged_cpm<-merge(annots, input, by.x="ENSEMBL", by.y=0)

#The combination of absolute age and gender works. 
#This is an example, not a final data processing
