#Classic edgeR with single factor design
library(edgeR)

#get data
setwd("~/Documents/NCBI_GEO/E-GEOD-70537")
files<-list.files(pattern="count.txt")
files
data<-readDGE(files, columns=c(1,3)) #column 4 is redundant
head(data$counts,5)
#define groups
cells<-c("Hep3B","Huh7","Hep3B","Huh7")
#rename columns
colnames(data$counts)<-paste(cells, 6:9, sep="_")

#the core lines
y <- DGEList(counts=data, group=cells)
y$samples
dim(y$counts)
keep <- filterByExpr(y,group=cells)
y <- y[keep, , keep.lib.sizes=FALSE]
design <- model.matrix(~ 0+group, y$samples)
design
y <- normLibSizes(y)

#make lognormed for visualization
lognormed<-cpm(y, log=T)
boxplot(lognormed, col=c("red","blue","red","blue"))
#check samples similarity
library(corrplot)
M<-cor(lognormed)
corrplot(M, method="pie", 
         col.lim = c(0.5,1), 
         order= "alphabet",
         col = rainbow(40), 
         bg="gray") 
plotMDS(lognormed, col=c("red","blue","red","blue"))

#estimate dispersion as recommended
y <- estimateCommonDisp(y, design= design) #design changes results
plotBCV(y)
#find DE genes
et <- exactTest(y, pair=1:2)
tT<-topTags(et, n=Inf)
#END of recommendations

#my part of extracting
tops<-tT[tT$table$FDR<0.01,]
dim(tops$table)
tops$table
#add gene names
#https://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdfinstructions in 
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
columns(org.Hs.eg.db)
library(AnnotationDbi)
annots <- select(org.Hs.eg.db, keys=rownames(tops),
                 columns=c("SYMBOL","GENENAME"), keytype="ENTREZID")
merged<-merge(annots, tops, by.y=0, by.x="ENTREZID")
merged<-na.omit(merged)
