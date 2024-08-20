#Testing Glimma with example data
#BiocManager::install("Glimma")
library(Glimma)
library(edgeR)

dge <- readRDS(system.file("RNAseq123/dge.rds", package = "Glimma"))
glimmaMDS(dge) #to the viewer
dge$samples
glMDSPlot(dge, 
          groups=as.numeric(as.factor(dge$samples$group))) #to html
dge$samples
design <- model.matrix(~0+group, dge$samples)
dge <- estimateDisp(dge, design=design)
gfit <- glmFit(dge, design)
gfit$design
cont.matrix <- makeContrasts(groupLP - groupML, 
                             levels=design)
glrt <- glmLRT(gfit, contrast = cont.matrix)
glimmaMA(glrt)
glimmaMD(glrt)
glMDPlot(glrt)
#glMDRmd(glrt) not

glimmaVolcano(glrt, dge = dge)
#no volcano for html

#DESeq2 version
library(DESeq2)
dds <- DESeqDataSetFromMatrix(
  countData = dge$counts,
  colData = dge$samples,
  rowData = dge$genes,
  design = ~group
)
glimmaMDS(dds)
glMDSPlot(dds)
dds <- DESeq(dds)#, quiet=TRUE)
glimmaMA(dds)
glMDPlot(dds) #to the html, not well
glimmaVolcano(dds)
results(dds)
