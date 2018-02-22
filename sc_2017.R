---
title: "sc_2017"
output: html_document
---

```{r,include=FALSE, cache=FALSE}
#library(reshape2)
library(ggplot2)
library(scater)
library(stringr)
#library(scRNA.seq.funcs)
library(scran)
#library(DESeq2)
#library(edgeR)
#library(RUVSeq)
options(stringsAsFactors = FALSE)
#library(limma)
library(readxl)
```


```{r read.data}
#bulk.NT<-read.table("/Users/yingy_adm/Downloads/T.txt",sep="\t",header = T)
#bulk.GT<-read.table("/Users/yingy_adm/Downloads/GT.txt",sep="\t",header = T)
# bulk.T<-data.frame(Name=bulk.GT$Name,
#                    bulk.count=bulk.GT$NumReads+bulk.NT$NumReads,
#                    bulk.tpm=bulk.GT$TPM+bulk.NT$TPM)
# single.T<-data.frame(Name=rownames(genecount),single.count=rowSums(genecount),
#                      single.tpm=rowSums(genetpm))
# merged.T<-merge(bulk.T,single.T,by=c("Name","Name"),all=F)
# rownames(merged.T)<-merged.T[,1]
# merged.T<-merged.T[,-1]
# merged.T<-merged.T[which(rowSums(merged.T)>0),]
# 
# plot(merged.T$single.count, merged.T$bulk.count, log="xy",pch=16, cex=0.6, xlab="summrized counts from single cell",ylab="counts from bulk data")
# corr=cor(merged.T$single.tpm/288, merged.T$bulk.tpm)
# 
# plot(merged.T$single.tpm/288, merged.T$bulk.tpm, log="xy",pch=16, cex=0.6, xlab="tpm from single cell", ylab="tpm from bulk data")

celltype<-read.table("~/Documents/sc_2017/celltype.txt",sep="\t",header = T)
celltype<-celltype[1:288,]
celltype$clone_exp<-NA
celltype$clone_exp[which(celltype$clonetype=="B"|celltype$clonetype=="C")]<-"3"
celltype$clone_exp[which(celltype$clonetype%in%c("D","E","F","G","H","I"))]<-"2"
celltype$clone_exp[which(celltype$clonetype=="S")]<-"1"

genecount <- read.delim("~/Documents/sc_2017/agg_count.txt",header=T,sep="\t",
                        row.names=1,check.names=FALSE)
genecount<-genecount[,1:288]
genetpm <- read.delim("~/Documents/sc_2017/agg_tpm.txt",,header=T,sep="\t",
                        row.names=1,check.names=FALSE)
genetpm<-genetpm[,1:288]

gene_biotype <- read.delim("/Users/yingy_adm/Documents/sc_2017/agg_gene_biotype.txt",sep = " ")
rownames(gene_biotype)<-gene_biotype$gene_name
gene_biotype<-gene_biotype[match(rownames(genecount),gene_biotype$gene_name),]
genecount<-data.matrix(genecount)
genetpm<-data.matrix(genetpm)
#rownames(genecount)=rownames(genetpm)=rownames(gene_biotype)=make.unique(as.character(gene_biotype$gene_name), sep = "_")
mapping<-read.table(file = "/Users/yingy_adm/Documents/sc_2017/meta_mapping.txt",header = F, sep = "")
rownames(mapping)<-mapping[,1]
mapping<-mapping[,-1]
colnames(mapping)<-celltype$well
mapping<-as.data.frame(t(mapping))
mapping<-mapping[1:288,]
```

#plot for mapping info
```{r mapping plot, echo=FALSE}
#library(reshape2)
plot_mapinfo(mapping=mapping[1:96,])
plot_mapinfo(mapping=mapping[97:192,])
plot_mapinfo(mapping=mapping[193:288,])
#plot_mapinfo(mapping=mapping[289:384,])
```

```{r}
scdata <- SingleCellExperiment(
    assays = list(counts = genecount), colData = celltype,rowData = gene_biotype)
#logcounts(scdata) <- log2(
#    calculateCPM(scdata, use.size.factors = FALSE) + 1)
assay(scdata, "logcounts_raw") <- log2(counts(scdata) + 1)
tpm(scdata)<-genetpm
```
#remove genes with 0 counts
#keep_feature <- rowSums(counts(scdata) > 0) > 0
#scdata <- scdata[keep_feature, ]
```{r}
#get index and list of IG,TR,MT,ERCC gene
ig_gene<-grepl("^IG_", gene_biotype$biotype, perl=TRUE)#a vector of 652 gene ID
tr_gene<-grepl("^TR_", gene_biotype$biotype, perl=TRUE)#a vector of 314 gene ID
mt_gene<-grepl("^MT-", gene_biotype$gene_name, perl=TRUE) &
  (gene_biotype$biotype=="protein_coding")#a vector of 13 gene ID
isSpike(scdata, type="ERCC")<-grepl("^ERCC-", rownames(genecount))
#blank_well<-celltype$well[grepl("^0",celltype$cell.type)]
#bulk_well<-celltype$well[grepl("^50",celltype$cell.type)]
```

```{r}
scdata <- calculateQCMetrics(
  scdata
  ,feature_controls = list(ERCC = isSpike(scdata, type="ERCC"),MT = mt_gene
  , IG = ig_gene, TR = tr_gene))

endog_genes <- !rowData(scdata)$is_feature_control
plotPCA(
    scdata[endog_genes, ],
    exprs_values = "logcounts_raw",
    shape_by = "source",
    colour_by = "plate",
    size_by = "total_features"
)

plotPCA(
    scdata[endog_genes, ],
    exprs_values = "tpm",
    shape_by = "source",
    colour_by = "plate",
    size_by = "total_features"
)#2D9,2C5 out after QC
```
scater_gui(scdata)"ENSG00000010610"CD4"ENSG00000153563"CD8A"ENSG00000172116"CD8B"ENSG00000111537"IFNG"ENSG00000169194"IL13"ENSG00000109471"IL2
ENSG00000010610 ENSG00000153563 ENSG00000172116 ENSG00000111537 ENSG00000169194 ENSG00000109471
ENSG00000168685"IL7R"
#check & filter on cells
```{r}
p1<-plotPhenoData(
     scdata,
     aes_string(x = "pct_counts_IG",
                y = "pct_counts_TR",
                colour = "cell_type",
                shape = "well"))
ggplotly(p1)

p2<-plotPhenoData(
    QC,
     aes_string(x = "pct_counts_IG",
                y = "pct_counts_TR",
                colour = "cell_type",
                shape = "well"))
ggplotly(p2)#2H10,2H6

plotPhenoData(
    scdata,
    aes_string(x = "total_features",
               y = "pct_counts_MT",
               shape = "plate",
               colour = "cell_specificity")
)
plotPhenoData(
    scdata,
    aes_string(x = "total_features",
               y = "pct_counts_ERCC",
               shape = "plate",
               colour = "cell_specificity")
)
plotPhenoData(
    scdata,
    aes_string(x = "total_counts",
               y = "total_features",
               shape = "plate",
               colour = "cell_specificity")
)
plotPhenoData(
    scdata,
    aes_string(x = "total_counts",
               y = "total_counts_ERCC",
               shape = "plate",
               colour = "cell_specificity")
)
```

```{r}
par(mfrow=c(2,2), mar=c(5.1, 4.1, 0.1, 0.1))
hist(scdata$total_counts/1e6, xlab="Library sizes (millions)", main="", 
    breaks=20, col="grey80", ylab="Number of cells")
hist(scdata$total_features, xlab="Number of expressed genes", main="", 
    breaks=20, col="grey80", ylab="Number of cells")
hist(scdata$pct_counts_ERCC, xlab="ERCC proportion (%)", 
    ylab="Number of cells", breaks=20, main="", col="grey80")
hist(scdata$pct_counts_MT, xlab="Mitochondrial proportion (%)", 
    ylab="Number of cells", breaks=20, main="", col="grey80")

ave.counts <- calcAverage(scdata)
hist(log10(ave.counts), breaks=100, main="", col="grey80", 
    xlab=expression(Log[10]~"average count"))

ave.counts <- calcAverage(QC)
hist(log10(ave.counts), breaks=100, main="", col="grey80", 
    xlab=expression(Log[10]~"average count"))

num.cells <- nexprs(QC, byrow=TRUE)
smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells", 
    xlab=expression(Log[10]~"average count"))
```
 # scater::plotPhenoData(
 #     scdata,
 #     aes_string(x = "total_features",
 #                 y = "pct_counts_MT",
 #                colour = "plate"))
```{r}
drop.pct.mapped<-(mapping$percent<30)
drop.n.reads<-(scdata$total_counts<200000)
drop.n.genes<-(scdata$total_features<2000)
drop.pct.mt<-(scdata$pct_counts_MT>30)
drop.pct.ercc<-(scdata$pct_counts_ERCC>50)
drop.control<-celltype$cell_type=="0"|celltype$cell_type=="50"

scdata$drop <- (drop.pct.mapped | drop.n.reads | drop.n.genes | drop.pct.mt | drop.pct.ercc|drop.control)
#scdata$clone_drop <-celltype[,10]=="N"|celltype[,10]=="S"|scdata$drop=="TRUE"

#GENE filtering
filter_genes <- apply(counts(scdata[ , !colData(scdata)$drop]), 1, 
                      function(x) length(x[x > 1]) >= 4)
rowData(scdata)$use <- filter_genes

QC <- scdata[rowData(scdata)$use, !colData(scdata)$drop]
data.frame(ByLibSize=sum(drop.n.reads),ByMapping=sum(drop.pct.mapped),
           ByFeature=sum(drop.n.genes),ByMito=sum(drop.pct.mt), 
           BySpike=sum(drop.pct.ercc), Remaining=ncol(QC))
scone_norm<-read.table("/Users/yingy_adm/scone_norm.txt",header = T, sep ="\t",check.names = F)
assay(QC, "lognorm") <- as.matrix(log2(scone_norm+1))
#QCclone <- scdata[rowData(scdata)$use, !colData(scdata)$fclone_drop ]
#visheatmap(cov(drop))
```

#Identifying confounding factors
```{r}
endog_genes <- !rowData(QC)$is_feature_control
plotPCA(
    QC[endog_genes, ],
    exprs_values = "logcounts_raw",
    shape_by = "specificity",
    colour_by = "source",
    size_by = "total_features"
)

plotPCA(
    QC[endog_genes, ],
    exprs_values = "lognorm",
    shape_by = "specificity",
    colour_by = "source",
    size_by = "total_features"
)

plotQC(
    QC[endog_genes, ],
    type = "find-pcs",
    exprs_values = "logcounts_raw",
    variable = "total_features"
)

plotQC(
    QC[endog_genes, ],
    type = "expl",
    exprs_values = "logcounts_raw",
    variables = c(
        "total_features",
        "total_counts",
        "plate",
        "source",
        "cell_specificity",
        "pct_counts_ERCC",
        "pct_counts_MT"
    )
)
plotQC(
    QC[endog_genes, ],
    type = "expl",
    exprs_values = "lognorm",
    variables = c(
        "total_features",
        "total_counts",
        "plate",
        "source",
        "cell_specificity",
        "pct_counts_ERCC",
        "pct_counts_MT"
    )
)
```
scater_gui(QC)
#normalization using scran all genes SF
```{r}
#before normalization
endog_genes <- !rowData(QC)$is_feature_control
plotPCA(
    QC[endog_genes, ],
    exprs_values = "logcounts_raw",
    colour_by = "plate",
    size_by = "total_features",
    shape_by = "source"
)
qclust <- quickCluster(QC, min.size = 15)
QC <- computeSumFactors(QC, sizes = 15, clusters = qclust,positive=TRUE)
QC <- normalize(QC)
#after normalization
plotPCA(
    QC[endog_genes, ],
    exprs_values = "logcounts",
    colour_by = "plate",
    size_by = "total_features",
    shape_by = "source"
)
plotRLE(
    QC[endog_genes, ], 
    exprs_mats = list(Raw = "logcounts_raw", scran = "logcounts"),
    exprs_logged = c(TRUE, TRUE),
    colour_by = "plate"
)

```

#normalization using ERCC_SF
```{r}
#QC<-computeSpikeFactors(QC, type="ERCC", general.use=FALSE)
#QC <- normalize(QC,return_log = FALSE) #normcounts

norm_sf<-function(data)
{ind_ercc<-grep("^ERCC-", rownames(data), perl=TRUE, value=FALSE)
sf_tech <- estimateSizeFactorsForMatrix(data[ind_ercc,])
sf_bio <- estimateSizeFactorsForMatrix(data[-ind_ercc,])
nm_tech <- t(t(data[ind_ercc,]) / sf_tech )
nm_bio <- t(t(data[-ind_ercc,]) / sf_bio)
nm_count<-rbind(nm_tech,nm_bio)
return(nm_count)}

ERCCnorm<-norm_sf(counts(QC))
ERCCnorm<-ERCCnorm[match(rownames(QC), rownames(ERCCnorm)),]
assay(QC,"ERCCnorm")<-ERCCnorm
assay(QC,"log_ERCCnorm")<-log2(ERCCnorm+1)

#PCA after normalization
endog_genes <- !rowData(QC)$is_feature_control
plotPCA(
    QC[endog_genes, ],
    exprs_values = "log_ERCCnorm",
    colour_by = "plate",
    size_by = "total_features",
    shape_by = "source"
)
plotRLE(
    QC[endog_genes, ], 
    exprs_mats = list(logRaw = "logcounts_raw", ERCCnorm = "log_ERCCnorm"),
    exprs_logged = c(TRUE, TRUE),
    colour_by = "plate"
)
```

```{r}
#rowMeans and rowVars to plot the relationship between mean expression and variance for all genes in this dataset. (Hint: use log=“xy” to plot on a log-scale).
plot(rowMeans(counts(QC)),rowVars(counts(QC)),log="xy")
#fit a mean-dependent trend to the variances of the log-expression values for the spike-in transcripts
var.fit <- trendVar(QC, parametric=TRUE, span=0.2,assay.type="counts")
var.out <- decomposeVar(QC, var.fit)
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
    ylab="Variance of log-expression")
curve(var.fit$trend(x), col="dodgerblue", lwd=2, add=TRUE)
cur.spike <- isSpike(QC)
points(var.out$mean[cur.spike], var.out$total[cur.spike], col="red", pch=16)

chosen.genes <- order(var.out$bio, decreasing=TRUE)[1:20]
plotExpression(QC, features=rownames(var.out)[chosen.genes],colour_by="source")

bulk_sig.genes_up<-c("PDCD1","ENTPD1","TNFRSF4","IL21","CXCL13")
ul<-c("PDCD1","ENTPD1","TNFRSF4","IL21","CXCL13")%in%rownames(QC)
plotExpression(QC[,QC$source=="Gut"],features=bulk_sig.genes_up[l],colour_by="specificity")
bulk_sig.genes_down<-c("NT5E","SELL","TNFRSF9","CCR4","KLRG1","CCL22")
dl<-c("NT5E","SELL","TNFRSF9","CCR4","KLRG1","CCL22")%in%rownames(QC)
plotExpression(QC[,QC$source=="Gut"], features=bulk_sig.genes_down[l],colour_by="specificity")

Gup_gene_mat<-melt(tpm(QC)[bulk_sig.genes_up[ul],QC$source=="Gut"])
Gup_gene_mat$log2tpm<-log2(Gup_gene_mat$value)
colnames(Gup_gene_mat)[1]<-"gene_name"
Gup_gene_mat<-cbind(Gup_gene_mat,specificity=QC$specificity[QC$source=="Gut"])

ggplot(Gup_gene_mat, aes(x=gene_name, y=log2tpm, fill=specificity)) +
  geom_violin()

Gdown_gene_mat<-melt(tpm(QC)[bulk_sig.genes_down[dl],QC$source=="Gut"])
Gdown_gene_mat$log2tpm<-log2(Gdown_gene_mat$value)
colnames(Gdown_gene_mat)[1]<-"gene_name"
Gdown_gene_mat<-cbind(Gdown_gene_mat,specificity=QC$specificity[QC$source=="Gut"])

ggplot(Gdown_gene_mat, aes(x=gene_name, y=log2tpm, col=specificity)) +
  geom_violin()+geom_jitter(shape=16, position = "stack" )
```
#M3drop and Brennecke for selecting genes
```{r}
library(M3Drop)
c_label<-paste(colData(QC)$source,colData(QC)$specificity,sep="-")
#filter and norm using CPM
Normalized_data<-M3DropCleanData(counts(QC),labels = c_label, is.counts = T)
#fits <- M3DropDropoutModels(Normalized_data$data)
#pdf("/Users/yingy_adm/Desktop/M3fit.pdf",height=15,width=15)
#dev.off()
#pdf("/Users/yingy_adm/Desktop/M3DE.pdf",height=15,width=15)
#M3_DE_genes <- M3DropDifferentialExpression(Normalized_data$data,mt_method="fdr",mt_threshold=0.005)
M3_DE_genes <- M3DropDifferentialExpression(ERCCnorm,mt_method="fdr",mt_threshold=0.05)
#expr_matrix <- Normalized_data$data # Normalized & filtered expression matrix
#celltype_labs <- factor(Normalized_data$labels) # filtered cell-type labels
#cell_colors <- brewer.pal(max(3,length(unique(celltype_labs))), "Set3")
#pdf("/Users/yingy_adm/Desktop/Br_DE.pdf",height=15,width=15)
Brennecke_HVG <- BrenneckeGetVariableGenes(
    #Normalized_data$data,# scran ERCC sf normalized data
    ERCCnorm,
    fdr = 0.05,
    minBiolDisp = 0.5
)
```
#Correlated Expression for selecting genes
```{r}
cor_mat <- cor(t(expr_matrix), method="spearman") #Gene-gene correlations
cor_mat <- cor(t(normcounts(QC)), method="spearman")
diag(cor_mat) <- rep(0, times=nrow(expr_matrix))
score <- apply(cor_mat, 1, function(x) {max(abs(x))}) #Correlation of highest magnitude
names(score) <- rownames(expr_matrix);
score <- score[order(-score)]
Cor_genes = names(score[1:1500])
```

#PCA for selecting genes
```{r}
reducedDim(QC, "PCA") <- prcomp(t(logcounts(QC)), scale. = TRUE)$x
pc1 <- reducedDim(QC, "PCA")[,1]
design <- model.matrix(~pc1)
library(limma)
fit <- lmFit(logcounts(QC), design)
fit <- eBayes(fit, trend=TRUE, robust=TRUE)
topTable(fit)
```
#PCA on selected genes and qc cells
plotPCA(
    QC[DE_genes[,1],],
    exprs_values = "logcounts",
    colour_by = "plate",
    size_by = "total_features",
    shape_by = "source"
)
red.pca<-prcomp(t(normcounts(QC[Brennecke_HVG,])),retx = TRUE)
#select significant components
sigPC<-permutationPA(logcounts(QC[DE_genes[,1],]), B = 100, threshold = 0.05, verbose = TRUE) # 6significant PCs
#score from the 6 PCs
pc6<-t(red.pca$x[,1:6])

plotReducedDim(QC, "PCA",ncomponents=4,colour_by="plate",shape_by = "source")
#check selected genes
```{r}
pdf("/Users/yingy_adm/Desktop/M3_gene_heatmap.pdf",height=15,width=15)
M3DropExpressionHeatmap(
    M3_DE_genes[,1],
    ERCCnorm,
    cell_labels = paste(QC$source,QC$specificity)
)
dev.off()

pdf("/Users/yingy_adm/Desktop/Brennecke_HVG_cran_heatmap.pdf",height=15,width=15)
M3DropExpressionHeatmap(
    Brennecke_HVG,
    log2(ERCCnorm+1),
    cell_labels = paste(QC$source,QC$specificity,sep = "-")
)
dev.off()

pdf("/Users/yingy_adm/Desktop/Cor_genes_cran_erccnorm.pdf",height=15,width=15)
M3DropExpressionHeatmap(
    Cor_genes,
    expr_matrix,
    cell_labels = celltype_labs
)
dev.off()
```

```{r}
endog_genes <- !rowData(QC)$is_feature_control
plotPCA(
    QC[endog_genes, ],
    exprs_values ="logcounts",
    colour_by = "plate",
    shape_by = "cell_specificity",
    size_by = "total_features"
)
plotPCA(
    QC[endog_genes, ],
    exprs_values ="logtpm",
    colour_by = "plate",
    shape_by = "cell_specificity",
    size_by = "total_features"
)
```
#add M3drop normalized DATA
```{r}
filter_genes<-rownames(scdata)%in%rownames(Normalized_data$data)
rowData(scdata)$use <- filter_genes
QC <- scdata[rowData(scdata)$use, !colData(scdata)$drop]
norm_exprs(QC) <-log2(Normalized_data$data+0.1)
endog_genes <- !rowData(QC)$is_feature_control
plotPCA(
    QC[endog_genes, ],
    exprs_values ="norm_exprs",
    colour_by = "plate",
    shape_by = "source",
    size_by = "total_features"
)
norm_exprs(QC) <-Normalized_data$data
```
# calculate error models using scde
```{r}
library(scde)
#t.ind<-which(celltype[colnames(Normalized_data$data),]$cell_type=="T")#176 "T"out of 181
#groups<-factor(celltype$sepcificity[!colData(scdata)$drop][t.ind],levels = c("Specific", "Non-specific"))
#names(groups) <- rownames(celltype)[!colData(scdata)$drop][t.ind]
#groups<-celltype_labs
groups<-factor(colData(QC)$source,levels = c("Blood", "Gut"))
names(groups)<-colData(QC)$well
table(groups)
counts<-apply(counts(QC),2,function(x) {storage.mode(x) <- 'integer'; x}) 
#o.ifm <- scde.error.models(counts = counts[,t.ind], groups = groups, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
o.ifm <- scde.error.models(counts = counts, groups = groups, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
# estimate gene expression prior
o.prior <- scde.expression.prior(models = o.ifm, counts = counts[,t.ind], length.out = 400, show.plot = FALSE)
```

#DE_source by MAST
```{r}
library(MAST)
#counts <- out_norm
counts <- logcounts(QC)
fData = data.frame(primerid=rownames(counts),name=rowData(QC)$gene_name)
rownames(fData) = rownames(counts)

groups<-factor(QC$source)
cData = data.frame(wellKey=QC$well,cond=groups)
rownames(cData) = colnames(counts)

obj <- FromMatrix(as.matrix(counts), cData, fData)
colData(obj)$cngeneson <- scale(colSums(assay(obj)>0))
cond <- factor(colData(obj)$cond)

# Model expression as function of condition & number of detected genes
zlmCond <- zlm.SingleCellAssay(~cond + cngeneson, obj,method = "bayesglm", 
                           ebayes = TRUE, ebayesControl = list(method = "MLE", model = "H1"))
slotNames(zlmCond)
lrt<-lrTest(zlmCond, "cond")
logfc_source<-getLogFC(zlmCond)
summaryCond <- summary(zlmCond, doLRT="condGut")
summaryDt <- summaryCond$datatable

summaryDt <- as.data.frame(summaryDt)
pVals <- unlist(summaryDt[summaryDt$component == "H",4]) # H = hurdle model
names(pVals) <- unlist(summaryDt[summaryDt$component == "H",1])
pVals <- p.adjust(pVals, method = "fdr")
sig_source_gene<-names(pVals[pVals<0.01])

# Still need to merge on logFC
res_gene <- data.table(melt(lrt))
res_gene[,primerid:=gsub("\\.\\d+","",primerid)]
res_gene <- merge(res_gene, fData(obj), by="primerid")
res_gene_hurdle <- res_gene[metric=="Pr(>Chisq)" & test.type=="hurdle"]
res_gene_hurdle[,adj:=p.adjust(value,"fdr")]
lfc<-getLogFC(zlmCond)[contrast=="condGut"]
setkey(lfc,primerid)
setkey(res_gene_hurdle,primerid)
res_gene_hurdle<-merge(lfc,res_gene_hurdle)
nrow(res_gene_hurdle[adj<0.01])
```
#DE specific vs unspecific in Blood_T
```{r}
#Bcounts<-counts[,QC$source=="Blood"]
Bcounts<-out_norm[,QC$source=="Blood"]
groups<-factor(QC$specificity[QC$source=="Blood"])
BcData = data.frame(wellKey=QC$well[QC$source=="Blood"],cond=groups)
rownames(BcData) = colnames(counts)[QC$source=="Blood"]
Bobj <- FromMatrix(as.matrix(Bcounts), BcData, fData)
colData(Bobj)$cngeneson <- scale(colSums(assay(Bobj)>0))
Bcond <- factor(colData(Bobj)$cond)
BzlmCond <- zlm.SingleCellAssay(~Bcond + cngeneson, Bobj) 

BsummaryCond <- summary(BzlmCond, doLRT="BcondSpecific")
BsummaryDt <- BsummaryCond$datatable
BsummaryDt <- as.data.frame(BsummaryDt)
BpVals <- unlist(BsummaryDt[BsummaryDt$component == "H",4]) # H = hurdle model
names(BpVals) <- unlist(BsummaryDt[BsummaryDt$component == "H",1])
BpVals <- p.adjust(BpVals, method = "fdr")
sig_gene_B<-names(BpVals[BpVals<0.01])

#
# Still need to merge on logFC
Blrt<-lrTest(BzlmCond, "Bcond")
Bres_gene <- data.table(melt(Blrt))
Bres_gene[,primerid:=gsub("\\.\\d+","",primerid)]
#Bres_gene <- merge(Bres_gene, fData(Bobj), by="primerid")
Bres_gene_hurdle <- Bres_gene[Bres_gene$metric=="Pr(>Chisq)" & Bres_gene$test.type=="hurdle",]
Bres_gene_hurdle[,adj:=p.adjust(value,"fdr")]
Blfc<-getLogFC(BzlmCond)[contrast=="BcondSpecific"]
setkey(Blfc,primerid)
setkey(Bres_gene_hurdle,primerid)
Bres_gene_hurdle<-merge(Blfc,Bres_gene_hurdle)
nrow(Bres_gene_hurdle[adj<0.01])#24
#write.table(Bres_gene_hurdle,file.path("~/Desktop","Bspe.sig_hurdle.txt"),row.names = T,col.names = T)

Blogfc<-Blogfc[Blogfc$contrast=="BcondSpecific"&!is.na(Blogfc$logFC),]
Bp<-data.frame(geneid=names(BpVals),pvalue=BpVals)
B_logfc_p<-merge(Blogfc, Bp, by.x="primerid",by.y="geneid")

logfc_source<-getLogFC(BzlmCond)

#volcano plot 
with(B_logfc_p, plot(logFC, -log10(pvalue), pch=20, main="Volcano plot"))

# Add colored points: red if pvalue<0.01, orange of log2FC>1.5, green if both)
with(subset(B_logfc_p, pvalue<.01), points(logFC, -log10(pvalue), pch=20, col="red"))
with(subset(B_logfc_p, abs(logFC)>1.5), points(logFC, -log10(pvalue), pch=20, col="orange"))
with(subset(B_logfc_p, pvalue<.01 & abs(logFC)>1.5), points(logFC, -log10(pvalue), pch=20, col="green"))

# Label points with the textxy function from the calibrate plot
library(calibrate)
with(subset(B_logfc_p, pvalue<.01 | abs(logFC)>1.5), textxy(logFC, -log10(pvalue), labs=primerid, cex=.6))
```

#DE specific vs unspecific in Gut_T
```{r}
#Gcounts<-counts[,QC$source=="Gut"]
Gcounts<-out_norm[,QC$source=="Gut"]
groups<-factor(QC$specificity[QC$source=="Gut"])
GcData = data.frame(wellKey=QC$well[QC$source=="Gut"],cond=groups)
rownames(GcData) = colnames(counts)[QC$source=="Gut"]
Gobj <- FromMatrix(as.matrix(Gcounts), GcData, fData)
colData(Gobj)$cngeneson <- scale(colSums(assay(Gobj)>0))
Gcond <- factor(colData(Gobj)$cond)
GzlmCond <- zlm.SingleCellAssay(~Gcond + cngeneson, Gobj) 

# Still need to merge on logFC
Glrt<-lrTest(GzlmCond, "Gcond")
Gres_gene <- data.table(melt(Glrt))
Gres_gene[,primerid:=gsub("\\.\\d+","",primerid)]
#Gres_gene <- merge(Gres_gene, fData(Gobj), by="primerid")
Gres_gene_hurdle <- Gres_gene[Gres_gene$metric=="Pr(>Chisq)" & Gres_gene$test.type=="hurdle",]
Gres_gene_hurdle[,adj:=p.adjust(value,"fdr")]
Glfc<-getLogFC(GzlmCond)[contrast=="GcondSpecific"]
setkey(Glfc,primerid)
setkey(Gres_gene_hurdle,primerid)
Gres_gene_hurdle<-merge(Glfc,Gres_gene_hurdle)
nrow(Gres_gene_hurdle[adj<0.01])#14
#write.table(Gres_gene_hurdle,file.path("~/Desktop","Gspe.sig_hurdle.txt"),row.names = T,col.names = T)

Glogfc<-getLogFC(GzlmCond)
GsummaryCond <- summary(GzlmCond, doLRT="GcondSpecific")
GsummaryDt <- GsummaryCond$datatable
GsummaryDt <- as.data.frame(GsummaryDt)
GpVals <- unlist(GsummaryDt[GsummaryDt$component == "H",4]) # H = hurdle model
names(GpVals) <- unlist(GsummaryDt[GsummaryDt$component == "H",1])
GpVals <- p.adjust(GpVals, method = "fdr")
sig_gene_G<-names(GpVals[GpVals<0.01])

Glogfc<-Glogfc[Glogfc$contrast=="GcondSpecific"&!is.na(Glogfc$logFC),]
Gp<-data.frame(geneid=names(GpVals),pvalue=GpVals)
G_logfc_p<-merge(Glogfc, Gp, by.x="primerid",by.y="geneid")

#volcano plot 
with(G_logfc_p, plot(logFC, -log10(pvalue), pch=20, main="Volcano plot"))

# Add colored points: red if pvalue<0.01, orange of log2FC>1.5, green if both)
with(subset(G_logfc_p, pvalue<.01), points(logFC, -log10(pvalue), pch=20, col="red"))
with(subset(G_logfc_p, abs(logFC)>1.5), points(logFC, -log10(pvalue), pch=20, col="orange"))
with(subset(G_logfc_p, pvalue<.01 & abs(logFC)>1.5), points(logFC, -log10(pvalue), pch=20, col="green"))

down<-bulk_down$gene_name[bulk_down$gene_name%in%Gres_gene_hurdle$primerid]#130
up<-bulk_up$gene_name[bulk_up$gene_name%in%Gres_gene_hurdle$primerid]#45
with(Gres_gene_hurdle, plot(logFC, -log10(adj), pch=20, main="Volcano plot in gut"))
with(Gres_gene_hurdle[Gres_gene_hurdle$primerid%in%up,], points(logFC, -log10(adj), pch=20, col="red"))
with(Gres_gene_hurdle[Gres_gene_hurdle$primerid%in%down,], points(logFC, -log10(adj), pch=20, col="green"))
# Label points with the textxy function from the calibrate plot
library(calibrate)
with(subset(G_logfc_p, pvalue<.01 | abs(logFC)>1.5), textxy(logFC, -log10(pvalue), labs=primerid, cex=.6))
```

```{r}
logfc_source<-logfc_source[logfc_source$contrast=="condGut"&!is.na(logfc_source$logFC),]
p<-data.frame(geneid=names(pVals),pvalue=pVals)
logfc_p_source<-merge(logfc_source, p, by.x="primerid",by.y="geneid")

#volcano plot 
with(logfc_p_source, plot(logFC, -log10(pvalue), pch=20, main="Volcano plot"))

# Add colored points: red if pvalue<0.01, orange of log2FC>2, green if both)
with(subset(logfc_p_source, pvalue<.0001), points(logFC, -log10(pvalue), pch=20, col="red"))
with(subset(logfc_p_source, abs(logFC)>2), points(logFC, -log10(pvalue), pch=20, col="orange"))
with(subset(logfc_p_source, pvalue<.0001 & abs(logFC)>2000), points(logFC, -log10(pvalue), pch=20, col="green"))

# Label points with the textxy function from the calibrate plot
library(calibrate)
with(subset(logfc_p_source, pvalue<.0001 & abs(logFC)>), textxy(logFC, -log10(pvalue), labs=primerid, cex=.6))
```
#DE by wilcox.test
#keep <- anno[,1] == "NA19101" | anno[,1] == "NA19239"
#data <- molecules[,keep]
group <- colData(QC)$source
lib_size = colSums(counts(QC))
norm <- t(t(counts(QC))/lib_size * median(lib_size)) 
wilcox <- apply(
    norm, 1, function(x) {
        wilcox.test(
            x[group == "Blood"], 
            x[group == "Gut"]
        )    })
wilcox.pVals<-wilcox$p.value

# multiple testing correction
pVals <- p.adjust(pVals, method = "fdr")
DE_Quality_AUC(pVals)
###1 QC PLOT FOR TRYING PARAMETERS
plotQC(scdata, type = "exprs-freq-vs-mean")

quality_control_plot<-function(scdata,min_total_reads, min_total_features)
  {
  #Histogram of library sizes for all cells
  hist(
    scdata$total_counts,
    breaks = 100, main = "Histogram of library sizes for all cells")
  abline(v = min_total_reads, col = "red")
  
  filter_by_total_counts <- (scdata$total_counts > min_total_reads)
  print(as.data.frame(table(filter_by_total_counts)))
  
  hist(scdata$total_features,
    breaks = 100,
    main = "Histograms of number of expressed genes for all cells" )
  abline(v = min_total_features, col = "red")
  
  filter_by_expr_features <- (scdata$total_features > min_total_features)
   print (as.data.frame(table(filter_by_expr_features)))
}

quality_control_plot(scdata=QCmet,min_total_reads=30000, min_total_features=2000)

#more QC plot
hist(pData(QCmet)$pct_counts_feature_controls_MT, xlab="Mitochondrial proportion (%)",
           ylab="Number of cells", breaks=20, main="", col="grey80")

hist(pData(QCmet)$pct_counts_feature_controls_ERCC, xlab="ERCC proportion (%)",
           ylab="Number of cells", breaks=20, main="", col="grey80")

hist(log10(fData(QCmet)$mean_exprs), breaks=100, main="Histogram of log-average counts for all genes", col="grey80",xlab=expression(Log[10]~"average count"))


#execute manual filter 
filter_by_total_counts <- (QCmet$total_counts > 30000)
filter_by_expr_features <- (QCmet$total_features > 2000)
QCmet$use <- (
  # sufficient genes
  filter_by_expr_features &
    # sufficient counted
    filter_by_total_counts )
QCmet$use[87]=TRUE

knitr::kable(
  as.data.frame(table(QCmet$use)),
  booktabs = TRUE,
  row.names = FALSE,
  caption = 'The number of cells removed by manual filter (FALSE)'
  )
  
scater::plotPhenoData(
  QCmet,
  aes_string(x = "total_features",
             y = "pct_counts_feature_controls_MT")
)

scater::plotPhenoData(
  QCmet,
  aes_string(x = "total_features",
             y = "pct_counts_feature_controls_ERCC")
)

scater::plotPhenoData(
  QCmet,
  aes_string(x = "tpm_endogenous_features",
             y = "tpm_feature_controls_ERCC")
  )  
#set_exprs(QCmet,"norm_counts")<-nm_count
#set_exprs(QCmet,"norm_tpm")<-nm_tpm

#plot to check IG genes
#plotPhenoData(QCmet, aesth = aes_string(x = "pct_counts_feature_controls_IG",
#                                        y = "pct_counts_feature_controls_TR", colour = "cell.type"))

#PCA on raw count/endog_genes using exprs_values="exprs/counts/tpm"
endog_genes <- !fData(QCmet)$is_feature_control
scater::plotPCA(QCmet[endog_genes, ], ncomponents = 2,
                colour_by = "cell_type",
                size_by = "total_features",
                shape_by = "plate",
                exprs_values = "tpm")

#PCA on raw count/endog_genes using exprs_values="exprs" #exprs=log2(counts)
endog_genes <- !fData(QCmet)$is_feature_control
scater::plotPCA(QCmet[endog_genes, ], ncomponents = 3,
                colour_by = "cell.type",
                size_by = "total_features",
                shape_by = "batch",
                exprs_values = "exprs")

#TSNE on raw count/endog_genes using exprs_values="exprs/counts/tpm"
plotTSNE(QCmet[endog_genes, ], ncomponents = 2,
         exprs_values = "exprs", colour_by = "cell.type", shape_by ="batch",
         size_by = "total_features", 
         perplexity = 20)

#filter on cells
# libsize.drop <- isOutlier(QCmet$total_counts, nmads=3, type="lower", log=TRUE)
# feature.drop <- isOutlier(QCmet$total_features, nmads=3, type="lower", log=TRUE)
# mito.drop <- isOutlier(QCmet$pct_counts_feature_controls_MT, nmads=3, type="higher") 
# spike.drop <- isOutlier(QCmet$pct_counts_feature_controls_ERCC, nmads=3, type="higher")
# QC <- QCmet[,!(libsize.drop | feature.drop | mito.drop | spike.drop)]
# data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
#            ByMito=sum(mito.drop), BySpike=sum(spike.drop),Remaining=ncol(QC))

#try Random forest??
map_stats<-cbind(sep2016_mapping,nov2016_mapping)
map_stats<-t(map_stats)
map_stats<-data.frame(map_stats)
drop.pct.mapped<-(map_stats$percent_mapped<20)
drop.n.reads<-(QCmet$total_counts<2^15)
drop.n.genes<-(QCmet$total_features<2000)
drop.pct.mt<-(QCmet$pct_counts_feature_controls_MT>6)
drop.pct.ercc<-(QCmet$pct_counts_feature_controls_ERCC>30)
QCmet$use<-!(drop.pct.mapped|drop.n.reads|drop.n.genes|drop.pct.mt|drop.pct.ercc)
data.frame(ByMapPerc=sum(drop.pct.mapped),ByLibSize=sum(drop.n.reads),
           ByFeature=sum(drop.n.genes),ByMito=sum(drop.pct.mt),
           BySpike=sum(drop.pct.ercc), Remaining=sum(QCmet$use))

#Filtering out low-abundance genes
filter_genes <- apply(counts(QCmet[, pData(QCmet)$use]), 1, 
                      function(x) length(x[x > 1]) >= 2)
fData(QCmet)$use <- filter_genes

dim(QCmet[fData(QCmet)$use, pData(QCmet)$use])
# ave.counts <- rowMeans(counts(QC))
# keep <- ave.counts >= 1
# sum(keep)#12119
# ave.tpm <- rowMeans(tpm(QC))
# keep1 <- ave.tpm >= 1
# sum(keep1)#11598
# numcells <- nexprs(QC, byrow=TRUE)
# keep2 <- numcells >= 4
# sum(keep2)#24101
# QC<-QC[keep1,]
# endog <- !fData(QC)$is_feature_control

#excute filter to object
QC <- QCmet[fData(QCmet)$use, pData(QCmet)$use]
endog_genes <- !fData(QC)$is_feature_control

#after QC:PCA on raw count/endog_genes using exprs_values="exprs",i.e. log2(counts)
scater::plotPCA(QC[endog_genes, ], ncomponents = 2,
                colour_by = "cell.type",
                size_by = "total_counts",
                shape_by = "batch",
                exprs_values = "exprs")

#identify principle components that correlate with total_features
scater::plotQC(QC[endog_genes, ],
               type = "find-pcs",
               variable = "cell.type",
               exprs_values = "exprs")

#normalize on libary size 
#function for evaluation
calc_cell_RLE <-
  function (expr_mat, spikes = NULL) 
  {
    RLE_gene <- function(x) {
      if (median(unlist(x)) > 0) {
        log((x + 1)/(median(unlist(x)) + 1))/log(2)
      }
      else {
        rep(NA, times = length(x))
      }
    }
    if (!is.null(spikes)) {
      RLE_matrix <- t(apply(expr_mat[-spikes, ], 1, RLE_gene))
    }
    else {
      RLE_matrix <- t(apply(expr_mat, 1, RLE_gene))
    }
    cell_RLE <- apply(RLE_matrix, 2, median, na.rm = T)
    return(cell_RLE)
  }

#Cell-wise RLE(relative log expression) on raw count
boxplot(calc_cell_RLE(counts(QC[endog_genes, ])),
        col = "grey50",
        ylab = "RLE",
        main = "", ylim=c(-1,1))

norm_sf<-function(data)
{ind_ercc<-grep("^ERCC", rownames(data), perl=TRUE, value=FALSE)
sf_tech <- estimateSizeFactorsForMatrix(data[ind_ercc,])
sf_bio <- estimateSizeFactorsForMatrix(data[-ind_ercc,])
nm_tech <- t(t(data[ind_ercc,]) / sf_tech )
nm_bio <- t(t(data[-ind_ercc,]) / sf_bio)
nm_count<-rbind(nm_tech,nm_bio)
return(nm_count)}

nm_count<-norm_sf(counts(QC))
nm_tpm<-norm_sf(tpm(QC))
nm_count<-data.matrix(nm_count)
nm_tpm<-data.matrix(nm_tpm)
nm_count<-nm_count[rownames(counts(QC)),]
nm_tpm<-nm_count[rownames(tpm(QC)),]
#write.table(nm_count, file = "~/Documents/sc_count/Nnm_count.txt",sep="\t")
#write.table(nm_tpm, file = "~/Documents/sc_count/Nnm_tpm.txt",sep="\t")
set_exprs(QC,"norm_counts")<-nm_count
set_exprs(QC,"norm_tpm")<-nm_tpm

#check after normalization
scater::plotPCA(QC[endog_genes, ],
                shape_by = "batch",
                size_by = "total_features",
                colour_by = "cell.type",
                exprs_values = "norm_counts")
boxplot(calc_cell_RLE(norm_counts(QC[endog_genes, ])),
        col = "grey50",
        ylab = "RLE",
        main = "", ylim=c(-1,1))

#RLE on raw counts
#endog_genes <- !fData(QC)$is_feature_control
QC <- scater::normaliseExprs(QC,method = "RLE",feature_set = endog_genes,
                             exprs_values = "counts",return_norm_as_exprs = TRUE)
scater::plotPCA(QC[endog_genes, ],
                shape_by = "batch",
                size_by = "total_features",
                colour_by = "cell.type",
                exprs_values = "norm_counts")
boxplot(calc_cell_RLE(norm_counts(QC[endog_genes, ])),
        col = "grey50",
        ylab = "RLE",
        main = "", ylim=c(-1,1))
#RLE on tpm???
QC <- scater::normaliseExprs(QC,method = "RLE",feature_set = endog_genes,
                             exprs_values = "tpm",return_norm_as_exprs = TRUE)
scater::plotPCA(QC[endog_genes, ],
                shape_by = "batch",
                size_by = "total_features",
                colour_by = "cell.type",
                exprs_values = "norm_tpm")
boxplot(calc_cell_RLE(norm_tpm(QC[endog_genes, ])),
        col = "grey50",
        ylab = "RLE",
        main = "", ylim=c(-1,1))

#weighted trimmed mean of M-values (TMM)
QC <- scater::normaliseExprs(QC,method = "TMM",feature_set = endog_genes,
                             exprs_values = "counts",return_norm_as_exprs = TRUE)
scater::plotPCA(QC[endog_genes, ],
                shape_by = "batch",
                size_by = "total_features",
                colour_by = "cell.type",
                exprs_values = "norm_counts")
boxplot(calc_cell_RLE(norm_counts(QC[endog_genes, ])),
        col = "grey50",
        ylab = "RLE",
        main = "", ylim=c(-1,1))

QC <- scater::normaliseExprs(QC,method = "TMM",feature_set = endog_genes,
                             exprs_values = "exprs",return_norm_as_exprs = TRUE)
scater::plotPCA(QC[endog_genes, ],
                shape_by = "batch",
                size_by = "total_features",
                colour_by = "cell.type",
                exprs_values = "norm_ex")
boxplot(calc_cell_RLE(norm_tpm(QC[endog_genes, ])),
        col = "grey50",
        ylab = "RLE",
        main = "", ylim=c(-1,1))
#Checking for important technical factors
plotExplanatoryVariables(QC, variables=c("total_counts_ERCC",
                                          "log10_total_counts_ERCC"))

library("SC3")
endog_genes <- !rowData(QC)$is_feature_control
k=sc3_estimate_k(QCclone[endog_genes,])
metadata(k)$sc3$k_estimation
QCclone<- sc3(QCclone[endog_genes, ], ks = 2, biology = TRUE)
sc3_plot_consensus(QCclone, k = 2, show_pdata = c("source","cell_specificity","clonetype"))
sc3_plot_expression(QCclone, k = 2, show_pdata = c("source","cell_specificity","clonetype"))

k=sc3_estimate_k(QC[Brennecke_HVG,])
metadata(k)$sc3$k_estimation #3
QC<- sc3(QC[Brennecke_HVG, ], ks = 3, biology = TRUE)
sc3_plot_consensus(QC, k = 3, show_pdata = c("source","cell_specificity"))
sc3_plot_expression(QC, k = 3, show_pdata = c("source","cell_specificity"))
TQC<-QC[,which()DE_genes]
QC<-sc3_calc_consens(QC)
QC <- sc3(QC[Brennecke_HVG,], ks = 3, biology = TRUE)
sc3_interactive(QC)

a<-data.frame(specificity=colData(QC)$specificity,tpm=tpm(QC)["ENSG00000168685",])
p <- ggplot(a, aes(x=specificity, y=tpm))
p + geom_violin() + geom_boxplot(width=.1, fill="black", outlier.colour=NA) +stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=2.5)

ERCC_count<-log(counts(QC)[grepl("^ERCC-", rownames(QC)),]+1)
EN_count<-log(counts(QC)[!grepl("^ERCC-", rownames(QC)),]+1)

plot(rowMeans(EN_count),rowVars(EN_count), pch=16, cex=0.6, xlab="Mean log-expression",  ylab="Variance of log-expression")
points(rowMeans(ERCC_count),rowVars(ERCC_count),col="red")

T<-as.vector(which(as.numeric(dge[which(rownames(dge)=="CD3D"),])>0|
         as.numeric(dge[which(rownames(dge)=="CD3E"),])>0|
         as.numeric(dge[which(rownames(dge)=="CD3G"),])>0|
         as.numeric(dge[which(rownames(dge)=="CD4"),])>0|
         as.numeric(dge[which(rownames(dge)=="CD8A"),])>0|
         as.numeric(dge[which(rownames(dge)=="CD8B"),])>0)
         
ax = sns.violinplot(x="day", y="total_bill", hue="smoker",
                    data=tips, palette="muted")
