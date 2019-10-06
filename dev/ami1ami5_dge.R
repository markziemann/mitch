library("getDEE2")
library("limma")
library("edgeR")
library("DESeq2")
library("ABSSeq")
library("topconfects")
source("nDrich.R")

mdat<-getDee2Metadata("hsapiens")
samples<-mdat[which(mdat$SRP_accession=="SRP096178"),]
samples<-samples[order(samples$SRR_accession),]

samples$group<-factor(c("Ctrl","Ctrl","Ctrl","Ami1","Ami1","Ami1","Ami5","Ami5","Ami5"),levels=c("Ctrl","Ami1","Ami5"))
samples$label<-c("Ctrl_1","Ctrl_2","Ctrl_3","Ami1_1","Ami1_2","Ami1_3","Ami5_1","Ami5_2","Ami5_3")
SRRlist<-as.vector(samples$SRR_accession)
x<-getDEE2("hsapiens",SRRlist)
x<-Tx2Gene(x)
y<-x$Tx2Gene

#here attach genenames
gt<-unique(x$TxInfo[,1:2])

#filter very low genes
y<-y[which(rowSums(y)/ncol(y)>=(10)),]

###########################################
#EdgeR analysis
###########################################
design<-model.matrix(~samples$group)
rownames(design)<-samples$SRR_accession

#DE for AMI1
des<-as.matrix(design[1:6,1:2])
counts<-y[,1:6]
z<-DGEList(counts=counts)
z<-calcNormFactors(z)
z <- estimateDisp(z, des,robust=TRUE)
fit <- glmFit(z, des)
lrt <- glmLRT(fit)
dge<-as.data.frame(topTags(lrt,n=Inf))
dge$dispersion<-lrt$dispersion
dge<-merge(dge,lrt$fitted.values,by='row.names')
rownames(dge)=dge$Row.names
dge$Row.names=NULL
dge<-merge(dge,z$counts,by='row.names')
ami1_edger<-dge[order(dge$PValue),]

#DE for AMI5
des<-as.matrix(design[c(1:3,7:9),c(1,3)])
counts<-y[,c(1:3,7:9)]
z<-DGEList(counts=counts)
z <- estimateDisp(z, des,robust=TRUE)
fit <- glmFit(z, des)
lrt <- glmLRT(fit)
dge<-as.data.frame(topTags(lrt,n=Inf))
dge$dispersion<-lrt$dispersion
dge<-merge(dge,lrt$fitted.values,by='row.names')
rownames(dge)=dge$Row.names
dge$Row.names=NULL
dge<-merge(dge,z$counts,by='row.names')
ami5_edger<-dge[order(dge$PValue),]

save.image("ami1ami5_dge.RData")

###########################################
#now do limma 
###########################################

#AMI1
des<-as.matrix(design[1:6,2])
counts<-y[,1:6]
z<-DGEList(counts=counts)
z <- calcNormFactors(z)
v <- voom(z,des,plot=FALSE)
fit <- lmFit(v, des)
fit.de <- eBayes(fit, robust=TRUE)
dge<-topTable(fit.de,n=Inf)
ami1_limma<-dge[order(dge$P.Value),]

#AMI5
des<-as.matrix(design[c(1:3,7:9),3])
counts<-y[,c(1:3,7:9)]
z<-DGEList(counts=counts)
z <- calcNormFactors(z)
v <- voom(z,des,plot=FALSE)
fit <- lmFit(v, des)
fit.de <- eBayes(fit, robust=TRUE)
dge<-topTable(fit.de,n=Inf)
ami5_limma<-dge[order(dge$P.Value),]

save.image("ami1ami5_dge.RData")

###########################################
#now do ABSSeq
###########################################

#AMI1
des<-design[1:6,1:2]
counts<-y[,1:6]
obj<-ABSDataSet(counts, factor(des[,2]))  #default normalisation is qtotal
obj<-ABSSeq(obj)
dge<- as.data.frame(cbind(obj$Amean,obj$Bmean,obj$foldChange,obj$pvalue,obj$adj.pvalue))
colnames(dge)=c("Amean","Bmean","foldChange","pvalue","adj.pvalue")
ami1_absseq<-dge[order(dge$pvalue),]

#AMI5
des<-design[c(1:3,7:9),c(1,3)]
counts<-y[,c(1:3,7:9)]
obj<-ABSDataSet(counts, factor(des[,2]))  #default normalisation is qtotal
obj<-ABSSeq(obj)
dge<- as.data.frame(cbind(obj$Amean,obj$Bmean,obj$foldChange,obj$pvalue,obj$adj.pvalue))
colnames(dge)=c("Amean","Bmean","foldChange","pvalue","adj.pvalue")
ami5_absseq<-dge[order(dge$pvalue),]

save.image("ami1ami5_dge.RData")


###########################################
#now do DESeq2
###########################################
y<-round(y)

#ami1
des<-samples[1:6,]
des$grp<-as.numeric(grepl("Ami",des$group))
dds <- DESeqDataSetFromMatrix(countData =y[,1:6], colData = des, design = ~ grp)
dds <- DESeq(dds)
z<- DESeq2::results(dds)
vsd <- vst(dds, blind=FALSE)
#stick on the normalised expression values to the table
zz<-cbind(z,assay(vsd))
#sort by adjusted p-value
ami1_deseq2<-as.data.frame(zz[order(zz$pvalue),])

#AMI5
des<-samples[c(1:3,7:9),]
des$grp<-as.numeric(grepl("Ami",des$group))
dds <- DESeqDataSetFromMatrix(countData =y[,c(1:3,7:9)], colData = des, design = ~ grp)
dds <- DESeq(dds)
z<- DESeq2::results(dds)
vsd <- vst(dds, blind=FALSE)
#stick on the normalised expression values to the table
zz<-cbind(z,assay(vsd))
#sort by adjusted p-value
ami5_deseq2<-as.data.frame(zz[order(zz$pvalue),])

save.image("ami1ami5_dge.RData")

###########################################
# now do TopConfects of DeSeq2
###########################################

#ami1
des<-samples[1:6,]
des$grp<-as.numeric(grepl("Ami",des$group))
dds <- DESeqDataSetFromMatrix(countData =y[,1:6], colData = des, design = ~ grp)
dds <- DESeq(dds)
confects <- deseq2_confects(dds, fdr=0.5)
ami1_confects<-confects$table

des<-samples[c(1:3,7:9),]
des$grp<-as.numeric(grepl("Ami",des$group))
dds <- DESeqDataSetFromMatrix(countData =y[,c(1:3,7:9)], colData = des, design = ~ grp)
dds <- DESeq(dds)
confects <- deseq2_confects(dds, fdr=0.5)
ami5_confects<-confects$table

#########################################################
# nDrich analysis
#########################################################

#read in GMT file
genesets<-gmt_import("ReactomePathways.gmt")

#edger analysis
x1<-list("ami1_edger"=ami1_edger,"ami5_edger"=ami5_edger)
y1<-mitch_import(x1,"edger",geneIDcol="Row.names",geneTable=gt)
res<-mitch_calc(y1,genesets,resrows=50)
mitch_plots(res,outfile="ami1ami5_edger.pdf")
mitch_report(res,"ami1ami5_edger.html")

#limma analysis
x2<-list("ami1_limma"=ami1_limma,"ami5_limma"=ami5_limma)
y2<-mitch_import(x2,"limma",geneTable=gt)
res<-mitch_calc(y2,genesets,resrows=50)
mitch_plots(res,outfile="ami1ami5_limma.pdf")
mitch_report(res,"ami1ami5_limma.html")

#absseq analysis
x3<-list("ami1_absseq"=ami1_absseq,"ami5_absseq"=ami5_absseq)
y3<-mitch_import(x3,"absseq",geneTable=gt)
res<-mitch_calc(y3,genesets,resrows=50)
mitch_plots(res,outfile="ami1ami5_absseq.pdf")
mitch_report(res3,"ami1ami5_absseq.html")

#deseq2 analysis
x4<-list("ami1_deseq2"=ami1_deseq2,"ami5_deseq2"=ami5_deseq2)
y4<-mitch_import(x4,"deseq2",geneTable=gt)
res<-mitch_calc(y4,genesets,resrows=50)
mitch_plots(res,outfile="ami1ami5_deseq.pdf")
mitch_report(res,"ami1ami5_deseq.html")

#topconfects analysis
x5<-list("ami1_confects"=ami1_confects,"ami5_confects"=ami5_confects)
y5<-mitch_import(x5,"topconfects",geneIDcol="name",geneTable=gt)
res<-mitch_calc(y5,genesets,resrows=50)
mitch_plots(res,outfile="ami1ami5_confects.pdf")
mitch_report(res,"ami1ami5_confects.html")


