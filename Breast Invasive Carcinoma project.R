####load the needed libraries ########  
library(readr)
library(org.Hs.eg.db)
library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(enrichplot)
library(clusterProfiler)
library(ggplot2)
##########  Preprocessing ##############
#Define your working directory here
setwd("your working directory")
getwd()

#Loading expression and pheno data

data <- read_tsv(file = 'TCGA-BRCA.GDC_phenotype.tsv')
data2 <- read_tsv(file = 'TCGA-BRCA.htseq_counts.tsv')
write.table(data, "file2.csv", 
            sep = ",", # use tab as separator 
            row.names = FALSE) # no row names
write.table(data2, "file.csv", 
            sep = ",", # use tab as separator 
            row.names = FALSE) # no row names

count_data <- as.data.frame( read.csv("file.csv", row.names = 1))
colnames(count_data)
colnames(count_data) <- gsub("[.]", "-", colnames(count_data))
pheno = read.csv("file2.csv")
rownames(pheno) = pheno[,1]

#Filtration of samples with pheno data and ordereing samples in the same order in count and pheno matrices
samples = intersect(pheno$submitter_id.samples , colnames(count_data))
pheno_filt = pheno[samples,]

dim(pheno_filt)
count_filt = count_data[,samples]
dim(count_filt) # Herein the number of samples with phenodata is ??? samples out of ??


View(head(count_filt))
table(pheno_filt$sample_type.samples) #Checking for the number of tumor and normal samples


##Now, we need to convert the data to raw read counts since the DESeq needs raw counts
count_raw=2^(count_filt)-1  ## get the read counts as the data is log2 transformed
View(head(count_raw))

save(count_raw, pheno_filt,file="raw_count_pheno.RDATA")

######### Differential Expression ###############
#Defining the conditions
table(pheno_filt$sample_type.samples)

cond1="Primary Tumor" 
cond2="Solid Tissue Normal"

#converting data to integers and loading it to DESeq followed by differential gene expression
exp.mtx = apply(count_raw,2,as.integer)
rownames(exp.mtx) = rownames(count_raw)

dds = DESeqDataSetFromMatrix( countData = exp.mtx, colData = pheno_filt , design = ~ sample_type.samples)
dds.run = DESeq(dds)
res=results(dds.run, contrast = c("sample_type.samples",cond1 ,cond2) )

# remove nulls
res=res[complete.cases(res), ]
summary(res)


res.df=as.data.frame(res)
res.df.ordered = res.df[order(res.df$padj, decreasing = F),]

write.table(res.df.ordered, file = "your_file_output.csv", sep = ",")

plotMA(res, ylim=c(-2,2)) 


res.degs_half=res.df.ordered[res.df.ordered$padj< 0.05 & abs(res.df.ordered$log2FoldChange)>1.5,] # results at log fold change of 1
res.degs_1 =res.df.ordered[res.df.ordered$padj< 0.05 & abs(res.df.ordered$log2FoldChange)>1,] #results at log fold change 1
dim(res.degs_half) 
dim(res.degs_1) 

degs = rownames(res.degs_half)
write.table(res.degs_half,"degs.csv", sep = ",")
write.csv(degs, file =  "DEGs.csv")

############ Normalization and Heatmap #########################
#### get the normalized and loggedtransformed values of all exp data
ntd=normTransform(dds)
exp.norm= assay(ntd)

#Getting names of the top  20differentially expressed genes
top20genes.names = rownames(res.degs_half[1:20,])

#top100genes = exp.norm[top100genes.names,]
top20genes = exp.norm[top20genes.names,]

##Complex Heatmap
column_ha = HeatmapAnnotation(Sample_Type = pheno_filt$type, Tissue_Site = pheno_filt$tissue_source_site, col = list(Sample_Type = c("Normal" = "red" , "Tumor" = "blue")) )
column_ha = HeatmapAnnotation(Sample_Type = pheno_filt$sample_type.samples,  col = list(Sample_Type = c("Primary Tumor" = "red" , "Solid Tissue Normal" = "blue", "Metastatic" = "black")) )

Heatmap(top20genes, col = colorRamp2(c(min(top20genes), mean(top20genes), max(top20genes)),c("blue", "white", "red")),name = "Exp", top_annotation = column_ha, cluster_columns = F, row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7))

