# Load files
extra_metadata = "sample_metadata_batch.tsv"
project_summary = "project-summary.csv"
counts_file = "combined.counts"
summarydata = data.frame(read.table(project_summary, header=TRUE, sep=","), row.names="Name", check.rows=FALSE)
summarydata$Name = rownames(summarydata)
summarydata = summarydata[order(summarydata$Name),]
summarydata$time = as.factor(summarydata$time)
extra_metadata = read.table(extra_metadata, header=T)
summarydata$batch = as.factor(extra_metadata$batch)

counts = read.table(counts_file, header=TRUE, row.names="id", check.names=FALSE)
counts = counts[, order(colnames(counts))]
colnames(counts) = gsub(".counts", "", colnames(counts))
counts <- counts[rowSums(counts>0)>3,]

library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(DESeq2)
library(reshape)
library(ggplot2)
library(CHBUtils)
library(vcd)
library(dplyr)
library(corrplot)
library(pheatmap)
library(gridExtra)
library(tidyr)


#PCA showing batch
dds = DESeqDataSetFromMatrix(counts, colData=summarydata, design=~1)
vst = varianceStabilizingTransformation(dds)
mds(assay(vst), condition = as.factor(summarydata$batch))

# Calculate correlation between metadata
# this  will plot a correlation matrix to
# detect any correlation between metadata and batch
variables = combn(c("treatment", "time", "batch", "concentration"),2,simplify = F)
cor_list = lapply(variables, function(var){
  dd = data.frame(x=summarydata[,var[1]], y=summarydata[,var[2]])
  sm = chisq.test(dd$x, dd$y) # get a p-value
  r= assocstats(xtabs(~x+y,data=dd))$cramer # get a R^2 value for categorical vars
  pval = sm$p.value
  c(r,pval)
})

cor = data.frame(cbind( data.frame(do.call(rbind, variables)),
                        data.frame(do.call(rbind, cor_list)) ))
names(cor) = c("var1", "var2", "cor", "pvalue")

cor_mat = cor[,c(1:3)] %>% tidyr::spread(var1,cor)
rownames = cor_mat[,1]
cor_mat = cor_mat[,2:ncol(cor_mat)]
cor_mat[is.na(cor_mat)] = 0
cor_mat = as.matrix(cor_mat)
row.names(cor_mat) = rownames

pval_mat = cor[,c(1,2,4)] %>% tidyr::spread(var1,pvalue)
rownames = pval_mat[,1]
pval_mat = pval_mat[,2:ncol(pval_mat)]
pval_mat[is.na(pval_mat)] = 1
pval_mat = as.matrix(pval_mat)
row.names(pval_mat) = rownames

corrplot(cor_mat, p.mat = pval_mat, insig = "n",
         is.corr=TRUE)
 ## ??? what does the color and circle size mean in the graph?? why is treatment and concentration so well correlated? Do you correlate the data or just the metadata??

# remove 24, 38, 56 samples
remove = which(summarydata$Name %in%
                 c("24_AGCGATAG-ATAGAGGC", "38_TCCGGAGA-GGCTCTGA", "56_TAATGCGC-AGGCGAAG")
               )

row.names(summarydata) = paste0(summarydata$treatment,"_",summarydata$time,
                                   "_", summarydata$concentration,"_",extra_metadata$sample_number)
names(counts) = row.names(summarydata)

# change gene names and column names. commented since needed to do once only.
# counts_exp = cbind(gene=row.names(counts), counts)
# map_gene = annotate_df(counts_exp, "gene", 'hsapiens_gene_ensembl', "ensembl_gene_id", "external_gene_name")
# map_gene$external_gene_name[is.na(map_gene$external_gene_name)] = map_gene$gene
# map_gene$external_gene_name = paste0(map_gene$external_gene_name,"_",1:nrow(map_gene))
# counts_exp = map_gene[ ,c( ncol(counts_exp)+1, 2:ncol(counts_exp) ) ]
# write.table(counts_exp, "combined.counts.symbol", row.names = F, quote=F, sep="\t")
counts_exp = read.table("combined.counts.symbol", row.names=1, sep="\t", header = T, check.names = F)

colnames(counts_exp)<-rownames(summarydata)
counts_exp<-t(counts_exp)

counts_exp<-transform(counts_exp,sample_number=summarydata$sample_number, time=summarydata$time, concentration=summarydata$concentration)


counts_exp2name<-counts_exp2name[as.matrix(sort(x))]
 
metadata = summarydata[-remove,]

counts_clean = counts_exp[,row.names(metadata)]

write.table(counts_exp,file = 'STARcounts.csv', sep = ',', row.names = TRUE)


### Differential expression by time only for DMSO
samples = row.names(metadata[metadata$treatment=="DMSO",])
# check if correlation between metadata and batch effect
# All conditions have samples from different batch
table(metadata[samples, c("time", "batch")])

dds = DESeqDataSetFromMatrix(countData=counts_clean[,samples],
                             colData=metadata[samples,], design = ~ time)
dds = DESeq(dds, reduced = ~ 1, test = "LRT")
# In case you want to add batch effect
# design(dds) = ~ batch  + time
# dds = DESeq(dds, full = ~ batch +  time, reduced = ~ batch , test = "LRT")
res = results(dds)
# to give a summary of the DE analysis
summary(res)
# plot 10 gene to see if it makes sense
head(res[order(res$padj),])
pp = lapply(10:20, function(i){
  dd = plotCounts(dds, row.names(res[order(res$padj),])[i],
                  intgroup="time", returnData = T)
  dd$treatment = metadata[row.names(dd), "treatment"]
  dd$batch = metadata[row.names(dd), "batch"]
  p=ggplot(dd, aes(x=time,y=count,color=batch,shape=time)) +
    geom_point()
  p
})
do.call(grid.arrange, pp)

##?? can we put gene names on the plots? you seem to use =, instead of <-, are there differences?

# plot MDS
rlog = rlog(dds)
sign = row.names(res[res$padj<0.01 & !is.na(res$padj),])
mds(assay(rlog)[sign,], condition = metadata[samples,"batch"])
# plot heatmap
pheatmap(assay(rlog)[sign,], clustering_distance_cols = "correlation",
         show_rownames = FALSE,
         clustering_method = "ward.D2",
         annotation_col = metadata[,c("time","treatment", "batch")])

# GO enrichment
library(org.Hs.eg.db)
library(clusterProfiler)
.accnum = convertIDs(sign, "ENSEMBL", "ENTREZID", org.Hs.eg.db, "useFirst")
ego <- enrichGO(gene = .accnum[!is.na(.accnum)],
                OrgDb = "org.Hs.eg.db", ont = "BP", pAdjustMethod = "BH",
                pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = TRUE)
# reduce list, you can print the whole list typing summary(ego)
ego = simplify(ego)
summary(ego)[,1:7]

# clustering over time
library(cluster)
rlogMat = assay(rlog)
m = (1-cor(t(rlogMat[sign,]), method = "kendall"))
m[m<0] = 0
d = as.dist(m^2)
c = diana(d, diss = TRUE, stand = FALSE)

groups = cutree(as.hclust(c), h = c$dc)

norm_sign = t(apply(rlogMat[sign,], 1, function(e){
  m = sapply(metadata$treatment, function(i){
    idx = which(metadata$treatment == i)
    mean(e[idx], na.rm=TRUE)
  })
  (e - min(m))/(max(m) - min(m))
}))

# here functions to plot
plot_cluster  = function(norm_sign, g_in_c, metadata, title) {
  ma = as.data.frame(norm_sign)[g_in_c,]
  ma_long = suppressMessages(melt(cbind(gene=row.names(ma), ma),
                                  variable_name = "sample"))
  ma_long$group = metadata[as.character(ma_long$sample), "time"]
  ma_sum = ma_long %>% group_by(gene, group) %>% summarise(average=mean(value)) %>%
    ungroup()

  ma_group = ma_sum %>% group_by(group) %>% summarise(average=median(average)) %>% ungroup()

  ggplot(ma_sum, aes(x=group, y=average)) +
    geom_boxplot() +
    stat_smooth(data=ma_group, aes(x=group, y=average, group=1),se=F,method = "lm",formula = y~poly(x,3)) +
    ggtitle(paste("Group:", title, "(", length(g_in_c), " genes )")) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


to_plot = names(table(groups))[table(groups) > 15]
plots = lapply(to_plot, function(x){
  plot_cluster(norm_sign, as.character(names(groups[groups==x])), metadata[samples,], x)
})

do.call(grid.arrange, plots)

# here Go enrichment of each cluster and put all together in a table

##?? KEGG.db contains mappings based on older data because the original resource was removed from the the public domain before the most recent update
#was produced. This package should now be considered deprecated and future versions of Bioconductor may not have it available.  Users who want
#more current data are encouraged to look at the KEGGREST or reactome.db packages

library(KEGG.db)
cat("\n\n\n")
.void = lapply(to_plot, function(x){
  cat("\n### group:", x)
  cat("\n\n")
  .g = as.character(names(groups[groups==x]))
  .accnum = convertIDs(.g,
                       "ENSEMBL", "ENTREZID", org.Hs.eg.db, "useFirst")
  .symbol = convertIDs(.g,
                       "ENSEMBL", "SYMBOL", org.Hs.eg.db, "useFirst")
  ego <- enrichKEGG(gene = .accnum[!is.na(.accnum)],
                    organism = "human",use_internal_data = TRUE)
  if ( is.data.frame(summary(ego)) )
    if ( nrow(summary(ego)) >0 )  print(kable(summary(ego)[,1:7]))

  cat("\n\n")
  ego <- enrichGO(gene = .accnum[!is.na(.accnum)],
                  OrgDb = "org.Hs.eg.db", ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff = 0.01, qvalueCutoff = 0.1, readable = TRUE)
  print(print_enrichGO(summary(ego), 30))
  cat("\n\n")
  cbind(ensembl=.g, symbol=.symbol, group=x)
})

cluster_table = data.frame(do.call(rbind, .void))

### Differential expression by time for Lapatinib not seen in DMSO
samples = row.names(metadata[(metadata$concentration==3 & metadata$treatment=="Lap") |
                             metadata$treatment=="DMSO",])
# check if correlation between metadata and batch effect
table(metadata[samples, c("treatment", "time", "batch")])

coldata = droplevels(metadata[samples,c("time", "treatment","concentration", "batch")])

full <- model.matrix(~ time + treatment + time:treatment, coldata)
reduced <- full[,1:4]

dds = DESeqDataSetFromMatrix(countData=counts_clean[,samples],
                             colData=metadata[samples,], design = ~ time + treatment)
# if you want to take into account batch effect
# design(dds) = ~ batch + time + treatment + batch:time + batch:treatment + time:treatment
# full <- model.matrix(~ batch + time + treatment + batch:time + batch:treatment + time:treatment, coldata)
# reduced <- full[,grepl("batch", colnames(full)) |
#                 grepl("batch[0-9]:time", colnames(full)) |
#                 grepl("Intercept", colnames(full))]

dds = DESeq(dds, full=full, reduced = reduced, test = "LRT")
res = results(dds)
# to give a summary of the DE analysis
summary(res)    ##?? of all cases, there were a lot of low counts genes, 60%. is this expected?
# plot 10 gene to see if it makes sense
head(res[order(res$padj),])
pp = lapply(1:10, function(i){
  dd = plotCounts(dds, row.names(res[order(res$padj),])[i],
                  intgroup="time", returnData = T)
  dd$treatment = metadata[row.names(dd), "treatment"]
  dd$batch = metadata[row.names(dd), "batch"]
  p=ggplot(dd, aes(x=time,y=count,color=treatment,shape=batch)) +
    geom_point()
  p
})
do.call(grid.arrange, pp)

# plot MDS
rlog = rlog(dds)
sign = row.names(res[res$padj<0.01 & !is.na(res$padj),])
mds_lap<-mds(assay(rlog)[sign,], condition = metadata[samples,"time"])
png(filename = "Lap_DMSO_PCA.png",width=600,height=400,pointsize = 12)
plot(mds_lap)
dev.off()

# plot heatmap
pheatmap(assay(rlog)[sign,], clustering_distance_cols = "correlation",
         show_rownames = FALSE,
         clustering_method = "ward.D2",
         annotation_col = metadata[,c("time","treatment", "batch")])

### Differential expression by time for Sorafenib treatment not seen in DMSO
samples = row.names(metadata[(metadata$concentration==3 & metadata$treatment=="Sora") & metadata$time!=168 |
                               metadata$treatment=="DMSO",])
# check if correlation between metadata and batch effect
table(metadata[samples, c("treatment", "time", "batch")])

coldata = droplevels(metadata[samples,c("time", "treatment","concentration", "batch")])

full <- model.matrix(~ time + treatment + time:treatment, coldata)
reduced <- full[,1:4]

dds = DESeqDataSetFromMatrix(countData=counts_clean[,samples],
                             colData=metadata[samples,], design = ~ time + treatment)
# if you want to take into account batch effect
#design(dds) = ~ batch + time + treatment + batch:time + batch:treatment + time:treatment
#full <- model.matrix(~ batch + time + treatment + batch:time + batch:treatment + time:treatment, coldata)
#reduced <- full[,grepl("batch", colnames(full)) |
#grepl("batch[0-9]:time", colnames(full)) |
#grepl("Intercept", colnames(full))]

dds = DESeq(dds, full=full, reduced = reduced, test = "LRT")
res = results(dds)
# to give a summary of the DE analysis
summary(res)    ##?? of all cases, there were a lot of low counts genes, 60%. is this expected?
# plot 10 gene to see if it makes sense
head(res[order(res$padj),])
pp = lapply(1:10, function(i){
  dd = plotCounts(dds, row.names(res[order(res$padj),])[i],
                  intgroup="time", returnData = T)
  dd$treatment = metadata[row.names(dd), "treatment"]
  dd$batch = metadata[row.names(dd), "batch"]
  p=ggplot(dd, aes(x=time,y=count,color=treatment,shape=batch)) +
    geom_point()
  p
})
do.call(grid.arrange, pp)

# plot MDS
rlog = rlog(dds)
sign = row.names(res[res$padj<0.01 & !is.na(res$padj),])
mds_sora<-mds(assay(rlog)[sign,], condition = metadata[samples,"time"])
png(filename = "Sorafenib_DMSO_PCA_no168.png",width=600,height=400,pointsize = 12)
plot(mds_sora)
dev.off()

# plot heatmap
heatmap_sora<-pheatmap(assay(rlog)[sign,], clustering_distance_cols = "correlation",
         show_rownames = FALSE,
         clustering_method = "ward.D2",
         annotation_col = metadata[,c("time","treatment", "batch")])


### Differential expression by time for Sunitinib treatment not seen in DMSO
samples = row.names(metadata[(metadata$concentration==3 & metadata$treatment=="Sun") |
                               metadata$treatment=="DMSO",])
# check if correlation between metadata and batch effect
table(metadata[samples, c("treatment", "time", "batch")])

coldata = droplevels(metadata[samples,c("time", "treatment","concentration", "batch")])

full <- model.matrix(~ time + treatment + time:treatment, coldata)
reduced <- full[,1:4]

dds = DESeqDataSetFromMatrix(countData=counts_clean[,samples],
                             colData=metadata[samples,], design = ~ time + treatment)
# if you want to take into account batch effect
#design(dds) = ~ batch + time + treatment + batch:time + batch:treatment + time:treatment
#full <- model.matrix(~ batch + time + treatment + batch:time + batch:treatment + time:treatment, coldata)
#reduced <- full[,grepl("batch", colnames(full)) |
#grepl("batch[0-9]:time", colnames(full)) |
#grepl("Intercept", colnames(full))]

dds = DESeq(dds, full=full, reduced = reduced, test = "LRT")
res = results(dds)
# to give a summary of the DE analysis
summary(res)    ##?? of all cases, there were a lot of low counts genes, 60%. is this expected?
# plot 10 gene to see if it makes sense
head(res[order(res$padj),])
pp = lapply(1:10, function(i){
  dd = plotCounts(dds, row.names(res[order(res$padj),])[i],
                  intgroup="time", returnData = T)
  dd$treatment = metadata[row.names(dd), "treatment"]
  dd$batch = metadata[row.names(dd), "batch"]
  p=ggplot(dd, aes(x=time,y=count,color=treatment,shape=batch)) +
    geom_point()
  p
})
do.call(grid.arrange, pp)

# plot MDS
rlog = rlog(dds)
sign = row.names(res[res$padj<0.01 & !is.na(res$padj),])
mds_sun<-mds(assay(rlog)[sign,], condition = metadata[samples,"time"])
png(filename = "Sunitinib_DMSO_PCA.png",width=600,height=400,pointsize = 12)
plot(mds_sun)
dev.off()

# plot heatmap
pheatmap(assay(rlog)[sign,], clustering_distance_cols = "correlation",
         show_rownames = FALSE,
         clustering_method = "ward.D2",
         annotation_col = metadata[,c("time","treatment", "batch")])


### Differential expression by time for Erlotinib treatment not seen in DMSO
samples = row.names(metadata[(metadata$concentration==3 & metadata$treatment=="Erl") |
                               metadata$treatment=="DMSO",])
# check if correlation between metadata and batch effect
table(metadata[samples, c("treatment", "time", "batch")])

coldata = droplevels(metadata[samples,c("time", "treatment","concentration", "batch")])

full <- model.matrix(~ time + treatment + time:treatment, coldata)
reduced <- full[,1:4]

dds = DESeqDataSetFromMatrix(countData=counts_clean[,samples],
                             colData=metadata[samples,], design = ~ time + treatment)
# if you want to take into account batch effect
#design(dds) = ~ batch + time + treatment + batch:time + batch:treatment + time:treatment
#full <- model.matrix(~ batch + time + treatment + batch:time + batch:treatment + time:treatment, coldata)
#reduced <- full[,grepl("batch", colnames(full)) |
#grepl("batch[0-9]:time", colnames(full)) |
#grepl("Intercept", colnames(full))]

dds = DESeq(dds, full=full, reduced = reduced, test = "LRT")
res = results(dds)
# to give a summary of the DE analysis
summary(res)    ##?? of all cases, there were a lot of low counts genes, 60%. is this expected?
# plot 10 gene to see if it makes sense
head(res[order(res$padj),])
pp = lapply(1:10, function(i){
  dd = plotCounts(dds, row.names(res[order(res$padj),])[i],
                  intgroup="time", returnData = T)
  dd$treatment = metadata[row.names(dd), "treatment"]
  dd$batch = metadata[row.names(dd), "batch"]
  p=ggplot(dd, aes(x=time,y=count,color=treatment,shape=batch)) +
    geom_point()
  p
})
do.call(grid.arrange, pp)

# plot MDS
rlog = rlog(dds)
sign = row.names(res[res$padj<0.01 & !is.na(res$padj),])
mds(assay(rlog)[sign,], condition = metadata[samples,"time"])
# plot heatmap
pheatmap(assay(rlog)[sign,], clustering_distance_cols = "correlation",
         show_rownames = FALSE,
         clustering_method = "ward.D2",
         annotation_col = metadata[,c("time","treatment", "batch")])