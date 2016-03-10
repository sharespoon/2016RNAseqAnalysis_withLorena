# Load files
WORKDIR="/Volumes/sysbio/SORGER PROJECTS/people/Sharon Wang/Data/2015-05 RNAseq/2016-02HSPH_RNAseq_analysis/makeCountsTable_forAdam/INPUT/"
setwd(WORKDIR)

extra_metadata = "sample_metadata_batch.tsv"
project_summary = "project-summary.csv"
counts_file = "combined.counts"

summarydata <- data.frame(read.table(project_summary, header=TRUE, sep=","), row.names="Name", check.rows=FALSE)
summarydata$Name = rownames(summarydata)
summarydata = summarydata[order(summarydata$Name),]
extra_metadata = read.table(extra_metadata, header=T)
summarydata$batch = extra_metadata$batch
row.names(summarydata) = paste0(summarydata$treatment,"_",summarydata$time,
                                "_", summarydata$concentration,"_",extra_metadata$sample_number)


counts = read.table(counts_file, header=TRUE, row.names="id", check.names=FALSE)
counts = counts[, order(colnames(counts))]
colnames(counts) = gsub(".counts", "", colnames(counts))
counts <- counts[rowSums(counts>0)>3,]
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

#fix the sample_number, as it is a factor and it has non-numeric values, such as 20C and 80_. Also, fixed treatment, which is a factor. 
summarydata$sample_number=as.character(summarydata$sample_number)
summarydata$sample_number[13]="20"
summarydata$sample_number[78]="80"
summarydata$sample_number=as.numeric(as.matrix(summarydata$sample_number))
summarydata$treatment=as.character(summarydata$treatment)

#add additional columns to counts_exp so that we can re-arrange rows based on the different columns. 
counts_exp<-transform(counts_exp,treatment=summarydata$treatment, sample_number=summarydata$sample_number, time=summarydata$time, concentration=summarydata$concentration)

counts_exp2<-counts_exp[with(counts_exp,order(treatment,concentration,time,sample_number)),]
# this is an alternative way to do shuffling of rows. counts_exp2<-counts_exp[order(counts_exp$treatment,counts_exp$concentration,counts_exp$time, counts_exp$sample_number),]

write.table(counts_exp2,"../OUTPUT/STARcounts.csv",sep = ',',row.names = TRUE)


