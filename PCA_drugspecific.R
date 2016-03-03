# remove 24, 38, 56 samples
remove = which(summarydata$Name %in%
                 c("24_AGCGATAG-ATAGAGGC", "38_TCCGGAGA-GGCTCTGA", "56_TAATGCGC-AGGCGAAG")
)
metadata_clean = summarydata[-remove,]
counts_clean = counts[,-remove]


dds = DESeqDataSetFromMatrix(counts_clean, colData=metadata_clean, design=~1)
vst = varianceStabilizingTransformation(dds)
mds(assay(vst), condition = as.factor(metadata_clean$batch))

condition=which(metadata_clean$treatment %in% c('DMSO'))
metadata_DMSO=metadata_clean[condition,]
counts_DMSO=counts_clean[,condition]
dds_DMSO= DESeqDataSetFromMatrix(counts_DMSO, colData=metadata_DMSO, design=~1)
vst_DMSO = varianceStabilizingTransformation(dds_DMSO)
mds(assay(vst_DMSO), condition = as.factor(metadata_DMSO$batch))

condition=which(metadata_clean$treatment %in% c('DMSO','Sora'))
metadata_DMSO_Sora=metadata_clean[condition,]
counts_DMSO_Sora=counts_clean[,condition]
dds_sora= DESeqDataSetFromMatrix(counts_DMSO_Sora, colData=metadata_DMSO_Sora, design=~1)
vst_sora = varianceStabilizingTransformation(dds_sora)
mds(assay(vst_sora), condition = as.factor(metadata_DMSO_Sora$batch))

condition=which(metadata_clean$treatment %in% c('DMSO','Sun'))
metadata_DMSO_Sun=metadata_clean[condition,]
counts_DMSO_Sun=counts_clean[,condition]
dds_Sun= DESeqDataSetFromMatrix(counts_DMSO_Sun, colData=metadata_DMSO_Sun, design=~1)
vst_Sun = varianceStabilizingTransformation(dds_Sun)
mds(assay(vst_Sun), condition = as.factor(metadata_DMSO_Sun$batch))

condition=which(metadata_clean$treatment %in% c('DMSO','Erl'))
metadata_DMSO_Erl=metadata_clean[condition,]
counts_DMSO_Erl=counts_clean[,condition]
dds_Erl= DESeqDataSetFromMatrix(counts_DMSO_Erl, colData=metadata_DMSO_Erl, design=~1)
vst_Erl = varianceStabilizingTransformation(dds_Erl)
mds(assay(vst_Erl), condition = as.factor(metadata_DMSO_Erl$batch))

condition=which(metadata_clean$treatment %in% c('DMSO','Lap'))
metadata_DMSO_Lap=metadata_clean[condition,]
counts_DMSO_Lap=counts_clean[,condition]
dds_Lap= DESeqDataSetFromMatrix(counts_DMSO_Lap, colData=metadata_DMSO_Lap, design=~1)
vst_Lap = varianceStabilizingTransformation(dds_Lap)
mds(assay(vst_Lap), condition = as.factor(metadata_DMSO_Lap$batch))