
library(countsimQC)
library(DESeq2)
set.seed(123)

genes <- nrow(cnts)
samples <- ncol(cnts)
ndiff<- 0.05

deseq_cnts <- as.data.frame(matrix(as.numeric(unlist(cnts)), ncol = ncol(cnts)))
rownames(cnts) <- sprintf("g%d", 1:genes)
deseq_cnts <- cbind(rownames(cnts), deseq_cnts)
colnames(deseq_cnts)[1] <- "genes"
  
deseq_compcodeR <- read.table("Huntley_compcodeR_set.txt", header = TRUE, sep = "\t")
deseq_compcodeR <- cbind(rownames(deseq_compcodeR), deseq_compcodeR)
colnames(deseq_compcodeR)[1] <- "genes"
  
deseq_seqgendiff <- read.table("Huntley_seqgendiff_set.txt", header = TRUE, sep = "\t")
deseq_seqgendiff <- cbind(rownames(deseq_seqgendiff), deseq_seqgendiff)
colnames(deseq_seqgendiff)[1] <- "genes"
  
deseq_powsimR <- read.table("Huntley_powsimR_set.txt", header = TRUE, sep = "\t")
deseq_powsimR <- cbind(rownames(deseq_powsimR), deseq_powsimR)
colnames(deseq_powsimR)[1] <- "genes"

deseq_simseq <- read.table("Huntley_simseq_set.txt", header = TRUE, sep = "\t")
deseq_simseq <- cbind(rownames(deseq_simseq), deseq_simseq)
colnames(deseq_simseq)[1] <- "genes"

deseq_splat <- read.table("Huntley_splat_set.txt", header = TRUE, sep = "\t")
deseq_splat <- cbind(rownames(deseq_splat), deseq_splat)
colnames(deseq_splat)[1] <- "genes"

deseq_SPsimSeq <- read.table("Huntley_SPsimSeq_set.txt", header = TRUE, sep = "\t")
deseq_SPsimSeq <- cbind(rownames(deseq_SPsimSeq), deseq_SPsimSeq)
colnames(deseq_SPsimSeq)[1] <- "genes"

samp <- sprintf("sample%d", 1:samples)
dex <- as.factor(c(rep("class1", samples/2),rep("class2", samples/2)))
metaData <- cbind(samp, dex)

dds_cnts <- DESeqDataSetFromMatrix(countData = deseq_cnts, colData = metaData, design = ~dex, tidy = TRUE)
dds_c <- DESeqDataSetFromMatrix(countData = deseq_compcodeR, colData = metaData, design = ~dex, tidy = TRUE)
dds_p <- DESeqDataSetFromMatrix(countData = deseq_powsimR, colData = metaData, design = ~dex, tidy = TRUE)
dds_sp <- DESeqDataSetFromMatrix(countData = deseq_splat, colData = metaData, design = ~dex, tidy = TRUE)
dds_SPsimSeq <- DESeqDataSetFromMatrix(countData = deseq_SPsimSeq, colData = metaData, design = ~dex, tidy = TRUE)

samp <- sprintf("sample%d", 1:(samples/2))
dex <- as.factor(c(rep("class1", samples/4),rep("class2", samples/4)))
metaData <- cbind(samp, dex)

dds_s <- DESeqDataSetFromMatrix(countData = deseq_seqgendiff, colData = metaData, design = ~dex, tidy = TRUE)
dds_sim <- DESeqDataSetFromMatrix(countData = deseq_simseq, colData = metaData, design = ~dex, tidy = TRUE)

ddsList <- list(Huntley = dds_cnts,
                compcodeR = dds_c,
                powsimR = dds_p,
                splat = dds_sp,
                SPsimSeq = dds_SPsimSeq,
                seqgendiff = dds_s,
                simseq = dds_sim)


countsimQCReport(ddsList = ddsList,
                 outputFile = "countsim_allH.html",
                 outputDir = "D:/Synthetic data review scripts",
                 outputFormat = "html_document", 
                 showCode = FALSE,
                 forceOverwrite = TRUE,
                 savePlots = TRUE,
                 description = "This is the report for Huntley and all synthetic sets.", 
                 maxNForCorr = 25,
                 maxNForDisp = Inf, 
                 calculateStatistics = TRUE,
                 subsampleSize = 25,
                 kfrac = 0.01,
                 kmin = 5, 
                 permutationPvalues = FALSE,
                 nPermutations = NULL)


### DE analysis
DE_results <- data.frame(DataType = character(),
                 SumNAvalues = numeric(), 
                 SumPadj05 = numeric(), 
                 stringsAsFactors = FALSE) 

dds_cnts <- DESeq(dds_cnts)
dds_c <- DESeq(dds_c)
dds_p <- DESeq(dds_p)
dds_sp <- DESeq(dds_sp)
dds_s <- DESeq(dds_s)
dds_sim <- DESeq(dds_sim)
dds_SPsimSeq <- DESeq(dds_SPsimSeq)


res_cnts <- results(dds_cnts, independentFiltering=FALSE)
res_cnts <- res_cnts[order(res_cnts$padj),]
DE_results[1,1] <- "Huntley"
DE_results[1,2] <- sum(is.na(res_cnts$padj))
DE_results[1,3] <- sum(res_cnts$padj <= 0.05, na.rm = TRUE)

res_c <- results(dds_c, independentFiltering=FALSE)
res_c <- res_c[order(res_c$padj),]
DE_results[2,1] <- "compcodeR"
DE_results[2,2] <- sum(is.na(res_c$padj))
DE_results[2,3] <- sum(res_c$padj <= 0.05, na.rm = TRUE)

res_p <- results(dds_p, independentFiltering=FALSE)
res_p <- res_p[order(res_p$padj),]
DE_results[3,1] <- "powsimR"
DE_results[3,2] <- sum(is.na(res_p$padj))
DE_results[3,3] <- sum(res_p$padj <= 0.05, na.rm = TRUE)

res_sp <- results(dds_sp, independentFiltering=FALSE)
res_sp <- res_sp[order(res_sp$padj),]
DE_results[4,1] <- "splat"
DE_results[4,2] <- sum(is.na(res_sp$padj))
DE_results[4,3] <- sum(res_sp$padj <= 0.05, na.rm = TRUE)

res_s <- results(dds_s, independentFiltering=FALSE)
res_s <- res_s[order(res_s$padj),]
DE_results[5,1] <- "seqgendiff"
DE_results[5,2] <- sum(is.na(res_s$padj))
DE_results[5,3] <- sum(res_s$padj <= 0.05, na.rm = TRUE)

res_sim <- results(dds_sim, independentFiltering=FALSE)
res_sim <- res_sim[order(res_sim$padj),]
DE_results[6,1] <- "simseq"
DE_results[6,2] <- sum(is.na(res_sim$padj))
DE_results[6,3] <- sum(res_sim$padj <= 0.05, na.rm = TRUE)

res_SPsimSeq <- results(dds_SPsimSeq, independentFiltering=FALSE)
res_SPsimSeq <- res_SPsimSeq[order(res_SPsimSeq$padj),]
DE_results[7,1] <- "SPsimSeq"
DE_results[7,2] <- sum(is.na(res_SPsimSeq$padj))
DE_results[7,3] <- sum(res_SPsimSeq$padj <= 0.05, na.rm = TRUE)

write.table(DE_results, "DE_allH.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
