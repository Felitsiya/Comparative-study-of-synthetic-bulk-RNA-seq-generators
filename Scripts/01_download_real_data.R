rm(list = ls())

library(recount3)
library(stringr)
library(matrixStats)

# SRP181886 Huntley data set

human_projects <- available_projects()
proj_info <- subset(human_projects,
  project == "SRP181886" & project_type == "data_sources")

rse_gene <- create_rse(proj_info)

class <- as.data.frame(rse_gene@colData@listData$sra.sample_attributes)
class <- str_split_fixed(class[,1], ";;", 10)
class <- str_split_fixed(class[,4], "[|]", 2)
class <- class[,1]
class <- gsub("'", "", class)

counts <- as.data.frame(rse_gene@assays@data@listData$raw_counts)
colnames(counts) <- class
length(base::which(class == "control"))
length(base::which(class == "Alzheimers disease"))

counts <- counts[rowSums(counts)>0,]
counts <- as.matrix(counts)
counts <- counts[rowMedians(counts)>0,]
attr(counts, "dimnames") <- NULL
anyNA(counts)

AD <- which(class == "Alzheimers disease")
ctrl <- which(class == "control")
AD <- counts[,AD]
ctrl <- counts[,ctrl]
AD <- AD[,1:50]
ctrl <- ctrl[,1:50]
cnts <- cbind(AD, ctrl)

rm(AD, ctrl, counts, class)
length(which(cnts > 10^7))
length(which(cnts > 10^6))
length(as.vector(as.matrix(cnts)))
write.table(cnts, "Huntley_counts.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

#hist(as.vector(as.matrix(cnts)), breaks = 10000)
#hist(as.vector(as.matrix(cnts)), breaks = 10000, xlim = c(1, 10^7))

