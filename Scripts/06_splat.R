
library(splatter)
#browseVignettes("splatter")

splat_name <- "TAMU2_splat_set.txt"
genes <- nrow(cnts)
samples <- ncol(cnts)
ndiff<- 0.05

params <- splatEstimate(cnts)
sim <- splatSimulate(params,
                     method = "groups",
                     group.prob = c(0.5, 0.5),
                     de.prob = ndiff,
                     dropout.type = "none")

splat_cnts <- counts(sim)
splat_groups <- sim$Group
class0 <- which(splat_groups == "Group1")
class1 <- which(splat_groups == "Group2")
c0 <- splat_cnts[, class0]
c1 <- splat_cnts[, class1]
splat_cnts <- cbind(c0, c1)

write.table(splat_cnts, splat_name, col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

library(scater)
sim <- logNormCounts(sim)
sim <- runPCA(sim)
plotPCA(sim, colour_by = "Group")
