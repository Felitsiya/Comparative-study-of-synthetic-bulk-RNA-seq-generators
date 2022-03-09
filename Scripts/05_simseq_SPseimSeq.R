
# SimSeq
library(SimSeq)
library(SPsimSeq)

simseq_name <- "TAMU2_simseq_set.txt"
SPsimSeq_name <- "TAMU1_SPsimSeq_set.txt"

genes <- nrow(cnts)
samples <- ncol(cnts)
ndiff<- 0.05

treatment <- as.factor(c(rep("class0", 50), rep("class1", 50)))
simseq_df <- SimData(counts = cnts,
                     treatment = treatment,
                     sort.method = "unpaired",
                     k.ind  = samples/4,
                     n.genes = genes,
                     n.diff = round(ndiff*genes))

simseq_cnts <- simseq_df$counts
simseq_classes <- simseq_df$treatment

write.table(simseq_cnts, simseq_name, col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)


# SPsimSeq

rownames(cnts) <- sprintf("g%d", 1:genes)
colnames(cnts) <- sprintf("s%d", 1:samples)
treatment <- c(rep(0, 50), rep(1, 50))
gc()
set.seed(123)
sim.data.bulk <- SPsimSeq(n.sim = 1,
                          s.data = cnts,
                          group = treatment,
                          n.genes = genes,
                          batch.config = 1,
                          group.config = c(0.5, 0.5),
                          tot.samples = samples, 
                          pDE = ndiff,
                          genewiseCor = FALSE,
                          lfc.thrld = 0.5,
                          result.format = "list",
                          return.details = FALSE)

sim.data.counts <- sim.data.bulk[[1]]$counts
sim.data.groups <- sim.data.bulk[[1]]$colData$Group

write.table(sim.data.counts, SPsimSeq_name, col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
