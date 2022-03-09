
library(powsimR)
#browseVignettes("powsimR")

tableName <- "Huntley_powsimR_set.txt"

genes <- nrow(cnts)
samples <- ncol(cnts[,51:100])
ndiff<- 0.05

estparam_real <- estimateParam(countData = cnts[,51:100],
                               readData = NULL,
                               batchData = NULL,
                               spikeData = NULL,
                               spikeInfo = NULL,
                               Lengths = NULL,
                               MeanFragLengths = NULL,
                               RNAseq = "bulk",
                               Protocol = "Read",
                               Distribution = "NB",
                               Normalisation = "MR",
                               GeneFilter = 0.1,
                               SampleFilter = 10,
                               sigma = 1.96,
                               verbose = TRUE)

png(file = "powsimR_params1.png", width = 1700, height = 1000)
plotParam(estParamRes = estparam_real, Annot = TRUE)
dev.off()

setupres <- Setup(ngenes = NULL,
                  nsims = 1,
                  p.DE = ndiff,
                  n1 = samples,
                  n2 = samples,
                  Thinning = NULL,
                  LibSize = "given",
                  estParamRes = estparam_real,
                  estSpikeRes = NULL,
                  DropGenes = FALSE,
                  setup.seed = 5299,
                  verbose = TRUE)

simres <- simulateDE(SetupRes = setupres,
                     Prefilter = NULL,
                     Imputation = NULL,
                     Normalisation = "MR",
                     DEmethod = "DESeq2",
                     DEFilter = FALSE,
                     Counts = TRUE,
                     verbose = TRUE)

powsim_cnts <- simres$Counts
powsim_cnts <- powsim_cnts$`50vs50`
powsim_cnts <- powsim_cnts[[1]]

write.table(powsim_cnts, tableName, col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

