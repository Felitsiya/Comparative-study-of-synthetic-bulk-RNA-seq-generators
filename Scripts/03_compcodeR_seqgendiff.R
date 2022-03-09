
library(compcodeR)
library(seqgendiff)
library(matrixStats)
source("D:/compcodeR_test/thin_2group_mod.R")

cnts_compcodeR_name <- "TAMU2_compcodeR_set.txt"
thout2_mat_name <- "TAMU2_seqgendiff_set.txt"

depth_cnts <- mean(rowMeans(cnts))
depth_cnts

genes <- nrow(cnts)
samples <- ncol(cnts)
ndiff<- 0.05

# compcodeR

depth_param <- sum(cnts)/samples

cnts_means <- rowMeans(cnts[,1:samples/2])
disp1_calc <- (rowVars(cnts[,1:samples/2])-
                 rowMeans(cnts[,1:samples/2]))/(rowMeans(cnts[,1:samples/2])^2)
disp2_calc <- (rowVars(cnts[,(samples/2+1):samples])-
                 rowMeans(cnts[,(samples/2+1):samples]))/(rowMeans(cnts[,(samples/2+1):samples])^2)
cnts_disp <- as.matrix(cbind(disp1_calc, disp2_calc))
  
compcodeR_data.obj <- generateSyntheticData(dataset = "mydata",
                                            n.vars = genes,
                                            n.diffexp = round(ndiff*genes), # preserving ratio in TAMU data (500/10000)
                                            samples.per.cond = samples/2,
                                            dispersions = cnts_disp,
                                            fraction.upregulated = 1, #default
                                            seqdepth = depth_param,
                                            minfact = 0.7, #default
                                            maxfact = 1.4, #default
                                            relmeans = cnts_means,
                                            effect.size = 1.5 #default
                                            )
  
cnts_compcodeR <- as.data.frame(compcodeR_data.obj@count.matrix)
write.table(cnts_compcodeR, cnts_compcodeR_name, col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

depthcompcodeR <- mean(rowMeans(cnts_compcodeR))
print(depthcompcodeR)
print(depth_cnts)

# seqgendiff
exp_fun <- function(n,rate, es) {
  stats::rexp(n = n, rate = rate) + es
}
prop_null = 1 - ndiff
n = round(genes * (1 - prop_null))
  
thout2 <- thin_2group_mod(mat = as.matrix(cnts[,51:100]), 
                          prop_null = prop_null, 
                          signal_fun = exp_fun,
                          signal_params = list(n = n, rate = 0.5, es = 1.5), #in compcodeR the effect sizes will be obtained by simulating numbers from an exponential distribution (with rate 1) and adding the results to the effect.size
                          group_prop = 0.5)
print(mean(abs(thout2$coef) < 10^-6))
  
thout2_mat <- thout2$mat
designmat <- thout2$designmat
class0 <- which(thout2$designmat == 0)
class1 <- which(thout2$designmat == 1)
c0 <- thout2_mat[, class0]
c1 <- thout2_mat[, class1]
thout2_mat <- cbind(c0, c1)
write.table(thout2_mat, thout2_mat_name, col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

