
library(ggplot2)
library(ggpubr)

cnts_compcodeR <- read.table("TAMU2_compcodeR_set.txt", header = TRUE, sep = "\t")
powsim_cnts <- read.table("TAMU2_powsimR_set.txt", header = TRUE, sep = "\t")
seqgendiff_cnts <- read.table("TAMU2_seqgendiff_set.txt", header = TRUE, sep = "\t")
simseq_cnts <- read.table("TAMU2_simseq_set.txt", header = TRUE, sep = "\t")
SPsimSeq_cnts <- read.table("TAMU2_SPsimSeq_set.txt", header = TRUE, sep = "\t")
splat_cnts <- read.table("TAMU2_splat_set.txt", header = TRUE, sep = "\t")
  
dim(cnts)
dim(cnts_compcodeR)
dim(powsim_cnts)
dim(splat_cnts)
dim(seqgendiff_cnts)
dim(simseq_cnts)
dim(SPsimSeq_cnts)

a1 <- as.vector(as.matrix(cnts[,1:5]))
a1b <- as.vector(as.matrix(cnts[,51:55]))
a2 <- as.vector(as.matrix(cnts_compcodeR[,1:5]))
a3 <- as.vector(as.matrix(powsim_cnts[,1:5]))
a4 <- as.vector(as.matrix(seqgendiff_cnts[,1:5]))
a5 <- as.vector(as.matrix(simseq_cnts[,1:5]))
a6 <- as.vector(as.matrix(splat_cnts[,1:5]))
a7 <- as.vector(as.matrix(SPsimSeq_cnts[,1:5]))

p1 <- as.data.frame(qqplot(a1, a2, plot.it=FALSE))
p2 <- as.data.frame(qqplot(a1b, a3, plot.it=FALSE)) 
p3 <- as.data.frame(qqplot(a1b, a4, plot.it=FALSE))
p4 <- as.data.frame(qqplot(a1, a5, plot.it=FALSE))
p5 <- as.data.frame(qqplot(a1, a6, plot.it=FALSE))
p6 <- as.data.frame(qqplot(a1, a7, plot.it=FALSE))

a <- ggplot(p1) + geom_point(aes(x, y)) + 
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlab("NGSSPPG2 5 samp") + ylab("compcodeR 5 samp") + 
  geom_smooth(aes(x, y), method='lm', se=FALSE, formula= y~x, color = "blue") +
  theme(axis.title.x = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold")) +
#  xlim(0, 1.6e8) + ylim(0, 3e8) #H
#  xlim(0, 2500) + ylim(0, 3100) #1
  xlim(0, 7900) + ylim(0, 8100) #2

b <- ggplot(p2) + geom_point(aes(x, y)) + 
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlab("NGSSPPG2 5 samp") + ylab("powsimR 5 samp") +
  geom_smooth(aes(x, y), method='lm', se=FALSE, formula= y~x, color = "blue") +
  theme(axis.title.x = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold")) +
#  xlim(0, 1.6e8) + ylim(0, 3e8) #H
#  xlim(0, 2500) + ylim(0, 3100) #1
  xlim(0, 7900) + ylim(0, 8100) #2

c <- ggplot(p3) + geom_point(aes(x, y)) + 
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlab("NGSSPPG2 5 samp") + ylab("seqgendiff 5 samp") +
  geom_smooth(aes(x, y), method='lm', se=FALSE, formula= y~x, color = "blue") +
  theme(axis.title.x = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold")) +
#  xlim(0, 1.6e8) + ylim(0, 3e8) #H
#  xlim(0, 2500) + ylim(0, 3100) #1
  xlim(0, 7900) + ylim(0, 8100) #2

d <- ggplot(p4) + geom_point(aes(x, y)) + 
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlab("NGSSPPG2 5 samp") + ylab("SimSeq 5 samp") +
  geom_smooth(aes(x, y), method='lm', se=FALSE, formula= y~x, color = "blue") +
  theme(axis.title.x = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold")) +
#  xlim(0, 1.6e8) + ylim(0, 3e8) #H
#  xlim(0, 2500) + ylim(0, 3100) #1
  xlim(0, 7900) + ylim(0, 8100) #2

e <- ggplot(p5) + geom_point(aes(x, y)) + 
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlab("NGSSPPG2 5 samp") + ylab("Splat 5 samp") +
  geom_smooth(aes(x, y), method='lm', se=FALSE, formula= y~x, color = "blue") +
  theme(axis.title.x = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold")) +
#  xlim(0, 1.6e8) + ylim(0, 3e8) #H
# xlim(0, 2500) + ylim(0, 3100) #1
  xlim(0, 7900) + ylim(0, 8100) #2

f <- ggplot(p6) + geom_point(aes(x, y)) + 
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlab("NGSSPPG2 5 samp") + ylab("SPsimSeq 5 samp") +
  geom_smooth(aes(x, y), method='lm', se=FALSE, formula= y~x, color = "blue") +
  theme(axis.title.x = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold")) +
#  xlim(0, 1.6e8) + ylim(0, 3e8) #H
#  xlim(0, 2500) + ylim(0, 3100) #1
  xlim(0, 7900) + ylim(0, 8100) #2

theme_set(theme_pubr())
figure <- ggarrange(a, b, e, f, d, c, ncol = 3, nrow = 2,
                    labels = c("A", "B", "C", "D", "E", "F"))
ggexport(figure, filename = "QQ_all2b.png",
         res =300,
         width = 3000, height = 2000)

