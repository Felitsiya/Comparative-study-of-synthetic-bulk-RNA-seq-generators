
# Get TAMU data
cnts <- read.table("mcmcparams2_trn_mod.txt", skip = 1, sep = ",", as.is = TRUE)
cnts <- t(cnts)

# Huntley
cnts <- read.table("Huntley_counts.txt", skip = 1, sep = "\t")

# Modification for powsimR data generation
cnts[cnts > 10^7] <- 10^7
