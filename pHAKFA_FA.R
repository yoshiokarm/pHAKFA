#  pHAKFA DATA ANALYSES --------------------------------------------------------
#  code by Reyn Yoshioka

# Data analyses accompanying the article Yoshioka and Galloway et al. 2024,
# Ecosphere, Benthic marine invertebrate herbivores diversify their algal
# diets in winter (DOI: )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

# 1: Getting Started -----------------------------------------------------------
# > 1.1: Clear everything ----- 
# clear workspace
# rm(list=ls()) # commented out so folks don't clear their environments!
# clear console
cat("\014")

# > 1.2: Load packages -----
library(here)
library(vegan)
library(MASS)
library(ggplot2)
library(grid)
library(scales)
library(ggrepel)
library(gridExtra)
library(purrr)
library(viridis)
library(DHARMa) # glm diagnostics
library(ggpubr) # ggarrange plots
library(easyCODA)
library(dplyr)
library(tidyr)
library(patchwork)

# > 1.2: here::i_am()
# point to current script's location
here::i_am("pHAKFA_FA.R")

# > 1.3: Read in data -----
df_FA = read.csv("pHAKFA_FA_conc_pub.csv")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

# 2: Logratio Selection --------------------------------------------------------

# > 2.i: Data Prep -----
# Drop unreliable macro sample
df_FA = df_FA[df_FA$sampleID != "126_pHAKFA_MACRO.qgd", ]

# > 2.ii: Aggregate FA -----
# Odd and branched chain FA
df_FA$OBCFA = 
  df_FA$i15_0 +
  df_FA$a15_0 +
  df_FA$C15_0 +
  df_FA$i16_0 +
  df_FA$i17_0 +
  df_FA$a17_0 +
  df_FA$C17_0 +
  df_FA$C17_1 +
  df_FA$i18_0 +
  df_FA$i19_0 +
  df_FA$C19_1 

df_FA$C14_1nX = 
  df_FA$C14_1nX_1 +
  df_FA$C14_1nX_2

df_FA$C16_1nX =
  df_FA$C16_1nX_1 +
  df_FA$C16_1nX_2 +
  df_FA$C16_1nX_3 +
  df_FA$C16_1nX_4

df_FA$C16_2nX = 
  df_FA$C16_2nX_1 +
  df_FA$C16_2nX_2

df_FA$C18_1nX = 
  df_FA$C18_1nX_1 +
  df_FA$C18_1nX_2

df_FA$C18_2nX = 
  df_FA$C18_2

df_FA$C20_1nX =
  df_FA$C20_1nX_1 +
  df_FA$C20_1nX_2 +
  df_FA$C20_1nX_3
  
df_FA$C20_2nX = 
  df_FA$C20_2nX_1 +
  df_FA$C20_2nX_2 +
  df_FA$C20_2nX_3

df_FA$C22_2nX = 
  df_FA$C22_2nX

# Remove columns
# Remove C19.0 standard, above OBCFA, and subsumed FAs
df_FA = subset(df_FA, select = -c(C19_0,
                                  i15_0,
                                  a15_0,
                                  C15_0,
                                  i16_0,
                                  i17_0,
                                  a17_0,
                                  C17_0, 
                                  C17_1, 
                                  i18_0,
                                  i19_0,
                                  C19_1,
                                  C14_1nX_1,
                                  C14_1nX_2,
                                  C16_1nX_1,
                                  C16_1nX_2,
                                  C16_1nX_3,
                                  C16_1nX_4,
                                  C16_2nX_1,
                                  C16_2nX_2,
                                  C18_1nX_1,
                                  C18_1nX_2,
                                  C18_2,
                                  C20_1nX_1,
                                  C20_1nX_2,
                                  C20_1nX_3,
                                  C20_2nX_1,
                                  C20_2nX_2,
                                  C20_2nX_3))

# > 2.iii: Matrix prep and proportions -----
ma_FA = as.matrix(df_FA[ , 12:length(df_FA)])

class(ma_FA) = "numeric" 

# Replace NA with 0
ma_FA[is.na(ma_FA)] = 0

# replace <0 with 0
ma_FA[ma_FA < 0] = 0

# Calculate percentages by dividing by row sums
ma_FA = ma_FA / rowSums(ma_FA) * 100

# Attach it to metadata
df_FA_perc = cbind(df_FA[ , 1:11],
                   data.frame(ma_FA))

# Reshape data to long format
df_FA_percL = pivot_longer(df_FA_perc,
                           C12_0:length(df_FA_perc),
                           names_to = "FA",
                           values_to = "percent")

# Summarize FA percentages
df_FA_percM = df_FA_percL %>%
  group_by(season, spCode, FA) %>%
  summarise(means = mean(percent),
            SD = sd(percent))

# Filter FA with at least 0.5%  and 5% by season x species combination
df_FA_percM0.5 = df_FA_percM[df_FA_percM$means >= .5,]
df_FA_percM5 = df_FA_percM[df_FA_percM$means >= 5,]

# Create list of those FA
FAlist0.5 = unique(df_FA_percM0.5$FA)
FAlist5 = unique(df_FA_percM5$FA)

# Subset df_FA_percL by that list
df_FA_percL0.5 = df_FA_percL[df_FA_percL$FA %in% FAlist0.5,]

# Reshape data to wide format
df_FA_perc0.5 = pivot_wider(df_FA_percL0.5,
                            names_from = FA,
                            values_from = percent)

# Replace zeroes
zerorep = min(df_FA_percL[df_FA_percL$percent > 0, ]$percent) / 2
# Note that Graeve & Greenacre do this FA-wise, rather than choose a global min
# I like a global min since a non-zero value could still be very large if 
# a given species has really no material amounts of that FA. However, that 
# may inflate/deflate the apparent log-ratios by using a really small number,
# which is maybe why G&G do it FA-wise... Think about this...
# ma_FA[ma_FA == 0] = zerorep 
# OR
ma_FA_0NA = ma_FA
ma_FA_0NA[ma_FA_0NA == 0] = NA
FAmin = apply(ma_FA_0NA, 2, min, na.rm = T)

for(j in 1:ncol(ma_FA)) {
  for(i in 1:nrow(ma_FA)) {
    if(ma_FA[i, j] == 0) ma_FA[i, j] <- 0.5 * FAmin[j]
  }
}

ma_FA <- ma_FA / rowSums(ma_FA) * 100

# > 2.iv: Subset by type -----
ma_FA_prod = ma_FA[df_FA$type == "Algae",]
ma_FA_prod_head = df_FA[df_FA$type == "Algae", 1:11]
ma_FA_cons = ma_FA[df_FA$type == "Invert",]
ma_FA_cons_head = df_FA[df_FA$type == "Invert", 1:11]

unique(df_FA$spCode)

# > 2.v: Name keys -----
df_divphy = data.frame(spCode = c("IDOTEA",
                                   "MESOCENT",
                                   "STRONGYD",
                                   "TEGULA",
                                   "HALKAM",
                                   "PUGETTIA",
                                   "AGFIM",
                                   "MACRO",
                                   "OPUNTIELLA",
                                   "CRYPTO",
                                   "OSMUNDIA",
                                   "SACCH"),
                       divphy = c("Arthropoda",
                                   "Echinodermata",
                                   "Echinodermata",
                                   "Mollusca",
                                   "Mollusca",
                                   "Arthropoda",
                                   "Ochrophyta",
                                   "Ochrophyta",
                                   "Rhodophyta",
                                   "Rhodophyta",
                                   "Rhodophyta",
                                   "Ochrophyta"))

spCode_add = data.frame(spCode = c("IDOTEA",
                                     "MESOCENT",
                                     "STRONGYD",
                                     "TEGULA",
                                     "HALKAM",
                                     "PUGETTIA",
                                     "AGFIM",
                                     "MACRO",
                                     "OPUNTIELLA",
                                     "CRYPTO",
                                     "OSMUNDIA",
                                     "SACCH"),
                         binom = c("P. resecata",
                                   "M. franciscanus",
                                   "S. droebachiensis",
                                   "T. pulligo",
                                   "H. kamtchatkana",
                                   "P. producta",
                                   "N. fimbriatum",
                                   "M. pyrifera",
                                   "O. californica",
                                   "C. ruprechtiana",
                                   "O. spectabilis",
                                   "H. nigripes"),
                         symb = c("\u260B",
                                  "\u26ED",
                                  "\u26EF",
                                  "\u265F",
                                  "\u2659",
                                  "\u260A",
                                  "\u2665",
                                  "\u2663",
                                  "\u2667",
                                  "\u2664",
                                  "\u2661",
                                  "\u2660"))

ma_FA_prod_head = left_join(ma_FA_prod_head,
                            df_divphy,
                            by = join_by(spCode))
ma_FA_cons_head = left_join(ma_FA_cons_head,
                            df_divphy,
                            by = join_by(spCode))

# > 2.vi: G&G2020's output function (verbatim): -----
print.ratios <- function(rationames, R2, procr=NA, N=10) {
  # function prints ratios and the corresponding R2, optionally Procrustes correlations
  # prints first 10 ratios by default  
  # split names of parts in ratio into FA1 and FA2
  # notice that the ratios can be reported as FA1/FA2 or FA2/FA1  
  foo    <- strsplit(rationames, "/")
  parts  <- matrix(unlist(foo), ncol=2, byrow=TRUE)
  df   <- as.data.frame(parts)[1:N,]
  if(is.na(procr)) {
    df <- cbind(df, round(100*R2[1:N], 2))
    colnames(df) <- c("FA1", "FA2","R2")
  }
  if(!is.na(procr)) {
    df <- cbind(df, round(100*R2[1:N], 2), round(procr[1:N], 3))
    colnames(df) <- c("FA1", "FA2", "R2","Procrustes")
  }
  print(df[1:N,])
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

# > 2.1 Algae -----
# > > 2.1.1: Logratio selection ----
# Step 1:
FA_prod_step1 = STEP(data = ma_FA_prod, nsteps = 1, top = 10)
print.ratios(FA_prod_step1$names.top, FA_prod_step1$R2.top)
# FA1     FA2    R2
# 1    C16_0 C18_4n3 46.73 <- palmitic is abundant, stearidonic for browns
# 2  C18_4n3 C24_1n9 45.41
# 3  C18_4n3 C22_4n3 45.41
# 4  C18_4n3 C22_4n6 45.41
# 5  C18_4n3   C22_0 45.41
# 6  C18_4n3 C22_2nX 45.41
# 7  C18_4n3 C22_5n3 45.41
# 8  C18_4n3 C22_3n6 45.41
# 9  C18_4n3 C22_6n3 45.31
# 10 C18_4n3   OBCFA 44.83
FA_prod_LR1 = FA_prod_step1$logratios.top[,1]
FA_prod_NR1 = FA_prod_step1$ratios.top[1,]

# Step 2
FA_prod_step2 = STEP(data = ma_FA_prod, nsteps = 1, top = 10, previous = FA_prod_LR1)
print.ratios(FA_prod_step2$names.top, FA_prod_step2$R2.top)
# FA1     FA2    R2
# 1  C18_2n6c C20_5n3 63.12 <- LIN is good for greens (not here), EPA for reds/diatoms
# 2   C18_1n9 C20_4n6 62.98
# 3   C20_5n3 C20_1nX 62.82
# 4  C18_2n6c C20_4n6 62.69
# 5   C18_1n9 C20_5n3 62.60
# 6     C18_0 C20_5n3 62.40
# 7   C18_3n3 C20_5n3 62.24
# 8  C18_2n6c C20_3n6 62.18
# 9  C18_2n6t C20_5n3 62.04
# 10    C20_0 C20_5n3 61.79
FA_prod_LR2 = cbind(FA_prod_LR1, FA_prod_step2$logratios.top[,1])
FA_prod_NR2 = rbind(FA_prod_NR1, FA_prod_step2$ratios.top[1,])

# Step 3
FA_prod_step3 = STEP(data = ma_FA_prod, nsteps = 1, top = 10, previous = FA_prod_LR2)
print.ratios(FA_prod_step3$names.top, FA_prod_step3$R2.top)
# FA1     FA2    R2
# 1  C20_4n6 C22_3n6 75.72
# 2  C20_4n6 C24_1n9 75.72
# 3  C20_4n6   C22_0 75.72
# 4  C20_4n6 C22_2nX 75.72
# 5  C20_4n6 C22_4n3 75.72
# 6  C20_4n6 C22_4n6 75.72
# 7  C20_4n6 C22_5n3 75.72
# 8  C20_4n6 C22_6n3 75.61 <- ARA and DHA for browns and dinos(??), resp.
# 9    C16_0 C20_4n6 75.42
# 10 C18_4n3 C20_4n6 75.42
FA_prod_LR3 = cbind(FA_prod_LR2, FA_prod_step3$logratios.top[,8])
FA_prod_NR3 = rbind(FA_prod_NR2, FA_prod_step3$ratios.top[8,])

# Step 4
FA_prod_step4 = STEP(data = ma_FA_prod, nsteps = 1, top = 10, previous = FA_prod_LR3)
print.ratios(FA_prod_step4$names.top, FA_prod_step4$R2.top)
# FA1      FA2    R2
# 1   C18_1n9  C20_5n3 83.87 <- detritus vs diatoms/reds
# 2   C18_1n9 C18_2n6c 83.87 
# 3   C16_1n7  C18_1n9 83.51
# 4   C16_1n7  C20_1n9 83.08
# 5  C18_2n6c  C22_3n6 83.03
# 6  C18_2n6c  C24_1n9 83.03
# 7   C20_5n3  C24_1n9 83.03
# 8   C20_5n3    C22_0 83.03
# 9   C20_5n3  C22_4n3 83.03
# 10  C20_5n3  C22_3n6 83.03
FA_prod_LR4 = cbind(FA_prod_LR3, FA_prod_step4$logratios.top[,1])
FA_prod_NR4 = rbind(FA_prod_NR3, FA_prod_step4$ratios.top[1,])

# Step 5
FA_prod_step5 = STEP(data = ma_FA_prod, nsteps = 1, top = 10, previous = FA_prod_LR4)
print.ratios(FA_prod_step5$names.top, FA_prod_step5$R2.top)
# FA1     FA2    R2
# 1  C16_3n4 C18_1nX 88.90
# 2  C14_1nX C18_1nX 88.90
# 3  C18_1nX C18_2nX 88.88
# 4  C20_1n9 C18_1nX 88.84
# 5  C20_2n6 C18_1nX 88.78
# 6  C20_4n6 C18_1nX 88.76 <- ARA and unIDd 18 MUFAs; not great but ok
# 7  C22_6n3 C18_1nX 88.76
# 8  C18_1nX C20_2nX 88.74
# 9    C14_0 C18_1nX 88.73
# 10 C20_5n3 C18_1nX 88.73
FA_prod_LR5 = cbind(FA_prod_LR4, FA_prod_step5$logratios.top[,6])
FA_prod_NR5 = rbind(FA_prod_NR4, FA_prod_step5$ratios.top[6,])

# Step 6
FA_prod_step6 = STEP(data = ma_FA_prod, nsteps = 1, top = 10, previous = FA_prod_LR5)
print.ratios(FA_prod_step6$names.top, FA_prod_step6$R2.top)
# FA1      FA2    R2
# 1  C18_2n6c  C18_4n3 92.35
# 2   C18_4n3  C20_5n3 92.35 <- browns vs. reds
# 3     C16_0 C18_2n6c 92.35
# 4     C16_0  C20_5n3 92.35
# 5     C16_0  C18_1n9 92.35
# 6   C18_1n9  C18_4n3 92.35
# 7     C14_0  C20_5n3 92.31
# 8     C14_0  C18_1n9 92.31
# 9     C14_0 C18_2n6c 92.31
# 10  C18_4n3  C22_4n6 92.08
FA_prod_LR6 = cbind(FA_prod_LR5, FA_prod_step6$logratios.top[,2])
FA_prod_NR6 = rbind(FA_prod_NR5, FA_prod_step6$ratios.top[2,])

# Step 7
FA_prod_step7 = STEP(data = ma_FA_prod, nsteps = 1, top = 10, previous = FA_prod_LR6)
print.ratios(FA_prod_step7$names.top, FA_prod_step7$R2.top)
# FA1     FA2    R2
# 1  C16_1n7 C22_2nX 94.77
# 2  C16_1n7 C22_4n6 94.77
# 3  C16_1n7 C24_1n9 94.77
# 4  C16_1n7 C22_5n3 94.77
# 5  C16_1n7   C22_0 94.77
# 6  C16_1n7 C22_4n3 94.77
# 7  C16_1n7 C22_3n6 94.77
# 8  C16_1n7 C20_4n6 94.75 <- diatoms vs. browns
# 9  C16_1n7 C22_6n3 94.75
# 10 C16_1n7 C18_1nX 94.75
FA_prod_LR7 = cbind(FA_prod_LR6, FA_prod_step7$logratios.top[,8])
FA_prod_NR7 = rbind(FA_prod_NR6, FA_prod_step7$ratios.top[8,])

# Step 8
FA_prod_step8 = STEP(data = ma_FA_prod, nsteps = 1, top = 10, previous = FA_prod_LR7)
print.ratios(FA_prod_step8$names.top, FA_prod_step8$R2.top)
# FA1      FA2    R2
# 1  C14_0    C18_0 96.47
# 2  C18_0 C18_2n6c 96.40 <- stearic common, LIN for plants
# 3  C18_0  C18_1n9 96.40
# 4  C18_0  C20_5n3 96.40
# 5  C18_0  C18_4n3 96.40
# 6  C16_0    C18_0 96.40
# 7  C18_0  C20_1n9 96.40
# 8  C14_0  C20_1nX 96.30
# 9  C18_0  C24_1n9 96.25
# 10 C18_0  C22_4n3 96.25
FA_prod_LR8 = cbind(FA_prod_LR7, FA_prod_step8$logratios.top[,2])
FA_prod_NR8 = rbind(FA_prod_NR7, FA_prod_step8$ratios.top[2,])

# Step 9
FA_prod_step9 = STEP(data = ma_FA_prod, nsteps = 1, top = 10, previous = FA_prod_LR8)
print.ratios(FA_prod_step9$names.top, FA_prod_step9$R2.top)
# FA1     FA2    R2
# 1   C16_1nX C20_1nX 97.52
# 2  C18_2n6t C16_1nX 97.52
# 3   C18_3n3 C16_1nX 97.51
# 4   C16_1n7 C16_1nX 97.49
# 5   C20_4n6 C16_1nX 97.49
# 6   C16_1nX C18_1nX 97.49
# 7   C22_6n3 C16_1nX 97.49
# 8   C22_4n6 C16_1nX 97.49
# 9   C24_1n9 C16_1nX 97.49
# 10  C22_5n3 C16_1nX 97.49

# Stop here. Diminishing returns and these MUFAs have limited interpretability 
# Captured 96.40% of the variation!

rownames(FA_prod_NR8) = paste("Step", 1:8, sep="")
colnames(FA_prod_NR8) = c("FA1","FA2")
FA_prod_FR = as.data.frame(cbind(FA_prod_NR8,
                                 Ratio = paste(colnames(ma_FA_prod)[FA_prod_NR8[,1]],
                                               "/",
                                               colnames(ma_FA_prod)[FA_prod_NR8[,2]],
                                               sep="")))
FA_prod_FR

colnames(FA_prod_LR8) = FA_prod_FR[,3]

FA_prod_PiR = sort(unique(as.numeric(FA_prod_NR8))) # "parts in ratios"
colnames(ma_FA_prod)[FA_prod_PiR]

FA_prod_PCA = PCA(FA_prod_LR8, weight = FALSE)

PLOT.PCA(FA_prod_PCA,
         map = "contribution",
         axes.inv = c(1, 1),
         rescale = 2)

# Extract and calculate contribution vectors
axes_inv = c(1, 1)
dim = c(1, 2)
FA_prod_RPC = FA_prod_PCA$rowcoord[, dim] %*% diag(FA_prod_PCA$sv[dim] * axes_inv)
FA_prod_CSC = FA_prod_PCA$colcoord[, dim] %*% diag(axes_inv)
FA_prod_CPC = FA_prod_CSC %*% diag(FA_prod_PCA$sv[dim])
FA_prod_CCC = FA_prod_CSC * sqrt(FA_prod_PCA$colmass)
FA_prod_CRD = FA_prod_CCC

# rescale contribution vectors for visibility
vectorscale = 2
FA_prod_CRD_scale = FA_prod_CRD * vectorscale

# Extract points
df_FA_prod_RPC = cbind(ma_FA_prod_head, FA_prod_RPC)
colnames(df_FA_prod_RPC)=c(colnames(ma_FA_prod_head),'dim1','dim2')

# Aggregate by season and sp
df_FA_prod_RPC_MSE = df_FA_prod_RPC %>%
  group_by(season, spCode) %>%
  summarize(dim1M = mean(dim1),
            dim2M = mean(dim2),
            dim1SE = sd(dim1) / sqrt(length(dim1)),
            dim2SE = sd(dim2) / sqrt(length(dim2)))

# Calculate variance explained by PCA dimensions
FA_prod_PCA_perc_1 = 100 * FA_prod_PCA$sv[dim[1]]^2 / sum(FA_prod_PCA$sv^2)
FA_prod_PCA_perc_2 = 100 * FA_prod_PCA$sv[dim[2]]^2 / sum(FA_prod_PCA$sv^2)

df_FA_prod_RPC$ID2 = paste(df_FA_prod_RPC$spCode,
                           rownames(df_FA_prod_RPC),
                           sep = "_")

# Format arrows
df_FA_prod_arrows = data.frame(dim1 = FA_prod_CRD_scale[ ,1],
                               dim2 = FA_prod_CRD_scale[, 2],
                               LR = FA_prod_PCA$colnames)

df_FA_prod_arrows$angle = 90 - ((atan2(df_FA_prod_arrows$dim1,
                                       df_FA_prod_arrows$dim2) *
                                   180) /
                                  pi)
df_FA_prod_arrows$angle = ifelse(df_FA_prod_arrows$angle > 90 & df_FA_prod_arrows$angle < 270,
                                 df_FA_prod_arrows$angle + 180,
                                 df_FA_prod_arrows$angle)
df_FA_prod_arrows$hjust = ifelse(df_FA_prod_arrows$angle > 90,
                                 1,
                                 0)

df_FA_prod_arrows$LR = gsub("C", "", df_FA_prod_arrows$LR)
df_FA_prod_arrows$LR = sub("c", "", df_FA_prod_arrows$LR)
df_FA_prod_arrows$LR = gsub("_", ":", df_FA_prod_arrows$LR)
df_FA_prod_arrows$LR = gsub("n", "*omega*-", df_FA_prod_arrows$LR)

df_FA_prod_RPC = left_join(df_FA_prod_RPC,
                           spCode_add[,c("spCode","binom")],
                           by = join_by(spCode))
df_FA_prod_RPC_MSE = left_join(df_FA_prod_RPC_MSE,
                               spCode_add[,c("spCode","binom")],
                               by = join_by(spCode))

# > > 2.1.2: Plot PCA ----
FA_prod_PCA_plot = ggplot() +
  geom_point(aes(x = dim1,
                 y = dim2,
                 color = season,
                 shape = binom),
             stroke = 0.5,
             alpha = 0.5,
             data = df_FA_prod_RPC) +
  # geom_text(aes(x = dim1,
  #               y = dim2,
  #               color = season,
  #               label = ID2),
  #           size = 2,
  #           data = df_FA_prod_RPC) +
  geom_segment(aes(x = 0,
                   y = 0,
                   xend = FA_prod_CRD_scale[ ,1],
                   yend = FA_prod_CRD_scale[ ,2]),
               color = "grey50",
               linewidth = 0.5,
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_errorbar(aes(x = dim1M,
                    ymin = dim2M - 1.96 * dim2SE,
                    ymax = dim2M + 1.96 * dim2SE,
                    color = season),
                width = 0,
                linewidth = 0.5,
                data = df_FA_prod_RPC_MSE) +
  geom_errorbarh(aes(y = dim2M,
                     xmin = dim1M - 1.96 * dim1SE,
                     xmax = dim1M + 1.96 * dim1SE,
                     color = season),
                 height = 0,
                 linewidth = 0.5,
                 data = df_FA_prod_RPC_MSE) +
  geom_text(aes(x = dim1 * 1.05,
                y = dim2 * 1.05,
                label = LR),
            colour = "grey50",family = "sans",
            size = 2,
            hjust = df_FA_prod_arrows$hjust,
            angle = df_FA_prod_arrows$angle,
            parse = TRUE,
            data = df_FA_prod_arrows) +
  geom_point(aes(x = dim1M,
                 y = dim2M,
                 color = season,
                 shape = binom),
             size = 3,
             stroke = 1,
             data = df_FA_prod_RPC_MSE) +
  scale_shape_manual(values = c(0, 1, 2, 4, 5, 6)) +
  scale_colour_manual(values = c("#9B1F1A", "#5FA1F7"), 
                      labels = c("Summer", "Winter")) + 
  scale_y_continuous(sec.axis = ~./vectorscale) +
  scale_x_continuous(sec.axis = ~./vectorscale) +
  labs(x = paste("PC 1 (",
                 round(FA_prod_PCA_perc_1,1),
                 "%)",
                 sep=""),
       y = paste("PC 2 (",
                 round(FA_prod_PCA_perc_2,1),
                 "%)",
                 sep=""),
       color = "Season",
       shape = "Species") +
  coord_fixed() +
  theme_classic() +
  guides(shape = guide_legend(label.theme = element_text(face = "italic",
                                                         size = unit(10, "pt"))),
         color = guide_legend(label.theme = element_text(size = unit(10, "pt")))) +
  theme(axis.text.y.right = element_text(colour="grey50"),
        axis.text.x.top = element_text(colour="grey50"),
        axis.ticks.y.right = element_line(colour="grey50"),
        axis.ticks.x.top = element_line(colour="grey50"),
        axis.line.y.right = element_line(colour="grey50"),
        axis.line.x.top = element_line(colour="grey50"))

FA_prod_PCA_plot

ggsave("pHAKFA_Fig_01.pdf",
       device = "pdf",
       height = 12,
       width = 18,
       units = "cm")

# > > 2.1.3: Univariate plot -----
FA_prod_LR8_long = cbind(ma_FA_prod_head, FA_prod_LR8)
FA_prod_LR8_long = pivot_longer(FA_prod_LR8_long,
                                cols = (length(FA_prod_LR8_long)-7):length(FA_prod_LR8_long),
                                names_to = "FAs",
                                values_to = "logratio")

FA_prod_LR8_long = left_join(FA_prod_LR8_long,
                             spCode_add[,c("spCode","binom")],
                             by = join_by(spCode))

FA_prod_LR8_long$binom = factor(FA_prod_LR8_long$binom,
                                  levels = c("N. fimbriatum",
                                             "M. pyrifera",
                                             "H. nigripes",
                                             "C. ruprechtiana",
                                             "O. californica",
                                             "O. spectabilis"))

FA_prod_LR8_MSE = FA_prod_LR8_long %>%
  group_by(season, binom, FAs) %>%
  summarize(M = mean(logratio),
            SE = sd(logratio)/sqrt(length(logratio)))

FA_prod_LR8_MSE2 = pivot_wider(FA_prod_LR8_MSE,
                               id_cols = c(FAs, binom),
                               names_from = season,
                               values_from = M)
FA_prod_LR8_MSE2$diff = FA_prod_LR8_MSE2$summer - FA_prod_LR8_MSE2$winter
FA_prod_LR8_MSE2$diff_abs = abs(FA_prod_LR8_MSE2$diff)

FA_prod_LR8_MSE2 = FA_prod_LR8_MSE2 %>%
  group_by(FAs) %>%
  summarize(diff_abs = mean(diff_abs))

FA_prod_LR8_long$FAs = gsub("C", "", FA_prod_LR8_long$FAs)
FA_prod_LR8_long$FAs = sub("c", "", FA_prod_LR8_long$FAs)
FA_prod_LR8_long$FAs = gsub("_", ":", FA_prod_LR8_long$FAs)
FA_prod_LR8_long$FAs = gsub("w", "n-", FA_prod_LR8_long$FAs)

FA_prod_LR8_MSE$FAs = gsub("C", "", FA_prod_LR8_MSE$FAs)
FA_prod_LR8_MSE$FAs = sub("c", "", FA_prod_LR8_MSE$FAs)
FA_prod_LR8_MSE$FAs = gsub("_", ":", FA_prod_LR8_MSE$FAs)
FA_prod_LR8_MSE$FAs = gsub("w", "n-", FA_prod_LR8_MSE$FAs)

FA_prod_LR8_MSE2$FAs = gsub("C", "", FA_prod_LR8_MSE2$FAs)
FA_prod_LR8_MSE2$FAs = sub("c", "", FA_prod_LR8_MSE2$FAs)
FA_prod_LR8_MSE2$FAs = gsub("_", ":", FA_prod_LR8_MSE2$FAs)
FA_prod_LR8_MSE2$FAs = gsub("w", "n-", FA_prod_LR8_MSE2$FAs)

FA_prod_LR8_long$FAs = factor(FA_prod_LR8_long$FAs,
                              levels = FA_prod_LR8_MSE2$FAs[order(FA_prod_LR8_MSE2$diff_abs,
                                                                  decreasing = TRUE)])
FA_prod_LR8_MSE$FAs = factor(FA_prod_LR8_MSE$FAs,
                              levels = FA_prod_LR8_MSE2$FAs[order(FA_prod_LR8_MSE2$diff_abs,
                                                                  decreasing = TRUE)])

FA_prod_univ_plot = ggplot(aes(x = binom,
                               y = logratio,
                               color = season),
                           data = FA_prod_LR8_long) +
  geom_point(size = 1,
             alpha = 0.2,
             position = position_dodge(.9)) +
  geom_point(aes(x = binom,
                 y = M,
                 color = season),
             size = 1,
             position = position_dodge(.9),
             data = FA_prod_LR8_MSE) +
  geom_errorbar(aes(x = binom,
                    y = M,
                    ymin = M - SE * 1.96,
                    ymax = M + SE * 1.96,
                    color = season),
                position = position_dodge(0.9),
                data = FA_prod_LR8_MSE) +
  scale_colour_manual(values = c("#9B1F1A", "#5FA1F7"), 
                      labels = c("Summer", "Winter")) + 
  labs(x = "Species",
       y = "Logratio",
       color = "Season") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   size = unit(6, "pt"),
                                   face = "italic"),
        axis.text.y = element_text(size = unit(6, "pt")),
        axis.title = element_text(size = unit(8, "pt")),
        strip.text = element_text(size = unit(5, "pt")),
        legend.title = element_text(size = unit(8, "pt")),
        legend.text = element_text(size = unit(6, "pt")),
        legend.position = "bottom") +
  facet_grid(. ~ FAs)
FA_prod_univ_plot

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

# > 2.2: Herbivores w/ Algal logratios -----
# > > 2.2.1: Use same logatios (no selection) ----
FA_cons_LR8P = data.frame(ma_FA_cons)
FA_cons_LR8P$"C16_0/C18_4n3" = log(FA_cons_LR8P$C16_0 /
                                     FA_cons_LR8P$C18_4n3)
FA_cons_LR8P$"C18_2n6c/C20_5n3" = log(FA_cons_LR8P$C18_2n6c /
                                        FA_cons_LR8P$C20_5n3)
FA_cons_LR8P$"C20_4n6/C22_6n3" = log(FA_cons_LR8P$C20_4n6 /
                                       FA_cons_LR8P$C22_6n3)
FA_cons_LR8P$"C18_1n9/C20_5n3" = log(FA_cons_LR8P$C18_1n9 /
                                       FA_cons_LR8P$C20_5n3)
FA_cons_LR8P$"C20_4n6/C18_1nX" = log(FA_cons_LR8P$C20_4n6 /
                                       FA_cons_LR8P$C18_1nX)
FA_cons_LR8P$"C18_4n3/C20_5n3" = log(FA_cons_LR8P$C18_4n3 /
                                       FA_cons_LR8P$C20_5n3)
FA_cons_LR8P$"C16_1n7/C20_4n6" = log(FA_cons_LR8P$C16_1n7 /
                                       FA_cons_LR8P$C20_4n6)
FA_cons_LR8P$"C18_0/C18_2n6c" = log(FA_cons_LR8P$C18_0 /
                                      FA_cons_LR8P$C18_2n6c)

FA_cons_LR8P = FA_cons_LR8P[ ,(length(FA_cons_LR8P) - 7):length(FA_cons_LR8P)]

FA_cons_LR8P = as.matrix(FA_cons_LR8P)

FA_cons_P_PCA = PCA(FA_cons_LR8P, weight = FALSE)

PLOT.PCA(FA_cons_P_PCA,
         map = "contribution",
         axes.inv = c(1, 1),
         rescale = 2)

FA_cons_P_RPC = FA_cons_P_PCA$rowcoord[, dim] %*% diag(FA_cons_P_PCA$sv[dim] * axes_inv)
FA_cons_P_CSC = FA_cons_P_PCA$colcoord[, dim] %*% diag(axes_inv)
FA_cons_P_CPC = FA_cons_P_CSC %*% diag(FA_cons_P_PCA$sv[dim])
FA_cons_P_CCC = FA_cons_P_CSC * sqrt(FA_cons_P_PCA$colmass)
FA_cons_P_CRD = FA_cons_P_CCC

# rescale contribution vectors for visibility
vectorscale = 2
FA_cons_P_CRD_scale = FA_cons_P_CRD * vectorscale

# Extract points
df_FA_cons_P_RPC = cbind(ma_FA_cons_head, FA_cons_P_RPC)
colnames(df_FA_cons_P_RPC) = c(colnames(ma_FA_cons_head),
                               'dim1',
                               'dim2')

# Aggregate by season and sp
df_FA_cons_P_RPC_MSE = df_FA_cons_P_RPC %>%
  group_by(season, spCode) %>%
  summarize(dim1M = mean(dim1),
            dim2M = mean(dim2),
            dim1SE = sd(dim1) / sqrt(length(dim1)),
            dim2SE = sd(dim2) / sqrt(length(dim2)))

# Calculate variance explained by PCA dimensions
FA_cons_P_PCA_perc_1 = 100 * FA_cons_P_PCA$sv[dim[1]]^2 / sum(FA_cons_P_PCA$sv^2)
FA_cons_P_PCA_perc_2 = 100 * FA_cons_P_PCA$sv[dim[2]]^2 / sum(FA_cons_P_PCA$sv^2)

df_FA_cons_P_RPC$ID2 = paste(df_FA_cons_P_RPC$spCode,
                           rownames(df_FA_cons_P_RPC),
                           sep = "_")

# Format arrows
df_FA_cons_P_arrows = data.frame(dim1 = FA_cons_P_CRD_scale[ ,1],
                               dim2 = FA_cons_P_CRD_scale[, 2],
                               LR = FA_cons_P_PCA$colnames)

df_FA_cons_P_arrows$angle = 90 - ((atan2(df_FA_cons_P_arrows$dim1,
                                       df_FA_cons_P_arrows$dim2) *
                                   180) /
                                  pi)
df_FA_cons_P_arrows$angle = ifelse(df_FA_cons_P_arrows$angle > 90 & df_FA_cons_P_arrows$angle < 270,
                                 df_FA_cons_P_arrows$angle + 180,
                                 df_FA_cons_P_arrows$angle)
df_FA_cons_P_arrows$hjust = ifelse(df_FA_cons_P_arrows$angle > 90,
                                 1,
                                 0)

df_FA_cons_P_arrows$LR = gsub("C", "", df_FA_cons_P_arrows$LR)
df_FA_cons_P_arrows$LR = sub("c", "", df_FA_cons_P_arrows$LR)
df_FA_cons_P_arrows$LR = gsub("_", ":", df_FA_cons_P_arrows$LR)
df_FA_cons_P_arrows$LR = gsub("n", "*omega*-", df_FA_cons_P_arrows$LR)

df_FA_cons_P_RPC = left_join(df_FA_cons_P_RPC,
                           spCode_add[,c("spCode","binom")],
                           by = join_by(spCode))
df_FA_cons_P_RPC_MSE = left_join(df_FA_cons_P_RPC_MSE,
                               spCode_add[,c("spCode","binom")],
                               by = join_by(spCode))

# > > 2.2.2: Plot PCA -----
FA_cons_P_PCA_plot = ggplot() +
  geom_point(aes(x = dim1,
                 y = dim2,
                 color = season,
                 shape = binom),
             stroke = 0.5,
             alpha = 0.5,
             data = df_FA_cons_P_RPC) +
  # geom_text(aes(x = dim1,
  #               y = dim2,
  #               color = season,
  #               label = ID2),
  #           size = 2,
  #           data = df_FA_cons_P_RPC) +
  geom_segment(aes(x = 0,
                   y = 0,
                   xend = FA_cons_P_CRD_scale[ ,1],
                   yend = FA_cons_P_CRD_scale[ ,2]),
               color = "grey50",
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(aes(x = dim1 * 1.05,
                y = dim2 * 1.05,
                label = LR),
            colour = "grey50",
            size = 2,
            parse = TRUE,
            hjust = df_FA_cons_P_arrows$hjust,
            angle = df_FA_cons_P_arrows$angle,
            data = df_FA_cons_P_arrows)+
  geom_errorbar(aes(x = dim1M,
                    ymin = dim2M - 1.96 * dim2SE,
                    ymax = dim2M + 1.96 * dim2SE,
                    color = season),
                linewidth = 0.5,
                width = 0,
                data = df_FA_cons_P_RPC_MSE) +
  geom_errorbarh(aes(y = dim2M,
                     xmin = dim1M - 1.96 * dim1SE,
                     xmax = dim1M + 1.96 * dim1SE,
                     color = season),
                 linewidth = 0.5,
                 height = 0,
                 data = df_FA_cons_P_RPC_MSE) +
  geom_point(aes(x = dim1M,
                 y = dim2M,
                 color = season,
                 shape = binom),
             size = 3,
             stroke = 1,
             data = df_FA_cons_P_RPC_MSE) +
  scale_shape_manual(values = c(0, 1, 2, 4, 5, 6)) +
  scale_colour_manual(values = c("#9B1F1A", "#5FA1F7"), 
                      labels = c("Summer", "Winter")) + 
  scale_y_continuous(limits = c(-2.1, 1.1),
                     sec.axis = ~./vectorscale)+
  scale_x_continuous(sec.axis = ~./vectorscale)+
  labs(x = paste("PC 1 (",
                 round(FA_cons_P_PCA_perc_1, 1),
                 "%)",
                 sep=""),
       y = paste("PC 2 (",
                 round(FA_cons_P_PCA_perc_2, 1),
                 "%)",
                 sep=""),
       color = "Season",
       shape = "Species") +
  coord_fixed() +
  theme_classic() +
  guides(shape = guide_legend(label.theme = element_text(face = "italic",
                                                         size = unit(10, "pt"))),
         color = guide_legend(label.theme = element_text(size = unit(10, "pt")))) +
  theme(axis.text.y.right = element_text(colour="grey50"),
        axis.text.x.top = element_text(colour="grey50"),
        axis.ticks.y.right = element_line(colour="grey50"),
        axis.ticks.x.top = element_line(colour="grey50"),
        axis.line.y.right = element_line(colour="grey50"),
        axis.line.x.top = element_line(colour="grey50"))

FA_cons_P_PCA_plot

ggsave("pHAKFA_Fig_02.pdf",
       device = "pdf",
       height = 12,
       width = 18,
       units = "cm")

# > > 2.2.3: Univariate plot -----
FA_cons_LR8P_long = cbind(ma_FA_cons_head, FA_cons_LR8P)
FA_cons_LR8P_long = pivot_longer(FA_cons_LR8P_long,
                                cols = (length(FA_cons_LR8P_long)-7):length(FA_cons_LR8P_long),
                                names_to = "FAs",
                                values_to = "logratio")

FA_cons_LR8P_long = left_join(FA_cons_LR8P_long,
                             spCode_add[,c("spCode","binom")],
                             by = join_by(spCode))

FA_cons_LR8P_long$binom = factor(FA_cons_LR8P_long$binom,
                                levels = c("M. franciscanus",
                                           "S. droebachiensis",
                                           "H. kamtchatkana",
                                           "T. pulligo",
                                           "P. resecata",
                                           "P. producta"))

FA_cons_LR8P_MSE = FA_cons_LR8P_long %>%
  group_by(season, binom, FAs) %>%
  summarize(M = mean(logratio),
            SE = sd(logratio)/sqrt(length(logratio)))

FA_cons_LR8P_MSE2 = pivot_wider(FA_cons_LR8P_MSE,
                               id_cols = c(FAs, binom),
                               names_from = season,
                               values_from = M)
FA_cons_LR8P_MSE2$diff = FA_cons_LR8P_MSE2$summer - FA_cons_LR8P_MSE2$winter
FA_cons_LR8P_MSE2$diff_abs = abs(FA_cons_LR8P_MSE2$diff)

FA_cons_LR8P_MSE2 = FA_cons_LR8P_MSE2 %>%
  group_by(FAs) %>%
  summarize(diff_abs = mean(diff_abs))

FA_cons_LR8P_long$FAs = gsub("C", "", FA_cons_LR8P_long$FAs)
FA_cons_LR8P_long$FAs = sub("c", "", FA_cons_LR8P_long$FAs)
FA_cons_LR8P_long$FAs = gsub("_", ":", FA_cons_LR8P_long$FAs)
FA_cons_LR8P_long$FAs = gsub("n", "n-", FA_cons_LR8P_long$FAs)

FA_cons_LR8P_MSE$FAs = gsub("C", "", FA_cons_LR8P_MSE$FAs)
FA_cons_LR8P_MSE$FAs = sub("c", "", FA_cons_LR8P_MSE$FAs)
FA_cons_LR8P_MSE$FAs = gsub("_", ":", FA_cons_LR8P_MSE$FAs)
FA_cons_LR8P_MSE$FAs = gsub("n", "n-", FA_cons_LR8P_MSE$FAs)

FA_cons_LR8P_MSE2$FAs = gsub("C", "", FA_cons_LR8P_MSE2$FAs)
FA_cons_LR8P_MSE2$FAs = sub("c", "", FA_cons_LR8P_MSE2$FAs)
FA_cons_LR8P_MSE2$FAs = gsub("_", ":", FA_cons_LR8P_MSE2$FAs)
FA_cons_LR8P_MSE2$FAs = gsub("n", "n-", FA_cons_LR8P_MSE2$FAs)

FA_cons_LR8P_long$FAs = factor(FA_cons_LR8P_long$FAs,
                               levels = FA_cons_LR8P_MSE2$FAs[order(FA_cons_LR8P_MSE2$diff_abs,
                                                                    decreasing = TRUE)])
FA_cons_LR8P_MSE$FAs = factor(FA_cons_LR8P_MSE$FAs,
                              levels = FA_cons_LR8P_MSE2$FAs[order(FA_cons_LR8P_MSE2$diff_abs,
                                                                   decreasing = TRUE)])

FA_consP_univ_plot = ggplot(aes(x = binom,
                                y = logratio,
                                color = season),
                            data = FA_cons_LR8P_long) +
  geom_point(size = 1,
             alpha = 0.2,
             position = position_dodge(.9)) +
  geom_point(aes(x = binom,
                 y = M,
                 color = season),
             size = 1,
             position = position_dodge(.9),
             data = FA_cons_LR8P_MSE) +
  geom_errorbar(aes(x = binom,
                    y = M,
                    ymin = M - SE * 1.96,
                    ymax = M + SE * 1.96,
                    color = season),
                position = position_dodge(0.9),
                data = FA_cons_LR8P_MSE) +
  scale_colour_manual(values = c("#9B1F1A", "#5FA1F7"), 
                      labels = c("Summer", "Winter")) + 
  theme_bw() +
  labs(x = "Species",
       y = "Logratio",
       color = "Season") +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   size = unit(6, "pt"),
                                   face = "italic"),
        axis.text.y = element_text(size = unit(6, "pt")),
        axis.title = element_text(size = unit(8, "pt")),
        strip.text = element_text(size = unit(5, "pt")),
        legend.title = element_text(size = unit(8, "pt")),
        legend.text = element_text(size = unit(6, "pt")),
        legend.position = "bottom") +
  facet_grid(. ~ FAs)
FA_consP_univ_plot


# > > 2.2.4: Combined PCA with algae and herbivores using algal LRs -----
ma_FA_head = rbind(ma_FA_prod_head,
                   ma_FA_cons_head)
FA_comb_LR8 = rbind(FA_prod_LR8,
                    FA_cons_LR8P)

FA_comb_PCA = PCA(FA_comb_LR8, weight = FALSE)

PLOT.PCA(FA_comb_PCA,
         map = "contribution",
         axes.inv = c(1, 1),
         rescale = 2)

# Extract and calculate contribution vectors
axes_inv = c(1, 1)
dim = c(1, 2)
FA_comb_RPC = FA_comb_PCA$rowcoord[, dim] %*% diag(FA_comb_PCA$sv[dim] * axes_inv)
FA_comb_CSC = FA_comb_PCA$colcoord[, dim] %*% diag(axes_inv)
FA_comb_CPC = FA_comb_CSC %*% diag(FA_comb_PCA$sv[dim])
FA_comb_CCC = FA_comb_CSC * sqrt(FA_comb_PCA$colmass)
FA_comb_CRD = FA_comb_CCC

# rescale contribution vectors for visibility
vectorscale = 2
FA_comb_CRD_scale = FA_comb_CRD * vectorscale

# Extract points
df_FA_comb_RPC = cbind(ma_FA_head, FA_comb_RPC)
colnames(df_FA_comb_RPC) = c(colnames(ma_FA_head),'dim1','dim2')

# Aggregate by season and sp
df_FA_comb_RPC_MSE = df_FA_comb_RPC %>%
  group_by(season, spCode) %>%
  summarize(dim1M = mean(dim1),
            dim2M = mean(dim2),
            dim1SE = sd(dim1) / sqrt(length(dim1)),
            dim2SE = sd(dim2) / sqrt(length(dim2)))

# Calculate variance explained by PCA dimensions
FA_comb_PCA_perc_1 = 100 * FA_comb_PCA$sv[dim[1]]^2 / sum(FA_comb_PCA$sv^2)
FA_comb_PCA_perc_2 = 100 * FA_comb_PCA$sv[dim[2]]^2 / sum(FA_comb_PCA$sv^2)

df_FA_comb_RPC$ID2 = paste(df_FA_comb_RPC$spCode,
                           rownames(df_FA_comb_RPC),
                           sep = "_")

# Format arrows
df_FA_comb_arrows = data.frame(dim1 = FA_comb_CRD_scale[ ,1],
                               dim2 = FA_comb_CRD_scale[, 2],
                               LR = FA_comb_PCA$colnames)

df_FA_comb_arrows$angle = 90 - ((atan2(df_FA_comb_arrows$dim1,
                                       df_FA_comb_arrows$dim2) *
                                   180) /
                                  pi)
df_FA_comb_arrows$angle = ifelse(df_FA_comb_arrows$angle > 90 & df_FA_comb_arrows$angle < 270,
                                 df_FA_comb_arrows$angle + 180,
                                 df_FA_comb_arrows$angle)
df_FA_comb_arrows$hjust = ifelse(df_FA_comb_arrows$angle > 90,
                                 1,
                                 0)

df_FA_comb_arrows$LR = gsub("C", "", df_FA_comb_arrows$LR)
df_FA_comb_arrows$LR = sub("c", "", df_FA_comb_arrows$LR)
df_FA_comb_arrows$LR = gsub("_", ":", df_FA_comb_arrows$LR)
df_FA_comb_arrows$LR = gsub("n", "*omega*-", df_FA_comb_arrows$LR)

df_FA_comb_RPC = left_join(df_FA_comb_RPC,
                           spCode_add,
                           by = join_by(spCode))

df_FA_comb_RPC$binom = factor(df_FA_comb_RPC$binom,
                              levels = c("N. fimbriatum",
                                         "M. pyrifera",
                                         "H. nigripes",
                                         "C. ruprechtiana",
                                         "O. californica",
                                         "O. spectabilis",
                                         "M. franciscanus",
                                         "S. droebachiensis",
                                         "T. pulligo",
                                         "H. kamtchatkana",
                                         "P. resecata",
                                         "P. producta"))

# plot
FA_comb_PCA_plot = ggplot() +
  geom_point(aes(x = dim1,
                 y = dim2,
                 color = binom,
                 shape = interaction(season, type)),
             alpha = 1.0,
             size = 3,
             data = df_FA_comb_RPC) +
  # geom_text(aes(x = dim1,
  #               y = dim2,
  #               color = season,
  #               label = symb),
  #           size = 5,
  #           family = "Apple Symbols",
  #           data = df_FA_comb_RPC) +
  geom_segment(aes(x = 0,
                   y = 0,
                   xend = FA_comb_CRD_scale[ ,1],
                   yend = FA_comb_CRD_scale[ ,2]),
               color = "grey70",
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(aes(x = dim1 * 1.1,
                y = dim2 * 1.1,
                label = LR),
            colour = "grey70",
            size = 2,
            parse = TRUE,
            hjust = df_FA_comb_arrows$hjust,
            angle = df_FA_comb_arrows$angle,
            data = df_FA_comb_arrows)+
  # geom_errorbar(aes(x = dim1M,
  #                   ymin = dim2M - dim2SE,
  #                   ymax = dim2M + dim2SE,
  #                   color = season),
  #               width = 0,
  #               data = df_FA_comb_RPC_MSE) +
  # geom_errorbarh(aes(y = dim2M,
  #                    xmin = dim1M - dim1SE,
  #                    xmax = dim1M + dim1SE,
  #                    color = season),
  #                height = 0,
  #                data = df_FA_comb_RPC_MSE) +
  # geom_point(aes(x = dim1M,
  #                y = dim2M,
  #                color = season,
  #                shape = spCode),
  #            size = 3,
  #            data = df_FA_comb_RPC_MSE) +
  scale_shape_manual(values = c(2,
                                17,
                                1,
                                19),
                     labels = c("producers, summer",
                                "producers, winter",
                                "consumers, summer",
                                "consumers, winter")) +
  scale_color_viridis_d(option = "turbo") +
  scale_y_continuous(sec.axis = ~./vectorscale) +
  scale_x_continuous(sec.axis = ~./vectorscale) +
  labs(x = paste("PC 1 (",
                 round(FA_comb_PCA_perc_1,1),
                 "%)",
                 sep=""),
       y = paste("PC 2 (",
                 round(FA_comb_PCA_perc_2,1),
                 "%)",
                 sep=""),
       color = "Species",
       shape = "Trophic Level, Season") +
  coord_fixed() +
  theme_classic() +
  guides(shape = guide_legend(label.theme = element_text(size = unit(10, "pt"))),
         color = guide_legend(label.theme = element_text(face = "italic",
                                                         size = unit(10, "pt")))) +
  theme(axis.text.y.right = element_text(colour="grey70"),
        axis.text.x.top = element_text(colour="grey70"),
        axis.ticks.y.right = element_line(colour="grey70"),
        axis.ticks.x.top = element_line(colour="grey70"),
        axis.line.y.right = element_line(colour="grey70"),
        axis.line.x.top = element_line(colour="grey70"))

FA_comb_PCA_plot

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

# > 2.3: Herbivores -----
# > > 2.3.1: Logratio selection -----
# Step 1:
FA_cons_step1 = STEP(data = ma_FA_cons, nsteps = 1, top = 10)
print.ratios(FA_cons_step1$names.top, FA_cons_step1$R2.top)
# FA1     FA2    R2
# 1   C22_5n3 C20_2nX 49.08 <- not sure how to interpret the 20:2, though
# 2     C18_0 C18_1nX 47.28
# 3     OBCFA C20_2nX 47.17
# 4   C22_5n3 C16_1nX 47.04
# 5   C14_1n5 C22_5n3 46.70
# 6   C14_1n5   C18_0 46.68
# 7   C14_1n5 C22_4n3 46.42
# 8   C22_5n3 C20_1nX 46.35
# 9   C14_1n5   OBCFA 46.00
# 10 C22_1n9c C22_6n3 45.70
FA_cons_LR1 = FA_cons_step1$logratios.top[,1]
FA_cons_NR1 = FA_cons_step1$ratios.top[1,]

# Step 2
FA_cons_step2 = STEP(data = ma_FA_cons, nsteps = 1, top = 10, previous = FA_cons_LR1)
print.ratios(FA_cons_step2$names.top, FA_cons_step2$R2.top)
# FA1     FA2    R2
# 1   C20_3n3 C22_5n3 75.80
# 2   C20_3n3 C20_2nX 75.80
# 3   C20_5n3 C22_5n3 75.68 <- EPA vs DPA slightly more interpretable?
# 4   C20_5n3 C20_2nX 75.68
# 5   C18_3n6 C22_5n3 75.37
# 6   C18_3n6 C20_2nX 75.37
# 7   C20_3n3 C20_1nX 75.28
# 8   C20_3n3   OBCFA 75.09
# 9  C18_1n7c C20_3n3 75.08
# 10 C18_2n6c C20_2nX 75.02
FA_cons_LR2 = cbind(FA_cons_LR1, FA_cons_step2$logratios.top[,3])
FA_cons_NR2 = rbind(FA_cons_NR1, FA_cons_step2$ratios.top[3,])

# Step 3
FA_cons_step3 = STEP(data = ma_FA_cons, nsteps = 1, top = 10, previous = FA_cons_LR2)
print.ratios(FA_cons_step3$names.top, FA_cons_step3$R2.top)
# FA1     FA2    R2
# 1   C18_4n3 C22_3n6 82.41 <- SDA and this n-6?
# 2     C14_0 C22_3n6 82.35
# 3   C18_3n3 C22_3n6 82.34
# 4   C20_4n3 C22_3n6 82.26
# 5  C18_2n6t C22_3n6 82.23
# 6     C16_0 C22_3n6 82.22
# 7   C22_3n6 C24_1n9 82.21
# 8   C22_3n6 C18_2nX 82.21
# 9   C20_5n3 C22_3n6 82.20
# 10  C22_3n6 C20_2nX 82.20
FA_cons_LR3 = cbind(FA_cons_LR2, FA_cons_step3$logratios.top[,1])
FA_cons_NR3 = rbind(FA_cons_NR2, FA_cons_step3$ratios.top[1,])

# Step 4
FA_cons_step4 = STEP(data = ma_FA_cons, nsteps = 1, top = 10, previous = FA_cons_LR3)
print.ratios(FA_cons_step4$names.top, FA_cons_step4$R2.top)
# FA1     FA2    R2
# 1  C22_3n6 C22_6n3 85.95
# 2  C18_4n3 C22_6n3 85.95 <- SDA for browns, DHA for dinos?
# 3  C18_3n3 C22_6n3 85.89
# 4    C14_0 C22_6n3 85.81
# 5  C20_4n3 C22_6n3 85.77
# 6    C16_0 C22_6n3 85.65
# 7  C22_6n3 C20_1nX 85.63
# 8  C18_1n9 C22_6n3 85.59
# 9  C18_3n6 C22_6n3 85.57
# 10 C20_5n3 C22_6n3 85.55
FA_cons_LR4 = cbind(FA_cons_LR3, FA_cons_step4$logratios.top[,2])
FA_cons_NR4 = rbind(FA_cons_NR3, FA_cons_step4$ratios.top[2,])

# Step 5
FA_cons_step5 = STEP(data = ma_FA_cons, nsteps = 1, top = 10, previous = FA_cons_LR4)
print.ratios(FA_cons_step5$names.top, FA_cons_step5$R2.top)
# FA1      FA2    R2
# 1  C22_1n9c  C16_1nX 88.46
# 2     C14_0  C16_1nX 88.40
# 3   C16_1nX  C20_1nX 88.37 <- 16 MUFAs have been used for microbes/bacteria, 20 MUFAs for zoops/copepods
# 4   C20_1n9  C16_1nX 88.30
# 5   C20_3n3  C16_1nX 88.27
# 6   C16_1n7  C16_1nX 88.23
# 7     C12_0  C20_3n3 88.22
# 8     C12_0  C20_1nX 88.19
# 9     C12_0 C22_1n9c 88.19
# 10  C18_1n9  C20_1nX 88.18
FA_cons_LR5 = cbind(FA_cons_LR4, FA_cons_step5$logratios.top[,3])
FA_cons_NR5 = rbind(FA_cons_NR4, FA_cons_step5$ratios.top[3,])

# Step 6
FA_cons_step6 = STEP(data = ma_FA_cons, nsteps = 1, top = 10, previous = FA_cons_LR5)
print.ratios(FA_cons_step6$names.top, FA_cons_step6$R2.top)
# FA1     FA2    R2
# 1  C18_1n9 C20_4n6 90.46 <- Oleic for detritus, ARA for browns (and reds-ish)
# 2  C18_1n9 C18_3n6 90.31
# 3    C14_0 C20_4n6 90.28
# 4    C14_0 C22_5n3 90.27
# 5    C14_0 C20_2nX 90.27
# 6    C14_0 C20_5n3 90.27
# 7  C18_1n9 C20_2nX 90.27
# 8  C18_1n9 C20_5n3 90.27
# 9  C18_1n9 C22_5n3 90.27
# 10   C16_0 C20_4n6 90.24
FA_cons_LR6 = cbind(FA_cons_LR5, FA_cons_step6$logratios.top[,1])
FA_cons_NR6 = rbind(FA_cons_NR5, FA_cons_step6$ratios.top[1,])

# Step 7
FA_cons_step7 = STEP(data = ma_FA_cons, nsteps = 1, top = 10, previous = FA_cons_LR6)
print.ratios(FA_cons_step7$names.top, FA_cons_step7$R2.top)
# FA1     FA2    R2
# 1  C22_1n9c C22_4n6 92.12
# 2   C20_1n9 C22_4n6 92.08
# 3  C18_2n6c C22_4n6 92.07 <- LIN for greens, docosatetraenoic for...?
# 4   C20_3n3 C22_4n6 92.05
# 5   C22_4n6 C18_1nX 92.04
# 6   C18_3n6 C22_4n6 92.04
# 7   C22_4n6 C20_1nX 92.03
# 8   C22_4n6 C16_1nX 92.03
# 9     C14_0 C22_4n6 92.02
# 10  C20_3n6 C22_4n6 92.01
FA_cons_LR7 = cbind(FA_cons_LR6, FA_cons_step7$logratios.top[,3])
FA_cons_NR7 = rbind(FA_cons_NR6, FA_cons_step7$ratios.top[3,])

# Step 8
FA_cons_step8 = STEP(data = ma_FA_cons, nsteps = 1, top = 10, previous = FA_cons_LR7)
print.ratios(FA_cons_step8$names.top, FA_cons_step8$R2.top)
# FA1     FA2    R2
# 1   C18_1n9 C20_1n9 93.34 <- oleic for detritus/animal material, 20:1n9 suggestive of copes/zoops
# 2   C20_1n9 C20_4n6 93.34
# 3   C20_1n9 C22_4n3 93.32
# 4   C18_3n3 C20_1n9 93.28
# 5   C20_1n9   OBCFA 93.27
# 6     C14_0 C20_4n6 93.26
# 7     C14_0 C18_1n9 93.26
# 8     C16_0 C20_1n9 93.25
# 9  C18_2n6c C20_1n9 93.25
# 10  C20_1n9 C22_4n6 93.25
FA_cons_LR8 = cbind(FA_cons_LR7, FA_cons_step8$logratios.top[,1])
FA_cons_NR8 = rbind(FA_cons_NR7, FA_cons_step8$ratios.top[1,])

# Step 9
FA_cons_step9 = STEP(data = ma_FA_cons, nsteps = 1, top = 10, previous = FA_cons_LR8)
print.ratios(FA_cons_step9$names.top, FA_cons_step9$R2.top)
# FA1      FA2    R2
# 1  C22_1n9c  C22_2nX 94.41
# 2     C14_0    C18_0 94.38
# 3  C18_2n6t C22_1n9c 94.36
# 4   C20_3n3  C22_2nX 94.36
# 5  C18_2n6c  C22_2nX 94.35
# 6   C22_2nX  C22_4n6 94.35
# 7     C14_0  C20_1n9 94.35
# 8     C14_0  C18_1n9 94.35 <- myristic is common enough SAFA, oleic as common MUFA
# 9     C14_0  C20_4n6 94.35
# 10    C18_0 C18_2n6c 94.35
FA_cons_LR9 = cbind(FA_cons_LR8, FA_cons_step9$logratios.top[,8])
FA_cons_NR9 = rbind(FA_cons_NR8, FA_cons_step9$ratios.top[8,])

# Step 10
FA_cons_step10 = STEP(data = ma_FA_cons, nsteps = 1, top = 10, previous = FA_cons_LR9)
print.ratios(FA_cons_step10$names.top, FA_cons_step10$R2.top)
# FA1      FA2    R2
# 1  C22_1n9c  C22_2nX 95.39
# 2   C22_2nX  C22_4n6 95.37
# 3  C18_2n6c  C22_2nX 95.37
# 4   C20_3n3  C22_2nX 95.37
# 5   C20_3n6  C22_2nX 95.37
# 6     C18_0 C18_2n6c 95.34
# 7     C18_0  C22_4n6 95.34
# 8   C18_3n6  C22_2nX 95.34
# 9   C22_2nX  C22_3n6 95.34
# 10  C22_2nX  C22_6n3 95.34
# Stop here. Diminishing returns and this 22:2 TwoFA has limited interpretability 
# Captured 94.35% of the variation!

rownames(FA_cons_NR9) = paste("Step", 1:9, sep="")
colnames(FA_cons_NR9) = c("FA1","FA2")
FA_cons_FR = as.data.frame(cbind(FA_cons_NR9,
                                 Ratio = paste(colnames(ma_FA_cons)[FA_cons_NR9[,1]],
                                               "/",
                                               colnames(ma_FA_cons)[FA_cons_NR9[,2]],
                                               sep="")))
FA_cons_FR

colnames(FA_cons_LR9) = FA_cons_FR[,3]

FA_cons_PiR = sort(unique(as.numeric(FA_cons_NR9))) # "parts in ratios"
colnames(ma_FA_cons)[FA_cons_PiR]

FA_cons_PCA = PCA(FA_cons_LR9, weight = FALSE)

PLOT.PCA(FA_cons_PCA,
         map = "contribution",
         axes.inv = c(1, 1),
         rescale = 2)

# Extract and calculate contribution vectors
axes_inv = c(1, 1)
dim = c(1, 2)
FA_cons_RPC = FA_cons_PCA$rowcoord[, dim] %*% diag(FA_cons_PCA$sv[dim] * axes_inv)
FA_cons_CSC = FA_cons_PCA$colcoord[, dim] %*% diag(axes_inv)
FA_cons_CPC = FA_cons_CSC %*% diag(FA_cons_PCA$sv[dim])
FA_cons_CCC = FA_cons_CSC * sqrt(FA_cons_PCA$colmass)
FA_cons_CRD = FA_cons_CCC

# rescale contribution vectors for visibility
vectorscale = 2
FA_cons_CRD_scale = FA_cons_CRD * vectorscale

# Extract points
df_FA_cons_RPC = cbind(ma_FA_cons_head, FA_cons_RPC)
colnames(df_FA_cons_RPC)=c(colnames(ma_FA_cons_head),'dim1','dim2')

# Aggregate by season and sp
df_FA_cons_RPC_MSE = df_FA_cons_RPC %>%
  group_by(season, spCode) %>%
  summarize(dim1M = mean(dim1),
            dim2M = mean(dim2),
            dim1SE = sd(dim1) / sqrt(length(dim1)),
            dim2SE = sd(dim2) / sqrt(length(dim2)))

# Calculate variance explained by PCA dimensions
FA_cons_PCA_perc_1 = 100 * FA_cons_PCA$sv[dim[1]]^2 / sum(FA_cons_PCA$sv^2)
FA_cons_PCA_perc_2 = 100 * FA_cons_PCA$sv[dim[2]]^2 / sum(FA_cons_PCA$sv^2)

df_FA_cons_RPC$ID2 = paste(df_FA_cons_RPC$spCode,
                           rownames(df_FA_cons_RPC),
                           sep = "_")

# Format arrows
df_FA_cons_arrows = data.frame(dim1 = FA_cons_CRD_scale[ ,1],
                               dim2 = FA_cons_CRD_scale[, 2],
                               LR = FA_cons_PCA$colnames)

df_FA_cons_arrows$angle = 90 - ((atan2(df_FA_cons_arrows$dim1,
                                       df_FA_cons_arrows$dim2) *
                                   180) /
                                  pi)
df_FA_cons_arrows$angle = ifelse(df_FA_cons_arrows$angle > 90 & df_FA_cons_arrows$angle < 270,
                                 df_FA_cons_arrows$angle + 180,
                                 df_FA_cons_arrows$angle)
df_FA_cons_arrows$hjust = ifelse(df_FA_cons_arrows$angle > 90,
                                 1,
                                 0)

df_FA_cons_arrows$LR = gsub("C", "", df_FA_cons_arrows$LR)
df_FA_cons_arrows$LR = sub("c", "", df_FA_cons_arrows$LR)
df_FA_cons_arrows$LR = gsub("_", ":", df_FA_cons_arrows$LR)
df_FA_cons_arrows$LR = gsub("n", "*omega*-", df_FA_cons_arrows$LR)

df_FA_cons_RPC = left_join(df_FA_cons_RPC,
                           spCode_add[,c("spCode","binom")],
                           by = join_by(spCode))
df_FA_cons_RPC_MSE = left_join(df_FA_cons_RPC_MSE,
                               spCode_add[,c("spCode","binom")],
                               by = join_by(spCode))

# > > 2.3.2: Plot PCA -----
FA_cons_PCA_plot = ggplot() +
  geom_point(aes(x = dim1,
                 y = dim2,
                 color = season,
                 shape = binom),
             stroke = 0.5,
             alpha = 0.5,
             data = df_FA_cons_RPC) +
  # geom_text(aes(x = dim1,
  #               y = dim2,
  #               color = season,
  #               label = ID2),
  #           size = 2,
  #           data = df_FA_cons_RPC) +
  geom_segment(aes(x = 0,
                   y = 0,
                   xend = FA_cons_CRD_scale[ ,1],
                   yend = FA_cons_CRD_scale[ ,2]),
               color = "grey50",
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(aes(x = dim1 * 1.1,
                y = dim2 * 1.1,
                label = LR),
            colour = "grey50",
            size = 2,
            parse = TRUE,
            hjust = df_FA_cons_arrows$hjust,
            angle = df_FA_cons_arrows$angle,
            data = df_FA_cons_arrows)+
  geom_errorbar(aes(x = dim1M,
                    ymin = dim2M - 1.96 * dim2SE,
                    ymax = dim2M + 1.96 * dim2SE,
                    color = season),
                width = 0,
                linewidth = 0.5,
                data = df_FA_cons_RPC_MSE) +
  geom_errorbarh(aes(y = dim2M,
                     xmin = dim1M - 1.96 * dim1SE,
                     xmax = dim1M + 1.96 * dim1SE,
                     color = season),
                 height = 0,
                 linewidth = 0.5,
                 data = df_FA_cons_RPC_MSE) +
  geom_point(aes(x = dim1M,
                 y = dim2M,
                 color = season,
                 shape = binom),
             size = 3,
             stroke = 1,
             data = df_FA_cons_RPC_MSE) +
  scale_shape_manual(values = c(0, 1, 2, 4, 5, 6)) +
  scale_colour_manual(values = c("#9B1F1A", "#5FA1F7"), 
                      labels = c("Summer", "Winter")) + 
  scale_y_continuous(sec.axis = ~./vectorscale) +
  scale_x_continuous(sec.axis = ~./vectorscale) +
  labs(x = paste("PC 1 (",
                 round(FA_cons_PCA_perc_1,1),
                 "%)",
                 sep=""),
       y = paste("PC 2 (",
                 round(FA_cons_PCA_perc_2,1),
                 "%)",
                 sep=""),
       color = "Season",
       shape = "Species") +
  coord_fixed() +
  theme_classic() +
  guides(shape = guide_legend(label.theme = element_text(face = "italic",
                                                         size = unit(10, "pt"))),
         color = guide_legend(label.theme = element_text(size = unit(10, "pt")))) +
  theme(axis.text.y.right = element_text(colour="grey50"),
        axis.text.x.top = element_text(colour="grey50"),
        axis.ticks.y.right = element_line(colour="grey50"),
        axis.ticks.x.top = element_line(colour="grey50"),
        axis.line.y.right = element_line(colour="grey50"),
        axis.line.x.top = element_line(colour="grey50"))

FA_cons_PCA_plot

ggsave("pHAKFA_Fig_03.pdf",
       device = "pdf",
       height = 12,
       width = 18,
       units = "cm")

# > > 2.3.3: Univariate plot -----
FA_cons_LR9_long = cbind(ma_FA_cons_head, FA_cons_LR9)
FA_cons_LR9_long = pivot_longer(FA_cons_LR9_long,
                                cols = (length(FA_cons_LR9_long)-8):length(FA_cons_LR9_long),
                                names_to = "FAs",
                                values_to = "logratio")

FA_cons_LR9_long = left_join(FA_cons_LR9_long,
                              spCode_add[,c("spCode","binom")],
                              by = join_by(spCode))

FA_cons_LR9_long$binom = factor(FA_cons_LR9_long$binom,
                                 levels = c("M. franciscanus",
                                            "S. droebachiensis",
                                            "H. kamtchatkana",
                                            "T. pulligo",
                                            "P. resecata",
                                            "P. producta"))

FA_cons_LR9_MSE = FA_cons_LR9_long %>%
  group_by(season, binom, FAs) %>%
  summarize(M = mean(logratio),
            SE = sd(logratio)/sqrt(length(logratio)))

FA_cons_LR9_MSE2 = pivot_wider(FA_cons_LR9_MSE,
                               id_cols = c(FAs, binom),
                               names_from = season,
                               values_from = M)
FA_cons_LR9_MSE2$diff = FA_cons_LR9_MSE2$summer - FA_cons_LR9_MSE2$winter
FA_cons_LR9_MSE2$diff_abs = abs(FA_cons_LR9_MSE2$diff)

FA_cons_LR9_MSE2 = FA_cons_LR9_MSE2 %>%
  group_by(FAs) %>%
  summarize(diff_abs = mean(diff_abs))

FA_cons_LR9_long$FAs = gsub("C", "", FA_cons_LR9_long$FAs)
FA_cons_LR9_long$FAs = sub("c", "", FA_cons_LR9_long$FAs)
FA_cons_LR9_long$FAs = gsub("_", ":", FA_cons_LR9_long$FAs)
FA_cons_LR9_long$FAs = gsub("n", "n-", FA_cons_LR9_long$FAs)

FA_cons_LR9_MSE$FAs = gsub("C", "", FA_cons_LR9_MSE$FAs)
FA_cons_LR9_MSE$FAs = sub("c", "", FA_cons_LR9_MSE$FAs)
FA_cons_LR9_MSE$FAs = gsub("_", ":", FA_cons_LR9_MSE$FAs)
FA_cons_LR9_MSE$FAs = gsub("n", "n-", FA_cons_LR9_MSE$FAs)

FA_cons_LR9_MSE2$FAs = gsub("C", "", FA_cons_LR9_MSE2$FAs)
FA_cons_LR9_MSE2$FAs = sub("c", "", FA_cons_LR9_MSE2$FAs)
FA_cons_LR9_MSE2$FAs = gsub("_", ":", FA_cons_LR9_MSE2$FAs)
FA_cons_LR9_MSE2$FAs = gsub("n", "n-", FA_cons_LR9_MSE2$FAs)

FA_cons_LR9_long$FAs = factor(FA_cons_LR9_long$FAs,
                              levels = FA_cons_LR9_MSE2$FAs[order(FA_cons_LR9_MSE2$diff_abs,
                                                                  decreasing = TRUE)])
FA_cons_LR9_MSE$FAs = factor(FA_cons_LR9_MSE$FAs,
                             levels = FA_cons_LR9_MSE2$FAs[order(FA_cons_LR9_MSE2$diff_abs,
                                                                 decreasing = TRUE)])

FA_cons_univ_plot = ggplot(aes(x = binom,
                               y = logratio,
                               color = season),
                           data = FA_cons_LR9_long) +
  geom_point(size = 1,
             alpha = 0.2,
             position = position_dodge(.9)) +
  geom_point(aes(x = binom,
                 y = M,
                 color = season),
             size = 1,
             position = position_dodge(.9),
             data = FA_cons_LR9_MSE) +
  geom_errorbar(aes(x = binom,
                    y = M,
                    ymin = M - SE * 1.96,
                    ymax = M + SE * 1.96,
                    color = season),
                position = position_dodge(0.9),
                data = FA_cons_LR9_MSE) +
  scale_colour_manual(values = c("#9B1F1A", "#5FA1F7"), 
                      labels = c("Summer", "Winter")) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   size = unit(6, "pt"),
                                   face = "italic"),
        axis.text.y = element_text(size = unit(6, "pt")),
        axis.title = element_text(size = unit(8, "pt")),
        strip.text = element_text(size = unit(5, "pt")),
        legend.title = element_text(size = unit(8, "pt")),
        legend.text = element_text(size = unit(6, "pt")),
        legend.position = "bottom") +
  facet_grid(. ~ FAs)
FA_cons_univ_plot

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

# 3: Analyses ------------------------------------------------------------------

# > 3.1: Producers -----
ma_FA_prod_head %>%
  group_by(spCode, season) %>%
  summarize(count = length(spCode))

# two of the 12 have only 3 samples, so we sample 3 from each remaining
# number of ways to leave one out of each group: 4^10 = 1048576

ma_FA_prod_head$ind = 1:length(ma_FA_prod_head$spCode)

FA_prod_LR8_full = cbind(subset(ma_FA_prod_head,
                                select = -c(genus,
                                            species,
                                            lengthDia,
                                            height,
                                            lengthHeightUnits,
                                            wetWeightG)),
                         FA_prod_LR8)

# > > 3.1.1: PERMANOVA -----

iterations = 1000 
PERMANOVA_prod_LR8_iter = data.frame(Df = rep(NA, iterations),
                                     SumOfSqs = rep(NA, iterations),
                                     R2 = rep(NA, iterations),
                                     pseudoF = rep(NA, iterations),
                                     P = rep(NA, iterations))

for (i in 1:iterations) {
  FA_prod_LR8_sub = FA_prod_LR8_full %>%
    group_by(season, spCode) %>%
    slice_sample(n = 3,
                 replace = FALSE)
  
  FA_prod_LR8_head = FA_prod_LR8_sub[,1:7]
  FA_prod_LR8_mat = as.matrix(FA_prod_LR8_sub[,8:length(FA_prod_LR8_sub)])
  
  permdes_prod_LR8 = how(within = Within(type = "free"), # samples are free to shuffle across seasons
                         plots = Plots(strata = FA_prod_LR8_head$spCode, type = "free"), # species within phyla are allowed to shuffle, but samples within species cannot shuffle
                         blocks = FA_prod_LR8_head$divphy, # samples are not allowed to shuffle across phyla
                         nperm = 999) 
  
  PERMANOVA_prod_LR8 = adonis2(FA_prod_LR8_mat ~ divphy + spCode + season,
                               method = "euclidean",
                               by = "terms", # sequential term significance vs. marginal significance (by = "margin")
                               data = FA_prod_LR8_head,
                               permutations = permdes_prod_LR8) # Uses our defined permutation structure
  
  PERMANOVA_prod_LR8_iter[i,]$Df = list(PERMANOVA_prod_LR8$Df)
  PERMANOVA_prod_LR8_iter[i,]$SumOfSqs = list(PERMANOVA_prod_LR8$SumOfSqs)
  PERMANOVA_prod_LR8_iter[i,]$R2 = list(PERMANOVA_prod_LR8$R2)
  PERMANOVA_prod_LR8_iter[i,]$pseudoF = list(PERMANOVA_prod_LR8$F)
  PERMANOVA_prod_LR8_iter[i,]$P = list(PERMANOVA_prod_LR8$`Pr(>F)`)
  
}

# Df for residual (4) and total (5)
median(sapply(PERMANOVA_prod_LR8_iter$Df, "[[", 4))
median(sapply(PERMANOVA_prod_LR8_iter$Df, "[[", 5))

# Sum of squares for residual (4) and total (5)
median(sapply(PERMANOVA_prod_LR8_iter$SumOfSqs, "[[", 4))
median(sapply(PERMANOVA_prod_LR8_iter$SumOfSqs, "[[", 5))

# First element is divphy
hist(sapply(PERMANOVA_prod_LR8_iter$R2, "[[", 1))
# Second element is spCode
hist(sapply(PERMANOVA_prod_LR8_iter$R2, "[[", 2))
# Third element is season
hist(sapply(PERMANOVA_prod_LR8_iter$R2, "[[", 3))

df_PERMANOVA_prod_LR8_iter = data.frame(variable = rep(c("divphy",
                                                         "spCode",
                                                         "season"),
                                                       each = iterations),
                                        Df = c(sapply(PERMANOVA_prod_LR8_iter$Df, "[[", 1),
                                               sapply(PERMANOVA_prod_LR8_iter$Df, "[[", 2),
                                               sapply(PERMANOVA_prod_LR8_iter$Df, "[[", 3)),
                                        SumOfSqs = c(sapply(PERMANOVA_prod_LR8_iter$SumOfSqs, "[[", 1),
                                                     sapply(PERMANOVA_prod_LR8_iter$SumOfSqs, "[[", 2),
                                                     sapply(PERMANOVA_prod_LR8_iter$SumOfSqs, "[[", 3)),
                                        R2 = c(sapply(PERMANOVA_prod_LR8_iter$R2, "[[", 1),
                                               sapply(PERMANOVA_prod_LR8_iter$R2, "[[", 2),
                                               sapply(PERMANOVA_prod_LR8_iter$R2, "[[", 3)),
                                        pseudoF = c(sapply(PERMANOVA_prod_LR8_iter$pseudoF, "[[", 1),
                                                    sapply(PERMANOVA_prod_LR8_iter$pseudoF, "[[", 2),
                                                    sapply(PERMANOVA_prod_LR8_iter$pseudoF, "[[", 3)),
                                        P = c(sapply(PERMANOVA_prod_LR8_iter$P, "[[", 1),
                                              sapply(PERMANOVA_prod_LR8_iter$P, "[[", 2),
                                              sapply(PERMANOVA_prod_LR8_iter$P, "[[", 3)))

df_PERMANOVA_prod_LR8_iter$variable = factor(df_PERMANOVA_prod_LR8_iter$variable,
                                             levels = c("divphy",
                                                        "spCode",
                                                        "season"))

summ_PERMANOVA_prod_LR8_iter = df_PERMANOVA_prod_LR8_iter %>%
  group_by(variable) %>%
  summarize(Df = unique(Df),
            mean_SumOfSqs = mean(SumOfSqs),
            med_SumOfSqs = median(SumOfSqs),
            mean_R2 = mean(R2),
            med_R2 = median(R2),
            mean_pseudoF = mean(pseudoF),
            med_pseudoF = median(pseudoF),
            mean_P = mean(P),
            med_P = median(P))

ggplot() +
  geom_histogram(aes(x = P),
                 data = df_PERMANOVA_prod_LR8_iter) +
  geom_vline(aes(xintercept = med_P),
             linetype = 2,
             data = summ_PERMANOVA_prod_LR8_iter) +
  labs(x = "p-value",
       y = "number of iterations") +
  theme_bw() +
  theme() +
  facet_grid(. ~ variable,
             scales = "free_x")

plot_R2_dens = ggplot() +
  geom_density(aes(x = R2,
                   fill = variable,
                   color = variable),
               alpha = 0.5,
               data = df_PERMANOVA_prod_LR8_iter) +
  geom_vline(aes(xintercept = med_R2,
                 color = variable),
             linetype = 2,
             data = summ_PERMANOVA_prod_LR8_iter) +
  geom_text(aes(x = med_R2 + 0.05,
                y = 30,
                label = round(med_R2, 2),
                color = variable),
            data = summ_PERMANOVA_prod_LR8_iter) +
  scale_color_viridis_d(option = "mako",
                        begin = 0.0,
                        end = 0.8,
                        labels = c("division",
                                   "species",
                                   "season")) +
  scale_fill_viridis_d(option = "mako",
                       begin = 0.0,
                       end = 0.8,
                       labels = c("division",
                                  "species",
                                  "season")) +
  labs(x = "R2",
       y = "number of iterations") +
  theme_bw() +
  theme(legend.position = "none")
plot_R2_dens

plot_p_dens = ggplot() +
  geom_density(aes(x = P,
                   fill = variable,
                   color = variable),
               alpha = 0.5,
               data = df_PERMANOVA_prod_LR8_iter) +
  geom_vline(aes(xintercept = med_P),
             color = "black",
             linetype = 2,
             data = summ_PERMANOVA_prod_LR8_iter) +
  geom_text(aes(x = med_P + 0.0025,
                y = 300,
                label = round(med_P, 3)),
            color = "black",
            data = summ_PERMANOVA_prod_LR8_iter) +
  scale_color_viridis_d(option = "mako",
                        begin = 0.0,
                        end = 0.8,
                        labels = c("division",
                                   "species",
                                   "season")) +
  scale_fill_viridis_d(option = "mako",
                       begin = 0.0,
                       end = 0.8,
                       labels = c("division",
                                  "species",
                                  "season")) +
  labs(x = "p-value",
       y = NULL) +
  guides(fill = guide_legend(position = "inside"),
         color = guide_legend(position = "inside")) +
  theme_bw() +
  theme(legend.position.inside = c(0.80, 0.80))
plot_p_dens

plot_R2_dens | plot_p_dens

# png("plot_pHAKFA_LR_prod_PERMANOVA.png",
#     height = 700 * 4,
#     width = 700 * 7,
#     res = 700)
# plot_R2_dens | plot_p_dens
# dev.off()
# 
# tiff("plot_pHAKFA_LR_prod_PERMANOVA.tiff",
#     height = 700 * 4,
#     width = 700 * 7,
#     res = 700)
# plot_R2_dens | plot_p_dens
# dev.off()

# # > > 3.1.2: PERMDISP2 -----

# Use all the algae data (no sampling down)
FA_prod_LR8_head = FA_prod_LR8_full[,1:7]
FA_prod_LR8_mat = as.matrix(FA_prod_LR8_full[,8:length(FA_prod_LR8_full)])


PERMDISP_prod_LR8_spCode = betadisper(vegdist(FA_prod_LR8_mat,
                                               method = "euclidean"),
                                       bias.adjust = TRUE,
                                       FA_prod_LR8_head$spCode)
PERMDISP_prod_LR8_spCode
anova(PERMDISP_prod_LR8_spCode)


PERMDISP_prod_LR8_season = betadisper(vegdist(FA_prod_LR8_mat,
                                              method = "euclidean"),
                                      bias.adjust = TRUE,
                                      FA_prod_LR8_head$season)
PERMDISP_prod_LR8_season
anova(PERMDISP_prod_LR8_season)

PERMDISP_prod_LR8_int = betadisper(vegdist(FA_prod_LR8_mat,
                                           method = "euclidean"),
                                   bias.adjust = TRUE,
                                   interaction(FA_prod_LR8_head$season,
                                               FA_prod_LR8_head$spCode))
PERMDISP_prod_LR8_int
anova(PERMDISP_prod_LR8_int)

# > 3.2: Herbivores w/ Algal logratios -----
ma_FA_cons_head %>%
  group_by(spCode, season) %>%
  summarize(count = length(spCode))
# The arthropods are the problems, so our best route is to omit them from analysis
# and discuss them qualitatively.

FA_cons_LR8P_full = cbind(subset(ma_FA_cons_head,
                                 select = -c(genus,
                                             species,
                                             lengthDia,
                                             height,
                                             lengthHeightUnits,
                                             wetWeightG)),
                          FA_cons_LR8P)

# isolate Pugettia
FA_cons_LR8P_Puge = FA_cons_LR8P_full[FA_cons_LR8P_full$spCode == "PUGETTIA",]

# drop the Arthropods
FA_cons_LR8P_full = FA_cons_LR8P_full[FA_cons_LR8P_full$divphy != "Arthropoda",]

# Check data without arthropods
FA_cons_LR8P_full %>%
  group_by(spCode, season) %>%
  summarize(count = length(spCode))

FA_cons_LR8P_head = FA_cons_LR8P_full[, 1:6]
FA_cons_LR8P_mat = as.matrix(FA_cons_LR8P_full[, 7:length(FA_cons_LR8P_full)])

FA_cons_LR8P_Puge_head = FA_cons_LR8P_Puge[, 1:6]
FA_cons_LR8P_Puge_mat = as.matrix(FA_cons_LR8P_Puge[, 7:length(FA_cons_LR8P_Puge)])

# > > 3.2.1: PERMANOVAs -----

# Molluscs and echinoderms
permdes_cons_LR8P = how(within = Within(type = "free"), # samples are free to shuffle across seasons
                       plots = Plots(strata = FA_cons_LR8P_head$spCode, type = "free"), # species within phyla are allowed to shuffle, but samples within species cannot shuffle
                       blocks = FA_cons_LR8P_head$divphy, # samples are not allowed to shuffle across phyla
                       nperm = 999) 

PERMANOVA_cons_LR8P = adonis2(FA_cons_LR8P_mat ~ divphy + spCode + season,
                             method = "euclidean",
                             by = "terms", # sequential term significance vs. marginal significance (by = "margin")
                             data = FA_cons_LR8P_head,
                             permutations = permdes_cons_LR8P) # Uses our defined permutation structure
PERMANOVA_cons_LR8P

# Pugettia
# season only
PERMANOVA_cons_LR8P_Puge = adonis2(FA_cons_LR8P_Puge_mat ~ season,
                                   method = "euclidean",
                                   by = "terms", # sequential term significance vs. marginal significance (by = "margin")
                                   data = FA_cons_LR8P_Puge_head)

PERMANOVA_cons_LR8P_Puge

# > > 3.2.2: PERMDISP2 -----

# Molluscs and echinoderms
PERMDISP_cons_LR8P_spCode = betadisper(vegdist(FA_cons_LR8P_mat,
                                                method = "euclidean"),
                                        FA_cons_LR8P_head$spCode)
PERMDISP_cons_LR8P_spCode
anova(PERMDISP_cons_LR8P_spCode)


PERMDISP_cons_LR8P_season = betadisper(vegdist(FA_cons_LR8P_mat,
                                               method = "euclidean"),
                                       FA_cons_LR8P_head$season)
PERMDISP_cons_LR8P_season
anova(PERMDISP_cons_LR8P_season)

PERMDISP_cons_LR8P_int = betadisper(vegdist(FA_cons_LR8P_mat,
                                            method = "euclidean"),
                                    interaction(FA_cons_LR8P_head$season,
                                                FA_cons_LR8P_head$spCode),
                                    bias.adjust = TRUE)
PERMDISP_cons_LR8P_int
anova(PERMDISP_cons_LR8P_int)

# Pugettia
# season only
PERMDISP_cons_LR8P_Puge_season = betadisper(vegdist(FA_cons_LR8P_Puge_mat,
                                                    method = "euclidean"),
                                            FA_cons_LR8P_Puge_head$season)
PERMDISP_cons_LR8P_Puge_season
anova(PERMDISP_cons_LR8P_Puge_season)

# > 3.3: Herbivores -----
FA_cons_LR9_full = cbind(subset(ma_FA_cons_head,
                                 select = -c(genus,
                                             species,
                                             lengthDia,
                                             height,
                                             lengthHeightUnits,
                                             wetWeightG)),
                          FA_cons_LR9)
# isolate Pugettia 
FA_cons_LR9_Puge = FA_cons_LR9_full[FA_cons_LR9_full$spCode == "PUGETTIA",]

# drop the Arthropods
FA_cons_LR9_full = FA_cons_LR9_full[FA_cons_LR9_full$divphy != "Arthropoda",]

# check sample sizes without arthropods
FA_cons_LR9_full %>%
  group_by(spCode, season) %>%
  summarize(count = length(spCode))

FA_cons_LR9_head = FA_cons_LR9_full[, 1:6]
FA_cons_LR9_mat = as.matrix(FA_cons_LR9_full[, 7:length(FA_cons_LR9_full)])

FA_cons_LR9_Puge_head = FA_cons_LR9_Puge[, 1:6]
FA_cons_LR9_Puge_mat = as.matrix(FA_cons_LR9_Puge[, 7:length(FA_cons_LR9_Puge)])

# > > 3.3.1: PERMANOVAs -----

# Molluscs and echinoderms
permdes_cons_LR9 = how(within = Within(type = "free"), # samples are free to shuffle across seasons
                        plots = Plots(strata = FA_cons_LR9_head$spCode, type = "free"), # species within phyla are allowed to shuffle, but samples within species cannot shuffle
                        blocks = FA_cons_LR9_head$divphy, # samples are not allowed to shuffle across phyla
                        nperm = 999) 

PERMANOVA_cons_LR9 = adonis2(FA_cons_LR9_mat ~ divphy + spCode + season,
                              method = "euclidean",
                              by = "terms", # sequential term significance vs. marginal significance (by = "margin")
                              data = FA_cons_LR9_head,
                              permutations = permdes_cons_LR9) # Uses our defined permutation structure

PERMANOVA_cons_LR9

# Pugettia
# season only
PERMANOVA_cons_LR9_Puge = adonis2(FA_cons_LR9_Puge_mat ~ season,
                                  method = "euclidean",
                                  by = "terms", # sequential term significance vs. marginal significance (by = "margin")
                                  data = FA_cons_LR9_Puge_head)

PERMANOVA_cons_LR9_Puge

# > > 3.3.2: PERMDISP2 -----
PERMDISP_cons_LR9_spCode = betadisper(vegdist(FA_cons_LR9_mat,
                                               method = "euclidean"),
                                       FA_cons_LR9_head$spCode)
PERMDISP_cons_LR9_spCode
anova(PERMDISP_cons_LR9_spCode)


PERMDISP_cons_LR9_season = betadisper(vegdist(FA_cons_LR9_mat,
                                              method = "euclidean"),
                                      FA_cons_LR9_head$season)
PERMDISP_cons_LR9_season
anova(PERMDISP_cons_LR9_season)

PERMDISP_cons_LR9_int = betadisper(vegdist(FA_cons_LR9_mat,
                                           method = "euclidean"),
                                   interaction(FA_cons_LR9_head$season,
                                               FA_cons_LR9_head$spCode))
PERMDISP_cons_LR9_int
anova(PERMDISP_cons_LR9_int)

# Pugettia
# season only
PERMDISP_cons_LR9_Puge_season = betadisper(vegdist(FA_cons_LR9_Puge_mat,
                                                   method = "euclidean"),
                                           FA_cons_LR9_Puge_head$season)
PERMDISP_cons_LR9_Puge_season
anova(PERMDISP_cons_LR9_Puge_season)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####
# END OF SCRIPT ----------------------------------------------------------------
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####
