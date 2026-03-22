# ▷▷▷▷▷▷▷▷▷▷▷▷▷▷
# ▷ Load packages   ▷
# ▷▷▷▷▷▷▷▷▷▷▷▷▷▷
lapply(c("bnlearn", "dplyr", "progress", "tidyverse"), library, character.only = TRUE)

# Read data from GSEA results
lnc2ssg2 <- read.table("lnc2ssg_result.txt", header = TRUE)
lnc2imm2 <- read.table("lnc2imm_result.txt", header = TRUE)
ssg2imm2 <- read.table("ssg2imm_result.txt", header = TRUE)

load("SSG_exp.RData")
load("lnc_exp.RData")
load("IMM_exp.RData")

# Prepare gene pairs (lncRNA-SSG)
gene_pair <- lnc2ssg2[, c(1, 2)]

# Function to count immune genes associated with a given lncRNA
imm_table <- function(x) {
  imm <- lnc2imm2[which(lnc2imm2$lnc %in% x[1]), 2]
  return(length(imm))
}

system.time({
  lnc2ssg_imm <- apply(gene_pair, 1, imm_table)
})
lnc2ssg_imm <- as.data.frame(lnc2ssg_imm)
lnc2ssg_imm <- cbind(gene_pair, lnc2ssg_imm)
colnames(lnc2ssg_imm) <- c("lnc", "ssg", "imm_num")

####  Extract qualified lnc-ssg pairs containing immune genes based on Bayesian inference ####
permutation_LS <- read.csv("permutation_result.csv", header = TRUE)
permutation_LS$P <- permutation_LS$num / 1000
# FDR correction
permutation_LS$FDR <- p.adjust(permutation_LS$P, method = "BH")
permutation_LS$imm_num <- lnc2ssg_imm$imm_num

LIS <- permutation_LS[which(permutation_LS$FDR < 0.05), c(1, 2)]  # Pairs containing immune genes
LS <- permutation_LS[which(permutation_LS$FDR >= 0.05), c(1, 2)]   # Pairs not containing immune genes

write.table(LIS, paste0("./", cancer_abbr, "/11.包含免疫基因/LIS.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(LS, paste0("./", cancer_abbr, "/11.包含免疫基因/LS.txt"), quote = FALSE, row.names = FALSE, sep = "\t")

# Update lnc2ssg pairs
lnc2ssg <- read.table(paste0("./", cancer_abbr, "/11.包含免疫基因/LIS.txt"), sep = "\t", header = TRUE)

# Build triplets
lnc2imm <- lnc2imm2[, 1:2]
lnc2imm <- lnc2imm[lnc2imm$lnc %in% LIS$lnc, ]
unique(lnc2imm$lnc)

ssg2imm <- ssg2imm2[, 1:2]
ssg2imm <- ssg2imm[ssg2imm$ssg %in% LIS$ssg, ]
unique(ssg2imm$ssg)

lnc_ssg_imm <- merge(lnc2ssg, ssg2imm, by = "ssg")
lnc_imm_ssg <- merge(lnc2imm, ssg2imm, by = "imm")
ssg_lnc_imm <- merge(lnc2ssg, lnc2imm, by = "lnc")

a <- c("lnc", "imm", "ssg")
a <- sort(a)

lnc_ssg_imm <- lnc_ssg_imm[, a]
lnc_imm_ssg <- lnc_imm_ssg[, a]
ssg_lnc_imm <- ssg_lnc_imm[, a]

triplet <- rbind(lnc_ssg_imm, lnc_imm_ssg, ssg_lnc_imm)
triplet <- unique(triplet)  # Ensure uniqueness

#### Bayesian inference ####
modeSelection <- function(lnc, ssg, imm) {
  learning.test <- data.frame(M = lnc, E = ssg, G = imm)  # Data for building Bayesian network
  
  # Define models
  models <- list(
    LIS = model2network("[M][G|M][E|G]"),
    LSI = model2network("[M][E|M][G|E]"),
    IR = model2network("[M][G|M][E|M]"),
    CR = model2network("[M][G][E|M:G]")
  )
  
  # Calculate BIC scores
  BICS <- sapply(models, function(model) score(model, learning.test, type = "bic-g"))
  
  # Determine best model
  BICS <- BICS[order(BICS, decreasing = FALSE)]
  BIC.min <- min(BICS)
  
  # Calculate relative likelihood and model weights
  dBIC <- BICS - BIC.min
  w <- exp(-0.5 * dBIC) / sum(exp(-0.5 * dBIC))
  
  model <- names(which.max(w))
  r <- c(
    model = model,
    BICS,
    RL_BIC = exp(0.5 * (min(BICS) - BICS[2])),  # Minimum minus second smallest
    dBIC = dBIC,
    weight = w,
    ME_mi = ci.test(lnc, ssg, test = "mi-g")$statistic,  # lnc vs ssg
    ME_p = ci.test(lnc, ssg, test = "mi-g")$p.value,
    MG_mi = ci.test(lnc, imm, test = "mi-g")$statistic,  # lnc vs Gene
    MG_p = ci.test(lnc, imm, test = "mi-g")$p.value,
    EG_mi = ci.test(ssg, imm, test = "mi-g")$statistic,  # ssg vs Gene
    EG_p = ci.test(ssg, imm, test = "mi-g")$p.value,
    CIT_ME_mi = ci.test(lnc, ssg, imm, test = "mi-g")$statistic,  # lnc, ssg | imm
    CIT_ME_p = ci.test(lnc, ssg, imm, test = "mi-g")$p.value,
    CIT_MG_mi = ci.test(lnc, imm, ssg, test = "mi-g")$statistic,  # lnc, imm | ssg
    CIT_MG_p = ci.test(lnc, imm, ssg, test = "mi-g")$p.value,
    CIT_EG_mi = ci.test(ssg, imm, lnc, test = "mi-g")$statistic,  # ssg, imm | lnc
    CIT_EG_p = ci.test(ssg, imm, lnc, test = "mi-g")$p.value
  )
  
  return(r)
}

# Model selection
pb <- progress_bar$new(total = nrow(triplet))

system.time({
  for (n in 1:nrow(triplet)) {
    pb$tick()
    
    # Extract expression data for the triplet
    lnc1 <- lnc_exp[triplet[n, "lnc"], ]
    ssg1 <- SSG_exp[triplet[n, "ssg"], ]
    imm1 <- IMM_exp[triplet[n, "imm"], ]
    
    # Convert to numeric and handle missing values
    lnc1 <- as.double(na.omit(lnc1))
    ssg1 <- as.double(na.omit(ssg1))
    imm1 <- as.double(na.omit(imm1))
    
    # Skip if any missing values
    if (length(lnc1) == 0 || length(ssg1) == 0 || length(imm1) == 0) next
    
    # Model selection
    aaa <- modeSelection(lnc = lnc1, ssg = ssg1, imm = imm1)
    
    # Store results in triplet data frame
    triplet[n, 4:29] <- aaa
  }
})

# 15. Result statistics
colnames(triplet)[4:29] <- c("model", "BIC_LIS", "BIC_LSI", "BIC_IR", "BIC_CR", "RL_BIC",
                             "dBIC_LIS", "dBIC_LSI", "dBIC_IR", "dBIC_CR",
                             "weight_LIS", "weight_LSI", "weight_IR", "weight_CR",
                             "MI_lncssg", "p_lncssg", "MI_lncimm", "p_lncimm", "MI_ssgimm", "p_ssgimm",
                             "MI_lncssg|imm", "p_lncssg|imm", "MI_lncimm|ssg", "p_lncimm|ssg",
                             "MI_ssgimm|lnc", "p_ssgimm|lnc")
write.table(triplet, "model_result.txt", quote = FALSE, row.names = FALSE, sep = "\t")

# Convert columns to numeric
triplet <- read.table("model_result.txt", sep = "\t", header = TRUE, check.names = FALSE)
triplet[, c(5:29)] <- apply(triplet[, c(5:29)], 2, as.numeric)

# ♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦
# ♦ Filter triplets by model type ♦
# ♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦
# IR model
IR_triplet <- triplet[which(triplet[, 4] == "IR"), ]  # IR model triplets
IR_triplet <- IR_triplet[which(IR_triplet$`p_ssgimm|lnc` < 0.05), ]  # Filter by p_ssgimm|lnc
a1 <- as.numeric(nrow(IR_triplet))  # Number of qualified triplets under IR model

# LIS model
LIS_triplet <- triplet[which(triplet[, 4] == "LIS"), ]  # LIS model triplets
LIS_triplet <- LIS_triplet[which(LIS_triplet$`p_lncssg|imm` < 0.05), ]  # Filter by p_lncssg|imm
a2 <- as.numeric(nrow(LIS_triplet))  # Number of qualified triplets under LIS model

# LSI model
LSI_triplet <- triplet[which(triplet[, 4] == "LSI"), ]  # LSI model triplets
LSI_triplet <- LSI_triplet[which(LSI_triplet$`p_lncimm|ssg` < 0.05), ]  # Filter by p_lncimm|ssg
a3 <- as.numeric(nrow(LSI_triplet))  # Number of qualified triplets under LSI model

# CR model
CR_triplet <- triplet[which(triplet[, 4] == "CR"), ]  # CR model triplets
CR_triplet <- CR_triplet[which(CR_triplet$`p_lncimm|ssg` < 0.05), ]  # Filter by p_lncimm|ssg
a4 <- as.numeric(nrow(CR_triplet))  # Number of qualified triplets under CR model

################
qualified_genepair <- rbind(CR_triplet, IR_triplet)
qualified_genepair <- rbind(qualified_genepair, LIS_triplet)
qualified_genepair <- rbind(qualified_genepair, LSI_triplet)  # Combine all qualified triplets
b <- as.numeric(nrow(qualified_genepair))  # Total number of qualified triplets

IR_percent <- a1 / b
LIS_percent <- a2 / b
LSI_percent <- a3 / b
CR_percent <- a4 / b

# Summarize filtering results
static_result <- data.frame(matrix(data = NA, nrow = 4, ncol = 5))
colnames(static_result) <- c("model", "qualified", "count", "pass_rate", "cancer_type")
static_result[, 1] <- c("IR", "LIS", "LSI", "CR")
static_result[, 2] <- c(a1, a2, a3, a4)
static_result[, 3] <- rep(b, 4)
static_result[, 4] <- c(IR_percent, LIS_percent, LSI_percent, CR_percent)
static_result[, 5] <- rep(cancer_abbr, 4)

write.table(static_result, "static_result.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Qualified gene pairs
qualified_genepair <- rbind(CR_triplet, IR_triplet)
qualified_genepair <- rbind(qualified_genepair, LIS_triplet)
qualified_genepair <- rbind(qualified_genepair, LSI_triplet)  # Combine all qualified triplets
write.table(qualified_genepair, "qualified_genepair.txt", quote = FALSE, row.names = FALSE, sep = "\t")

# Extract triplets qualified by both Bayesian and permutation
qualified_genepair <- read.table("qualified_genepair.txt", sep = "\t", header = TRUE)
LIS_genepair <- merge(LIS, qualified_genepair, by = c("lnc", "ssg"))

# Keep unique lnc-ssg pairs
LIS_genepair_unique <- LIS_genepair %>%
  select(lnc, ssg) %>%
  distinct()
LIS_genepair <- drop_na(LIS_genepair_unique)

# Count number of immune genes per lnc-ssg pair
LIS_gene <- LIS_genepair[, c(2, 1)]
LIS_gene <- LIS_gene[!duplicated(LIS_gene), ]
LIS_gene <- merge(permutation_LS, LIS_gene, all.y = TRUE)
LIS_gene <- drop_na(LIS_gene)

###
write.table(LIS_gene, "LIS_gene.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(LIS_genepair, "LIS_genepair.txt", sep = "\t", row.names = FALSE, quote = FALSE)