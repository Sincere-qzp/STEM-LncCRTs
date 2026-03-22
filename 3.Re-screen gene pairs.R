# Update interaction pairs
lnc2ssg2 <- lnc2ssg2[which(lnc2ssg2$lnc %in% lnc_gene$ENSG), ]
lnc2ssg2 <- lnc2ssg2[which(lnc2ssg2$ssg %in% SSG_gene$ENSG), ]

lnc2imm2 <- lnc2imm2[which(lnc2imm2$lnc %in% lnc_gene$ENSG), ]
lnc2imm2 <- lnc2imm2[which(lnc2imm2$imm %in% imm_gene$ENSG), ]

ssg2imm2 <- ssg2imm2[which(ssg2imm2$ssg %in% SSG_gene$ENSG), ]
ssg2imm2 <- ssg2imm2[which(ssg2imm2$imm %in% imm_gene$ENSG), ]

# ⛓⛓⛓⛓⛓⛓⛓⛓⛓⛓⛓⛓⛓⛓⛓⛓⛓
# ⛓ Use tumor purity score as covariate ⛓
# ⛓ Calculate partial correlation between lnc and SSG ⛓
# ⛓⛓⛓⛓⛓⛓⛓⛓⛓⛓⛓⛓⛓⛓⛓⛓⛓
gene_pair <- lnc2ssg2[, c(1, 2)]

#### Calculate tumor purity
# rforge <- "http://r-forge.r-project.org"
# install.packages ("estimate", repos = rforge, dependencies = TRUE)
exp <- read.table(file_to_read, row.names = 1, check.names = FALSE, header = TRUE, sep = "\t")
id <- rownames(exp)
exp <- cbind(id, exp)
write.table(exp, "ESTIMATE-EM_exp.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Calculate tumor purity score
filterCommonGenes(input.f = "ESTIMATE-EM_exp.txt",
                  output.f = "ESTIMATE-EM.gct",
                  id = "GeneSymbol")
estimateScore(input.ds = "ESTIMATE-EM.gct",
              output.ds = "estimate_score.gct",
              platform = "illumina")

# View results
scores <- read.table("estimate_score.gct",
                     skip = 2, header = TRUE)
rownames(scores) <- scores[, 1]
scores <- t(scores[, 3:ncol(scores)])
TumorPurity <- cos(0.6049872018 + 0.0001467884 * scores[, 3])
scores <- as.data.frame(scores)
scores$TumorPurity <- TumorPurity
write.table(scores, "TumorPurity.txt",
            quote = FALSE, sep = "\t")

### 1. Read tumor purity score file
purity_score <- read.table("TumorPurity.txt",
                           header = TRUE, sep = "\t")
purity_score <- data.frame(cbind(rownames(purity_score), purity_score$TumorPurity))
colnames(purity_score) <- c("V1", "TumorPurity")

### 2. Preprocess tumor purity score file
# Filter out normal samples
purity_score <- purity_score[-which(substr(purity_score$V1, start = 14, stop = 15) > 10), ]
# Keep sample names consistent with expression matrix format
purity_score$V1 <- gsub("\\.", "-", purity_score$V1)

### 3. Align sample order
lnc_exp <- as.data.frame(t(lnc_exp))  # rows: genes, columns: samples
lnc_exp <- lnc_exp[, purity_score$V1]

SSG_exp <- as.data.frame(t(SSG_exp))  # rows: genes, columns: samples
SSG_exp <- SSG_exp[, purity_score$V1]

IMM_exp <- as.data.frame(t(IMM_exp))  # rows: genes, columns: samples
IMM_exp <- IMM_exp[, purity_score$V1]

### 4. Calculate PCC coefficients and p-values
PCOR_LI <- function(x) {
  PCC <- pcor.test(as.numeric(lnc_exp[x[1], ]),
                   as.numeric(SSG_exp[x[2], ]),
                   as.numeric(purity_score$TumorPurity),
                   method = "pearson")
  r <- c(estimate = PCC$estimate, p = PCC$p.value)
  return(r)
}
PCC_LI <- apply(gene_pair, 1, PCOR_LI)
PCC_LI <- as.data.frame(t(PCC_LI))
PCC_LI <- cbind(gene_pair, PCC_LI)
PCC_LI$RS <- (-log10(PCC_LI$p)) * (sign(PCC_LI$estimate))
rm(PCOR_LI)

index <- which(PCC_LI$estimate > 0.5 & PCC_LI$p < 0.05)
PCC_LI <- PCC_LI[index, ]
lnc2ssg2 <- lnc2ssg2[index, ]  # Update lnc2ssg interaction pairs
rm(index, gene_pair)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>> Perform GSEA enrichment analysis, score lnc-ssg pairs >
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
### 1. Read pathway file related to immune genes
gmt <- read_excel("imm_genelist.xlsx")
colnames(gmt) <- c("ont", "gene")  # ont is biological pathway
length(table(gmt$ont))

### 2. Read IMM gene expression matrix
load("after_group_exp.RData")
gmtlast <- data.frame()

### 3. Convert Ensembl IDs to Symbol IDs in expression matrix
name <- as.data.frame(rownames(exp))
colnames(name) <- "name"
SYMBOL_ENSEMBL <- read.table("gene table.txt",
                             sep = "\t", header = TRUE)
name$name <- name$symbol

### 4. Filter IMM pathways present in expression matrix
for (i in 1:nrow(gmt)) {
  if (gmt$gene[i] %in% name$symbol) {
    gmtlast <- rbind(gmtlast, gmt[i, ])
  }
}
GO_term <- gmtlast[which(gmtlast$ont == "GO"), ]
GO_term <- GO_term[which(GO_term$gene %in% gene_top_200$ENSG), ]
gmtlast <- gmtlast[which(!(gmtlast$ont == "GO")), ]
gmtlast <- rbind(gmtlast, GO_term)
gmt <- gmtlast  # Genes present in GeneList
rm(gmtlast, GO_term)

### 4. Differential gene analysis
index <- which(substr(colnames(exp), start = 14, stop = 15) > 10)
group <- c(rep("normal", length(index)), rep("disease", ncol(exp) - length(index)))
targets <- cbind(colnames(exp), group)  # first column sampleID, second column group
targets <- as.data.frame(targets)
colnames(targets) <- c("FileName", "Target")
lev <- unique(targets$Target)
design <- model.matrix(~ 0 + factor(targets$Target, levels = lev))  # sample matrix
colnames(design) <- lev

# Comparison between groups
cont.wt <- makeContrasts(disease - normal, levels = design)
fit <- lmFit(exp, design)
fit2 <- contrasts.fit(fit, cont.wt)
fit2 <- eBayes(fit2)
tT <- topTable(fit2, adjust = "fdr", n = Inf)
tT <- subset(tT, select = c("adj.P.Val", "P.Value", "logFC"))
colnames(tT) <- c("FDR", "P.Value", "logFC")

### 5. Prepare gene set for GSEA
GSEAdata <- data.frame(gene = rownames(tT), logFC = tT$logFC)
geneList <- GSEAdata$logFC  # second column is logFC
names(geneList) <- GSEAdata$gene
geneList <- sort(geneList, decreasing = TRUE)  # sort by log2FC descending

### 6. GSEA enrichment
gene_enrichment_w <- data.frame()

for (n in 1:length(table(gmt$ont))) {
  gmt1 <- gmt[which(gmt$ont == unique(gmt$ont)[n]), ]  # extract gene set for the nth pathway
  if (nrow(gmt1) < 2) {
    next  # skip if less than 2 genes
  }
  gmts <- gmt1
  for (m in 1:nrow(gmt1)) {
    if (nrow(gmt1) > 1) {  # ensure at least two rows
      gmt2 <- data.frame(ont = m, gene = gmt1[-m, 2])
      gmts <- rbind(gmts, gmt2)
    }
  }
  set.seed(1234)
  GSEAresult2 <- GSEA(geneList, TERM2GENE = gmts, pvalueCutoff = 1)
  GSEAresult2 <- as.data.frame(GSEAresult2)
  if (nrow(GSEAresult2) == 0) {
    next
  }
  ES <- GSEAresult2$enrichmentScore[which(GSEAresult2$ID == unique(gmt$ont)[n])]
  if (length(ES) == 0) {
    next
  }
  GSEAresult2 <- GSEAresult2[-which(GSEAresult2$ID == unique(gmt$ont)[n]), ]
  GSEAresult2$ID <- as.numeric(GSEAresult2$ID)
  GSEAresult2 <- GSEAresult2[order(GSEAresult2$ID), ]
  coef <- ES / GSEAresult2$enrichmentScore
  w <- data.frame("gene" = gmt1$gene, "w" = coef)
  gene_enrichment_w <- rbind(gene_enrichment_w, w)
}

colnames(gene_enrichment_w) <- c("SYMBOL", "w")
imm_gene$ENSG <- imm_gene$SYMBOL
imm_p_w <- merge(imm_gene, gene_enrichment_w, by = "SYMBOL")  # includes cases where one gene may correspond to multiple pathways
keys <- colnames(imm_p_w)[!grepl("w", colnames(imm_p_w))]
imm_p_w <- as.data.table(imm_p_w)
imm_p_w <- imm_p_w[, list(w = mean(w)), keys]  # take mean for multiple
imm_p_w <- as.data.frame(imm_p_w)
colnames(imm_p_w)[c(1, 2, 5)] <- c("SYMBOL", "ENSEMBL", "w")

write.table(imm_p_w, "imm_pagerank_and_enrichment.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)
rm(coef, cont.wt, design, ES, fit, fit2, gene_enrichment_w, geneList,
   gmt1, gmt2, gmts, group, GSEAdata, GSEAresult2, keys, lev, m, n, i, tT,
   w, targets, index, gmt)

## 7. Filter IMM genes
lnc2imm2 <- lnc2imm2[which(lnc2imm2$imm %in% imm_p_w$ENSEMBL), ]
ssg2imm2 <- ssg2imm2[which(ssg2imm2$imm %in% imm_p_w$ENSEMBL), ]
lnc2ssg2 <- lnc2ssg2[which(lnc2ssg2$lnc %in% lnc2imm2$lnc), ]

write.table(lnc2ssg2, "lnc2ssg_result.txt",
            quote = FALSE, row.names = FALSE, sep = "\t")
write.table(ssg2imm2, "ssg2imm_result.txt",
            quote = FALSE, row.names = FALSE, sep = "\t")
write.table(lnc2imm2, "lnc2imm_result.txt",
            quote = FALSE, row.names = FALSE, sep = "\t")
write.table(purity_score, "purity_score.txt",
            sep = "\t", quote = FALSE)
save(lnc_exp, "lnc_exp.RData")
save(SSG_exp, "SSG_exp.RData")
save(IMM_exp, "IMM_exp.RData")

##########################
SYMBOL_ENSEMBL <- read.table("gene table.txt",
                             sep = "\t", header = TRUE)
lnc2ssg2 <- read.table("lnc2ssg_result.txt", header = TRUE)
lnc2imm2 <- read.table("lnc2imm_result.txt", header = TRUE)
ssg2imm2 <- read.table("ssg2imm_result.txt", header = TRUE)
purity_score <- read.table("purity_score.txt", sep = "\t")
load("SSG_exp.RData")
load("lnc_exp.RData")
load("IMM_exp.RData")
imm_p_w <- read.table("imm_pagerank_and_enrichment.txt",
                      sep = "\t", header = TRUE)

#### 11. Score lncRNA-SSG interaction pairs, WPC score ####
gene_pair <- lnc2ssg2[, c(1, 2)]

### 1. Calculate partial correlation for lnc2imm and ssg2imm
PCOR_LI <- function(x) {
  PCC <- pcor.test(as.numeric(lnc_exp[x[1], ]),
                   as.numeric(IMM_exp[x[2], ]),
                   as.numeric(purity_score$TumorPurity),
                   method = "pearson")
  r <- c(estimate = PCC$estimate, p = PCC$p.value)
  return(r)
}
PCC_LM <- apply(lnc2imm2, 1, PCOR_LI)
PCC_LM <- as.data.frame(t(PCC_LM))
PCC_LM <- cbind(lnc2imm2[, c(1, 2)], PCC_LM)

# Ensure p-value column is numeric
PCC_LM$p <- as.numeric(PCC_LM$p)
# Find rows where p-value is 0
zero_p_indices <- which(PCC_LM$p == 0)
# If there are zero p-values
if (length(zero_p_indices) > 0) {
  # Find non-zero p-values
  non_zero_p_values <- PCC_LM$p[PCC_LM$p != 0]
  # Find the minimum non-zero p-value
  if (length(non_zero_p_values) > 0) {
    min_non_zero_p <- min(non_zero_p_values)
    # Assign minimum non-zero p-value to all zero entries
    PCC_LM$p[zero_p_indices] <- min_non_zero_p
  }
}
rm(PCOR_LI)

###
PCOR_SI <- function(x) {
  PCC <- pcor.test(as.numeric(SSG_exp[x[1], ]),
                   as.numeric(IMM_exp[x[2], ]),
                   as.numeric(purity_score$TumorPurity),
                   method = "pearson")
  r <- c(estimate = PCC$estimate, p = PCC$p.value)
  return(r)
}
PCC_IM <- apply(ssg2imm2, 1, PCOR_SI)
PCC_IM <- as.data.frame(t(PCC_IM))
PCC_IM <- cbind(ssg2imm2[, c(1, 2)], PCC_IM)

# Ensure p-value column is numeric
PCC_IM$p <- as.numeric(PCC_IM$p)
# Find rows where p-value is 0
zero_p_indices <- which(PCC_IM$p == 0)
# If there are zero p-values
if (length(zero_p_indices) > 0) {
  # Find non-zero p-values
  non_zero_p_values <- PCC_IM$p[PCC_IM$p != 0]
  # Find the minimum non-zero p-value
  if (length(non_zero_p_values) > 0) {
    min_non_zero_p <- min(non_zero_p_values)
    # Assign minimum non-zero p-value to all zero entries
    PCC_IM$p[zero_p_indices] <- min_non_zero_p
  }
}
rm(PCOR_SI)

### 2. Count number of immune genes for each lnc-ssg pair
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
rm(imm_table)
gc()

WPC_calculate <- function(x) {  # x is lnc gene
  # 1. Find all immune genes associated with the lnc gene
  imm <- lnc2imm2[which(lnc2imm2$lnc %in% x[1]), 2]
  weight <- imm_p_w[which(imm_p_w$ENSEMBL %in% imm), "w"]
  sums <- 0
  for (i in seq_along(imm)) {
    cor1 <- PCC_LM[which(PCC_LM$lnc %in% x[1] & PCC_LM$imm %in% imm[i]), 3]
    if (length(cor1) == 0) next
    p1 <- PCC_LM[which(PCC_LM$lnc %in% x[1] & PCC_LM$imm %in% imm[i]), 4]
    # if(length(p1)==0){next}
    cor2 <- PCC_IM[which(PCC_IM$imm %in% imm[i] & PCC_IM$ssg %in% x[2]), 3]
    if (length(cor2) == 0) next
    p2 <- PCC_IM[which(PCC_IM$imm %in% imm[i] & PCC_IM$ssg %in% x[2]), 4]
    
    sums <- sums + weight[i] * ((-log10(p1)) * sign(cor1) + (-log10(p2)) * sign(cor2))
  }
  return(sums)
}
system.time({
  WPC <- apply(gene_pair, 1, WPC_calculate)
})
WPC <- as.data.frame(WPC)  # scores for lnc-ssg interaction pairs
WPC <- cbind(gene_pair, WPC)
write.csv(WPC, "WPC score.csv", row.names = FALSE)

## permutation
library(ppcor)
library(progress)
WPC <- read.csv("WPC score.csv")
gene_pair <- lnc2ssg2[, c(1, 2)]

#### 12. Permutation ####
permutation_LS <- lnc2ssg2[, c(1, 2)]
permutation_LS[, "num"] <- 0

################################################################################
### 1. Function 1
PCOR_LI <- function(x) {
  PCC <- pcor.test(as.numeric(lnc_exp[x[1], ]),
                   as.numeric(IMM_exp2[x[2], ]),
                   as.numeric(purity_score$TumorPurity),
                   method = "pearson")
  r <- c(estimate = PCC$estimate, p = PCC$p.value)
  return(r)
}

################################################################################
### 2. Function 2
PCOR_SI <- function(x) {
  PCC <- pcor.test(as.numeric(SSG_exp[x[1], ]),
                   as.numeric(IMM_exp2[x[2], ]),
                   as.numeric(purity_score$TumorPurity),
                   method = "pearson")
  r <- c(estimate = PCC$estimate, p = PCC$p.value)
  return(r)
}

################################################################################
### 3. Function 3
WPC_calculate <- function(x) {   # x is lnc gene
  # 1. First, find how many imm genes are associated with each lnc gene
  imm <- lnc2imm2[which(lnc2imm2$lnc %in% x[1]), 2]
  weight <- imm_p_w[which(imm_p_w$ENSEMBL %in% imm), "w"]
  sums <- 0
  for (i in seq_along(imm)) {
    cor1 <- PCC_LM2[which(PCC_LM2$lnc %in% x[1] & PCC_LM2$imm %in% imm[i]), 3]
    if (length(cor1) == 0) next
    p1 <- PCC_LM2[which(PCC_LM2$lnc %in% x[1] & PCC_LM2$imm %in% imm[i]), 4]
    # if (length(p1) == 0) next
    cor2 <- PCC_IM2[which(PCC_IM2$imm %in% imm[i] & PCC_IM2$ssg %in% x[2]), 3]
    if (length(cor2) == 0) next
    p2 <- PCC_IM2[which(PCC_IM2$imm %in% imm[i] & PCC_IM2$ssg %in% x[2]), 4]
    # if (length(p2) == 0) next
    sums <- sums + weight[i] * ((-log10(p1)) * sign(cor1) + (-log10(p2)) * sign(cor2))
  }
  return(sums)
}

################################################################################
pb <- progress_bar$new(total = 1000)
system.time({
  for (i in 1:1000) {   # 1000 times
    pb$tick()
    # 1. Randomly permute sample labels of IMM genes
    col <- colnames(IMM_exp)
    IMM_exp2 <- IMM_exp[, sample(col, length(col))]
    colnames(IMM_exp2) <- col
    
    # 2. Calculate partial correlation coefficients for lnc2imm and ssg2imm
    PCC_LM2 <- apply(lnc2imm2, 1, PCOR_LI)
    PCC_LM2 <- as.data.frame(t(PCC_LM2))
    PCC_LM2 <- cbind(lnc2imm2[, c(1, 2)], PCC_LM2)
    
    PCC_IM2 <- apply(ssg2imm2, 1, PCOR_SI)
    PCC_IM2 <- as.data.frame(t(PCC_IM2))
    PCC_IM2 <- cbind(ssg2imm2[, c(1, 2)], PCC_IM2)
    
    WPC2 <- apply(gene_pair, 1, WPC_calculate)
    WPC2 <- as.data.frame(WPC2)   # Scores for lnc-ssg pairs
    WPC2 <- cbind(gene_pair, WPC2)
    
    # 3. Compare values to calculate p-value
    for (j in seq_len(nrow(WPC))) {
      if (WPC[j, 3] < WPC2[j, 3]) {
        permutation_LS[j, 3] <- permutation_LS[j, 3] + 1
      } else {
        permutation_LS[j, 3] <- permutation_LS[j, 3] + 0
      }
    }
    ### Permutation results are stored in permutation_LS
  }
})   # Timing finished

write.csv(permutation_LS,
          "permutation_result.csv",
          row.names = FALSE, quote = FALSE)