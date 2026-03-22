# ▷▷▷▷▷▷▷▷▷▷▷▷▷▷
# ▷ 加载多个包   ▷
# ▷▷▷▷▷▷▷▷▷▷▷▷▷▷
lapply(c("dplyr", "readr", "clusterProfiler", "org.Hs.eg.db", "stringr", 
         "edgeR", "limma", "purrr", "ppcor", "progress", "utils","estimate",
         "data.table", "tidyr", "readxl","igraph","splitstackshape"), library, character.only = TRUE)

### Read expression profile
exp <- read.table(file_to_read, row.names = 1, check.names = FALSE, header = TRUE, sep = "\t")

### Remove genes with zero expression in >70% of samples
## Calculate number of zeros per row
fc_zero <- function(f) {
  zero_num <- sum(f == 0)
  return(zero_num)
}
zero <- as.data.frame(apply(exp, 1, fc_zero))
colnames(zero) <- "num"
zero[, "ratio"] <- zero$num / length(colnames(exp))
exp <- exp[which(zero$ratio < 0.7), ]

# +-------------------------------+
# | Separate normal and tumor samples, normal first, tumor after |
# +-------------------------------+
index <- which(substr(colnames(exp), start = 14, stop = 15) > 10)
normal_exp <- exp[, index]
exp <- exp[, -index]
# Sort samples
nt <- data.frame(normal_exp, exp, check.names = FALSE)
exp <- nt
remove(nt, fc_zero, index, normal_exp, zero)

# Save file
save(exp, "after_group_exp.RData")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# > Differential Expression Analysis >
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
load("after_group_exp.RData")
index <- which(substr(colnames(exp), start = 14, stop = 15) > 10)
group <- c(rep("normal", length(index)), rep("tumor", length(colnames(exp)) - length(index)))
targets <- cbind(colnames(exp), group)
targets <- as.data.frame(targets)
colnames(targets) <- c("Sample", "Group")
lev <- unique(targets$Group)
design <- model.matrix(~ 0 + factor(targets$Group, levels = lev))  # sample matrix
colnames(design) <- lev  # rename columns to level names

# Pairwise comparison, use topTable to obtain differentially expressed genes
cont.wt <- makeContrasts(tumor - normal, levels = design)
fit <- lmFit(exp, design)
fit2 <- contrasts.fit(fit, cont.wt)
fit2 <- eBayes(fit2)
tT <- topTable(fit2, adjust = "fdr", n = Inf)
tT <- subset(tT, select = c("adj.P.Val", "P.Value", "logFC"))
colnames(tT) <- c("FDR", "P.Value", "logFC")
tT <- tT[which(tT$FDR < 0.05), ]  
exp <- exp[rownames(tT), ]
exp <- exp[, -c(1:length(index))]  # remove normal samples
save(exp, file = "tumor_exp.RData") # Save expression profile of differentially expressed genes only
remove(cont.wt, design, fit, fit2, group, lev, targets, tT, index)

# ♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦
# ♦ Split expression profiles for three gene types ♦
# ♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦♦
### 1. Extract genes from expression profile
load("tumor_exp.RData")
name <- as.data.frame(rownames(exp))
colnames(name) <- "ENSG"

### 2. Read annotation files
SYMBOL_ENSEMBL <- read.table("gene table.txt", sep = "\t", header = TRUE)

##
lncRNA <- SYMBOL_ENSEMBL[SYMBOL_ENSEMBL$gene_biotype == "lncRNA", ]
lncRNA <- lncRNA[lncRNA$hgnc_symbol != "", ]
colnames(lncRNA) <- c("ENSG", "SYMBOL", "TYPE")  # extract lncRNA
##
SSG <- read.csv("SSG.csv")
IMM <- read_excel("imm_genelist.xlsx")

### 3. Partition into three gene types
lncRNA_list <- name[which(name$ENSG %in% lncRNA$SYMBOL), 1]  # number of differentially expressed lncRNAs in the expression profile
SSG_list <- name[which(name$ENSG %in% SSG$common_elements), 1]  # number of differentially expressed SSGs in the expression profile
IMM_list <- name[which(name$ENSG %in% IMM$Symbol), 1]  # number of differentially expressed immune genes in the expression profile

### 4. Separate expression profiles for the three types
name[which(name$ENSG %in% lncRNA_list), "class"] <- "lnc"
name[which(name$ENSG %in% IMM_list), "class"] <- "IMM"
name[which(name$ENSG %in% SSG_list), "class"] <- "SSG"
table(name$class)

lnc_exp <- exp[lncRNA_list, ]
SSG_exp <- exp[SSG_list, ]
IMM_exp <- exp[IMM_list, ]

# Transpose expression profiles so that rows are samples and columns are genes (for subsequent map2)
SSG_exp <- as.data.frame(t(SSG_exp))
lnc_exp <- as.data.frame(t(lnc_exp))
IMM_exp <- as.data.frame(t(IMM_exp))

# Save the three expression profiles
save(lnc_exp, file = "lnc_exp.RData")  # columns: genes, rows: samples
save(SSG_exp, file = "SSG_exp.RData")  # columns: genes, rows: samples
save(IMM_exp, file = "IMM_exp.RData")  # columns: genes, rows: samples

# .................................
# . Calculate Pearson partial correlation .
# .................................

### 1. Compute lnc2ssg
k <- 1
lnc2ssg <- matrix(ncol = 2, nrow = length(SSG_list) * length(lncRNA_list))
for (i in 1:length(lncRNA_list)) {
  for (j in 1:length(SSG_list)) {
    lnc2ssg[k, ] <- c(lncRNA_list[i], SSG_list[j])
    k <- k + 1
  }
}
lnc2ssg <- as.data.frame(lnc2ssg)
lnc2ssg_data1 <- lnc_exp[, lnc2ssg[, 1]]
lnc2ssg_data2 <- SSG_exp[, lnc2ssg[, 2]]

# Compute Pearson correlation coefficient
PC_lnc2ssg <- function(x, y) {
  PC <- cor.test(as.numeric(x), as.numeric(y), method = "pearson")
  PC_estimate <- PC$estimate
  PC_pvalue <- PC$p.value
  r <- c(R = PC_estimate, p = PC_pvalue)
  return(r)
}
system.time({
  lnc2ssg_PCS <- map2_df(lnc2ssg_data1, lnc2ssg_data2, PC_lnc2ssg)
})
lnc2ssg <- cbind(lnc2ssg, lnc2ssg_PCS)
colnames(lnc2ssg) <- c("lnc", "ssg", "R", "p")
lnc2ssg2 <- lnc2ssg
lnc2ssg2 <- lnc2ssg2[which(lnc2ssg2$p < 0.05), ]
lnc2ssg2 <- lnc2ssg2[which(abs(lnc2ssg2$R) > 0.5), ]  #
write.table(lnc2ssg2, "lnc2ssg_result.txt", quote = FALSE, row.names = FALSE, sep = "\t")
remove(PC_lnc2ssg, lnc2ssg_PCS, lnc2ssg_data1, lnc2ssg_data2, lnc2ssg, k)

### 2. Compute ssg2imm
k <- 1
ssg2imm <- matrix(ncol = 2, nrow = length(SSG_list) * length(IMM_list))
for (i in 1:length(SSG_list)) {
  for (j in 1:length(IMM_list)) {
    ssg2imm[k, ] <- c(SSG_list[i], IMM_list[j])
    k <- k + 1
  }
}
ssg2imm <- as.data.frame(ssg2imm)
ssg2imm_data1 <- SSG_exp[, ssg2imm[, 1]]
ssg2imm_data2 <- IMM_exp[, ssg2imm[, 2]]

# Compute Pearson correlation coefficient
PC_ssg2imm <- function(x, y) {
  PC <- cor.test(as.numeric(x), as.numeric(y), method = "pearson")
  PC_estimate <- PC$estimate
  PC_pvalue <- PC$p.value
  r <- c(R = PC_estimate, p = PC_pvalue)
  return(r)
}

system.time({
  ssg2imm_PCS <- map2_df(ssg2imm_data1, ssg2imm_data2, PC_ssg2imm)
})
ssg2imm <- cbind(ssg2imm, ssg2imm_PCS)
colnames(ssg2imm) <- c("ssg", "imm", "R", "p")
ssg2imm2 <- ssg2imm
ssg2imm2 <- ssg2imm2[which(ssg2imm2$p < 0.05), ]
ssg2imm2 <- ssg2imm2[which(abs(ssg2imm2$R) > 0.5), ]  #
write.table(ssg2imm2, "ssg2imm_result.txt", quote = FALSE, row.names = FALSE, sep = "\t")
remove(PC_ssg2imm, ssg2imm_PCS, ssg2imm_data1, ssg2imm_data2, ssg2imm, k)

### 3. Compute lnc2imm
k <- 1
lnc2imm <- matrix(ncol = 2, nrow = length(lncRNA_list) * length(IMM_list))
for (i in 1:length(lncRNA_list)) {
  for (j in 1:length(IMM_list)) {
    lnc2imm[k, ] <- c(lncRNA_list[i], IMM_list[j])
    k <- k + 1
  }
}
lnc2imm <- as.data.frame(lnc2imm)
lnc2imm_data1 <- lnc_exp[, lnc2imm[, 1]]
lnc2imm_data2 <- IMM_exp[, lnc2imm[, 2]]

# Compute Pearson correlation
PC_lnc2imm <- function(x, y) {
  PC <- cor.test(as.numeric(x), as.numeric(y), method = "pearson")
  PC_estimate <- PC$estimate
  PC_pvalue <- PC$p.value
  r <- c(R = PC_estimate, p = PC_pvalue)
  return(r)
}
system.time({
  lnc2imm_PCS <- map2_df(lnc2imm_data1, lnc2imm_data2, PC_lnc2imm)
})
lnc2imm <- cbind(lnc2imm, lnc2imm_PCS)
colnames(lnc2imm) <- c("lnc", "imm", "R", "p")
lnc2imm2 <- lnc2imm
lnc2imm2 <- lnc2imm2[which(lnc2imm2$p < 0.05), ]
lnc2imm2 <- lnc2imm2[which(abs(lnc2imm2$R) > 0.5), ]  #
write.table(lnc2imm2, "lnc2imm_result.txt", quote = FALSE, row.names = FALSE, sep = "\t")
remove(PC_lnc2imm, lnc2imm_PCS, lnc2imm_data1, lnc2imm_data2, lnc2imm, k, i, j)
