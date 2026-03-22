# *********************************
# * Build co-expression network after Pearson filtering *
# *********************************
lnc2ssg <- lnc2ssg2[, -c(4)]
ssg2imm <- ssg2imm2[, -c(4)]
lnc2imm <- lnc2imm2[, -c(4)]
net <- data.frame()
net <- rbind(net, lnc2ssg)
colnames(net) <- colnames(ssg2imm)
net <- rbind(net, ssg2imm)
colnames(net) <- colnames(lnc2imm)
net <- rbind(net, lnc2imm)  # net is the co-expression network
colnames(net) <- c("gene1", "gene2", "cor")
write.table(net, "net.txt", row.names = FALSE, quote = FALSE, sep = "\t")
remove(lnc2ssg, ssg2imm, lnc2imm, net)

# =================================
# ========== page rank =============
# =================================
ssg <- read.csv("SSG.csv")

# Read file and create graph with edge weights
G_weighted_df <- read.table("net.txt",
                            header = TRUE, stringsAsFactors = FALSE)
G_weighted <- graph_from_data_frame(d = G_weighted_df, directed = FALSE)

# Set edge weight attribute
E(G_weighted)$weight <- G_weighted_df$cor  

# Check if edges and weights are correctly set
print(E(G_weighted)$weight)
print(length(E(G_weighted)$weight))  # Ensure length matches number of edges

# Initialize propagate_input vector
pIDs <- ssg$common_elements
network_nodes <- V(G_weighted)$name
propagate_input <- setNames(rep(0, length(network_nodes)), network_nodes)
propagate_input[names(propagate_input) %in% pIDs] <- 1

# Compute weighted personalized PageRank
weighted_personalized_pagerank <- page_rank(
  G_weighted,
  personalized = propagate_input,
  weights = E(G_weighted)$weight,  # Use edge weights
  damping = 0.85
)$vector

# Create data frame and save results
df_metrics <- data.frame(ENSG = names(weighted_personalized_pagerank),
                         pagerank = weighted_personalized_pagerank)
write.table(df_metrics, "result.txt",
            sep = "\t", row.names = FALSE)

#####################
#### File preparation ###
#####################
SYMBOL_ENSEMBL <- read.table("gene table.txt",
                             sep = "\t", header = TRUE)
lncRNA <- SYMBOL_ENSEMBL[SYMBOL_ENSEMBL$gene_biotype == "lncRNA", ]
lncRNA <- lncRNA[lncRNA$hgnc_symbol != "", ]
colnames(lncRNA) <- c("ENSG", "SYMBOL", "TYPE")

SSG <- read.csv("SSG.csv")
IMM <- read_excel("imm_genelist.xlsx")

lnc2ssg2 <- read.table("lnc2ssg_result.txt",
                       header = TRUE)
lnc2imm2 <- read.table("lnc2imm_result.txt",
                       header = TRUE)
ssg2imm2 <- read.table("ssg2imm_result.txt",
                       header = TRUE)

load("SSG_exp.RData")
load("lnc_exp.RData")
load("IMM_exp.RData")

# ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
# * Select highly influential genes closely related to SSG *
# ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
result <- read.table("result.txt", header = TRUE)
selection <- result[which(!(result$ENSG %in% SSG$common_elements)), ]

# Extract rows with non-zero probability
node_0 <- subset(selection, pagerank > 0)

top_200 <- 200
# Sort by pagerank descending and select top 200
gene_top_200 <- head(node_0[order(-node_0$pagerank), ], top_200)

# Classify selected genes
gene_top_200[which(gene_top_200$ENSG %in% lncRNA$SYMBOL), "class"] <- "lnc"
gene_top_200[which(gene_top_200$ENSG %in% IMM$Symbol), "class"] <- "IMM"

table(gene_top_200$class)

imm_gene <- gene_top_200[which(gene_top_200$class == "IMM"), ]
lnc_gene <- gene_top_200[which(gene_top_200$class == "lnc"), ]
SSG_gene <- result[which(result$ENSG %in% SSG$common_elements), ]

rm(result, selection)
