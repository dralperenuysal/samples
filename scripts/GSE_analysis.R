library(GEOquery)
library(glue)
library(dplyr)
library(DESeq2)
library(annotables)
library(tidyr)

# Define the data
GEO <- "GSE------"
output <- glue("output/rnaseq/{GEO}/")

# Get the sample info
gse <- getGEO(GEO = GEO, GSEMatrix = T, AnnotGPL = T)[[1]]
sample_info <- pData(phenoData(gse))

# Get the supplementary file if it was not downloaded
getGEOSuppFiles(GEO = GEO, makeDirectory = T, baseDir = glue("resource/rnaseq/{GEO}"))

# Manipulate the sample_info for further differentially expressed gene analysis
sample_info <- sample_info %>%
  select(c(1, 2, 47)) # Select the required columns

sample_info <- sample_info %>%
  rename("condition:ch1" = "state") %>%  # Rename the column
  mutate("state" = case_when(
    state == "PC" ~ "control",
    TRUE ~ state
  )) %>%
  mutate(title = gsub("RNA-seq_", "", title))

# Read the counts if already downladed
counts <- read.delim(file = "")
rownames(counts) <- counts$id # Id column should be rownames for DEG analysis
counts <- counts %>% select(-id)

# Clean the fused the columns
counts <- counts %>%
  rename("A27.21.A27.22" = "A27.21") %>%
  rename("A26.18.A26.19" = "A26.18")

# Fix the index and columns (AI idea, worked well)
index <- match(sample_info$title, colnames(counts))
index <- index[!is.na(index)]
counts <- counts %>% select(c(index))
rownames(sample_info) <- sample_info$title

# Find rownames in sample_info that do not exist in counts columns
rows_to_delete <- rownames(sample_info)[!(rownames(sample_info) %in% colnames(counts))]

# Remove rows from sample_info
sample_info <- sample_info[!(sample_info$title %in% rows_to_delete), ]

# Checking the design
all(colnames(counts) %in% rownames(sample_info)) # Should be TRUE
all(colnames(counts) == rownames(sample_info)) # Should be TRUE

# Differentially Expressed Gene Analysis with DESeq2
counts <- round(counts) # This is sometimes necessary to avoid errors.

# Make the data numeric.
original_rownames <- rownames(counts) # Keep original rownames
counts_numeric <- as.matrix(sapply(counts, as.numeric))
rownames(counts_numeric) <- original_rownames
counts <- counts_numeric

# DEG analysis
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_info,
                              design = ~ state)

# Drop the low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Relevel the factors
dds$state <- relevel(dds$state, ref = "control")

# Run Deseq2
dds <- DESeq(dds)
res <- results(dds)

# Get the summary
summary(res)

# Make the results dataframe to see important values including logFC and adjusted p value.
df <- as.data.frame(res)
df <- df %>%
  filter(padj <0.05) %>%
  filter(log2FoldChange > 1 | log2FoldChange < -1) # No significant results were found for GSE53697

# Map the entrez IDs of the DEGs
# For this purpose, use the DAVID.
# Check if the directory exists, if not, create it
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

# Get the gene symbols with annotables package
grch38 <- grch38

# Make the gene column for gene symbol conversion
df['ensgene'] <- rownames(df)

# map with annotables
map <- grch38 %>%
  filter(ensgene %in% df$ensgene)

# Save the mapped dataframe
mapped_df <- left_join(df, map, by = "ensgene")
write.table(mapped_df, file = paste0(output, GEO, "_mapped.tsv"), sep = '\t', row.names = F, col.names = T) # Save the data

# Significant DEGs were found.


