# ----------- Part 1: ROCK1/2 expression analysis on normal and tumor thyroid samples  ----------
library(dplyr)
library(ggplot2)
library(pROC)
library(stringr)
library(tidyr)
library(ggpattern)
library(ggridges)
library(ggrepel)
library(jsonlite)

# Load metadata from JSON file
raw.metadata <- fromJSON("/Volumes/Extreme SSD/thyroid-ml/metadata.repository.2024-07-16.json")

# Filter metadata for RNA-Seq experimental strategy
raw.metadata <- filter(raw.metadata, experimental_strategy == "RNA-Seq")

# Extract case_id information from JSON files
case_ids <- character(nrow(raw.metadata))
for (i in 1:nrow(raw.metadata)) {
        case_ids[i] <- raw.metadata[[3]][[i]]$case_id
}

# Convert case_ids to a character vector and add to metadata
case_ids <- unlist(case_ids)
raw.metadata$case_id <- case_ids

# Keep only the 4th and last column of metadata
raw.metadata <- raw.metadata[,c(4,ncol(raw.metadata))]

# Load sample metadata from TSV file
sample.metadata <- read.delim("/Volumes/Extreme SSD/thyroid-ml/biospecimen.cart.2024-07-16/sample.tsv", sep = "\t")
THCA.metadata <- subset(sample.metadata, project_id == "TCGA-THCA")

# Filter sample metadata for specific sample types
sample.metadata <- filter(THCA.metadata, sample_type %in% c("Solid Tissue Normal", "Primary Tumor"))

# Merge sample metadata with raw metadata on case_id
sample.metadata <- merge(sample.metadata, raw.metadata, by = 'case_id')

# Remove duplicate file names, keeping all other information
sample.metadata <- distinct(sample.metadata, file_name, .keep_all = TRUE)

# Separate tumor and control tissues
control.metadata <- filter(sample.metadata, sample_type == "Solid Tissue Normal")
tumor.metadata <- filter(sample.metadata, sample_type == "Primary Tumor")

# Download RNA-seq data for control samples
file_names <- control.metadata$file_name
control_data <- data.frame()
path <- "/Volumes/Extreme SSD/thyroid-ml/rna-seq/"
for (i in file_names) {
        file_path <- paste0(path, i)
        if (file.exists(file_path)) {
                file <- read.table(file_path, header = TRUE, sep = "\t")
                file <- filter(file, gene_name %in% c("ROCK1", "ROCK2", "MIR222", "MIR221"))
                file <- file[, c(2, 7)]
                colnames(file)[2] <- i
                if (ncol(control_data) == 0) {
                        control_data <- file
                } else {
                        control_data <- bind_cols(control_data, file[, 2, drop = FALSE])
                }
                print(paste(i, "was downloaded"))
        } else {
                print(paste(i, "was not found"))
        }
}

# Transpose control data and set column names
control_data <- t(control_data[,2:ncol(control_data)])
colnames(control_data) <- c("ROCK1", "ROCK2", "MIR222", "MIR221")

# Convert control data to data frame and set all columns to numeric
control_data <- as.data.frame(control_data)
control_data <- control_data %>% mutate_all(as.numeric)

# Download RNA-seq data for tumor samples
file_names <- tumor.metadata$file_name
tumor_data <- data.frame()
path <- "/Volumes/Extreme SSD/thyroid-ml/rna-seq/"
for (i in file_names) {
        file_path <- paste0(path, i)
        if (file.exists(file_path)) {
                file <- read.table(file_path, header = TRUE, sep = "\t")
                file <- filter(file, gene_name %in% c("ROCK1", "ROCK2", "MIR222", "MIR221"))
                file <- file[, c(2, 7)]
                colnames(file)[2] <- i
                if (ncol(tumor_data) == 0) {
                        tumor_data <- file
                } else {
                        tumor_data <- bind_cols(tumor_data, file[, 2, drop = FALSE])
                }
                print(paste(i, "was downloaded"))
        } else {
                print(paste(i, "was not found"))
        }
}

# Transpose tumor data and set column names
tumor_data <- t(tumor_data[,2:ncol(tumor_data)])
colnames(tumor_data) <- c("ROCK1", "ROCK2", "MIR222", "MIR221")
tumor_data <- as.data.frame(tumor_data)

# Calculate mean ROCK1 and ROCK2 expression for control and tumor samples
t.test(tumor_data$ROCK1, control_data$ROCK1)
t.test(tumor_data$ROCK2, control_data$ROCK2)
t.test(tumor_data$MIR222, control_data$MIR222)
t.test(tumor_data$MIR221, control_data$MIR221)

"From these initial results, we observe that after the ontogenic transformation in thyroid tumors, there is an increase in the expression of ROCK1 and ROCK2. In the next step, I will check the correlation of ROCK1 and ROCK2 expression with the project's target microRNAs: MIR221 and MIR222."


# ---------------- Part 2: Correlation analysis between ROCK1/2 and MIR221/222 on thyroid tumor samples ----------------

###### Density plot

# Assuming tumor_data is your original dataframe
tumor_data_long <- tumor_data %>%
        pivot_longer(cols = c(MIR222, MIR221, ROCK1, ROCK2), names_to = "Gene", values_to = "Expression")

# Create the density plots side by side
ggplot(tumor_data_long, aes(x = Expression)) +
        geom_density(color = "black", size = 1) +  # Add density line
        labs(title = NULL, x = NULL, y = NULL) +
        theme(
                panel.background = element_rect(size = 0.5, linetype = "solid"),
                panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"), 
                panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "white"),
                axis.text.x = element_text(size = 12),  # Increase x-axis text size
                axis.text.y = element_text(size = 12)   # Increase y-axis text size
        ) +
        facet_wrap(~ Gene, scales = "free_x") +  # Facet by Gene with free x scales
        coord_flip()  # Flip coordinates

#### Scatterplot
# Reorganize the data to long format
library(ggpubr)
tumor_data_long <- tumor_data %>%
        pivot_longer(cols = c(MIR222, MIR221), names_to = "MIR", values_to = "MIR_Expression") %>%
        pivot_longer(cols = c(ROCK1, ROCK2), names_to = "ROCK", values_to = "ROCK_Expression")

# Create the scatter plots side by side with correlation values
ggplot(tumor_data_long, aes(x = ROCK_Expression, y = MIR_Expression)) +
        geom_point(size = 2) +  # Add points to the plot
        geom_smooth(method = "lm", color = "blue") +  # Add linear regression line
        labs(x = NULL, y = NULL) +  # Remove axis titles
        theme(
                panel.background = element_rect(size = 0.5, linetype = "solid"),
                panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"), 
                panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "white"),
                axis.text.x = element_text(size = 12),  # Increase x-axis text size
                axis.text.y = element_text(size = 12)   # Increase y-axis text size
        ) +
        stat_cor(method = "pearson", label.x = 45, label.y = 250, size = 4) +  # Add Pearson correlation values
        facet_grid(MIR ~ ROCK)  # Create facets for each combination of MIR and ROCK

"For ROCK2 we find a very low correlation. This may suggest that some samples have a positive MIR/ROCK relationship, 
but others may have a negative relationship. I will try to separate these samples."

# ----------- Part 3: Clustering samples based on the expression of ROCK/MIR222 MIR221 ----------
# Define the number of clusters
# First, I will try to remove some outliers
tumor_samples <- tumor_data
tumor_data <- tumor_samples

tumor_data <- tumor_data %>%
        filter(MIR221 <= 40 & MIR222 <= 6)

set.seed(123)  # For reproducibility
num_clusters <- 6
data_selected <- tumor_data %>% select(MIR222, MIR221, ROCK2)
data_scaled <- scale(data_selected)
# Perform K-means
kmeans_result <- kmeans(data_scaled, centers = num_clusters)

# Add cluster labels to the original dataset
tumor_data$cluster <- kmeans_result$cluster

ggplot(tumor_data, aes(x = ROCK2, y = MIR221)) +
        geom_point(size = 2) +  # Add points to the plot
        geom_smooth(method = "lm", color = "blue") +  # Add linear regression line
        labs(x = NULL, y = NULL) +
        facet_wrap(~cluster)

# After identifying the groups, I will separate different dataframes containing the filenames of each sample. 
# After defining the dataframes with the 8 groups, I will download the RNA-seq data again and compare the DEG between the groups.

# miR-222
group_1 = tumor_data %>% # mirHIGH/rockHIGH - miR222
        filter(cluster %in% c(1,2,4,3))
group_1$group <- "1"
group_1 <- group_1 %>%
        filter(MIR221 > 15 & ROCK2 > 15)
group_1 <- group_1[,c(2,4,6)]

group_2 = tumor_data %>% # mirLOW/rockLOW - miR222
        filter(cluster %in% c(6))
group_2$group <- "2"
group_2 <- group_2 %>%
        filter(MIR221 < 10 & ROCK2 < 10)
group_2 <- group_2[,c(2,4,6)]

group_3 = tumor_data %>% # mirHIGH/rockLOW - miR222
        filter(cluster %in% c(2,4,5))
group_3$group <- "3"
group_3 <- group_3 %>%
        filter(MIR221 > 15 & ROCK2 < 15)
group_3 <- group_3[,c(2,4,6)]

group_4 = tumor_data %>% # mirLOW/rockHIGH - miR222
        filter(cluster %in% c(3))
group_4$group <- "4"
group_4 <- group_4 %>%
        filter(MIR221 < 15 & ROCK2 > 20)
group_4 <- group_4[,c(2,4,6)]


# ------------- Part 4: What are the clinical implications of ROCK2 modulation by miRs? ----------------------
# Initially, I will select the groups of patients where there is a regulation of ROCK triggered by miRs.
group_mir221 <- rbind(group_3, group_4)
group_mir221$file_name <- rownames(group_mir221)

# Read pathology metadata
pathology_metadata <- read.delim("/Volumes/Extreme SSD/thyroid-ml/clinical.cart.2024-07-16/clinical.tsv", sep = "\t")
pathology_metadata <- pathology_metadata[,c(1,29:32)] # Select relevant columns
pathology_metadata <- lapply(pathology_metadata, function(x) gsub("'--", NA, x)) # Replace '--' with NA
pathology_metadata <- as.data.frame(pathology_metadata) # Convert to data frame

# Segregate ROCK1/2 high and low groups based on expression levels
filtered.metadata <- tumor.metadata[,c(1,36,40)] # Select relevant columns from tumor metadata
filtered.metadata <- merge(filtered.metadata, pathology_metadata, by = "case_id") # Merge metadata with expression data
filtered.metadata <- distinct(filtered.metadata, file_name, .keep_all = TRUE)

group_mir221 <- merge(filtered.metadata, group_mir221, by = "file_name")

######## Count the number of cases in each group and pathology stage
table(group_mir221$group, group_mir221$ajcc_pathologic_n)

# At the end, the statistical analysis was performed on graphpad


# ------------ Part 5: What is the survival of patients based on the modulation of ROCK2 by miRs? ---------------

# Load the Bioconductor installer.
if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
BiocManager::install(version = "3.18")

# Install the main RTCGA package and its dependencies
BiocManager::install("RTCGA")
BiocManager::install("RTCGA.clinical")
BiocManager::install("RTCGA.mRNA")
install.packages("survminer")

# Load necessary libraries for survival analysis
library(survminer)
library(survival)
library(RTCGA)
infoTCGA()
library(RTCGA.clinical)

# Merge barcodes with ROCK expression data
bar_code <- sample.metadata[,c(1,3)]
bar_code <- merge(bar_code, group_mir221, by = "case_id")
bar_code <- unique(bar_code)
colnames(bar_code)[2] <- "bcr_patient_barcode"

# Create the clinical data
clin <- survivalTCGA(THCA.clinical, extract.cols="admin.disease_code")
clin <- merge(bar_code, clin, by = "bcr_patient_barcode")

# Tabulate survival data by outcome for ROCK1
xtabs(~group + patient.vital_status, data = clin) %>% addmargins()
coxph(Surv(times, patient.vital_status) ~ group, data = clin)
sfit <- survfit(Surv(times, patient.vital_status) ~ group, data = clin)

# Plot survival analysis for ROCK1
ggsurvplot(sfit, conf.int = TRUE, pval = TRUE, risk.table = TRUE, 
           palette = c("dodgerblue2", "#AFABAB"))

# No difference between groups.

# ------------- Part 6: What are the clinical implications of miRs when they fail to modulate ROCK? ----------------------
# Now I will select the patient groups where there is NO regulation of ROCK triggered by miRs. They are: 1 and 2 // 5 and 6

group_mir221 <- rbind(group_1, group_2)
group_mir221$file_name <- rownames(group_mir221)

# Read pathology metadata
pathology_metadata <- read.delim("/Volumes/Extreme SSD/thyroid-ml/clinical.cart.2024-07-16/clinical.tsv", sep = "\t")
pathology_metadata <- pathology_metadata[,c(1,29:32)] # Select relevant columns
pathology_metadata <- lapply(pathology_metadata, function(x) gsub("'--", NA, x)) # Replace '--' with NA
pathology_metadata <- as.data.frame(pathology_metadata) # Convert to data frame

# Segregate ROCK1/2 high and low groups based on expression levels
filtered.metadata <- tumor.metadata[,c(1,36,40)] # Select relevant columns from tumor metadata
filtered.metadata <- merge(filtered.metadata, pathology_metadata, by = "case_id") # Merge metadata with expression data
filtered.metadata <- distinct(filtered.metadata, file_name, .keep_all = TRUE)

group_mir221 <- merge(filtered.metadata, group_mir221, by = "file_name")

######## Count the number of cases in each group and pathology stage
table(group_mir221$group, group_mir221$ajcc_pathologic_n)


# ------ Part 7: Analysis of differentially expressed genes between groups ----------
# ------ Part 7.1: Analysis of patient groups that can OR cannot modulate the ROCK-miR axis ------
############ DESeq Analysis
library(dplyr)
library(ggplot2)
library(pROC)
library(stringr)
library(tidyr)
library(ggpattern)
library(ggridges)
library(DESeq2)
library(clusterProfiler)
library(circlize)
library(org.Hs.eg.db)

# After separating the groups, I will download the RNA-seq data again and compare ALL DEG between groups.
# The groups that can modulate are 'pos', and those that cannot are 'neg'. My groups are:
pos_mir <- rbind(group_3, group_4)
neg_mir <- rbind(group_1, group_2)

pos_mir$file_names <- rownames(pos_mir)
neg_mir$file_names <- rownames(neg_mir)

groups <- list(pos_mir = pos_mir, neg_mir = neg_mir)

path <- "/Volumes/Extreme SSD/thyroid-ml/rna-seq/"

# Loop to download the full RNA-seq data for samples corresponding to each group
for (group_name in names(groups)) {
        group <- groups[[group_name]]
        file_names <- group$file_name
        data_j <- data.frame()
        
        for (file_name in file_names) {
                file_path <- paste0(path, file_name)
                
                if (file.exists(file_path)) {
                        file <- read.table(file_path, header = TRUE, sep = "\t")
                        file <- file %>%
                                dplyr::select(gene_id, tpm_unstranded, gene_type) %>%
                                dplyr::filter(gene_type %in% c("miRNA", "protein_coding"))
                        colnames(file)[2] <- file_name
                        
                        if (ncol(data_j) == 0) {
                                data_j <- file
                        } else {
                                data_j <- cbind(data_j, file[file_name])
                        }
                        print(paste(file_name, "was downloaded"))
                } else {
                        print(paste(file_name, "was not found"))
                }
        }
        
        # Ensuring the first column is gene_id
        gene_ids <- data_j$gene_id
        data_j <- data_j[ , !(names(data_j) %in% c("gene_id", "gene_type"))]
        data_j <- cbind(gene_id = gene_ids, data_j)
        
        # Saving the data frame in the global environment
        assign(paste0(group_name, "_cts"), data_j, envir = .GlobalEnv)
}

# First, check if the column names of cts match the row names of metadata
cts <- c("pos_mir_cts", "neg_mir_cts")
metadata <- c("pos_mir", "neg_mir")
for (i in cts) {
        for (j in metadata) {
                df <- get(i)
                df_2 <- get(j)
                rownames(df) <- df$gene_id  # Assuming 'gene_id' is the column to use as row names
                df <- df[, -1]  # Remove the 'gene_id' column from the data frame
                
                if (all(rownames(df_2) == colnames(df))) {
                        print(paste(i, "is ready for DESeq2 analysis"))
                        assign(i, df, envir = .GlobalEnv)
                }
        }
}

# After checking and ensuring everything is correct, I will run DESeq now in a loop
cts <- c("pos_mir_cts", "neg_mir_cts")
metadata <- c("pos_mir", "neg_mir")

for (idx in seq_along(cts)) {
        df_name <- cts[idx]
        meta_name <- metadata[idx]
        
        df <- get(df_name)
        df_2 <- get(meta_name)
        
        dds <- DESeqDataSetFromMatrix(countData = round(df),
                                      colData = df_2,
                                      design = ~ group)
        
        dds <- DESeq(dds)
        res <- results(dds)
        res <- as.data.frame(res)
        
        file <- read.table(file_path, header = TRUE, sep = "\t")
        file <- as_tibble(file) %>%
                dplyr::select(gene_id, tpm_unstranded, gene_type, gene_name) %>%
                dplyr::filter(gene_type %in% c("miRNA", "protein_coding"))
        file <- file[5:nrow(file),c(1,3,4)]
        
        res$gene_id <- rownames(res)
        res <- merge(file, res, by = "gene_id")
        res <- res %>%
                filter(padj <= 0.05, log2FoldChange >= 1 | log2FoldChange <= -1 )
        
        assign(paste0(df_name, "_res"), res, envir = .GlobalEnv)
}

group_modulation_res <- pos_mir_cts_res
group_modulation_cts <- pos_mir_cts

group_not_modulation_res <- neg_mir_cts_res
group_not_modulation_cts <- neg_mir_cts

# Now I will perform over-representation analysis with the differentially expressed genes
pos_genes <- pos_mir_cts_res$gene_name
pos_genes <- enrichGO(gene = pos_genes,
                      OrgDb = org.Hs.eg.db, keyType = "SYMBOL", readable = TRUE,
                      ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.05)

neg_genes <- neg_mir_cts_res$gene_name
neg_genes <- enrichGO(gene = neg_genes,
                      OrgDb = org.Hs.eg.db, keyType = "SYMBOL", readable = TRUE,
                      ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.05)

# Found multiple signaling pathways involved with different cellular processes.
# Among the observed pathways in genes, I selected the most representative pathways that could be influencing migration and metastasis.
# They are:
metastasis <- c("tissue migration",
                "ameboidal-type cell migration",
                "epithelial cell migration",
                "extracellular matrix disassembly",
                "extracellular matrix organization",
                "positive regulation of cell adhesion",
                "mesenchymal cell differentiation")

# Now I will select in genes the data containing the regulated genes involved with these signaling pathways
library(enrichplot)
selected_genes <- pos_genes@result[pos_genes$Description %in% metastasis, ] # Select genes and pathways involved with metastasis
go_result <- pos_genes # Copy raw ORA result to a new variable
go_result@result <- selected_genes # Add only the pathway of interest to the result dataframe

# Now that I have the list of differentially expressed genes between my groups, and involved with metastasis, I will collect the log-2 fold change data to plot in the cnetplot. In this part of the article, we will state that "as expected, we observed regulation of pathways involved with cell migration, such as 'blablablabla'"
fold_change <- pos_mir_cts_res[,c(3,5)] # Select fold-change and padj data from deseq
core_genes <- str_split(as.data.frame(selected_genes)[,"geneID"], "/") # Extract gene names present in the pathways involved with metastasis
core_genes <- stringi::stri_list2matrix(core_genes)
core_genes <- as.data.frame(core_genes)
colnames(core_genes) <- metastasis

gene_list <- c(core_genes[[1]], core_genes[[2]],
               core_genes[[3]], core_genes[[4]],
               core_genes[[5]], core_genes[[6]],
               core_genes[[7]])

# Here I will have a string containing the filtered genes of interest
fold_change <- fold_change[fold_change$gene_name %in% gene_list, ] # Select genes that were previously filtered
fold_change <- fold_change$gene_name

# For loop to select my genes
core_genes # name of df
fold_change # list of genes

# For loop to select genes
for (col in 1:ncol(core_genes)) {
        for (row in 1:nrow(core_genes)) {
                if (core_genes[row, col] %in% fold_change) {
                        print(paste("Gene found:", core_genes[row, col]))
                } else {
                        print(paste("Gene removed:", core_genes[row, col]))
                        core_genes[row, col] <- NA
                }
        }
}

# Re-transform core_genes into a list
core_genes <- as.list(as.data.frame(core_genes))
listed_genes <- function(core_genes) {
        # Apply the collapse function to each vector in the list
        collapsed_genes <- sapply(core_genes, function(genes) {
                # Remove NA before collapsing
                genes <- genes[!is.na(genes)]
                # Join genes with "/"
                paste(genes, collapse = "/")
        })
        return(collapsed_genes)
}

# Adjust the dataframes containing the results of enrichplot
core_genes <- listed_genes(core_genes)
core_genes <- as.data.frame(core_genes)
core_genes$Description <- metastasis
colnames(core_genes)[1] <- "geneID"
go_filtered <- go_result

# Adjust the dataframes containing the enrichplot results
core_genes <- listed_genes(core_genes)
core_genes <- as.data.frame(core_genes)
core_genes$Description <- metastasis
colnames(core_genes)[1] <- "geneID"
go_filtered <- go_result
go_filtered@result$geneID <- core_genes$geneID
go_filtered@result <- go_filtered@result[-c(7), ]

# Select fold-change values for genes
fold_change_value <- pos_mir_cts_res[,c(3,5)]
fold_change_value <- fold_change_value[fold_change_value$gene_name %in% fold_change, ]
fold_change_value <- fold_change_value$log2FoldChange
names(fold_change_value) <- fold_change

# Now, create a cnetplot of filtered up-regulated genes
cnetplot(go_filtered, categorySize = "none", foldChange = fold_change_value,
         showCategory = 7, node_label = "gene") +
        scale_colour_gradient2(name = "Log2 Fold-Change", low = "darkblue", mid = "white", high = "darkred")

# ------ Part 7.2: Comparison of heatmap between the 4 groups ------
# After analyzing DEG between samples that cannot regulate ROCK via miR, we discovered that
# there is no modulation of specific migration pathways.
# Now, we will compare via heatmap the expression of migration genes between the groups that regulate ROCK by miR (group 1) or do not regulate ROCK by miR.

##################################### Migration Genes
# List of differentially expressed genes in group 1
core_genes_migration <- core_genes[3:4,]  # Select only migration genes
core_genes_migration <- str_split(as.data.frame(core_genes_migration)[,"geneID"], "/")  # Remove slash
core_genes_migration <- c(core_genes_migration[[1]], core_genes_migration[[2]])  # Filter only gene names
core_genes_migration <- unique(core_genes_migration)

# Next, select the genes from this list with higher log2 fold-change
core_migration_top <- pos_mir_cts_res[pos_mir_cts_res$gene_name %in% core_genes_migration, ]  # Select genes and pathways involved in metastasis
core_migration_top <- core_migration_top %>%
        dplyr::filter(log2FoldChange > 1)

core_migration_top <- core_migration_top$gene_name

# After selecting the genes with higher fold-change, now let's select the raw counts expression across different groups
group_modulation_cts 
group_not_modulation_cts

# Read gene IDs and names
gene_ids <- read.table(file_path, header = TRUE, sep = "\t")
gene_ids <- gene_ids %>%
        dplyr::select(gene_id, gene_name, gene_type) %>%
        dplyr::filter(gene_type %in% c("miRNA", "protein_coding")) %>%
        dplyr::select(gene_id, gene_name)

group_modulation_cts <- cbind(gene_ids, group_modulation_cts)
group_not_modulation_cts <- cbind(gene_ids, group_not_modulation_cts) 

group_modulation_migration <- merge(group_modulation_cts[,2:ncol(group_modulation_cts)], group_not_modulation_cts[,2:ncol(group_not_modulation_cts)], by = "gene_name")

group_modulation_migration <- group_modulation_migration %>% 
        filter(gene_name %in% core_migration_top)

rownames(group_modulation_migration) <- group_modulation_migration$gene_name
group_modulation_migration <- group_modulation_migration[,2:ncol(group_modulation_migration)]

# Annotation for heatmap
annotation_not_modulation <- rbind(group_1, group_2)
annotation_not_modulation$group <- sub("1", "High/High", annotation_not_modulation$group)
annotation_not_modulation$group <- sub("2", "Low/Low", annotation_not_modulation$group)
annotation_not_modulation$file_name <- rownames(annotation_not_modulation)
annotation_not_modulation <- annotation_not_modulation %>%
        dplyr::arrange(desc(ROCK2))
annotation_not_modulation <- annotation_not_modulation[,3:4]
annotation_not_modulation$type <- "Uncontrolled"

annotation_modulation <- rbind(group_3, group_4)
annotation_modulation$group <- sub("3", "High/Low", annotation_modulation$group)
annotation_modulation$group <- sub("4", "Low/High", annotation_modulation$group)
annotation_modulation$file_name <- rownames(annotation_modulation)
annotation_modulation <- annotation_modulation %>%
        dplyr::arrange(desc(ROCK2))
annotation_modulation <- annotation_modulation[,3:4]
annotation_modulation$type <- "Controlled"

annotation = rbind(annotation_not_modulation, annotation_modulation)
group_modulation_migration_filtered <- group_modulation_migration[, annotation$file_name]  # Select columns from counts that are present in annotation
annotation = annotation[,c(1,3)]
annotation_colors <- list(
        group = c("High/High" = "#212529", "Low/Low" = "#495057", "High/Low" = "#adb5bd", "Low/High" = "#dee2e6"))

# Heatmap
library(pheatmap)
pheatmap(group_modulation_migration_filtered, 
         cluster_rows = TRUE,  # Cluster rows
         cluster_cols = FALSE, scale = "row", # Do not cluster columns
         color = colorRampPalette(c("white", "#FBF2F2", "#AF0000"))(20),  # Color scheme
         show_rownames = TRUE,  # Show row names
         show_colnames = FALSE,
         clustering_distance_rows = "canberra",
         fontsize_row = 8,
         annotation_col = annotation,
         annotation_colors = annotation_colors)

###################################### Adhesion Loss
# Migration Genes
# List of differentially expressed genes in group 1
core_genes_adesion <- core_genes[6,]  # Select only adhesion loss genes
core_genes_adesion <- str_split(as.data.frame(core_genes_adesion)[,"geneID"], "/")  # Remove slash
core_genes_adesion <- c(core_genes_adesion[[1]])  # Filter only gene names
core_genes_adesion <- unique(core_genes_adesion)

# Next, select the genes from this list with lower log2 fold-change
core_adesion_top <- pos_mir_cts_res[pos_mir_cts_res$gene_name %in% core_genes_adesion, ]  # Select genes and pathways involved in metastasis
core_adesion_top <- core_adesion_top %>%
        dplyr::filter(log2FoldChange < 1)

core_adesion_top <- core_adesion_top$gene_name
# Uncomment if specific genes are known
# core_adesion_top <- c("SH3BP1", "CD40", "CORO1A", "CTSH", "HSPB1", "IFNG", "RRAS", "APOE", "AGT", "S100A9", "LRG1", "PTP4A3", "MIR126", "MIR200C", "TNF", "GPX1")

# After selecting the genes with lower fold-change, now let's select the raw counts expression across different groups
group_modulation_res <- pos_mir_cts_res
group_modulation_cts <- pos_mir_cts

group_not_modulation_res <- neg_mir_cts_res
group_not_modulation_cts <- neg_mir_cts

# Read gene IDs and names
gene_ids <- read.table(file_path, header = TRUE, sep = "\t")
gene_ids <- gene_ids %>%
        dplyr::select(gene_id, gene_name, gene_type) %>%
        dplyr::filter(gene_type %in% c("miRNA", "protein_coding")) %>%
        dplyr::select(gene_id, gene_name)

group_modulation_cts <- cbind(gene_ids, group_modulation_cts)
group_not_modulation_cts <- cbind(gene_ids, group_not_modulation_cts) 

group_modulation_adesion <- merge(group_modulation_cts[,2:ncol(group_modulation_cts)], group_not_modulation_cts[,2:ncol(group_not_modulation_cts)], by = "gene_name")

group_modulation_adesion <- group_modulation_adesion %>% 
        filter(gene_name %in% core_adesion_top)

rownames(group_modulation_adesion) <- group_modulation_adesion$gene_name
group_modulation_adesion <- group_modulation_adesion[,2:ncol(group_modulation_adesion)]

# Annotation for heatmap
annotation_not_modulation <- rbind(group_1, group_2)
annotation_not_modulation$group <- sub("1", "High/High", annotation_not_modulation$group)
annotation_not_modulation$group <- sub("2", "Low/Low", annotation_not_modulation$group)
annotation_not_modulation$file_name <- rownames(annotation_not_modulation)
annotation_not_modulation <- annotation_not_modulation %>%
        dplyr::arrange(desc(ROCK2))
annotation_not_modulation <- annotation_not_modulation[,3:4]
annotation_not_modulation$type <- "Uncontrolled"

annotation_modulation <- rbind(group_3, group_4)
annotation_modulation$group <- sub("3", "High/Low", annotation_modulation$group)
annotation_modulation$group <- sub("4", "Low/High", annotation_modulation$group)
annotation_modulation$file_name <- rownames(annotation_modulation)
annotation_modulation <- annotation_modulation %>%
        dplyr::arrange(desc(ROCK2))
annotation_modulation <- annotation_modulation[,3:4]
annotation_modulation$type <- "Controlled"

annotation = rbind(annotation_not_modulation, annotation_modulation)
group_modulation_adesion_filtered <- group_modulation_adesion[, annotation$file_name]  # Select columns from counts that are present in annotation
annotation = annotation[,c(1,3)]
annotation_colors <- list(
        group = c("High/High" = "#212529", "Low/Low" = "#495057", "High/Low" = "#adb5bd", "Low/High" = "#dee2e6"))

# Heatmap
library(pheatmap)
pheatmap(group_modulation_adesion_filtered, 
         cluster_rows = TRUE,  # Cluster rows
         cluster_cols = FALSE, scale = "row", # Do not cluster columns
         color = colorRampPalette(c("white", "#FBF2F2", "#AF0000"))(20),  # Color scheme
         show_rownames = TRUE,  # Show row names
         show_colnames = FALSE,
         clustering_distance_rows = "euclidean",
         fontsize_row = 8,
         annotation_col = annotation,
         annotation_colors = annotation_colors)
