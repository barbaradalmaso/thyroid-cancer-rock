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

# Supondo que tumor_data seja o seu dataframe original
tumor_data_long <- tumor_data %>%
        pivot_longer(cols = c(MIR222, MIR221, ROCK1, ROCK2), names_to = "Gene", values_to = "Expression")

# Criar os density plots lado a lado
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
# Reorganizar os dados para o formato longo
library(ggpubr)
tumor_data_long <- tumor_data %>%
        pivot_longer(cols = c(MIR222, MIR221), names_to = "MIR", values_to = "MIR_Expression") %>%
        pivot_longer(cols = c(ROCK1, ROCK2), names_to = "ROCK", values_to = "ROCK_Expression")

# Criar os scatter plots lado a lado com valores de correlação
ggplot(tumor_data_long, aes(x = ROCK_Expression, y = MIR_Expression)) +
        geom_point(size = 2) +  # Adicionar pontos ao gráfico
        geom_smooth(method = "lm", color = "blue") +  # Adicionar linha de regressão linear
        labs(x = NULL, y = NULL) +  # Remover títulos dos eixos
        theme(
                panel.background = element_rect(size = 0.5, linetype = "solid"),
                panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"), 
                panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "white"),
                axis.text.x = element_text(size = 12),  # Aumentar o tamanho do texto do eixo x
                axis.text.y = element_text(size = 12)   # Aumentar o tamanho do texto do eixo y
        ) +
        stat_cor(method = "pearson", label.x = 45, label.y = 250, size = 4) +  # Adicionar valores de correlação de Pearson
        facet_grid(MIR ~ ROCK)  # Criar facetas para cada combinação de MIR e ROCK

"Para ROCK2 encontramos uma correlação muito baixa. Isso pode sugerir que algumas amostras possuem relação positiva MIR/ROCK,
mas outras podem possuir uma relação negativa. Vou tentar separar essas amostras"


# ----------- Parte 3: Agrupamento de amostras com base na expressao de ROCK/MIR222 MIR221 ----------
# Definir o número de clusters
# Primeiro vou tentar tirar alguns outliers
tumor_data <- tumor_data %>%
        filter(MIR221 <= 60 & MIR222 <= 10)


set.seed(123)  # Para reprodutibilidade
num_clusters <- 6 
data_selected <- tumor_data %>% select(MIR222, MIR221, ROCK2)
data_scaled <- scale(data_selected)
# Executar K-means
kmeans_result <- kmeans(data_scaled, centers = num_clusters)

# Adicionar os rótulos dos clusters ao conjunto de dados original
tumor_data$cluster <- kmeans_result$cluster

ggplot(tumor_data, aes(x = ROCK2, y = MIR222)) +
        geom_point(size = 2) +  # Adicionar pontos ao gráfico
        geom_smooth(method = "lm", color = "blue") +  # Adicionar linha de regressão linear
        labs(x = NULL, y = NULL) +
        facet_wrap(~cluster)

# Apos conseguir identificar os grupos, vou separar diferentes dataframes contendo os files-names de cada amostra. Apos definir os dataframes com os 8 grupos, vou fazer o download do RNA-seq novamente e comparar os DEG entre os grupos.

# miR-222
group_1 = tumor_data %>% #mirHIGH/rockHIGH - miR222
        filter(cluster %in% c(4))
group_1 = group_1[,2:3]
cor.test(group_1$MIR222, group_1$ROCK2)
group_1$group <- "group_1"

group_2 = tumor_data %>% #mirLOW/rockLOW - miR222
        filter(cluster %in% c(3))
group_2 = group_2[,2:3]
cor.test(group_2$MIR222, group_2$ROCK2)
group_2$group <- "group_2"


group_3 = tumor_data %>% #mirHIGH/rockLOW - miR222
        filter(cluster %in% c(1))
group_3 = group_3[,2:3]
cor.test(group_3$MIR222, group_3$ROCK2)
group_3$group <- "group_3"


group_4 = tumor_data %>% #mirlOW/rockHIGH - miR222
        filter(cluster %in% c(6))
group_4 = group_4[,2:3]
cor.test(group_4$MIR222, group_4$ROCK2)
group_4$group <- "group_4"

# miR-221
group_5 = tumor_data %>% #mirHIGH/rockHIGH - miR221
        filter(cluster %in% c(4))
group_5 = group_5[,c(2,4)]
cor.test(group_5$MIR221, group_5$ROCK2)
group_5$group <- "group_5"

group_6 = tumor_data %>% #mirLOW/rockLOW - miR221
        filter(cluster %in% c(3))
group_6 = group_6[,c(2,4)]
cor.test(group_6$MIR221, group_6$ROCK2)
group_6$group <- "group_6"

group_7 = tumor_data %>% #mirHIGH/rockLOW - miR221
        filter(cluster %in% c(1))
group_7 = group_7[,c(2,4)]
cor.test(group_7$MIR221, group_3$ROCK2)
group_7$group <- "group_7"

group_8 = tumor_data %>% #mirlOW/rockHIGH - miR221
        filter(cluster %in% c(6))
group_8 = group_8[,c(2,4)]
cor.test(group_8$MIR221, group_8$ROCK2)
group_8$group <- "group_8"

metadata_mir222_1 <- rbind(group_2, group_3)
metadata_mir222_1$file_names <- rownames(metadata_mir222_1)
metadata_mir222_2 <- rbind(group_1, group_4)
metadata_mir222_2$file_names <- rownames(metadata_mir222_2)


metadata_mir221_1 <- rbind(group_6, group_7)
metadata_mir221_1$file_names <- rownames(metadata_mir221_1)
metadata_mir221_2 <- rbind(group_5, group_8)
metadata_mir221_2$file_names <- rownames(metadata_mir221_2)

# ------------- Parte 4: Quais as implicações clinicas da modulação de ROCK2 pelos miRs? ----------------------
# Inicialmente, vou selecionar os grupos de pacientes onde existe uma regulacao de ROCK desencadeada pelos miRs. Eles sao: 3 e 4 // 7 e 8
group_mir222 <- rbind(group_3, group_4)
group_mir221 <- rbind(group_7, group_8)

group_mir222$file_name <- rownames(group_mir222)
group_mir221$file_name <- rownames(group_mir221)

# Read pathology metadata
pathology_metadata <- read.delim("/Volumes/Extreme SSD/thyroid-ml/clinical.cart.2024-07-16/clinical.tsv", sep = "\t")
pathology_metadata <- pathology_metadata[,c(1,29:32)] # Select relevant columns
pathology_metadata <- lapply(pathology_metadata, function(x) gsub("'--", NA, x)) # Replace '--' with NA
pathology_metadata <- as.data.frame(pathology_metadata) # Convert to data frame

# Segregate ROCK1/2 high and low groups based on expression levels
filtered.metadata = tumor.metadata[,c(1,36,40)] # Select relevant columns from tumor metadata
filtered.metadata = merge(filtered.metadata, pathology_metadata, by = "case_id") # Merge metadata with expression data
filtered.metadata = distinct(filtered.metadata, file_name, .keep_all = TRUE)

group_mir222 = merge(filtered.metadata, group_mir222, by = "file_name")
group_mir221 = merge(filtered.metadata, group_mir221, by = "file_name")

######## Count the number of cases in each group and pathology stage
table(group_mir221$group, group_mir221$ajcc_pathologic_n)

# At the end, the statistical analysis was performed on graphpad"

# ------------ Parte 5: Qual a sobrevivencia de pacientes a partir da da modulação de ROCK2 pelos miRs? ---------------

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
bar_code <- merge(bar_code, group_mir222, by = "case_id")
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

#------------- Parte 5: Quais as implicações clinicas dos miRs quando nao conseguem modular ROCK? ----------------------
# Agora vou selecionar os grupos de pacientes onde NAO existe uma regulacao de ROCK desencadeada pelos miRs. Eles sao: 1 e 2 // 5 e 6

group_mir222 <- rbind(group_1, group_2)
group_mir221 <- rbind(group_5, group_6)

group_mir222$file_name <- rownames(group_mir222)
group_mir221$file_name <- rownames(group_mir221)

# Read pathology metadata
pathology_metadata <- read.delim("/Volumes/Extreme SSD/thyroid-ml/clinical.cart.2024-07-16/clinical.tsv", sep = "\t")
pathology_metadata <- pathology_metadata[,c(1,29:32)] # Select relevant columns
pathology_metadata <- lapply(pathology_metadata, function(x) gsub("'--", NA, x)) # Replace '--' with NA
pathology_metadata <- as.data.frame(pathology_metadata) # Convert to data frame

# Segregate ROCK1/2 high and low groups based on expression levels
filtered.metadata = tumor.metadata[,c(1,36,40)] # Select relevant columns from tumor metadata
filtered.metadata = merge(filtered.metadata, pathology_metadata, by = "case_id") # Merge metadata with expression data
filtered.metadata = distinct(filtered.metadata, file_name, .keep_all = TRUE)

group_mir222 = merge(filtered.metadata, group_mir222, by = "file_name")
group_mir221 = merge(filtered.metadata, group_mir221, by = "file_name")

######## Count the number of cases in each group and pathology stage
table(group_mir221$group, group_mir221$ajcc_pathologic_n)









# ------ Parte 6: Analise de genes diferencialmente expressos entre os grupos ----------
############ Deseq Analysis
library(dplyr)
library(ggplot2)
library(pROC)
library(stringr)
library(tidyr)
library(ggpattern)
library(ggridges)
library(DESeq2)

# Apos conseguir separar os grupos, vou fazer o download do RNA-seq novamente e comparar TODOS os DEG entre os grupos.
# Meus grupos sao:
groups <- list(metadata_mir221_1 = metadata_mir221_1, metadata_mir221_2 = metadata_mir221_2,
               metadata_mir222_1 = metadata_mir222_1, metadata_mir222_2 = metadata_mir222_2)

path <- "/Volumes/Extreme SSD/thyroid-ml/rna-seq/"

# Loop sobre os grupos
for (group_name in names(groups)) {
        group <- groups[[group_name]]
        file_names <- group$file_name
        data_j <- data.frame()
        
        for (file_name in file_names) {
                file_path <- paste0(path, file_name)
                
                if (file.exists(file_path)) {
                        file <- read.table(file_path, header = TRUE, sep = "\t")
                        file <- file %>%
                                select(gene_id, tpm_unstranded, gene_type) %>%
                                filter(gene_type %in% c("miRNA", "protein_coding"))
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
        
        # Garantindo que a primeira coluna seja gene_id
        gene_ids <- data_j$gene_id
        data_j <- data_j[ , !(names(data_j) %in% c("gene_id", "gene_type"))]
        data_j <- cbind(gene_id = gene_ids, data_j)
        
        # Salvando o data frame no ambiente global
        assign(paste0(group_name, "_cts"), data_j, envir = .GlobalEnv)
}



# Primeiro preciso checar se as colunas de cts tem os mesmos nomes dos rownames de matadata
cts <- c("metadata_mir221_1_cts", "metadata_mir221_2_cts", "metadata_mir222_1_cts", "metadata_mir222_2_cts")
metadata <- c("metadata_mir221_1", "metadata_mir221_2", "metadata_mir222_1", "metadata_mir222_2")

for (i in cts) {
        for (j in metadata) {
                df <- get(i)
                df_2 <- get(j)
                rownames(df) <- df$gene_id  # Supondo que 'gene_id' seja a coluna que você deseja usar como nomes das linhas
                df <- df[, -1]  # Remover a coluna 'gene_id' do data frame
                
                if (all(rownames(df_2) == colnames(df))) {
                        print(paste(i, "is ready for DESEQ2 analysis"))
                        assign(i, df, envir = .GlobalEnv)
                }
        }
}

# Depois de checar e ver que esta tudo certo, vou correr o DESEQ agora em forma de loop
cts <- c("metadata_mir221_1_cts", "metadata_mir221_2_cts", "metadata_mir222_1_cts", "metadata_mir222_2_cts")
metadata <- c("metadata_mir221_1", "metadata_mir221_2", "metadata_mir222_1", "metadata_mir222_2")

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
                file <- file %>%
                        select(gene_id, tpm_unstranded, gene_type, gene_name) %>%
                        filter(gene_type %in% c("miRNA", "protein_coding"))
                file <- file[5:nrow(file),c(1,3,4)]
                
                res$gene_id <- rownames(res)
                res <- merge(file, res, by = "gene_id")
                res <- res %>%
                        filter(padj <= 0.05)
                
                
                assign(paste0(df_name, "res"), res, envir = .GlobalEnv)
                
                
        }
        
# Agora vou tentar observar quais sao os overlapping genes entre os grupos:
# " miR High/ ROCK High + miR Low / Rock Low x  miR Low/ ROCK High + miR Low / Rock High"

