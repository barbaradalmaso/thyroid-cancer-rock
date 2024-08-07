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
tumor_samples <- tumor_data
tumor_data <- tumor_samples

tumor_data <- tumor_data %>%
        filter(MIR221 <= 40 & MIR222 <= 6)

set.seed(123)  # Para reprodutibilidade
num_clusters <- 6 
data_selected <- tumor_data %>% select(MIR222, MIR221, ROCK2)
data_scaled <- scale(data_selected)
# Executar K-means
kmeans_result <- kmeans(data_scaled, centers = num_clusters)

# Adicionar os rótulos dos clusters ao conjunto de dados original
tumor_data$cluster <- kmeans_result$cluster

ggplot(tumor_data, aes(x = ROCK2, y = MIR221)) +
        geom_point(size = 2) +  # Adicionar pontos ao gráfico
        geom_smooth(method = "lm", color = "blue") +  # Adicionar linha de regressão linear
        labs(x = NULL, y = NULL) +
        facet_wrap(~cluster)

# Apos conseguir identificar os grupos, vou separar diferentes dataframes contendo os files-names de cada amostra. Apos definir os dataframes com os 8 grupos, vou fazer o download do RNA-seq novamente e comparar os DEG entre os grupos.

# miR-222
group_1 = tumor_data %>% #mirHIGH/rockHIGH - miR222
        filter(cluster %in% c(1,2,4,3))
group_1$group <- "1"
group_1 <- group_1 %>%
        filter(MIR221 > 15 & ROCK2 > 15)
group_1 <- group_1[,c(2,4,6)]


group_2 = tumor_data %>% #mirLOW/rockLOW - miR222
        filter(cluster %in% c(6))
group_2$group <- "2"
group_2 <- group_2 %>%
        filter(MIR221 < 10 & ROCK2 < 10)
group_2 <- group_2[,c(2,4,6)]

group_3 = tumor_data %>% #mirHIGH/rockLOW - miR222
        filter(cluster %in% c(2,4,5))
group_3$group <- "3"
group_3 <- group_3 %>%
        filter(MIR221 > 15 & ROCK2 < 15)
group_3 <- group_3[,c(2,4,6)]

group_4 = tumor_data %>% #mirlOW/rockHIGH - miR222
        filter(cluster %in% c(3))
group_4$group <- "4"
group_4 <- group_4 %>%
        filter(MIR221 < 15 & ROCK2 > 20)
group_4 <- group_4[,c(2,4,6)]


# ------------- Parte 4: Quais as implicações clinicas da modulação de ROCK2 pelos miRs? ----------------------
# Inicialmente, vou selecionar os grupos de pacientes onde existe uma regulacao de ROCK desencadeada pelos miRs. 
group_mir221 <- rbind(group_3, group_4)
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

#------------- Parte 6: Quais as implicações clinicas dos miRs quando nao conseguem modular ROCK? ----------------------
# Agora vou selecionar os grupos de pacientes onde NAO existe uma regulacao de ROCK desencadeada pelos miRs. Eles sao: 1 e 2 // 5 e 6

group_mir221 <- rbind(group_1, group_2)
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

group_mir221 = merge(filtered.metadata, group_mir221, by = "file_name")

######## Count the number of cases in each group and pathology stage
table(group_mir221$group, group_mir221$ajcc_pathologic_n)

# ------ Parte 7: Analise de genes diferencialmente expressos entre os grupos ----------
# ------ Parte 7.1: Analise de grupos de pacientes que conseguem OU NAO conseguem modular o eixo ROCK-miR ------
############ Deseq Analysis
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

# Apos conseguir separar os grupos, vou fazer o download do RNA-seq novamente e comparar TODOS os DEG entre os grupos.
# Os grupos que conseguem modular sao 'pos', e os que nao conseguem sao 'neg'. Meus grupos sao:
pos_mir <- rbind(group_3, group_4)
neg_mir <- rbind(group_1, group_2)

pos_mir$file_names <- rownames(pos_mir)
neg_mir$file_names <- rownames(neg_mir)

groups <- list(pos_mir = pos_mir, neg_mir = neg_mir)

path <- "/Volumes/Extreme SSD/thyroid-ml/rna-seq/"

# Loop para fazer o download dos dados completos de RNA-seq de amostras referentes a cada grupo
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
        
        # Garantindo que a primeira coluna seja gene_id
        gene_ids <- data_j$gene_id
        data_j <- data_j[ , !(names(data_j) %in% c("gene_id", "gene_type"))]
        data_j <- cbind(gene_id = gene_ids, data_j)
        
        # Salvando o data frame no ambiente global
        assign(paste0(group_name, "_cts"), data_j, envir = .GlobalEnv)
}

# Primeiro preciso checar se as colunas de cts tem os mesmos nomes dos rownames de matadata
cts <- c("pos_mir_cts", "neg_mir_cts")
metadata <- c("pos_mir", "neg_mir")
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

# Agora vou fazer a over-representation analysis com os genes diferencialmente expressos
pos_genes <- pos_mir_cts_res$gene_name
pos_genes <- enrichGO(gene = pos_genes,
                  OrgDb = org.Hs.eg.db, keyType = "SYMBOL", readable = TRUE,
                  ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.05)

neg_genes <- neg_mir_cts_res$gene_name
neg_genes <- enrichGO(gene = neg_genes,
                      OrgDb = org.Hs.eg.db, keyType = "SYMBOL", readable = TRUE,
                      ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.05)


# Encontrei multiplas vias de sinalizacao envolvidas com diferentes processos celulares.
# Dentre as vias observadas em genes, selecionei as vias mais representativas que poderiam estar influenciando migracao e metastase.
# Sao elas:
metastasis <- c("tissue migration",
                "ameboidal-type cell migration",
                "epithelial cell migration",
                "extracellular matrix disassembly",
                "extracellular matrix organization",
                "positive regulation of cell adhesion",
                "mesenchymal cell differentiation")

# Agora vou selecionar em genes os dados contendo os genes regulados e envolvidos com essas vias de sinalizacao
library(enrichplot)
selected_genes <- pos_genes@result[pos_genes$Description %in% metastasis, ] # Selecionar genes e vias de sinalizacao envolvidas com metastase
go_result <- pos_genes # Copiar resultado bruto da analise ORA para uma nova variavel
go_result@result <- selected_genes # Adicionar no dataframe de resultado, dentro do resultado bruto, somente a via de sinalizacao de interesse

# Agora que tenho a lista de genes diferencialmente expressos entre os meus grupos, e envolvidos com metastase, vou coletar o dado de log-2 fold change para conseguir plotar no cnetplot. Nessa parte do artigo, vamos colocar que "como esperado, observamos uma regulacao de vias envolvidas com migracao celular, como "blablablabla"
fold_change <- pos_mir_cts_res[,c(3,5)] # Selecionar dados de fold-change e padj pelo deseq
core_genes <- str_split(as.data.frame(selected_genes)[,"geneID"], "/") # Extrair nome dos genes presentes nas vias envolvidas com metastase
core_genes <- stringi::stri_list2matrix(core_genes)
core_genes <- as.data.frame(core_genes)
colnames(core_genes) <- metastasis

gene_list <- c(core_genes[[1]], core_genes[[2]],
               core_genes[[3]], core_genes[[4]],
               core_genes[[5]], core_genes[[6]],
               core_genes[[7]])

# Aqui eu vou ter um sring contendo os genes filtrados de meu interesse
fold_change <- fold_change[fold_change$gene_name %in% gene_list, ] # Selecionar genes que filtrei previamente
fold_change <- fold_change$gene_name

# for loop pra selecionar meus genes
core_genes #nome do df
fold_change #lista de genes

# For loop para selecionar genes
for (col in 1:ncol(core_genes)) {
        for (row in 1:nrow(core_genes)) {
                if (core_genes[row, col] %in% fold_change) {
                        print(paste("Gene encontrado:", core_genes[row, col]))
                } else {
                        print(paste("Gene apagado:", core_genes[row, col]))
                        core_genes[row, col] <- NA
                }
        }
}

# Retransformar core_genes em uma lista
core_genes <- as.list(as.data.frame(core_genes))
listed_genes <- function(core_genes) {
        # Aplica a função de colapso em cada vetor da lista
        collapsed_genes <- sapply(core_genes, function(genes) {
                # Remove NA antes de colapsar
                genes <- genes[!is.na(genes)]
                # Junta genes com "/"
                paste(genes, collapse = "/")
        })
        return(collapsed_genes)
}

# Ajustar os dataframes contendo os resultados de enrichplot
core_genes <- listed_genes(core_genes)
core_genes <- as.data.frame(core_genes)
core_genes$Description <- metastasis
colnames(core_genes)[1] <- "geneID"
go_filtered <- go_result
go_filtered@result$geneID <- core_genes$geneID
go_filtered@result <- go_filtered@result[-c(7), ]

# Selecionar fg de genes
fold_change_value <- pos_mir_cts_res[,c(3,5)]
fold_change_value <- fold_change_value[fold_change_value$gene_name %in% fold_change, ]
fold_change_value <- fold_change_value$log2FoldChange
names(fold_change_value) <- fold_change

# Agora vou fazer cnetplot de filtered genes up regulated
cnetplot(go_filtered, categorySize = "none", foldChange = fold_change_value,
         showCategory = 7, node_label = "gene") +
         scale_colour_gradient2(name = "Log2 Fold-Change", low = "darkblue", mid = "white", high = "darkred"
                                        )
        
# ------ Parte 7.2: Analise de comparacao de heatmap entre os 4 grupos ------
# Depois de fazer as analises de DEG entre amostras que nao conseguem regular ROCK via miR, descobrimos que
# nao ocorre modulacao de vias especificas de migracao.
# entao agora vou comparar via heatmap a expressao de genes de migracao entre os grupos que regulam ROCK por miR (grupo 1) ou nao regulam ROCK por mir

##################################### Genes de migracao
core_genes # lista de genes que sao diferencialmente expressos no grupo 1
core_genes_migration <- core_genes[3:4,] # Selecinar so genes de migracao
core_genes_migration <- str_split(as.data.frame(core_genes_migration)[,"geneID"], "/") # Tirar barra
core_genes_migration <- c(core_genes_migration[[1]], core_genes_migration[[2]]) # Filtrar somente o nome dos genes
core_genes_migration <- unique(core_genes_migration)

# Em seguida selecionar os genes dessa lista com maior log2fg
core_migration_top <- pos_mir_cts_res[pos_mir_cts_res$gene_name %in% core_genes_migration, ] # Selecionar genes e vias de sinalizacao envolvidas com metastase
core_migration_top <- core_migration_top %>%
        dplyr::filter(log2FoldChange > 1)

core_migration_top <- core_migration_top$gene_name

# Apos selecoinar os genes com maior  fd, agora vamos selecionar a expressao bruta por cts nos diferentes grupos
group_modulation_cts 
group_not_modulation_cts

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
group_modulation_migration_filtered <- group_modulation_migration[, annotation$file_name] # Selecionar colunas de cts que estejam presentes no annotation
annotation = annotation[,c(1,3)]
annotation_colors <- list(
        group = c("High/High" = "#212529", "Low/Low" = "#495057", "High/Low" = "#adb5bd", "Low/High" = "#dee2e6"))

# Heatmap
library(pheatmap)
pheatmap(group_modulation_migration_filtered, 
         cluster_rows = TRUE,  # Agrupar as linhas
         cluster_cols = FALSE, scale = "row", # Agrupar as colunas
         color = colorRampPalette(c("white", "#FBF2F2", "#AF0000"))(20),  # Esquema de cores
         show_rownames = TRUE,  # N??o mostrar os nomes das linhas
         show_colnames = FALSE,
         clustering_distance_rows = "canberra",
         fontsize_row = 8,
         annotation_col = annotation,
         annotation_colors = annotation_colors)

###################################### Perda de adesao
# Genes de migracao
core_genes # lista de genes que sao diferencialmente expressos no grupo 1
core_genes_adesion <- core_genes[6,] # Selecinar so genes de migracao
core_genes_adesion <- str_split(as.data.frame(core_genes_adesion)[,"geneID"], "/") # Tirar barra
core_genes_adesion <- c(core_genes_adesion[[1]]) # Filtrar somente o nome dos genes
core_genes_adesion <- unique(core_genes_adesion)

# Em seguida selecionar os genes dessa lista com maior log2fg
core_adesion_top <- pos_mir_cts_res[pos_mir_cts_res$gene_name %in% core_genes_adesion, ] # Selecionar genes e vias de sinalizacao envolvidas com metastase
core_adesion_top <- core_adesion_top %>%
        dplyr::filter(log2FoldChange < 1)

core_adesion_top <- core_adesion_top$gene_name
# core_adesion_top <- c("SH3BP1", "CD40", "CORO1A", "CTSH", "HSPB1", "IFNG", "RRAS", "APOE", "AGT", "S100A9",  "LRG1", "PTP4A3",  "MIR126", "MIR200C", "TNF", "GPX1")

# Apos selecoinar os genes com maior  fd, agora vamos selecionar a expressao bruta por cts nos diferentes grupos
group_modulation_res <- pos_mir_cts_res
group_modulation_cts <- pos_mir_cts

group_not_modulation_res <- neg_mir_cts_res
group_not_modulation_cts <- neg_mir_cts

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
group_modulation_adesion_filtered <- group_modulation_adesion[, annotation$file_name] # Selecionar colunas de cts que estejam presentes no annotation
annotation = annotation[,c(1,3)]
annotation_colors <- list(
        group = c("High/High" = "#212529", "Low/Low" = "#495057", "High/Low" = "#adb5bd", "Low/High" = "#dee2e6"))

# Heatmap
library(pheatmap)
pheatmap(group_modulation_adesion_filtered, 
         cluster_rows = TRUE,  # Agrupar as linhas
         cluster_cols = FALSE, scale = "row", # Agrupar as colunas
         color = colorRampPalette(c("white", "#FBF2F2", "#AF0000"))(20),  # Esquema de cores
         show_rownames = TRUE,  # N??o mostrar os nomes das linhas
         show_colnames = FALSE,
         clustering_distance_rows = "euclidean",
         fontsize_row = 8,
         annotation_col = annotation,
         annotation_colors = annotation_colors)

