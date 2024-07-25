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

# Filter sample metadata for specific sample types
sample.metadata <- filter(sample.metadata, sample_type %in% c("Solid Tissue Normal", "Primary Tumor"))

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
expr_control_rock1 <- t(control_data[1021,2:ncol(control_data)])
expr_tumor_rock1 <- t(tumor_data[1021,2:ncol(tumor_data)])

expr_control_rock2 <- t(control_data[6789,2:ncol(control_data)])
expr_tumor_rock2 <- t(tumor_data[6789,2:ncol(tumor_data)])

"From these initial results, we observe that after the ontogenic transformation in thyroid tumors, there is an increase in the expression of ROCK1 and ROCK2. In the next step, I will check the correlation of ROCK1 and ROCK2 expression with the project's target microRNAs: MIR221 and MIR222."


# ---------------- Part 2: Correlation analysis between ROCK1/2 and MIR221/222 on thyroid tumor samples ----------------

# Calculate Pearson correlation test between ROCK2 and MIR222 in tumor data
cor_test <- cor.test(tumor_data$ROCK2, tumor_data$MIR222, method = "pearson")

# Create a string with the correlation test results
cor_results <- paste("r =", round(cor_test$estimate, 2), 
                     "\np-value =", format.pval(cor_test$p.value, digits = 2))

# Create a scatter plot with ggplot2 and add correlation results
ggplot(tumor_data, aes(x = ROCK2, y = MIR222)) +
        geom_point(size = 2) +  # Add points to the plot
        geom_smooth(method = "lm", color = "blue") +  # Add linear regression line
        labs(x = NULL, y = NULL) +  # Remove axis titles
        theme(
                panel.background = element_rect(size = 0.5, linetype = "solid"),
                panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"), 
                panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "white"),
                axis.text.x = element_text(size = 12),  # Increase x-axis text size
                axis.text.y = element_text(size = 12)   # Increase y-axis text size
        )

# Create density plot for MIR222 expression in tumor data
ggplot(tumor_data, aes(x = MIR222)) +
        geom_density(color = "black", size = 1) +  # Add density line
        labs(title = NULL, x = NULL, y = NULL) +
        theme(
                panel.background = element_rect(size = 0.5, linetype = "solid"),
                panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"), 
                panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "white"),
                axis.text.x = element_text(size = 12),  # Increase x-axis text size
                axis.text.y = element_text(size = 12)   # Increase y-axis text size
        ) +
        coord_flip()  # Flip coordinates


"For each gene/MIR analyzed, change the arguments in the correlation analysis functions and graph generation. 
As a result, we observed a moderate positive correlation between ROCK1 and MIR221/222, and ROCK2 and MIR222 (approximately ~0.5). Thus, although the relationship is not extremely strong, we see that when the expression of ROCK1 increases, the expression of MIR221/222 tends to increase as well, as does the expression of ROCK2 and MIR222. The correlation between ROCK2 and MIR221 expression is very low and not significant (~0.1).

We also observed the overall distribution of ROCK1/2 and MIR221/222 expression in thyroid tumor samples. The histogram of both ROCK1/2 suggests that there are distinct populations of patients with its high and low expression."


# ----------- ROCK1/2 expression and tumor AJCC pathology ----------

# Select only tumor samples and combine ROCK1 and ROCK2 expression
expr_tumor_rock = cbind(expr_tumor_rock2, expr_tumor_rock1)
colnames(expr_tumor_rock) = c("ROCK2", "ROCK1")

# Read pathology metadata
pathology_metadata <- read.delim("/Volumes/Extreme SSD/thyroid-ml/clinical.cart.2024-07-16/clinical.tsv", sep = "\t")
pathology_metadata <- pathology_metadata[,c(1,29:32)] # Select relevant columns
pathology_metadata <- lapply(pathology_metadata, function(x) gsub("'--", NA, x)) # Replace '--' with NA
pathology_metadata <- as.data.frame(pathology_metadata) # Convert to data frame

# Segregate ROCK1/2 high and low groups based on expression levels
filtered.metadata = tumor.metadata[,c(1,36,40)] # Select relevant columns from tumor metadata
tumor_data_filtered = tumor_data
tumor_data_filtered$file_name = rownames(tumor_data) # Add file names to the filtered data
filtered.metadata = merge(filtered.metadata, tumor_data_filtered, by = "file_name") # Merge metadata with expression data

# Configure ROCK1/2 expression table by calculating quartiles
q1.ROCK1 = quantile(tumor_data$ROCK1, 0.25, na.rm = TRUE)
q4.ROCK1 = quantile(tumor_data$ROCK1, 0.75, na.rm = TRUE)
q1.ROCK2 = quantile(tumor_data$ROCK2, 0.25, na.rm = TRUE)
q4.ROCK2 = quantile(tumor_data$ROCK2, 0.75, na.rm = TRUE)

# Segregate samples into low, high, and control groups for ROCK1
ROCK1_neg <- subset(filtered.metadata, ROCK1 <= q1.ROCK1)
ROCK1_neg$type <- c(rep("Low", nrow(ROCK1_neg)))
ROCK1_pos <- subset(filtered.metadata, ROCK1 >= q4.ROCK1)
ROCK1_pos$type <- c(rep("High", nrow(ROCK1_pos)))
ROCK1_ct <- subset(filtered.metadata, ROCK1 <= q4.ROCK1 & ROCK1 >= q1.ROCK1)
ROCK1_ct$type <- c(rep("Control", nrow(ROCK1_ct)))
ROCK1_data <- rbind(ROCK1_neg, ROCK1_pos, ROCK1_ct)

# Segregate samples into low, high, and control groups for ROCK2
ROCK2_neg <- subset(filtered.metadata, ROCK2 <= q1.ROCK2)
ROCK2_neg$type <- c(rep("Low", nrow(ROCK2_neg)))
ROCK2_pos <- subset(filtered.metadata, ROCK2 >= q4.ROCK2)
ROCK2_pos$type <- c(rep("High", nrow(ROCK2_pos)))
ROCK2_ct <- subset(filtered.metadata, ROCK2 <= q4.ROCK2 & ROCK2 >= q1.ROCK2)
ROCK2_ct$type <- c(rep("Control", nrow(ROCK2_ct)))
ROCK2_data <- rbind(ROCK2_neg, ROCK2_pos, ROCK2_ct)

# Create a combined expression data frame for ROCK1 and ROCK2
ROCK_expression <- ROCK1_data[,c(1,2,ncol(ROCK1_data))]
colnames(ROCK_expression)[3] <- "ROCK1" 
ROCK2_data = ROCK2_data[,c(1, ncol(ROCK2_data))]
colnames(ROCK2_data)[2] <- "ROCK2" 
ROCK_expression = merge(ROCK_expression, ROCK2_data, by = "file_name")

# Merge ROCK expression data with pathology metadata
metadata.expression = merge(ROCK_expression, pathology_metadata, by = "case_id")
metadata.expression = unique(metadata.expression)
metadata.expression.ROCK1 = metadata.expression[metadata.expression$ROCK1 != "High", ]
metadata.expression.ROCK2 = metadata.expression[metadata.expression$ROCK2 != "High", ]

######## Perform statistical analysis using chi-squared test
contingency_table <- table(metadata.expression$ROCK1, metadata.expression$ajcc_pathologic_n)
chisq.test(contingency_table)

######## Count the number of cases in each group and pathology stage
metadata.expression.ROCK2 %>%
        group_by(ROCK2, ajcc_pathologic_t) %>%
        summarise(count = n(), .groups = 'drop')

# Export data and construct graphs using graph pad prism
metadata.expression <- lapply(metadata.expression, function(x) gsub("High", NA, x)) # Replace High group with NA
metadata.expression = as.data.frame(metadata.expression)

write.table(metadata.expression, 
            file = "/Volumes/Extreme SSD/thyroid-ml/thyroid.cancer/metadata.ROCK.csv", 
            row.names = FALSE, 
            col.names = TRUE, 
            sep = ",", 
            quote = FALSE)

"Using the chi-squared test, we observed that there is a trend of decreased AJCC lymph node metastasis in patients with low ROCK1 expression. We did not find any differences in other pathological factors, such as stage and tumor size."

# ------------------- Part 3: ROCK1/2 expression and patient survival analysis --------------
# ------------------- Part 3: ROCK1/2 expression and patient survival analysis --------------
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
bar_code <- merge(bar_code, ROCK_expression, by = "case_id")
bar_code <- unique(bar_code)
colnames(bar_code)[2] <- "bcr_patient_barcode"

# Create the clinical data
clin <- survivalTCGA(THCA.clinical, extract.cols="admin.disease_code")
clin <- merge(bar_code, clin, by = "bcr_patient_barcode")
clin.ROCK1 = clin[clin$ROCK1 != "High", ]
clin.ROCK2 = clin[clin$ROCK2 != "High", ]

# Tabulate survival data by outcome for ROCK1
xtabs(~ROCK1 + patient.vital_status, data = clin.ROCK1) %>% addmargins()
coxph(Surv(times, patient.vital_status) ~ ROCK1, data = clin.ROCK1)
sfit <- survfit(Surv(times, patient.vital_status) ~ ROCK1, data = clin.ROCK1)

# Plot survival analysis for ROCK1
ggsurvplot(sfit, conf.int = TRUE, pval = TRUE, risk.table = TRUE, 
           legend.labs = c("Control", "ROCK1 Low"),  
           palette = c("dodgerblue2", "#AFABAB"))

# Tabulate survival data by outcome for ROCK2
xtabs(~ROCK2 + patient.vital_status, data = clin.ROCK2) %>% addmargins()
coxph(Surv(times, patient.vital_status) ~ ROCK2, data = clin.ROCK2)
sfit <- survfit(Surv(times, patient.vital_status) ~ ROCK2, data = clin.ROCK2)

# Plot survival analysis for ROCK2
ggsurvplot(sfit, conf.int = TRUE, pval = TRUE, risk.table = TRUE, 
           legend.labs = c("Control", "ROCK2 Low"),  
           palette = c("dodgerblue2", "#AFABAB"))

"We observed that the decrease in ROCK1 expression tends to increase the survival probability of thyroid cancer patients compared to controls. No difference was observed between the ROCK2 Low group and controls."

