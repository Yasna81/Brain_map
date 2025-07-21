#amygdal
library(readr)
file <- read.delim("~/project/brain_map/tsv_files/amygdala_tpm.tsv", header = TRUE , sep = "\t")
library(biomaRt)
ensembl <- useMart("ensembl",dataset = "hsapiens_gene_ensembl")
ensembl_ids <- file$Name

anot <- getBM(
    attributes = c("ensembl_gene_id" , "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = ensembl
)

exp_anot_amyg <- merge(file, anot ,by.x = "Name",by.y = "ensembl_gene_id",all.x = TRUE)
exp_anot_amyg$Name <- sub("\\..*","",exp_anot_amyg$Name)
write.csv(exp_anot_amyg,"~/project/brain_map/csv_files/exp_anot_amyg.csv", row.names = FALSE)
amg <- read.csv("~/project/brain_map/csv_files/exp_anot_amyg.csv")
library(dplyr)
amg$amygdal <- rowMeans(amg[,-c(1,2)],na.rm = TRUE)
mean_tpm_amg <- amg[,c(1,2, ncol(amg))]
#write.csv(mean_tpm_amg,"~/project/brain_map/mean_tpm/mean_tpm_amg_1.csv",row.names = FALSE)
colnames(mean_tpm_amg)[1] <- "gene_id"
colnames(mean_tpm_amg)[2] <- "gene_symbol"


# hipocomp
library(readr)
file_h <- read.delim("~/project/brain_map/tsv_files/gene_tpm_v10_brain_hippocampus.tsv", header = TRUE , sep = "\t")
library(biomaRt)
ensembl <- useMart("ensembl",dataset = "hsapiens_gene_ensembl")
ensembl_ids <- file_h$Name

anot <- getBM(
    attributes = c("ensembl_gene_id" , "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = ensembl
)

exp_anot_hipc <- merge(file_h, anot ,by.x = "Name",by.y = "ensembl_gene_id",all.x = TRUE)
exp_anot_hipc$Name <- sub("\\..*","",exp_anot_hipc$Name)
#write.csv(exp_anot_amyg,"~/project/brain_map/csv_files/exp_anot_amyg.csv", row.names = FALSE)
#amg <- read.csv("~/project/brain_map/csv_files/exp_anot_amyg.csv")
library(dplyr)
exp_anot_hipc$hipocomp <- rowMeans(exp_anot_hipc[,-c(1,2)],na.rm = TRUE)
mean_tpm_hipc <- exp_anot_hipc[,c(1,2, ncol(exp_anot_hipc))]
#write.csv(mean_tpm_amg,"~/project/brain_map/mean_tpm/mean_tpm_amg_1.csv",row.names = FALSE)
colnames(mean_tpm_hipc)[1] <- "gene_id"
colnames(mean_tpm_hipc)[2] <- "gene_symbol"

#hipotal
#library(readr)
file_hp <- read.delim("~/project/brain_map/tsv_files/gene_tpm_v10_brain_hypothalamus.tsv", header = TRUE , sep = "\t")
#library(biomaRt)
ensembl <- useMart("ensembl",dataset = "hsapiens_gene_ensembl")
ensembl_ids <- file_hp$Name

anot <- getBM(
    attributes = c("ensembl_gene_id" , "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = ensembl
)

exp_anot_hypo <- merge(file_hp, anot ,by.x = "Name",by.y = "ensembl_gene_id",all.x = TRUE)
exp_anot_hypo$Name <- sub("\\..*","",exp_anot_hypo$Name)
#write.csv(exp_anot_amyg,"~/project/brain_map/csv_files/exp_anot_amyg.csv", row.names = FALSE)
#amg <- read.csv("~/project/brain_map/csv_files/exp_anot_amyg.csv")
#library(dplyr)
exp_anot_hypo$hypotalamus <- rowMeans(exp_anot_hypo[,-c(1,2)],na.rm = TRUE)
mean_tpm_hypo <- exp_anot_hypo[,c(1,2, ncol(exp_anot_hypo))]
#write.csv(mean_tpm_amg,"~/project/brain_map/mean_tpm/mean_tpm_amg_1.csv",row.names = FALSE)
colnames(mean_tpm_hypo)[1] <- "gene_id"
colnames(mean_tpm_hypo)[2] <- "gene_symbol"

#frontal cortex

#library(readr)
file_c_f <- read.delim("~/project/brain_map/tsv_files/gene_tpm_v10_brain_frontal_cortex_ba9.tsv", header = TRUE , sep = "\t")
#library(biomaRt)
ensembl <- useMart("ensembl",dataset = "hsapiens_gene_ensembl")
ensembl_ids <- file_c_f$Name

anot <- getBM(
    attributes = c("ensembl_gene_id" , "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = ensembl
)

exp_anot_c_f <- merge(file_c_f, anot ,by.x = "Name",by.y = "ensembl_gene_id",all.x = TRUE)
exp_anot_c_f$Name <- sub("\\..*","",exp_anot_c_f$Name)
#write.csv(exp_anot_amyg,"~/project/brain_map/csv_files/exp_anot_amyg.csv", row.names = FALSE)
#amg <- read.csv("~/project/brain_map/csv_files/exp_anot_amyg.csv")
#library(dplyr)
exp_anot_c_f$frontal_cortex <- rowMeans(exp_anot_c_f[,-c(1,2)],na.rm = TRUE)
mean_tpm_c_f <- exp_anot_c_f[,c(1,2, ncol(exp_anot_c_f))]
#write.csv(mean_tpm_amg,"~/project/brain_map/mean_tpm/mean_tpm_amg_1.csv",row.names = FALSE)
colnames(mean_tpm_c_f)[1] <- "gene_id"
colnames(mean_tpm_c_f)[2] <- "gene_symbol"

#antrior cortex 


#library(readr)
file_c_a <- read.delim("~/project/brain_map/tsv_files/gene_tpm_v10_brain_anterior_cingulate_cortex_ba24.tsv", header = TRUE , sep = "\t")
#library(biomaRt)
ensembl <- useMart("ensembl",dataset = "hsapiens_gene_ensembl")
ensembl_ids <- file_c_a$Name

anot <- getBM(
    attributes = c("ensembl_gene_id" , "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = ensembl
)

exp_anot_c_a <- merge(file_c_a, anot ,by.x = "Name",by.y = "ensembl_gene_id",all.x = TRUE)
exp_anot_c_a$Name <- sub("\\..*","",exp_anot_c_a$Name)
#write.csv(exp_anot_amyg,"~/project/brain_map/csv_files/exp_anot_amyg.csv", row.names = FALSE)
#amg <- read.csv("~/project/brain_map/csv_files/exp_anot_amyg.csv")
#library(dplyr)
exp_anot_c_a$anterior_cortex <- rowMeans(exp_anot_c_a[,-c(1,2)],na.rm = TRUE)
mean_tpm_c_a <- exp_anot_c_a[,c(1,2, ncol(exp_anot_c_a))]
#write.csv(mean_tpm_amg,"~/project/brain_map/mean_tpm/mean_tpm_amg_1.csv",row.names = FALSE)
colnames(mean_tpm_c_a)[1] <- "gene_id"
colnames(mean_tpm_c_a)[2] <- "gene_symbol"
#mean tpms are ready !
tissue_df <- list(
    amygdal = mean_tpm_amg,
    hippocampus = mean_tpm_hipc,
    hypotalamus = mean_tpm_hypo,
    frontal_cortex = mean_tpm_c_f,
    anterior_cortex = mean_tpm_c_a
)

for ( tissue in names(tissue_df)){
    colnames(tissue_df[[tissue]]) [1:3] <- c("gene_id","gene_symbol", tissue)
}
merged_df <- Reduce(function(x, y) merge(x,y,by= c("gene_id","gene_symbol"),all = TRUE) , tissue_df)
# filtering unanotated genes 
filter_df <- merged_df[merged_df$gene_symbol != merged_df$gene_id,]
#handling duplicated genes with mean value 
library(dplyr)
collapseed_df <- filter_df %>% 
    group_by(gene_symbol) %>%
    summarise(across(where(is.numeric), mean ,na.rm = TRUE))

#save the matrix- log transfor also
collapsed_matrix <- as.data.frame(collapseed_df)
rownames(collapsed_matrix) <-collapsed_matrix$gene_symbol
collapsed_matrix$gene_symbol <- NULL
log_matrix <-log2(collapsed_matrix + 1)
saveRDS(log_matrix,file = "exp.mat.rds")
