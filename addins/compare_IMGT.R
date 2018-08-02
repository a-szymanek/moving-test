args = commandArgs(trailingOnly = TRUE)
# Script for preparing comparison bewteen ensemb and IMGT data

# ==================================================  load required packages  ================================================
library(tigger, quietly = T, verbose = F)
library(doParallel, quietly = T, verbose = F)
library(msa, quietly = T, verbose = F)
library(dplyr, quietly = T, verbose = F, warn.conflicts = F)
library(utils)
source("../main/download_IMGT.R")
source("../main/filter_IMGT.R")
options(stringsAsFactors = F)
# =================================================== args ===================================================================
args <- c()
IMGT_file_1 <- "NO"
Ensembl_file_1 <- "data/2018_07_31_biomart_IG_nucleotide.txt"
gene_to_compare <- "IGLV"
args[1] <- IMGT_file_1
args[2] <- Ensembl_file_1
args[3] <- gene_to_compare

if (length(args) != 3 ) {
  stop("You have to provide 3 parameters: 
        1. path to IMGT  data or NO
        2. path to Ensembl nucleotide data (this files should be downloaded yet)
        3. name of the gene for comparing or gene type eg: IGHV", call. = FALSE)
} 

path_to_IMGT    <- args[1]
path_to_Ensembl <- args[2]
gene_to_compare <- args[3]

# ================================================= read data =================================================================
if (!file.exists(path_to_Ensembl)) {
  stop("Ensembl file not exist")
}

Ensembl_data <- readIgFasta(path_to_Ensembl, strip_down_name = F)
component_nb_ensembl <- Ensembl_data %>% strsplit("") %>% unlist(use.names = F) %>% table() %>% names() %>% grep("[ACTGN]",., invert = T) %>% length()

if (path_to_IMGT %in% "NO") {
  if (component_nb_ensembl == 0) {
    type_seq <- "nucleotide" 
  } else if (component_nb_ensembl > 0) {
    type_seq <- "protein" 
  }
  path_to_IMGT <- download_IMGT(type_seq, out_file = "../data")
  IMGT_data <- filter_IMGT(input = path_to_IMGT, chain_type = c("heavy", "light"), region = c("V", "D", "J"))
} else {
  if (!file.exists(path_to_IMGT)) {
    stop("File not exist")
  } else {
    IMGT_data <- readIgFasta(path_to_IMGT, strip_down_name = F)
  }
}

component_nb_imgt <- IMGT_data %>% strsplit("") %>% unlist(use.names = F) %>% table() %>% names()  %>% grep("[ACTGN/.]",., invert = T) %>% length()

if (component_nb_ensembl != component_nb_imgt) {
  stop("You input data in different types nucleotide/protein")
}

# ============================================= filter genes ===================================================================
IMGT_data <- IMGT_data[grep(gene_to_compare, names(IMGT_data), value = T, fixed = T)]
if (grepl("*", gene_to_compare, fixed = T)) {
  gene_to_compare <- gsub("*", "", gene_to_compare)
}
Ensembl_data <- Ensembl_data[grep(gene_to_compare, names(Ensembl_data), value = T, fixed = T)]

# ========================================================= make MSA ==========================================================
if (component_nb_ensembl == 0) {
  seq_object <- DNAStringSet(c(Ensembl_data, IMGT_data))
  type_seq <- "nucleotide" 
} else {
  seq_object <- AAStringSet(c(Ensembl_data, IMGT_data))
  type_seq <- "protein" 
}

sequences <- msa(seq_object, order = "input")
#print(sequences, show = "complete")
sequences <- msaConvert(sequences, type = "seqinr::alignment")
names_seq <- sequences$nam
sequences <- sequences$seq
names(sequences) <- names_seq

#sel <- grep("[A-Z]", strsplit(sequences[length(Ensembl_data) + 1], split = "") %>% unlist()) %>% min()
#out <- substr(sequences, sel, nchar(sequences[1]))
out <- sequences
# ======================================================= prepare comparison ===================================================
cl <- makeCluster(7)
registerDoParallel(cl)
out_tbl <- foreach(x = 1:length(out), .packages = "dplyr") %dopar% {
  split_x <- strsplit(out[x], split = "") %>% unlist()
  lapply(1:length(out), function(y) { 
    result <- out[x] == out[y]
    differences <- split_x != strsplit(out[y], split = "") %>% unlist()
    nucl <- unlist(strsplit(out[x], split = ""))[differences]
    position <- which(differences == TRUE)
    data.frame(first = names(out[x]), 
               second = names(out[y]), 
               result = result, 
               diff = nucl %>% unname() %>% paste(., collapse = ", "), 
               positions = position %>% unname() %>% paste(., collapse = ", "),
               len_x = nchar(gsub("-", "", out[[x]])), 
               len_y = nchar(gsub("-", "", out[[y]])), 
               nb_diff = length(nucl), row.names = NULL, stringsAsFactors = F
    )
  })
}

out_tbl <- do.call(rbind, unlist(out_tbl, recursive = F))
out_tbl <- out_tbl[out_tbl$first != out_tbl$second, ]

out_file_path <- file.path("../results", paste("ensembl_IMGT", gene_to_compare, type_seq, Sys.Date(), ".csv", sep = "_"))
write.table(out_tbl, out_file_path, row.names = F, quote = F, sep = "\t")
# ========================================== summary ================================
# identical_seq <- out_tbl[out_tbl$result == T, c(1,2)]
# length_diff <- out_tbl[strsplit(out_tbl$diff, ",") %>% lapply(., grep, pattern = "[A-Z]") %>% lapply(., function(x) length(x) == 0 ) %>% unlist(), ]
# df <- out_tbl[(grepl("^ENSG", out_tbl$first) & grepl("^[^ENSG]", out_tbl$second)) | (grepl("^ENSG", out_tbl$second) & grepl("^[^ENSG]", out_tbl$first)), ]  
# df_ <- split.data.frame(df, df$first) %>% lapply(., function(x) x[order(x$nb_diff, decreasing = F), ] %>% head(4)) %>% unname() %>% do.call(rbind, .)
# cat("Identical sequences:") 
# identical_seq
# cat("Different in length:")
# length_diff[!(rownames(length_diff) %in% rownames(identical_seq)), ]
# cat("Enseble sequences are the most similar to:")
# df
