args = commandArgs(trailingOnly = TRUE)
# Script for preparing comparison bewteen ensemb and IMGT data

# ==================================================  load required packages  ================================================
library(tigger, quietly = T, verbose = F)
library(msa, quietly = T, verbose = F)
library(dplyr, quietly = T, verbose = F, warn.conflicts = F)
source("main/download_IMGT.R")
source("main/filter_IMGT.R")

# =================================================== args ===================================================================
args <- c()
#IMGT_file_1 <- "data/2018-08-01_IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP_2018-08-01"
IMGT_file_1 <- NA
#IMGT_file_1 <- "data/2018-08-01_IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F+ORF+inframeP"
Ensembl_file_1 <- "data/2018_07_31_biomart_IG_protein.txt"
#Ensembl_file_1 <- "data/2018_07_31_biomart_IG_protein.txt"
gene_to_compare <- "IGHV2-5"
args[1] <- IMGT_file_1
args[2] <- Ensembl_file_1
args[3] <- gene_to_compare

if (length(args) != 3 ) {
  stop("You have to provide 3 parameters: 
        1. path to IMGT  data or NA
        2. path to Ensembl nucleotide data (this files should be downloaded yet)
        3. name of the gene for comparing", call. = FALSE)
} 

path_to_IMGT    <- args[1]
path_to_Ensembl <- args[2]
gene_to_compare <- args[3]

# ================================================= read data =================================================================
if (!file.exists(path_to_Ensembl)) {
  stop("File not exist")
}

Ensembl_data <- readIgFasta(path_to_Ensembl, strip_down_name = F)
component_nb_ensembl <- Ensembl_data %>% strsplit("") %>% unlist(use.names = F) %>% table() %>% names() %>% grep("[ACTGN]",., invert = T) %>% length()

if (is.na(path_to_IMGT)) {
  if (component_nb_ensembl == 0 & is.na(path_to_IMGT)) {
    type_seq <- "nucleotide" 
  } else if (component_nb_ensembl > 0 & is.na(path_to_IMGT)) {
    type_seq <- "protein" 
  }
  path_to_IMGT <- download_IMGT(type_seq)
  #IMGT_data <- filter_IMGT(input = path_to_IMGT, chain_type = substr(gene_to_compare, 3,3), region = substr(gene_to_compare, 4, 4))
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
IMGT_data <- IMGT_data[grep(paste0(args[3], "*"), names(IMGT_data), value = T, fixed = T)]
Ensembl_data <- Ensembl_data[grep(args[3], names(Ensembl_data), value = T, fixed = T)]

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

sel <- grep("[A-Z]", strsplit(sequences[length(Ensembl_data) + 1], split = "") %>% unlist()) %>% min()
out <- substr(sequences, sel, nchar(sequences[1]))

# ======================================================= prepare comparison ===================================================
n <- length(sequences) * length(sequences)
comparison <- data.frame(first = rep(NA, times = n), 
                    second = rep(NA, times = n), 
                    result = rep(NA, times = n), 
                    diff = rep(NA, times = n), 
                    positions = rep(NA, times = n), 
                    len_x = rep(NA, times = n), 
                    len_y = rep(NA, times = n), 
                    nb_diff = rep(NA, times = n))
counter <- 1
for (x in 1:length(out)) {
  for (y in 1:length(out)) {
    result <- out[x] == out[y]
    differences <- strsplit(out[x], split = "") %>% unlist() != strsplit(out[y], split = "") %>% unlist()
    nucl <- unlist(strsplit(out[x], split = ""))[differences]
    position <- which(differences == TRUE)
    comparison[counter, ] <- c(names(out[x]), 
                          names(out[y]), 
                          result, 
                          nucl %>% unname() %>% paste(., collapse = ", "), 
                          position %>% unname() %>% paste(., collapse = ", "),
                          nchar(out[x]), 
                          nchar(out[y]), 
                          length(nucl)
    )
    counter <- counter + 1
  }
}

comparison[, 6:8] <- apply(comparison[, 6:8], 2, as.numeric)
comparison <- comparison[comparison$first != comparison$second, ]

out_file_path <- file.path("../results", paste("ensembl_IMGT", gene_to_compare, type_seq, Sys.Date(), sep = "_"))
write.csv(comparison, out_file_path, row.names = F, quote = F)
# ========================================== summary ================================
identical_seq <- comparison[comparison$result == T, c(1,2)]
length_diff <- comparison[strsplit(comparison$diff, ",") %>% lapply(., grep, pattern = "[A-Z]") %>% lapply(., function(x) length(x) == 0 ) %>% unlist(), ]
df <- comparison[(grepl("^ENSG", comparison$first) & grepl("^[^ENSG]", comparison$second)) | (grepl("^ENSG", comparison$second) & grepl("^[^ENSG]", comparison$first)), ]  
df <- df[order(df$nb_diff, decreasing = F), ] %>% head(length(Ensembl_data)*2)

cat("Identical sequences:") 
identical_seq
cat("Different in length:")
length_diff[!(rownames(length_diff) %in% rownames(identical_seq)), ]
cat("Enseble sequences are the most similar to:")
df
