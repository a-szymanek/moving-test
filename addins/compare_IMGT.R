# Compare IMGT sequences and Ensembl sequences.
# Required input info: 
#  * Ensembl file in a fasta format
#  * gene name with a star at the end eg: IGHV2-5* or chain type and region eg: "IGKV"

# Run script in a console eg: Rscript compare_IMGT.R NO ../data/2018_07_31_biomart_IG_nucleotide.txt IGHV2-5*

# An output table includes following columns: 
#  first and second - names of comparing sequences
#  result (logical) TRUE if sequences are identical, otherwise FALSE
#  diff  (character) nucleotides/amino acids differentiate comparing sequences
#  positions - positions of differences
#  len_x - orginal length of the first sequence (not relative to differences positions)
#  len_y - orginal length of the second sequence (not relative to differences positions)
#  nb_diff - number of different positions bewteen first and second sequence
#  add_into - additional information about comparing sequences

args = commandArgs(trailingOnly = TRUE)

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
  dir.create("../data", showWarnings = F)
  if (component_nb_ensembl == 0) {
    type_seq <- "nucleotide" 
    path_to_IMGT <- paste0("../data/", Sys.Date(), "_IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP", Sys.Date(), "_filter")
  } else if (component_nb_ensembl > 0) {
    type_seq <- "protein" 
    path_to_IMGT <- paste0("../data/", Sys.Date(), "_IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F+ORF+inframeP", Sys.Date(), "_filter")
  }
  check_if_exist <- file.exists(path_to_IMGT)
  if (check_if_exist == F) {
    path_to_IMGT <- download_IMGT(type_seq, out_file = "../data")
    IMGT_data <- filter_IMGT(input = path_to_IMGT, chain_type = c("heavy", "light"), region = c("V", "D", "J"))
  } else {
    IMGT_data <- readIgFasta(path_to_IMGT, strip_down_name = F)
  }
} else {
  if (!file.exists(path_to_IMGT)) {
    stop("IMGT file not exist /or name in incorrect")
  } else {
    IMGT_data <- readIgFasta(path_to_IMGT, strip_down_name = F)
  }
}

# component_nb_imgt <- IMGT_data %>% strsplit("") %>% unlist(use.names = F) %>% table() %>% names()  %>% grep("[ACTGN/.]",., invert = T) %>% length()
# 
# if (component_nb_ensembl != component_nb_imgt) {
#   stop("You input data in different types nucleotide/protein")
# }

# ============================================= filter genes ===================================================================
IMGT_data <- IMGT_data[grep(gene_to_compare, names(IMGT_data), value = T, fixed = T)]
if (grepl("*", gene_to_compare, fixed = T)) {
  gene_to_compare <- gsub("*", "", gene_to_compare, fixed = T)
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
sequences <- msaConvert(sequences, type = "seqinr::alignment")
names_seq <- sequences$nam
sequences <- sequences$seq
names(sequences) <- names_seq

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
# diff in length
out_tbl$add_info[strsplit(out_tbl$diff, ",") %>% lapply(., grep, pattern = "[A-Z]") %>% lapply(., function(x) length(x) == 0 ) %>% unlist()] <- "identical/or only different in length"

# most similar to
df <- out_tbl[(grepl("^ENSG", out_tbl$first) & grepl("^[^ENSG]", out_tbl$second)) | (grepl("^ENSG", out_tbl$second) & grepl("^[^ENSG]", out_tbl$first)), ]
min_diff_values <- tapply(df$nb_diff, df$first, min)

selected_most_similar <- sapply(1:length(min_diff_values), function(x) {
  (out_tbl$nb_diff %in% min_diff_values[x] & out_tbl$first %in% names(min_diff_values[x]))
})

selected_most_similar <- ifelse(rowSums(selected_most_similar) == 1, T, F)

out_tbl$add_info[selected_most_similar] <- paste(out_tbl$add_info[selected_most_similar], "the most similar sequences", sep = "; ") %>% gsub("NA;", "", .)

#write table
out_file_path <- file.path("../results", paste("ensembl_IMGT", gene_to_compare, type_seq, Sys.Date(), ".csv", sep = "_"))
write.table(out_tbl, out_file_path, row.names = F, quote = F, sep = "\t")
