# August 01st; 2018
# Script for filtering data from IMGT database
require(tigger, quietly = T,  warn.conflicts = F)
require(dplyr, quietly = T,  warn.conflicts = F)
require(alakazam, quietly = T,  warn.conflicts = F)
require(seqinr, quietly = T,  warn.conflicts = F)

#' @description Function for filtering raw IMGT data (selecting organism, chain type, and region)
#' @param input - path to raw IMGT data
#' @param organism - strig; default "Homosapiens"
#' @param chain_type - string; heavy/light/kappa/lambda or both
#' @param region - select only V, D or J region (you can select also more than one region at once)
filter_IMGT <- function(input, organism = "Homosapiens", chain_type, region) {
  ref_seq <- readIgFasta(input, strip_down_name = F)
  #select organism
  ref_seq <- ref_seq[grep(organism, names(ref_seq), ignore.case = T)]
  chains <- c()
  if ("heavy" %in% chain_type | "h" %in% chain_type | "H" %in% chain_type) {
    chains <- c(chains, "H")
  } 
  if ("light" %in% chain_type) {
    chains <- c(chains, "L", "K")
  } 
  if ("kappa" %in% chain_type & !("K" %in% chains) | "K" %in% chain_type & !("K" %in% chains)) {
    chains <- c(chains, "K")
  }
  if ("lambda" %in% chain_type & !("L" %in% chains) | "L" %in% chain_type & !("L" %in% chains)) {
    chains <- c(chains, "L") 
  }
  if (("L" %in% chains | "K" %in% chains) & !("H" %in% chains) & "D" %in% region) {
    cat("D region not exist in light chains! Selected only V and J regions.")
    region <- region[region != "D"]
  }
  #select chains and regions
  query <- outer(paste0("IG", chains), region, FUN = "paste0")  %>% paste(., collapse = "|")
  ref_seq <- ref_seq[grep(query, names(ref_seq), value = T, ignore.case = T)]
  #convert headers into only allele names
  conv_names <- strsplit(names(ref_seq), "[|]")
  names(ref_seq) <- sapply(conv_names, "[", 2)
  write.fasta(sequences = as.list(ref_seq), names = names(ref_seq), file.out = paste(input, Sys.Date(), "filter", sep = "_"))
  return(ref_seq)
}