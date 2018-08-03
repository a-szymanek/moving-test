# August 01st; 2018
# Script for downloading data from IMGT database

#' @description Function for downloading data from the IMGT database
#' @param sequence_type - nucleotide/protein
#' @param gapped_format - logical, if the sequence have to be in IMGT gapped format
#' @param out_file - name of output directory (default saved in 'data/'; named as in IMGT database)
download_IMGT <- function(sequence_type, gapped_format = T, out_file = "data") {
  if (sequence_type == "nucleotide" & gapped_format == T) {
    file_name <- "IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP"
  } else if (sequence_type == "nucleotide" & gapped_format == F) {
    file_name <- "IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+inframeP"
  } else if (sequence_type == "protein" & gapped_format == T) {
    file_name <- "IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F+ORF+inframeP"
  } else if (sequence_type == "protein" & gapped_format == T) {
    file_name <- "IMGTGENEDB-ReferenceSequences.fasta-AA-WithoutGaps-F+ORF+inframeP"
  }
  download.file(file.path("http://www.imgt.org/download/GENE-DB", file_name),
                destfile = file.path(out_file, paste(Sys.Date(), file_name, sep = "_")),
                method = "wget")
  cat("Downloading", sequence_type, "sequence completetd\n")
  return(file.path(out_file, paste(Sys.Date(), file_name, sep = "_")))
}

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
  if ("heavy" %in% chain_type) {
    chains <- c(chains, "H")
  } else if ("light" %in% chain_type) {
    chains <- c(chains, "L", "K")
  } else if ("kappa" %in% chain_type & !("K" %in% chains)) {
    chains <- c(chains, "K")
  } else if ("lambda" %in% chain_type & !("L" %in% chains)) {
    chains <- c(chains, "K") 
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
  return(ref_seq)
}
