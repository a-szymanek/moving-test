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
  
  return(cat("Downloading", sequence_type, "sequence completetd"))
}

