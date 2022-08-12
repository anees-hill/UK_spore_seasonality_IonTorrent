# blastn_to_qiime.R - Samuel Anees-Hill
# The following script is designed for UNIX-based execution and takes blastn (NCBI) output (tsv. In blastn output
# format (outfmt) 6 style, with nident included as a column) and converts it into a .tsv format importable into QIIME2 via the function 'qiime tools import'
# (type 'FeatureData[Taxonomy]'/input-format HeaderlessTSVTaxonomyFormat).
# An optional additional argument (4th argument) can take the fasta file used in the blastn search and locate ASVs/sequences that
# are not present in the BLASTn output (note: blastn outfmt 6 style does not
# include ASVs with no matches - hence the use of this feature). These blank entries are labelled with the 5th argument (i.e. "k__").

# Check packages are installed
list.of.packages <- c("tidyverse", "taxizedb", "Biostrings")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(tidyverse)
library(taxizedb)
library(Biostrings)

# Note: package taxizedb requires a local database of NCBI entries. This should be automatic but can be downloaded with this function:
# db_download_ncbi()

# UNIX inputs
input_postblast <- NULL; nident_col_num <- NULL; output_filename <- NULL; input_preblast <- NULL; missing_replace <- "k__"
input_postblast <- commandArgs(TRUE)[1]
nident_col_num <- as.double(commandArgs(TRUE)[2])
output_filename <- commandArgs(TRUE)[3]
input_preblast <- commandArgs(TRUE)[4]
missing_replace <- commandArgs(TRUE)[5]
args = commandArgs(trailingOnly=TRUE)


# Pre-execution checks
if((is.null(input_postblast)==TRUE)){stop("ERROR: input BLASTn output as the first argument")}
if((is.null(nident_col_num)==TRUE)){stop("ERROR: 'nident' column number not input (second argument)")}
if((is.null(output_filename)==TRUE)){stop("ERROR: output filename missing (third argument)")}

# Define functions:
retrieve_single_entry <- function(entry, output_list_selection){
  output_list <- output_list_selection
  single_output <- pluck(output_list, entry)
  k <- single_output[which(single_output$rank == "kingdom"),"name"][1]
  p <- single_output[which(single_output$rank == "phylum"),"name"][1]
  c <- single_output[which(single_output$rank == "class"),"name"][1]
  o <- single_output[which(single_output$rank == "order"),"name"][1]
  f <- single_output[which(single_output$rank == "family"),"name"][1]
  g <- single_output[which(single_output$rank == "genus"),"name"][1]
  s <- single_output[which(single_output$rank == "species"),"name"][1]
  # Remove genus from species entry
  s_corrected <- substr(s, (unlist(gregexpr(' ', s))[1])+1, nchar(s)) 
  # Concatenate
  total <- paste0("k__",k,"; p__",p,"; c__",c,
                  "; o__",o,"; f__",f,"; g__",g,"; s__",s_corrected)
  return(total)
}

retrieve_full_taxonomy <- function(input_selection){
  # keep top ncbi result only
  input <- input_selection %>% 
    group_by(X1) %>% 
    dplyr::slice(1)
  # adding id col
  input$id <- seq(1:nrow(input))
  # find double named entries
  doubles_ind <- grepl(";", input$nident, fixed = TRUE)
  doubles <- input[doubles_ind,]
  singles <- input[!doubles_ind,]
  # correct double named entries
  doubles <- doubles %>% 
    mutate(nident = trimws(str_sub(nident, 
                                      start = 1, 
                                      end = (str_locate(nident, ";")[, 1])-1)))
  doubles$nident <- as.double(doubles$nident)
  singles$nident <- as.double(singles$nident)
  input_corrected <- bind_rows(singles, doubles)
  input_corrected <- input_corrected %>% 
    arrange(id) %>% select(-id)
  # retrieve taxonomy
  print("[1/2] Retrieving raw taxonomy from local SQLite database...(may take a few minutes)")
  output_list <- classification(input_corrected$nident, db = "ncbi")
  print("[2/2] Parsing format for QIIME2 import")
  output <- map(1:length(output_list), retrieve_single_entry, output_list_selection = output_list)
  output_tax <- do.call(rbind.data.frame, output); colnames(output_tax) <- "taxonomy"
  output_df <- data.frame("ASV" = input$X1, "taxonomy" = output_tax)
  return(output_df)
}

check_for_missing_entries <- function(){
  preblastfasta = readDNAStringSet(input_preblast)
  seq_name = names(preblastfasta)
  sequence = paste(preblastfasta)
  ref <- data.frame(seq_name, sequence)
  pre_blast <- unique(ref$seq_name)
  post_blast <- unique(input$X1)
  check <- pre_blast %in% post_blast
  nonmatches <- pre_blast[!pre_blast %in% post_blast]
  print(paste0("Comparing blastn output with blastn fasta input: ",length(nonmatches),"/",length(pre_blast)," entries had no mactches via BLASTn (",
               (round(length(nonmatches)/length(pre_blast)*100,1)),"%)"))
  print(paste0("Any missing entries will be classified '",missing_replace,"'"))
  return(nonmatches)
}

# Process ========================================================================

# Import postblast results query
input <- read_tsv(input_postblast, col_names = FALSE)[,c(1,nident_col_num)]
colnames(input) <- c("X1","nident")

# Process
if(length(args)>3){
  nonmatches <- check_for_missing_entries()
  } else {
  nonmatches <- NULL
}
output <- retrieve_full_taxonomy(input)
if(((is_character(nonmatches))==TRUE) & (length(args)>3) & length(nonmatches>0)){
  nonmatches <- data.frame("ASV" = nonmatches, "taxonomy" = missing_replace)
  output <- output %>% 
    bind_rows(nonmatches)
}

# Write-out
write_tsv(output, file = output_filename, col_names = FALSE)
print(paste0("File saved as ",output_filename))



