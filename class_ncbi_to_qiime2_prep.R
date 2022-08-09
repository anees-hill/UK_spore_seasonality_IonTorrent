# Check packages are installed
list.of.packages <- c("tidyverse", "taxizedb")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(tidyverse)
library(taxizedb)

# Note: taxizedb requires a local database of NCBI entries. This can be downloaded with this function
# db_download_ncbi()

# filepath_refseqs
input_filepath <- commandArgs(TRUE)

op_f1 <- as.data.frame(str_locate_all(input_filepath, "/"))
op_f2 <- op_f1[nrow(op_f1),1]
output_filepath <- str_sub(input_filepath, op_f2+1, nchar(input_filepath))

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
    slice(1)
  # adding id col
  input$id <- seq(1:nrow(input))
  # find double named entries
  doubles_ind <- grepl(";", input$X9, fixed = TRUE)
  doubles <- input[doubles_ind,]
  singles <- input[!doubles_ind,]
  # correct double named entries
  doubles <- doubles %>% 
    mutate(X9 = trimws(str_sub(X9, 
                                      start = 1, 
                                      end = (str_locate(X9, ";")[, 1])-1)))
  input_corrected <- bind_rows(singles, doubles)
  input_corrected <- input_corrected %>% 
    arrange(id) %>% select(-id)
  # retrieve taxonomy
  print("[1/2] Retrieving raw taxonomy from local SQLite database...")
  output_list <- classification(input_corrected$X9, db = "ncbi")
  print("[2/2] Parsing format...")
  output <- map(1:length(output_list), retrieve_single_entry, output_list_selection = output_list)
  output_tax <- do.call(rbind.data.frame, output); colnames(output_tax) <- "taxonomy"
  output_df <- data.frame("ASV" = input$X1, "taxonomy" = output_tax)
  return(output_df)
}

# Process ========================================================================

input <- read_tsv(input_filepath, col_names = FALSE)[,c(1,9)]
output <- retrieve_full_taxonomy(input)

write_tsv(output, file = paste0(output_filepath,"_parse.tsv"), col_names = FALSE)
























