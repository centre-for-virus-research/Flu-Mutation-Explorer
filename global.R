library(shiny)
library(shiny.react)
library(jsonlite)
library(tidyverse) # TODO individual packages
library(magrittr)
# library(gt)
library(rBLAST)
library(DT)
library(janitor)
library(ape)
library(Biostrings)
library(pwalign)
library(tools)
library(readxl)
library(ggiraph)
library(bslib)
library(gdtools)
library(config)
library(logger)

# Initialize logger
log_threshold(INFO)
log_layout(layout_glue_generator(format = "{time} [{level}] {msg}"))

# Load configuration
conf <- config::get()

register_gfont("Open Sans")

tree_files <- list.files(path=conf$paths$data$trees,
                         pattern = "\\.nwk$",
                         full.names=TRUE) # %T>% message("Tree files: ", .)

tree_names <-
  tree_files %>%
  basename %>%
  map_vec(tools::file_path_sans_ext) %>%
  map_vec(tools::file_path_sans_ext) %>%
  map_vec(tools::file_path_sans_ext) %>%
  str_replace("_cluster_rep", "")

log_info("Tree names: {paste(tree_names, collapse=', ')}")

# read tree file text into list of strings and set tree names
tree_data <<- 
  tree_files %>% 
  map(read_file) %>% 
  setNames(tree_names) # %T>% message("Trees: ", .) # list of tree names and Newick format strings

# extract tree tip names from Newick format strings as Ape MultiPhylo object
trees_multi <<- 
  tree_data %>% 
  paste(collapse = "") %>%  # concatenate Newick strings
  read.tree(text = ., tree.names = tree_data %>% names) # read Newick strings into Ape MultiPhylo object

full_metadata <<- 
  read_rds(conf$paths$data$metadata) %>% 
  # select(-isolate) %>% # exclude strain and isolate from tree metadata
  select(parsed_strain, primary_accession, h_subtype, n_subtype, host_group, host_order, segment) %>% 
  dplyr::rename(strain = parsed_strain) %>%  # rename parsed_strain column to strain
  mutate(host_group = case_when( # update mammals to human and other
    host_group == "Mammalia" & host_order == "Primates" ~ "Human",
    host_group == "Mammalia" & host_order != "Primates" ~ "Other Mammals",
    host_group == "Aves" ~ "Birds", # rename
    .default = host_group
  )) #%>% 
  # mutate(strain = str_replace_all(strain, " ", "_")) # replace spaces with underscores in strain column

discarded_data <- read_rds(conf$paths$data$discarded) # amino acid sequence data with discarded gaps
names(discarded_data) <- tree_names
discarded <<- discarded_data

transposed <<- read_rds(conf$paths$data$transposed) # amino acid sequence data with gapless positions

# Pre-process adaptation mutations data into a list
adaptation_mutations_all <<- list(
  seg1 = read_csv(conf$paths$adaptation_mutations$seg1, show_col_types = FALSE) %>% 
    select(-1) %>% clean_names() %>% 
    dplyr::select(-segment, -further_notes, -experimentally_verified) %>% 
    dplyr::rename(doi = x, new_amino_acid = mutant),
  
  seg2 = read_csv(conf$paths$adaptation_mutations$seg2, show_col_types = FALSE) %>% 
    select(-1) %>% clean_names() %>% 
    dplyr::select(-segment, -further_notes, -experimentally_verified, -starts_with("x")) %>% 
    dplyr::rename(new_amino_acid = mutant),
  
  seg3 = read_csv(conf$paths$adaptation_mutations$seg3, show_col_types = FALSE) %>% 
    select(-1) %>% clean_names() %>% 
    dplyr::select(-segment, -extra_notes, -experimentally_verified) %>% 
    dplyr::rename(new_amino_acid = mutant),
  
  seg4 = read_csv(conf$paths$adaptation_mutations$seg4, show_col_types = FALSE) %>% 
    select(-1) %>% clean_names() %>% 
    mutate(position_h5 = as.integer(str_extract(mutation_h5_numbering, "-?\\d+")), .after = position) %>% 
    dplyr::select(-segment, -experimentally_verified) %>% 
    dplyr::rename(new_amino_acid = mutation, mutation = mutation_h3_numbering),
    
  seg5 = read_csv(conf$paths$adaptation_mutations$seg5, show_col_types = FALSE) %>% 
    select(-1) %>% clean_names() %>% 
    dplyr::select(-segment, -extra_notes, -experimentally_verified) %>% 
    dplyr::rename(new_amino_acid = mutant),
    
  seg6 = read_csv(conf$paths$adaptation_mutations$seg6, show_col_types = FALSE) %>% 
    select(-1) %>% clean_names() %>% 
    dplyr::select(-segment, -x, -experimentally_verified) %>% 
    dplyr::rename(new_amino_acid = mutation_1),
    
  seg7 = read_csv(conf$paths$adaptation_mutations$seg7, show_col_types = FALSE) %>% 
    select(-1) %>% clean_names() %>% 
    dplyr::select(-segment, -extra_notes, -experimentally_verified) %>% 
    dplyr::rename(new_amino_acid = mutant),
    
  seg8 = read_csv(conf$paths$adaptation_mutations$seg8, show_col_types = FALSE) %>% 
    select(-1) %>% clean_names() %>% 
    dplyr::select(-segment, -extra_notes, -experimentally_verified) %>% 
    dplyr::rename(new_amino_acid = mutant)
)

# Pre-calculate initial metadata CSV string for faster session startup
initial_metadata_csv <<- 
  full_metadata %>% 
  select(-primary_accession, -segment) %>% 
  format_csv()

# read reference sequences for each segment and subtype
# pre-sorted by segment, H and N numbers 
# TODO resolve NA subtypes
cluster_glue <<- 
  read_rds("cluster_reference.rds") %>% 
  drop_na 

# Influenza A virus (A/turkey/England/50-92/91(H5N1)) reference sequences
ref_seqs <<- 
  list.files(pattern = "^seg\\d\\.fasta$") %>%  # reference sequence files 
  purrr::set_names() %>% # set names to file paths
  map(function(x) {
    readAAStringSet(x)   # get reference sequence (singles sequence in each file) %>% 
  }) %>% 
  purrr::set_names(file_path_sans_ext(names(.))) # set names removing file extension

ref_set <<- read_tsv("BLAST_segment_recognizer/ref_set.tsv", col_names = c("accession_version", "segment")) # reference accessions for each segment
status <<- 
  read_tsv("current_version/IAV_DB_summary.log", col_names = FALSE, n_max = 2) #%>%  # status of each tree file
  #select(-X3)

status_segments <<- 
  read_tsv("current_version/IAV_DB_summary.log", col_names = TRUE, skip = 2, 
           col_types = list(segment = col_character(),
                            total = col_integer(),
                            clustered = col_integer())) %>%  # status of each tree file 
  mutate(segment = case_when(
    segment == "Segment_1" ~ "PB2",
    segment == "Segment_2" ~ "PB1",
    segment == "Segment_3" ~ "PA",
    segment == "Segment_4" ~ "HA",
    segment == "Segment_5" ~ "NP",
    segment == "Segment_6" ~ "NA",
    segment == "Segment_7" ~ "M",
    segment == "Segment_8" ~ "NS",
    TRUE ~ segment
  )) %>% 
  rename_with(~ str_to_title(.), everything()) # rename columns to title case

# Source all files in R directory
sapply(list.files("R", full.names = TRUE), source)

# Run validation checks
tryCatch({
  validate_data()
}, error = function(e) {
  log_error("Validation failed: {e$message}")
})
