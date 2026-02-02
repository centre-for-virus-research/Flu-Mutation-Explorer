# function to take a sequence name and a gapless position for that individual sequence
# and return a list of amino acid at the consensus position
# TODO put function in separate file
consensus <- function(segment, sequence_id, gapless_position) {
  message("Segment: ", segment)
  message("Sequence ID: ", sequence_id)
  message("Gapless position: ", gapless_position)
  
  # get consensus position from gapless sequence position
  consensus_position <- gapless_to_consensus(segment, sequence_id, gapless_position)
  
  message("Consensus position: ", consensus_position)
  if(is.na(consensus_position)) {
    return(NULL) # return null if consensus position is NA
  }
  
  # get amino acid at consensus position for all sequences in segment
  transposed %>% 
    pluck(segment, consensus_position) %>% 
    unlist %>% # convert to vector
    enframe(name = "Node", value = "amino_acid")
}

# Function to convert a gapless position to a consensus position
gapless_to_consensus <- function(segment, sequence_id, gapless_position){
  # Check if sequence exists for this segment
  if (is.null(discarded[[segment]]) || !sequence_id %in% names(discarded[[segment]])) {
    return(NA_integer_)
  }

  seq_data <- discarded[[segment]][[sequence_id]]
  
  if (gapless_position > length(seq_data) || gapless_position < 1) {
    return(NA_integer_)
  }

  consensus_position <- 
    seq_data[gapless_position] %>% # get amino acid at gapless position
    names %>% # get consensus position as a string
    as.integer # convert to index
  
  if (length(consensus_position) == 0) {
    return(NA_integer_)
  }
  
  consensus_position
}

# Function to convert multiple papers and DOIs to HTML links with line breaks
papers_dois <- function(paper, doi) {
  # Check if paper and doi are not empty
  if (!is.na(paper) & 
      !is.na(doi) & 
      is.character(paper) & 
      is.character(doi) &
      length(paper) != 0 & 
      length(doi) != 0) { 
    paper_list <- str_split_1(paper, ";") %>% str_trim()
    doi_list <- str_split_1(doi, ";") %>% str_trim()
    
    map2_chr(paper_list, doi_list, function(paper_s, doi){
      str_c("<a href='https://doi.org/", doi, "'target='_blank'>", paper_s,"</a>")
    }) %>% str_flatten(collapse = "<br>")
  } else {
    # return empty string if paper or doi is NA or empty
    ""
  }
}
