Chunk_by_2000=function(fasta_name, prefix = "chunk")
{
#Requires seqinr and Biostrings.  Give a fasta file name and get it chunked into 2000 sequences by segment name.
  #Uses a fasta file from a GISAID download, where the segment name is the 2nd "|" field of each header.
  #fasta_name is the path to that file and prefix names the output files, e.g.
  #Chunk_by_2000("cow24.fasta", "cow2024") writes cow2024_<segment>_<n>.fasta
  #to the current working directory.
  library(seqinr)
  cow2024=read.fasta(fasta_name)
  
library(Biostrings)
cow2024 <- DNAStringSet(sapply(cow2024, function(x) paste(x, collapse = "")))
# Extract the segment name (2nd field between the "|" separators)
segments <- sapply(strsplit(names(cow2024), "|", fixed = TRUE), `[`, 2)

# Check what was found
print(table(segments))

max_per_file <- 2000

for (seg in unique(na.omit(segments))) {
  
  # Sequences belonging to this segment, in original file order
  seg_idx <- which(segments == seg)
  
  # Assign each sequence to a chunk of up to 2000
  chunk_id <- ceiling(seq_along(seg_idx) / max_per_file)
  
  for (i in unique(chunk_id)) {
    chunk_seqs <- cow2024[seg_idx[chunk_id == i]]
    out_file   <- paste0(prefix, "_", seg, "_", i, ".fasta")
    
    writeXStringSet(chunk_seqs, filepath = out_file)
    message("Wrote ", length(chunk_seqs), " sequences to ", out_file)
  }
}
}