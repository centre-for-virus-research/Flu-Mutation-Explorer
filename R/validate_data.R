#' Validate Data Integrity on Startup
#' 
#' Checks for existence of required files and validates that reference accessions
#' exist in the loaded data.
#'
#' @return NULL or stops on error
#' @export
validate_data <- function() {
  log_info("Starting data validation...")
  
  # 1. File Existence Checks
  # Iterate over all paths in config (recursive)
  # But simple checks first
  required_files <- c(
    conf$paths$data$discarded,
    conf$paths$data$transposed,
    conf$paths$data$metadata
  )
  
  missing_files <- required_files[!file.exists(required_files)]
  if(length(missing_files) > 0) {
    log_warn("Missing core data files: {paste(missing_files, collapse=', ')}")
    # We might want to stop, or allow partial load. Strict mode: Stop.
    stop("Critical data files missing. Check config and file system.")
  }
  
  # 2. Reference ID Integrity
  # Check if Turkey references exist in 'discarded' object
  # adaptationValues$segment logic maps to keys 1..8 in global.R (renamed from seg1..seg8 in config?)
  # Wait, in global.R we renamed discarded keys to tree_names (seg1, seg2...).
  # And config has seg1: ID.
  
  turkey_refs <- conf$references$turkey
  
  # Iterate names of logic (seg1, seg2...)
  for (seg_key in names(turkey_refs)) {
    ref_id <- turkey_refs[[seg_key]]
    
    # Check if segment key exists in discarded
    if (is.null(discarded[[seg_key]])) {
      log_warn("Segment key '{seg_key}' missing from 'discarded' object.")
      next
    }
    
    # Check if sequence ID exists in that segment
    if (!ref_id %in% names(discarded[[seg_key]])) {
      log_error("CRITICAL: Reference ID '{ref_id}' for '{seg_key}' NOT FOUND in sequence data.")
      # Note: We log error but maybe not stop app entirely, or we should?
      # The plan said "warn/error immediately".
      # Given this is the exact bug the user faced, we should make it very visible.
      # User might want to fix it in config. 
    } else {
      log_success("Verified reference '{ref_id}' for '{seg_key}'.")
    }
  }
  
  log_info("Data validation complete.")
}
