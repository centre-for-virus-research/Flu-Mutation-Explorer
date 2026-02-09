# Amino Acid Color Palette
# Standardized colors based on chemical properties with ranges for each type.

AA_PALETTE <- c(
  # Hydrophobic (Blues)
  "A" = "#C8D7FF",
  "V" = "#A0C0FF",
  "I" = "#78A8FF",
  "L" = "#5090FF",
  "M" = "#2878FF",
  "F" = "#0060FF",
  "W" = "#0048D1",
  
  # Polar (Greens)
  "S" = "#C8FFC8",
  "T" = "#96FF96",
  "N" = "#64FF64",
  "Q" = "#32CD32",
  
  # Positive Charge (Reds)
  "R" = "#FFB4B4",
  "K" = "#FF6E6E",
  "H" = "#FF0000",
  
  # Negative Charge (Purples/Magentas)
  "D" = "#E6BEFF",
  "E" = "#9B30FF",
  
  # Special Cases
  "C" = "#F08080", # Cysteines (Pink)
  "G" = "#F09048", # Glycines (Orange)
  "P" = "#C0C000", # Prolines (Yellow)
  "Y" = "#15A4A4", # Tyrosines (Teal)
  
  # Gaps and Unknowns (Greys)
  "-" = "#E0E0E0",
  "X" = "#B0B0B0",
  "?" = "#D0D0D0",
  "*" = "#A0A0A0",
  "B" = "#C0C0C0",
  "Z" = "#C0C0C0",
  "J" = "#C0C0C0"
)

# Function to get palette for JS (RGB arrays for Taxonium)
get_aa_palette_js <- function() {
  # Convert hex to RGB list
  palette_rgb <- lapply(AA_PALETTE, function(hex) {
    as.vector(col2rgb(hex))
  })
  jsonlite::toJSON(palette_rgb, auto_unbox = TRUE)
}

# Function to get sorted amino acid levels
get_aa_levels <- function() {
  names(AA_PALETTE)
}
