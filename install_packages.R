#!/usr/bin/env Rscript
# Install required R packages for Flu-GDB app

cat("Installing CRAN packages...\n")
install.packages(c(
  "shiny",
  "shiny.react",
  "jsonlite",
  "tidyverse",
  "magrittr",
  "DT",
  "janitor",
  "bslib",
  "gdtools",
  "ggiraph",
  "config",
  "logger",
  "ape",
  "rBLAST",
  "tools",
  "readxl"
), repos = "https://cloud.r-project.org")

cat("\nInstalling Bioconductor packages...\n")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")

BiocManager::install(c(
  "Biostrings",
  "pwalign"
))

cat("\nPackage installation complete!\n")
