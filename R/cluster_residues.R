#' What a cluster carries at one alignment column
#'
#' A tree tip is a cluster representative standing for up to tens of thousands of
#' members, so after a Position search its own residue says no more about the cluster
#' than its own host does. This is the residue equivalent of `cluster_hosts`, and it is
#' formatted to match: shares as whole percents, "; " between them, most common first.
#'
#' The denominator is stated rather than left implicit, because it is not the cluster
#' size. A member gapped at the column carries no residue there and is not counted, so
#' the tally can cover fewer sequences than the cluster holds - and the pop-up shows the
#' cluster size a line above, where the difference would otherwise look like an error.
#'
#' No commas anywhere in the output. Taxonium parses the metadata CSV with a bare
#' split(",") that ignores quoting, so one comma would shift every following column and
#' the pop-up would show the wrong field values with no error at all.
#'
#' @param product product key, "seg1" or "seg7_M2"
#' @param column alignment column, as resolved by gapless_to_consensus()
#' @param residues the per-cluster tally, defaulting to the global cluster_residues
#' @return tibble(representative, cluster_residues), or NULL when no tally is loaded
#' @export
cluster_residue_summary <- function(product, column, residues = cluster_residues) {
  if(is.null(residues) || is.null(column) || is.na(column)) return(NULL)

  wanted_product <- product
  wanted_column  <- as.integer(column)

  residues %>%
    dplyr::filter(product == wanted_product, column == wanted_column) %>%
    dplyr::arrange(representative, dplyr::desc(n), amino_acid) %>% # ties resolved by residue, so the string is stable
    dplyr::summarise(
      cluster_residues = str_c(str_c(amino_acid, " ", residue_share(n, sum(n))),
                               collapse = "; "),
      covered = sum(n),
      .by = representative) %>%
    dplyr::mutate(cluster_residues = str_c(cluster_residues, " (n = ", covered, ")")) %>%
    dplyr::select(representative, cluster_residues)
}

#' Whole-percent shares, bounded away from 0% and 100%
#'
#' Anything present rounds to at least "<1%", as cluster_host_composition.R does for
#' hosts: a residue carried by one sequence in fifty thousand should not read as "0%".
#'
#' The other end matters more here than it does for hosts. Clusters are built at 95%
#' identity, so a site is usually near-monomorphic - 55,771 of 55,775 members carrying V
#' is typical - and a plain round() renders that as "V 100%" followed by three other
#' residues, which contradicts itself. Only an unmixed site gets "100%".
#'
#' @export
residue_share <- function(n, total) {
  percent <- 100 * n / total

  dplyr::case_when(
    percent > 0   & percent < 0.5   ~ "<1%",
    percent < 100 & percent > 99.5  ~ ">99%",
    .default = str_c(round(percent), "%"))
}
