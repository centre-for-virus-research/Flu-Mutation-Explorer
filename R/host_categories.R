#' Host categories that never reach a bar
#'
#' Shared by the app and by the pipeline stage that builds the all-sequences counts, so
#' the two modes of the Adaptation Mutations plot filter identically rather than drifting
#' apart as separate lists.
#'
#' "Unknown" and "Environment" are real values assigned by the taxonomy stage - a
#' sequence with no host recorded, and one sampled from the environment rather than from
#' an animal. "Unclassified" is the marker for a host name the taxonomy could not place
#' at all, kept distinct from "Unknown" so an incomplete taxonomy run is visible as
#' itself rather than merged into the sequences that genuinely record no host.
#'
#' None of the three names a taxon, so none belongs on a host axis. They are excluded
#' when the counts are built and again when they are plotted; the second filter is
#' redundant today and is kept so the two cannot drift apart if the first ever changes.
EXCLUDED_HOSTS <- c("Unknown", "Environment", "Unclassified")
