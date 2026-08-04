# @Project: RandomSampling
# @File: utils.R
# @Author: Kai Ma
# Date: 2026-08-04
# Description: Shared utility functions for community construction scripts.
#   Extracted from cons_seedcomm.R and cons_mockcomm.R to eliminate duplication.


#' Construct metacommunity files from individual community CSVs
#'
#' Joins individual community CSV files (species x sites), partitions into
#' abundant/rare fractions at 80% cumulative abundance, computes pairwise
#' beta diversity, and writes metacomm_*.csv and beta.csv outputs.
#'
#' @param comm_files Character vector. Paths to community CSV files
#'   (each with a "species" column followed by site abundance columns).
#' @param output_path Character. Directory path for output files.
#' @param dist_method Character. Beta diversity metric passed to vegan::vegdist.
#'   Default is "bray".
#' @return Invisible NULL. Called for side effects (writes files).
construct_metacomm <- function(comm_files, output_path, dist_method = "bray") {
  require(tidyfst)
  require(vegan)

  metadatac <- lapply(comm_files, fread, integer64 = "numeric") %>%
    purrr::reduce(full_join_dt, by = "species")
  metadatac[is.na(metadatac)] <- 0
  metadatac <- metadatac[order(rowSums(metadatac[, -1]), decreasing = T), ]

  # whole community beta diversity
  dist_whole <- vegdist(t(metadatac[, -1]), method = dist_method) %>% as.vector()
  fwrite(metadatac, paste0(output_path, "/metacomm_0whole.csv"))

  # partition into abundant / rare at 80% cumulative abundance
  cutoff <- which(cumsum(rowSums(metadatac[, -1]) / sum(metadatac[, -1])) > 0.8)[1]
  metadataAT <- metadatac[1:cutoff, ]
  metadataRT <- metadatac[(cutoff + 1):nrow(metadatac), ]

  dist_AT <- vegdist(t(metadataAT[, -1]), method = dist_method) %>% as.vector()
  fwrite(metadataAT, paste0(output_path, "/metacomm_abundant.csv"))

  dist_RT <- vegdist(t(metadataRT[, -1]), method = dist_method) %>% as.vector()
  fwrite(metadataRT, paste0(output_path, "/metacomm_rare.csv"))

  # write beta diversity summary
  n_whole <- ncol(metadatac) - 1
  n_AT <- ncol(metadataAT) - 1
  n_RT <- ncol(metadataRT) - 1

  dist <- data.table(
    beta = c(dist_whole, dist_AT, dist_RT),
    Taxa = c(
      rep("whole", n_whole * (n_whole - 1) / 2),
      rep("abundant", n_AT * (n_AT - 1) / 2),
      rep("rare", n_RT * (n_RT - 1) / 2)
    )
  )
  fwrite(dist, paste0(output_path, "/beta.csv"))

  invisible(NULL)
}
