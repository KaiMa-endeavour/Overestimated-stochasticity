# @Project: RandomSampling
# @File: cons_mockcomm.R
# @Author: Kai Ma
# Date: 25/11/2020


cons_mockcomm <- function(input_path, seq_dep, dist_method = "bray", output_path) {
  require(vegan)
  require(tidyfst)
  library(MicroEcoTk)
  source("./utils.R")

  if (!is.null(output_path)) {
    if (!dir.exists(output_path)) {
      dir.create(output_path)
    }
  }
  seed_comm_file <- Sys.glob(paste0(input_path, "/*_seed_0whole.csv"))
  otunum <- matrix(NA, 3 * length(seed_comm_file), 2)

  filename_whole <- dir(input_path, pattern = "/*_seed_0whole.csv") %>%
    strsplit("_") %>%
    as.data.table() %>%
    t()
  filename_whole[, 2] <- "mock"

  filename_AT <- dir(input_path, pattern = "/*_seed_AT.csv") %>%
    strsplit("_") %>%
    as.data.table() %>%
    t()
  filename_AT[, 2] <- "mock"

  filename_RT <- dir(input_path, pattern = "/*_seed_RT.csv") %>%
    strsplit("_") %>%
    as.data.table() %>%
    t()
  filename_RT[, 2] <- "mock"

  mock_comm_files <- character(0)
  for (i in 1:length(seed_comm_file)) {
    cat(paste0("Construct ", i, "-th(st/nd/rd) mock community, sequencing depth = ", formatC(seq_dep, format = "d", big.mark = ","), "."), "\n")
    seed_comm <- fread(seed_comm_file[i])
    mock_comm <- data.table(seed_comm$species, MicroEcoTk::rarefy_vt(t(seed_comm[, -1]), depth = seq_dep, replace = FALSE))
    colnames(mock_comm) <- colnames(seed_comm)

    mock_comm <- mock_comm[which(mock_comm[, 2] > 0), ]
    mock_comm <- mock_comm[order(mock_comm[, -1], decreasing = T), ]
    mock_file <- paste0(output_path, "/", filename_whole[i, ] %>% paste0(collapse = "_"))
    fwrite(mock_comm, mock_file)
    mock_comm_files <- c(mock_comm_files, mock_file)

    cutoff_m <- which(cumsum(mock_comm[, 2]) / sum(mock_comm[, 2]) > 0.8)[1]
    # abundant
    mock_abundant <- mock_comm[1:cutoff_m, ]
    fwrite(mock_abundant, paste0(output_path, "/", filename_AT[i, ] %>% paste0(collapse = "_")))
    # rare
    mock_rare <- mock_comm[(cutoff_m + 1):nrow(mock_comm), ]
    fwrite(mock_rare, paste0(output_path, "/", filename_RT[i, ] %>% paste0(collapse = "_")))
    ## result
    otunum[3 * i - 2, 1] <- nrow(mock_comm)
    otunum[3 * i - 2, 2] <- "whole"
    otunum[3 * i - 1, 1] <- nrow(mock_abundant)
    otunum[3 * i - 1, 2] <- "abundant"
    otunum[3 * i, 1] <- nrow(mock_rare)
    otunum[3 * i, 2] <- "rare"
  }
  colnames(otunum) <- c("OTUnum", "Taxa")
  fwrite(as.data.table(otunum), paste0(output_path, "/OTUnum.csv"))

  # construct metacommunity files (shared logic from utils.R)
  construct_metacomm(mock_comm_files, output_path, dist_method)
}
