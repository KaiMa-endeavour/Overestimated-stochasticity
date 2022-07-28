# @Project: RandomSampling
# @File: cons_mockcomm.R
# @Author: Kai Ma
# Date: 2020/11/25 10:02


# input_path <- 'C:/Users/KaiMa/Desktop/111'
# output_path <- 'C:/Users/KaiMa/Desktop/222'
# seq_dep <- 30000
cons_mockcomm <- function(input_path, seq_dep, dist_method = "bray", output_path) {
  require(vegan)
  require(tidyfst)

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
  rarefy_vt <- function(x, depth, prob = NULL, replace = F) {
    depth <- round(as.numeric(depth), 0)
    stopifnot(is.numeric(x))
    if (depth > sum(x) & replace == F) {
      replace <- T
      warning("replace = TRUE")
    }
    if (is.null(prob)) {
      prob <- prob
    } else {
      prob <- prob[match(x, sort(x, decreasing = T))]
      prob <- rep(prob, x)
    }
    y <- sample(x = rep(1:length(x), x), size = depth, prob = prob, replace = replace)
    y_tab <- table(y)
    z <- numeric(length(x))
    z[as.numeric(names(y_tab))] <- y_tab
    z
  }

  for (i in 1:length(seed_comm_file)) {
    cat(paste0("Construct ", i, "-th(st/nd/rd) mock community, sequencing depth = ", formatC(seq_dep, format = "d", big.mark = ","), "."), "\n")
    seed_comm <- fread(seed_comm_file[i])
    mock_comm <- data.table(seed_comm$species, rarefy_vt(t(seed_comm[, -1]), depth = seq_dep, replace = F))
    colnames(mock_comm) <- colnames(seed_comm)

    mock_comm <- mock_comm[which(mock_comm[, 2] > 0), ]
    mock_comm <- mock_comm[order(mock_comm[, -1], decreasing = T), ]
    fwrite(mock_comm, paste0(output_path, "/", filename_whole[i, ] %>% paste0(collapse = "_")))

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

  metadatac <- list.files(output_path, "/*_mock_0whole.csv", full.names = T) %>%
    lapply(fread, integer64 = "numeric") %>%
    purrr::reduce(full_join_dt, by = "species")
  metadatac[is.na(metadatac)] <- 0
  metadatac <- metadatac[order(rowSums(metadatac[, -1]), decreasing = T), ]
  (dist_whole <- vegdist(t(metadatac[, -1]), method = dist_method) %>% as.vector())
  fwrite(metadatac, paste0(output_path, "/metacomm_0whole.csv"))

  cutoff_m <- which(cumsum(rowSums(metadatac[, -1]) / sum(metadatac[, -1])) > 0.8)[1]
  metadataAT <- metadatac[1:cutoff_m, ]
  metadataRT <- metadatac[(cutoff_m + 1):nrow(metadatac), ]

  (dist_AT <- vegdist(t(metadataAT[, -1]), method = dist_method) %>% as.vector())
  fwrite(metadataAT, paste0(output_path, "/metacomm_abundant.csv"))

  (dist_RT <- vegdist(t(metadataRT[, -1]), method = dist_method) %>% as.vector())
  fwrite(metadataRT, paste0(output_path, "/metacomm_rare.csv"))

  dist <- data.table(beta = c(dist_whole, dist_AT, dist_RT), Taxa = c(rep("whole", (ncol(metadatac) - 1) * (ncol(metadatac) - 2) / 2), rep("abundant", (ncol(metadataAT) - 1) * (ncol(metadataAT) - 2) / 2), rep("rare", (ncol(metadataRT) - 1) * (ncol(metadataRT) - 2) / 2)))
  fwrite(dist, paste0(output_path, "/beta.csv"))
}
# cons_mockcomm(input_path = input_path, seq_dep = seq_dep, output_path = output_path)
