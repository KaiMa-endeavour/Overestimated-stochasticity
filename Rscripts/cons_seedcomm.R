# @Project: RandomSampling
# @File: cons_seedcomm.R
# @Author: Kai Ma
# Date: 24/11/2020
# Ref: May, F., Gerstner, K., McGlinn, D.J., Xiao, X., Chase, J.M., 2018. mobsim: An r package for the simulation and measurement of biodiversity across spatial scales. Methods in Ecology and Evolution 9, 1401–1408. https://doi.org/10.1111/2041-210x.12986


cons_seedcomm <- function(s_pool = 1e4, n_sim = 1e8, sad_type = "lnorm", drift, dispersal, n_comm = 50, dist_method = "bray", output_path) {
  require(mobsim)
  require(tidyfst)
  require(vegan)
  source("./utils.R")

  s_pool <- as.integer(s_pool)
  n_sim <- as.integer(n_sim)

  if (!dir.exists(output_path)) {
    dir.create(output_path)
  }
  if (!dir.exists(paste0(output_path, "/pseudo/"))) {
    dir.create(paste0(output_path, "/pseudo/"))
  }

  otunum <- matrix(NA, 3 * n_comm, 2)
  seed_comm_files <- character(0)
  for (n in 1:n_comm) {
    cat(paste0("Construction of the ", n, "-th(st/nd/rd) seed community ", "(Total: ", n_comm, ") ", date()), "\n")
    pseudo <- sim_sad(s_pool = s_pool, n_sim = n_sim, sad_type = sad_type, sad_coef = list("cv_abund" = 11), fix_s_sim = T)
    pseudo_dt <- rn_col(as.data.frame(pseudo))
    names(pseudo_dt) <- c("species", paste("site", n, sep = ""))
    fwrite(pseudo_dt, paste0(output_path, "/pseudo/", 10 + n, "_seed_0whole.csv"))

    sp_diff <- sample(pseudo_dt$species, round(drift * s_pool))
    sp_constant <- setdiff(pseudo_dt$species, sp_diff)
    seed_comm <- pseudo_dt[match(c(sp_diff, sp_constant), pseudo_dt$species), ]
    seed_comm$species <- paste0("species", seq(1, nrow(seed_comm)))

    rename <- sample(seed_comm$species, round(dispersal * s_pool))
    seed_comm$species[match(rename, seed_comm$species)] <- paste0(rename, "_1")
    seed_comm <- seed_comm[order(seed_comm[, -1], decreasing = T), ]
    seed_file <- paste0(output_path, "/", 10 + n, "_seed_0whole.csv")
    fwrite(seed_comm, seed_file)
    seed_comm_files <- c(seed_comm_files, seed_file)

    cutoff_s <- which(cumsum(seed_comm[, 2]) / sum(seed_comm[, 2]) > 0.8)[1]
    # abundant
    seed_abundant <- seed_comm[1:cutoff_s, ]
    fwrite(seed_abundant, paste0(output_path, "/", 10 + n, "_seed_AT.csv"))
    # rare
    seed_rare <- seed_comm[(cutoff_s + 1):nrow(seed_comm), ]
    fwrite(seed_rare, paste0(output_path, "/", 10 + n, "_seed_RT.csv"))
    ## result
    otunum[3 * n - 2, 1] <- nrow(seed_comm)
    otunum[3 * n - 2, 2] <- "whole"
    otunum[3 * n - 1, 1] <- nrow(seed_abundant)
    otunum[3 * n - 1, 2] <- "abundant"
    otunum[3 * n, 1] <- nrow(seed_rare)
    otunum[3 * n, 2] <- "rare"
  }
  colnames(otunum) <- c("OTUnum", "Taxa")
  fwrite(as.data.table(otunum), paste0(output_path, "/OTUnum.csv"))

  # construct metacommunity files (shared logic from utils.R)
  construct_metacomm(seed_comm_files, output_path, dist_method)
}
