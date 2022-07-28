# @Project: RandomSampling
# @File: cons_seedcomm.R
# @Author: Kai Ma
# Date: 2020/11/24 21:23
# Ref: May, F., Gerstner, K., McGlinn, D.J., Xiao, X., Chase, J.M., 2018. mobsim: An r package for the simulation and measurement of biodiversity across spatial scales. Methods in Ecology and Evolution 9, 1401–1408. https://doi.org/10.1111/2041-210x.12986


# drift <- 0.35
# dispersal <- 0.35
# output_path <- 'C:/Users/KaiMa/Desktop/111'
cons_seedcomm <- function(s_pool = 1e4, n_sim = 1e8, sad_type = "lnorm", drift, dispersal, n_comm = 50, dist_method = "bray", output_path) {
  require(mobsim)
  require(tidyfst)
  require(vegan)

  s_pool <- as.integer(s_pool)
  n_sim <- as.integer(n_sim)

  if (!dir.exists(output_path)) {
    dir.create(output_path)
  }
  if (!dir.exists(paste0(output_path, "/pseudo/"))) {
    dir.create(paste0(output_path, "/pseudo/"))
  }

  otunum <- matrix(NA, 3 * n_comm, 2)
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
    fwrite(seed_comm, paste0(output_path, "/", 10 + n, "_seed_0whole.csv"))

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

  metadatac <- list.files(output_path, "*_seed_0whole_*", full.names = T) %>%
    lapply(fread, integer64 = "numeric") %>%
    purrr::reduce(full_join_dt, by = "species")
  metadatac[is.na(metadatac)] <- 0
  metadatac <- metadatac[order(rowSums(metadatac[, -1]), decreasing = T), ]
  (dist_whole <- vegdist(t(metadatac[, -1]), method = dist_method) %>% as.vector())
  if (is.null(output_path)) {
    fwrite(metadatac, "/metacomm_0whole.csv")
  } else {
    fwrite(metadatac, paste0(output_path, "/metacomm_0whole.csv"))
  }

  cutoff_s <- which(cumsum(rowSums(metadatac[, -1]) / sum(metadatac[, -1])) > 0.8)[1]
  metadataAT <- metadatac[1:cutoff_s, ]
  metadataRT <- metadatac[(cutoff_s + 1):nrow(metadatac), ]

  (dist_AT <- vegdist(t(metadataAT[, -1]), method = dist_method) %>% as.vector())
  if (is.null(output_path)) {
    fwrite(metadataAT, "/metacomm_abundant.csv")
  } else {
    fwrite(metadataAT, paste0(output_path, "/metacomm_abundant.csv"))
  }

  (dist_RT <- vegdist(t(metadataRT[, -1]), method = dist_method) %>% as.vector())
  if (is.null(output_path)) {
    fwrite(metadataRT, "/metacomm_rare.csv")
  } else {
    fwrite(metadataRT, paste0(output_path, "/metacomm_rare.csv"))
  }

  dist <- data.table(beta = c(dist_whole, dist_AT, dist_RT), Taxa = c(rep("whole", (ncol(metadatac) - 1) * (ncol(metadatac) - 2) / 2), rep("abundant", (ncol(metadataAT) - 1) * (ncol(metadataAT) - 2) / 2), rep("rare", (ncol(metadataRT) - 1) * (ncol(metadataRT) - 2) / 2)))
  if (is.null(output_path)) {
    fwrite(dist, "/beta.csv")
  } else {
    fwrite(dist, paste0(output_path, "/beta.csv"))
  }
}
# cons_seedcomm(drift = drift, dispersal = dispersal, output_path = output_path, n_comm = 3)
