# @Project: RandomSampling
# @File: null_stoc4.R
# @Author: Kai Ma
# Date: 20/12/2020


null_stoc4 <- function(input_path, output_path, rand = 1000, nworker = 100) {
  require(NST)
  require(tidyfst)
  myfiles <<- Sys.glob(paste0(input_path, "/metacomm_*.csv"))
  # source("C:/Users/KaiMa/Desktop/nullstoc/tNSTmod.R")
  source("./tNSTmod.R")
  res <- lapply(1:length(myfiles), function(i) {
    comm <- fread(myfiles[i])
    comm <- as.data.table(t(comm[, -1]))
    nr <- nrow(comm)
    ecosphere <- tNSTmod(comm = comm, sp.freq = "prop", samp.rich = "fix", swap.method = "not", abundance = "shuffle", rand = rand, nworker = nworker)
    beta_ecosphere <- ecosphere$details$rand.mean
    ST_ecosphere <- ecosphere$index.pair$ST.ij.bray
    
    ecosim <- tNSTmod(comm, sp.freq = "fix", samp.rich = "fix", swap.method = "tswap", abundance = "shuffle", rand = rand, nworker = nworker)
    beta_ecosim <- ecosim$details$rand.mean
    ST_ecosim <- ecosim$index.pair$ST.ij.bray
    
    ecosphere_stegen <- tNSTmod(comm, sp.freq = "prop", samp.rich = "fix", swap.method = "not", abundance = "region", rand = rand, nworker = nworker)
    beta_ecosphere_stegen <- ecosphere_stegen$details$rand.mean
    ST_ecosphere_stegen <- ecosphere_stegen$index.pair$ST.ij.bray
    
    ecosim_stegen <- tNSTmod(comm, sp.freq = "fix", samp.rich = "fix", swap.method = "tswap", abundance = "region", rand = rand, nworker = nworker)
    beta_ecosim_stegen <- ecosim_stegen$details$rand.mean
    ST_ecosim_stegen <- ecosim_stegen$index.pair$ST.ij.bray
    
    beta <- data.table(ecosphere = beta_ecosphere, ecosim = beta_ecosim, ecosphere_stegen = beta_ecosphere_stegen, ecosim_stegen = beta_ecosim_stegen)
    ST <- data.table(ecosphere = ST_ecosphere, ecosim = ST_ecosim, ecosphere_stegen = ST_ecosphere_stegen, ecosim_stegen = ST_ecosim_stegen)
    list(beta, ST)
  })
  beta <- lapply(1:length(res), function(i) {
    res[[i]][[1]]
  }) %>% rbindlist() %>% mutate_dt(Taxa = c(rep('whole', nr*(nr-1)/2), rep('abundant', nr*(nr-1)/2), rep('rare', nr*(nr-1)/2)))
  ST <- lapply(1:length(res), function(i) {
    res[[i]][[2]]
  }) %>% rbindlist() %>% mutate_dt(Taxa = c(rep('whole', nr*(nr-1)/2), rep('abundant', nr*(nr-1)/2), rep('rare', nr*(nr-1)/2)))
  
  if (!dir.exists(output_path)) dir.create(output_path)
  fwrite(beta, paste0(output_path, "/nullstoc_beta.csv"))
  fwrite(ST, paste0(output_path, "/nullstoc_ST.csv"))
}

