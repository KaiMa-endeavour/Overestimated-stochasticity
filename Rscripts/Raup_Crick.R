# @Project: RandomSampling
# @File: Raup_Crick.R
# @Author: Kai Ma
# Date: 22/1/2021


Raup_Crick <- function(input_path, output_path, null_model = "shuffle", reps = 1000, nworker = 100) {
  require(tidyfst)
  require(ecodist)
  
  myfiles <- Sys.glob(paste0(input_path, "/metacomm_*.csv"))
  # source("C:/Users/KaiMa/Desktop/Raup_Crick_Abundance.r")
  source('./Raup_Crick_Abundance.r')
  
  res <- sapply(1:length(myfiles), function(i) {
    comm <- fread(myfiles[i])
    comm <- t(comm[, -1])
    comm <- apply(comm, 2, as.numeric)
    rc <- raup_crick_abundance(comm, plot_names_in_col1 = F, null_model = null_model, reps = reps, nworker = nworker) %>% as.vector()
  })
  colnames(res) <- c('whole', 'abundant', 'rare')
  if (!dir.exists(output_path)) dir.create(output_path)
  fwrite(res, paste0(output_path, paste0("/RC_", null_model, ".csv")))
}
