# @Project: RandomSampling
# @File: tNSTmod.R
# @Author: Daliang Ning
# Date: 18/12/2020
# Ref: Ning, D., Deng, Y., Tiedje, J.M., Zhou, J., 2019. A general framework for quantitatively assessing ecological stochasticity. Proc. Natl. Acad. Sci. U.S.A. 116, 16892–16898. https://doi.org/10.1073/pnas.1904623116
# @Reviser: Kai Ma


tNSTmod <- function(comm, group = data.frame(x = rep(1, nrow(comm))), meta.group = NULL, meta.com = NULL,
                    dist.method = "bray", abundance.weighted = T, sp.freq, samp.rich, swap.method, abundance,
                    rand = 1000, output.rand = TRUE, nworker = 4,
                    between.group = FALSE, SES = TRUE, RC = TRUE) {
  if (sum(is.na(comm)) > 0) {
    comm[is.na(comm)] <- 0
    warning("NA is not allowed in comm. automatically filled zero.")
  }
  groupck <- function(group) {
    grp.tab <- table(as.vector(group[, 1]))
    if (sum(grp.tab < 2) > 0) {
      stop("some group(s) has only one sample. impossible to perform beta diversity analysis.")
    }
    if (sum(grp.tab < 6) > 0) {
      warning("some groups have less than 6 samples, for which NST can be calculated but not recommened.")
    }
    invisible()
  }
  if (is.null(meta.group)) {
    sampc <- match.name(rn.list = list(comm = comm, group = group))
    comm <- sampc$comm
    group <- sampc$group
    groupck(group)
    if (is.null(meta.com)) {
      meta.abs <- NULL
    } else {
      meta.abs <- rep(0, ncol(comm))
      meta.com <- meta.com[, which(colnames(meta.com) %in% colnames(comm)), drop = FALSE]
      meta.abs[match(colnames(meta.com), colnames(comm))] <- colSums(meta.com)
      names(meta.abs) <- colnames(comm)
    }
  } else {
    sampc <- match.name(rn.list = list(comm = comm, group = group, meta.group = meta.group))
    comm <- sampc$comm
    group <- sampc$group
    groupck(group)
    meta.group <- sampc$meta.group
    meta.lev <- unique(as.vector(meta.group[, 1]))
    if (is.null(meta.com)) {
      meta.abs <- lapply(1:length(meta.lev), function(i) {
        NULL
      })
    } else {
      meta.com <- meta.com[, which(colnames(meta.com) %in% colnames(comm)), drop = FALSE]
      if (sum(!(rownames(comm) %in% rownames(meta.com))) == 0) {
        meta.abs <- lapply(
          1:length(meta.lev),
          function(i) {
            sampi <- rownames(meta.group)[which(meta.group[, 1] == meta.lev[i])]
            outi <- rep(0, ncol(comm))
            outi[match(colnames(meta.com), colnames(comm))] <- colSums(meta.com[which(rownames(meta.com) %in% sampi), , drop = FALSE])
            names(outi) <- colnames(comm)
            outi
          }
        )
      } else {
        if (sum(!(meta.lev %in% rownames(meta.com))) == 0) {
          meta.abs <- lapply(
            1:length(meta.lev),
            function(i) {
              outi <- rep(0, ncol(comm))
              outi[match(colnames(meta.com), colnames(comm))] <- meta.com[which(rownames(meta.com) == meta.lev[i]), ]
              names(outi) <- colnames(comm)
              outi
            }
          )
        } else {
          stop("meta.group and meta.com setting are not fit.")
        }
      }
    }
    names(meta.abs) <- meta.lev
  }

  obs3 <- as.matrix(beta.g(comm,
    dist.method = dist.method,
    abundance.weighted = abundance.weighted,
    as.3col = TRUE, out.list = FALSE
  ))
  obs <- as.numeric(obs3[, 3])

  # # Null models definition
  # null.models=data.frame(sp.freq=c("equip", "equip", "equip",
  #                                  "prop", "prop", "prop",
  #                                  "fix", "fix", "fix",
  #                                  "fix","fix", "fix", "fix"),
  #                        samp.rich=c("equip", "prop", "fix",
  #                                    "equip", "prop", "fix",
  #                                    "equip", "prop", "fix",
  #                                    "fix", "fix", "fix", "fix"),
  #                        swap.method=c("not", "not", "not",
  #                                      "not", "not", "not",
  #                                      "not", "not", "swap",
  #                                      "swap", "tswap", "quasiswap", "backtrack"),
  #                        stringsAsFactors = FALSE)
  # rownames(null.models)=c("EE", "EP", "EF",
  #                         "PE", "PP", "PF",
  #                         "FE", "FP", "FF",
  #                         "FF.swap", "FF.tswap",
  #                         "FF.quasiswap", "FF.backtrack")
  #
  # sp.freq=null.models[null.model,"sp.freq"]
  # samp.rich=null.models[null.model,"samp.rich"]
  # swap.method=null.models[null.model,"swap.method"]

  dist.rand <- function(i, comm, meta.group = meta.group, dist.method = dist.method, abundance.weighted = abundance.weighted, sp.freq = sp.freq, samp.rich = samp.rich, swap.method = swap.method, abundance = abundance, meta.abs = meta.abs) {
    message("Now randomizing i=", i, ". ", date())
    # if(abundance.weighted){abundance.null="region"}else{abundance.null="not"}
    if (is.null(meta.group)) {
      comr <- NST::taxo.null(
        comm = comm, sp.freq = sp.freq, samp.rich = samp.rich,
        swap.method = swap.method, burnin = i,
        abundance = abundance, region.meta = meta.abs
      )
    } else {
      meta.lev <- unique(as.vector(meta.group[, 1]))
      meta.abs <- meta.abs[match(meta.lev, names(meta.abs))]
      comrl <- Reduce(rbind, lapply(
        1:length(meta.lev),
        function(v) {
          comv <- comm[which(meta.group[, 1] == meta.lev[v]), , drop = FALSE]
          NST::taxo.null(
            comm = comv, sp.freq = sp.freq, samp.rich = samp.rich,
            swap.method = swap.method, burnin = i,
            abundance = abundance, region.meta = meta.abs[[v]]
          )
        }
      ))
      comr <- comrl[match(rownames(comm), rownames(comrl)), , drop = FALSE]
    }
    out <- as.matrix(beta.g(comr,
      dist.method = dist.method,
      abundance.weighted = abundance.weighted,
      as.3col = TRUE, out.list = FALSE
    ))
    as.numeric(out[, 3])
  }
  comm <- apply(comm, 2, as.numeric)
  if (nworker == 1) {
    dist.ran <- sapply(1:rand, dist.rand, comm, meta.group, dist.method, abundance.weighted, sp.freq, samp.rich, swap.method, abundance, meta.abs)
  } else {
    require(foreach)
    require(parallel)
    c1 <- makeCluster(nworker, type = "PSOCK")
    cat("Now randomizing by parallel computing. Begin at ", date(), ". Please wait...", "\n")
    clusterEvalQ(c1, library(NST))
    dist.ran <- parSapply(c1, 1:rand, dist.rand, comm, meta.group, dist.method, abundance.weighted, sp.freq, samp.rich, swap.method, abundance, meta.abs)
    stopCluster(c1)
  }
  if (length(dim(dist.ran)) == 1) {
    rand.mean <- mean(dist.ran)
  } else {
    rand.mean <- apply(dist.ran, 1, mean, na.rm = TRUE)
  }
  if (SES) {
    rand.sd <- apply(dist.ran, 1, stats::sd, na.rm = TRUE)
    ses <- obs3
    ses[, 3] <- (obs - rand.mean) / rand.sd
    ses[which(obs == rand.mean), 3] <- 0
    colnames(ses)[3] <- paste0("SES.", colnames(obs3)[3])
  }
  if (RC) {
    obsmm <- matrix(obs, nrow = nrow(dist.ran), ncol = ncol(dist.ran))
    a1 <- rowSums(dist.ran > obsmm, na.rm = TRUE) / ncol(dist.ran)
    a2 <- rowSums(dist.ran == obsmm, na.rm = TRUE) / ncol(dist.ran)
    rc <- obs3
    rc[, 3] <- (0.5 - (a1 + (a2 / 2))) * 2
    colnames(rc)[3] <- paste0("RC.", colnames(obs3)[3])
  }

  # upper limit of different beta diversity metrics
  beta.limit <- data.frame(
    Dmax.in = c(
      NA, NA, 1, 1,
      1, 1, 1, 1,
      1, 1, 1, NA,
      1, NA, NA, 1.4954,
      NA, 1, 1, 1
    ),
    Dmax.ab = c(
      NA, NA, 1, 1,
      1, 1, 1, 1,
      1, NA, NA, 1,
      1, NA, 1, NA,
      NA, NA, 1, 1
    ),
    stringsAsFactors = FALSE
  )
  rownames(beta.limit) <- c(
    "manhattan", "euclidean", "canberra", "bray",
    "sorensen", "jaccard", "ruzicka", "kulczynski",
    "gower", "altGower", "mGower", "morisita",
    "horn", "binomial", "chao", "cao",
    "mEuclidean", "mManhattan", "chao.jaccard", "chao.sorensen"
  )

  if (abundance.weighted) {
    Dmax <- beta.limit[dist.method, "Dmax.ab"]
  } else {
    Dmax <- beta.limit[dist.method, "Dmax.in"]
  }

  if (is.na(Dmax)) {
    rand.amax <- apply(dist.ran, 1, max, na.rm = TRUE)
    rand.bmax <- apply(dist.ran, 1, function(v) {
      max(stats::density(v, na.rm = TRUE)$x)
    })
    d.up1 <- (sum((obs <= 1) * (rand.amax <= 1)) == length(obs)) * (rand.bmax > 1)
    rand.max <- d.up1 + ((1 - d.up1) * rand.bmax)
    Dmax <- max(obs, rand.max)
  }

  GD <- rand.mean / cbind(dist.ran, obs) # G / D and G / Gk
  EC <- (Dmax - rand.mean) / (Dmax - cbind(dist.ran, obs)) # E / C and E / Ek
  EC[is.nan(EC)] <- 1
  ECGD <- EC
  ECGD[which(EC > 1)] <- GD[which(EC > 1)]

  Cij <- (Dmax - obs) / Dmax
  Eij <- (Dmax - rand.mean) / Dmax
  Dsij <- obs / Dmax
  Gsij <- rand.mean / Dmax
  MSTij <- (Eij / Cij) * (Dsij / Gsij)
  MSTij[which(MSTij > 1)] <- ((Cij / Eij) * (Gsij / Dsij))[which(MSTij > 1)]
  indexs <- data.frame(obs3[, 1:2],
    D.ij = obs, G.ij = rand.mean, Ds.ij = Dsij,
    Gs.ij = Gsij, C.ij = Cij, E.ij = Eij, ST.ij = ECGD[, ncol(ECGD)], MST.ij = MSTij
  )
  colnames(indexs)[3:ncol(indexs)] <- paste0(colnames(indexs)[3:ncol(indexs)], ".", colnames(obs3)[3])
  STmin <- Eij
  STmin[which(obs > rand.mean)] <- (1 - Eij)[which(obs > rand.mean)]

  grp.lev <- unique(group[, 1])
  outi <- outij <- list()
  if (between.group) {
    outik <- outikj <- list()
    bgi <- 1
  }

  for (i in 1:length(grp.lev)) {
    sampi <- rownames(group)[which(group[, 1] == grp.lev[i])]
    ni <- length(sampi)
    idi <- which((obs3[, 1] %in% sampi) & (obs3[, 2] %in% sampi))
    ECGDi <- ECGD[idi, , drop = FALSE]
    ECGD.mi <- colMeans(ECGDi, na.rm = TRUE)
    STi.max <- max(ECGD.mi, na.rm = TRUE)
    STi.min <- mean(STmin[idi])
    STi <- ECGD.mi[[length(ECGD.mi)]]
    if (STi == STi.min) {
      NSTi <- 0
    } else {
      NSTi <- (STi - STi.min) / (STi.max - STi.min)
    }
    MSTi <- mean(MSTij[idi], na.rm = TRUE)
    outi[[i]] <- data.frame(group = grp.lev[i], size = length(idi), ST.i = STi, NST.i = NSTi, MST.i = MSTi, stringsAsFactors = FALSE)

    STij <- ECGD[idi, ncol(ECGD)]
    NSTij <- (STij - STi.min) / (STi.max - STi.min)
    MSTiji <- MSTij[idi]
    outij[[i]] <- data.frame(
      group = rep(grp.lev[i], length(idi)),
      obs3[idi, 1:2, drop = FALSE], C.ij = Cij[idi], E.ij = Eij[idi],
      ST.ij = STij, NST.ij = NSTij, MST.ij = MSTiji, stringsAsFactors = FALSE
    )
    if (between.group) {
      if (i < length(grp.lev)) {
        for (k in (i + 1):length(grp.lev))
        {
          sampk <- rownames(group)[which(group[, 1] == grp.lev[k])]
          nk <- length(sampk)
          idik <- which(((obs3[, 1] %in% sampi) & (obs3[, 2] %in% sampk)) | ((obs3[, 1] %in% sampk) & (obs3[, 2] %in% sampi)))

          ECGDik <- ECGD[idik, , drop = FALSE]
          ECGD.mik <- colMeans(ECGDik)
          STik.max <- max(ECGD.mik)
          STik.min <- mean(STmin[idik])
          STik <- ECGD.mik[[length(ECGD.mik)]]
          if (STik == STik.min) {
            NSTik <- 0
          } else {
            NSTik <- (STik - STik.min) / (STik.max - STik.min)
          }
          MSTik <- mean(MSTij[idik])
          outik[[bgi]] <- data.frame(
            group = paste0(grp.lev[i], ".vs.", grp.lev[k]), size = length(idik),
            ST.i = STik, NST.i = NSTik, MST.i = MSTik, stringsAsFactors = FALSE
          )

          STikj <- ECGD[idik, ncol(ECGD)]
          NSTikj <- (STikj - STik.min) / (STik.max - STik.min)
          MSTikj <- MSTij[idik]
          outikj[[bgi]] <- data.frame(
            group = rep(paste0(grp.lev[i], ".vs.", grp.lev[k]), length(idik)),
            obs3[idik, 1:2, drop = FALSE], C.ij = Cij[idik], E.ij = Eij[idik],
            ST.ij = STikj, NST.ij = NSTikj, MST.ij = MSTikj, stringsAsFactors = FALSE
          )
          bgi <- bgi + 1
        }
      }
    }
  }
  outi <- Reduce(rbind, outi)
  outij <- Reduce(rbind, outij)
  colnames(outi)[3:ncol(outi)] <- paste0(colnames(outi)[3:ncol(outi)], ".", colnames(obs3)[3])
  colnames(outij)[4:ncol(outij)] <- paste0(colnames(outij)[4:ncol(outij)], ".", colnames(obs3)[3])
  output <- list(index.pair = indexs, index.grp = outi, index.pair.grp = outij, Dmax = Dmax, dist.method = colnames(obs3)[3])
  if (between.group) {
    outik <- Reduce(rbind, outik)
    outikj <- Reduce(rbind, outikj)
    colnames(outik)[3:ncol(outik)] <- paste0(colnames(outik)[3:ncol(outik)], ".", colnames(obs3)[3])
    colnames(outikj)[4:ncol(outikj)] <- paste0(colnames(outikj)[4:ncol(outikj)], ".", colnames(obs3)[3])
    output <- c(output, list(index.between = outik, index.pair.between = outikj))
  }
  if (SES) {
    output$index.pair <- data.frame(output$index.pair, ses[, 3, drop = FALSE], stringsAsFactors = FALSE)
  }
  if (RC) {
    output$index.pair <- data.frame(output$index.pair, rc[, 3, drop = FALSE], stringsAsFactors = FALSE)
  }
  if (output.rand) {
    details <- list(rand.mean = rand.mean, Dmax = Dmax, obs3 = obs3, dist.ran = dist.ran, group = group, meta.group = meta.group)
    output <- c(output, list(details = details))
  }
  output
}
