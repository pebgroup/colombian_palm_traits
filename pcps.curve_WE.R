# modified from PCPS

pcps.curve_WE <- function (comm, phylodist, trait, checkdata = TRUE, method = "bray", 
          squareroot = TRUE, ranks = TRUE, null.model.ts = FALSE, null.model.bm = FALSE, 
          tree, runs = 99, progressbar = FALSE, parallel = NULL) 
{
  res <- list(call = match.call())
  if (inherits(comm, "metacommunity.data")) {
    if (!missing(phylodist) | !missing(trait)) {
      stop("\n When you use an object of class metacommunity.data the arguments phylodist and trait must not be specified. \n")
    }
    phylodist <- comm$phylodist
    trait <- comm$traits
    comm <- comm$community
  }
  list.warning <- list()
  if (checkdata) {
    organize.temp <- organize.pcps(comm, phylodist = phylodist, 
                                   trait = trait, check.comm = TRUE)
    if (!is.null(organize.temp$stop)) {
      organize.temp$call <- match.call()
      return(organize.temp)
    }
    list.warning <- organize.temp$list.warning
    comm <- organize.temp$community
    phylodist <- organize.temp$phylodist
    trait <- organize.temp$traits
  }
  if (length(list.warning) > 0) {
    res$list.warning <- list.warning
  }
  if (ncol(trait) != 1) {
    stop("\n Only one trait is allowed\n")
  }
  MT <- SYNCSA::matrix.t(comm, trait, scale = FALSE, ranks = ranks, 
                         notification = FALSE)$matrix.T
  res.pcps <- pcps(comm, phylodist, method = method, squareroot = squareroot, 
                   correlations = FALSE)
  res.values <- res.pcps$values
  res.vectors <- res.pcps$vectors
  curve.obs <- pcpc.curve.calc(res.values, res.vectors, MT)
  rownames(curve.obs) <- rownames(res.values)
  res$curve.obs <- curve.obs
  if (progressbar) {
    if (null.model.ts & null.model.bm) {
      BarRuns <- runs * 2
    }
    else {
      BarRuns <- runs
    }
  }
  newClusters <- FALSE
  if (is.numeric(parallel)) {
    parallel <- parallel::makeCluster(parallel, type = "PSOCK")
    newClusters <- TRUE
  }
  ptest.ts <- function(samp, comm, phylodist, method, squareroot, 
                       mt) {
    phylodist_rand <- phylodist
    rownames(phylodist_rand) <- colnames(phylodist_rand) <- rownames(phylodist)[samp]
    pcps.null <- PCPS::pcps(comm, phylodist_rand, 
                            method = method, squareroot = squareroot, correlations = FALSE)
    res <- PCPS::pcpc.curve.calc(pcps.null$values, pcps.null$vectors, 
                                 mt)
    return(res)
  }
  ptest.bm <- function(samp, tree, comm, values, vectors, ranks) {
    trait.null <- cbind(ape::rTraitCont(phy = tree, model = "BM"))
    match.names <- match(colnames(comm), rownames(trait.null))
    MT.null <- SYNCSA::matrix.t(comm, trait.null[match.names, 
                                                 , drop = FALSE], scale = FALSE, ranks = ranks, notification = FALSE)$matrix.T
    res <- PCPS::pcpc.curve.calc(values, vectors, MT.null)
    return(res)
  }
  if (null.model.ts) {
    seqpermutation <- SYNCSA::permut.vector(ncol(phylodist), 
                                            nset = runs)
    seqpermutation <- lapply(seq_len(nrow(seqpermutation)), 
                             function(i) seqpermutation[i, ])
    if (!inherits(parallel, "cluster")) {
      res.curve.null.ts <- vector("list", runs)
      for (i in 1:runs) {
        res.curve.null.ts[[i]] <- ptest.ts(samp = seqpermutation[[i]], # modified this line 24/9/2021, WLE  
                                           comm = comm, phylodist = phylodist, method = method, 
                                           squareroot = squareroot, mt = MT)
        if (progressbar) {
          SYNCSA::ProgressBAR(i, BarRuns, style = 3)
        }
      }
    }
    else {
      res.curve.null.ts <- parallel::clusterApply(parallel, 
                                                  seqpermutation, ptest.ts, comm = comm, phylodist = phylodist, 
                                                  method = method, squareroot = squareroot, mt = MT)
    }
    res$curve.null.ts <- res.curve.null.ts
  }
  if (null.model.bm) {
    seqpermutation <- vector("list", runs)
    if (!inherits(parallel, "cluster")) {
      res.curve.null.bm <- vector("list", runs)
      for (i in 1:runs) {
        res.curve.null.bm[[i]] <- ptest.bm(NULL, tree, 
                                           comm, res.values, res.vectors, ranks = ranks)
        if (progressbar) {
          SYNCSA::ProgressBAR(i + runs, BarRuns, style = 3)
        }
      }
    }
    else {
      res.curve.null.bm <- parallel::clusterApply(parallel, 
                                                  seqpermutation, ptest.bm, tree = tree, comm = comm, 
                                                  values = res.values, vectors = res.vectors, ranks = ranks)
    }
    res$curve.null.bm <- res.curve.null.bm
  }
  if (newClusters) {
    parallel::stopCluster(parallel)
  }
  class(res) <- "pcpscurve"
  return(res)
}
