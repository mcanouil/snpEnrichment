mclapply2 <- function(X, FUN, ..., mc.preschedule = TRUE, mc.set.seed = TRUE, mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L), mc.cleanup = TRUE, mc.allow.recursive = TRUE) {
  if (Sys.info()[["sysname"]] != "Linux") {
    mc.cores <- 1
  } else {
    mc.cores <- min(parallel::detectCores(), mc.cores)
  }
  return(parallel::mclapply(
    X = X, FUN = FUN, ...,
    mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed, mc.silent = mc.silent,
    mc.cores = maxCores(mc.cores), mc.cleanup = mc.cleanup, mc.allow.recursive = mc.allow.recursive
  ))
}
