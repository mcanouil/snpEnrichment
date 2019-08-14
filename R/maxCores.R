maxCores <- function(mc.cores = 1) {
  if (Sys.info()[["sysname"]] == "Linux") {
    nbCores <- parallel::detectCores()
    mc.cores.old <- mc.cores
    if (file.exists("/proc/meminfo")) {
      memInfo <- readLines("/proc/meminfo")
      sysMemFree <- memInfo[grep("^MemFree:", memInfo)]
      sysMemCached <- memInfo[grep("^Cached:", memInfo)]
      sysMemAvailable <- 0.95 * (as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", sysMemFree)) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", sysMemCached)))
      sysProc <- as.numeric(unlist(strsplit(system(paste("ps v", Sys.getpid()), intern = TRUE)[2], " +"), use.names = FALSE)[8])
      mc.cores <- max(min(as.numeric(mc.cores), floor(sysMemAvailable / sysProc)), 1)
      if (mc.cores > nbCores) mc.cores <- nbCores
      if (mc.cores != mc.cores.old) {
        warning(paste0('To avoid memory overload "mc.cores" was decreased to "', mc.cores, '".'), call. = FALSE)
      }
    } else {
      mc.cores <- ifelse(mc.cores.old > nbCores, nbCores, mc.cores.old)
    }
  } else {
    mc.cores <- 1
  }
  mc.cores
}
