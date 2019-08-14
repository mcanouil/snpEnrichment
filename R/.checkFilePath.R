.checkFilePath <- function(path) {
  gsub("/*$", "/", path)
}
