.readTranscript <- function(transcriptFile) {
  if (all(class(try(close(file(transcriptFile)), silent = TRUE)) != "try-error")) {
    transcript <- utils::read.delim(file = transcriptFile, header = TRUE, stringsAsFactors = FALSE, na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE)
    transcript <- transcript[, 1:4]
    colnames(transcript) <- c("TRANSCRIPT", "CHR", "START", "END")
  } else {
    if (class(transcriptFile) %in% c("matrix", "data.frame")) {
      transcript <- transcriptFile[, 1:4]
      colnames(transcript) <- c("TRANSCRIPT", "CHR", "START", "END")
    } else {
      stop('[Enrichment:readEnrichment] "transcriptFile" required "matrix", "data.frame" or "txt file".', call. = FALSE)
    }
  }
  transcript
}
