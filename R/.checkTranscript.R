#' .checkTranscript
#'
#' @param distThresh
#' @keywords internal
.checkTranscript <- function(data, transcriptFile, distThresh) {
  if (any(transcriptFile != FALSE)) {
    transcript <- .readTranscript(transcriptFile = transcriptFile)
    transcriptCHR <- stats::na.exclude(transcript[transcript[, "CHR"] == unique(data[, "CHR"]), c("START", "END")])

    cisFunc <- function(line, distThresh, dataTranscript) {
      position <- as.numeric(line[4])
      CIS <- any(position > (dataTranscript[, "START"] - distThresh) & (position < (dataTranscript[, "END"] + distThresh)))
      c(line, CIS = as.numeric(CIS))
    }

    temp <- unique(as.data.frame(t(apply(data, 1, cisFunc, distThresh = distThresh * 10^3, dataTranscript = transcriptCHR)), stringsAsFactors = FALSE))
    temp[temp[, "CIS"] == 1, -grep("CIS", colnames(temp))]
  } else {
    unique(data)
  }
}
