# ============================================================
# Record the R and package versions used for the analysis
#
# Run this file after the main analysis in the same R session:
#   source("record_session_info.R")
# ============================================================

output_file <- "session-info.txt"

capture.output(
  {
    cat("Analysis environment\n")
    cat("====================\n\n")
    cat("Recorded on: ", format(Sys.time(), tz = Sys.timezone()), "\n", sep = "")
    cat("Working directory: ", normalizePath(getwd()), "\n\n", sep = "")
    print(sessionInfo())
  },
  file = output_file
)

message("Session information saved to ", normalizePath(output_file))
