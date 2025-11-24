# .Rprofile - Clean startup without conda
# DO NOT source any files that contain conda commands

# Set environment variables
Sys.setenv(
  MKL_NUM_THREADS = "16",
  OMP_NUM_THREADS = "16",
  MKL_DYNAMIC = "TRUE"
)

# Display startup message
cat("=====================================\n")
cat("ESTARFM Project Environment\n")
cat("R Version:", R.version.string, "\n")
cat("Working Directory:", getwd(), "\n")
cat("Threads:", Sys.getenv("OMP_NUM_THREADS"), "\n")
cat("=====================================\n\n")

# Do NOT source these files:
# source("force_mkl_integration.R")  # COMMENTED OUT
# source("complete_intel_arc_solution.R")  # COMMENTED OUT