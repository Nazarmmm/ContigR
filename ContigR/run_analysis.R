#!/usr/bin/env Rscript

# Load required libraries
# library(ContigR)

Rcpp::sourceCpp("./src/contig_analysis.cpp")

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript run_analysis.R <input_fasta_file> [maxk] [mink] [output_dir]")
}

# Set default values
input_file <- args[1]
maxk <- if (length(args) > 1) as.integer(args[2]) else 100
mink <- if (length(args) > 2) as.integer(args[3]) else 50
output_dir <- if (length(args) > 3) args[4] else "output"

# Run the analysis
cat("Starting contig analysis...\n")
cat("Input file:", input_file, "\n")
cat("Parameters:\n")
cat("  maxk:", maxk, "\n")
cat("  mink:", mink, "\n")
cat("  output directory:", output_dir, "\n")

# Run the analysis
result <- analyze_contigs_cpp(
  filepath = input_file,
  maxk = maxk,
  mink = mink,
  output_dir = output_dir
)

# Print summary
cat("\nAnalysis completed successfully!\n")
cat("Results are saved in:", output_dir, "\n")
cat("Adjacency table:", result$adjacency_table_path, "\n")
cat("Execution times are saved in:", result$execution_times_path, "\n") 