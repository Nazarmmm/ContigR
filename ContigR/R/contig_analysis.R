#' Analyze contigs from a FASTA file
#' 
#' @param filepath Path to the input FASTA file
#' @param maxk Maximum k-mer size for analysis
#' @param mink Minimum k-mer size for analysis
#' @param output_dir Directory to save output files
#' @param num_iterations Number of iterations for performance analysis
#' @return A list containing analysis results and execution times
#' @export
analyze_contigs <- function(filepath, maxk = 50, mink = 5, output_dir = "Output", num_iterations = 100) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Run the analysis
  results <- analyze_contigs_cpp(filepath, maxk, mink, output_dir, num_iterations)
  
  # Create visualizations
  create_visualizations(results$adjacency_table_path, output_dir)
  
  return(results)
}

#' Create visualizations from contig analysis results
#' 
#' @param data_path Path to the adjacency table file
#' @param output_dir Directory to save visualizations
#' @return NULL
#' @export
create_visualizations <- function(data_path, output_dir) {
  # Read the data
  data <- read.delim(data_path, header = TRUE, comment.char = "#")
  
  # Convert column names to match the data
  colnames(data) <- c("Index", "Contig_name", "Length", "Coverage", "GC", 
                     "Multiplicity", "Annotation", "Start", "End")
  
  # Create visualizations
  create_length_distribution(data, output_dir)
  create_gc_distribution(data, output_dir)
  create_coverage_distribution(data, output_dir)
  create_multiplicity_distribution(data, output_dir)
  create_length_vs_gc_plot(data, output_dir)
  
  # Save summary statistics
  save_summary_statistics(data, output_dir)
}

#' Create length distribution plot
#' @keywords internal
create_length_distribution <- function(data, output_dir) {
  p <- ggplot(data, aes(x = Length)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black") +
    theme_minimal() +
    labs(title = "Distribution of Contig Lengths",
         x = "Length (bp)",
         y = "Count") +
    scale_x_log10()
  
  ggsave(file.path(output_dir, "contig_lengths_distribution.png"), 
         p, width = 10, height = 6, dpi = 300)
}

#' Create GC content distribution plot
#' @keywords internal
create_gc_distribution <- function(data, output_dir) {
  p <- ggplot(data, aes(x = GC)) +
    geom_histogram(bins = 30, fill = "darkgreen", color = "black") +
    theme_minimal() +
    labs(title = "Distribution of GC Content",
         x = "GC Content (%)",
         y = "Count")
  
  ggsave(file.path(output_dir, "gc_content_distribution.png"), 
         p, width = 10, height = 6, dpi = 300)
}

#' Create coverage distribution plot
#' @keywords internal
create_coverage_distribution <- function(data, output_dir) {
  p <- ggplot(data[data$Coverage != "0.000000",], 
             aes(x = as.numeric(Coverage))) +
    geom_histogram(bins = 30, fill = "coral", color = "black") +
    theme_minimal() +
    labs(title = "Distribution of Coverage",
         x = "Coverage",
         y = "Count")
  
  ggsave(file.path(output_dir, "coverage_distribution.png"), 
         p, width = 10, height = 6, dpi = 300)
}

#' Create multiplicity distribution plot
#' @keywords internal
create_multiplicity_distribution <- function(data, output_dir) {
  p <- ggplot(data, aes(x = factor(Multiplicity))) +
    geom_bar(fill = "purple", color = "black") +
    theme_minimal() +
    labs(title = "Distribution of Multiplicity",
         x = "Multiplicity",
         y = "Count")
  
  ggsave(file.path(output_dir, "multiplicity_distribution.png"), 
         p, width = 10, height = 6, dpi = 300)
}

#' Create length vs GC content scatter plot
#' @keywords internal
create_length_vs_gc_plot <- function(data, output_dir) {
  p <- ggplot(data, aes(x = Length, y = GC)) +
    geom_point(alpha = 0.6, color = "blue") +
    theme_minimal() +
    labs(title = "Length vs GC Content",
         x = "Length (bp)",
         y = "GC Content (%)") +
    scale_x_log10()
  
  ggsave(file.path(output_dir, "length_vs_gc_content.png"), 
         p, width = 10, height = 6, dpi = 300)
}

#' Save summary statistics to a file
#' @keywords internal
save_summary_statistics <- function(data, output_dir) {
  sink(file.path(output_dir, "contig_statistics.txt"))
  cat("=== Contig Analysis Summary ===\n\n")
  cat("Total number of contigs:", nrow(data), "\n")
  cat("\nLength Statistics:\n")
  print(summary(data$Length))
  cat("\nGC Content Statistics:\n")
  print(summary(data$GC))
  cat("\nMultiplicity Statistics:\n")
  print(table(data$Multiplicity))
  sink()
} 