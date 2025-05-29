library(ggplot2)
library(tidyr)
library(dplyr)

# Read the data
data <- read.delim("output__adjacent_contigs.tsv", header=TRUE, comment.char="#")
# Convert column names to match the data
colnames(data) <- c("Index", "Contig_name", "Length", "Coverage", "GC", "Multiplicity", "Annotation", "Start", "End")

# Basic statistics
cat("Summary Statistics:\n")
print(summary(data))

# 1. Contig Length Distribution
p1 <- ggplot(data, aes(x=Length)) +
  geom_histogram(bins=30, fill="steelblue", color="black") +
  theme_minimal() +
  labs(title="Distribution of Contig Lengths",
       x="Length (bp)",
       y="Count") +
  scale_x_log10()
ggsave("contig_lengths_distribution.png", p1, width=10, height=6, dpi=300)

# 2. GC Content Distribution
p2 <- ggplot(data, aes(x=GC)) +
  geom_histogram(bins=30, fill="darkgreen", color="black") +
  theme_minimal() +
  labs(title="Distribution of GC Content",
       x="GC Content (%)",
       y="Count")
ggsave("gc_content_distribution.png", p2, width=10, height=6, dpi=300)

# 3. Coverage Distribution (excluding missing values)
p3 <- ggplot(data[data$Coverage != "0.000000",], aes(x=as.numeric(Coverage))) +
  geom_histogram(bins=30, fill="coral", color="black") +
  theme_minimal() +
  labs(title="Distribution of Coverage",
       x="Coverage",
       y="Count")
ggsave("coverage_distribution.png", p3, width=10, height=6, dpi=300)

# 4. Multiplicity Distribution
p4 <- ggplot(data, aes(x=factor(Multiplicity))) +
  geom_bar(fill="purple", color="black") +
  theme_minimal() +
  labs(title="Distribution of Multiplicity",
       x="Multiplicity",
       y="Count")
ggsave("multiplicity_distribution.png", p4, width=10, height=6, dpi=300)

# 5. Scatter plot: Length vs GC Content
p5 <- ggplot(data, aes(x=Length, y=GC)) +
  geom_point(alpha=0.6, color="blue") +
  theme_minimal() +
  labs(title="Length vs GC Content",
       x="Length (bp)",
       y="GC Content (%)") +
  scale_x_log10()
ggsave("length_vs_gc_content.png", p5, width=10, height=6, dpi=300)

# Additional analysis: Calculate correlations
cat("\nCorrelations:\n")
numeric_data <- data %>%
  mutate(Coverage = as.numeric(Coverage)) %>%
  select(Length, GC, Coverage, Multiplicity)
print(cor(numeric_data, use="pairwise.complete.obs"))

# Save summary statistics to a file
sink("contig_statistics.txt")
cat("=== Contig Analysis Summary ===\n\n")
cat("Total number of contigs:", nrow(data), "\n")
cat("\nLength Statistics:\n")
print(summary(data$Length))
cat("\nGC Content Statistics:\n")
print(summary(data$GC))
cat("\nMultiplicity Statistics:\n")
print(table(data$Multiplicity))
sink() 