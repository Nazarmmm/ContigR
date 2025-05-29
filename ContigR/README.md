# ContigR

ContigR is an R package for analyzing and visualizing contigs from FASTA files. It provides functionality for detecting overlaps between contigs, assigning multiplicity, and generating various visualizations and reports.

## Installation

To install the package, you'll need to have R and Rcpp installed. Then you can install ContigR using:

```R
install.packages("devtools")
devtools::install_github("yourusername/ContigR")
```

## Dependencies

- R (>= 3.6.0)
- Rcpp (>= 1.0.7)
- ggplot2
- tidyr
- dplyr

## Usage

The package can be used either through R functions or via a command-line script.

### Using R Functions

```R
library(ContigR)

# Analyze contigs from a FASTA file
result <- analyze_contigs(
  filepath = "path/to/your/contigs.fasta",
  maxk = 100,  # maximum k-mer size
  mink = 50,   # minimum k-mer size
  output_dir = "output",
  num_iterations = 1
)
```

### Using Command Line Script

```bash
Rscript run_analysis.R path/to/your/contigs.fasta [maxk] [mink] [output_dir] [num_iterations]
```

Parameters:
- `input_fasta_file`: Path to the input FASTA file (required)
- `maxk`: Maximum k-mer size (default: 100)
- `mink`: Minimum k-mer size (default: 50)
- `output_dir`: Output directory (default: "output")
- `num_iterations`: Number of iterations to run (default: 1)

## Output Files

The analysis generates several output files in the specified output directory:

1. `__adjacent_contigs.tsv`: A tab-separated file containing information about adjacent contigs
2. `_full_matching_log.txt`: Detailed log of contig matches
3. `_summary.txt`: Summary statistics of the analysis
4. `_annotated_genbank.gtf`: GenBank format file with contig annotations
5. Various visualization plots:
   - Length distribution
   - GC content distribution
   - Coverage distribution
   - Multiplicity distribution
   - Length vs GC content scatter plot

## Functions

### Main Functions

- `analyze_contigs()`: Main function to run the contig analysis
- `create_visualizations()`: Generate visualizations from the analysis results

### Visualization Functions

- `create_length_distribution()`: Create histogram of contig lengths
- `create_gc_distribution()`: Create histogram of GC content
- `create_coverage_distribution()`: Create histogram of coverage
- `create_multiplicity_distribution()`: Create histogram of multiplicity
- `create_length_vs_gc_plot()`: Create scatter plot of length vs GC content

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. 