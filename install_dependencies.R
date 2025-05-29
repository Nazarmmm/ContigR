# Install required packages
if (!require("devtools")) install.packages("devtools")
if (!require("Rcpp")) install.packages("Rcpp")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("tidyr")) install.packages("tidyr")
if (!require("dplyr")) install.packages("dplyr")

# Install ContigR package
devtools::install("ContigR") 