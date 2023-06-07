#! /usr/bin/env Rscript

# Load required libraries
library(rmarkdown)
library(optparse)
library(knitr)
library(dplyr)
arg_list <- list(
  make_option(
    c('--csv_files'), action = 'store', type = 'character',
    help = 'Path to the directory containing CSV files'))

opt <- parse_args(OptionParser(option_list=arg_list))

# Get the list of CSV files in the directory
csv_files <- list.files(opt$csv_files, pattern = "\\.csv$", full.names = TRUE)

# Define a list to store the loaded data frames
data_frames <- list()

# Load the CSV files into separate data frames
for (file in csv_files) {
  df <- read.csv(file)
  data_frames[[length(data_frames) + 1]] <- df
}

# Merge the data frames based on the "scRNA_Metrics" column
merged_df <- Reduce(function(x, y) merge(x, y, by = "scRNA_Metrics"), data_frames)

# Remove the "X" prefix from column names
colnames(merged_df) <- gsub("^X", "", colnames(merged_df))

# Save the final result
write.csv(merged_df, "experiment_result.csv", row.names = FALSE)

rmd_content <- paste(
  "---",
  "title: Single-cell RNA-seq Report",
  "output: html_document",
  "---",
  "",
  "# QC Report per Experiment",
  "",
  "[Download CSV file](", "experiment_result.csv", ")",
  "",
  "```{r, echo=FALSE}",
  "knitr::kable(merged_df)",
  "```",
  sep = "\n"
)

# Save R Markdown template file
writeLines(rmd_content, "experiment_report.Rmd")

# Generate HTML report
rmarkdown::render("experiment_report.Rmd", output_file = "experiment_report.html")

