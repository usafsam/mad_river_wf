#!/usr/bin/env Rscript

library(argparser, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(readr, quietly = TRUE)
library(readxl, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(tools, quietly = TRUE)
library(writexl, quietly = TRUE)

my_parser <-
  arg_parser(
    "Return plot of a coverage for a gene, given a nextstrain tabular file and the gene's boundaries."
  ) %>%
  add_argument("in_file",
               "Tabular file from nextalign. Must have a column labeled 'missing'.",
               type = "character") %>%
  add_argument("gene_start", "Start of gene. (e.g. 21563)", type = "numeric") %>%
  add_argument("gene_end", "End of gene. (e.g. 25384)", type = "numeric") %>%
  add_argument("--sample_id", "Column to use for the sample ID.", default = "SampleID") %>%
  add_argument(
    "--out_type",
    "Type of output file to write tabular data to. Defaults to extension of input file.",
    type = "character"
  )

my_args <- parse_args(my_parser)

in_file_extension <- file_ext(my_args[["in_file"]])
csv_delimiter <- ","

if (in_file_extension == "xls" || in_file_extension == "xlsx") {
  my_data <-
    read_excel(my_args[["in_file"]],
               col_types = cols(missing = col_character()),
               trim_ws = TRUE)
} else if (in_file_extension == "tsv") {
  my_data <- read_delim(
    my_args[["in_file"]],
    "\t",
    escape_double = FALSE,
    col_types = cols(missing = col_character()),
    trim_ws = TRUE
  )
} else if (in_file_extension == "csv") {
  my_data <- read_csv(my_args[["in_file"]],
                      col_types = cols(missing = col_character()),
                      trim_ws = TRUE)
  if (ncol(my_data) == 1) {
    my_data <- read_delim(
      my_args[["in_file"]],
      ";",
      escape_double = FALSE,
      col_types = cols(missing = col_character()),
      trim_ws = TRUE
    )
    csv_delimiter <- ";"
  }
} else {
  stop(paste("Unrecognized extension:", in_file_extension), call. = FALSE)
}

if (!my_args[["sample_id"]] %in% colnames(my_data)) {
  print(my_parser)
  stop(paste("Column `", my_args[["sample_id"]], "` doesn't exist in `", my_args[["in_file"]], "`!", sep =
               ""),
       call. = FALSE)
}

gene_range <- seq(my_args[["gene_start"]], my_args[["gene_end"]])

get_coverage_gaps <- function(gap_string) {
  coverage_gaps <- c()
  gap_vector <- gap_string %>% strsplit(",") %>% unlist()
  for (gap in gap_vector) {
    gap_start_end <- gap %>% strsplit("-") %>% unlist()
    gap_start <- as.numeric(gap_start_end[1])
    if (length(gap_start_end) == 2) {
      gap_end <- as.numeric(gap_start_end[2]) - 1
    } else {
      gap_end <- as.numeric(gap_start_end[1])
    }
    if (is.na(gap_start)) {
      return(NULL)
    }
    coverage_gaps <- c(coverage_gaps, seq(gap_start, gap_end))
  }
  return(sort(unique(coverage_gaps)))
}

get_gene_coverage <- function(gap_string, gene_range) {
  coverage_gaps <- get_coverage_gaps(gap_string)
  if (is.null(coverage_gaps)) {
    covered_bases <- c()
  } else {
    covered_bases <- setdiff(gene_range, coverage_gaps)
  }
  return(length(covered_bases) / length(gene_range))
}

get_missing_bases <- function(gap_string, gene_range) {
  coverage_gaps <- get_coverage_gaps(gap_string)
  if (is.null(coverage_gaps)) {
    return(as.numeric(length(gene_range)))
  } else {
    return(as.numeric(length(
      intersect(coverage_gaps, gene_range)
    )))
  }
}

my_data <-
  my_data %>%
  mutate(missing_no = replace_na(missing, paste(my_args[["gene_start"]], my_args[["gene_end"]] +
                                                  1, sep = "-"))) %>%
  mutate(gene_coverage = as.numeric(lapply(
    missing_no, get_gene_coverage, gene_range
  ))) %>%
  mutate(missing_bases = as.numeric(lapply(
    missing_no, get_missing_bases, gene_range
  )))

out_data <- my_data %>% select(my_args[["sample_id"]])

out_data[[paste("CoverageProportion(", my_args[["gene_start"]], "-", my_args[["gene_end"]], ")", sep =
                  "")]] <- my_data$gene_coverage
out_data[[paste("NumMissingBases(", my_args[["gene_start"]], "-", my_args[["gene_end"]], ")", sep =
                  "")]] <- my_data$missing_bases

out_file_extension <- in_file_extension

if (!is.na(my_args[["out_type"]])) {
  out_file_extension <- my_args[["out_type"]]
}

if (out_file_extension == "xls" || out_file_extension == "xlsx") {
  write_xlsx(out_data,
             path = paste(
               file_path_sans_ext(my_args[["in_file"]]),
               "_",
               my_args[["gene_start"]],
               "-",
               my_args[["gene_end"]],
               ".xlsx",
               sep = ""
             ))
} else if (out_file_extension == "tsv") {
  write_delim(
    out_data,
    paste(
      file_path_sans_ext(my_args[["in_file"]]),
      "_",
      my_args[["gene_start"]],
      "-",
      my_args[["gene_end"]],
      ".tsv",
      sep = ""
    ),
    "\t"
  )
} else if (out_file_extension == "csv") {
  if (csv_delimiter == ",") {
    write_csv(
      out_data,
      paste(
        file_path_sans_ext(my_args[["in_file"]]),
        "_",
        my_args[["gene_start"]],
        "-",
        my_args[["gene_end"]],
        ".csv",
        sep = ""
      )
    )
  } else {
    write_delim(
      out_data,
      paste(
        file_path_sans_ext(my_args[["in_file"]]),
        "_",
        my_args[["gene_start"]],
        "-",
        my_args[["gene_end"]],
        ".csv",
        sep = ""
      ),
      ";"
    )
  }
} else {
  stop(paste("Unrecognized extension:", out_file_extension),
       call. = FALSE)
}
