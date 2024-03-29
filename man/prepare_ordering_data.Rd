% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare_enriched_regions_for_ordering.R
\name{prepare_ordering_data}
\alias{prepare_ordering_data}
\title{Process the enriched data into the format required for the ordering step with the Plackett-Luce model}
\usage{
prepare_ordering_data(
  annotated_segments_file,
  tumour_type,
  enriched_dir,
  hg19_chr_coordinates,
  output_dir
)
}
\arguments{
\item{annotated_segments_file}{Full path to the file with the annotated CNA}

\item{tumour_type}{String that represents the type of tumour we work with}

\item{enriched_dir}{Full path the directory with the enriched regions}

\item{hg19_chr_coordinates}{Full path to the file with the chromosome coordinates of the hg19 build genome}

\item{output_dir}{Full path to the directory where the files with the enriched regions will be output}
}
\value{
files with enriched regions (separate files for LOH, HD and gain) in the format for the ordering script
}
\description{
Process the enriched data into the format required for the ordering step with the Plackett-Luce model
}
