% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare_enriched_regions_for_ordering.R
\name{remove_artefacts}
\alias{remove_artefacts}
\title{Remove artefacts from the enriched regions, eg. segments near telomere, centromere, HLA region, in a few samples}
\usage{
remove_artefacts(
  annotated_segments_file,
  tumour_type,
  output_dir,
  pvalues_dir,
  hg19_chr_coordinates,
  min_region = 10000,
  minN = 3
)
}
\arguments{
\item{annotated_segments_file}{Full path to the file with the annotated CNA}

\item{tumour_type}{String that represents the type of tumour we work with}

\item{output_dir}{Full path to the directory where the files with the enriched regions (and the p-values) will be output}

\item{pvalues_dir}{Full path to the directory with the FDR-corrected p-values and the enriched regions}

\item{hg19_chr_coordinates}{Full path to the file with the chromosome coordinates of the hg19 build genome}

\item{min_region}{A numeric representing the minimum length of enriched region; default = 1e4}

\item{minN}{A numeric representing the minimum number of events per segment to be used; default = 3}
}
\value{
files with enriched regions (separate files for LOH, HD and gain), with any artefacts removed
}
\description{
Remove artefacts from the enriched regions, eg. segments near telomere, centromere, HLA region, in a few samples
}
