% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare_subclones_for_timing.R
\name{CNA_annotation}
\alias{CNA_annotation}
\title{Collate the subclones data from Battenberg into a single file and add the ploidy of the sample.}
\usage{
CNA_annotation(allsegs_file, tumour_type, output_dir)
}
\arguments{
\item{allsegs_file}{Full path to the directory with the allsegs file with the subclones information}

\item{tumour_type}{String that represents the type of tumour we work with}

\item{output_dir}{Full path to the directory with the annotated allsegments.txt output from this function}
}
\value{
\emph{tumour_type}_allsegments.txt file with information about the samplenames, copy number aberrations (startpos, endpos, copy number state of major/minor allele), ploidy and type of CNA
}
\description{
Collate the subclones data from Battenberg into a single file and add the ploidy of the sample.
}
