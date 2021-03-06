% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bic_plots.R, R/tmp_bic_plots.R
\name{bic.plot.read.distribution}
\alias{bic.plot.read.distribution}
\alias{bic.plot.read.distribution}
\title{Plot RSeQC read distribution stats}
\usage{
bic.plot.read.distribution(dat, file = NULL, stack = TRUE, pct = FALSE,
  col.pal = "Set3")

bic.plot.read.distribution(dat, file = NULL, stack = TRUE, pct = FALSE,
  col.pal = "Set3")
}
\arguments{
\item{dat}{data frame containing merged output from multiple runs of
RSeQC's read_distribution.py script, where rows are
samples and columns are metrics; must contain "Samples" slot
and may contain any other slots}

\item{file}{PDF file to which plot should be saved (optional)}

\item{stack}{logical indicating that bar chart should be stacked; Default: TRUE}

\item{pct}{logical indicating that plot should show percentages}

\item{col.pal}{name of color palette; must be from list of RColorBrewer palettes
Default: "Set3"}

\item{horizontal}{logical indicating that bars should be horizontal; Default: TRUE}

\item{dat}{data frame containing merged output from multiple runs of
RSeQC's read_distribution.py script, where rows are
samples and columns are metrics; must contain "Samples" slot
and may contain any other slots}

\item{pct}{logical indicating that plot should show percentages}

\item{stack}{logical indicating that bar chart should be stacked; Default: TRUE}

\item{col.pal}{name of color palette; must be from list of RColorBrewer palettes
Default: "Set3"}

\item{file}{PDF file to which plot should be saved (optional)}
}
\description{
Plot distribution of reads across different genomic features
like exons, introns, TSS, etc. for all samples in a project

Plot distribution of reads across different genomic features
like exons, introns, TSS, etc. for all samples in a project
}
