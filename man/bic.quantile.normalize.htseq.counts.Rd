% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bic_analyze_counts.R
\name{bic.quantile.normalize.htseq.counts}
\alias{bic.quantile.normalize.htseq.counts}
\title{Use quantile normalization to normalize raw HTSeq counts}
\usage{
bic.quantile.normalize.htseq.counts(htseq.counts, key = NULL)
}
\arguments{
\item{htseq.counts}{matrix containing raw counts from HTSeq. One column per 
sample, one row per gene}

\item{key}{matrix where column one contains sample IDs that matches 
those in counts file, and column two contains the condition 
groups to which they belong. Only those samples that are in 
the key matrix will be included in the normalization}
}
\value{
list containing two matrices (one of raw counts and one of quantile 
        normalized counts) and a vector of gene symbols if one was included in 
        the raw counts matrix
}
\description{
Run \code{normalizeBetweenArrays()} function to normalize
raw HTSeq counts and return a list containing both raw
and normalized counts
}
