% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare_gnomAD.R
\name{prepareGnomAdClp}
\alias{prepareGnomAdClp}
\title{Prepare gnomAD variants table for xipher::prepareVcfClp.}
\usage{
prepareGnomAdClp(outGnomADPath, inGnomADVariantsTablePath, maf.min = 0.001)
}
\arguments{
\item{outGnomADPath}{path to write gnomAD variants table.  You can then pass this to prepareVcfClp() with gnomAdPath=}

\item{inGnomADVariantsTablePath}{path to gnomAD variants table produced by GATK VariantsToTable}

\item{maf.min}{minimum allele frequency}
}
\description{
Prepare gnomAD variants table for xipher::prepareVcfClp.
}
