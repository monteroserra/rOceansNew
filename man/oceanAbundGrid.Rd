% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions_rOceans_2_global patterns
%   biodiversity v.0.0.6.r
\name{oceanAbundGrid}
\alias{oceanAbundGrid}
\title{Function 2.2 oceanAbundGrid: compute spatial patterns of abundance of occurrences at 
multiple scales
compute spatial patterns of abundance from larga occurrences datasets at multiple scales}
\usage{
oceanAbundGrid(occurrences, lat_name = "decimalLatitude",
  long_name = "decimalLongitude", extent = "global", cell_size = 5,
  min_long = -180, max_long = 180, min_lat = -90, max_lat = 90)
}
\arguments{
\item{occurrences}{dataset  of occurrences downloaded from OBIS or GBIF}

\item{extent}{spatial extent of the analysis:"global", "Mediterranean"}

\item{cell_size}{size of the grid cells in degrees (?)}

\item{min_long}{minimum longitude of the analysis}

\item{max_long}{maximum longitude of the analysis}

\item{max_lat}{maximum latitude of the analysis}

\item{from}{source data, it can be "OBIS" or"GBIF"}

\item{latitude_name}{variable name indicating latitude in the occurrences dataset}

\item{longitude_name}{variable name indicating longitude in the occurrences dataset}

\item{max_long}{minimum latitute of the analysis}
}
\value{
The function returns a an object of class raster with abundance per cell
}
\description{
computes global grid of abundance of occurrences at multiple scales
}
\details{
The function creates downloads occurrences data from the specific sources
and stores them in a data frame.
}
\author{
I. Montero-Serra,  E. Aspillaga, V. Barve, N Barve & K. Kaplan,
}
