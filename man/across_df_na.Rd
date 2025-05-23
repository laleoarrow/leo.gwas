% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/0.basic_tools.R
\name{across_df_na}
\alias{across_df_na}
\title{Across a df to count na}
\usage{
across_df_na(df)
}
\arguments{
\item{df}{a data frame}
}
\value{
a data frame with the number of NA values in each column
}
\description{
This function summarize NA values in each column of a data frame.
}
\examples{
df <- data.frame(a = c(1, 2, NA, 4), b = c(NA, 2, 3, 4))
summarize_na(df)
}
