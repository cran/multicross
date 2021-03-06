% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multicross.R
\name{split_mat}
\alias{split_mat}
\title{Split a data frame or matrix into subsets based on a particular categorical variable}
\usage{
split_mat(obj, by)
}
\arguments{
\item{obj}{is a data frame or matrix to be split into subsets, divided by the categorical variable}

\item{by}{is a character-string that specifies the columns that need to be subsetted}
}
\value{
A list containing the subsetted data sets. The names of the list corresponds to the value of the subsetted list
}
\description{
Split a data frame or matrix into subsets based on a particular categorical variable
}
