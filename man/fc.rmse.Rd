% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Performance-metrics.R
\name{fc.rmse}
\alias{fc.rmse}
\title{Forecast RMSE}
\usage{
fc.rmse(fc, o)
}
\arguments{
\item{fc}{Vector or array of forecasts to evaluate}

\item{o}{Vector or array of observations. Dimensions must match those of fc}
}
\value{
Calculated RMSE
}
\description{
Root Mean Squared Error, calculated over forecasts and observations
}

