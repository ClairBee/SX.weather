
#' Load .rda directly into variable
#' 
#' Wrapper function to load .rda object directly into a named variable, rather than loadng & reassigning
#' @param filenm Name of file (including path) to be loaded
#' @return Loaded object
#' @export
#' 
load.data <- function(filenm) {
    get(load(filenm))
}



#' Import forecasts & obs
#' 
#' Function to neatly load required observation & forecast data (keeps code tidy), truncating forecast to cover same date range as observations
#' @export
#' 
load.all <- function() {
    ecmwf <<- readRDS("./Data/ECMWF-forecasts.rds")
    ncep <<- readRDS("./Data/NCEP-forecasts.rds")
    ukmo <<- readRDS("./Data/UKMO-forecasts.rds")
    
    obs <<- readRDS("./Data/Observations.rds")
}


#' Calculate ensemble spread
#' 
#' For a given forecast ensemble, calculate the spread (standard deviation).
#' @param ens Array of ensemble data, variables · days · years · leadtimes · ensemble members
#' @return Array of spreads for each leadtime, variables · leadtimes
#' @export
#' 
ensemble.spread <- function(ens) {
    
    d <- dim(ens)
    # check against array calculation
    daily.var <- array(apply(array(ens, dim = c(d[1], d[2] * d[3], d[4], d[5])), 1:3, var),
                       dim = d[1:4], 
                       dimnames = list(dimnames(ens)[[1]], NULL, dimnames(ens)[[3]], NULL))
    
    sqrt(apply(daily.var, c(1, 4), mean))
}


#' Calculate forecast errors
#' 
#' Get forecast errors against observed values, accounting for forecast leadtime
#' @param fc Array of forecasts, variables · days · years · leadtimes · ensemble members
#' @param actual Array of observations, variables · days · years. Default is object \code{obs}.
#' @return Array of errors, with same dimensions as fc
#' @export
#' 
forecast.errors <- function(fc, actual = obs) {
    
    err <- array(dim = dim(fc[,16:105,,,]), dimnames = dimnames(fc[,16:105,,,]))
    
    invisible(sapply(0:14, function(lt) {
        err[,,,toString(lt),] <<- sweep(fc[,(16:105) - lt,,toString(lt),], 1:3, actual, "-")
    }))
    
    err
}


#' Forecast RMSE over full forecast period
#' 
#' Calculate forecast Root Mean Squared Error against observed values
#' @param fc Array of forecasts, variables · days · years · leadtimes · ensemble members
#' @param actual Array of observations, variables · days · years. Default is object \code{obs}.
#' @return Array of RMSE values, variables · leadtimes · ensemble members
#' @export
#' 
forecast.rmse <- function(fc, actual = obs) {
    err <- forecast.errors(fc, actual = actual)
    sqrt(apply(err^2, c(1, 4, 5), mean))
}


