
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
    
    daily.var <- array(dim = dim(ens[, 16:105, , , 1]), dimnames = dimnames(ens[, 16:105, , , 1]))
    
    invisible(sapply(0:14, function(lt) {
        daily.var[,,,toString(lt)] <<- apply(ens[,(16:105)-lt,,toString(lt),], 1:3, var)
    }))
    
    sqrt(apply(daily.var, c(1, 4), mean))
}


#' Adjust forecast by leadtime
#' 
#' Takes forecasts at multiple leadtimes and adjusts them so that forecasts for the same day are aligned (eg. look ahead 1 day for t-1, look ahead 2 days for t-2, and so on); rather than aligning records by date on which forecast was produced.
#' @param fc Forecast array
#' @return Adjusted forecast array
#' @export
#' 
offset.forecast <- function(fc) {
    os <- array(dim = dim(fc[, 16:105, , , ]), dimnames = dimnames(fc[, 16:105, , , ]))
    
    invisible(sapply(0:14, function(lt) {
        os[, , , toString(lt), ] <<- fc[, (16:105) - lt, , toString(lt), ]
    }))
    os
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
    
    fc <- pad.to.5d(fc)
    err <- pad.to.5d(array(dim = dim(fc), dimnames = dimnames(fc)))[,16:105,,,, drop = F]
    
    invisible(sapply(0:14, function(lt) {
        err[,,,toString(lt),] <<- sweep(fc[,(16:105) - lt,,toString(lt),], 1:3, actual, "-")
    }))
    
    err <- err[,,,,,drop = T]
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
    
    err <- pad.to.5d(forecast.errors(fc, actual = actual))
    sqrt(apply(err^2, c(1, 4, 5), mean))[,,,drop = T]
}


#' Create superensemble
#' 
#' Join ensemble forecasts into superensemble, ordered with control forecasts first, followed by perturbations.
#' @param data.list List of model names to be included in superensemble
#' @return Single array containing combined data, with attribute \code{m} a numeric vector labelling source ensemble
#' @export
#' 
superensemble <- function(data.list = list(ecmwf, ncep, ukmo)) {
    
    ctrl <- abind(lapply(data.list, "[",,,,,1), along = 5)
    pert <- abind(lapply(data.list, "[",,,,,-1), along = 5)
    
    super <- abind(ctrl, pert, along = 5)
    attr(super, "m") <- c(rep(0, length(data.list)), 
                         unlist(sapply(1:length(data.list), 
                                       function(n) rep(n, dim(data.list[[n]])[5]-1))))
    return(super)
}


#' Pad 4d array to 5d
#' 
#' Support function - pad an array forecast containing only one ensemble member to behave like an ensemble forecast.
#' @param arr Array to be padded
#' @return 5d array to be used in ensemble processing
#' @export
#' 
pad.to.5d <- function(arr) {
    if(length(dim(arr)) == 4) {
        arr <- array(arr, dim = c(dim(arr), 1), dimnames = append(dimnames(arr), "1"))
    }
    return(arr)
}





#' Square matrix
#' 
#' Support function - re-square an array after applying eg. cov
#' @param arr Array whose first dimension is to be rearranged into a square matrix
#' @return Rearranged array (if possible)
#' @export
square.mat <- function(arr) {
    
    d <- dim(arr)
    
    if(sqrt(d[1]) %% 1 != 0) {
        cat("Cannot produce square array.", "\n")
        return(arr)
    } else {
        d1 <- sqrt(d[1])
        return(array(arr, dim = c(d1, d1, d[-1])))
    }
}







