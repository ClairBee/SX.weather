
#' Calculate fitted RMSE
#' 
#' Calculate RMSE for a fitted model
#' @param Tau Vectors of means for all dimensions
#' @param actual Array of observations, against which fitted means are to be assessed. Default is to use built-in \link{obs} data.
#' @export
model.rmse <- function(Tau, actual = obs) {
    
    actual <- actual[dimnames(Tau)[[1]],,]
    diff.mat <- abind(actual, Tau, along = 4)
    
    # subtract obs for each leadtime
    diff <- apply(diff.mat, 2:3, 
                  function(arr) sweep(arr[,-1], 1, arr[,1], "-"))
    
    res <- array(sqrt(apply(diff^2, 1, mean)),
                 dim = dim(Tau[,1,1,]),
                 dimnames = dimnames(Tau[,1,1,]))
    return(res)
}


#' Calculate fitted CRPS
#' 
#' Calculate CRPS for fitted model
#' @param Tau Vectors of means for all dimensions
#' @param S Variance/covariance matrices for all dimensions
#' @return Array containing CRPS scores for each variable, day, year and leadtime.
#' @export
#' 
model.crps <- function(Tau, S) {
    
    actual <- obs[dimnames(Tau)[[1]],,]
    
    array(verification::crps(c(rep(actual, dim(Tau)[[4]])),
               cbind(c(Tau),
                     sqrt(apply(S, 3:5, diag))))$crps, 
          dim = dim(Tau), dimnames = dimnames(Tau))
}


#' Calculate spread, RMSE and CRPS for a fitted model
#' 
#' @param Tau Vectors of means for all dimensions
#' @param S Variance/covariance matrices for all dimensions
#' @param actual Array of observations, against which fitted means are to be assessed. Default is to use built-in \link{obs} data.
#' @return List containing arrays of calculated spread, RMSE and CRPS
#' @export
#' 
model.performance <- function(fitted, actual = obs) {
    
    Tau <- fitted[[1]]
    S <- fitted[[2]]
    
    list(spread = sqrt(apply(apply(S, 3:5, diag), c(1,4), mean)),
         rmse = model.rmse(Tau, actual),
         crps = model.crps(Tau, S))
}


#' Ensemble CRPS
#' 
#' Calculate CRPS for ensemble prediction, for each variable & each leadtime.
#' @param model Array of predictions for which CRPS is to be calculated: variables · days · years · leadtimes · ensemble members.
#' @return Array of CRPS scores for each variable and leadtime.
#' @export
#' 
ensemble.crps <- function(model) {
    
    # bind observations & forecasts to allow easier array processing
    fudge <- abind("obs" = array(rep(obs, dim(model)[[4]]), dim = c(dim(obs), dim(model)[[4]])),
                   offset.forecast(model), along = 5)
    
    # collapse days & years into single dimension - note change in dimension order!
    fudge <- apply(fudge, c(1,4:5), rbind)
    
    # apply crpsDecomposition across all variables & leadtimes
    apply(fudge, c(2:3), function(arr) crpsDecomposition(obs = arr[,1], eps = arr[,-1])$CRPS)
}


####################################################################################################

# FUNCTIONS RUNNING OVER FORECAST (NO MODEL CONSTRAINTS)                                        ####

#' Forecast MAE
#' 
#' Median Absolute Error calculated over forecasts and observations
#' @param fc Vector or array of forecasts to evaluate
#' @param o Vector or array of observations. Dimensions must match those of fc
#' @return Calculated MAE
#' @export
#' 
fc.mae <- function(fc, o) {
    mean(abs(fc - o))
}


#' Forecast RMSE
#' 
#' Root Mean Squared Error, calculated over forecasts and observations
#' @param fc Vector or array of forecasts to evaluate
#' @param o Vector or array of observations. Dimensions must match those of fc
#' @return Calculated RMSE
#' @export
#' 
fc.rmse <- function(fc, o) {
    sqrt(mean((fc - o)^2))
}


####################################################################################################
