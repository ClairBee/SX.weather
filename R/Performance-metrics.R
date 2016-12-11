
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

# MULTIVARIATE METRICS                                                                          ####

#' CRPS for multivariate-normal predictive density
#' 
#' CRPS (energy score) for fitted multivariate-normal predictive density.
#' @details EITHER parameter 'efc' (for an ensemble forecast) OR parameters 'mu', 'sig' and 'k' (for a fitted univariate or multivariate normal density) must be supplied
#' @param o Vector of observations
#' @param mu Vector of means of fitted model
#' @param sig Covariance matrix or variance of fitted model - NOT the standard deviation
#' @param efc Matrix of multivariate forecasts, each column an ensemble member.
#' @param k Number of samples to use in Monte Carlo approximation. Default is 1000.
#' @return Energy score (multivariate CRPS)
#' @export
#' 
es.crps <- function(o, mu, sig, efc, k = 1000) {
    
    if(!missing(efc)) {
        
        # ENSEMBLE FORECAST
        
        # if univariate, pad dimensions of ensemble
        if(length(o) == 1) {efc <- array(efc, dim = c(1, length(efc)))}
            
        m <- ncol(efc)
        
        norm.1 <- mean(apply(efc-o, 2, function(v) sqrt(sum(v^2))))
        norm.2 <- sum(colSums(apply(efc, 2, function(i) apply(efc-i, 2, function(v) sqrt(sum(v^2)))))) / (2 * m^2)
    
    } else {
        
        if (length(mu) == 1) {
            
            # UNIVARIATE NORMAL FITTED DENSITY
            
            x <- rnorm(k, mean = mu, sd = sqrt(sig))
            
            norm.1 <- mean(apply(t(x)-o, 2, function(v) sqrt(sum(v^2))))
            norm.2 <- sum(sapply(x[1:k-1] - x[2:k], function(v) sqrt(sum(v^2)))) / (2 * (k-1))
            
        } else {
            
            # MULTIVARIATE NORMAL PREDICTIVE DENSITY
            require(mvtnorm)            # needed to simulate from mimic
            
            x <- rmvnorm(k, mean = mu, sigma = sig)
            
            norm.1 <- mean(apply(t(x)-o, 2, function(v) sqrt(sum(v^2))))
            norm.2 <- sum(apply(x[1:k-1,] - x[2:k,], 1, function(v) sqrt(sum(v^2)))) / (2 * (k-1))
        }
    }
    norm.1 - norm.2
}



#' Verification ranks
#' 
#' Rank of verifying observation against an ensemble forecast
#' @details Treats multivariate ensemble forecast as is univariate, by collapsing all dimensions except ensemble size.
#' @param o Vector or matrix of observations
#' @param ens Matrix of forecasts from ensemble. Last dimension gives ensemble size, all preceding dimensions must match observation dimensions.
#' @return Ranks of observation within ensemble spread.
#' @export
#' 
verif.ranks <- function(o, ens) {
    
    # univariate case: don't distinguish between variables
    # collapse all dimensions into single vector (except for ensemble members)
    
    l <- length(dim(ens))
    
    apply(abind(o, ens, along = l), -l, function(v) which(sort(v) == v[1]))
}



#' Verification rank histogram
#'
#' Plots verification rank histogram by calling \code{\link{verif.ranks}}. Treats multivariates forecast as if univariate.
#' @param o Vector or matrix of observations
#' @param ens Matrix of forecasts from ensemble. Last dimension gives ensemble size, all preceding dimensions must match observation dimensions.
#' @export
#' 
vr.hist <- function(o, ens, breaks = c(0:10)/10, main = "", col = "skyblue", ...) {
    
    m <- dim(ens)[length(dim(ens))] + 1
    vr <- verif.ranks(o, ens)
    hist(vr/m, breaks = breaks, col = col, main = main, xlab = "", ylab = "", prob = T, ...)
    abline(h = 1, col = "red3", lty = 2)
}

#' Box density ordinate transform
#' 
#' Run Box density ordinate transform over multivariate normal model to assess calibration.
#' @details Returns a centre-outward ordering: outliers tend to have low values, inliers tend to have high values. (Tells us nothing about direction of any bias, though)
#' @param o Vector of observations
#' @param mu Vector of fitted means
#' @param s Fitted covariance matrix
#' @return Ordinate-transformed value
#' @export
#' 
box.dot <- function(o, mu, s) {
    
    1 - pchisq(t(o - mu) %*% solve(s) %*% (o-mu), length(mu))
    
}


#' Determinant sharpness
#' 
#' Calculate determinant sharpness of a multivariate predictive density.
#' @param sig Covariance matrix of fitted model
#' @return Calculated determinant sharpness
#' @export
#' 
det.sharpness <- function(sig) {
    det(sig)^(1 / (2 * dim(sig)[1]))
}


#' Deviation from uniformity
#' 
#' Calculate deviation from uniformity of a vector of verification ranks
#' @param vr Vector of verification ranks, as returned by \code{\link{verif.ranks}}.
#' @param l Maximum possible rank. Default is \code{max(vr)}, but if observation is never higher than all ensemble members, this should be manually set to ensemble size + 1.
#' @return Score evaluating degree of deviation from uniformity.
#' @export
#' 
u.dev <- function(vr, l = max(vr)) {
    sum(sapply(1:l, function(j) abs(mean(vr == j) - 1/l)))
}