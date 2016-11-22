
#' Apply 'solve' over an array
#' 
#' Applies 'solve' to find a %*% x = b for x, over a multidimensional array, and returns the result as an array of the original dimensions.
#' @param arr Array over which 'solve' is to be applied
#' @param solve.over Indices to retain in output
#' @export
array.solve <- function(arr, solve.over) {
    
    sol <- apply(arr, solve.over, solve)
    array(sol, dim = c(rep(sqrt(dim(sol)[1]),2), dim(sol)[-1]), 
          dimnames = dimnames(arr)[c(c(1:length(dim(arr)))[-solve.over], solve.over)])
}


#' Calculate covariance matrices
#' 
#' For a superensemble of offset-corrected forecasts, return an array of covariances between perturbed ensemble members for all variables across a single model.
#' @export
se.covariances <- function(fc) {
    array(apply(aperm(fc[,,,,-1], c(5,1,2,3,4)), 3:5, cov), dim = c(5,5,90,7,15))
}


#' Calculate sigma matrices
#' 
#' For a superensemble of offset-corrected forecasts, return an array of covariances between means of all perturbed ensemble members for all variables across all models.
#' @export
se.sigma <- function(data.list = list(ecmwf, ncep, ukmo)) {
    
    mean.fc <- abind(lapply(data.list,
                            function (model) {
                                apply(offset.forecast(model)[,,,,-1], 1:4, mean)
                            }), along = 0)
    
    sig <- apply(mean.fc, 3:5, cov)
    sig <- array(sig, dim = c(rep(sqrt(dim(sig)[[1]]),2), dim(sig[1,,,])),
                 dimnames = dimnames(mean.fc))
    dimnames(sig)[[1]] <- dimnames(sig)[[2]]
    
    return(sig)
}


#' Calculate lambda matrix
#' 
#' For a superensemble of offset-corrected forecasts, return an array of covariances between means of all perturbed ensemble members for all variables across all models.
#' @export
se.lambda <- function() {
    super.error <- forecast.errors(superensemble())
    ens.mean.error <- apply(super.error[,,,,-(1:3)], 1:4, mean)
    
    array(apply(aperm(ens.mean.error, c(3,1,2,4)), c(3:4), cov), dim = c(5,5,90,15),
          dimnames = append(dimnames(ens.mean.error[,,1,]), 
                            list(dimnames(ens.mean.error)[[1]]), after = 1))
} 


#' Calculate precision matrix
#' 
#' For a superensemble of offset-corrected forecasts, return array of precision matrices.
#' @export
se.precision <- function(D, L, vars = 1:dim(D)[[2]]) {
    
    # truncate lambda & D to include only required variables
    D <- D[,vars, vars,,,]
    L <- L[vars, vars,,]
    
    array.solve(sweep(array.solve(apply(array.solve(D,
                                                    c(1, 4:6)),          # solve per model/day/lt
                                        c(1:2, 4:6), sum),               # sum over models
                                  c(3:5)),                               # solve for sum
                      c(1:3, 5), L, "+"),                                # add lambda to solution
                3:5)   
}

#' Calculate tau
#' 
#' Calculate superensemble tau using Sichun's method 1
#' @export
se.tau <- function(D, L, Y, S, Eta, vars = dimnames(S)[[1]]) {
    
    # truncate data to match S-matrix
    D <- D[,vars, vars,,,]
    L <- L[vars, vars,,]
    Y <- Y[,vars,,,]
    Eta <- Eta[vars,,]

    # calculate
    d.inv.sum <-  apply(array.solve(D, c(1, 4:6)), c(1:2, 4:6), sum)
    
    yy <- dim(D)[[5]]
    nvar <- length(vars)
    
    # first term of product: (I + sum(D^-1) * lambda)^-1
    dsl.mat <- abind(d.inv.sum,
                     "lambda" = aperm(array(rep(L, yy), dim = c(dim(L), yy)), c(1:3, 5, 4)),
                     along = 2)
    d.inv.sum.lambda <- array(apply(dsl.mat, 3:5, 
                                    function(arr) arr[,1:nvar] %*% arr[,-(1:nvar)]), 
                              dim = dim(d.inv.sum))
    
    I.plus.d.l <- array.solve(sweep(d.inv.sum.lambda, 1:2, diag(nvar), "+"), solve.over = 3:5)
    
    # second term of product: sum(D^-1 * y_bar)
    dy.mat <- abind("y.bar" = aperm(Y, c(2,1,3:5)),
                    array.solve(D, c(1,4:6)),
                    along = 1)
    d.y.sum <- apply(apply(dy.mat, 3:6, function(arr) arr[-1,] %*% arr[1,]), c(1, 3:5), sum)
    
    
    # multiply two elements
    Idl.dy.mat <- abind(d.y.sum, I.plus.d.l, along = 1)
    Idl.dy <- apply(Idl.dy.mat, 3:5, function(arr) arr[-1,] %*% arr[1,])
    
    # multiply by S
    S.Idl.dy.mat <- abind(Idl.dy, S, along = 1)
    S.Idl.dy <- apply(S.Idl.dy.mat, 3:5, function(arr) arr[1,] %*% arr[-1,])
    
    # subtract eta
    c.tau <- sweep(S.Idl.dy, c(1:2,4), Eta, "-")
    dimnames(c.tau) <- dimnames(D[1,1,,,,])
    return(c.tau)
}



#' Fit model
#'
#' Fit a multivariate normal model to all available data
#' @param varbls Vector of variables to be included. Default is to include all variables.
#' @return List containing vector of means tau and variance/covariance matrix S
#' @export
#' 
run.model <- function(varbls = dimnames(obs)[[1]]) {
    
    ecmwf.cov <- se.covariances(offset.forecast(ecmwf))
    ncep.cov <- se.covariances(offset.forecast(ncep))
    ukmo.cov <- se.covariances(offset.forecast(ukmo))
    
    c.sigma <- se.sigma()
    c.lambda <- se.lambda()
    
    ens.mean.error <- apply(forecast.errors(superensemble())[,,,,-(1:3)], 1:4, mean)
    c.eta <- apply(ens.mean.error, c(1,2,4), mean)
    
    y.bar <- abind("ecmwf" = apply(offset.forecast(ecmwf[,,,,-1]), 1:4, mean),
                   "ncep" = apply(offset.forecast(ncep[,,,,-1]), 1:4, mean),
                   "ukmo" = apply(offset.forecast(ukmo[,,,,-1]), 1:4, mean),
                   along = 0)
    
    d <- abind("ecmwf" = ecmwf.cov / (dim(ecmwf)[[5]] - 1) + c.sigma,
               "ncep" = ncep.cov / (dim(ncep)[[5]] - 1) + c.sigma,
               "ukmo" = ukmo.cov / (dim(ukmo)[[5]] - 1) + c.sigma,
               along = 0)
    
    m.precision <- se.precision(d, c.lambda, vars = varbls)
    m.s <- array.solve(m.precision, 3:5)
    m.tau <- se.tau(d, c.lambda, y.bar, m.s, c.eta)
    
    return(list(tau = m.tau, s = m.s))
}