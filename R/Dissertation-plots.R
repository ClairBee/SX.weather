
#' Plot forecast errors
#' 
#' Produce a plot of errors for each forecast at each leadtime, for all variables
#' @param model Name of dataset to use (use \code{data{package = "SX.weather"}} to list available options)
#' @export
#' 
plot.forecast.errors <- function(model) {
    
    org.par <- par()
    
    # add ensemble mean to model
    m.plus <- abind("em" = apply(model[,,,,-1], 1:4, mean), model, along = 5)
    
    errors <- forecast.errors(m.plus)
    mean.errors <- apply(errors, c(1, 4, 5), mean)
    
    par(mfrow = c(2,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    invisible(sapply(dimnames(mean.errors)[[1]][c(3:5, 1:2)], function(varbl) {
        
        matplot(mean.errors[varbl,,dim(mean.errors)[[3]]:1], type = "l", lty = 1,
                col = c(rep(adjustcolor("grey", alpha = 0.5), dim(mean.errors)[[3]]-2), "blue", "black"),
                xlab = "", ylab = "", main = varbl)
    }))
    
    # add legend in final box
    plot.new()
    legend("left", lty = 1, col = c("blue", "black", adjustcolor("grey", alpha = 0.5)), bty = "n", cex = 1.1,
           legend = c("Control forecast", "Perturbed mean", "Perturbed members"))
    
    # add overall title
    mtext(paste0(toupper(toString(as.list(match.call())$model)), " mean error at each forecast lead time"),
          outer = TRUE, cex = 1)
    
    # reset device parameters to original values
    par(mfrow = org.par$mfrow, mar = org.par$mar, oma = org.par$oma)
}


#' Plot RMSE
#' 
#' Produce a plot of RMSE for each forecast at each leadtime, for all variables
#' @param model Name of dataset to use (use \code{data{package = "SX.weather"}} to list available options)
#' @export
#' 
plot.forecast.rmse <- function(model) {
    
    org.par <- par()
    
    # add ensemble mean to model
    m.plus <- abind("em" = apply(model[,,,,-1], 1:4, mean), model, along = 5)
    
    rmse <- forecast.rmse(m.plus)
    rmse <- abind("em2" = apply(rmse[,,-c(1:2)], 1:2, mean), rmse)
    
    par(mfrow = c(2,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    invisible(sapply(dimnames(rmse)[[1]][c(3:5, 1:2)], function(varbl) {
        
        matplot(rmse[varbl,,dim(rmse)[[3]]:1], type = "l", lty = 1,
                col = c(rep(adjustcolor("grey", alpha = 0.5), dim(rmse)[[3]]-3), "blue", "black", "red"),
                xlab = "", ylab = "", main = varbl)
    }))
    
    # add legend in final box
    plot.new()
    legend("left", lty = 1, col = c("red", "blue", "black", adjustcolor("grey", alpha = 0.5)), bty = "n", cex = 1.1,
           legend = c("Mean of perturbed RMSE", "Control forecast", "RMSE of perturbed mean", "Perturbed members"))
    
    # add overall title
    mtext(paste0(toupper(toString(as.list(match.call())$model)), " RMSE at each forecast lead time"),
          outer = TRUE, cex = 1)
    
    # reset device parameters to original values
    par(mfrow = org.par$mfrow, mar = org.par$mar, oma = org.par$oma)
}


#' Plot RMSE vs spread
#' 
#' Produce a plot of RMSE vs spread for each forecast at each leadtime, for all variables
#' @param model Name of dataset to use (use \code{data{package = "SX.weather"}} to list available options)
#' @export
#' 
plot.rmse.spread <- function(model) {
    
    org.par <- par()
    
    spread <- ensemble.spread(model[,,,,-1])
    
    # add ensemble mean to model
    m.plus <- abind("em" = apply(model[,,,,-1], 1:4, mean), model, along = 5)
    rmse <- forecast.rmse(m.plus)
    rmse <- abind("em2" = apply(rmse[,,-c(1:2)], 1:2, mean), rmse)
    
    dat <- abind("spread" = spread, rmse, along = 3)
    
    par(mfrow = c(2,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    invisible(sapply(dimnames(dat)[[1]][c(3:5, 1:2)], function(varbl) {
        
        matplot(dat[varbl,,dim(dat)[[3]]:1], type = "l", lty = c(rep(1, dim(dat)[[3]]-1), 2),
                col = c(rep(adjustcolor("grey", alpha = 0.5), dim(dat)[[3]]-4), "blue", "black", "red", "green3"),
                xlab = "", ylab = "", main = varbl)
    }))
    
    # add legend in final box
    plot.new()
    legend("left", lty = c(2, rep(1, dim(dat)[[3]]-1)), col = c("green3", "red", "blue", "black", adjustcolor("grey", alpha = 0.5)),
           bty = "n", cex = 1.1,
           legend = c("Ensemble spread", "Mean of perturbed RMSE", "Control forecast RMSE",
                      "RMSE of perturbed mean", "Perturbed member RMSE"))
    
    # add overall title
    mtext(paste0(toupper(toString(as.list(match.call())$model)), " RMSE vs spread at each forecast lead time"),
          outer = TRUE, cex = 1)
    
    # reset device parameters to original values
    par(mfrow = org.par$mfrow, mar = org.par$mar, oma = org.par$oma)
    
}