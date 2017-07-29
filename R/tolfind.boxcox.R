## A grid search over tol






































#' Grid search over tol for NPPML estimation of random effect and variance component models
#' 
#' A grid search over the parameter \code{tol}, to set the initial values of
#' the EM algorithm.
#' 
#' A grid search over \code{tol} can be performed using \code{tolfind.boxcox()}
#' function, which works for \code{np.boxcoxmix()} to find the
#' optimal solution.
#' 
#' @param formula a formula describing the transformed response and the fixed
#' effect model (e.g. y ~ x).
#' @param groups the random effects. To fit overdispersion models , set \code{groups} = 1.
#' @param data a data frame containing variables used in the fixed and random
#' effect models.
#' @param K the number of mass points.
#' @param lambda a transformation parameter, setting \code{lambda}=1 means 'no
#' transformation'.
#' @param EMdev.change a small scalar, with default 0.0001, used to determine
#' when to stop EM algorithm.
#' @param plot.opt Set \code{plot.opt=2}, to plot the disparity against
#' iteration number. 
#' @param s number of points in the grid search of \code{tol}.
#' @param steps maximum number of iterations for the EM algorithm.
#' @param find.in.range search in a range of \code{tol}, with default (0,2) in
#' step of 0.1 .
#' @param start a description of the initial values to be used in the fitted
#' model, Quantile-based version "quantile" or Gaussian Quadrature "gq" can be
#' set.
#' @param verbose If set to FALSE, no printed output on progress.
#' @param noformat Set \code{noformat = TRUE}, to change the formatting of the plots.
#' @param \dots extra arguments will be ignored.
#' @return \item{MinDisparity}{the minimum disparity found.} \item{Mintol}{the
#' value of \code{tol} corresponding to \code{MinDisparity}.}
#' \item{AllDisparities }{a vector containing all disparities calculated on the
#' grid.} \item{Alltol }{list of \code{tol} values used in the grid.}
#' \item{AllEMconverged }{1 is TRUE, means the EM algorithm converged.}
#' @author Amani Almohaimeed and Jochen Einbeck
#' @seealso \code{\link{np.boxcoxmix}}.
#' @keywords tolfind boxcox
#' @examples
#' # The Pennsylvanian Hospital Stay Data
#' data(hosp, package = "npmlreg")
#' test1 <- tolfind.boxcox(duration ~ age , data = hosp, K = 2, lambda = 0, 
#'     s = 10, steps = 600, start = "quantile", plot.opt = 0)
#' # Minimal Disparity: 137.8364 at tol= 1.4  
#' # Minimal Disparity with EM converged: 137.8364 at tol= 1.4 
#' 
#' # Effect of Phenylbiguanide on Blood Pressure
#' data(PBG, package = "nlme")
#' test3 <- tolfind.boxcox(deltaBP ~ dose , groups = PBG$Rabbit, 
#'     data = PBG, K = 2, lambda = -1, s = 15, steps = 500, start = "quantile", plot.opt = 0)
#' # Minimal Disparity: 449.5876 at tol= 1.6 
#' # Minimal Disparity with EM converged: 449.5876 at tol= 1.6
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' @export tolfind.boxcox
tolfind.boxcox<- function (formula, groups=1, data, K=3,  lambda=1, EMdev.change = 1e-04, plot.opt = 2, s = 20, steps=500,
                           find.in.range = c(0, 2), start="gq", verbose = FALSE, noformat = FALSE, ...) 
{
  call <- match.call()
  if (K == 1) {
    stop("Please choose K > 1.")
  }
  all.Disparities <-  all.converged <- rep(0, s)
  min.Disp <- min.Disp.conv <- 10^8
  s.min <- step.min.conv <- 1
  if (!noformat) {
    if (steps > 8) 
      graphics::par(mfrow = c(4, 4), cex = 0.5)
    else graphics::par(mfrow = c(3, 3), cex = 0.5, cex.axis = 1.1)
  }
  for (t in 1:s) {
    tol <- find.in.range[1] + (find.in.range[2] - find.in.range[1]) * 
      t/s
    tol0 <- tol
      fit <- try(np.boxcoxmix(formula=formula, groups= groups, data=data, K=K, lambda=lambda,steps= steps, 
                           tol = tol, start=start, EMdev.change = EMdev.change, 
                           plot.opt = plot.opt, verbose = verbose))
      if (class(fit) == "try-error") {
        cat("boxcox.tolfind failed using tol=", tol, ". Hint:  specify another range of tol values and try again. ")
        return()
      }
      all.Disparities[t] <-fit$disparity
      all.converged[t] <- fit$EMconverged
      iter <- fit$EMiteration
    if (all.Disparities[t] < min.Disp) {
      min.Disp <- all.Disparities[t]
      s.min <- t
    }
    if (all.Disparities[t] < min.Disp.conv && all.converged[t]) {
      min.Disp.conv <- all.Disparities[t]
      step.min.conv <- t
      iter <- fit$EMiteration
    }
  } 
  tol.min <- find.in.range[1] + (find.in.range[2] - find.in.range[1]) * 
    s.min/s
  tol.min.conv <- find.in.range[1] + (find.in.range[2] - find.in.range[1]) * 
    step.min.conv/s
  npcolors <- 2 + all.converged
  if (plot.opt==2){
    graphics::plot(find.in.range[1] + (find.in.range[2] - find.in.range[1]) * 
                     (1:s)/s, all.Disparities, type = "o", xlab = "tol", 
                   ylab = "Disparity", col= npcolors)
    graphics::segments(tol.min, min.Disp, tol.min, 1.1 * min.Disp, col = 4)
  }
  cat("Minimal Disparity:", min.Disp, "at tol=", tol.min, "\n")
  if (max(all.converged) == 1) {
    cat("Minimal Disparity with EM converged:", min.Disp.conv, 
        "at tol=", tol.min.conv, "\n")
  }
  else {
    cat(" No convergence achieved for any choice of tol.", 
        "\n")
  }
  result <- list( "call"=call, MinDisparity = min.Disp.conv, Mintol = tol.min.conv, 
                  AllDisparities = all.Disparities, Alltol = find.in.range[1] + 
                    (find.in.range[2] - find.in.range[1]) * (1:s)/s, 
                  AllEMconverged = all.converged == TRUE, "EMiteration"=iter,
                  "kind"=2, "npcolors"=npcolors, "tol.min"=tol.min, "min.Disp"=min.Disp)
  class(result)<-"boxcoxmix"
  return(result)
}
