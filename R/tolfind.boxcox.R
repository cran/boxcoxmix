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
#' @param plot.opt Set \code{plot.opt=2}, to plot the EM trajectories and the development of the disparity over
#' iteration number.  And \code{plot.opt=0}, for none of them.
#' @param s number of points in the grid search of \code{tol}.
#' @param steps maximum number of iterations for the EM algorithm.
#' @param find.in.range search in a range of \code{tol}, with default (0,1.5) in
#' step of 0.1 .
#' @param start a description of the initial values to be used in the fitted
#' model, Quantile-based version "quantile" or Gaussian Quadrature "gq" can be
#' set.
#' @param verbose If set to FALSE, no printed output on progress.
#' @param noformat Set \code{noformat = TRUE}, to change the formatting of the plots.
#' @param \dots extra arguments will be ignored.
#' @return 
#' List with class \code{boxcoxmix} containing:
#' \item{MinDisparity}{the minimum disparity found.} \item{Mintol}{the
#' value of \code{tol} corresponding to \code{MinDisparity}.}
#' \item{AllDisparities }{a vector containing all disparities calculated on the
#' grid.} \item{Alltol }{list of \code{tol} values used in the grid.}
#' \item{AllEMconverged }{1 is TRUE, means the EM algorithm converged.}
#' \item{aic}{the Akaike information criterion of the fitted regression model.}
#' \item{bic}{the Bayesian information criterion of the fitted regression model.}
#' @author Amani Almohaimeed and Jochen Einbeck
#' @seealso \code{\link{np.boxcoxmix}}.
#' @keywords tolfind boxcox
#' @examples
#' # The Pennsylvanian Hospital Stay Data
#' data(hosp, package = "npmlreg")
#' test1 <- tolfind.boxcox(duration ~ age , data = hosp, K = 2, lambda = 0, 
#'            find.in.range = c(0, 2), s = 10,  start = "gq")
#' # Minimal Disparity: 137.8368 at tol= 2 
#' # Minimal Disparity with EM converged: 137.8368 at tol= 2
#' 
#' # Effect of Phenylbiguanide on Blood Pressure
#' \donttest{data(PBG, package = "nlme")
#' test2 <- tolfind.boxcox(deltaBP ~ dose , groups = PBG$Rabbit, find.in.range = c(0, 2),
#'     data = PBG, K = 2, lambda = -1, s = 15,  start = "quantile", plot.opt = 0)
#' test2$Mintol
#' # [1] 1.6
#' test2$MinDisparity
#' # [1] 449.5876}
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' @export tolfind.boxcox
tolfind.boxcox<- function (formula, groups=1, data, K=3,  lambda=1, EMdev.change = 1e-04, plot.opt = 2, s = 15, steps=500,
                           find.in.range = c(0, 1.5), start="gq", verbose = FALSE, noformat = FALSE, ...) 
{
	call <- match.call()
	if (K == 1) {
		stop("Please choose K > 1.")
	}
	tol <- find.in.range[1] + (find.in.range[2] - find.in.range[1]) * 
		1:s/s
	all.Disparities <-  all.converged <- rep(0, s)
	min.Disp <- min.Disp.conv <- 10^8
	s.min <- step.min.conv <- 1
	if (!noformat) {
		oldpar <- graphics::par(no.readonly = TRUE)
		on.exit(graphics::par(oldpar))
		if (steps > 8) 
			graphics::par(mfrow = c(4, 4), cex = 0.5)
		else graphics::par(mfrow = c(3, 3), cex = 0.5, cex.axis = 1.1)
	}
	for (t in 1:s) {
		# tol <- find.in.range[1] + (find.in.range[2] - find.in.range[1]) * 
		#   t/s
		# tol0 <- tol
		message("Executing NPML estimation of random effect and variance component model for tol = ", tol[t], " in range (", find.in.range[1], ", ", find.in.range[2], "].")
		fit <- np.boxcoxmix(formula=formula, groups= groups, data=data, K=K, lambda=lambda,steps= steps, 
				    tol = tol[t], start=start, EMdev.change = EMdev.change, 
				    plot.opt = plot.opt, verbose = verbose)
		message("Step for tol = ", tol[t], ": done!")
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
	fit <- np.boxcoxmix(formula=formula, groups= groups, data=data, K=K, lambda=lambda,steps= steps, 
			    tol = tol.min, start=start, EMdev.change = EMdev.change, 
			    plot.opt = 0, verbose = verbose)
	aic<- fit$aic
	bic<- fit$bic
	npcolors <- 2 + all.converged
	if (plot.opt==2){
		graphics::plot(tol, all.Disparities, type = "o", xlab = "tol", 
			       ylab = "Disparity", col= npcolors)
		graphics::segments(tol.min, min.Disp, tol.min, 1.1 * min.Disp, col = 4)
		message("Minimal Disparity: ", min.Disp, " at tol = ", tol.min, "\n")
		if (max(all.converged) == 1) {
			message("Minimal Disparity with EM converged: ", min.Disp.conv, " at tol = ", tol.min.conv, "\n")
		}
		else {
			message("No convergence achieved for any choice of tol.", "\n")
		}
		result <- list("call"=call, fit=fit, aic=aic,bic=bic, MinDisparity = min.Disp.conv, Mintol = tol.min.conv, 
			       AllDisparities = all.Disparities, Alltol = tol, 
			       AllEMconverged = all.converged == TRUE, "EMiteration"=iter,
			       "kind"=2, "npcolors"=npcolors, "tol.min"=tol.min, "min.Disp"=min.Disp)
		class(result)<-"boxcoxmix"
	}else{
		result <- list("call"=call, fit=fit, aic=aic, bic=bic, MinDisparity = min.Disp.conv, Mintol = tol.min.conv, 
			       AllDisparities = all.Disparities, Alltol = tol, 
			       AllEMconverged = all.converged == TRUE, "EMiteration"=iter,
			       "kind"=2, "npcolors"=npcolors, "tol.min"=tol.min, "min.Disp"=min.Disp)
		class(result)<-"boxcoxmix"
	}
	return(result)
}
