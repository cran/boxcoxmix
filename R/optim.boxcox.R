## A grid search over lambda






































#' Response Transformations for Random Effect and Variance Component Models
#' 
#' @rdname optim.boxcox 
#' @description The \code{optim.boxcox()} performs a grid search over the parameter \code{lambda}
#' for overdispersed generalized linear models and variance component models and 
#' then optimizes over this grid, to calculate the maximum likelihood estimator of
#' the transformation. 
#' 
#' 
#'
#' @details The Box-Cox transformation (Box & Cox, 1964) is applied to the overdispersed
#' generalized linear models and variance component models with an unspecified
#' mixing distribution. The NPML estimate of the mixing distribution is known
#' to be a discrete distribution involving a finite number of mass-points and corresponding
#' masses (Aitkin et al., 2009). An Expectation-Maximization (EM) algorithm is
#' used for fitting the finite mixture distribution, one needs to specify the
#' number of components \code{K} of the finite mixture in advance. To stop the EM-algorithm when it reached its convergence point,
#' we need to defined the convergence criteria that is the absolute change in 
#' the successive log-likelihood function values being less than an arbitrary 
#' parameter such as \code{EMdev.change} = 0.0001 (Einbeck et
#' at., 2014). This algorithm can be implemented using
#' the function \code{np.boxcoxmix()}, which is designed to account for overdispersed generalized
#' linear models and variance component models using the non-parametric
#' profile maximum likelihood (NPPML) estimation.
#' 
#'  
#' The ability of the EM algorithm to locate the global maximum in fewer iterations
#' can be affected by the choice of initial values, the function \code{optim.boxcox()} 
#' allows us to choose from two different methods to set the initial value of the mass
#' points. When option "gq" is set, then Gauss-Hermite masses and mass points are used
#' as starting points in the EM algorithm, while setting start= "quantile" uses the 
#' Quantile-based version to select the starting points. 
#' 
#' 
#' \code{optim.boxcox()} performs a grid search over the parameter \code{lambda} and then
#' optimizes over this grid, to calculate the maximum likelihood estimator of the transformation.
#' It produces a plot of the non-parametric profile likelihood function that summarises information
#' concerning \code{lambda}, including a vertical line indicating the best value of \code{lambda}
#' that maximizes the non-parametric profile log-likelihood.
#' 
#' @param formula a formula describing the transformed response and the fixed
#' effect model (e.g. y ~ x).
#' @param groups the random effects. To fit overdispersion models, set \code{groups} = 1.
#' @param data a data frame containing variables used in the fixed and random
#' effect models.
#' @param K the number of mass points.
#' @param steps maximum number of iterations for the EM algorithm.
#' @param tol a positive scalar (usually, 0<\code{tol} <= 2)
#' @param start a description of the initial values to be used in the fitted
#' model, Quantile-based version "quantile" or Gaussian Quadrature "gq" can be
#' set.
#' @param EMdev.change a small scalar, with default 0.0001, used to determine
#' when to stop EM algorithm.
#' @param find.in.range search in a range of \code{lambda}, with default (-3,3)
#' in step of 0.1.
#' @param s number of points in the grid search of \code{lambda}.
#' @param plot.opt Set \code{plot.opt=3}, to plot the disparity against
#' iteration number and the profile log-likelihood against \code{lambda}. 
#' Use \code{plot.opt=0}, to only plot the profile log-likelihood against \code{lambda}.
#' @param verbose If set to FALSE, no printed output on progress.
#' @param noformat Set \code{noformat = TRUE}, to change the formatting of the plots.
#' @param \dots extra arguments will be ignored.
#' @return 
#' \item{All.lambda}{list of \code{lambda} values used in the grid.}
#' \item{Maximum}{the best estimate of \code{lambda} found.} \item{objective}{the value of the profile
#' log-likelihood corresponding to Maximum.} 
#' \item{EMconverged}{1 is TRUE, means the EM algorithm converged.} \item{EMiteration}{provides the number of iterations of the EM algorithm.} 
#' \item{mass.point}{the fitted mass points.} \item{p}{the masses corresponding to the mixing proportions.} \item{beta}{the
#' vector of coefficients.} \item{sigma}{the standard deviation of the mixing distribution (the square root of the variance).} \item{se}{the standard error of the estimate.}
#' \item{w}{a matrix of posterior probabilities that element i comes from cluster k.}
#' \item{loglik}{the profile log-likelihood of the fitted regression model.}
#' \item{profile.loglik}{the profile complete log-likelihood of the fitted regression model.}
#' \item{disparity}{the disparity of the fitted regression model.} 
#' \item{call}{the matched call.} \item{formula}{the formula provided.}
#' \item{data}{the data argument.} \item{aic}{the Akaike information criterion of the fitted regression model.}
#' \item{fitted}{the fitted values for the individual observations.} \item{fitted.transformed}{the fitted values for
#' the individual transformed observations.} \item{residuals}{the difference between the observed values and the fitted values.}
#' \item{residuals.transformed}{the difference between the transformed observed values and the transformed fitted values.}
#' \item{predicted.re}{a vector of predicted residuals.}
#' 
#' The other outcomes are not relevant to users and they are intended for internal use only.
#'  
#' @author Amani Almohaimeed and Jochen Einbeck
#' @seealso \code{\link{np.boxcoxmix}},
#' \code{\link{tolfind.boxcox}}.
#' @references 
#' Box G. and Cox D. (1964). An analysis of transformations. Journal of
#' the Royal Statistical Society. Series B (Methodological), pages 211-252.
#' 
#' Aitkin, M. A., Francis, B., Hinde, J., and Darnell, R. (2009). Statistical
#' modelling in R. Oxford University Press Oxford.
#' 
#' Jochen Einbeck, Ross Darnell and John Hinde (2014). npmlreg:
#' Nonparametric maximum likelihood estimation for random effect
#' models. R package version 0.46-1.
#' 
#' @keywords optim boxcox 
#' @examples
#' # The strength Data
#' data(strength, package = "mdscore")
#' maxlam <- optim.boxcox(y ~ cut*lot, data = strength, K = 3,  
#'            start = "gq" ,  find.in.range = c(-2, 2), s = 5)
#' # Maximum profile log-likelihood: 33.6795 at lambda= -0.4  
#' 
#' \donttest{data(Oxboys, package = "nlme")
#' Oxboys$boy <- gl(26,9)
#' maxlamvc <- optim.boxcox(height ~  age, groups = Oxboys$boy,
#'                          data = Oxboys,   K = 2,  start = "gq",
#'                          find.in.range=c(-1.2,1), s=6, plot.opt = 0) 
#' maxlamvc$Maximum
#' #[1] -0.8333333
#' plot(maxlamvc,8)}
#' 
#' 
#' 
#' 
#' 
#' @export 
optim.boxcox <- function (formula, groups = 1, data, K = 3, steps = 500, tol = 0.5, 
                          start = "gq", EMdev.change = 1e-04, find.in.range = c(-3, 
                                                                                3), s = 60, plot.opt = 3, verbose = FALSE, noformat = FALSE, 
                          ...) 
{
  call <- match.call()
  mform <- strsplit(as.character(groups), "\\|")
  mform <- gsub(" ", "", mform)
  if (!noformat) {
    if (steps > 8) 
      graphics::par(mfrow = c(4, 4), cex = 0.5)
    else graphics::par(mfrow = c(3, 3), cex = 0.5, cex.axis = 1.1)
  }
  result <- disp <- loglik <- EMconverged <- rep(0, s + 1)
  S <- 0:s
  lambda <- find.in.range[1] + (find.in.range[2] - find.in.range[1]) * 
    S/s
  for (t in 1:(s + 1)) {
    fit <- try(np.boxcoxmix(formula = formula, groups = groups, 
                            data = data, K = K, lambda = lambda[t], steps = steps, 
                            tol = tol, start = start, EMdev.change = EMdev.change, 
                            plot.opt = plot.opt, verbose = verbose))
    if (class(fit) == "try-error") {
      cat("optim.boxcox failed using lambda=", lambda[t], 
          ". Hint:  specify another range of lambda values and try again.")
      return()
    }
    EMconverged[t] <- fit$EMconverged
    result[t] <- fit$loglik
    if (!all(is.finite(result[t]))) {
      print.plot <- FALSE
    }
  }
  s.max <- which.max(result)
  max.result <- result[s.max]
  lambda.max <- lambda[s.max]
  fit <- np.boxcoxmix(formula = formula, groups = groups, data = data, 
                      K = K, lambda = lambda.max, steps = steps, tol = tol, 
                      start = start, EMdev.change = EMdev.change, plot.opt = 0, 
                      verbose = verbose)
  W <- fit$w
  P <- fit$p
  se <- fit$se
  iter <- fit$EMiteration
  names(P) <- paste("MASS", 1:K, sep = "")
  Z <- fit$mass.point
  names(Z) <- paste("MASS", 1:K, sep = "")
  Beta <- fit$beta
  Sigma <- fit$sigma
  Disp <- fit$disparity
  Disparities <- fit$Disparities
  n <- NROW(data)
  if (fit$model == "pure") {
    if(K==1){
      aic <- Disp + 2 * 2 # 2 here is the intercept + (c=1)
      bic <- Disp + log(n) * 2 # 2 here is the intercept + (c=1)
    }
    else {
      aic <- Disp + 2 * (2 * K)
      bic <- Disp + log(n) * (2 * K)
    }}
  else {if(K==1){
    aic <- Disp + 2 * (length(Beta) +1) # 1 here is (c=1)
    bic <- Disp + log(n) * (length(Beta) +1) # 1 here is (c=1)
  }
    else {
      aic <- Disp + 2 * (length(Beta) + 2 * K)
      bic <- Disp + log(n) * (length(Beta) + 2 * K)
    }}
  y <- fit$y
  yt <- fit$yt
  fitted <- fit$fitted
  fitted.transformed <- fit$fitted.transformed
  masses <- fit$masses
  ylim <- fit$ylim
  residuals <- fit$residuals
  residuals.transformed <- fit$residuals.transformed
  predicted.re <- fit$predicted.re
  Class <- fit$Class
  xx <- fit$xx
  model <- fit$model
  Disp <- fit$disparity
  Disparities <- fit$Disparities
  Loglik <- fit$loglik
  npcolors <- "green"
  ylim1 = range(result, max.result)
  maxl <- paste("Maximum profile log-likelihood:", round(max.result, 
                                                         digits = 3), "at lambda=", round(lambda.max, digits = 2), 
                "\n")
  if (plot.opt == 3) {
    graphics::plot(lambda, result, type = "l", xlab = expression(lambda), 
                   ylab = "Profile log-likelihood", ylim = ylim1, col = "green")
    plims <- graphics::par("usr")
    y0 <- plims[3]
    graphics::segments(lambda.max, y0, lambda.max, max.result, 
                       lty = 1, col = "red", lwd = 2)
    cat("Maximum profile log-likelihood:", max.result, "at lambda=", 
        lambda.max, "\n")
    result <- list(call = call, y0 = y0, p = P, mass.point = Z, 
                   beta = Beta, sigma = Sigma, se = se, w = W, Disparities = Disparities, 
                   formula = formula, data = data, loglik = Loglik, 
                   aic = aic, bic = bic, masses = masses, y = y, yt = yt, 
                   All.lambda = lambda, profile.loglik = result, disparity = Disp, 
                   EMconverged = EMconverged, Maximum = lambda.max, 
                   mform = length(mform), ylim = ylim, fitted = fitted, 
                   Class = Class, fitted.transformed = fitted.transformed, 
                   predicted.re = predicted.re, residuals = residuals, 
                   residuals.transformed = residuals.transformed, objective = max.result, 
                   kind = 3, EMiteration = iter, ss = s, s.max = s.max, 
                   npcolor = npcolors, ylim1 = ylim1, xx = xx, maxl = maxl, 
                   model = model)
    class(result) <- "boxcoxmix"
  }
  else {
    result <- list(call = call, p = P, mass.point = Z, beta = Beta, 
                   sigma = Sigma, se = se, w = W, Disparities = Disparities, 
                   formula = formula, data = data, loglik = Loglik, 
                   aic = aic, bic = bic, masses = masses, y = y, yt = yt, 
                   All.lambda = lambda, profile.loglik = result, disparity = Disp, 
                   EMconverged = EMconverged, Maximum = lambda.max, 
                   mform = length(mform), ylim = ylim, fitted = fitted, 
                   Class = Class, fitted.transformed = fitted.transformed, 
                   predicted.re = predicted.re, residuals = residuals, 
                   residuals.transformed = residuals.transformed, objective = max.result, 
                   kind = 3, EMiteration = iter, ss = s, s.max = s.max, 
                   npcolor = npcolors, ylim1 = ylim1, xx = xx, maxl = maxl, 
                   model = model)
    class(result) <- "boxcoxmix"
  }
  return(result)
}