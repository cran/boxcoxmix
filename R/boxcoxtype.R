# Box-Cox-type link function


#' Box-Cox-type link function for logistic mixed-effects Models
#' @rdname boxcoxtype 
#' @import "npmlreg"
#' 
#' 
#' 
#' @description The \code{boxcoxtype()} performs a grid search over the parameter \code{Lambda}
#' for logistic mixed-effects models and then optimizes over this grid,
#' to calculate the maximum likelihood estimator of the transformation. 
#' 
#' 
#' 
#'
#' @details The Box-Cox transformation (Box & Cox, 1964) is applied to the
#' logistic mixed-effects models with an unspecified
#' mixing distribution. The NPML estimate of the mixing distribution is known
#' to be a discrete distribution involving a finite number of mass-points and corresponding
#' masses (Aitkin et al., 2009). An Expectation-Maximization (EM) algorithm is
#' used for fitting the finite mixture distribution, one needs to specify the
#' number of components \code{k} of the finite mixture in advance.
#' This algorithm can be implemented using the npmlreg function \code{\link{alldist}}
#' for the logistic-type overdispersion model and the npmlreg function \code{\link{allvc}} for the 
#' two-level logistic-type model, setting \code{family = binomial(link = boxcoxpower(Lambda))} where 
#' \code{Lambda} is the value of the power transformation. When \code{k}=1, the npmlreg function \code{alldist()} 
#' fits the logistic regression model without random effects. 
#'  
#' 
#' 
#' \code{boxcoxtype()} performs a grid search over the parameter \code{Lambda} and then
#' optimizes over this grid, to calculate the maximum likelihood estimator of the transformation.
#' It produces a plot of the profile likelihood function that summarises information
#' concerning \code{Lambda}, including a vertical line indicating the best value of \code{Lambda}
#' that maximizes the profile log-likelihood.
#' 
#' @param formula a formula describing the transformed response and the fixed
#' effect model (e.g. y ~ x).
#' @param random a formula defining the random model. Set \code{random= ~1} to model
#' logistic-type overdispersion model. For a two-level logistic-type model,
#' set \code{random= ~1|groups}, where groups are at the upper level.
#' @param data a data frame containing variables used in the fixed and random
#' effect models.
#' @param k the number of mass points.
#' @param trials optional prior weights for the data. For Bernoulli distribution, set trials=1.
#' @param find.in.range search in a range of \code{Lambda}, with default (-2,2)
#' in step of 0.1.
#' @param s number of points in the grid search of \code{Lambda}.
#' @param random.distribution the mixing distribution, Gaussian Quadrature (gq) or NPML (np) can be set.
#' @param plot.opt Set \code{plot.opt=1}, to plot the profile log-likelihood against \code{Lambda}.
#' if \code{plot.opt=0}, no plot is printed.
#' @param \dots extra arguments will be ignored.
#' @return
#' List with class \code{boxcoxmix} containing:
#' \item{Maximum}{the best estimate of \code{Lambda} found.} \item{objective}{the value of the profile
#' log-likelihood corresponding to Maximum.} 
#' \item{coef}{the vector of coefficients.} 
#' \item{profile.loglik}{the profile log-likelihood of the fitted regression model.}
#' \item{fit}{the fitted alldist object from the last EM iteration.}
#' \item{aic}{the Akaike information criterion of the fitted regression model.} 
#' \item{bic}{the Bayesian information criterion of the fitted regression model.}
#' 
#' The other outcomes are not relevant to users and they are intended for internal use only.
#'  
#' @author Amani Almohaimeed and Jochen Einbeck
#' @seealso \code{\link{np.boxcoxmix}}, \code{\link{optim.boxcox}},
#' \code{\link{tolfind.boxcox}}, \code{\link{Kfind.boxcox}}.
#' 
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
#' @keywords boxcoxtype 
#' @examples
#' #Beta blockers data
#' data("betablocker", package = "flexmix")
#' library(npmlreg)
#' betavc <-allvc(cbind(Deaths, Total - Deaths) ~ Treatment, data = betablocker,random=~1|Center,
#'  k=3,random.distribution='np',family = binomial(link = boxcoxpower(0)))
#' betavc$disparity
#' #[1] 318.7211
#' betavc3 <-boxcoxtype(cbind(Deaths, Total - Deaths) ~ Treatment,random=~1|Center, 
#' data = betablocker, find.in.range = c(-2,0.4),s=40,k=3,random.distribution='np')
#' #Maximum Profile Log-likelihood: -158.6025 at lambda= -0.56
#' betavc3$fit$disparity
#' #[1] 317.2049
#' betavc3$aic
#' #[1] 331.2049
#' betavc3$bic
#' #[1] 343.6942
#' 
#' @export 
boxcoxtype <- function(formula,random = ~1, k=3,trials=1, data,find.in.range = c(-2,2), s = 20, plot.opt=1,random.distribution='np',...) {
	call <- match.call()
	N <- NROW(data)
	weights= rep(trials, N)
	data$weights <- weights
	S <- 0:s
	lam <- find.in.range[1] + (find.in.range[2] - find.in.range[1]) * 
		S/s
	loglik <- rep(0, s+1 )
	#if (plot.opt==1){graphics::par(mfrow=c(1,1))}
	mform <- strsplit(as.character(random)[2], "\\|")[[1]]
	mform <- gsub(" ", "", mform)
	for (t in 1:(s+1)) {
		if (length(mform) == 1) {
			message("Executing NPML estimation of random effect model accounting for overdispersion for lambda = ", lam[t], " in range ", find.in.range[1], "; ", find.in.range[2], ".")
			fit <- alldist(formula = formula,random = formula(random), k=k,weights = weights,data = data, family=binomial(link=boxcoxpower(lam[t])),
				       verbose=FALSE, plot.opt = 0)
			message("Step for lambda = ", lam[t], ": done!")
			loglik[t] <- -1/2*(fit$disparity)
		}else {
			message("Executing NPML estimation of random effect variance component model for lambda = ", lam[t], " in range ", find.in.range[1], "; ", find.in.range[2], ".")
			fit <- allvc(formula = formula,random = formula(random), k=k,weights = weights,data = data, family=binomial(link=boxcoxpower(lam[t])),
				     verbose=FALSE, plot.opt = 0)
			message("Step for lambda = ", lam[t], ": done!")
			loglik[t] <- -1/2*(fit$disparity)
		}
	}
	max.result <- which.max(loglik)
	maxloglik<-max(loglik)
	lambda.max <- lam[max.result]
	if(length(mform) == 1){
		fit <- alldist(formula = formula,random = formula(random), k=k, weights=weights, data = data, 
			       family = binomial(link=boxcoxpower(lambda.max )), verbose=FALSE, random.distribution = random.distribution, plot.opt = 0)
	}else{
		fit <- allvc(formula =formula,random = formula(random), k=k, weights=weights, data = data,
			     family=binomial(link=boxcoxpower(lambda.max )), verbose=FALSE, random.distribution = random.distribution, plot.opt = 0)
	}
	if(k==1){
		coef<-fit$coefficients
	}else{
		coef<-fit$coefficients[1:abs(length(fit$coefficients)-length(fit$mass.points))]
	}
	ylim1 = range(loglik, maxloglik)
	ml <- paste("Maximum profile log-likelihood:",round(maxloglik) , "at lambda=", lambda.max, "\n")
	if (plot.opt==1){
		oldpar <- graphics::par(no.readonly = TRUE)
		on.exit(graphics::par(oldpar))
		graphics::plot(lam, loglik, type = "l", xlab = expression(lambda), 
			       ylab = "Profile log-likelihood", ylim = ylim1,col="green", main=ml)
		plims <- graphics::par("usr")
		y0 <- plims[3]
		graphics::segments(lambda.max, y0, lambda.max, maxloglik, lty = 1,col="red")
		message("Maximum Profile Log-likelihood: ", maxloglik, " at lambda = ", lambda.max, "\n")
	}
	aic<- fit$disparity+2*(length(coef)+2*k)
	bic<-fit$disparity+log(N)*(length(coef)+2*k)
	result <- list(coef=coef,fit=fit, All.lam=lam, profile.loglik = loglik, lastfit = list(fit), objective = maxloglik,
		       Maximum = lambda.max,y0=y0,ml=ml, "kind"=5, ylim1=ylim1, aic = aic, bic = bic)
	class(result)<-"boxcoxmix"
	return(result)
}

### Here's the Box-Cox Transformation is applied to the logistic model as a link function:
#' @rdname boxcoxtype
#' @param Lambda the power of the transformation
#' @export
boxcoxpower <- function(Lambda=0)
{
	linkfun <- function(mu){   if(Lambda==0) log(mu/(1-mu)) 
	else  (((mu/(1-mu))^Lambda)-1)/Lambda}
	linkinv <- function(eta) {  if(Lambda==0)  plogis(eta) 
	else (((1+Lambda*eta)^(-1/Lambda))+1)^(-1)}
	mu.eta<- function(eta) { if(Lambda==0)  ifelse(abs(eta)>30,.Machine$double.eps,
						       exp(eta)/(1+exp(eta))^2) else
							       pmax(((1+Lambda*eta)^((1/Lambda)-1))/((1+Lambda*eta)^(1/Lambda)+1)^2, 
								    .Machine$double.eps)}
	valideta <- function(eta) TRUE
	link <-paste("boxcoxpower(",Lambda,")", sep="")
	structure(list(linkfun = linkfun, linkinv = linkinv,
		       mu.eta = mu.eta, valideta = valideta, 
		       name = link),
		  class = "link-glm")
}



# Modified binomial family function (that allows boxcoxbower link function)
#' @rdname boxcoxtype
#' @param link the link function to be used.
#' @import "stats"
#' 
#' @export
binomial <- function (link = boxcoxpower(0))
{
	linktemp <- substitute(link)
	if (!is.character(linktemp)) linktemp <- deparse(linktemp)
	okLinks <- c("logit", "probit", "cloglog", "cauchit", "log", "boxcoxpower")
	if (linktemp %in% okLinks)
		stats <- make.link(linktemp)
	else if (is.character(link)) {
		stats <- make.link(link)
		linktemp <- link
	} else {
		## what else shall we allow?  At least objects of class link-glm.
		if(inherits(link, "link-glm")) {
			stats <- link
			if(!is.null(stats$name)) linktemp <- stats$name
		} else {
			stop(gettextf('link "%s" not available for binomial family; available links are %s',
				      linktemp, paste(sQuote(okLinks), collapse =", ")),
			     domain = NA)
		}
	}
	variance <- function(mu) mu * (1 - mu)
	validmu <- function(mu) all(is.finite(mu)) && all(mu>0 &mu<1)
	dev.resids <- stats::binomial()$dev.resids #function(y, mu, wt) .Call(stats:::C_binomial_dev_resids, y, mu, wt)
	#partialAIC <- function(y, n, mu, wt, dev)NA
	#finalAIC <- function(paic, n, dev)NA
	aic <- function(y, n, mu, wt, dev) {
		m <- if(any(n > 1)) n else wt
		-2*sum(ifelse(m > 0, (wt/m), 0)*
		       dbinom(round(m*y), round(m), mu, log=TRUE))
	}
	initialize <- expression({
					 if (NCOL(y) == 1) {
						 ## allow factors as responses
						 ## added BDR 29/5/98
						 if (is.factor(y)) y <- y != levels(y)[1L]
						 n <- rep.int(1, nobs)
						 ## anything, e.g. NA/NaN, for cases with zero weight is OK.
						 y[weights == 0] <- 0
						 if (any(y < 0 | y > 1))
							 stop("y values must be 0 <= y <= 1")
						 mustart <- (weights * y + 0.5)/(weights + 1)
						 m <- weights * y
						 if(any(abs(m - round(m)) > 1e-3))
							 warning("non-integer #successes in a binomial glm!")
					 }
					 else if (NCOL(y) == 2) {
						 if(any(abs(y - round(y)) > 1e-3))
							 warning("non-integer counts in a binomial glm!")
						 n <- y[, 1] + y[, 2]
						 y <- ifelse(n == 0, 0, y[, 1]/n)
						 weights <- weights * n
						 mustart <- (n * y + 0.5)/(n + 1)
					 }
					 else 
						 stop("for the 'binomial' family, y must be a vector of 0 and 1\'s\nor a 2 column 
						   matrix where col 1 is no. successes and col 2 is no. failures")
	})
						 
	simfun <- function(object, nsim){
		ftd <- fitted(object)
		n <- length(ftd)
		ntot <- n*nsim
		wts <- object$prior.weights
		if (any(wts %% 1 != 0))
			stop("cannot simulate from non-integer prior.weights")
		## Try to fathom out if the original data were
		## proportions, a factor or a two-column matrix
		if (!is.null(m <- object$model)) {
			y <- model.response(m)
			if(is.factor(y)) {
				## ignote weights
				yy <- factor(1+rbinom(ntot, size = 1, prob = ftd),
					     labels = levels(y))
				split(yy, rep(seq_len(nsim), each = n))
			} else if(is.matrix(y) && ncol(y) == 2) {
				yy <- vector("list", nsim)
				for (i in seq_len(nsim)) {
					Y <- rbinom(n, size = wts, prob = ftd)
					YY <- cbind(Y, wts - Y)
					colnames(YY) <- colnames(y)
					yy[[i]] <- YY
				}
				yy
			} else
				rbinom(ntot, size = wts, prob = ftd)/wts
		} else rbinom(ntot, size = wts, prob = ftd)/wts
	}
	structure(list(family = "binomial",
		       link = linktemp,
		       linkfun = stats$linkfun,
		       linkinv = stats$linkinv,
		       variance = variance,
		       dev.resids = dev.resids,
		       aic = aic,
		       mu.eta = stats$mu.eta,
		       initialize = initialize,
		       validmu = validmu,
		       valideta = stats$valideta,
		       simulate = simfun),
		  class = "family")
  }

