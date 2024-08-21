#' Summary of boxcoxmix functions
#' 
#'
#' @description \code{summary()} and \code{print()} are generic functions used to produce the results of the functions:
#' \code{np.boxcoxmix()}, \code{optim.boxcox()} and \code{tolfind.boxcox()}.
#' @aliases summary.boxcoxmix print.boxcoxmix print.boxcoxmixpure
#' summary.boxcoxmixpure
#' @rdname summary
#' @method print boxcoxmix
#' @return Return invisibly the main object provided, while print a summary of its content. 
#' @export
print.boxcoxmix <- function(x,digits=max(3,getOption('digits')-3),na.print='', ...){
	if (x$kind == '1'){
		np <- length(x$beta)
		m <- seq(1,np)[substr(attr(x$beta,'names'),1,4)=='MASS']
		mass.points <- x$beta[m]
		cat('\nCall: ',deparse(x$call),'\n\n', fill = TRUE)
		cat('Coefficients', fill = TRUE)
		cat(":\n")
		if (is.na(x$beta)[np]){np<-np-1}
		K <- length(x$mass.point)
		names(x$beta)<- if(K==1) names(x$beta) else x$xx
		print.default(format(x$beta[1:np],digits = digits), print.gap = 2,quote = FALSE);cat('\n')
		cat('MLE of sigma:\t  ',
		    format(signif(x$sigma,digits)),'\n', fill = TRUE)
		p <- x$p
		cat('Mixture proportions')
		cat(":\n")
		print.default(format(p,digits),print.gap=2,quote=FALSE)
		cat('-2 log L:\t   ',
		    format(round(x$disparity,digits=1)),'and AIC = ',
		    format(round(x$aic,digits)),"\n", fill = TRUE)
	}
	if (x$kind == '2'){
		cat("Minimal Disparity with EM converged:",x$MinDisparity, 
		    "at tol=", x$Mintol, "\n", fill = TRUE)
		cat('MinDisparity')
		cat(":\n")
		print.default(format(x$MinDisparity,digits),quote=FALSE) 
		cat('Mintol')
		cat(":\n")
		print.default(format(x$Mintol,digits),quote=FALSE) 
		cat('AllDisparities')
		cat(":\n")
		print.default(format(x$AllDisparities,digits),print.gap=2,quote=FALSE) 
		cat('Alltol')
		cat(":\n")
		print.default(format(x$Alltol,digits),print.gap=2,quote=FALSE) 
		cat('AllEMconverged')
		cat(":\n")
		print.default(format(x$AllEMconverged,digits),print.gap=2,quote=FALSE) 
	}
	if (x$kind == '3'){
		cat("Maximum profile log-likelihood:", x$objective, "at lambda=", x$Maximum, "\n", fill = TRUE)
		np <- length(x$beta)
		m <- seq(1,np)[substr(attr(x$beta,'names'),1,4)=='MASS']
		mass.points <- x$beta[m]
		cat('\nCall: ',deparse(x$call),'\n\n', fill = TRUE)
		cat('Coefficients')
		cat(":\n")
		if (is.na(x$beta)[np]){np<-np-1}
		K <- length(x$mass.point)
		names(x$beta)<- if(K==1) names(x$beta) else x$xx
		print.default(format(x$beta[1:np],digits = digits), print.gap = 2,quote = FALSE);cat('\n')
		cat('MLE of sigma:\t  ',
		    format(signif(x$sigma,digits)),'\n', fill = TRUE)
		p <- x$p
		cat('Mixture proportions')
		cat(":\n")
		print.default(format(p,digits),print.gap=2,quote=FALSE)
		cat('-2 log L:\t   ',
		    format(round(x$disparity,digits=1)), 'and AIC = ',
		    format(round(x$aic,digits)),"\n", fill = TRUE)
	}
	invisible(x)
}
#' @rdname summary
#' @method print boxcoxmixpure
#' @export
print.boxcoxmixpure <- function(x, digits = max(3, getOption('digits') - 3), na.print = '', ...){
	if (x$kind == '1'){
		np <- length(x$mass.point)
		cat('\nCall: ',deparse(x$call),'\n\n', fill = TRUE)
		cat('Coefficients')
		cat(":\n")
		print.default(format(x$mass.point[1:np],digits = digits), print.gap = 2,quote = FALSE);cat('\n')
		cat('MLE of sigma:\t  ',
		    format(signif(x$sigma,digits)),'\n', fill = TRUE)
		p <- x$p
		cat('Mixture proportions')
		cat(":\n")
		print.default(format(p,digits),print.gap=2,quote=FALSE)
		cat('-2 log L:\t   ',
		    format(round(x$disparity,digits=1)),'and AIC = ',
		    format(round(x$aic,digits)),"\n", fill = TRUE)
	}
	if (x$kind == '2'){
		cat("Minimal Disparity with EM converged:",x$MinDisparity, 
		    "at tol=", x$Mintol, "\n", fill = TRUE)
		cat('MinDisparity')
		cat(":\n")
		print.default(format(x$MinDisparity,digits),quote=FALSE) 
		cat('Mintol')
		cat(":\n")
		print.default(format(x$Mintol,digits),quote=FALSE) 
		cat('AllDisparities')
		cat(":\n")
		print.default(format(x$AllDisparities,digits),print.gap=2,quote=FALSE) 
		cat('Alltol')
		cat(":\n")
		print.default(format(x$Alltol,digits),print.gap=2,quote=FALSE) 
		cat('AllEMconverged')
		cat(":\n")
		print.default(format(x$AllEMconverged,digits),print.gap=2,quote=FALSE) 
	}
	if (x$kind == '3'){
		cat("Maximum profile log-likelihood:", x$objective, "at lambda=", x$Maximum, "\n", fill = TRUE)
		np <- length(x$beta)
		m <- seq(1,np)[substr(attr(x$beta,'names'),1,4)=='MASS']
		mass.points <- x$beta[m]
		cat('\nCall: ',deparse(x$call),'\n\n', fill = TRUE)
		cat('Coefficients')
		cat(":\n")
		print.default(format(x$beta[1:np],digits = digits), print.gap = 2,quote = FALSE);cat('\n')
		cat('MLE of sigma:\t  ',
		    format(signif(x$sigma,digits)),'\n', fill = TRUE)
		p <- x$p
		cat('Mixture proportions')
		cat(":\n")
		print.default(format(p,digits),print.gap=2,quote=FALSE)
		cat('-2 log L:\t   ',
		    format(round(x$disparity,digits=1)), 'and AIC = ',
		    format(round(x$aic,digits)),"\n", fill = TRUE)
	}
	invisible(x)
}
#' @rdname summary
#' @method summary boxcoxmix
#' @param x an object for which a summary is desired. 
#' @param object an object for which a summary is desired.
#' @param digits an integer number format.
#' @param na.print a character string which is used to indicate NA values output format.
#' @param \dots additional arguments.
#' @export
summary.boxcoxmix <- function(object,digits=max(3,getOption('digits')-3), ...){
	if (object$kind == '1'){
		np <- length(object$beta)
		K <- length(object$mass.point)
		m <- seq(1,np)[substr(attr(object$beta,'names'),1,4)=='MASS']
		mass.points <- object$beta[m]
		MASS <- object$mass.point
		Std.Error <-object$se
		tvalue <- object$beta/Std.Error
		coef.table <- matrix(c(object$beta, Std.Error, tvalue), np, ncol=3, byrow = FALSE)
		rownames(coef.table) <-if(K==1)names(object$beta) else object$xx
		colnames(coef.table) <- c("Estimate", "Std. Error", "t value")
		cat('\nCall: ',deparse(object$call),'\n\n', fill = TRUE)
		cat('Coefficients')
		cat(":\n")
		print(coef.table)
		cat('\nMass points:\n')
		names(MASS) <- paste('MASS',seq(1,ncol(object$w)),sep='')
		print.default(format(MASS,digits),print.gap=2,quote=FALSE)
		p <- object$p
		cat('\nMixture proportions:\n')
		names(p) <- paste('MASS',seq(1,ncol(object$w)),sep='')
		print.default(format(p,digits),print.gap=2,quote=FALSE)
		cat('\nMLE of sigma:\t  ',format(signif(object$sigma,digits)),'\n', fill = TRUE)
		cat('-2 log L:\t ',format(round(object$disparity,digits=1)),  'and AIC =  ', format(round(object$aic,digits)), fill = TRUE)
		if (!is.null(object$w)) cat('     Convergence at iteration ',round(object$EMiteration,0), fill = TRUE)
		cat('\n')
	}
	if (object$kind == '2'){
		cat("Minimal Disparity with EM converged:",object$MinDisparity, 
		    "at tol=", object$Mintol, "\n", fill = TRUE)
		cat('MinDisparity')
		cat(":\n")
		print.default(format(object$MinDisparity,digits),quote=FALSE) 
		cat('Mintol')
		cat(":\n")
		print.default(format(object$Mintol,digits),quote=FALSE) 
		cat('AllDisparities')
		cat(":\n")
		print.default(format(object$AllDisparities,digits),print.gap=2,quote=FALSE) 
		cat('Alltol')
		cat(":\n")
		print.default(format(object$Alltol,digits),print.gap=2,quote=FALSE) 
		cat('AllEMconverged')
		cat(":\n")
		print.default(format(object$AllEMconverged,digits),print.gap=2,quote=FALSE) 
	}
	if (object$kind == '3'){
		cat("Maximum profile log-likelihood:", object$objective, "at lambda=", object$Maximum, "\n", fill = TRUE)
		np <- length(object$beta)
		m <- seq(1,np)[substr(attr(object$beta,'names'),1,4)=='MASS']
		mass.points <- object$beta[m]
		Std.Error <-object$se
		tvalue <- object$beta/Std.Error
		MASS <- object$mass.point
		coef.table <- matrix(c(object$beta, Std.Error, tvalue), np, ncol=3, byrow = FALSE)
		colnames(coef.table) <- c("Estimate", "Std. Error", "t value")
		rownames(coef.table) <- object$xx
		cat('\nCall: ',deparse(object$call),'\n\n', fill = TRUE)
		cat('Coefficients')
		cat(":\n")
		print(coef.table)
		cat('\nMass points:\n')
		print.default(format(MASS,digits),print.gap=2,quote=FALSE)
		p <- object$p
		cat('\nMixture proportions:\n')
		print.default(format(p,digits),print.gap=2,quote=FALSE)
		cat('\nMLE of sigma:\t  ',format(signif(object$sigma,digits)),'\n', fill = TRUE)
		cat('-2 log L:\t   ',format(round(object$disparity,digits=1)),  'and AIC =  ', format(round(object$aic,digits)), fill = TRUE)
		if (!is.null(object$w)) cat('     Convergence at iteration ',round(object$EMiteration,0), fill = TRUE)
		cat('\n')
	}
	invisible(object)
}
#' @rdname summary
#' @method summary boxcoxmixpure
#' @export
summary.boxcoxmixpure <- function(object, digits = max(3, getOption('digits') - 3), ...){
	if (object$kind == '1'){
		k <- length(object$mass.point)
		mass.points <- object$mass.point
		Std.Error <-object$se
		tvalue <- object$mass.point/Std.Error
		coef.table <- matrix(c(mass.points, Std.Error, tvalue),k, ncol=3, byrow = FALSE)
		rownames(coef.table) <-names(object$mass.point)
		colnames(coef.table) <- c("Estimate", "Std. Error", "t value")
		cat('\nCall: ',deparse(object$call),'\n\n', fill = TRUE)
		cat('Coefficients')
		cat(":\n")
		print(coef.table)
		p <- object$p
		cat('\nMixture proportions:\n')
		print.default(format(p,digits),print.gap=2,quote=FALSE)
		cat('\nMLE of sigma:\t  ',format(signif(object$sigma,digits)),'\n', fill = TRUE)
		cat('-2 log L:\t ',format(round(object$disparity,digits=1)),  'and AIC =  ', format(round(object$aic,digits)), fill = TRUE)
		if (!is.null(object$w)) cat('     Convergence at iteration ',round(object$EMiteration,0), fill = TRUE)
		cat('\n')
	}
	if (object$kind == '2'){
		cat("Minimal Disparity with EM converged:",object$MinDisparity, 
		    "at tol=", object$Mintol, "\n", fill = TRUE)
		cat('MinDisparity')
		cat(":\n")
		print.default(format(object$MinDisparity,digits),quote=FALSE) 
		cat('Mintol')
		cat(":\n")
		print.default(format(object$Mintol,digits),quote=FALSE) 
		cat('AllDisparities')
		cat(":\n")
		print.default(format(object$AllDisparities,digits),print.gap=2,quote=FALSE) 
		cat('Alltol')
		cat(":\n")
		print.default(format(object$Alltol,digits),print.gap=2,quote=FALSE) 
		cat('AllEMconverged')
		cat(":\n")
		print.default(format(object$AllEMconverged,digits),print.gap=2,quote=FALSE) 
	}
	if (object$kind == '3'){
		cat("Maximum profile log-likelihood:", object$objective, "at lambda=", object$Maximum, "\n", fill = TRUE)
		np <- length(object$beta)
		m <- seq(1,np)[substr(attr(object$beta,'names'),1,4)=='MASS']
		mass.points <- object$beta[m]
		Std.Error <-object$se
		tvalue <- object$beta/Std.Error
		MASS <- object$mass.point
		coef.table <- matrix(c(object$beta, Std.Error, tvalue), np, ncol=3, byrow = FALSE)
		colnames(coef.table) <- c("Estimate", "Std. Error", "t value")
		rownames(coef.table) <-names(object$beta)
		cat('\nCall: ',deparse(object$call),'\n\n', fill = TRUE)
		cat('Coefficients')
		cat(":\n")
		print(coef.table)
		cat('\nMass points:\n')
		print.default(format(MASS,digits),print.gap=2,quote=FALSE)
		p <- object$p
		cat('\nMixture proportions:\n')
		print.default(format(p,digits),print.gap=2,quote=FALSE)
		cat('\nMLE of sigma:\t  ',format(signif(object$sigma,digits)),'\n', fill = TRUE)
		cat('-2 log L:\t   ',format(round(object$disparity,digits=1)),  'and AIC =  ', format(round(object$aic,digits)), fill = TRUE)
		if (!is.null(object$w)) cat('     Convergence at iteration ',round(object$EMiteration,0), fill = TRUE)
		cat('\n')
	}
	invisible(object)
}

#' Plot diagnostics for boxcoxmix functions
#'
#' @description \code{plot()} is a generic function used to produce some useful diagnostic plotting of the functions:
#' \code{np.boxcoxmix()}, \code{optim.boxcox()} and \code{tolfind.boxcox()}.
#' @aliases plot.boxcoxmix
#' 
#' @rdname plot
#' @param x an object for which a plot is desired. 
#' @param plot.opt an integer value between 1 and 8.
#' @param \dots additional arguments. 
#' @return The plots to be printed depend on the number given in \code{plot.opt},
#' for the \code{np.boxcoxmix()}, \code{optim.boxcox()} and \code{tolfind.boxcox()} functions:
#' \item{1}{the disparities with the iteration number against the mass points} 
#' \item{2}{the fitted value against the response of the original and the transformed Data.}
#' \item{3}{probability plot of residuals of the original against the transformed data.}
#' \item{4}{individual posterior probabilities.}
#' \item{5}{control charts of residuals of the original against the transformed data.}
#' \item{6}{The histograms of residuals of the original against the transformed data.}
#' \item{7}{works only for the \code{tolfind.boxcox()} function and plots the specified range of \code{tol} against the disparities}
#' \item{8}{works only for the \code{optim.boxcox()} function and gives the profile likelihood function that summarises information concerning \code{lambda}.}
#' \item{9}{works only for the \code{Kfind.boxcox()} function and plots the specified range of \code{K} against the AIC or BIC information criteria}
#' \item{10}{works only for the \code{boxcoxtype()} function and gives the profile likelihood function that summarises information concerning \code{lambda} for generalized linear Mixed-effects Models.}
#' 
#' @import "stats"
#' @import "graphics"
#' @import "statmod"
#' @import "qicharts"
#' @method plot boxcoxmix
#' @export
plot.boxcoxmix <- function(x, plot.opt = 1, ...){
  
	K<-length(x$mass.point) 
	oldpar <- graphics::par(no.readonly = TRUE)
	on.exit(graphics::par(oldpar))
	if (plot.opt %in% c(7,8,9,10)){graphics::par(mfrow=c(1,1))}
	else if (plot.opt%in% c(1,4,5)){
		graphics::par(mfrow=c(2,1),cex=0.5,cex.axis=1.5,cex.lab=1.5)
	}
	else {
		graphics::par(mfrow=c(1,2),cex=0.5,cex.axis=1.5,cex.lab=1.5)
	}   


	if (x$kind %in% c(1,3) && K > 1 && plot.opt == 1){ #Disparities
		graphics::plot(1:x$EMiteration,x$Disparities,, col=1,type="l",xlab='EM iterations',ylab='-2logL')      
	}
	if (x$kind %in% c(1,3) && K > 1 && plot.opt==1){#  EM Trajectories
		ylim<- x$ylim; followmass<-x$masses     
		graphics::plot(1:x$EMiteration,followmass[,1],col=1,type='l',ylim=ylim,ylab='mass points',xlab='EM iterations')
		for (k in 1:K){ graphics::lines(1:x$EMiteration, followmass[,k],col=k)
		}
	}

	if(x$kind %in% c(1,3) && !identical(x$fitted, "none") && plot.opt == 2 ){# true values vs fitted
		class.col<-masspoint.class(x)
		graphics::plot(x$y,x$fitted, xlab="True response", ylab="fitted" ,col=class.col)
		graphics::abline(0,1) 
		graphics::plot(x$yt,x$fitted.transformed, xlab="Transformed response", ylab="fitted.transformed" ,col=class.col)
		graphics::abline(0,1)
	}
	if(x$kind %in% c(1,3) && plot.opt == 3){# QQplot residuals
		R<- x$residuals; RT<- x$residuals.transformed
		graphics::par(mfrow=c(1,2))
		stats::qqnorm(R, main="Q-Q plot of Original data")
		stats::qqline(R,col = "red", lwd = 3)
		stats::qqnorm(RT, main="Q-Q plot of Transformed data")
		stats::qqline(RT,col = "red", lwd = 3)
	}
	if(x$kind %in% c(1,3) && x$mform == "1" && K > 1 && plot.opt == 4){ 

		if (is.infinite(abs(max(x$residuals)))){return(c("Posterior probabilty plot not applicable for this object."))}
		#wik
		pmax<-max(x$w)
		graphics::plot(x$residuals, x$w[,1],col=1,type="p",ylim=c(0,pmax),ylab="post.prob", xlab="Residuals")
		for (k in 2:K) {
			graphics::points(x$residuals, x$w[,k],col=k,pch=k,type="p")
		}
		graphics::plot(x$residuals.transformed, x$w[,1],col=1,type="p",ylim=c(0,pmax),ylab="post.prob", xlab="Residuals transformed")
		for (k in 2:K) {
			graphics::points(x$residuals.transformed, x$w[,k],col=k,pch=k,type="p")
		}
	}

	if(x$kind %in% c(1,3) && plot.opt == 5){# control chart of residuals
		qicharts::qic(x$residuals,chart = 'i',cex=1.3,xlab=NULL , ylab=NULL, main="Control Chart of Original Data")
		qicharts::qic(x$residuals.transformed,chart = 'i',cex=1.3,xlab=NULL , ylab=NULL, main="Control Chart of Transformed Data")
	}

	if(x$kind %in% c(1,3) && plot.opt == 6){# hist residuals
		r<-x$residuals
		rt<- x$residuals.transformed
		graphics::hist(r, main="Original data",xlab=NULL , ylab=NULL,col="lightblue")
		graphics::hist(rt, main="Transformed data",xlab=NULL , ylab=NULL,col="pink")
	}
	if (x$kind=='2' && plot.opt == 7){#tolfind
		npc<-x$npcolors
		graphics::plot(x$Alltol, x$AllDisparities, type = "o", xlab = "tol", 
			       ylab = "Disparity", col= npc)
		graphics::segments(x$tol.min, x$min.Disp, x$tol.min, 1.1 * x$min.Disp, col = 4)
	}
	if (x$kind=='3' && plot.opt == 8){#optim
		plims <- graphics::par("usr")
		y0 <- plims[3]
		npc0<-x$npcolor
		result<- x$ylim1
		graphics::plot(x$All.lambda, x$profile.loglik, type = "l", xlab = expression(lambda), 
			       ylab = "Profile log-likelihood", ylim = result, col= npc0, main = x$maxl)
		graphics::segments(x$Maximum, y0, x$Maximum, x$objective, lty = 1, col = "red", lwd = 2)
	}
	if (x$kind=='4' && plot.opt == 9){#kfind
		plims <- graphics::par("usr")
		y0 <- plims[3]
		if(x$mselect=='1'){
			graphics::plot(x$All.K, x$Allaic, type = "o", xlab = "K", 
				       ylab = "AIC", col= "green", main=x$md)
			graphics::segments(x$Best.K, y0, x$Best.K, x$MinAIC, lty = 2, col = "red", lwd = 2)
		}else{
			graphics::plot(x$All.K, x$Allbic, type = "o", xlab = "K", 
				       ylab = "BIC", col= "green", main=x$md)
			graphics::segments(x$Best.K, y0, x$Best.K, x$MinBIC, lty = 2, col = "red", lwd = 2)
		}
	}
	if (x$kind=='5' && plot.opt == 10){#boxcoxtype
		graphics::plot(x$All.lam, x$profile.loglik, type = "l", xlab = expression(lambda), 
			       ylab = "Profile log-likelihood", ylim = x$ylim1, col= "green", main=x$ml)
		graphics::segments(x$Maximum, x$y0, x$Maximum, x$objective, lty = 1, col = "red", lwd = 3)
	}
	invisible(x)
}
#' @rdname utils
#' @param object ..
masspoint.class<-function(object){ 
	N<- length(object$Class)
	K<- length(object$mass.point)
	classif<-rep(0,N)
	for (i in 1:N){
		classif[i]<-order(object$w[i,])[K]
	} 
	classif
}  
