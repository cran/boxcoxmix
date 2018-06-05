# Box-Cox Transformation with Random Effects
#'
#' @rdname vc.em
#' @param p ..
#' @param beta ..
#' @param sigma .. 
#' @param z .. 
#' @export
np.estep <- function(y, x, lambda, p, beta, z, sigma){
  n <- length(y)
  K <- length(z)
  w <-d <- s <- matrix(0, n,K)
  f <- fik(y, x, lambda,  beta, z, sigma)
  logpf <- t(apply(f, 1, "+", log(p)))
  Mi <- apply(logpf, 1, max)
  Sik <- logpf - Mi
  ifelse(Sik < -760, 0, Sik)
  eSik <- exp(Sik)
  SeSik <- as.vector(apply(eSik, 1, sum))
  w <- eSik/SeSik
  return(w) }

##Zk
#'
#' @rdname vc.em
#' @export
np.zk <- function(y, x, w, beta, lambda){
  n <- dim(w)[1]
  K <- dim(w)[2]
  z <- rep(0,K)
  if  (all(!is.na(x))){  
      xbeta <- x%*%beta
  } else {xbeta <- 0}    
  for(k in 1:K){
  z[k] <- sum(w[,k]*(ytrans(y, lambda) - xbeta))/sum(w[,k])
  }
  return(z)
}
## fik
#'
#' @rdname vc.em
#' @export
fik <- function(y, x, lambda,  beta, z, sigma){
  n <- length(y)
  K <- length(z)
  f <- matrix(0,n,K)
  theta <- np.theta(y, x, lambda,  beta, z)
  for(i in 1:n){
   for(k in 1:K){
   ytrx <- (1/(2*sigma^2))*theta[i,k]
   d <- ((-1/2)*log(2*pi))
   j <- log(sigma)
   g <- (lambda -1)*log(y[i])
   f[i,k] <- d-j-ytrx + g
   }
  }
  return(f)
}
##theta
#' @rdname vc.em
#' @export
np.theta <- function(y, x, lambda,  beta, z){
  n <- length(y)
  K <- length(z)
 theta <- matrix(0, n,K) 
for(i in 1:n){
  if  (all(!is.na(x))){  
    xb <- x[i,]%*%beta
  } else { xb<- 0 }    
  for(k in 1:K){
    theta[i,k] <- (ytrans(y[i],lambda) - xb - z[k])^2
   }
 }
  return(theta)
}
##yhat
#'
#' @rdname vc.em
#' @param v ..
#' @export
yhat <-  function(v, lambda=1){
  if (is.list(v)) { r<- length(v)
  for (i in 1:r){
    if(lambda == 0){ y[[i]] <- exp(v)
    }
    else {
      y[[i]]<- (lambda*v +1)^(1/lambda) 
    }
  }
  } else {
    if(lambda == 0){y <- exp(v)
    }
    else {
      y <- (lambda*v + 1)^(1/lambda)
    }
  }
  return(y)
}
#ytrans
#'
#' @rdname vc.em
#' @export
ytrans <-  function(y, lambda=1){
  if (is.list(y)) { r<- length(y)
  v <- list()
  for (i in 1:r){
    if(lambda == 0){v[[i]]<-  log(y[[i]])
    }
    else {
      v[[i]] <-  (y[[i]]^lambda - 1)/ lambda
    }
  }
  } else {
    if(lambda == 0){v<-  log(y)
    }
    else {
      v <-  (y^lambda - 1)/ lambda
    }
  }
  return(v)
}
## bhat
#' @importFrom stats lm
#' @importFrom stats vcov
#' @rdname vc.em
#' @export
np.bhat <- function(y, x, w, z, lambda){
  n <- dim(w)[1]
  Ynew <- ytrans(y,lambda)- w%*%z
  out <- lm(Ynew ~ -1 + x)
  se<- sqrt(diag(vcov(out)))
  beta<-  out$coef
  return(list("beta"=beta,"se"=se))
}


## M-step
#'
#' @rdname vc.em
#' @export
np.mstep<- function(y, x, beta, lambda, w){
  n <- dim(w)[1]
  K <- dim(w)[2]
  p <- apply(w,2,sum)/n
  z <- rep(0,K)
  for (S in 1:40){
    z <- np.zk(y, x, w, beta, lambda)
    if  (all(!is.na(x))){    
        bse <- np.bhat(y, x, w, z, lambda)
        beta<-bse$beta
        se<- bse$se
    } else {
        beta <- 0
        se<-NA
    }       
  }
  var <- 0
  for(i in 1:n){
  theta <- np.theta(y, x, lambda,  beta, z)
  for(k in 1:K){
   var<- var+ w[i,k]*theta[i,k]
   }
  }  
  sigma <- sqrt(var/n)
  return(list("p"=p, "z"=z, "beta"=beta, "se"=se, "sigma"=sigma))
} 
## EMstep
#' @importFrom stats coef
#' @importFrom stats sd
#' @rdname vc.em
#' @export
np.em <- function(y, x, K, lambda=1, steps= 500,tol=0.5, start="gq", EMdev.change = 1e-04, plot.opt = 1,verbose = TRUE, ...){
  n <- length(y)
  if (!is.matrix(x)){
    x<- matrix(x,n,1) }
  p    <- rep(1/K,K)
  yt <- ytrans(y, lambda)
  if (all(!is.na(x))){
    a <- lm(formula = yt ~ -1 + x)
    se <- sqrt(diag(vcov(a)))
    beta <- coef(a)
  } else {
    beta<- 0
    se <- NA
  }
  sigma <- sd(yt)
  tol0 <- tol
  lambda0 <- lambda
  if (K > 1 && start == "quantile") {
    z <- mean(yt)+tol*stats::quantile(yt-mean(yt), probs= (1:K)/K-1/(2*K))
  }
  if (K > 1 && start == "gq") {
    if (K > 60) {
      K <- 60
      cat("K was set equal to 60, since the number of mass points supported are only up to 60 mass points. \n")
    }
    if  (all(!is.na(x))){
      b <- lm(formula = yt ~ x)
    } else {
      b <- lm(formula = yt ~ 1)
    }           
    Z <- b$coef[1] +tol*summary(b)$s*gqz(K, minweight = 1e-50)$location
    z<-Z[order(Z)]
  }
  if (K == 1){
    out <- lm(yt ~ x) 
    s <- 1
    z <- out$coef[1]
    sizes <-"none"
    se<- sqrt(diag(vcov(out)))
    beta <- coef(out)
    w <- matrix(1, n,1)
    masses <- "none"
    sigma <- summary(out)$sigma
    lik <- 0
    for(i in 1:n){
    lik <- lik + ((y[i]^(lambda -1))/(2*pi*sigma)^(n/2))*exp(((ytrans(y[i],lambda) - x[i]%*%beta )^2)/2*sigma^2)
    }
    logLik <- (-n/2)*log(2*pi)-n*log(sigma) - n/2 + (lambda -1)*sum(log(y))
    Disp <- -2*logLik
    complete.loglik <- "none"
    converged <- "none"
    Disparities<- "none"
  } else {
    s<-1
    previous_loglik <- -(2^1000)
    converged <- FALSE
    loglik <- Disparities <- 0
    masses<- matrix(0,0,K)
    while (s <= steps && !converged){
      if (verbose) {
        cat(s, "..")
      }
      w   <- np.estep(y, x, lambda, p, beta, z, sigma)
      fit <- np.mstep(y, x, beta, lambda, w)
      p  <- fit$p
      z  <- fit$z
      beta <- fit$beta
      sigma <-fit$sigma
      se <- fit$se
      lik <-matrix(0,n, K)                          
      f <- fik(y, x, lambda,  beta, z, sigma) 
      f <- exp(f)
      for(k in 1:K){
      lik[,k] <- p[k]*f[,k]
      }
      loglik[s] <- sum(log(apply(lik,1,sum)))
      Disparities[s] <- -2*loglik[s]
      converged <- abs(loglik[s] - previous_loglik)< EMdev.change
      previous_loglik <- loglik[s]
      masses <- rbind(masses,z)
      logLik <- loglik[s]
      Disp <- Disparities[s]
      s<-s+1
    }
    if (verbose) {
      cat("\n")
      if (converged) {
        cat("EM algorithm met convergence criteria at iteration # ", 
            s - 1, "\n")
      }
      else {
        cat("EM algorithm failed to meet convergence criteria at iteration # ", 
            s - 1, "\n")
      }
    }
    if (plot.opt==1){
      graphics::par(mfrow=c(2,1),cex=0.5,cex.axis=1.5,cex.lab=1.5)
    }
    if ((plot.opt == 2 ) && graphics::par("mfrow")[1] > 
        2) {
      plot.main <- substitute("tol" == tol0, list(tol0 = tol0))
    }
    else if ( (plot.opt == 3) && graphics::par("mfrow")[1] >
              2 ) {
      plot.main <- substitute("lambda"== lambda0, list(lambda0 = lambda0))
    }
    else {
      plot.main <- c("")
    }
    if (plot.opt == 1 || plot.opt == 2 || plot.opt == 3 ) {
      graphics::plot(1:(s - 1), Disparities, col = 1, type = "l", xlab = "EM iterations", 
                     ylab = "-2logL", main = plot.main)
      if (verbose) {
        cat("Disparity trend plotted.\n")
      }
    }
    complete.lik <-matrix(0,n, K)
    k <- 1:K
    complete.lik[,k] <- (lik[,k])^w[,k]
    complete.loglik <- sum(log(complete.lik))
  }
  np.fit<- list("EMiteration"=s-1, "masses"=masses, "kind"=1,"p"=p, "mass.point"=z, "beta"=beta, "se"=se, "sigma"=sigma, "w" =w, "loglik"= logLik, "complete.loglik" = complete.loglik, "Disparities"=Disparities, "disparity"= Disp, "likelihood"= lik, "EMconverged" = converged)
  class(np.fit)<-"boxcoxmix"
  return(np.fit)
}
#np.boxcox
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats model.response
#' @rdname vc.em
#' @export np.boxcox
np.boxcox<- function(formula, groups=1, data, K = 3, tol = 0.5,  lambda = 1, steps=500, EMdev.change = 1e-04,
                     plot.opt = 1, verbose = TRUE, start="gq", ...){
  call <- match.call()
  ddim <- dim(data)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mf <- model.frame(formula=formula, data=data)
  model <- "mixed"
  x <- model.matrix(attr(mf, "terms"), data = mf)[, -1, drop = FALSE]
  if (dim(x)[2]== 0) {
    x <- model.matrix(attr(mf, "terms"), data = mf)
    model <- "pure"
  }
  y <- model.response(mf, "numeric")
  if (any(y <= 0)) 
   stop("response variable must be positive")
  n <- NROW(y)
  xnames <- dimnames(x)[[2]]
  sizes <- groups
  mform <- strsplit(as.character(sizes), "\\|")
  mform <- gsub(" ", "", mform)
  if (length(mform) >= 2) {
    stop("Please use function vc.boxcox for two-level models")
  }
  npfit <- try(np.em(y=y, x=x, sizes= sizes,  K = K, lambda = lambda, steps = steps,tol = tol,
                      start = start, EMdev.change = EMdev.change, 
                     plot.opt = plot.opt, verbose = verbose))
  if (class(npfit) == "try-error") {
    cat("Check model specification, change the number of components or specify another value of lambda or tol and try again.")
    return()
  }
  s <- npfit$EMiteration
  w <- npfit$w
  p  <- npfit$p
  names(p) <- paste('MASS',1:K,sep='')
  z  <- npfit$mass.point
  names(z) <- paste('MASS',1:K,sep='')
  beta <- npfit$beta
  sigma <-npfit$sigma
  se <- npfit$se
  disp <- npfit$disparity
  Disparities<- npfit$Disparities
  aic<- disp+2*(length(beta)+2*K)
  bic<-disp+log(n)*(length(beta)+2*K)
  loglik <- npfit$loglik
  complete.loglik <- npfit$complete.loglik
  masses <- npfit$masses
  EMconverged <- npfit$EMconverged
  yt <- ytrans(y, lambda)
  if (K == 1){ylim<- "none"}
  if (model=='mixed'){
    if (K == 1){
      predicted.re <-rep(1,n)
      fitted.transformed <- rep(0,n)
      x <- model.matrix(attr(mf, "terms"), data=mf)
      for(i in 1:n){
      fitted.transformed[i] <- x[i,]%*%beta
      }
    } else{
      predicted.re <- rep(0,n)
      fitted.transformed <- rep(0,n)
      for(i in 1:n){
        predicted.re[i] <- sum(w[i,]%*%z)
        fitted.transformed[i]  <- x[i,]%*%beta + predicted.re[i]
       }
     }
    fitted <- yhat(fitted.transformed, lambda = lambda)
    residuals <- y- fitted
    residuals.transformed <- yt - fitted.transformed
    if (K > 1){
      ylim <- c(min(masses[, ]), max(masses[, ]))
      if (plot.opt == 1) {
        graphics::plot(1:s, masses[,1], col=1, type = "l", ylim=ylim, ylab='mass points',xlab='EM iterations')
        for(k in 2:K){
        graphics::lines(1:s, masses[,k], type="l", lwd= 1.5,lty=1, col=k)}
        if (verbose==T){cat("EM Trajectories plotted.\n")}
      }
    }
    npresult<- list("call"=call, "yt"=yt, "aic"=aic, "bic"=bic,"xx"=xnames, "Class"=y,"mform"=mform,"ylim"=ylim, "masses"=masses,
                    "y"=y, "formula"=formula,"data"=data,"EMiteration"= s, "Disparities"=Disparities,
                    "p"=p, "mass.point"=z, "beta"=beta, "se"=se, "sigma"=sigma,
                    "w" =w, "loglik"= loglik, "complete.loglik" = complete.loglik,
                    "disparity"= disp, "EMconverged" = EMconverged, "fitted" = fitted,
                    "fitted.transformed"= fitted.transformed, "predicted.re"= predicted.re,
                    "residuals"=residuals, "residuals.transformed"=residuals.transformed,
                    "kind"=1, "model"=model)
    class(npresult)<-"boxcoxmix"
  } else  {
    if (K > 1){
      ylim <- c(min(masses[, ]), max(masses[, ]))
      if (plot.opt == 1) {
        graphics::plot(1:s, masses[,1], col=1, type = "l", ylim=ylim, ylab='mass points',xlab='EM iterations')
        for(k in 2:K){
        graphics::lines(1:s, masses[,k], type="l", lwd= 1.5,lty=1, col=k)}
        if (verbose==T){cat("EM Trajectories plotted.\n")}
      }
    }
    fitted <- "none"
    fitted.transformed <- "none"
    predicted.re <- "none"
    npresult<- list("call"=call, "yt"=yt, "aic"=aic, "bic"=bic,"xx"=xnames, "Class"=y,"mform"=mform,"ylim"=ylim, "masses"=masses,
                    "y"=y, "formula"=formula,"data"=data,"EMiteration"= s, "Disparities"=Disparities,
                    "p"=p, "mass.point"=z, "beta"=beta, "se"=se, "sigma"=sigma,
                    "w" =w, "loglik"= loglik, "complete.loglik" = complete.loglik,
                    "disparity"= disp, "EMconverged" = EMconverged,"fitted" = fitted,
                    "fitted.transformed"= fitted.transformed, "predicted.re"= predicted.re,
                    "residuals"=y, "residuals.transformed"=yt,
                     "kind"=1, "model"= model)
    class(npresult)<-"boxcoxmixpure"
  }
  return(npresult)
}

#Box-Cox Transformation with Variance Component model
#'
#' @rdname vc.em
#' @export
vc.estep <- function(Y, X, sizes=1, lambda, p, beta, z, sigma){
  K <- length(z)
  r <- length(sizes)
  w <-d <- s <- matrix(0, r,K)
  m <- mik(Y, X, sizes, lambda, beta, z, sigma)
  logpm <- t(apply(m, 1, "+", log(p)))
  Mi <- apply(logpm, 1, max)
  Sik <- logpm - Mi
  ifelse(Sik < -760, 0, Sik)
  eSik <- exp(Sik)
  SeSik <- as.vector(apply(eSik, 1, sum))
  w <- eSik/SeSik
  return(w) }

##Zk
#'
#' @rdname vc.em
#' @export
zk <- function(Y, X, sizes, w, beta, lambda){
  r <-  dim(w)[1]
  K <- dim(w)[2]
  z <- rep(0,K)
  yx <- rep(0,r)
  for(i in 1:r){
  if  (all(!is.na(X))){  
      Xbeta <- X[[i]]%*%beta
  } else {Xbeta <- 0} 
  yx[i] <- sum(Y[[i]] - Xbeta)}
  for(k in 1:K){
    z[k] <- sum(w[,k]*yx)/sum(sizes*w[,k])
  }
  return(z)
}
## bhat
#'
#' @rdname vc.em
#' @param Y ..
#' @param X ..
#' @param w ..
#' @export
bhat <- function(Y, X, sizes, w, z, lambda){
  r <- dim(w)[1]
  wz  <- w%*%z
  WZ <- rep(wz, sizes)
  Xlong<-matrix(0, 0, dim(X[[1]])[2])
  for(i in 1:r){
  Xlong<-rbind(Xlong, X[[i]])}
  All<-as.data.frame(cbind(Ynew=unlist(Y)-WZ, Xlong))
  out <-  lm(Ynew ~   -1 + Xlong, data=All)
  se<- sqrt(diag(vcov(out)))
  beta<-  out$coef
  return(list("beta"=beta,"se"=se))
}

##mik
#'
#' @rdname vc.em
#' @export
mik <- function(Y, X, sizes, lambda,  beta, z, sigma){
  K <- length(z)
  r <- length(sizes)
  m <- matrix(0,r,K)
  theta <- vc.theta(Y, X, sizes, lambda,  beta, z)
  for(i in 1:r){ 
  for(k in 1:K){
      m[i,k] <- (-sizes[i]/2)*log(2*pi)-sizes[i]*log(sigma)-(1/(2*sigma^2))*theta[i,k] + (lambda -1)*sum(log( Y[[i]]))
   }
  }
  return(m) 
}
##theta
#' @rdname vc.em
#' @export
vc.theta <- function(Y, X, sizes, lambda,  beta, z){
  K <- length(z)
  r <- length(sizes)
  theta <- matrix(0, r,K)
  Ytrans<-ytrans(Y, lambda)
  for(i in 1:r){
    if  (all(!is.na(X))){  
      Xb <- X[[i]]%*%beta
    } else { Xb<- 0 }     
    for(k in 1:K){
      theta[i,k] <- sum((Ytrans[[i]] - Xb - z[k])^2)
    }
  }
  return(theta)
}
## M-step
#'
#' @rdname vc.em
#' @export
vc.mstep<- function(Y, X, sizes=1, beta, lambda, w){
  r <- dim(w)[1]
  K <- dim(w)[2]
  p <- apply(w,2,sum)/r
  z <- rep(0,K)
  Ytrans <- ytrans(Y, lambda)
  for (S in 1:40){
    z <- zk(Ytrans, X, sizes, w, beta, lambda)
    if  (all(!is.na(X))){
    bse<-  bhat(Ytrans, X, sizes, w, z, lambda)
    beta<-bse$beta
    se<- bse$se
    } else {
      beta <- 0
      se<-NA
    }  
  }
  var <- 0
  theta <- vc.theta(Y, X, sizes, lambda,  beta, z)
  for(i in 1:r){
   for(k in 1:K){
   var<- var+ w[i,k]*theta[i,k]
    }
  }
  sigma <- sqrt(var/sum(sizes))
  return(list("p"=p, "z"=z, "beta"=beta, "se"=se, "sigma"=sigma))
}
## EMstep
































#' Internal boxcoxmix functions
#' 
#' @title Internal boxcoxmix functions
#' @description auxiliary functions are not intended to be directly called from the user.
#' 
#' @rdname vc.em
#' @aliases np.boxcox vc.boxcox np.estep np.zk fik yhat ytrans np.bhat 
#' nb.se np.mstep np.em vc.estep bhat vc.se mik vc.mstep vc.em gqz
#' @param y ..
#' @param x ..
#' @param sizes ..
#' @author Amani Almohaimeed and Jochen Einbeck
#' @keywords em
#' @export 
vc.em <- function (y, x, sizes = 1, K, lambda, steps = 500, tol = 0.5, 
                   start = "gq", EMdev.change = 1e-04, plot.opt = 1, verbose = TRUE, 
                   ...) 
{
  n <- length(y)
  if (length(sizes) == 1) {
    if (sizes == 1) 
      r <- n
    else r <- 1
  }
  else {
    r <- length(sizes)
  }
  if (!is.matrix(x)) {
    x <- matrix(x, n, 1)
  }
  p <- rep(1/K, K)
  yt <- ytrans(y, lambda)
  if (all(!is.na(x))) {
    a <- lm(formula = yt ~ -1 + x)
    se <- sqrt(diag(vcov(a)))
    beta <- coef(a)
  }
  else {
    beta <- 0
    se <- NA
  }
  sigma <- sd(yt)
  tol0 <- tol
  lambda0 <- lambda
  if (K > 1 && start == "quantile") {
    z <- mean(yt) + tol * stats::quantile(yt - mean(yt), 
                                          probs = (1:K)/K - 1/(2 * K))
  }
  if (K > 1 && start == "gq") {
    if (K > 60) {
      K <- 60
      cat("K was set equal to 60, since the number of mass points supported are only up to 60 mass points. \n")
    }
    if (all(!is.na(x))) {
      b <- lm(formula = yt ~ x)
    }
    else {
      b <- lm(formula = yt ~ 1)
    }
    Z <- b$coef[1] + tol * summary(b)$s * gqz(K, minweight = 1e-50)$location
    z <- Z[order(Z)]
  }
  cumsize <- cumsum(sizes)
  Y <- Ytrans <- list()
  Y[[1]] <- y[1:sizes[1]]
  for(i in 2:r){
    Y[[i]] <- y[(cumsize[i - 1] + 1):cumsize[i]]}
  Ytrans <- ytrans(Y, lambda)
  X <- list()
  X[[1]] <- x[1:sizes[1], ]
  X[[1]] <- matrix(X[[1]], ncol = length(beta))
  for(i in 2:r){
    X[[i]] <- x[(cumsize[i - 1] + 1):cumsize[i], ]
    X[[i]] <- matrix(X[[i]], ncol = length(beta))
  }
  if (K == 1) {
    out <- lm(yt ~ x)
    s <- 1
    z <- coef(out)[1]
    sizes <- "none"
    se <- sqrt(diag(vcov(out)))
    beta <- coef(out)
    w <- matrix(1, n, 1)
    masses <- "none"
    sigma <- summary(out)$sigma
    lik <- 0
    for(i in 1:n){
    lik <- lik + ((y[i]^(lambda - 1))/(2 * pi * sigma)^(n/2)) * 
        exp(((ytrans(y[i], lambda) - x[i] %*% beta)^2)/2 * 
              sigma^2)}
    logLik <- (-n/2) * log(2 * pi) - n * log(sigma) - n/2 + 
      (lambda - 1) * sum(log(y))
    Disp <- -2 * logLik
    complete.loglik <- "none"
    converged <- "none"
    Disparities<- "none"
  }
  else {
    s <- 1
    previous_loglik <- -(2^1000)
    converged <- FALSE
    loglik <- Disparities <- 0
    masses <- matrix(0, 0, K)
    while (s <= steps && !converged) {
      if (verbose) {
        cat(s, "..")
      }
      w <- vc.estep(Y, X, sizes, lambda, p, beta, z, sigma)
      fit <- vc.mstep(Y, X, sizes, beta, lambda, w)
      p <- fit$p
      z <- fit$z
      beta <- fit$beta
      sigma <- fit$sigma
      se <- fit$se
      lik <- matrix(0, r, K)
      m <- mik(Y, X, sizes, lambda, beta, z, sigma)
      m <- exp(m)
      for(k in 1:K){
      lik[, k] <- p[k] * m[, k]}
      loglik[s] <- sum(log(apply(lik, 1, sum)))
      Disparities[s] <- -2 * loglik[s]
      converged <- abs(loglik[s] - previous_loglik) < EMdev.change
      previous_loglik <- loglik[s]
      masses <- rbind(masses, z)
      logLik <- loglik[s]
      Disp <- Disparities[s]
      s <- s + 1
    }
    if (verbose) {
      cat("\n")
      if (converged) {
        cat("EM algorithm met convergence criteria at iteration # ", 
            s - 1, "\n")
      }
      else {
        cat("EM algorithm failed to meet convergence criteria at iteration # ", 
            s - 1, "\n")
      }
    }
    if (plot.opt == 1) {
      graphics::par(mfrow = c(2, 1), cex = 0.5, cex.axis = 1.5, 
                    cex.lab = 1.5)
    }
    if ((plot.opt == 2) && graphics::par("mfrow")[1] > 2) {
      plot.main <- substitute("tol" == tol0, list(tol0 = tol0))
    }
    else if ((plot.opt == 3) && graphics::par("mfrow")[1] > 
             2) {
      plot.main <- substitute("lambda" == lambda0, list(lambda0 = lambda0))
    }
    else {
      plot.main <- c("")
    }
    if (plot.opt == 1 || plot.opt == 2 || plot.opt == 3) {
      graphics::plot(1:(s - 1), Disparities, col = 1, type = "l", 
                     xlab = "EM iterations", ylab = "-2logL", main = plot.main)
      if (verbose) {
        cat("Disparity trend plotted.\n")
      }
    }
    complete.lik <- matrix(0, r, K)
    k <- 1:K
    complete.lik[, k] <- (lik[, k])^w[, k]
    complete.loglik <- sum(log(complete.lik))
  }
  fit <- list(EMiteration = s - 1, masses = masses, kind = 1, Disparities = Disparities,
              p = p, mass.point = z, beta = beta, se = se,  sigma = sigma, 
              w = w, loglik = logLik, complete.loglik = complete.loglik, 
              disparity = Disp, likelihood = lik, EMconverged = converged)
  class(fit) <- "boxcoxmix"
  return(fit)
}
# boxcoxmix






































#' Response Transformations for Random Effect and Variance Component Models
#' 
#' 
#' 
#' @rdname np.boxcoxmix
#' 
#' @description The function \code{np.boxcoxmix()} fits an overdispersed generalized linear model
#' and variance component models using nonparametric profile maximum likelihood.
#' 
#' 
#' @details The Box-Cox transformation (Box & Cox, 1964) is applied to the overdispersed
#' generalized linear models and variance component models with an unspecified
#' mixing distribution. The NPML estimate of the mixing distribution is known
#' to be a discrete distribution involving a finite number of mass-points and corresponding
#' masses (Aitkin et al., 2009). An Expectation-Maximization (EM) algorithm is
#' used for fitting the finite mixture distribution, one needs to specify the
#' number of components \code{K} of the finite mixture in advance. To stop the EM-algorithm
#' when it reached its convergence point, we need to defined the convergence criteria that is 
#' the absolute change in the successive log-likelihood function values being less than an arbitrary 
#' parameter such as \code{EMdev.change} = 0.0001 (Einbeck et at., 2014). This algorithm can be 
#' implemented using the function \code{np.boxcoxmix()}, which is designed to account for 
#' overdispersed generalized linear models and variance component models using the non-parametric
#' profile maximum likelihood (NPPML) estimation.
#' 
#' 
#' The ability of the EM algorithm to locate the global maximum in fewer iterations
#' can be affected by the choice of initial values, the function \code{np.boxcoxmix()}
#' allows us to choose from two different methods to set the initial value of the mass
#' points. When option "gq" is set, then Gauss-Hermite masses and mass points are used
#' as starting points in the EM algorithm, while setting start= "quantile" uses the 
#' Quantile-based version to select the starting points. 
#' 
#' @aliases np.boxcoxmix
#' @param formula a formula describing the transformed response and the fixed
#' effect model (e.g. y ~ x).
#' @param groups the random effects. To fit overdispersion models , set \code{groups} = 1.
#' @param data a data frame containing variables used in the fixed and random
#' effect models.
#' @param K the number of mass points.
#' @param tol a positive scalar (usually, 0< \code{tol} <= 2)
#' @param lambda a transformation parameter, setting \code{lambda}=1 means 'no
#' transformation'.
#' @param steps maximum number of iterations for the EM algorithm.
#' @param EMdev.change a small scalar, with default 0.0001, used to determine
#' when to stop EM algorithm.
#' @param plot.opt Set \code{plot.opt=1}, to plot the disparity against
#' iteration number. Use \code{plot.opt=2} for \code{tolfind.boxcox()} and \code{plot.opt=3}
#' for \code{optim.boxcox()}.
#' @param verbose If set to FALSE, no printed output on progress.
#' @param start a description of the initial values to be used in the fitted
#' model, Quantile-based version "quantile" or Gaussian Quadrature "gq" can be
#' set.
#' @param \dots extra arguments will be ignored.
#' @return 
#' \item{mass.point}{the fitted mass points.} \item{p}{the masses corresponding to the mixing proportions.} \item{beta}{the
#' vector of coefficients.} \item{sigma}{the standard deviation of the mixing distribution (the square root of the variance).} \item{se}{the standard error of the estimate.}
#' \item{w}{a matrix of posterior probabilities that element i comes from cluster k.}
#' \item{loglik}{the log-likelihood of the fitted regression model.}
#' \item{complete.loglik}{the complete log-likelihood of the fitted regression model.}
#' \item{disparity}{the disparity of the fitted regression model.} \item{EMiteration}{provides the number of iterations of the EM algorithm.} 
#' \item{EMconverged}{TRUE means the EM algorithm converged.} \item{call}{the matched call.} \item{formula}{the formula provided.}
#' \item{data}{the data argument.} \item{aic}{the Akaike information criterion of the fitted regression model.}
#'  \item{bic}{the Bayesian information criterion of the fitted regression model.}
#' \item{fitted}{the fitted values for the individual observations.} \item{fitted.transformed}{the fitted values for
#' the individual transformed observations.} \item{residuals}{the difference between the observed values and the fitted values.}
#' \item{residuals.transformed}{the difference between the transformed observed values and the transformed fitted values.}
#' \item{predicted.re}{a vector of predicted residuals.}
#' 
#' The other outcomes are not relevant to users and they are intended for internal use only.
#' 
#' @author Amani Almohaimeed and Jochen Einbeck
#' @seealso \code{\link{optim.boxcox}}, \code{\link{tolfind.boxcox}}.
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
#' @keywords boxcox random variance
#' @examples
#' #Pennsylvanian Hospital Stay Data
#' data(hosp, package = "npmlreg")
#' test1 <- np.boxcoxmix(duration ~ age + wbc1, data = hosp, K = 2, tol = 1, 
#'     start = "quantile", lambda = 1)
#' round(summary(test1)$w, digits = 3)
#' # [1,] 1.000 0.000
#' 
#' # Refinery yield of gasoline Data
#' data(Gasoline, package = "nlme")
#' test2.vc <- np.boxcoxmix(yield ~ endpoint + vapor, groups = Gasoline$Sample, 
#'       data = Gasoline, K = 3, tol = 1.7, start = "quantile", lambda = 0)
#' test2.vc$disparity
#' # [1] 176.9827
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
# vc.boxcox 
#' @export np.boxcoxmix
np.boxcoxmix<-function(formula, groups = 1, data, K = 3, tol = 0.5, lambda = 1, 
                    steps = 500, EMdev.change = 1e-04, plot.opt = 1, verbose = TRUE, 
                    start="gq", ...){
  call <- match.call()
  mform <- strsplit(as.character(groups), "\\|")
  mform <- gsub(" ", "", mform)
  if (length(mform) == 1) {
    fit <- try(np.boxcox(formula=formula, groups= groups, data=data, K=K, lambda=lambda,steps= steps, 
                         tol = tol, start=start, EMdev.change = EMdev.change, plot.opt = plot.opt, verbose = verbose))
    if (class(fit) == "try-error") {
      cat("Check model specification, change the number of components or specify another value of lambda or tol and try again.")
      return()
    }
  }
  else {
    fit <- try(vc.boxcox(formula=formula, groups= groups, data=data, K=K, lambda=lambda, steps= steps,
                         tol = tol, start=start, EMdev.change = EMdev.change, plot.opt = plot.opt, verbose = verbose))
    if (class(fit) == "try-error") {
      cat("Check model specification, change the number of components or specify another value of lambda or tol and try again.")
      return()
    }
  }
  W <- fit$w
  P <-  fit$p
  se <- fit$se
  iter <- fit$EMiteration
  names(P) <- paste('MASS',1:K,sep='')
  Z <- fit$mass.point
  names(Z) <- paste('MASS',1:K,sep='')
  Beta <- fit$beta
  Sigma<- fit$sigma
  aic<- fit$aic
  bic<- fit$bic
  y <- fit$y
  yt <- fit$yt
  fitted <- fit$fitted
  fitted.transformed <- fit$fitted.transformed
  masses<- fit$masses
  ylim<- fit$ylim
  residuals<- fit$residuals
  residuals.transformed<- fit$residuals.transformed
  predicted.re <- fit$predicted.re
  Class<-fit$Class
  xx <- fit$xx
  Disp <- fit$disparity
  Disparities <- fit$Disparities
  Loglik <- fit$loglik
  complete.loglik <- fit$complete.loglik
  kind <- fit$kind
  EMconverged <- fit$EMconverged
  model <- fit$model
  if (model == "mixed") {
    result<- list("call"=call,  "p"=P, "mass.point"=Z, "beta"=Beta, "sigma"=Sigma, "se"=se, "w" =W, "Disparities" = Disparities,
                  "formula" = formula, "data" = data, "loglik"= Loglik, "aic"=aic, "bic"=bic,"masses"=masses, "y"=y, "yt"=yt,
                  "complete.loglik" = complete.loglik, 
                  "disparity"= Disp, "EMconverged" = EMconverged, "mform"=length(mform),"ylim"=ylim, 
                  "fitted" = fitted, "Class"= Class, "fitted.transformed"= fitted.transformed, "predicted.re"= predicted.re,
                  "residuals"=residuals, "residuals.transformed"=residuals.transformed,"kind"=kind,
                  "EMiteration"= iter, "xx" = xx)
    class(result)<-"boxcoxmix"
  }else{
    result<- list("call"=call,  "p"=P, "mass.point"=Z, "beta"=Beta, "sigma"=Sigma, "se"=se, "w" =W, "Disparities" = Disparities,
                  "formula" = formula, "data" = data, "loglik"= Loglik, "aic"=aic, "bic"=bic, "masses"=masses, "y"=y, "yt"=yt,
                  "complete.loglik" = complete.loglik, 
                  "disparity"= Disp, "EMconverged" = EMconverged, "mform"=length(mform),"ylim"=ylim, 
                  "fitted" = fitted, "Class"= Class, "fitted.transformed"= fitted.transformed, "predicted.re"= predicted.re,
                  "residuals"=residuals, "residuals.transformed"=residuals.transformed,"kind"=kind,
                  "EMiteration"= iter, "xx" = xx)
    class(result)<-"boxcoxmixpure"
  }
  return(result)
}
#' @rdname vc.em
#' @param formula a formula describing the transformed response and the fixed
#' effect model (e.g. y ~ x).
#' @param groups the random effects. To fit overdispersion models , set \code{groups} = 1.
#' @param data a data frame containing variables used in the fixed and random
#' effect models.
#' @param K the number of mass points.
#' @param tol a positive scalar (usually, 0< \code{tol} <= 2)
#' @param lambda a transformation parameter, setting \code{lambda}=1 means 'no
#' transformation'.
#' @param steps maximum number of iterations for the EM algorithm.
#' @param EMdev.change a small scalar, with default 0.0001, used to determine
#' when to stop EM algorithm.
#' @param plot.opt Set plot.opt=1, to plot the disparity against
#' iteration number. Use \code{plot.opt=2} for \code{tolfind.boxcox} and \code{plot.opt=3}
#' for \code{optim.boxcox}.
#' @param verbose If set to FALSE, no printed output on progress.
#' @param start a description of the initial values to be used in the fitted
#' model, Quantile-based version "quantile" or Gaussian Quadrature "gq" can be
#' set.
#' @param \dots extra arguments will be ignored. 
#' @export  
vc.boxcox<- function (formula, groups = 1, data, K = 3, tol = 0.5, lambda = 1, 
                      steps = 500, EMdev.change = 1e-04, plot.opt = 1, verbose = TRUE, 
                      start="gq", ...) 
{
  call <- match.call()
  data <- as.data.frame(data)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  datay <- data
  groupy <- groups
  ordered <- data[with(data, order(groups)), ]
  mf <- model.frame(formula = formula, data = ordered)
  model <- "mixed"
  x <- model.matrix(attr(mf, "terms"), data = mf)[, -1, drop = FALSE]
  if (dim(x)[2]== 0) {
    x <- model.matrix(attr(mf, "terms"), data = mf)
    model <- "pure"
  }
  y <- model.response(mf)
  if (any(y <= 0)) 
    stop("response variable must be positive")
  N <- NROW(y)
  data <- if (is.matrix(y)) 
    data[dimnames(y)[[1]], ]
  else data[names(y), ]
  xnames <- dimnames(x)[[2]]
  groups <- groups[match(rownames(ordered), rownames(datay))]
  sizes <- table(groups)
  mform <- strsplit(as.character(sizes), "\\|")
  mform <- gsub(" ", "", mform)
  if (length(mform) == 1) {
    stop("Please use function np.boxcox for overdispersion models")
  }
  vfit <- try(vc.em(y = y, x = x, sizes = sizes, K = K, lambda = lambda, 
                    steps = steps, tol = tol, start = start, EMdev.change = EMdev.change, 
                    plot.opt = plot.opt, verbose = verbose))
   if (class(vfit) == "try-error") {
    cat("Check model specification, change the number of components or specify another value of lambda or tol and try again.")
    return()
   }
  EMiter <- vfit$EMiteration
  if (K == 1) {
    w <- vfit$w
  }
  else {
    w <- vfit$w
    w <- w[match(unique(groupy), unique(groups)), ]
    rownames(w) <- unique(groupy)
  }
  p <- vfit$p
  names(p) <- paste("MASS", 1:K, sep = "")
  z <- vfit$mass.point
  names(z) <- paste("MASS", 1:K, sep = "")
  beta <- vfit$beta
  sigma <- vfit$sigma
  se <- vfit$se
  disp <- vfit$disparity
  Disparities <- vfit$Disparities
  loglik <- vfit$loglik
  EMconverged <- vfit$EMconverged
  complete.loglik <- vfit$complete.loglik
  masses <- vfit$masses
  aic <- disp + 2 * (length(beta) + 2 * K)
  bic<-disp+log(N)*(length(beta)+2*K)
  yt <- ytrans(y, lambda)
  if (K == 1) {
    ylim <- "none"
  }
  if (model == "mixed") {
    if (K == 1) {
      predicted.re <- rep(1, N)
    }
    else {
      r <- dim(w)[1]
      predicted <- rep(0, r)
      predicted[1] <- sum(w[1, ] %*% z)
      i <- 2:r
      predicted[i] <- sum(w[i, ] %*% z)
      names(predicted) <- unique(groupy)
      predicted.re <- predicted
    }
    if (K == 1) {
      fitted.transformed <- rep(0, N)
      x <- model.matrix(attr(mf, "terms"), data = mf)
      for(i in 1:N){
        fitted.transformed[i] <- x[i, ] %*% beta
      }
    }
    else {
      predicted <- predicted[match(unique(groups), unique(groupy))]
      fitted.transformed <- matrix(0, N, 1)
      cumsize <- cumsum(sizes)
      for(j in 1:cumsize[1]){
        fitted.transformed[j, ] <- x[j, ] %*% beta + predicted[1]
      }
      for(i in 2:r){
      for(j in (cumsize[i - 1] + 1):cumsize[i]){
      fitted.transformed[j, ] <- x[j, ] %*% beta + predicted[i]
       }
      }
      rownames(fitted.transformed) <- rownames(ordered)
      fitted.transformed <- fitted.transformed[match(row.names(datay), 
                                                     row.names(ordered)), 1, drop = FALSE]
    }
    fitted <- yhat(fitted.transformed, lambda = lambda)
    residuals <- y - fitted
    residuals.transformed <- yt - fitted.transformed
    if (K > 1) {
      ylim <- c(min(masses[, ]), max(masses[, ]))
      if (plot.opt == 1) {
        graphics::plot(1:EMiter, masses[, 1], col = 1, type = "l", 
                       ylim = ylim, ylab = "mass points", xlab = "EM iterations")
        for(k in 2:K){
        graphics::lines(1:EMiter, masses[, k], type = "l", 
                          lwd = 1.5, lty = 1, col = k)}
        if (verbose == T) {
          cat("EM Trajectories plotted.\n")
        }
      }
    }
    result <- list(call = call, aic = aic, bic = bic, xx = xnames, yt = yt, Disparities = Disparities,
                   Class = sizes, mform = mform, ylim = ylim, masses = masses, 
                   y = y, formula = formula, data = data, EMiteration = EMiter, 
                   p = p, mass.point = z, beta = beta, se = se, sigma = sigma, 
                   w = w, loglik = loglik, complete.loglik = complete.loglik, 
                   disparity = disp, EMconverged = EMconverged, fitted = fitted, 
                   fitted.transformed = fitted.transformed, predicted.re = predicted.re, 
                   residuals = residuals, residuals.transformed = residuals.transformed, 
                   kind = 1,model=model)
    class(result) <- "boxcoxmix"
  }
  else {
    if (K > 1) {
      ylim <- c(min(masses[, ]), max(masses[, ]))
      if (plot.opt == 1) {
        graphics::plot(1:EMiter, masses[, 1], col = 1, type = "l", 
                       ylim = ylim, ylab = "mass points", xlab = "EM iterations")
        for(k in 2:K){
        graphics::lines(1:EMiter, masses[, k], type = "l", 
                          lwd = 1.5, lty = 1, col = k)}
        if (verbose == T) {
          cat("EM Trajectories plotted.\n")
        }
      }
    }
    fitted <- "none"
    fitted.transformed <- "none"
    predicted.re <- "none"
    result <- list(call = call, aic = aic, bic = bic, xx = xnames, yt = yt, Disparities = Disparities,
                   Class = sizes, mform = mform, ylim = ylim, masses = masses, 
                   y = y, formula = formula, data = data, EMiteration = EMiter, 
                   p = p, mass.point = z, beta = beta, se = se, sigma = sigma, 
                   w = w, loglik = loglik, complete.loglik = complete.loglik, fitted = fitted, 
                   fitted.transformed = fitted.transformed, predicted.re = predicted.re,
                   disparity = disp, EMconverged = EMconverged, residuals = y, 
                   residuals.transformed = yt, kind = 1, model= model)
    class(result) <- "boxcoxmixpure"
  }
  return(result)
}
##gqz
#' @import "statmod"
#' @rdname vc.em
#' @param numnodes ..
#' @param minweight ..
#' @export
gqz <- function (numnodes = 20, minweight = 1e-06) 
{
  out <- statmod::gauss.quad(numnodes, "hermite")
  h <- rbind(out$nodes * sqrt(2), out$weights/sum(out$weights))
  ord <- order(h[1, ], decreasing = TRUE)
  h <- h[, ord]
  h <- cbind(h[1, ], h[2, ])
  h <- subset(as.data.frame(h), (h[, 2] >= minweight))
  names(h) <- c("location", "weight")
  h
}