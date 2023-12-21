### R code from vignette source 'boxcoxmix.Rnw'

###################################################
### code chunk number 1: boxcoxmix.Rnw:333-338
###################################################
library(boxcoxmix)
data(strength, package="mdscore")
test.inv <- np.boxcoxmix(y ~ cut *lot, data = strength, K = 3,  
                         tol = 1.8,  start = "gq", lambda = -1, verbose=FALSE) 
test.inv


###################################################
### code chunk number 2: boxcoxmix.Rnw:344-347
###################################################
test.gauss <- np.boxcoxmix(y ~ cut *lot, data = strength, K = 3,  
                         tol = 1.8,  start = "gq", lambda = 1, verbose=FALSE) 
test.gauss


###################################################
### code chunk number 3: boxcoxmix.Rnw:354-358 (eval = FALSE)
###################################################
## test.optim <- optim.boxcox(y ~ cut*lot, data = strength,  K = 3,
##                            tol = 1.8, start = "gq", find.in.range = c(-3, 3),
##                            s = 60) 
## plot(test.optim, 8)


###################################################
### code chunk number 4: boxcoxmix.Rnw:375-378
###################################################
library(npmlreg)
inv.gauss <- alldist(y ~ cut*lot,  data = strength, k = 3,   tol = 0.45, verbose=FALSE,  family = "inverse.gaussian")
inv.gauss


###################################################
### code chunk number 5: boxcoxmix.Rnw:526-531
###################################################
data(Oxboys, package="nlme")
Oxboys$boy <- gl(26,9)
Oxboys$boy 
testox <- np.boxcoxmix(height ~  age, groups = Oxboys$boy,  data = Oxboys, 
                       K = 6, tol = 1, start = "gq", lambda=1, verbose=FALSE)


###################################################
### code chunk number 6: boxcoxmix.Rnw:533-534 (eval = FALSE)
###################################################
## plot(testox, 1)


###################################################
### code chunk number 7: boxcoxmix.Rnw:547-550 (eval = FALSE)
###################################################
## testo <- optim.boxcox(height ~ age, groups = Oxboys$boy, data = Oxboys, 
##          K = 6,  tol =1, start = "gq",  find.in.range=c( -1.2, 0.1),  s=15)
## plot(testo, 8)


