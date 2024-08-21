### R code from vignette source 'boxcoxmix.Rnw'

###################################################
### code chunk number 1: boxcoxmix.Rnw:57-58
###################################################
options(width=60)


###################################################
### code chunk number 2: boxcoxmix.Rnw:339-345
###################################################
library(boxcoxmix)
data(strength, package="mdscore")
test.inv <- np.boxcoxmix(y ~ cut *lot, data = strength, K = 3,  
                         tol = 1.8,  start = "gq", lambda = -1,
                         verbose=FALSE) 
test.inv


###################################################
### code chunk number 3: boxcoxmix.Rnw:351-355
###################################################
test.gauss <- np.boxcoxmix(y ~ cut *lot, data = strength, K = 3,  
                           tol = 1.8,  start = "gq", lambda = 1,
                           verbose=FALSE) 
test.gauss


###################################################
### code chunk number 4: boxcoxmix.Rnw:362-366 (eval = FALSE)
###################################################
## test.optim <- optim.boxcox(y ~ cut*lot, data = strength,  K = 3,
##                            tol = 1.8, start = "gq", find.in.range = c(-3, 3),
##                            s = 60) 
## plot(test.optim, 8)


###################################################
### code chunk number 5: boxcoxmix.Rnw:383-387
###################################################
library(npmlreg)
inv.gauss <- alldist(y ~ cut*lot, data = strength, k = 3, tol = 0.45,
                     verbose=FALSE,  family = "inverse.gaussian")
inv.gauss


###################################################
### code chunk number 6: boxcoxmix.Rnw:535-541
###################################################
data(Oxboys, package="nlme")
Oxboys$boy <- gl(26,9)
Oxboys$boy 
testox <- np.boxcoxmix(height ~  age, groups = Oxboys$boy,
                       data = Oxboys, K = 6, tol = 1, start = "gq",
                       lambda=1, verbose=FALSE)


###################################################
### code chunk number 7: boxcoxmix.Rnw:543-544 (eval = FALSE)
###################################################
## plot(testox, 1)


###################################################
### code chunk number 8: boxcoxmix.Rnw:557-560 (eval = FALSE)
###################################################
## testo <- optim.boxcox(height ~ age, groups = Oxboys$boy, data = Oxboys, 
##          K = 6, tol =1, start = "gq", find.in.range = c( -1.2, 0.1), s=15)
## plot(testo, 8)


