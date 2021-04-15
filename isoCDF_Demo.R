########################################################
### Isotonic distributional regression (total order) ###
### January 22, 2021                                 ###
########################################################

# This file illustrates the methods described in
# "Accelerating the pool-adjacent-violators algorithm
#  for isotonic distributional regression",
#  A. Henzi, A. Moesching, L. Duembgen (2021).
                     

source("isoCDF.R")

n <- 70
z <- ((n+1)/2 - (1:n))/n + rnorm(n)

PAVA_demo(z)

z2 <- round(z)
PAVA_demo(z2)
PAVA2_demo(z2)


#-------------- A first simulation ----------------------------------------#
#
# X ~ Unif(0, 10)
# Y ~ Gamma(shape = sqrt(X), scale = min(max(X, 1), 6))
# Sample size: n = 3000
# Weights: W ~ Unif(0, 1)

set.seed(231120)

n <- 1000
# X <- round(runif(n, 0, 10),2)
# Y <- round(rgamma(n, shape = sqrt(X), scale = pmin(pmax(X, 1), 6)),2)
X <- sort(runif(n, 0, 10))
Y <- rgamma(n,
	shape = sqrt(X),
	scale = 2 + (X - 5)/sqrt(2 + (X - 5)^2))

# # Later on:
# X <- c(0,X,10)
# Y <- c(max(Y),Y,min(Y))


# W <- runif(n)
par(mfrow=c(1,1))

plot(X, Y)

t1 <- Sys.time()
idr1 <- isoCdfNaive(X = X, Y = Y)
(t1 <- Sys.time() - t1)

t2 <- Sys.time()
idr2 <- isoCdfNaive2(X = X, Y = Y)
(t2 <- Sys.time() - t2)

t3 <- Sys.time()
idr3 <- isoCdf(X = X, Y = Y)
(t3 <- Sys.time() - t3)


round(as.numeric(t1)/as.numeric(t2),5)
round(as.numeric(t2)/as.numeric(t3),5)
round(as.numeric(t1)/as.numeric(t3),5)


# all.equal(idr1$cdf, idr2$cdf, tolerance = 1e-14)
# all.equal(idr2$cdf, idr3$cdf, tolerance = 1e-14)

# idr1 <- isoCdf(X,Y,beta=c(0.05,0.25,0.5,0.75,0.95))
idr2 <- isoQuant(X,Y,beta=c(0.05,0.25,0.5,0.75,0.95))
# all.equal(idr1$Q,idr2$Q,tolerance=1e-14)

x <- idr2$x
Qmin <- idr2$Qmin
Qmax <- idr2$Qmax
Q <- idr2$Q

par(cex=1.3,
	mai=c(0.5,0.5,0.02,0.02),
	mgp=c(0.5,0.6,0))
plot(X,Y,col='darkgrey',xlab='',ylab='')
for (j in 1:dim(Q)[2]){
  lines(x,Q[,j],lwd=2,col='blue')
}

# idr2 <- isoQuant(X,Y,beta=0.5)
# idr2$Q


# Illustrate isoCdf step by step:

set.seed(231120)

n <- 30
# X <- round(runif(n, 0, 10),2)
# Y <- round(rgamma(n, shape = sqrt(X), scale = pmin(pmax(X, 1), 6)),2)
X <- 1:n
Y <- X + rnorm(n,0,7)

plot(X,Y)

isoCdf_demo(Y)


#### Monte-Carlo experiment for computation times:

set.seed(231120)
n <- 1000

mcsim <- 500

tv1 <- rep(NA,mcsim)
tv2 <- rep(NA,mcsim)
tv3 <- rep(NA,mcsim)

for (s in 1:mcsim){
	print(paste('Simulation',s),quote=FALSE)
	X <- sort(runif(n, 0, 10))
	Y <- rgamma(n,
		shape = sqrt(X),
		scale = 2 + (X - 5)/sqrt(2 + (X - 5)^2))
	t1 <- Sys.time()
	idr1 <- isoCdfNaive(X = X, Y = Y)
	(t1 <- Sys.time() - t1)
	tv1[s] <- as.numeric(t1)
	t2 <- Sys.time()
	idr2 <- isoCdfNaive2(X = X, Y = Y)
	(t2 <- Sys.time() - t2)
	tv2[s] <- as.numeric(t2)
	t3 <- Sys.time()
	idr3 <- isoCdf(X = X, Y = Y)
	(t3 <- Sys.time() - t3)
	tv3[s] <- as.numeric(t3)
}


boxplot(tv1/tv2, tv2/tv3, tv1/tv3,
	names=c('stand./modif.','modif./abr.','stand./abr.'),
	lty=1,lwd=2,
	ylim=c(1,max(tv1/tv3)))
abline(h=1,col='blue')

boxplot(log10(tv1/tv2), log10(tv2/tv3), log10(tv1/tv3),
	names=c('stand./modif.','modif./abr.','stand./abr.'),
	lty=1,lwd=2,
	ylim=c(0,max(log10(tv1/tv3))))
abline(h=0,col='blue')

c(mean(tv1),sd(tv1))
c(mean(tv2),sd(tv2))
c(mean(tv3),sd(tv3))

c(mean(tv1/tv1),sd(tv1/tv2))
c(mean(tv1/tv3),sd(tv1/tv3))
c(mean(tv2/tv3),sd(tv2/tv3))
