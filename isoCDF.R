########################################################
### Isotonic distributional regression (total order) ###
### January 22, 2021                                 ###
########################################################

# This file contains implementations of the methods 
# described in
# "Accelerating the pool-adjacent-violators algorithm
#  for isotonic distributional regression",
#  A. Henzi, A. Moesching, L. Duembgen (2021).

PAVA <- function(z, w)
# Pool adjacent violators
#
# Computation of the antitonic least squares regression
# of a vector z with weights w. That means,
# the program computes an antitonic vector a (i.e.
# with non-increasing components) such that
#    sum(w*(a - z)^2)
# is minimial.
# 
# Input:
# - z:  Vector z to be fitted.
# - w:  Weight vector w with strictly positive components.
#
# Output:
# - d:  The size of the final partition.
# - PP: A vector containing the right-most indices bi, 0 <= i <= d,
#       of the partition P = {(a0,...,b0), (a1,...,b1), ..., (ad,...,bd)}
#       of {0,1,...,length(z)} into index intervals, where a0 = b0 = 0
#       and b0 + 1 = a1, b1 + 1 = a2, ... That means,
#       PP[1] = 0, PP[2] = b1, ..., PP[d+1] = bd = length(z).
# - WW: A vector containing the weights corresponding to each element of
#       the partition PP, where WW[1] = 0.
# - MM: A vector containing the value of the antitonic approximation A(z)
#       on each interval of the partition PP, where MM[1] = Inf.
# - CA: The number of times an interval was added to the partition.
# - CM: The number of times two adjacent intervals of the partition
#       have been merged.
#
# The vectors PP, WW, MM are predefined to have size m+1 for computation
# purposes. Only the coefficients 2 to d+1 of PP, WW and MM are relevant.
{
	m <- length(z)
	PP <- c(0, 1, integer(m - 1))
	WW <- c(0, w[1], numeric(m - 1))
	MM <- c(Inf, z[1], numeric(m - 1))
	
	if (m == 1){
		return(list(d = 1, PP = PP, WW = WW,
			MM = MM, CA = 1, CM = 0))
	}
	
	d <- 2
	# During the algorithm, d is the number of relevant
	# intervals in the partition plus 1
	for (i in 2:m) {
		d <- d + 1
		PP[d] <- i
		WW[d] <- w[i]
		MM[d] <- z[i]
		
		while (MM[d - 1] <= MM[d]) {
			dd <- c(d - 1, d)
			d <- d - 1
			MM[d] <- sum(WW[dd] * MM[dd])
			WW[d] <- sum(WW[dd])
			MM[d] <- MM[d]/WW[d]
			PP[d] <- PP[d + 1]
		}
	}
	d <- d - 1
	# In the end we return the number of relevant
	# intervals in the partition.
	return(list(d = d, PP = PP, WW = WW,
		MM = MM, CA = m, CM = m-d))
}

PAVA2 <- function(z, w)
# Same as PAVA, but checking z for constant blocks.
{
	m <- length(z)
	PP <- c(0, rep(0,m))
	WW <- c(0, rep(0,m))
	MM <- c(Inf, rep(0,m))
	
	if (m == 1){
		return(list(d = 1, PP = c(0,1), WW = c(0,w[1]),
			MM = c(Inf,z[1]), CA = 1, CM = 0))
	}
	
	CA <- 0
	CM <- 0
	d <- 1
	# During the algorithm, d is the number of relevant
	# intervals in the partition plus 1.
	while (PP[d] < m) {
		pp <- PP[d] + 1
		ww <- w[pp]
		d <- d+1
		MM[d] <- z[pp]
		while (pp < m && z[pp+1] == MM[d]){
			pp <- pp + 1
			ww <- ww + w[pp]
		}
		CA <- CA + 1
		PP[d] <- pp
		WW[d] <- ww
		
		while (MM[d - 1] <= MM[d]) {
			dd <- c(d - 1, d)
			d <- d - 1
			MM[d] <- sum(WW[dd] * MM[dd])
			WW[d] <- sum(WW[dd])
			MM[d] <- MM[d]/WW[d]
			PP[d] <- PP[d + 1]
			CM <- CM + 1
		}
	}
	d <- d - 1
	# In the end we return the number d of relevant
	# intervals in the partition
	
	return(list(d = d, PP = PP, WW = WW,
		MM = MM, CA = CA, CM = CM))
}

prepareData <- function(X, Y, W=rep(1,length(X)))
# Convert raw data into input for functions isoCdf, isoCdfNaive, isoCdf2.
# 
# Input:
# - X:       covariate vector (numeric vector)
# - Y:       response vector (numeric vector)
# - W:       weight vector (positive numeric vector)
#
# The inputs X, Y and W must all have the same length.
#
# Output:
# - x:       sorted distinct values of X (numeric vector)
# - y:       sorted distinct values of Y (numeric vector)
# - w:       aggregated weights for all entries of x (numeric vector)
# - W        weight vector W sorted according to order of Y (numeric vector)
# - Y        sorted vector Y (numeric vector)
# - pos.Y:   indices of entries of x to which sorted values of Y belong to 
#            (integer vector)
{
	N <- length(X)
	
	# Sort everything with respect to order of X
	order.X <- order(X, Y)
	X <- X[order.X]
	W <- W[order.X]
	Y <- Y[order.X]
	
	# Determine vector x of unique values of the sorted vector X,
	# the positions of the latter within x and the vector w of
	# weights of the components of x:
	x <- c(X[1], numeric(N - 1))
	pos.X <- integer(N)
	w <- numeric(N)
	m.X <- 1
	for (k in seq_len(N)) {
		if (X[k] > x[m.X]) {
			m.X <- m.X + 1
			x[m.X] <- X[k]
		}
		pos.X[k] <- m.X
		w[m.X] <- w[m.X] + W[k]
	}
	x <- x[seq_len(m.X)]
	w <- w[seq_len(m.X)]
	
	# Sort Y and keep its ordering permutation
	# (ordering permutation *after* sorting
	#  with respect to order of X)
	order.Y <- order(Y)
	Y <- Y[order.Y]
	
	# Determine vector y of unique values of Y:
	y <- c(Y[1], numeric(N - 1))
	m.Y <- 1
	for (k in seq_len(N)) {
		if (Y[k] > y[m.Y]) {
			m.Y <- m.Y + 1
			y[m.Y] <- Y[k]
		}
	}
	y <- y[seq_len(m.Y)]
	
	# Format for output
	pos.Y <- pos.X[order.Y]
	W <- W[order.Y]
	
	return(list(x=x, y=y, w=w, W=W, Y=Y, pos.Y=pos.Y))
}

isoCdf <- function(X, Y, W = rep(1, length(X)), beta=NULL)
# Compute the isotonic distributional regression of Y
# conditional on X using abridged PAVA.
# 
# Input:
# - X:       covariate vector (numeric vector)
# - Y:       response vector (numeric vector)
# - W:       weight vector (positive numeric vector)
# - beta:    optional vector of probabilities in (0,1)
#            for which the quantile curve should be returned.
#
# The inputs X, Y and W must all have the same length N.
#
# Output:
# - x       sorted distinct values of X
#           (numeric or ordered factor vector)
# - y       sorted distinct values of Y (numeric vector)
# - cdf     matrix containing conditional CDF of Y given X = x
#           evaluated at y
#           (numeric matrix; row i is the CDF conditional on
#            X = x[i])
# - CA:     The number of times an element was added to the
#           partition or modified.
# - CM:     The number of times two elements of the partition
#           were merged.
{
	# Prepare data
	input <- prepareData(X, Y, W)
	x <- input$x
	y <- input$y
	w <- input$w
	m <- length(x)
	W <- input$W
	Y <- input$Y
	pos.Y <- input$pos.Y
	N <- length(Y)
  
	# Output is returned as matrix 'cdf'
	cdf <- matrix(nrow = m, ncol = length(y), 1)
	
	# Initialization:
	z <- rep(0,m)
	PP <- c(0,m,rep(NA,m-1))
	WW <- c(0,sum(w),rep(NA,m-1))
	MM <- c(Inf,0,rep(NA,m-1))
	d <- 2
	CA <- 1
	CM <- 0
	
	# Induction:
	i <- 1 # next observation
	K <- 0 # number of columns of cdf which have been computed
	while (Y[i] < Y[N]){
		j0 <- pos.Y[i]
		# Update of z:
		z[j0] <- z[j0] + W[i]/w[j0]    
		# Update of PP, WW, MM:
		s0 <- 1
		while (PP[s0] < j0){
			s0 <- s0 + 1
		}
		b0 <- PP[s0]
		a0 <- PP[s0-1] + 1
		# Update of CA:
		CA <- CA + d - s0 + b0 - a0 + 1
		if (s0 < d){
			# remember the intervals to the right of P_s0:
			rem <- (s0+1):d
			PPrem <- PP[rem]
			WWrem <- WW[rem]
			MMrem <- MM[rem]
			drem <- d - s0
		}
		d <- s0
		Pdnew <- a0:j0
		PP[d] <- j0
		WW[d] <- sum(w[Pdnew])
		MM[d] <- sum(w[Pdnew]*z[Pdnew])/WW[d]
		while (MM[d-1] <= MM[d]){
			# PAV:
			dnew <- d-1
			PP[dnew] <- PP[d]
			MM[dnew] <- WW[dnew]*MM[dnew] + WW[d]*MM[d]
			WW[dnew] <- WW[dnew] + WW[d]
			MM[dnew] <- MM[dnew]/WW[dnew]
			d <- dnew
			CM <- CM + 1
		}
		while (PP[d] < b0){
			# add a singleton or more to
			# the partition PP:
			dnew <- d + 1
			pp <- PP[d] + 1
			ww <- w[pp]
			MM[dnew] <- z[pp]
			while (pp < b0 & z[pp + 1] == MM[dnew]){
				pp <- pp + 1
				ww <- ww + w[pp]
			}
			PP[dnew] <- pp
			WW[dnew] <- ww
			d <- dnew
			while (MM[d-1] <= MM[d]){
				# PAV:
				dnew <- d-1
				PP[dnew] <- PP[d]
				MM[dnew] <- WW[dnew]*MM[dnew] + WW[d]*MM[d]
				WW[dnew] <- WW[dnew] + WW[d]
				MM[dnew] <- MM[dnew]/WW[dnew]
				d <- dnew
				CM <- CM + 1
			}    	
		}
		if (PP[d] < m){
			# add the remaining intervals from the old partition:
			dnew <- d + drem
			srem <- (d+1):dnew
			PP[srem] <- PPrem
			WW[srem] <- WWrem
			MM[srem] <- MMrem
			d <- dnew
		}
		if (Y[i] < Y[i+1]){
			# update cdf:
			K <- K+1
			cdf[,K] <- rep(MM[2:d], times = diff(PP[seq_len(d)]))
		}
		i <- i+1
	}

	if (is.null(beta)){
	  Q <- NULL
	}else{
	  Q <- matrix(NA,m,length(beta))
	  dimnames(Q)[[1]] <- paste('x_',1:m,sep='')
	  dimnames(Q)[[2]] <- paste('Q',beta,sep='')
	  for (j in 1:length(beta)){
	    qL <- quant_L(y, m, cdf, beta[j])
	    qU <- quant_U(y, m, cdf, beta[j])
	    if (qL[m] <= qU[1]){
	      Q[,j] <- (qL[M] + qU[1])/2
	    }else{
	      Q[,j] <- TautString0(pmax(qL,qU[1]),pmin(qU,qL[m]),x)$Y
	    }
	  }
	}
	return(list(x=x, y=y, cdf=cdf, CA=CA, CM=CM, Q=Q))
}


isoCdfNaive <- function(X, Y, W=rep(1,length(X)))
# Compute the isotonic distributional regression of Y
# conditional on X without abridged PAVA.
# 
# Input:
# - X:       covariate vector (numeric vector)
# - Y:       response vector (numeric vector)
# - W:       weight vector (positive numeric vector)
#
# The inputs X, Y and W must all have the same length N.
#
# Output:
# - x       sorted distinct values of X
#           (numeric or ordered factor vector)
# - y       sorted distinct values of Y (numeric vector)
# - cdf     matrix containing conditional CDF of Y given X = x
#           evaluated at y
#           (numeric matrix; row i is the CDF conditional on
#            X = x[i])
# - CA:     The number of times an element was added to the
#           partition or modified.
# - CM:     The number of times two elements of the partition
#           were merged.
{
  input <- prepareData(X, Y, W)
  x <- input$x
  y <- input$y
  w <- input$w
  m <- length(x)
  W <- input$W
  Y <- input$Y
  pos.Y <- input$pos.Y
  N <- length(Y)
  
  # Output is returned as matrix 'cdf'
  cdf <- matrix(nrow = m, ncol = length(y), 1)
  CA <- 0
  CM <- 0
  
  # Induction:
  z <- rep(0,m)
  i <- 1 # next observation
  K <- 0 # number of columns of cdf which have been computed so far
  while (Y[i] < Y[N]){
    j0 <- pos.Y[i]
    # Update of z:
    z[j0] <- z[j0] + W[i]/w[j0]
    if (Y[i] < Y[i+1]){
      tmp <- PAVA(z=z,w=w)
      d <- tmp$d
      PP <- tmp$PP
      MM <- tmp$MM
      K <- K+1
      cdf[,K] <- rep(MM[2:(d + 1)], times = diff(PP[seq_len(d + 1)]))
      CA <- CA + tmp$CA
      CM <- CM + tmp$CM
    }
    i <- i+1
  }
  return(list(x=x, y=y, cdf=cdf, CA=CA, CM=CM))
}

isoCdfNaive2 <- function(X, Y, W=rep(1,length(X)))
# Compute the isotonic distributional regression of Y
# conditional on X without
# abridged PAVA.
# 
# Input:
# - X:       covariate vector (numeric vector)
# - Y:       response vector (numeric vector)
# - W:       weight vector (positive numeric vector)
#
# The inputs X, Y and W must all have the same length N.
#
# Output:
# - x       sorted distinct values of X (numeric or ordered factor vector)
# - y       sorted distinct values of Y (numeric vector)
# - cdf     matrix containing conditional CDF of Y given X = x evaluated at y
#           (numeric matrix; row i is the CDF conditional on X = x[i])
# - CA:     The number of times an element was added to the partition or
#           modified.
# - CM:     The number of times two elements of the partition were merged.
{
  input <- prepareData(X, Y, W)
  x <- input$x
  y <- input$y
  w <- input$w
  m <- length(x)
  W <- input$W
  Y <- input$Y
  pos.Y <- input$pos.Y
  N <- length(Y)
  
  # Output is returned as matrix 'cdf'
  cdf <- matrix(nrow = m, ncol = length(y), 1)
  CA <- 0
  CM <- 0
  
  # Induction:
  z <- rep(0,m)
  i <- 1 # next observation
  K <- 0 # number of columns of cdf which have been computed so far
  while (Y[i] < Y[N]){
    j0 <- pos.Y[i]
    # Update of z:
    z[j0] <- z[j0] + W[i]/w[j0]
    if (Y[i] < Y[i+1]){
      tmp <- PAVA2(z=z,w=w)
      d <- tmp$d
      PP <- tmp$PP
      MM <- tmp$MM
      K <- K+1
      cdf[,K] <- rep(MM[2:(d + 1)], times = diff(PP[seq_len(d + 1)]))
      CA <- CA + tmp$CA
      CM <- CM + tmp$CM
    }
    i <- i+1
  }
  return(list(x=x, y=y, cdf=cdf, CA=CA, CM=CM))
}

isoQuant <- function(X, Y, W = rep(1, length(X)),
	beta=c(0.1,0.25,0.5,0.75,0.9))
# Compute isotonic regression quantiles of Y conditional on X using
# abridged PAVA.
# The difference to isoCdf is that we don't compute and store all
# estimated conditional distribution functions in a huge matrix.
# 
# Input:
# - X:       covariate vector (numeric vector)
# - Y:       response vector (numeric vector)
# - W:       weight vector (positive numeric vector)
# - beta:    vector of strictly increasing values in (0,1)
#            for which the quantile curve(s) should be computed
#
# The inputs X, Y and W must all have the same length N.
#
# Output:
# - x       sorted distinct values of X (numeric or ordered factor vector)
# - y       sorted distinct values of Y (numeric vector)
# - Qmin, Qmax, Q
#           matrices containing the minimal, the maximal and the smoothest
#           beta-quantiles
{
	# Prepare data
	input <- prepareData(X, Y, W)
	x <- input$x
	y <- input$y
	w <- input$w
	m <- length(x)
	W <- input$W
	Y <- input$Y
	pos.Y <- input$pos.Y
	N <- length(Y)
	
	# Initialize output:
	Q <- matrix(NA, nrow = m, ncol = length(beta))
	dimnames(Q)[[1]] <- paste('x_',1:m,sep='')
	dimnames(Q)[[2]] <- paste('Q',beta,sep='')
	Qmin <- Q
	Qmin[] <- Y[1]
	Qmax <- Q
	Qmax[] <- Y[N]
	
	# Initialization:
	z <- rep(0,m)
	cdf <- rep(0,m)
	PP <- c(0,m,rep(NA,m-1))
	WW <- c(0,sum(w),rep(NA,m-1))
	MM <- c(Inf,0,rep(NA,m-1))
	d <- 2
	
	# Induction:
	betamax <- max(beta)
	i <- 1 # next observation
	while (Y[i] < Y[N] && min(cdf) < betamax){
		j0 <- pos.Y[i]
		# Update of z:
		z[j0] <- z[j0] + W[i]/w[j0]    
		# Update of PP, WW, MM:
		s0 <- 1
		while (PP[s0] < j0){
			s0 <- s0 + 1
		}
		b0 <- PP[s0]
		a0 <- PP[s0-1] + 1
		if (s0 < d){
			# remember the intervals to the right of P_s0:
			rem <- (s0+1):d
			PPrem <- PP[rem]
			WWrem <- WW[rem]
			MMrem <- MM[rem]
			drem <- d - s0
		}
		d <- s0
		Pdnew <- a0:j0
		PP[d] <- j0
		WW[d] <- sum(w[Pdnew])
		MM[d] <- sum(w[Pdnew]*z[Pdnew])/WW[d]
		while (MM[d-1] <= MM[d]){
			# PAV:
			dnew <- d-1
			PP[dnew] <- PP[d]
			MM[dnew] <- WW[dnew]*MM[dnew] + WW[d]*MM[d]
			WW[dnew] <- WW[dnew] + WW[d]
			MM[dnew] <- MM[dnew]/WW[dnew]
			d <- dnew
		}
		while (PP[d] < b0){
			# add a singleton to the partition PP:
			dnew <- d + 1
			pp <- PP[d] + 1
			ww <- w[pp]
			MM[dnew] <- z[pp]
			while (pp < b0 & z[pp + 1] == MM[dnew]){
				pp <- pp + 1
				ww <- ww + w[pp]
			}
			PP[dnew] <- pp
			WW[dnew] <- ww
			d <- dnew
			while (MM[d-1] <= MM[d]){
				# PAV:
				dnew <- d-1
				PP[dnew] <- PP[d]
				MM[dnew] <- WW[dnew]*MM[dnew] + WW[d]*MM[d]
				WW[dnew] <- WW[dnew] + WW[d]
				MM[dnew] <- MM[dnew]/WW[dnew]
				d <- dnew
			}    	
		}
		if (PP[d] < m){
			# add the remaining intervals from the
			# old partition:
			dnew <- d + drem
			srem <- (d+1):dnew
			PP[srem] <- PPrem
			WW[srem] <- WWrem
			MM[srem] <- MMrem
			d <- dnew
		}
		if (Y[i] < Y[i+1]){
			# update quantiles:
			cdf <- rep(MM[2:d], times = diff(PP[seq_len(d)]))
			for (j in 1:length(beta)){
				Qmin[cdf <  beta[j],j] <- Y[i+1]
				tmp <- (Qmax[,j] > Y[i] & cdf > beta[j])
				Qmax[tmp,j] <- Y[i]
			}
		}
		i <- i+1
	}
	
	for (j in 1:length(beta)){
		if (sum(Qmin[,j] < Qmax[,j]) == 0){
			Q[,j] <- Qmin[,j]
		}else{
			if (Qmax[1,j] >= Qmin[m,j]){
				Q[,j] <- (Qmin[m,j] + Qmax[1,j])/2
			}else{
				Q[,j] <- TautString0(
					Y.low=pmax(Qmax[1,j],Qmin[,j]),
					Y.upp=pmin(Qmax[,j],Qmin[m,j]),x=x)$Y
			}
		}
	}  
	return(list(x=x, y=y, Qmin=Qmin, Qmax=Qmax, Q=Q))
}


TautString0 <- function(Y.low,Y.upp,
                        x=1:length(Y.low))
# Input: Three vectors Y.low, Y.upp of length n > 2
# such that
#    Y.low[i] <= Y.upp[i]  for all i
# with equality for i in {1,n}, and
#    x[1] < x[2] < ... < x[n] .
# Output: The corresponding taut string Y.
# That means, Y is a vector of length n such that
#    Y.low[i] <= Y[i] <= Y.upp[i]  for all i ,
# for 1 < i < n,
#    Y[i] = Y.low[i]  if
#       dY[i-1]/dx[i-1] > dY[i]/dx[i] ,
#    Y[i] = Y.upp[i]  if
#       dY[i-1]/dx[i-1] < dY[i]/dx[i] ,
# where dz[j] := z[j+1] - z[j] for any vector z and
# 1 <= j < length(z).
{
	n <- length(Y.low)
	a <- 1
	# Largest index up to which lower and upper
	# string coincide.
	k.low <- c(1,2,rep(NA,n-2))
	# indices of knots of lower string
	slope.low <- c(Inf,
		(Y.low[2] - Y.low[1])/(x[2] - x[1]),
		rep(NA,n-2))
	#  left slopes of the lower string at knots
	b.low <- 2
	# index of last knot of lower string
	k.upp <- c(1,2,rep(NA,n-2))
	# indices of knots of upper string
	slope.upp <- c(-Inf,
		(Y.upp[2] - Y.upp[1])/(x[2] - x[1]),
		rep(NA,n-2))
	# left slopes of the upper string at knots
	b.upp <- 2
	# index of last know of upper string
	for (c in 3:n){
		# Extend lower string:
		b.low <- b.low + 1
		k.low[b.low] <- c
		slope.low[b.low] <-
			(Y.low[c] - Y.low[c-1])/(x[c] - x[c-1])
		# Pull lower string:
		while (b.low > a + 1 &&
			slope.low[b.low] >= slope.low[b.low-1]){
			slope.low[b.low-1] <-
				(Y.low[c] - Y.low[k.low[b.low-2]])/
				(x[c] - x[k.low[b.low-2]])
			k.low[b.low-1] <- c
			b.low <- b.low-1
		}
		
		# Extend upper string:
		b.upp <- b.upp + 1
		k.upp[b.upp] <- c
		slope.upp[b.upp] <-
			(Y.upp[c] - Y.upp[c-1])/(x[c] - x[c-1])
		# Pull upper string:
		while (b.upp > a + 1 &&
			slope.upp[b.upp] <= slope.upp[b.upp-1]){
			slope.upp[b.upp-1] <-
				(Y.upp[c] - Y.upp[k.upp[b.upp-2]])/
				(x[c] - x[k.upp[b.upp-2]])
			k.upp[b.upp-1] <- c
			b.upp <- b.upp-1
		}
		
		# Check whether lower and upper string are still
		# ordered and modifiy them, if necessary:
		if (b.low == a+1){
			# Bend lower string, if necessary:
			while (b.upp > a + 1 &&
				slope.low[b.low] >= slope.upp[a+1]){
				a <- a + 1
				b.low <- b.low + 1
				k.low[a] <- k.upp[a]
				slope.low[a] <- slope.upp[a]
				k.low[b.low] <- c
				slope.low[b.low] <-
					(Y.low[c] - Y.upp[k.upp[a]])/
					(x[c] - x[k.upp[a]])
				Y.low[k.low[a]] <- Y.upp[k.upp[a]]
			}
		}
		if (b.upp == a+1){
			# Bend upper string, if necessary:
			while (b.low > a + 1 &&
				slope.upp[b.upp] <= slope.low[a+1]){
				a <- a + 1
				b.upp <- b.upp + 1
				k.upp[a] <- k.low[a]
				slope.upp[a] <- slope.low[a]
				k.upp[b.upp] <- c
				slope.upp[b.upp] <-
					(Y.upp[c] - Y.low[k.low[a]])/
					(x[c] - x[k.low[a]])
				Y.upp[k.upp[a]] <- Y.low[k.low[a]]
			}
		}
	}
	
	# Transform slope.low (== slope.upp) into a
	# proper vector of length n:
	slope <- c(NA,rep(0,n-1))
	a <- 1
	for (b in 2:b.low){
		slope[(a+1):k.low[b]] <- slope.low[b]
		a <- k.low[b]
	}
	
	# Compute vector Y from its slopes:
	Y <- Y.low
	dx <- x[2:n] - x[1:(n-1)]
	Y[2:n] <- Y[1] + cumsum(dx * slope[2:n])
	
	return(list(knot.ind=k.low,Y=Y))
}

quant_L <- function(y, m, hat_F, beta)
# Lower quantile function
{
	q <- rep(0,m)
	for (j in 1:m){
		q[j] <- min(y[beta <= hat_F[j,]])
	}
	return(q)
}

quant_U <- function(y, m, hat_F, beta) 
# Upper quantile function
{
	q <- rep(0,m)
	for (j in 1:m){
		q[j] <- min(y[beta < hat_F[j,]])
	}
	return(q)
}


# Demo versions of isoCdf and isoCDFnaive:

isoCdf_demo <- function(Y)
# Illustrates isoCdf with X = 1:length(Y)
# and W = rep(1,length(Y)).
{
	Y0 <- Y
	# Prepare data
	input <- prepareData(1:length(Y),Y,
		rep(1,length(Y)))
	x <- input$x
	y <- input$y
	w <- input$w
	m <- length(x)
	W <- input$W
	Y <- input$Y
	pos.Y <- input$pos.Y
	N <- length(Y)
	
	par(mfrow=c(2,1),mai=c(0.9,0.9,0.01,0.01))
	
	# Output is returned as matrix 'cdf'
	cdf <- matrix(nrow = m, ncol = length(y), 1)
	
	# Initialization:
	z <- rep(0,m)
	PP <- c(0,m,rep(NA,m-1))
	WW <- c(0,sum(w),rep(NA,m-1))
	MM <- c(Inf,0,rep(NA,m-1))
	d <- 2
	CA <- 1
	CM <- 0
	
	tmp <- readline('Start: ')
	plot(1:N,Y0,xlim=c(0,N),xlab='X',ylab='Y')
	plot(cumsum(w),cumsum(w*z),
		xlab='i',ylab='S(i)',
		xlim=c(0,sum(w)),ylim=c(0,sum(w)))
	lines(PP[1:d],cumsum(WW[1:d]*c(0,MM[2:d])),
		col='blue')
	points(PP[1:d],cumsum(WW[1:d]*c(0,MM[2:d])),
		pch=16,col='blue')
	# Induction:
	i <- 1 # next observation
	K <- 0 # number of columns of cdf which have been computed
	while (Y[i] < Y[N]){
		j0 <- pos.Y[i]
		# Update of z:
		z[j0] <- z[j0] + W[i]/w[j0]    
		tmp <- readline('Increase y: ')
		plot(1:N,Y0,xlim=c(0,N),xlab='X',ylab='Y')
		abline(h=Y[i],col='magenta')
		abline(v=j0,col='red',lty=3)
		plot(cumsum(w),cumsum(w*z),
			xlab='i',ylab='S(i)',
			xlim=c(0,sum(w)),ylim=c(0,sum(w)))
		lines(PP[1:d],cumsum(WW[1:d]*c(0,MM[2:d])),
			col='red')
		points(PP[1:d],cumsum(WW[1:d]*c(0,MM[2:d])),
			pch=16,col='red')
		# Update of PP, WW, MM:
		s0 <- 1
		while (PP[s0] < j0){
			s0 <- s0 + 1
		}
		b0 <- PP[s0]
		abline(v=j0,col='red',lty=3)
		abline(v=b0, col='forestgreen',lty=3)
		a0 <- PP[s0-1] + 1
		# Update of CA:
		CA <- CA + d - s0 + b0 - a0 + 1
		if (s0 < d){
			# remember the intervals to the right of P_s0:
			rem <- (s0+1):d
			PPrem <- PP[rem]
			WWrem <- WW[rem]
			MMrem <- MM[rem]
			drem <- d - s0
		}
		d <- s0
		Pdnew <- a0:j0
		PP[d] <- j0
		WW[d] <- sum(w[Pdnew])
		MM[d] <- sum(w[Pdnew]*z[Pdnew])/WW[d]
		tmp <- readline('Update: ')
		plot(1:N,Y0,xlim=c(0,N),xlab='X',ylab='Y')
		abline(h=Y[i],col='magenta')
		abline(v=j0,col='red',lty=3)
		plot(cumsum(w),cumsum(w*z),
			xlab='i',ylab='S(i)',
			xlim=c(0,sum(w)),ylim=c(0,sum(w)))
		lines(PP[1:d],cumsum(WW[1:d]*c(0,MM[2:d])),
			col='blue')
		points(PP[1:d],cumsum(WW[1:d]*c(0,MM[2:d])),
			pch=16,col='blue')
		abline(v=j0,col='red',lty=3)
		abline(v=b0,col='forestgreen',lty=3)
		while (MM[d-1] <= MM[d]){
			# PAV:
			dnew <- d-1
			PP[dnew] <- PP[d]
			MM[dnew] <- WW[dnew]*MM[dnew] + WW[d]*MM[d]
			WW[dnew] <- WW[dnew] + WW[d]
			MM[dnew] <- MM[dnew]/WW[dnew]
			d <- dnew
			tmp <- readline('PAV: ')
			plot(1:N,Y0,xlim=c(0,N),xlab='X',ylab='Y')
			abline(h=Y[i],col='magenta')
			abline(v=j0,col='red',lty=3)
			plot(cumsum(w),cumsum(w*z),
				xlab='i',ylab='S(i)',
				xlim=c(0,sum(w)),ylim=c(0,sum(w)))
			lines(PP[1:d],cumsum(WW[1:d]*c(0,MM[2:d])),
				col='blue')
			points(PP[1:d],cumsum(WW[1:d]*c(0,MM[2:d])),
				pch=16,col='blue')
			abline(v=j0,col='red',lty=3)
			abline(v=b0,col='forestgreen',lty=3)
			CM <- CM + 1
		}
		while (PP[d] < b0){
			# add a singleton or more
			# to the partition PP:
			dnew <- d + 1
			pp <- PP[d] + 1
			ww <- w[pp]
			MM[dnew] <- z[pp]
			while (pp < b0 & z[pp + 1] == MM[dnew]){
				pp <- pp + 1
				ww <- ww + w[pp]
			}
			PP[dnew] <- pp
			WW[dnew] <- ww
			d <- dnew
			tmp <- readline('Update: ')
			plot(1:N,Y0,xlim=c(0,N),xlab='X',ylab='Y')
			abline(h=Y[i],col='magenta')
			abline(v=j0,col='red',lty=3)
			plot(cumsum(w),cumsum(w*z),
				xlab='i',ylab='S(i)',
				xlim=c(0,sum(w)),ylim=c(0,sum(w)))
			lines(PP[1:d],cumsum(WW[1:d]*c(0,MM[2:d])),
				col='blue')
			points(PP[1:d],cumsum(WW[1:d]*c(0,MM[2:d])),
				pch=16,col='blue')
			abline(v=j0,col='red',lty=3)
			abline(v=b0,col='forestgreen',lty=3)
			while (MM[d-1] <= MM[d]){
				# PAV:
				dnew <- d-1
				PP[dnew] <- PP[d]
				MM[dnew] <- WW[dnew]*MM[dnew] + WW[d]*MM[d]
				WW[dnew] <- WW[dnew] + WW[d]
				MM[dnew] <- MM[dnew]/WW[dnew]
				d <- dnew
				tmp <- readline('PAV: ')
				plot(1:N,Y0,xlim=c(0,N),xlab='X',ylab='Y')
				abline(h=Y[i],col='magenta')
				abline(v=j0,col='red',lty=3)
				plot(cumsum(w),cumsum(w*z),
					xlab='i',ylab='S(i)',
					xlim=c(0,sum(w)),ylim=c(0,sum(w)))
				lines(PP[1:d],cumsum(WW[1:d]*c(0,MM[2:d])),
					col='blue')
				points(PP[1:d],cumsum(WW[1:d]*c(0,MM[2:d])),
					pch=16,col='blue')
				abline(v=j0,col='red',lty=3)
				abline(v=b0,col='forestgreen',lty=3)
				CM <- CM + 1
			}    	
		}
		if (PP[d] < m){
			# add the remaining intervals from the old partition:
			dnew <- d + drem
			srem <- (d+1):dnew
			PP[srem] <- PPrem
			WW[srem] <- WWrem
			MM[srem] <- MMrem
			d <- dnew
		}

		tmp <- readline('Update: ')
		plot(1:N,Y0,xlim=c(0,N),xlab='X',ylab='Y')
		abline(h=Y[i],col='magenta')
		abline(v=j0,col='red',lty=3)
		plot(cumsum(w),cumsum(w*z),
			xlab='i',ylab='S(i)',
			xlim=c(0,sum(w)),ylim=c(0,sum(w)))
		lines(PP[1:d],cumsum(WW[1:d]*c(0,MM[2:d])),
			col='blue')
		points(PP[1:d],cumsum(WW[1:d]*c(0,MM[2:d])),
			pch=16,col='blue')
		
		if (Y[i] < Y[i+1]){
			# update cdf:
			K <- K+1
			cdf[,K] <- rep(MM[2:d], times = diff(PP[seq_len(d)]))
		}
		i <- i+1
	}
	tmp <- readline('End: ')
	plot(1:N,Y0,xlim=c(0,N),xlab='X',ylab='Y')
	abline(h=Y[N],col='magenta')
	plot(cumsum(w),cumsum(w),
		xlab='i',ylab='S(i)',
		xlim=c(0,sum(w)),ylim=c(0,sum(w)))
	lines(c(0,sum(w)),c(0,sum(w)),
		col='blue')
	points(c(0,sum(w)),c(0,sum(w)),
		pch=16,col='blue')

	return(list(x=x, y=y, cdf=cdf, CA=CA, CM=CM))
}


PAVA_demo <- function(z)
# Illustrate PAVA with w = rep(1,\length(z)).
{
	par(mfrow=c(2,1),mai=c(0.9,0.9,0.01,0.01))
	m <- length(z)
	PP <- c(0, rep(0,m))
	WW <- c(0, rep(0,m))
	MM <- c(Inf, rep(0,m))
	tmp <- readline('Start: ')
	plot(1:m,z,xlab='i',ylab='z(i)',xlim=c(0,m))
	plot(0:m,cumsum(c(0,z)),xlab='i',ylab='S(i)')
	points(0,0,pch=16,col='blue')

	if (m == 1){
		return(list(d = 1, PP = c(0,1), WW = c(0,1),
			MM = c(Inf,z[1]), CA = 1, CM = 0))
	}
	
	d <- 1
	# During the algorithm, d is the number of relevant elements in the
	# partition + 1
	for (i in 1:m) {
    d <- d + 1
    PP[d] <- i
    WW[d] <- 1
    MM[d] <- z[i]
	tmp <- readline('New knot: ')
	plot(1:m,z,xlab='i',ylab='z(i)',xlim=c(0,m))
	for (i in 2:d){
		lines(c(PP[i-1]+0.6,PP[i]+0.4),
			rep(MM[i],2),col='blue',lwd=2)
	}
	plot(0:m,cumsum(c(0,z)),xlab='i',ylab='S(i)')
	lines(PP[1:d],cumsum(c(0,WW[2:d]*MM[2:d])),
		col='blue')
	points(PP[1:d],cumsum(c(0,WW[2:d]*MM[2:d])),
		pch=16,col='blue')

    while (MM[d - 1] <= MM[d]) {
      dd <- c(d - 1, d)
      d <- d - 1
      MM[d] <- sum(WW[dd] * MM[dd])
      WW[d] <- sum(WW[dd])
      MM[d] <- MM[d]/WW[d]
      PP[d] <- PP[d + 1]
		tmp <- readline('Pull string: ')
		plot(1:m,z,xlab='i',ylab='z(i)',xlim=c(0,m))
		for (i in 2:d){
			lines(c(PP[i-1]+0.6,PP[i]+0.4),
				rep(MM[i],2),col='blue',lwd=2)
		}
		plot(0:m,cumsum(c(0,z)),xlab='i',ylab='S(i)')
		lines(PP[1:d],cumsum(c(0,WW[2:d]*MM[2:d])),
			col='blue')
		points(PP[1:d],cumsum(c(0,WW[2:d]*MM[2:d])),
			pch=16,col='blue')
	}
  }
  d <- d - 1
  # In the end we return the number of relevant intervals
  # in the partition
	
	return(list(d = d, PP = PP, WW = WW,
		MM = MM, CA = m, CM = m-d))
}

PAVA2_demo <- function(z)
# Illustrate PAVA with w = rep(1,\length(z)).
{
	par(mfrow=c(2,1),mai=c(0.9,0.9,0.01,0.01))
	m <- length(z)
	PP <- c(0, rep(0,m))
	WW <- c(0, rep(0,m))
	MM <- c(Inf, rep(0,m))
	tmp <- readline('Start: ')
	plot(1:m,z,xlab='i',ylab='z(i)',xlim=c(0,m))
	plot(0:m,cumsum(c(0,z)),xlab='i',ylab='S(i)')
	points(0,0,pch=16,col='blue')

	if (m == 1){
		return(list(d = 1, PP = c(0,1), WW = c(0,1),
			MM = c(Inf,z[1]), CA = 1, CM = 0))
	}
	
	CA <- 0
	CM <- 0
	d <- 1
	# During the algorithm, d is the number of relevant elements in the
	# partition + 1
	while (PP[d] < m) {
		pp <- PP[d] + 1
		ww <- 1
		d <- d+1
		MM[d] <- z[pp]
		while (pp < m && z[pp+1] == MM[d]){
			pp <- pp + 1
			ww <- ww + 1
		}
		PP[d] <- pp
		WW[d] <- ww
		CA <- CA + 1
		tmp <- readline('New knot: ')
		plot(1:m,z,xlab='i',ylab='z(i)',xlim=c(0,m))
		for (i in 2:d){
			lines(c(PP[i-1]+0.6,PP[i]+0.4),
				rep(MM[i],2),col='blue',lwd=2)
		}
		plot(0:m,cumsum(c(0,z)),xlab='i',ylab='S(i)')
		lines(PP[1:d],cumsum(c(0,WW[2:d]*MM[2:d])),
			col='blue')
		points(PP[1:d],cumsum(c(0,WW[2:d]*MM[2:d])),
			pch=16,col='blue')
		
		while (MM[d - 1] <= MM[d]) {
			dd <- c(d - 1, d)
			d <- d - 1
			MM[d] <- sum(WW[dd] * MM[dd])
			WW[d] <- sum(WW[dd])
			MM[d] <- MM[d]/WW[d]
			PP[d] <- PP[d + 1]
			CM <- CM + 1
			tmp <- readline('Pull string: ')
			plot(1:m,z,xlab='i',ylab='z(i)',xlim=c(0,m))
			for (i in 2:d){
				lines(c(PP[i-1]+0.6,PP[i]+0.4),
					rep(MM[i],2),col='blue',lwd=2)
			}
			plot(0:m,cumsum(c(0,z)),xlab='i',ylab='S(i)')
			lines(PP[1:d],cumsum(c(0,WW[2:d]*MM[2:d])),
				col='blue')
			points(PP[1:d],cumsum(c(0,WW[2:d]*MM[2:d])),
				pch=16,col='blue')
		}
	}
	d <- d - 1
	# In the end we return the number of relevant intervals
	# in the partition
	
	return(list(d = d, PP = PP, WW = WW,
		MM = MM, CA = CA, CM = CM))
}