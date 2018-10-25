source('Derivatives.R', chdir = TRUE)

#' Range Check
#'
#' Checks if a value is within the range [x,y]
#' @param x lower value of range
#' @param y upper value of range
#' @param a value to evaluate if falls within the range [x,y]

range.check	<-	function(x,y,a){
	
	if((x<=a)&(a<=y)){ret	<-	1}else{ret	<-	0}
	return(ret)
}



#' Expit transformation
#'
#' expit(x) = frac(1,(1+exp(-x)))
#' @param x  a single number or a vector to be transformed
#' @export
#' @examples
#' x	<-	rnorm(40)
#' expit(x)

expit	<-	function(x){1/(1+exp(-x))}

#' Logit transformation
#'
#' Logit transformation for data
#' logit(x) = log(x)-log(1-x)
#' @param x  a single number or a vector to be transformed
#' @export
#' @examples
#' x	<-	rnorm(40)
#' logit(x)
logit	<-	function(x){log(x)-log(1-x)}

# pgf of X_t
# This is the pgf og the number of LIZs in a specific quadrant in time t given the number of LIZs in the same
# quadrant in time t-1

phi.g <- function(z,a,b,x,gam){
	
	((1-a + a*z)^x)*((b+a*(1-z))/(b+1-z))^gam
}

#'	Llikelihood of X_t > 8 
#'
#'	Llikelihood of X_t>8 given X_{t-1}
#'	@param a = alpha
#'	@param b = beta
#'	@param x = x_{t-1}
#'	@param gam = gamma

phi.g.9p	<-	function(z,a,b,x,gam){
	
	len.x	<-	length(x)
	prob.x	<-	matrix(0,nrow=10,ncol=len.x)
	x.vals	<-	0:8		
	
	for(i in 1:length(x.vals)){
		
		ith.x		<-	x.vals[i]
		prob.x[i,]	<-	phi.g.list[[i]](z=z,a=a,b=b,x=x, gam=gam)/factorial(ith.x)
	
	}
	res	<-	1-apply(prob.x,2,sum)
	return(res)
}

phi.g.list	<-	list(phi.g.0 = phi.g, phi.g.1 = phi.g.1
					,phi.g.2 = phi.g.2,phi.g.3 = phi.g.3
					,phi.g.4 = phi.g.4,phi.g.5 = phi.g.5
					,phi.g.6 = phi.g.6,phi.g.7 = phi.g.7
					,phi.g.8 = phi.g.8,phi.g.9 = phi.g.9p)
	

#' Expected value of X_t
#'
#' Expected value of X_t based on the DSD process
#' \code{E[x_t]} = (\code{alpha}*\code{x_{t-1}}) + \code{frac((\code{gamma}*(1-\code{alpha})),\code{beta})}
#' @param a  \code{alpha}
#' @param b  \code{beta}
#' @param gam  \code{gamma}
#' @param psi  \code{psi} (Zero Inflation Parameter)
#' @param x  \code{x_{t-1}}
#' @export
#' @examples
#' a = 0.8
#' b = 3
#' gam = 1.7
#' psi	= 0.9
#' xtm1 = rpois(-gam*log(a))
#' x.t	=	mu.x(a,b,gam,psi,xtm1)

mu.x.nb	<-	function(a,b,gam,psi,x){
		
		psi*((a*x)+(gam*(1-a))/b)
	
	}

mu.x.pois	<-	function(a,psi,lambda,x){
	
	a*x+(1-a)*lambda
	
}

#' Variance of \code{X_t}
#'
#' Variance of X_t based on the prediction of the DSD process
#'	@param	a	\code{alpha}
#'	@param	b	\code{beta}
#'	@param	gam	\code{gamma}
#' @param psi  \code{psi} (Zero Inflation Parameter)
#'	@param	x	\code{x_{t-1}}
#'	@return	Values of the variance of \code{X_t}
#' @export
#' @examples
#' a = 0.8
#' b = 3
#' gam = 1.7
#' psi = 0.9
#' xtm1 = rpois(-gam*log(a))
#' x.t	=	mu.x(a,b,gam,psi,xtm1)
#' var.x.t	=	Var.x(a,b,gam,psi,xtm1)

Var.x.nb	<-	function(a,b,gam,psi,x){
	sigma2.y	<-	(a*x)+(gam*(1-a))/b -(a^2)*x+((gam*(1-(a^2)))/(b^2))
	mu.y2		<-	(((a*x)+(gam*(1-a))/b))^2
	mu.x2		<-	sigma2.y+mu.y2
	mu.x		<-	(a*x)+(gam*(1-a))/b
	var.psi		<-	psi*mu.x2+(mu.x)^2
	return(var.psi)
}

Var.x.pois	<-	function(a,psi,lambda,x){
	
	
	a*(1-a)*x + (1-a)*lambda
	
}


#' Time series data simulation
#'
#' Function for data simulation of LIZs according to the dsd model
#'	@param t.steps  length of the time series
#'	@param t0  An integer value that will be taken as the starting point for the simulation. If t0 is missing, 
#'			   the simulation will start with a random value with probability given by 
#'			   a negative binomial distribution with parameters gamma and beta/alpha
#'	@param alpha		Decay of the counts through time. It can be a single value for all of the 
#'						simulations a vector for varying alpha according to some covariate or a matrix varying alpha for each
#'						time step. 0=<alpha=<1
#'	@param beta		Beta parameter of the negative binomial process. Constant throughout the reps.
#'	@param gamma		gamma parameter of the DSD process. Constant throughout the reps.
#'	@param psi		Zero inflation parameter. It could be also a function of some covariates. Similarly to alpha,
#'					it could be a single number, a vector or a matrix. 0=<psi=<1
#'	@param reps 	 	replicates of the time series process.

#'	@return a matrix with ncol = t.steps, and nrow = reps.
#'	@export
#'	@examples
#'	quadrants <-	40
#'	t.steps		<-	10
#'	b <- 3
#'	gam	<-	1.7 
#'	t0	<-	2
#'	betas.a	<-	c(0.1,2)
#'	betas.psi	<-	c(2,0.9)

#'	# Covariates for alpha and psi
#'	ones	<-	rep(1,quadrants)
#'	set.seed(1234)
#'	x.alpha	<-	cbind(ones,rnorm(quadrants))
#'	alphas	<-	expit(x.alpha%*%betas.a)
#'	set.seed(4521)
#'	x.psi	<-	cbind(ones,rnorm(quadrants))
#'	psi		<-	expit(x.psi%*%betas.psi)

#'	cov		<-	data.frame(x.alpha=x.alpha[,2],x.psi=x.psi[,2])

#'	rdsd(t.steps=t.steps,t0=t0,alpha=alphas,beta=b,gamma=gam,psi=psi,reps=quadrants)	

rdsd	<-	function(t.steps,t0,alpha,beta,gamma,psi,reps){
	
	if(class(alpha)!="matrix"){stop("alpha must be a matrix")}
	if(class(psi)!="matrix"){stop("psi must be a matrix")}
	if(missing(reps)){reps	<-	nrow(alpha)}
	## Set up the parameters for the simulation
	
	if(nrow(alpha)==1)	{alpha		<-	matrix(alpha,ncol=reps,nrow=t.steps)}
	if(nrow(psi)==1)	{psi		<-	matrix(psi,ncol=reps,nrow=t.steps)}
	
	p			<-	beta/(beta+1)
	lam			<-	-gamma*log(alpha)	
	
	## Data matrix and first observation

	tsdata		<-	matrix(NA,ncol=reps,nrow=t.steps)
	
	if(!missing(t0)){
		
			tsdata[1,]	<-	t0
		
		}else{
				
				tsdata[1,]	<-	rnbinom(reps,gamma,prob=beta/(beta+1))
			
			}
	
	
	## The first loop will determine the time series for each quadrant
	#for(t in 1:reps){
	
	## The second loop will loop over the length of the time series for each quadrant
	
	for (i in 1:(t.steps-1)){
			
		## Obtain the number of LIZs in a quadrant given the decay rate alpha for quadrant t
		## In this case alpha can be constant over the quadrants or can be a function of some covariates
		## In the case in which alpha is a function of covariates, alpha should be calculated previous to
		## simulation of data
			
			
		Xt		<-	tsdata[i,]
		
		thinning	<-	rbinom(reps,Xt,alpha[i,])		
						
		## The third loop will determine the number of individuals that immigrate in each wave, 
		## and how many of them will die in the quadrant.
					
		epsilon.t	<-	rep(NA,reps)
		
		## This determines how many immigration waves there are to the quadrant which would add to the 
		## LIZs already present in the quadrant
	
		waves	<-	rpois(reps,lam[i,])
		
		for(j in 1:reps){
	
			Ui				<-	runif(waves[j],min=0,max=1)
			N.imm			<-	rnbinom(waves[j],1,prob=p)
			dead.tot		<-	rbinom(waves[j],N.imm,alpha[i,j]^Ui)
			epsilon.t[j]	<-	sum(dead.tot)
						
		}		
				
			tsdata[i+1,]	<-	(thinning+epsilon.t)
			
		## Now, the zeros could be determined by the suitablity of the quadrant, which might also 
		## be constant or be a function of some covariates. If your don't want to consider the heterogeneity
		## in the suitablity of the quadrants, then set psi = 1.
				
		suitable		<-	rbinom(reps,1,psi[i,])
		not.suitable	<-	which(suitable==0)
		
		tsdata[i+1,not.suitable]	<-	0
					
	}
			
	return(tsdata)
	
}

# Time series simulation under the poisson model

rdsd.pois	<-	function(t.steps,t0,alpha,psi,lambda,reps){
	
	if(class(alpha)!="matrix"){stop("alpha must be a matrix")}
	if(class(psi)!="matrix"){stop("psi must be a matrix")}
	if(missing(reps)){reps	<-	nrow(alpha)}
	## Set up the parameters for the simulation
	
	if(nrow(alpha)==1)	{alpha		<-	matrix(alpha,ncol=reps,nrow=t.steps)}
	if(nrow(psi)==1)	{psi		<-	matrix(psi,ncol=reps,nrow=t.steps)}
	
	## Data matrix and first observation

	tsdata		<-	matrix(NA,ncol=reps,nrow=t.steps)
	
	if(missing(t0)){tsdata[1,]	<-	rpois(reps,lambda)}else{tsdata[1,]	<-	t0}
		
	## The first loop will determine the time series for each quadrant
	#for(t in 1:reps){
	
	## The second loop will loop over the length of the time series for each quadrant
	
	for (i in 1:(t.steps-1)){
			
		## Obtain the number of LIZs in a quadrant given the decay rate alpha for quadrant t
		## In this case alpha can be constant over the quadrants or can be a function of some covariates
		## In the case in which alpha is a function of covariates, alpha should be calculated previous to
		## simulation of data
			
			
		Xt		<-	tsdata[i,]
		
		tsdata[i+1,]	<-	rbinom(reps,Xt,alpha[i,])+rpois(reps,lambda*(1-alpha[i,]))
									
	## Now, the zeros could be determined by the suitablity of the quadrant, which might also 
	## be constant or be a function of some covariates. If your don't want to consider the heterogeneity
	## in the suitablity of the quadrants, then set psi = 1.
				
		suitable		<-	rbinom(reps,1,psi[i,])
		not.suitable	<-	which(suitable==0)
		
		tsdata[i+1,not.suitable]	<-	0
					
	}
			
	return(tsdata)
}	



#' Likelihood function of DSD with covariates
#'
#' Likelihood function of DSD with covariates. This function is to be used in conjunction with the dsd.glm function
#' To be used by itself, alpha.matrix and psi.matrix have to be the design matrix
#' for the relationship of alpha and psi with its covariates, and
#'	@param guess = c(betas.alpha,betas.psi,b,gam)
#'	@param alpha.matrix =	the design matrix usually is the result of the function model.matrix() but can be 
#'						 	specified manually in which the first column should be 1 and the following columns should have the
#'							covariates values.
#'	@param	psi.matrix	=	design matrix of the relationship of psi and its covariates
#'	@return negative log likelihood of the parameters given the time series data.
#'	@export

dsd.like.NB <- function(guess,tsdata,alpha.desmat,psi.desmat){
	
	t.steps		<-	nrow(tsdata)
	quads		<-	ncol(tsdata)
	
	# Alpha		
	a.cov		<-	ncol(alpha.desmat)
	betas.alpha	<-	guess[1:a.cov]
	a.hat		<-	alpha.desmat%*%betas.alpha
	a.hat		<-	matrix(plogis(a.hat),ncol=quads,nrow=t.steps)
	
	# Psi
	psi.cov		<-	ncol(psi.desmat)
	betas.psi	<-	guess[(a.cov+1):(a.cov+psi.cov)]
	psi.hat		<-	psi.desmat%*%betas.psi
	psi.hat		<-	matrix(plogis(psi.hat),ncol=quads,nrow=t.steps)

	# Other Parameters
	b		<-	exp(guess[a.cov+psi.cov+1])
	gam		<-	exp(guess[a.cov+psi.cov+2])
	
	# Log Likelihood
	
	reps	<-	ncol(tsdata)
		
	q 		<- (nrow(tsdata)-1)
	qp1 	<- q+1
	
	lik.mat	<-	matrix(0,nrow=q,ncol=reps)

	tsdata[tsdata>8]	<-	9	
	
			# All Data
	x		<-	tsdata[2:qp1,]
	xm1		<-	tsdata[1:q,]
	
	a.t		<-	a.hat[1:q,]
	psi.t		<-	psi.hat[1:q,]
	
	x.vals	<-	sort(unique(c(x)))
			
			
	for(i in 1:length(x.vals)){
		
		ith.x	<-	x.vals[i]
		x.i		<-	x[x==ith.x]
		xm1.i	<-	xm1[x==ith.x]
		a.i		<-	a.t[x==ith.x]
		psi.i	<-	psi.t[x==ith.x]
		
	if(ith.x==0){lik.mat[x==ith.x]	<-	(1-psi.i)+(psi.i*phi.g.list[[i]](z=0,a=a.i,b=b,x=xm1.i, gam=gam)/factorial(ith.x))
		}else{lik.mat[x==ith.x]	<-	(psi.i*phi.g.list[[i]](z=0,a=a.i,b=b,x=xm1.i, gam=gam)/factorial(ith.x))}
			
	}			

		
	lik.mat[lik.mat<=0]	<-	.Machine$double.xmin
	llik.mat			<-	log(lik.mat)
	
	return(-sum(llik.mat))
	
}

dsd.like.Pois <- function(guess,tsdata,alpha.desmat,psi.desmat){
	
	t.steps		<-	nrow(tsdata)
	quads		<-	ncol(tsdata)
	
	# Alpha		
	a.cov		<-	ncol(alpha.desmat)
	betas.alpha	<-	guess[1:a.cov]
	a.hat		<-	alpha.desmat%*%betas.alpha
	a.hat		<-	matrix(plogis(a.hat),ncol=quads,nrow=t.steps)
	
	# Psi
	psi.cov		<-	ncol(psi.desmat)
	betas.psi	<-	guess[(a.cov+1):(a.cov+psi.cov)]
	psi.hat		<-	psi.desmat%*%betas.psi
	psi.hat		<-	matrix(plogis(psi.hat),ncol=quads,nrow=t.steps)

	# Other Parameters

	lambda	<-	exp(guess[a.cov+psi.cov+1])

	# Log Likelihood
	
	reps	<-	ncol(tsdata)
		
	q 		<- (nrow(tsdata)-1)
	qp1 	<- q+1
	
	lik.mat	<-	matrix(0,nrow=q,ncol=reps)
			
	# All Data
	x		<-	tsdata[2:qp1,]
	xm1		<-	tsdata[1:q,]
	
	a.t		<-	a.hat[1:q,]
	psi.t		<-	psi.hat[1:q,]
	
	super.x	<-	cbind(as.vector(x),as.vector(xm1))
	m.vec	<-	apply(super.x,1,min)
	unique.m<-	sort(unique(m.vec))
			
	for(i in 1:length(unique.m)){
		m.i		<-	unique.m[i]
		x.i		<-	x[m.vec==m.i]
		xm1.i	<-	xm1[m.vec==m.i]
		a.i		<-	a.t[m.vec==m.i]
		psi.i	<-	psi.t[m.vec==m.i]
				
		like1	<-	factorial(xm1.i)*exp(-((1-a.i)*lambda))
		like2.1	<-	like2.2	<-	matrix(0,ncol=m.i+1,nrow=length(like1))
		for(k in 0:m.i){
			like2.1[,k+1]	<- (a.i^k)*((1-a.i)^(xm1.i+x.i-2*k))*lambda^(x.i-k)
			like2.2[,k+1]	<-	factorial(k)*factorial((xm1.i-k))*factorial(x.i-k)
		}
		like2.mat			<-	like2.1/like2.2
		like2				<-	apply(like2.mat,1,sum)
		full.like			<-	numeric(length(x.i))
		full.lik			<-	psi.i*(like1*like2)
		full.lik[x.i==0]	<-	(1-psi.i[x.i==0])+full.lik[x.i==0]
		
		lik.mat[m.vec==m.i]	<-	full.lik
			
	}		
		
	lik.mat[lik.mat<=0]	<-	.Machine$double.xmin
	llik.mat			<-	log(lik.mat)
	
	return(-sum(llik.mat))
	
}




dsd.proflike <- function(guess,tsdata,alpha.desmat,psi.desmat,b,gam){
	
	t.steps		<-	nrow(tsdata)
	quads		<-	ncol(tsdata)
	
	# Alpha		
	a.cov		<-	ncol(alpha.desmat)
	betas.alpha	<-	guess[1:a.cov]
	a.hat		<-	alpha.desmat%*%betas.alpha
	a.hat		<-	matrix(plogis(a.hat),ncol=quads,nrow=t.steps)
	
	# Psi
	psi.cov		<-	ncol(psi.desmat)
	betas.psi	<-	guess[(a.cov+1):(a.cov+psi.cov)]
	psi.hat		<-	psi.desmat%*%betas.psi
	psi.hat		<-	matrix(plogis(psi.hat),ncol=quads,nrow=t.steps)
	
	# Log Likelihood
	tsdata[tsdata>8]	<-	9	
	
	reps	<-	ncol(tsdata)
		
	q 		<- (nrow(tsdata)-1)
	qp1 	<- q+1
	
	lik.mat	<-	matrix(0,nrow=q,ncol=reps)
	
	# All Data
	x		<-	tsdata[2:qp1,]
	xm1		<-	tsdata[1:q,]
	
	a.t		<-	a.hat[1:q,]
	psi.t		<-	psi.hat[1:q,]
	
	x.vals	<-	sort(unique(c(x)))
	
	for(i in 1:length(x.vals)){
		
		ith.x	<-	x.vals[i]
		x.i		<-	x[x==ith.x]
		xm1.i	<-	xm1[x==ith.x]
		a.i		<-	a.t[x==ith.x]
		psi.i	<-	psi.t[x==ith.x]
		
		if(ith.x==0){lik.mat[x==ith.x]	<-	(1-psi.i)+(psi.i*phi.g.list[[i]](z=0,a=a.i,b=b,x=xm1.i, gam=gam)/factorial(ith.x))
			}else{lik.mat[x==ith.x]	<-	(psi.i*phi.g.list[[i]](z=0,a=a.i,b=b,x=xm1.i, gam=gam)/factorial(ith.x))}
			
	}
	
	lik.mat[lik.mat<=0]	<-	.Machine$double.xmin
	llik.mat			<-	log(lik.mat)
	
	return(-sum(llik.mat))
	
}



#'	Discrete Self Decomposable Generalized model
#'
#'	Estimates of the parameters of the Discrete Self Decomposable model and allows covariates in 
#'	the decay and zero inflation process
#'	The input should be in the following way:
#'	@param 	Alpha.mod		formula for covariates of alpha (~x.alpha)
#'	@param	psi.mod 	 	formula for covariates of psi (~x.psi)
#'	@param	tsdata 			time series data in the form ncol = time steps and nrows = replicates of the time series. It has to be a Matrix.
#'	@param	covariates		a data frame (Not a Matrix!) with the covaiates data for both alpha and psi. 
#'							It allows for any type of covariates including factors. If the covariates are dynamic (i.e. they 
#'							change from tiem step to time step) then enter them	as a list with each element of the list being 
#'							a single data.frame with the values of a single covariate for each quadrant and each time step. 
#'							Then, the length of the list should be equal to the number of covariates you want to use and each element 
#'							should have dimensions identical to ts.data.
#'	@param	method			method for optimization. You can use any of the methods available in optim.
#'							It defaults to "Nelder-mead"
#'	@param CI				If confidence interval should be calculated and if so, which type. The options are
#'							"no", "wald", "boot". "wald" computes wald's confidence intervals of the parameters based
#'							on th computation of the hessian matrix by the optim function. "boot" computes bootstrap
#'							confidence intervals of the parameters. Bear in mind that the option "boot" requires many
#'							optimizations of the process and thus it can be slow depending on the size of your dataset.
#'	@param inits			specification of initial values of the parameters for the optimization process. inits can be
#'							left NULL and in such case, the function initiates the optimization by obtaining random values of 
#'							the parameters. If inits is specified, they have to be in the following order: 
#'							c(betas.alpha,betas.psi,beta,gamma).
#'	@param boot.reps		In case that you want bootstrap confidence intervals you can specify the number of bootstrap replicates to
#'							perform. It can be left NULL and in such case 2000 replicates will be performed and a warning will be issued.

#'	@return					a list with the following elements: 
#'	@return					MLEs 	a vector with the values of The Maximum likelihood Estimates of the parameters in the following order
#'									the intercept and coefficients of the alpha parameter relationship with the covariates
#'									the intercept and coefficients of the psi parameter relationship with the covariates
#'									beta and gamma.
#'	@return					LogLik	The log likelihood value of the data given the parameters
#'	@return					Predicted.psi  predicted values of psi for each quadrant
#'	@return					Predicted.alpha	predicted values of alpha for each quadrant
#'	@return					Expected.val  Expected value of the times series in each quadrant
#'	@return					If CI = "wald" or "boot" it also returns the confidence intervals of each parameters and
#'							the CI of psi and alpha for each quadrant. 
#'	@export
#'	@examples
#'	quadrants <-	40
#'	t.steps		<-	10
#'	b <- 3
#'	gam	<-	1.7 
#'	x0	<-	2
#'	betas.a	<-	c(0.1,2)
#'	betas.psi	<-	c(2,0.9)

#'	# Covariates for alpha and psi
#'	ones	<-	rep(1,quadrants)
#'	set.seed(1234)
#'	x.alpha	<-	cbind(ones,rnorm(quadrants))
#'	alphas	<-	expit(x.alpha%*%betas.a)
#'	set.seed(4521)
#'	x.psi	<-	cbind(ones,rnorm(quadrants))
#'	psi		<-	expit(x.psi%*%betas.psi)

#'	covars		<-	data.frame(x.alpha=x.alpha[,2],x.psi=x.psi[,2])

#'	tsdata	<-	rdsd(t.steps=t.steps,t0=x0,alpha=alphas,beta=b,gamma=gam,psi=psi,reps=quadrants)	
#'	mod		<-	dsd.glm(~x.alpha,~x.psi,tsdata=tsdata,covariates=covars,CI="wald")
dsd.glm	<-	function(alpha.mod,psi.mod,data,method="Nelder-Mead",inits, pboot,family){
	
	if(class(data)!="list"&class(data)!="data.frame"){stop("data must be a list or a data.frame")}
	if(class(data)=="list"){tsdata	<-	data[[1]]}else{tsdata<-data}
	
	############ Setting up of the data and parameters for the optimization:
	# No Covariates in any parameter
	
			quads	<-	ncol(tsdata)
			tsteps	<-	nrow(tsdata)
			
			newdat	<-	lapply(data,as.vector)
		
			alpha.frame	<-	model.frame(alpha.mod,newdat)
			alpha.desmat<-	model.matrix(alpha.mod,alpha.frame)
			alpha.ncoefs<-	ncol(alpha.desmat)
				
			psi.frame	<-	model.frame(psi.mod,newdat)
			psi.desmat	<-	model.matrix(psi.mod,psi.frame)
			psi.ncoefs	<-	ncol(psi.desmat)
			
	if(missing(inits)){
		
		switch(family
			,"NB"={
					
					guess	<-	c(rep(0,alpha.ncoefs),rep(0,psi.ncoefs),0,0)
					}
			,"Poisson"={
					lambda.init		<-	mean(tsdata)
					guess	<-	c(rep(0,alpha.ncoefs),rep(0,psi.ncoefs),log(lambda.init))
				})
		}else{guess	<-	inits}
		
	############ Optimization of the dsd.like.cov function
	
	switch(family
	,"NB"={
		mles		<-	optim(guess,dsd.like.NB,method=method,tsdata=tsdata,alpha.desmat=alpha.desmat,psi.desmat=psi.desmat)	
		loglike		<-	-mles$value
		mod.coefs	<-	mles$par
		alpha.coefs	<-	mod.coefs[1:alpha.ncoefs]
		psi.coefs	<-	mod.coefs[(alpha.ncoefs+1):(alpha.ncoefs+psi.ncoefs)]
	
		b.hat		<-	exp(mod.coefs[(alpha.ncoefs+psi.ncoefs+1)])
		gamma.hat	<-	exp(mod.coefs[length(mod.coefs)])
		}
	,"Poisson"={
		mles		<-	optim(guess,dsd.like.Pois,method=method,tsdata=tsdata,alpha.desmat=alpha.desmat,psi.desmat=psi.desmat)	
		loglike		<-	-mles$value
		mod.coefs	<-	mles$par
		alpha.coefs	<-	mod.coefs[1:alpha.ncoefs]
		psi.coefs	<-	mod.coefs[(alpha.ncoefs+1):(alpha.ncoefs+psi.ncoefs)]
		lambda.hat	<-	exp(mod.coefs[(alpha.ncoefs+psi.ncoefs+1)])
	}
	)

	if(missing(pboot)){
		
		switch(family
			,"NB"={
				aic			<-	2*(alpha.ncoefs+psi.ncoefs+2)-2*loglike
				results		<-	list(loglike=loglike,AIC=aic,alpha=alpha.coefs,psi=psi.coefs,b=b.hat,gamma=gamma.hat)			
				}
			,"Poisson"={
				aic			<-	2*(alpha.ncoefs+psi.ncoefs+1)-2*loglike
				results	<-	list(loglike=loglike,AIC=aic,alpha=alpha.coefs,psi=psi.coefs,lambda=lambda.hat)}
		)

		return(results)
		
	}else{
		
		#mle alpha matrix (values of alpha for each cell in raw data matrix)
	
		# Alpha		
		a.cov		<-	ncol(alpha.desmat)
		betas.alpha	<-	mod.coefs[1:a.cov]
		a.hat		<-	alpha.desmat%*%betas.alpha
		a.hat		<-	matrix(plogis(a.hat),ncol=quads,nrow=tsteps)
	
		#mle psi matrix (values of psi for each cell in raw data matrix)
		
		# Psi
		psi.cov		<-	ncol(psi.desmat)
		betas.psi	<-	mod.coefs[(a.cov+1):(a.cov+psi.cov)]
		psi.hat		<-	psi.desmat%*%betas.psi
		psi.hat		<-	matrix(plogis(psi.hat),ncol=quads,nrow=tsteps)
		
		boot.reps  <- pboot
		boot.mat <- matrix(0,nrow=boot.reps,ncol=length(mod.coefs))
		pb				<-	txtProgressBar(min=0,max=boot.reps,style=3)
		
		switch(family
			,"NB"={
				
				aic			<-	2*(alpha.ncoefs+psi.ncoefs+2)-2*loglike
				
				for(i in 1:boot.reps){

					## Simulate data based on estimated parameters
					ts.mat			<-	rdsd(t.steps=tsteps,alpha=a.hat,beta=b.hat,gamma=gamma.hat,psi=psi.hat,reps=quads)
					ith.mles		<-	optim(mod.coefs,dsd.like.NB,method=method,tsdata=ts.mat,alpha.desmat=alpha.desmat,psi.desmat=psi.desmat)
					boot.mat[i,]<-	ith.mles$par
					setTxtProgressBar(pb,i)	
				}
				close(pb)
		
				par.ci	<-	apply(boot.mat,2,quantile,prob=c(0.025,0.975))
		
				alpha.coefs	<-	matrix(c(par.ci[1,1:alpha.ncoefs]
										,mod.coefs[1:alpha.ncoefs]
										,par.ci[2,1:alpha.ncoefs])
										,ncol=3,nrow=alpha.ncoefs
										,dimnames=list(colnames(alpha.desmat),c("2.5%","mle","97.5%")))
				psi.coefs	<-	matrix(c(par.ci[1,(alpha.ncoefs+1):(alpha.ncoefs+psi.ncoefs)]
										,mod.coefs[(alpha.ncoefs+1):(alpha.ncoefs+psi.ncoefs)]
										,par.ci[2,(alpha.ncoefs+1):(alpha.ncoefs+psi.ncoefs)])
										,ncol=3,nrow=psi.ncoefs
										,dimnames=list(colnames(psi.desmat),c("2.5%","mle","97.5%")))
				b.hat		<-	matrix(exp(c(par.ci[1,(alpha.ncoefs+psi.ncoefs+1)]
											,mod.coefs[(alpha.ncoefs+psi.ncoefs+1)]
											,par.ci[2,(alpha.ncoefs+psi.ncoefs+1)]))
											,ncol=3,nrow=1
											,dimnames=list("beta",c("2.5%","mle","97.5%")))
				gamma.hat	<-	matrix(exp(c(par.ci[1,length(mod.coefs)]
											,mod.coefs[length(mod.coefs)]
											,par.ci[2,length(mod.coefs)]))
											,ncol=3,nrow=1
											,dimnames=list("gamma",c("2.5%","mle","97.5%")))
				results		<-	list(loglike=loglike,AIC=aic,alpha=alpha.coefs,psi=psi.coefs,b=b.hat,gamma=gamma.hat)
									}
			,"Poisson"={
				
				aic			<-	2*(alpha.ncoefs+psi.ncoefs+1)-2*loglike
				
				for(i in 1:boot.reps){

					## Simulate data based on estimated parameters
					ts.mat			<-	rdsd.pois(t.steps=tsteps,alpha=a.hat,psi=psi.hat,lambda=lambda.hat,reps=quads)
					ith.mles		<-	optim(mod.coefs,dsd.like.Pois,method=method,tsdata=ts.mat,alpha.desmat=alpha.desmat,psi.desmat=psi.desmat)
					boot.mat[i,]<-	ith.mles$par
					setTxtProgressBar(pb,i)	
				}
				close(pb)
		
				par.ci	<-	apply(boot.mat,2,quantile,prob=c(0.025,0.975))
		
				alpha.coefs	<-	matrix(c(par.ci[1,1:alpha.ncoefs]
										,mod.coefs[1:alpha.ncoefs]
										,par.ci[2,1:alpha.ncoefs])
										,ncol=3,nrow=alpha.ncoefs
										,dimnames=list(colnames(alpha.desmat),c("2.5%","mle","97.5%")))
				psi.coefs	<-	matrix(c(par.ci[1,(alpha.ncoefs+1):(alpha.ncoefs+psi.ncoefs)]
										,mod.coefs[(alpha.ncoefs+1):(alpha.ncoefs+psi.ncoefs)]
										,par.ci[2,(alpha.ncoefs+1):(alpha.ncoefs+psi.ncoefs)])
										,ncol=3,nrow=psi.ncoefs
										,dimnames=list(colnames(psi.desmat),c("2.5%","mle","97.5%")))
				lambda.hat		<-	matrix(exp(c(par.ci[1,(alpha.ncoefs+psi.ncoefs+1)]
											,mod.coefs[(alpha.ncoefs+psi.ncoefs+1)]
											,par.ci[2,(alpha.ncoefs+psi.ncoefs+1)]))
											,ncol=3,nrow=1
											,dimnames=list("beta",c("2.5%","mle","97.5%")))
				results		<-	list(loglike=loglike,AIC=aic,alpha=alpha.coefs,psi=psi.coefs,lambda=lambda.hat)
			
				}
				
		)

		return(results)
		
	}
	
		
}

dsd.profglm	<-	function(alpha.mod,psi.mod,b.vec,gam.vec,data,method="Nelder-Mead",inits){
	
	if(class(data)!="list"&class(data)!="data.frame"){stop("data must be a list or a data.frame")}
	if(class(data)=="list"){tsdata	<-	data[[1]]}else{tsdata<-data}
	
	############ Setting up of the data and parameters for the optimization:
	# No Covariates in any parameter
	
			quads	<-	ncol(tsdata)
			tsteps	<-	nrow(tsdata)
			
			newdat	<-	lapply(data,as.vector)
		
			alpha.frame	<-	model.frame(alpha.mod,newdat)
			alpha.desmat<-	model.matrix(alpha.mod,alpha.frame)
			alpha.ncoefs<-	ncol(alpha.desmat)
				
			psi.frame	<-	model.frame(psi.mod,newdat)
			psi.desmat	<-	model.matrix(psi.mod,psi.frame)
			psi.ncoefs	<-	ncol(psi.desmat)		
			
	if(missing(inits)){guess	<-	c(rep(0,alpha.ncoefs),rep(0,psi.ncoefs))
		}else{guess	<-	inits}
		
	############ Optimization of the dsd.like.cov function
	
	lik.mat	<-	matrix(0,ncol=length(b.vec),nrow=length(gam.vec),dimnames=list(gam.vec,b.vec))
	for(i in 1:length(b.vec)){
		
		for(j in 1:length(gam.vec)){
			
				mles		<-	optim(guess,dsd.proflike,method=method,tsdata=tsdata,alpha.desmat=alpha.desmat,psi.desmat=psi.desmat,b=b.vec[i],gam=gam.vec[j])
				lik.mat[j,i]		<-	-mles$value			
						
		}
		
	}
	
	return(lik.mat)
}



dredge.dsd	<-	function(alpha.mod,psi.mod,dataset,family){
	
	# All possible model combinations for alpha
	alpha.terms	 <-	terms(alpha.mod)	
	alpha.names	 <-	attr(alpha.terms,"term.labels")
	n.terms  	 <- length(alpha.names)
	nmods    	 <- 2^n.terms
	alpha.vec	<-	rep("NA",length=nmods)
	alpha.vec[nmods] 	 <- paste(all.vars(alpha.terms)[1],"1",sep="~")

	if(nmods>1){
		alpha.vec[(nmods-1)] <- paste(all.vars(alpha.terms)[1],paste(alpha.names, collapse="+"),sep="~") 
		#alpha.models	<-	list()

		nsubcombs <- rep(0,length((n.terms-1):1))
		start.p  <- 0
		end.p    <- 0

		for(i in 1:(n.terms-1)){
	
			m <- i 
			combs <- t(combn(alpha.names,m))
			ith.ncombs <- nrow(combs)
			nsubcombs[i] <- ith.ncombs
	
			for(j in 1:ith.ncombs){
		
				start.p		<-	end.p + j		
				x.s			<-	paste(combs[j,],collapse="+")
				jmname      <- paste(all.vars(alpha.terms)[1],x.s,sep="~")
				ith.form	<-	formula(jmname)
				alpha.vec[start.p]	<-	jmname
				#alpha.models[[start.p]]	<-	dsd.glm(alpha.mod=ith.form,psi.mod=y~1,data=dataset)	
			}
			end.p <- sum(nsubcombs)
	
		}
	
	}
	
	# All possible model combinations for psi
	psi.terms	 <-	terms(psi.mod)	
	psi.names	 <-	attr(psi.terms,"term.labels")
	n.terms  	 <- length(psi.names)
	nmods    	 <- 2^n.terms
	psi.vec	<-	rep("NA",length=nmods)
	psi.vec[nmods] 	 <- paste(all.vars(psi.terms)[1],"1",sep="~")
	
	if(nmods>1){

		psi.vec[(nmods-1)] <- paste(all.vars(psi.terms)[1],paste(psi.names, collapse="+"),sep="~") 
		#psi.models	<-	list()

		nsubcombs <- rep(0,length((n.terms-1):1))
		start.p  <- 0
		end.p    <- 0

		for(i in 1:(n.terms-1)){
	
			m <- i 
			combs <- t(combn(psi.names,m))
			ith.ncombs <- nrow(combs)
			nsubcombs[i] <- ith.ncombs
	
			for(j in 1:ith.ncombs){
		
				start.p		<-	end.p + j		
				x.s			<-	paste(combs[j,],collapse="+")
				jmname      <- paste(all.vars(psi.terms)[1],x.s,sep="~")
				ith.form		<-	formula(jmname)
				psi.vec[start.p]			<-	jmname
				#psi.models[[start.p]]	<-	dsd.glm(alpha.mod=y~1,psi.mod=ith.form,data=dataset)	
			}
			end.p <- sum(nsubcombs)
	
		}
	}
	
	final.fmls	<-	data.frame(alpha=rep(alpha.vec,each=length(psi.vec)),psi=rep(psi.vec,length(alpha.vec)))
	
	psi.coef.mat	<-	matrix(NA,nrow=nrow(final.fmls),ncol=length(psi.names)+1,dimnames=list(1:nrow(final.fmls),c("Intercept",psi.names)))
	alpha.coef.mat	<-	matrix(NA,nrow=nrow(final.fmls),ncol=length(alpha.names)+1,dimnames=list(1:nrow(final.fmls),c("Intercept",alpha.names)))
	
	logLikes.vec <- AIC.vec <- beta.hats <- gamma.hats <- lambda.hats <- rep(0,nrow(final.fmls));
		
	
	for(i in 1:nrow(final.fmls)){
		
		fmla.alpha			<-	formula(as.character(final.fmls$alpha[i]))
		fmla.psi           <-	formula(as.character(final.fmls$psi[i]))
		ith.mod             <-  dsd.glm(alpha.mod=fmla.alpha,psi.mod=fmla.psi,data=dataset,family=family)
		
		logLikes.vec[i] <- ith.mod$loglike
		AIC.vec[i]      <- ith.mod$AIC
		
		if(family=="NB"){
			beta.hats[i]    <- ith.mod$b
			gamma.hats[i]   <- ith.mod$gamma  
		}else{lambda.hats[i]<- ith.mod$lambda}
		
		ith.alpha.terms		<-	terms(fmla.alpha)
		ith.alpha.vars		<-	attr(ith.alpha.terms,"term.labels")
		ith.alpha.coefs		<-	ith.mod$alpha
		alpha.coef.mat[i,1]	<-	ith.alpha.coefs[1]
		alpha.coef.mat[i,match(ith.alpha.vars,colnames(alpha.coef.mat))]	<-	ith.alpha.coefs[2:length(ith.alpha.coefs)]
		
		ith.psi.terms		<-	terms(fmla.psi)
		ith.psi.vars			<-	attr(ith.psi.terms,"term.labels")
		ith.psi.coefs		<-	ith.mod$psi
		psi.coef.mat[i,1]	<-	ith.psi.coefs[1]
		psi.coef.mat[i,match(ith.psi.vars,colnames(psi.coef.mat))]	<-	ith.psi.coefs[2:length(ith.psi.coefs)]

		

		
	}
	
	colnames(psi.coef.mat)	<-	paste(colnames(psi.coef.mat),"psi",sep="_")
	colnames(alpha.coef.mat)<-	paste(colnames(alpha.coef.mat),"alpha",sep="_")
	if(family=="NB"){
		final.res	<-	data.frame(final.fmls,alpha.coef.mat,psi.coef.mat,Beta=beta.hats,Gamma=gamma.hats,Loglik=logLikes.vec,AIC=AIC.vec)
	}else{final.res	<-	data.frame(final.fmls,alpha.coef.mat,psi.coef.mat,lambda=lambda.hats,Loglik=logLikes.vec,AIC=AIC.vec)}
	final.res	<-	final.res[order(final.res$AIC),]
	return(final.res)	
	
}

# Likelihood of the normal approximation of the DSD

dsd.cls	<-	function(guess,alpha.desmat,tsdata,tausq=FALSE){

	qp1     <- nrow(tsdata)
	q       <- qp1-1
	nreps	<-	ncol(tsdata)
	npars	<-	length(guess)
	
	a.cov		<-	ncol(alpha.desmat)
	betas.alpha	<-	guess[1:a.cov]
	alpha		<-	alpha.desmat%*%betas.alpha
	alpha		<-	matrix(plogis(alpha),ncol=nreps,nrow=qp1)
	
	if(tausq){
		beta	<- exp(guess[(npars-2)])
		gamma	<-	exp(guess[(npars-1)])	
		tausq	<- exp(guess[(npars)])
		}else{
		beta	<- exp(guess[(npars-1)])
		gamma	<- exp(guess[npars])
		tausq	<-	0
		}

	psi		<- 1
	
	mu0		<- gamma/beta
	sigsq0  <- mu0 + mu0/(beta)
	mu.t    <- mu0*(1-alpha)
	sigsq.t <- mu.t*(1+ (1+alpha)/beta )
	mu.Es   <- rbind(rep(mu0,nreps),mu.t[2:qp1,])
	ofn		<-	numeric(nreps)
	for(i in 1:nreps){
		Sigma   <- diag(c(sigsq0,sigsq.t[2:qp1,i]))	
		A0      <- toeplitz(c(1,cumprod(alpha[2:qp1,i])))
		A       <- A0*lower.tri(A0,diag=TRUE)
		At      <- t(A)
		mu.vec  <- A%*%matrix(mu.Es[,i],nrow=qp1,ncol=1) 
		ASAt    <- A%*%Sigma%*%At
		
		#Variance of observation error
		Itausq       <- diag(tausq,nrow=qp1,ncol=qp1);
		V            <- ASAt + Itausq;
		Vinv         <- ginv(V)
		
		SSQ     <- t(tsdata[,i] - mu.vec)%*%Vinv%*%(tsdata[,i]-mu.vec)
		ofn[i]  <- sum((qp1/2)*log(2*pi) + (0.5*log(det(V))) + 0.5*SSQ)
	}
	return(sum(ofn))
	
}


dsd.glm2	<-	function(alpha.mod,data,method="Nelder-Mead",tausq=FALSE,inits, pboot){
	
	if(class(data)!="list"&class(data)!="data.frame"){stop("data must be a list or a data.frame")}
	if(class(data)=="list"){tsdata	<-	data[[1]]}else{tsdata<-data}
	
	############ Setting up of the data and parameters for the optimization:
	# No Covariates in any parameter
	
	# First the model for alpha, beta, gamma and tausq.
			quads	<-	ncol(tsdata)
			tsteps	<-	nrow(tsdata)
			
			newdat	<-	lapply(data,as.vector)
		
			alpha.frame	<-	model.frame(alpha.mod,newdat)
			alpha.desmat<-	model.matrix(alpha.mod,alpha.frame)
			alpha.ncoefs<-	ncol(alpha.desmat)
				
	if(missing(inits)){
		if(tausq){guess	<-	c(rep(0,alpha.ncoefs),0,0,0)}else{
			guess	<-	c(rep(0,alpha.ncoefs),0,0)
			}
			}else{guess	<-	inits}
		
	############ Optimization of the dsd.like.cov function
		mles		<-	optim(guess,dsd.cls,method=method,tsdata=tsdata,alpha.desmat=alpha.desmat,tausq=tausq)	
		loglike		<-	-mles$value
		mod.coefs	<-	mles$par
		alpha.coefs	<-	mod.coefs[1:alpha.ncoefs]
		npars		<-	length(mod.coefs)
		
		if(tausq){	
			b.hat		<-	exp(mod.coefs[(npars-2)])
			gamma.hat	<-	exp(mod.coefs[(npars-1)])
			tausq.hat	<-	exp(mod.coefs[npars])

		}else{
			b.hat		<-	exp(mod.coefs[(npars-1)])
			gamma.hat	<-	exp(mod.coefs[(npars)])
			tausq.hat	<-	0
			}
			
		aic			<-	2*(npars)-2*loglike	
		
			
	if(missing(pboot)){
		

			results		<-	list(loglike=loglike,AIC=aic,alpha=alpha.coefs,beta=b.hat,gamma=gamma.hat,tausq=tausq.hat)			
			
		return(results)
		
	}else{
		
		#mle alpha matrix (values of alpha for each cell in raw data matrix)
	
		# Alpha		
		a.cov		<-	ncol(alpha.desmat)
		betas.alpha	<-	mod.coefs[1:a.cov]
		a.hat		<-	alpha.desmat%*%betas.alpha
		a.hat		<-	matrix(plogis(a.hat),ncol=quads,nrow=tsteps)
			
		boot.reps  <- pboot
		boot.mat <- matrix(0,nrow=boot.reps,ncol=length(mod.coefs))
		pb				<-	txtProgressBar(min=0,max=boot.reps,style=3)
		
		for(i in 1:boot.reps){

			## Simulate data based on estimated parameters
			ts.mat			<-	rdsd(t.steps=tsteps,alpha=a.hat,beta=b.hat,gamma=gamma.hat
									,psi=matrix(1,ncol=quads,nrow=tsteps),reps=quads)
			ith.mles		<-	optim(mod.coefs,dsd.cls,method=method,tsdata=ts.mat
										,alpha.desmat=alpha.desmat,tausq=FALSE)
			hats			<-	ith.mles$par
			alpha.boot		<-	hats[1:alpha.ncoefs]
			beta.boot		<-	exp(hats[npars-1])
			gamma.boot		<-	exp(hats[npars])
			
			boot.mat[i,]<-	c(alpha.boot,beta.boot,gamma.boot)
			setTxtProgressBar(pb,i)	
		}
		close(pb)
		
		par.ci	<-	apply(boot.mat,2,quantile,prob=c(0.025,0.975))
		
		alpha.coefs	<-	matrix(c(par.ci[1,1:alpha.ncoefs]
								,alpha.coefs
								,par.ci[2,1:alpha.ncoefs])
								,ncol=3,nrow=alpha.ncoefs
								,dimnames=list(colnames(alpha.desmat),c("2.5%","mle","97.5%")))
		b.hat		<-	matrix(exp(c(par.ci[1,(npars-1)]
									,b.hat
									,par.ci[2,(npars-1)]))
									,ncol=3,nrow=1
									,dimnames=list("beta",c("2.5%","mle","97.5%")))
		gamma.hat	<-	matrix(exp(c(par.ci[1,npars]
									,gamma.hat
									,par.ci[2,npars]))
									,ncol=3,nrow=1
									,dimnames=list("gamma",c("2.5%","mle","97.5%")))
		results		<-	list(loglike=loglike,AIC=aic,alpha=alpha.coefs,b=b.hat,gamma=gamma.hat)
	}

	return(results)
		
}
	


