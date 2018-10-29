# Analysis of the chamois data.
# Load dsd functions and other packages
#library(snow)
source('~/Dropbox/JuanPabloPostdoc/DSD_Code/LIZstoolkit_v2.1.R', chdir = TRUE)

# Load data
# For 2010
files2010	<-	list.files("~/Documents/UFL/RESEARCH/MathieuBasille/2010",full.names=TRUE)
files2010	<-	files2010[-c(3,7,8)]

F2010	<-	list(t(read.csv("~/Documents/UFL/RESEARCH/MathieuBasille/2010/F2010.csv",row.names=1)))
MT2010	<-	list(t(read.csv("~/Documents/UFL/RESEARCH/MathieuBasille/2010/MT2010.csv",row.names=1)))
MNT2010	<-	list(t(read.csv("~/Documents/UFL/RESEARCH/MathieuBasille/2010/MNT2010.csv",row.names=1)))

# For 2011
files2011	<-	list.files("~/Documents/UFL/RESEARCH/MathieuBasille/2011",full.names=TRUE)
files2011	<-	files2011[-c(3,7,8)]

F2011	<-	list(t(read.csv("~/Documents/UFL/RESEARCH/MathieuBasille/2011/F2011.csv",row.names=1)))
MT2011	<-	list(t(read.csv("~/Documents/UFL/RESEARCH/MathieuBasille/2011/MT2011.csv",row.names=1)))
MNT2011	<-	list(t(read.csv("~/Documents/UFL/RESEARCH/MathieuBasille/2011/MNT2011.csv",row.names=1)))

# For 2012
files2012	<-	list.files("~/Documents/UFL/RESEARCH/MathieuBasille/2012",full.names=TRUE)
files2012	<-	files2012[-c(3,7,8)]

F2012	<-	list(t(read.csv("~/Documents/UFL/RESEARCH/MathieuBasille/2012/F2012.csv",row.names=1)))
MT2012	<-	list(t(read.csv("~/Documents/UFL/RESEARCH/MathieuBasille/2012/MT2012.csv",row.names=1)))
MNT2012	<-	list(t(read.csv("~/Documents/UFL/RESEARCH/MathieuBasille/2012/MNT2012.csv",row.names=1)))


for(i in 1:length(files2010)){
	
	cov2010			<-	t(read.csv(files2010[i],row.names=1))
	cov2011			<-	t(read.csv(files2011[i],row.names=1))
	cov2012			<-	t(read.csv(files2012[i],row.names=1))
	F2010[[i+1]]	<-	MT2010[[i+1]]	<-	MNT2010[[i+1]]	<-	cov2010
	F2011[[i+1]]		<-	MT2011[[i+1]]	<-	MNT2011[[i+1]]	<-	cov2011
	F2012[[i+1]]		<-	MT2012[[i+1]]	<-	MNT2012[[i+1]]	<-	cov2012
	
}

names.dat		<-	c("y","elev","elevxsnow","kF","kMNT","kMT")
names(F2010) <-	names(F2011)	<-	names(F2012)	<-	names.dat
names(MT2010)	<-	names(MT2011)	<-	names(MT2012)	<-	names.dat
names(MNT2010)	<-	names(MNT2011)	<-	names(MNT2012)	<-	names.dat

# Formula for models
mod.form	<-	formula(y~elev+elevxsnow+kMT+kMNT)



# Some Models for 2010

F2010.null	<-	dsd.glm(alpha.mod=y~1,psi.mod=y~1,data=F2010, family="NB")
F2010.mod	<-	dsd.glm(alpha.mod=mod.form,psi.mod=mod.form,data=F2010,pboot=10, family="NB")
MT2010.null	<-	dsd.glm(alpha.mod=y~1,psi.mod=y~1,data=MT2010)
MT2010.mod	<-	dsd.glm(alpha.mod=mod.form,psi.mod=mod.form,data=MT2010,pboot=1000)
MNT2010.null<-	dsd.glm(alpha.mod=y~1,psi.mod=y~1,data=MNT2010)
MNT2010.mod	<-	dsd.glm(alpha.mod=mod.form,psi.mod=mod.form,data=MNT2010,pboot=1000)

# All models for the Females 2010 data: output is a table of coefficients, log likelihood and AIC scores
# for each model
Mod.SelFem2010 <- dredge.dsd(alpha.mod= mod.form, psi.mod =y~1, dataset=F2010, family="NB")
# Copy the best model, re-run it using parametric bootstrap to calculate the confidence intervals
Best.modFem2010 <- dredge.dsd(alpha.mod= formula(as.character(Mod.SelFem2010$alpha[1])), psi.mod =formula(as.character(Mod.SelFem2010$psi[1])), dataset=F2010, pboot=1000)


# Now do the same for the MT and the MNT data...
Mod.SelMT2010 <- dredge.dsd(alpha.mod= mod.form, psi.mod =mod.form, dataset=MT2010)
Mod.SelMNT2010 <- dredge.dsd(alpha.mod= mod.form, psi.mod =mod.form, dataset=MNT2010)


# Some Models for 2011
F2011.null	<-	dsd.glm(alpha.mod=y~1,psi.mod=y~1,data=F2011)
F2011.mod	<-	dsd.glm(alpha.mod=mod.form,psi.mod=mod.form,data=F2011,pboot=1000)
MT2011.null	<-	dsd.glm(alpha.mod=y~1,psi.mod=y~1,data=MT2011)
MT2011.mod	<-	dsd.glm(alpha.mod=mod.form,psi.mod=mod.form,data=MT2011,pboot=1000)
MNT2011.null<-	dsd.glm(alpha.mod=y~1,psi.mod=y~1,data=MNT2011)
MNT2011.mod	<-	dsd.glm(alpha.mod=mod.form,psi.mod=mod.form,data=MNT2011,pboot=1000)

# All models for the 2011 data: output is a table of coefficients, log likelihood and AIC scores
# for each model
Mod.SelFem2011 <- dredge.dsd(alpha.mod= mod.form, psi.mod =mod.form, dataset=F2011)
Mod.SelMT2011 <- dredge.dsd(alpha.mod= mod.form, psi.mod =mod.form, dataset=MT2011)
Mod.SelMNT2011 <- dredge.dsd(alpha.mod= mod.form, psi.mod =mod.form, dataset=MNT2011)


# Some Models for 2012
F2012.null	<-	dsd.glm(alpha.mod=y~1,psi.mod=y~1,data=F2012)
F2012.mod	<-	dsd.glm(alpha.mod=mod.form,psi.mod=mod.form,data=F2012,pboot=1000)
MT2012.null	<-	dsd.glm(alpha.mod=y~1,psi.mod=y~1,data=MT2012)
MT2012.mod	<-	dsd.glm(alpha.mod=mod.form,psi.mod=mod.form,data=MT2012,pboot=1000)
MNT2012.null<-	dsd.glm(alpha.mod=y~1,psi.mod=y~1,data=MNT2012)
MNT2012.mod	<-	dsd.glm(alpha.mod=mod.form,psi.mod=mod.form,data=MNT2012,pboot=1000)

# All models for the 2012 data: output is a table of coefficients, log likelihood and AIC scores
# for each model

Mod.SelFem2012 <- dredge.dsd(alpha.mod= mod.form, psi.mod =mod.form, dataset=F2012)
Mod.SelMT2012 <- dredge.dsd(alpha.mod= mod.form, psi.mod =mod.form, dataset=MT2012)
Mod.SelMNT2012 <- dredge.dsd(alpha.mod= mod.form, psi.mod =mod.form, dataset=MNT2012)

save.image("F-MT-MNT-TRIALRUN.RData")


### An example of :  simulating data and then estimating model parameters and from these getting the model predictions
quadrants <- 100
t.steps <- 50
b <- 3
gam <- 1.7
x0 <- 2
betas.a <- c(4.1,4)
betas.psi <- c(2,0.9)

# Covariate for alpha and psi
# We will assume there is only one covariate (and intercept)
# so the model for the probability alpha is
# 1/(1+exp(-intercept + slope*x)); where intercept is 0.1 and slope is 2
ones	<-	rep(1,quadrants*t.steps)
x.alpha	<-	cbind(ones,cov=rnorm(quadrants*t.steps)) 
x.acovar <- t(matrix(x.alpha[,2], nrow=quadrants,ncol=t.steps)) # this is the simulated matrix for the covariate of alpha
alphas	<-	t(matrix(expit(x.alpha%*%betas.a), nrow=quadrants,ncol=t.steps)) 

x.psi	<-	cbind(ones,rnorm(quadrants*t.steps))
x.pcovar <- t(matrix(x.psi[,2], nrow=quadrants, ncol=t.steps))
psi		<-	t(matrix(expit(x.psi%*%betas.psi), nrow=quadrants, ncol=t.steps))

tsdata	<-	rdsd(t.steps=t.steps,t0=x0,alpha=alphas,beta=b,gamma=gam,psi=psi,reps=quadrants)	
list.4test <- list(y = tsdata, acovar = x.acovar, psicovar = x.pcovar)
mod		<-	dsd.glm(alpha.mod=y~acovar,psi.mod=y~psicovar,data=list.4test, family="NB")
print(mod)

### Now that we got estimates of alpha and psi for each cell and time step, use these to compute model predictions 
# using function mu.x.nb

dsd.pred <- function(data,model,cov.alpha,cov.psi) { # two options for data format:
  # opt 1: data is a list with the first object being the data, the second the covariate for alpha, and the third the covariate for psi
  # opt 2: data is a data frame with the data, and cov.alpha and cov.psi are data frames with the covariates
  
  # Set up 
  if(class(data)!="list"&class(data)!="matrix"){stop("data must be a list or a matrix")}
  if(class(data)=="list"){tsdata	<-	data[[1]]}else{tsdata<-data}
  
  if(class(data)=="list"){cov.alpha	<-	data[[2]]}else{cov.alpha<-cov.alpha}
  if(class(data)=="list"){cov.psi	<-	data[[3]]}else{cov.psi<-cov.psi}
  
  quads	<-	ncol(tsdata)
  tsteps	<-	nrow(tsdata)

  # Handle covariates
  ones	<-	rep(1,quads*tsteps)
  
  x.alpha	<-	cbind(ones,cov=as.vector(cov.alpha))
  x.psi	<-	cbind(ones,cov=as.vector(cov.psi))
  
  # Compute matrix of alpha and psi
  alphas <- t(matrix(expit(x.alpha%*%model$alpha), nrow=ncol(data),ncol=nrow(data)))
  psis <- t(matrix(expit(x.psi%*%model$psi), nrow=ncol(data),ncol=nrow(data)))
  
  preds <- matrix(nrow=tsteps, ncol=quads)
  preds[1,] <- tsdata[1,]
  
  for (i in 2:tsteps) {
    for (j in 1:quads) {
      preds[i,j] <- mu.x.nb(a=alphas[i,j],b=mod$b,gam=mod$gamma,psi=psis[i,j],x=preds[i-1,j])
    }
    }

return(preds)

}

mod.pred <- dsd.pred(data=tsdata, model=mod, cov.alpha=x.acovar, cov.psi=x.pcovar)
