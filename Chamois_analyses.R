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
names(F2010)	<-	names(F2011)	<-	names(F2012)	<-	names.dat
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

