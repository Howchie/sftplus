pkgname <- "sftplus"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('sftplus')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("assessment")
### * assessment

flush(stderr()); flush(stdout())

### Name: assessment
### Title: Workload Assessment Functions
### Aliases: assessment
### Keywords: sft

### ** Examples

c1c.12 <- rexp(10000, .015)
c1i.12 <- rexp(10000, .01)
c1c    <- rexp(10000, .015)
c1i    <- rexp(10000, .01)

c2c.12 <- rexp(10000, .014)
c2i.12 <- rexp(10000, .01)
c2c    <- rexp(10000, .014)
c2i    <- rexp(10000, .01)

RT.1 <- pmin(c1c, c1i)
CR.1 <- c1c < c1i
RT.2 <- pmin(c2c, c2i)
CR.2 <- c2c < c2i

c1Correct <- c1c.12 < c1i.12
c2Correct <- c2c.12 < c2i.12

# OR Detection
CR.12 <- c1Correct | c2Correct
RT.12 <- rep(NA, 10000)
RT.12[c1Correct & c2Correct] <- pmin(c1c.12, c2c.12)[c1Correct & c2Correct]
RT.12[c1Correct & !c2Correct] <- c1c.12[c1Correct & !c2Correct]
RT.12[!c1Correct & c2Correct] <- c2c.12[!c1Correct & c2Correct]
RT.12[!c1Correct & !c2Correct] <- pmax(c1i.12, c2i.12)[!c1Correct & !c2Correct]

RT <- list(RT.12, RT.1, RT.2)
CR <- list(CR.12, CR.1, CR.2)
a.or.cf <- assessment(RT, CR, stopping.rule="OR", correct=TRUE, fast=TRUE, detection=TRUE)
a.or.cs <- assessment(RT, CR, stopping.rule="OR", correct=TRUE, fast=FALSE, detection=TRUE)
a.or.if <- assessment(RT, CR, stopping.rule="OR", correct=FALSE, fast=TRUE, detection=TRUE)
a.or.is <- assessment(RT, CR, stopping.rule="OR", correct=FALSE, fast=FALSE, detection=TRUE)

par(mfrow=c(2,2))
plot(a.or.cf, ylim=c(0,2))
plot(a.or.cs, ylim=c(0,2))
plot(a.or.if, ylim=c(0,2))
plot(a.or.is, ylim=c(0,2))


# AND 
CR.12 <- c1Correct & c2Correct
RT.12 <- rep(NA, 10000)
RT.12[CR.12] <- pmax(c1c.12, c2c.12)[CR.12]
RT.12[c1Correct & !c2Correct] <- c2i.12[c1Correct & !c2Correct]
RT.12[!c1Correct & c2Correct] <- c1i.12[!c1Correct & c2Correct]
RT.12[!c1Correct & !c2Correct] <- pmin(c1i.12, c2i.12)[!c1Correct & !c2Correct]

RT <- list(RT.12, RT.1, RT.2)
CR <- list(CR.12, CR.1, CR.2)
a.and.cf <- assessment(RT, CR, stopping.rule="AND", correct=TRUE, fast=TRUE, detection=TRUE)
a.and.cs <- assessment(RT, CR, stopping.rule="AND", correct=TRUE, fast=FALSE, detection=TRUE)
a.and.if <- assessment(RT, CR, stopping.rule="AND", correct=FALSE, fast=TRUE, detection=TRUE)
a.and.is <- assessment(RT, CR, stopping.rule="AND", correct=FALSE, fast=FALSE, detection=TRUE)

par(mfrow=c(2,2))
plot(a.and.cf, ylim=c(0,2))
plot(a.and.cs, ylim=c(0,2))
plot(a.and.if, ylim=c(0,2))
plot(a.and.is, ylim=c(0,2))



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("assessmentGroup")
### * assessmentGroup

flush(stderr()); flush(stdout())

### Name: assessmentGroup
### Title: Assessment Functions
### Aliases: assessmentGroup
### Keywords: sft

### ** Examples

## Not run: 
##D data(dots)
##D assessmentGroup(subset(dots, Condition=="OR"), 
##D   stopping.rule="OR", correct=TRUE, fast=FALSE, 
##D   detection=TRUE)
##D assessmentGroup(subset(dots, Condition=="AND"), 
##D   stopping.rule="AND", correct=TRUE, fast=TRUE, )
## End(Not run)



cleanEx()
nameEx("capacity.and")
### * capacity.and

flush(stderr()); flush(stdout())

### Name: capacity.and
### Title: Capacity Coefficient for Exhaustive (AND) Processing
### Aliases: capacity.and
### Keywords: sft

### ** Examples

rate1 <- .35
rate2 <- .3
RT.pa <- rexp(100, rate1)
RT.ap <- rexp(100, rate2)
RT.pp.limited <- pmax( rexp(100, .5*rate1), rexp(100, .5*rate2))
RT.pp.unlimited <- pmax( rexp(100, rate1), rexp(100, rate2))
RT.pp.super <- pmax( rexp(100, 2*rate1), rexp(100, 2*rate2))
tvec <- sort(unique(c(RT.pa, RT.ap, RT.pp.limited, RT.pp.unlimited, RT.pp.super)))

cap.limited <- capacity.and(RT=list(RT.pp.limited, RT.pa, RT.ap))
print(cap.limited$Ctest)
cap.unlimited <- capacity.and(RT=list(RT.pp.unlimited, RT.pa, RT.ap))
cap.super <- capacity.and(RT=list(RT.pp.super, RT.pa, RT.ap))

matplot(tvec, cbind(cap.limited$Ct(tvec), cap.unlimited$Ct(tvec), cap.super$Ct(tvec)),
  type='l', lty=1, ylim=c(0,3), col=2:4, main="Example Capacity Functions", xlab="Time", 
  ylab="C(t)")
abline(1,0)
legend('topright', c("Limited", "Unlimited", "Super"), lty=1, col=2:4)




cleanEx()
nameEx("capacity.id")
### * capacity.id

flush(stderr()); flush(stdout())

### Name: capacity.id
### Title: Capacity Coefficient for Full Identification (ID) Exhaustive
###   Processing
### Aliases: capacity.id
### Keywords: sft

### ** Examples

rate1p <- .35
rate1a <- .25
rate2p <- .35
rate2a <- .25
RT.pa <- pmax(rexp(100, rate1p), rexp(100, rate2a))
RT.ap <- pmax(rexp(100, rate2p), rexp(100, rate1a))
RT.nt <- pmax(rexp(100, rate1a), rexp(100, rate2a))
RT.pp.limited <- pmax( rexp(100, .5*rate1p), rexp(100, .5*rate2p))
RT.pp.unlimited <- pmax( rexp(100, rate1p), rexp(100, rate2p))
RT.pp.super <- pmax( rexp(100, 2*rate1p), rexp(100, 2*rate2p))
tvec <- sort(unique(c(RT.pa, RT.ap, RT.pp.limited, RT.pp.unlimited, RT.pp.super)))

cap.limited <- capacity.id(dt.rt=RT.pp.limited, nt.rt=RT.nt, 
			   st.rts=list(RT.pa, RT.ap))
print(cap.limited$Ctest)
cap.unlimited <- capacity.id(dt.rt=RT.pp.unlimited, nt.rt=RT.nt, 
			     st.rts=list(RT.pa, RT.ap))
cap.super <- capacity.id(dt.rt=RT.pp.super, nt.rt=RT.nt, 
			 st.rts=list(RT.pa, RT.ap))

matplot(tvec, cbind(cap.limited$Ct(tvec), cap.unlimited$Ct(tvec), cap.super$Ct(tvec)),
  type='l', lty=1, ylim=c(0,3), col=2:4, main="Example Capacity Functions", xlab="Time", 
  ylab="C(t)")
abline(1,0)
legend('topright', c("Limited", "Unlimited", "Super"), lty=1, col=2:4)




cleanEx()
nameEx("capacity.or")
### * capacity.or

flush(stderr()); flush(stdout())

### Name: capacity.or
### Title: Capacity Coefficient for First-Terminating (OR) Processing
### Aliases: capacity.or
### Keywords: sft

### ** Examples

rate1 <- .35
rate2 <- .3
RT.pa <- rexp(100, rate1)
RT.ap <- rexp(100, rate2)
RT.pp.limited <- pmin( rexp(100, .5*rate1), rexp(100, .5*rate2))
RT.pp.unlimited <- pmin( rexp(100, rate1), rexp(100, rate2))
RT.pp.super <- pmin( rexp(100, 2*rate1), rexp(100, 2*rate2))
tvec <- sort(unique(c(RT.pa, RT.ap, RT.pp.limited, RT.pp.unlimited, RT.pp.super)))

cap.limited <- capacity.or(RT=list(RT.pp.limited, RT.pa, RT.ap))
print(cap.limited$Ctest)
cap.unlimited <- capacity.or(RT=list(RT.pp.unlimited, RT.pa, RT.ap))
cap.super <- capacity.or(list(RT=RT.pp.super, RT.pa, RT.ap))

matplot(tvec, cbind(cap.limited$Ct(tvec), cap.unlimited$Ct(tvec), cap.super$Ct(tvec)),
  type='l', lty=1, ylim=c(0,3), col=2:4, main="Example Capacity Functions", xlab="Time", 
  ylab="C(t)")
abline(1,0)
legend('topright', c("Limited", "Unlimited", "Super"), lty=1, col=2:4)



cleanEx()
nameEx("capacity.stst")
### * capacity.stst

flush(stderr()); flush(stdout())

### Name: capacity.stst
### Title: Capacity Coefficient for Single-Target Self-Terminating (STST)
###   Processing
### Aliases: capacity.stst
### Keywords: sft

### ** Examples

rate1 <- .35
RT.pa <- rexp(100, rate1)
RT.pp.limited <- rexp(100, .5*rate1)
RT.pp.unlimited <- rexp(100, rate1)
RT.pp.super <- rexp(100, 2*rate1)
tvec <- sort(unique(c(RT.pa, RT.pp.limited, RT.pp.unlimited, RT.pp.super)))

cap.limited <- capacity.stst(RT=list(RT.pp.limited, RT.pa))
print(cap.limited$Ctest)
cap.unlimited <- capacity.stst(RT=list(RT.pp.unlimited, RT.pa))
cap.super <- capacity.stst(RT=list(RT.pp.super, RT.pa))

matplot(tvec, cbind(cap.limited$Ct(tvec), cap.unlimited$Ct(tvec), cap.super$Ct(tvec)),
  type='l', lty=1, ylim=c(0,5), col=2:4, main="Example Capacity Functions", xlab="Time", 
  ylab="C(t)")
abline(1,0)
legend('topright', c("Limited", "Unlimited", "Super"), lty=1, col=2:4, bty="n")




cleanEx()
nameEx("capacityGroup")
### * capacityGroup

flush(stderr()); flush(stdout())

### Name: capacityGroup
### Title: Capacity Analysis
### Aliases: capacityGroup
### Keywords: sft

### ** Examples

## Not run: 
##D data(dots)
##D capacityGroup(subset(dots, Condition=="OR"), 
##D   stopping.rule="OR")
##D capacityGroup(subset(dots, Condition=="AND"), 
##D   stopping.rule="AND")
## End(Not run)



cleanEx()
nameEx("dots")
### * dots

flush(stderr()); flush(stdout())

### Name: dots
### Title: RT and and Accuracy from a Simple Detection Task
### Aliases: dots
### Keywords: datasets

### ** Examples

data(dots)
summary(dots)
## Not run: 
##D sicGroup(dots)
##D capacityGroup(dots)
## End(Not run)



cleanEx()
nameEx("estimate.bounds")
### * estimate.bounds

flush(stderr()); flush(stdout())

### Name: estimate.bounds
### Title: Bounds on Response Time Cumulative Distribution Functions for
###   Parallel Processing Models
### Aliases: estimate.bounds
### Keywords: sft

### ** Examples

#randomly generated data
rate1 <- .35
rate2 <- .3
rate3 <- .4
RT.paa <- rexp(100, rate1)
RT.apa <- rexp(100, rate2)
RT.aap <- rexp(100, rate3)
RT.or <- pmin(rexp(100, rate1), rexp(100, rate2), rexp(100, rate3))
RT.and <- pmax(rexp(100, rate1), rexp(100, rate2), rexp(100, rate3))
tvec <- sort(unique(c(RT.paa, RT.apa, RT.aap, RT.or, RT.and)))

or.bounds <- estimate.bounds(RT=list(RT.paa, RT.apa, RT.aap), CR=NULL, assume.ID=FALSE, 
  unified.space=FALSE)
and.bounds <- estimate.bounds(RT=list(RT.paa, RT.apa, RT.aap))

## Not run: 
##D #plot the or bounds together with a parallel OR model
##D matplot(tvec, 
##D   cbind(or.bounds$Upper.Bound(tvec), or.bounds$Lower.Bound(tvec), ecdf(RT.or)(tvec)),
##D   type='l', lty=1, ylim=c(0,1), col=2:4, main="Example OR Bounds", xlab="Time", 
##D   ylab="P(T<t)")
##D abline(1,0)
##D legend('topright', c("Upper Bound", "Lower Bound", "Parallel OR Model"), 
##D   lty=1, col=2:4, bty="n")
##D 
##D #using the dots data set in sft package
##D data(dots)
##D attach(dots)
##D RT.A <- dots[Subject=='S1' & Condition=='OR' & Channel1==2 & Channel2==0, 'RT']
##D RT.B <- dots[Subject=='S1' & Condition=='OR' & Channel1==0 & Channel2==2, 'RT']
##D RT.AB <- dots[Subject=='S1' & Condition=='OR' & Channel1==2 & Channel2==2, 'RT']
##D tvec <- sort(unique(c(RT.A, RT.B, RT.AB)))
##D Cor.A <- dots[Subject=='S1' & Condition=='OR' & Channel1==2 & Channel2==0, 'Correct']
##D Cor.B <- dots[Subject=='S1' & Condition=='OR' & Channel1==0 & Channel2==2, 'Correct']
##D Cor.AB <- dots[Subject=='S1' & Condition=='OR' & Channel1==2 & Channel2==2, 'Correct']
##D capacity <- capacity.or(list(RT.AB,RT.A,RT.B), list(Cor.AB,Cor.A,Cor.B), ratio=TRUE)
##D bounds <- estimate.bounds(list(RT.A,RT.B), list(Cor.A,Cor.B), unified.space=TRUE)
##D 
##D #plot unified capacity coefficient space
##D plot(tvec, capacity$Ct(tvec), type="l", lty=1, col="red", lwd=2)
##D lines(tvec, bounds$Upper.Bound(tvec), lty=2, col="blue", lwd=2)
##D lines(tvec, bounds$Lower.Bound(tvec), lty=4, col="blue", lwd=2)
##D abline(h=1, col="black", lty=1)
## End(Not run)



cleanEx()
nameEx("estimateNAH")
### * estimateNAH

flush(stderr()); flush(stdout())

### Name: estimateNAH
### Title: Neslon-Aalen Estimator of the Cumulative Hazard Function
### Aliases: estimateNAH
### Keywords: survival sft

### ** Examples

x <- rexp(50, rate=.5)
censoring <- runif(50) < .90
H.NA <- estimateNAH(x, censoring)

# Plot the estimated cumulative hazard function
plot(H.NA$H, 
  main="Cumulative Hazard Function\n X ~ Exp(.5)    n=50", 
  xlab="X", ylab="H(x)")

# Plot 95% Confidence intervals
times <- seq(0,10, length.out=100)
lines(times, H.NA$H(times) + sqrt(H.NA$Var(times))*qnorm(1-.05/2), lty=2)
lines(times, H.NA$H(times) - sqrt(H.NA$Var(times))*qnorm(1-.05/2), lty=2)

# Plot the true cumulative hazard function
abline(0,.5, col='red')



cleanEx()
nameEx("estimateNAK")
### * estimateNAK

flush(stderr()); flush(stdout())

### Name: estimateNAK
### Title: Neslon-Aalen Estimator of the Reverse Cumulative Hazard Function
### Aliases: estimateNAK
### Keywords: survival sft

### ** Examples

x <- rexp(50, rate=.5)
censoring <- runif(50) < .90
K.NA <- estimateNAK(x, censoring)

# Plot the estimated cumulative reverse hazard function
plot(K.NA$K, 
  main="Cumulative Reverse Hazard Function\n X ~ Exp(.5)    n=50", 
  xlab="X", ylab="K(x)")

# Plot 95% Confidence intervals
times <- seq(0,10, length.out=100)
lines(times, K.NA$K(times) + sqrt(K.NA$Var(times))*qnorm(1-.05/2), lty=2)
lines(times, K.NA$K(times) - sqrt(K.NA$Var(times))*qnorm(1-.05/2), lty=2)

# Plot the true cumulative reverse hazard function
lines(times, log(pexp(times, .5)), col='red')



cleanEx()
nameEx("estimateUCIPand")
### * estimateUCIPand

flush(stderr()); flush(stdout())

### Name: estimateUCIPand
### Title: UCIP Performance on AND Tasks
### Aliases: estimateUCIPand
### Keywords: sft

### ** Examples

# Channel completion times and accuracy
rt1 <- rexp(100, rate=.5)
cr1 <- runif(100) < .90
rt2 <- rexp(100, rate=.4)
cr2 <- runif(100) < .95
Kucip = estimateUCIPand(list(rt1, rt2), list(cr1, cr2))


# Plot the estimated UCIP cumulative reverse hazard function
plot(Kucip$K, do.p=FALSE, 
  main="Estimated UCIP Cumulative Reverse Hazard Function\n
    X~max(X1,X2)    X1~Exp(.5)    X2~Exp(.4)", 
  xlab="X", ylab="K_UCIP(x)")
# Plot 95% Confidence intervals
times <- seq(0,10, length.out=100)
lines(times, Kucip$K(times) + sqrt(Kucip$Var(times))*qnorm(1-.05/2), lty=2)
lines(times, Kucip$K(times) - sqrt(Kucip$Var(times))*qnorm(1-.05/2), lty=2)
# Plot true UCIP cumulative reverse hazard function
lines(times[-1], log(pexp(times[-1], .5)) + log(pexp(times[-1], .4)), col='red')




cleanEx()
nameEx("estimateUCIPor")
### * estimateUCIPor

flush(stderr()); flush(stdout())

### Name: estimateUCIPor
### Title: UCIP Performance on OR Tasks
### Aliases: estimateUCIPor
### Keywords: sft

### ** Examples

# Channel completion times and accuracy
rt1 <- rexp(100, rate=.5)
cr1 <- runif(100) < .90
rt2 <- rexp(100, rate=.4)
cr2 <- runif(100) < .95
Hucip = estimateUCIPor(list(rt1, rt2), list(cr1, cr2))


# Plot the estimated UCIP cumulative hazard function
plot(Hucip$H, do.p=FALSE, 
  main="Estimated UCIP Cumulative Hazard Function\n
    X~min(X1,X2)    X1~Exp(.5)    X2~Exp(.4)", 
  xlab="X", ylab="H_UCIP(t)")
# Plot 95% Confidence intervals
times <- seq(0,10, length.out=100)
lines(times, Hucip$H(times) + sqrt(Hucip$Var(times))*qnorm(1-.05/2), lty=2)
lines(times, Hucip$H(times) - sqrt(Hucip$Var(times))*qnorm(1-.05/2), lty=2)
#Plot true UCIP cumulative hazard function
abline(0,.9, col='red')



cleanEx()
nameEx("fPCAassessment")
### * fPCAassessment

flush(stderr()); flush(stdout())

### Name: fPCAassessment
### Title: Functional Principal Components Analysis for the Assessment
###   Functions
### Aliases: fPCAassessment
### Keywords: sft

### ** Examples

## Not run: 
##D data(dots)
##D fPCAassessment(dots, dimensions=2, stopping.rule="OR", register="median",
##D                correct=TRUE, fast=FALSE, detection=TRUE, plotPCs=TRUE)
## End(Not run)



cleanEx()
nameEx("fPCAcapacity")
### * fPCAcapacity

flush(stderr()); flush(stdout())

### Name: fPCAcapacity
### Title: Functional Principal Components Analysis for the Capacity
###   Coefficient
### Aliases: fPCAcapacity
### Keywords: sft

### ** Examples

## Not run: 
##D data(dots)
##D fPCAcapacity(dots, dimensions=2,stopping.rule="OR", 
##D   plotPCs=TRUE)
## End(Not run)



cleanEx()
nameEx("mic_test")
### * mic_test

flush(stderr()); flush(stdout())

### Name: mic.test
### Title: Test of the Mean Interaction Contrast
### Aliases: mic.test
### Keywords: sft

### ** Examples

T1.h <- rweibull(300, shape=2 , scale=400 )
T1.l <- rweibull(300, shape=2 , scale=800 )
T2.h <- rweibull(300, shape=2 , scale=400 )
T2.l <- rweibull(300, shape=2 , scale=800 )

Serial.hh <- T1.h + T2.h
Serial.hl <- T1.h + T2.l
Serial.lh <- T1.l + T2.h
Serial.ll <- T1.l + T2.l
mic.test(HH=Serial.hh, HL=Serial.hl, LH=Serial.lh, LL=Serial.ll)

Parallel.hh <- pmax(T1.h, T2.h)
Parallel.hl <- pmax(T1.h, T2.l)
Parallel.lh <- pmax(T1.l, T2.h)
Parallel.ll <- pmax(T1.l, T2.l)
mic.test(HH=Parallel.hh, HL=Parallel.hl, LH=Parallel.lh, LL=Parallel.ll, method="art")




cleanEx()
nameEx("semiparametricSFT.bayes")
### * semiparametricSFT.bayes

flush(stderr()); flush(stdout())

### Name: semiparametricSFT.bayes
### Title: OR-centred semiparametric hierarchical Bayesian SFT model
### Aliases: semiparametricSFT.bayes

### ** Examples

## Not run: 
##D   fit <- semiparametricSFT.bayes(sftData, Condition = "OR", n_bins = 20, seed = 801)
##D   head(fit$tidy_curves)
##D 
##D   # Full double-factorial fit with hierarchical SIC and MIC.
##D   dfp <- semiparametricSFT.bayes(sftData, Condition = "OR", salience_split = TRUE, seed = 801)
##D   mic.bayes(dfp)$summary   # group and subject MIC with sign tests
##D   sic.bayes(dfp)$summary   # SIC(t) at every level
## End(Not run)



cleanEx()
nameEx("sft_data_to_rt")
### * sft_data_to_rt

flush(stderr()); flush(stdout())

### Name: sft_data_to_rt
### Title: Convert row-wise SFT data to RT/CR list input
### Aliases: sft_data_to_rt
### Keywords: sft

### ** Examples

## Not run: 
##D   lists <- sft_data_to_rt(sftData, Condition = "OR", Subject = "P1")
##D   ucip.test(lists$RT, lists$CR, stopping.rule = "OR")
##D   capacityGroup.bayes(sftData, Condition = "OR")
## End(Not run)



cleanEx()
nameEx("siDominance")
### * siDominance

flush(stderr()); flush(stdout())

### Name: siDominance
### Title: Dominance Test for Selective Influence
### Aliases: siDominance
### Keywords: sft

### ** Examples

T1.h <- rexp(50, .2)
T1.l <- rexp(50, .1)
T2.h <- rexp(50, .21)
T2.l <- rexp(50, .11)

HH <- T1.h + T2.h
HL <- T1.h + T2.l
LH <- T1.l + T2.h
LL <- T1.l + T2.l
siDominance(HH, HL, LH, LL)



cleanEx()
nameEx("sic")
### * sic

flush(stderr()); flush(stdout())

### Name: sic
### Title: Calculate the Survivor Interaction Contrast
### Aliases: sic
### Keywords: sft

### ** Examples

T1.h <- rexp(50, .2)
T1.l <- rexp(50, .1)
T2.h <- rexp(50, .21)
T2.l <- rexp(50, .11)

SerialAND.hh <- T1.h + T2.h
SerialAND.hl <- T1.h + T2.l
SerialAND.lh <- T1.l + T2.h
SerialAND.ll <- T1.l + T2.l
SerialAND.sic <- sic(HH=SerialAND.hh, HL=SerialAND.hl, LH=SerialAND.lh, 
  LL=SerialAND.ll)
print(SerialAND.sic$Dvals)
plot(SerialAND.sic$SIC, do.p=FALSE, ylim=c(-1,1))

p1 <- runif(200) < .3
SerialOR.hh <- p1[1:50]    * T1.h + (1-p1[1:50]   )*T2.h
SerialOR.hl <- p1[51:100]  * T1.h + (1-p1[51:100] )*T2.l
SerialOR.lh <- p1[101:150] * T1.l + (1-p1[101:150])*T2.h
SerialOR.ll <- p1[151:200] * T1.l + (1-p1[151:200])*T2.l
SerialOR.sic <- sic(HH=SerialOR.hh, HL=SerialOR.hl, LH=SerialOR.lh, LL=SerialOR.ll)
print(SerialOR.sic$Dvals)
plot(SerialOR.sic$SIC, do.p=FALSE, ylim=c(-1,1))

ParallelAND.hh <- pmax(T1.h, T2.h)
ParallelAND.hl <- pmax(T1.h, T2.l)
ParallelAND.lh <- pmax(T1.l, T2.h)
ParallelAND.ll <- pmax(T1.l, T2.l)
ParallelAND.sic <- sic(HH=ParallelAND.hh, HL=ParallelAND.hl, LH=ParallelAND.lh, 
  LL=ParallelAND.ll)
print(ParallelAND.sic$Dvals)
plot(ParallelAND.sic$SIC, do.p=FALSE, ylim=c(-1,1))

ParallelOR.hh <- pmin(T1.h, T2.h)
ParallelOR.hl <- pmin(T1.h, T2.l)
ParallelOR.lh <- pmin(T1.l, T2.h)
ParallelOR.ll <- pmin(T1.l, T2.l)
ParallelOR.sic <- sic(HH=ParallelOR.hh, HL=ParallelOR.hl, LH=ParallelOR.lh, 
  LL=ParallelOR.ll)
print(ParallelOR.sic$Dvals)
plot(ParallelOR.sic$SIC, do.p=FALSE, ylim=c(-1,1))



cleanEx()
nameEx("sicGroup")
### * sicGroup

flush(stderr()); flush(stdout())

### Name: sicGroup
### Title: SIC Analysis for a Group
### Aliases: sicGroup
### Keywords: ~sft

### ** Examples

## Not run: 
##D data(dots)
##D sicGroup(dots)
## End(Not run)



cleanEx()
nameEx("sic_test")
### * sic_test

flush(stderr()); flush(stdout())

### Name: sic.test
### Title: Statistical test of the SIC.
### Aliases: sic.test
### Keywords: sft

### ** Examples

T1.h <- rexp(50, .2)
T1.l <- rexp(50, .1)
T2.h <- rexp(50, .21)
T2.l <- rexp(50, .11)

SerialAND.hh <- T1.h + T2.h
SerialAND.hl <- T1.h + T2.l
SerialAND.lh <- T1.l + T2.h
SerialAND.ll <- T1.l + T2.l
sic.test(HH=SerialAND.hh, HL=SerialAND.hl, LH=SerialAND.lh, LL=SerialAND.ll)

p1 <- runif(200) < .3
SerialOR.hh <- p1[1:50]    * T1.h + (1-p1[1:50]   )*T2.h
SerialOR.hl <- p1[51:100]  * T1.h + (1-p1[51:100] )*T2.l
SerialOR.lh <- p1[101:150] * T1.l + (1-p1[101:150])*T2.h
SerialOR.ll <- p1[151:200] * T1.l + (1-p1[151:200])*T2.l
sic.test(HH=SerialOR.hh, HL=SerialOR.hl, LH=SerialOR.lh, LL=SerialOR.ll)

ParallelAND.hh <- pmax(T1.h, T2.h)
ParallelAND.hl <- pmax(T1.h, T2.l)
ParallelAND.lh <- pmax(T1.l, T2.h)
ParallelAND.ll <- pmax(T1.l, T2.l)
sic.test(HH=ParallelAND.hh, HL=ParallelAND.hl, LH=ParallelAND.lh, LL=ParallelAND.ll)

ParallelOR.hh <- pmin(T1.h, T2.h)
ParallelOR.hl <- pmin(T1.h, T2.l)
ParallelOR.lh <- pmin(T1.l, T2.h)
ParallelOR.ll <- pmin(T1.l, T2.l)
sic.test(HH=ParallelOR.hh, HL=ParallelOR.hl, LH=ParallelOR.lh, LL=ParallelOR.ll)



cleanEx()
nameEx("ucip_id_test")
### * ucip_id_test

flush(stderr()); flush(stdout())

### Name: ucip.id.test
### Title: A Statistical Test for Super or Limited Capacity
### Aliases: ucip.id.test
### Keywords: ~sft

### ** Examples

rate1p <- .35
rate1a <- .25
rate2p <- .35
rate2a <- .25
RT.pa <- pmax(rexp(100, rate1p), rexp(100, rate2a))
RT.ap <- pmax(rexp(100, rate2p), rexp(100, rate1a))
RT.nt <- pmax(rexp(100, rate1a), rexp(100, rate2a))
RT.pp.limited <- pmax( rexp(100, .5*rate1p), rexp(100, .5*rate2p))
RT.pp.unlimited <- pmax( rexp(100, rate1p), rexp(100, rate2p))
RT.pp.super <- pmax( rexp(100, 2*rate1p), rexp(100, 2*rate2p))

z.limited   <- ucip.id.test(dt.rt=RT.pp.limited,   nt.rt=RT.nt, st.rts=list(RT.pa, RT.ap)) 
z.unlimited <- ucip.id.test(dt.rt=RT.pp.unlimited, nt.rt=RT.nt, st.rts=list(RT.pa, RT.ap)) 
z.super     <- ucip.id.test(dt.rt=RT.pp.super,     nt.rt=RT.nt, st.rts=list(RT.pa, RT.ap)) 



cleanEx()
nameEx("ucip_test")
### * ucip_test

flush(stderr()); flush(stdout())

### Name: ucip.test
### Title: A Statistical Test for Super or Limited Capacity
### Aliases: ucip.test
### Keywords: ~sft

### ** Examples

rate1 <- .35
rate2 <- .3
RT.pa <- rexp(100, rate1)
RT.ap <- rexp(100, rate2)

CR.pa <- runif(100) < .98
CR.ap <- runif(100) < .98
CR.pp <- runif(100) < .96
CRlist <- list(CR.pp, CR.pa, CR.ap)

#  OR Processing
RT.pp.limited <- pmin( rexp(100, .5*rate1), rexp(100, .5*rate2))
RT.pp.unlimited <- pmin( rexp(100, rate1), rexp(100, rate2))
RT.pp.super <- pmin( rexp(100, 2*rate1), rexp(100, 2*rate2))
z.limited   <- ucip.test(RT=list(RT.pp.limited, RT.pa, RT.ap), CR=CRlist, stopping.rule="OR")
z.unlimited <- ucip.test(RT=list(RT.pp.unlimited, RT.pa, RT.ap), CR=CRlist, stopping.rule="OR")
z.super     <- ucip.test(RT=list(RT.pp.super, RT.pa, RT.ap), CR=CRlist, stopping.rule="OR")

#  AND Processing
RT.pp.limited <- pmax( rexp(100, .5*rate1), rexp(100, .5*rate2))
RT.pp.unlimited <- pmax( rexp(100, rate1), rexp(100, rate2))
RT.pp.super <- pmax( rexp(100, 2*rate1), rexp(100, 2*rate2))
z.limited   <- ucip.test(RT=list(RT.pp.limited, RT.pa, RT.ap), CR=CRlist, stopping.rule="AND")
z.unlimited <- ucip.test(RT=list(RT.pp.unlimited, RT.pa, RT.ap), CR=CRlist, stopping.rule="AND")
z.super     <- ucip.test(RT=list(RT.pp.super, RT.pa, RT.ap), CR=CRlist, stopping.rule="AND")



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
