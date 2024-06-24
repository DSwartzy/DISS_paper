#This file contains the code for the simulation for the DISS paper submitted to PLOS1 in July 2024.

# required packages ====
library(rdrobust)
library(dbscan)
library(ks)
library(RDHonest)
library(rdlocrand)
library(parallel)

# mean functions ====

#this is the design based on the Indiana data with different derivatives
design3.fn <- function(x){
  (.15-.15*x+2.5*x^2-1.5*x^3)*(x>=0) + (.05+1.5*x+3.2*x^2+2.7*x^3)*(x<0) + e
}

#this is a modified version of desgin 3 from AK
design1.fn <- function(x){
  (x+1)^2-2*pmax(x+.2,0)^2+2*pmax(x-.2,0)^2-2*pmax(x-.4,0)^2+2*pmax(x-.7,0)^2-.92+0.1*(x>=0) + e
}

#this is design 3 from from AI
design2.fn <- function(x) {
  (0.52+0.84*x-3.0*x^2+7.99*x^3-9.01*x^4+3.56*x^5)*(x>=0)+(0.42+0.84*x-3.0*x^2+7.99*x^3-9.01*x^4+3.56*x^5)*(x<0) + e
}

#this is the design based on the Indiana data with the same derivatives on either side
#had to shift it from what is in the power simulation
design4.fn <- function(x){
  -0.06*x^3+0.04*x^2+.39*x+0.53-0.1*(x>=0)+e
}

#this is design 4 from AI
design5.fn <- function(x) {
  (0.09 + 5.76*x - 42.56*x^2 + 120.90*x^3 - 139.71*x^4 + 55.59*x^5)*(x>=0) + (0.03 - 2.26*x - 13.14*x^2 - 30.89*x^3 - 31.98*x^4 - 12.1*x^5)*(x<0) + e
}

#this is a function that is flat at the cutoff but has curvature elsewhere
design6.fn <- function(x) {
  (-1*(x+0.2)^2)*(x<(-0.2))+0.1*(0<=x)*(x<=0.2)+((x-0.2)^2+0.1)*(x>0.2) + e
}

#this is a function that is flat everywhere
design7.fn <- function(x) {
  0.1*(x>=0) + e
}

#Functions for m-hat ====
#code in NPR_MROT.fit
NPR_MROT.fit <- function (d) 
{
  if (inherits(d, "RDData")) {
    max(NPR_MROT.fit(LPPData(data.frame(Y = d$Yp, X = d$Xp), 
                             0)), NPR_MROT.fit(LPPData(data.frame(Y = d$Ym, X = d$Xm), 
                                                       0)))
  }
  else if (inherits(d, "FRDData")) {
    c(M1 = max(NPR_MROT.fit(LPPData(data.frame(Y = d$Yp[, 
                                                        1], X = d$Xp), 0)), NPR_MROT.fit(LPPData(data.frame(Y = d$Ym[, 
                                                                                                                     1], X = d$Xm), 0))), M2 = max(NPR_MROT.fit(LPPData(data.frame(Y = d$Yp[, 
                                                                                                                                                                                            2], X = d$Xp), 0)), NPR_MROT.fit(LPPData(data.frame(Y = d$Ym[, 
                                                                                                                                                                                                                                                         2], X = d$Xm), 0))))
  }
  else if (inherits(d, "LPPData")) {
    r1 <- unname(stats::lm(d$Y ~ 0 + outer(d$X, 0:4, "^"))$coefficients)
    f2 <- function(x) abs(2 * r1[3] + 6 * x * r1[4] + 12 * 
                            x^2 * r1[5])
    f2e <- if (abs(r1[5]) <= 1e-10) 
      Inf
    else -r1[4]/(4 * r1[5])
    M <- max(f2(min(d$X)), f2(max(d$X)))
    if (min(d$X) < f2e & max(d$X) > f2e) 
      M <- max(f2(f2e), M)
    M
  }
}

#code in LPPData
LPPData <- function (d, point) 
{
  if (is.unsorted(d[[2]])) 
    d <- d[sort(unAsIs(d[[2]]), index.return = TRUE)$ix, 
    ]
  df <- list(Y = d[[1]], X = d[[2]] - point, orig.point = point, 
             var.names = names(d)[1:2])
  df$sigma2 <- d$sigma2
  df$w <- if (is.null(d$weights)) 
    rep(1L, length(df$X))
  else d$weights
  structure(df, class = "LPPData")
}

#code in unAsIs
unAsIs <- function (X) 
{
  if ("AsIs" %in% class(X)) {
    class(X) <- class(X)[-match("AsIs", class(X))]
  }
  X
}

#code in RDData
RDData <- function (d, cutoff) 
{
  if (is.unsorted(d[[2]])) 
    d <- d[sort(unAsIs(d[[2]]), index.return = TRUE)$ix, 
    ]
  X <- d[[2]] - cutoff
  df <- list(Ym = d[[1]][X < 0], Yp = d[[1]][X >= 0], Xm = X[X < 
                                                               0], Xp = X[X >= 0], orig.cutoff = cutoff, var.names = names(d)[1:2])
  df$sigma2m <- d$sigma2[X < 0]
  df$sigma2p <- d$sigma2[X >= 0]
  if (is.null(d$weights)) 
    d$weights <- rep(1L, length(X))
  df$wm <- d$weights[X < 0]
  df$wp <- d$weights[X >= 0]
  structure(df, class = "RDData")
}


#functions for ple ====

#Epanechnikov function
epanechnikov.fn <- function(x) {
  0.75*(1-x^2)*(abs(x)<=1)
}

#Note that I also added a valid input check.
#Note that you need to be very careful about using the "vector" option, as it may interact surprisingly with "apply" functions.
kernelweight.fn <- function(x.i, x.j, h, kernelname = "epanechnikov.fn", standardize.over="x.j") {
  #x.i is a vector of running variable values (or similar, like a single cutoff)
  #x.j is a vector of running variable values (or similar, like a single cutoff)
  #h is the scalar bandwidth
  #kernelname is the string name of the kernel function
  #standardize.over is an indicator of which vector the weights should be standardized over; options are:
  #"x.j" returns columns that sum to one
  #"x.i" returns rows that sum to one
  #"vector" returns a vector (not a matrix) that sums to one
  #"none" returns no normalization
  
  #this function takes two vectors and returns a **matrix** of normalized weights (even if x.i and x.j are scalars).
  #UNLESS the standardize.over option requests a vector, then it'll return a vector.
  
  #check to make sure some inputs are valid
  if(!(standardize.over) %in% c("x.j", "x.i", "vector", "none")) stop("you have asked for an invalid standardization.")
  if(standardize.over == "vector" & length(x.i)>1 & length(x.j)>1) stop("the vector option is only valid if at least one inputs are scalar.")
  if((standardize.over =="x.j" & length(x.i)==1) | (standardize.over =="x.i" & length(x.j)==1)) warning("I think you have chosen the incorrect standardization margin.") 
  
  #first, calculate the distance between each pair of elements in the two vectors.
  mydistance <- outer(x.i/h, x.j/h, "-") #this is a matrix where rows correspond to x.i an columns correspond to x.j
  #here we calculate unscaled weights using the generic name for the kernel function -- thus we can swap it out easily later.
  myunscaledweight <- do.call(kernelname, list(x=mydistance)) #this is a matrix of the same dimensions as mydistance
  
  #now we need to standardize
  if(standardize.over == "x.j")  myweights <- prop.table(myunscaledweight, margin = 2)
  if(standardize.over == "x.i")  myweights <- prop.table(myunscaledweight, margin = 1)
  if(standardize.over == "vector") myweights <- as.vector(myunscaledweight)/sum(as.vector(myunscaledweight))
  if(standardize.over == "none")  myweights <- myunscaledweight
  
  return(myweights)
}


#Weight function for local linear regression
LL.fn <- function(i, x, h, n){
  #requires function kernelweight.fn
  
  if(!("MASS" %in% .packages(all.available=TRUE))) stop("You must install the MASS package.")
  
  w.xi <- as.vector(kernelweight.fn(x[i], x, h, standardize.over="none"))
  
  #create the x design matrix
  x.xi <- cbind(rep(1,n), x-x[i])
  
  #create the e vector
  e.vec <- c(1,0)
  
  #get the results
  return((w.xi*x.xi)%*%MASS::ginv(t(x.xi)%*%(w.xi*x.xi))%*%e.vec)
}


#calculates the ple0 and ple1 point estimates and jackknifed variance estimates
ple.jack.fn <- function(y.vec, x.vec, c, h, mykernel="epanechnikov.fn"){
  
  #first getting the original estimate without anything deleted
  n <- length(x.vec)
  l.ll.matrix <- sapply(1:n, LL.fn, x=x.vec, h=h)
  l.kern.matrix <- kernelweight.fn(x.i=x.vec, x.j=x.vec, h=h, kernelname = mykernel, standardize.over = "x.j")
  f.vec <- as.numeric(x.vec>c)
  
  f.fitted.ple0 <- t(l.kern.matrix)%*%f.vec
  f.resid.ple0 <- f.vec-f.fitted.ple0
  
  y.fitted.ple0 <- t(l.kern.matrix)%*%y.vec
  y.resid.ple0 <- y.vec-y.fitted.ple0
  
  f.fitted.ple1 <- t(l.ll.matrix)%*%f.vec
  f.resid.ple1 <- f.vec-f.fitted.ple1
  
  y.fitted.ple1 <- t(l.ll.matrix)%*%y.vec
  y.resid.ple1 <- y.vec-y.fitted.ple1
  
  tau.ple0 <- sum(f.resid.ple0^2)^(-1)*sum(f.resid.ple0*y.resid.ple0)
  tau.ple1 <- sum(f.resid.ple1^2)^(-1)*sum(f.resid.ple1*y.resid.ple1)
  
  df <- data.frame(f.ple1=f.resid.ple1, y.ple1=y.resid.ple1, 
                   f.ple0=f.resid.ple0, y.ple0=y.resid.ple0, index=seq(1,n))
  
  #getting Wu's jackknife with formulas that don't require a loop
  
  hat.ple0 <- (f.resid.ple0)%*%MASS::ginv(t(f.resid.ple0)%*%(f.resid.ple0))%*%t(f.resid.ple0)
  r.ple0 <- (diag(n) - hat.ple0)%*%(y.resid.ple0)
  
  hat.ple1 <- (f.resid.ple1)%*%MASS::ginv(t(f.resid.ple1)%*%(f.resid.ple1))%*%t(f.resid.ple1)
  r.ple1 <- (diag(n) - hat.ple1)%*%(y.resid.ple1)
  
  se.ple0 <- sqrt(sum(df$f.ple0^2)^(-1)*sum(r.ple0^2/(1-diag(hat.ple0))*f.resid.ple0^2)*sum(df$f.ple0^2)^(-1))
  se.ple1 <- sqrt(sum(df$f.ple1^2)^(-1)*sum(r.ple1^2/(1-diag(hat.ple1))*f.resid.ple1^2)*sum(df$f.ple1^2)^(-1))
  
  #calculating the half-length of the interval to be more compatible with AK
  hl.ple0 <- 1.959963986*se.ple0
  hl.ple1 <- 1.959963986*se.ple1
  return(c(tau.ple0=tau.ple0, hl.ple0=hl.ple0, tau.ple1=tau.ple1, hl.ple1=hl.ple1))
  
}


#constants for the SKA bandwidth function
v.k <- 1792/121
a.p2 <- c(5/4, 15/7, 35/4)
b.p2 <- c(0,3/7,0)
a.p3 <- c(5/4, 525/44, 35/4, 1575/44)
b.p3 <- c(-1/21, 0, 2/3, 0)
a.p4 <- c(1575/832, 525/44, 19845/208, 1575/44, 121275/832)
b.p4 <- c(0,-5/33,0,10/11,0)

#SKA bandwidth function
skaband.fn <- function(y, x, c) {
  
  n <- length(x)
  
  #using options from the ks package
  
  h0.pilot1 <- ks::hpi(x, deriv.order=0)
  h1.pilot1 <- ks::hpi(x, deriv.order=1)
  h2.pilot1 <- ks::hpi(x, deriv.order=2)
  
  f0 <- ks::kdde(x, h=h0.pilot1, deriv.order=0, eval.points=c(c))$estimate
  f1 <- ks::kdde(x, h=h1.pilot1, deriv.order=1, eval.points=c(c))$estimate
  f2 <- ks::kdde(x, h=h2.pilot1, deriv.order=2, eval.points=c(c))$estimate
  
  #getting the variance estimates based on the nearest neighbor idea from CCT, with 3 neighbors
  
  x.u <- unique(x)
  y.u <- y[duplicated(x)==FALSE]
  x.right <- x.u[x.u>c]
  y.right <- y.u[x.u>c]
  right.index <- order(x.right)[1:3]
  x.right.nn <- x.right[right.index]
  y.right.nn <- y.right[right.index]
  x.left <- x[x<c]
  y.left <- y[x<c]
  left.index <- order(x.left, decreasing=TRUE)[1:3]
  x.left.nn <- x.left[left.index]
  y.left.nn <- y.left[left.index]
  sigma2.hat.right <- sum((y.right.nn-mean(y.right.nn))^2)/(length(y.right.nn)-1)
  sigma2.hat.left <- sum((y.left.nn-mean(y.left.nn))^2)/(length(y.left.nn)-1)
  
  #get the mean function derivative estimates
  
  #using a global quartic function to estimate the fourth derivative of the mean and sigma (here we assume equal variance)
  pilot.global4.lm <- lm(y ~ I(as.numeric((x-c)>0)) + I(x-c) + I((x-c)^2) + I((x-c)^3) + I((x-c)^4))
  m4.hat <- 24*pilot.global4.lm$coef[[6]]
  sigma.pilot4.hat <- summary(pilot.global4.lm)$sigma
  
  #using a global quintic function to estimate the fifth derivative of the mean and sigma (here we assume equal variance)
  pilot.global5.lm <- lm(y ~ I(as.numeric((x-c)>0)) + I(x-c) + I((x-c)^2) + I((x-c)^3) + I((x-c)^4) + I((x-c)^5))
  m5.hat <- 120*pilot.global5.lm$coef[[7]]
  sigma.pilot5.hat <- summary(pilot.global5.lm)$sigma
  
  #using a global cubic function to estimate the third derivative of the mean and sigma (here we assume equal variance)
  pilot.global3.lm <- lm(y ~ I(as.numeric((x-c)>0)) + I(x-c) + I((x-c)^2) + I((x-c)^3))
  m3.hat.1 <- 6*pilot.global3.lm$coef[[5]]
  sigma.pilot3.hat <- summary(pilot.global3.lm)$sigma
  
  #plugging in our estimates to the pilot bandwidth formula
  #when v=2 I can use the quartic (p=3)
  v <- 2
  h.pilot.v2 <- (((2*v+1)*(factorial(3+1))^2*(a.p3[v+1])*sigma.pilot4.hat^2)/(2*(4-v)*(b.p3[v+1]^2)*m4.hat^2*n*f0))^(1/9)
  
  #when v=3 I can use the quintic (p=4)
  v <- 3
  h.pilot.v3 <- (((2*v+1)*(factorial(4+1))^2*(a.p4[v+1])*sigma.pilot5.hat^2)/(2*(5-v)*(b.p4[v+1]^2)*m5.hat^2*n*f0))^(1/11)
  
  #when v=1 I can use the cubic (p=2)
  v <- 1
  h.pilot.v1 <- (((2*v+1)*(factorial(2+1))^2*(a.p2[v+1])*sigma.pilot3.hat^2)/(2*(3-v)*(b.p2[v+1]^2)*m3.hat.1^2*n*f0))^(1/7)
  
  #keeping only values within these bandwidths of the cutoff
  x.keep.v1 <- x[x>(c-h.pilot.v1) & x<(c+h.pilot.v1)]
  x.keep.v2 <- x[x>(c-h.pilot.v2) & x<(c+h.pilot.v2)]
  x.keep.v3 <- x[x>(c-h.pilot.v3) & x<(c+h.pilot.v3)]
  y.keep.v1 <- y[x>(c-h.pilot.v1) & x<(c+h.pilot.v1)]
  y.keep.v2 <- y[x>(c-h.pilot.v2) & x<(c+h.pilot.v2)]
  y.keep.v3 <- y[x>(c-h.pilot.v3) & x<(c+h.pilot.v3)]
  
  #running a v+1 degree polynomial jump for v=1,2,3
  local.v2.lm <- lm(y.keep.v2 ~ I(as.numeric((x.keep.v2-c)>0)) + I(x.keep.v2-c) + I((x.keep.v2-c)^2) + I((x.keep.v2-c)^3))
  local.v1.lm <- lm(y.keep.v1 ~ I(as.numeric((x.keep.v1-c)>0)) + I(x.keep.v1-c) + I((x.keep.v1-c)^2))
  local.v3.lm <- lm(y.keep.v3 ~ I(as.numeric((x.keep.v3-c)>0)) + I(x.keep.v3-c) + I((x.keep.v3-c)^2) + I((x.keep.v3-c)^3) + I((x.keep.v3-c)^4))
  
  m1.hat <- local.v1.lm$coefficient[[3]]
  m2.hat <- 2*local.v2.lm$coefficient[[4]]
  m3.hat <- 6*local.v3.lm$coefficient[[5]]
  
  #calculating b.b (will change with a different kernel)
  g2 <- m1.hat*f1+m2.hat*f0/2
  g2.prime <- m1.hat*f2+m2.hat*f1+(m2.hat*f1+m3.hat*f0)/2
  b.b <- 2*(1/10)*(f0*33/560)^(-1)*(f1/f0*g2*(1/10)-g2.prime*(1/20))
  
  #calculating the final bandwidth
  h.ska <- (v.k*(sigma2.hat.left+sigma2.hat.right)/n/b.b^2/24/f0)^(1/7)
  
  
  return(h.ska)
  
}

# Functions for estimation ====
tau.hat.one.y.fn <- function(y.vec, x.vec, c, dbound, cctkernel = "triangular"){
  
  n <- length(x.vec)
  
  #getting the value of mhat
  mhat <- NA
  try({
    df <- data.frame(y.vec, x.vec)
    df2 <- RDData(df, c=0)
    mhat <- NPR_MROT.fit(df2)
  }, silent=TRUE)
  
  #getting values for all bandwidths
  #first setting initial NAs for each bandwidth
  h.ska0 <- NA
  h.ik1 <- NA
  h.cct1 <- NA
  h.ak1 <- NA
  b.cct1 <- NA
  
  #now attempting the bandwidth values from rdrobust, rdhonest, and skaband.fn
  
  try({
    invisible(capture.output(fit.ik1 <- rdbwselect_2014(y.vec, x.vec, c=c, p=1, kernel=cctkernel, bwselect="IK")))
    h.ik1 <- unname(fit.ik1$bws[1,1]) 
  }, silent=TRUE)
  
  
  try({
    invisible(capture.output(fit.cct1 <- rdbwselect(y.vec, x.vec, c=c, p=1, kernel=cctkernel)))
    h.cct1 <- fit.cct1$bws[1,1] 
    b.cct1 <- fit.cct1$bws[1,3]
  }, silent=TRUE)
  
  try({
    suppressMessages(fit.akh1 <- RDHonest(y.vec~x.vec, cutoff=c, M=mhat, opt.criterion = "MSE"))
    h.ak1 <- unname(fit.akh1$coefficients$bandwidth)
  }, silent=TRUE)
  
  try({
    h.ska0 <- skaband.fn(y=y.vec, x=x.vec, c=c)
  }, silent=TRUE)
  
  #getting estimators for each of the bandwidth values above, with default NAs if the try block has an error
  
  rob.ik1 <- t(as.matrix(c(tau.us=NA_real_, tau.bc=NA_real_, se.us=NA_real_, se.rb=NA_real_)))
  rob.cct1 <- t(as.matrix(c(tau.us=NA_real_, tau.bc=NA_real_, se.us=NA_real_, se.rb=NA_real_)))
  rob.ak1 <- t(as.matrix(c(tau.us=NA_real_, tau.bc=NA_real_, se.us=NA_real_, se.rb=NA_real_)))
  rob.ska0 <- t(as.matrix(c(tau.us=NA_real_, tau.bc=NA_real_, se.us=NA_real_, se.rb=NA_real_)))
  
  hmod.ik1 <- list(coefficients=data.frame(estimates=NA_real_, std.error=NA_real_, conf.high=NA_real_))
  hmod.cct1 <-list(coefficients=data.frame(estimates=NA_real_, std.error=NA_real_, conf.high=NA_real_))
  hmod.ak1 <- list(coefficients=data.frame(estimates=NA_real_, std.error=NA_real_, conf.high=NA_real_))
  hmod.ska0 <- list(coefficients=data.frame(estimates=NA_real_, std.error=NA_real_, conf.high=NA_real_))
  
  plejack.ik1 <- c(tau.ple0=NA_real_, hl.ple0=NA_real_, tau.ple1=NA_real_, hl.ple1=NA_real_)
  plejack.cct1 <- c(tau.ple0=NA_real_, hl.ple0=NA_real_, tau.ple1=NA_real_, hl.ple1=NA_real_)
  plejack.ak1 <- c(tau.ple0=NA_real_, hl.ple0=NA_real_, tau.ple1=NA_real_, hl.ple1=NA_real_)
  plejack.ska0 <- c(tau.ple0=NA_real_, hl.ple0=NA_real_, tau.ple1=NA_real_, hl.ple1=NA_real_)
  
  #first the rdrobust estimators
  try({rob.ik1 <- rdrobust(y=y.vec, x=x.vec, c=c, p=1, h=h.ik1, kernel=cctkernel)$Estimate}, silent=TRUE)
  try({rob.cct1 <- rdrobust(y=y.vec, x=x.vec, c=c, p=1, h=h.cct1, b=b.cct1, kernel=cctkernel)$Estimate}, silent=TRUE)
  try({rob.ak1 <- rdrobust(y=y.vec, x=x.vec, c=c, p=1, h=h.ak1, kernel=cctkernel)$Estimate}, silent=TRUE)
  try({rob.ska0 <- rdrobust(y=y.vec, x=x.vec, c=c, p=1, h=h.ska0, kernel=cctkernel)$Estimate}, silent=TRUE)
  
  #then for the rdhonest estimators
  try({hmod.ik1 <- suppressMessages(RDHonest(y.vec~x.vec, cutoff=c, M=mhat, kern=cctkernel, h=h.ik1))}, silent=TRUE)
  try({hmod.cct1 <- suppressMessages(RDHonest(y.vec~x.vec, cutoff=c, M=mhat, kern=cctkernel, h=h.cct1))}, silent=TRUE)
  try({hmod.ak1 <- suppressMessages(RDHonest(y.vec~x.vec, cutoff=c, M=mhat, kern=cctkernel, h=h.ak1))}, silent=TRUE)
  try({hmod.ska0 <- suppressMessages(RDHonest(y.vec~x.vec, cutoff=c, M=mhat, kern=cctkernel, h=h.ska0))}, silent=TRUE)
  
  #then for the plejack estimators
  try({plejack.ik1 <- ple.jack.fn(y.vec=y.vec, x.vec=x.vec, c=c, h=h.ik1)}, silent=TRUE)
  try({plejack.cct1 <- ple.jack.fn(y.vec=y.vec, x.vec=x.vec, c=c, h=h.cct1)}, silent=TRUE)
  try({plejack.ak1 <- ple.jack.fn(y.vec=y.vec, x.vec=x.vec, c=c, h=h.ak1)}, silent=TRUE)
  try({plejack.ska0 <- ple.jack.fn(y.vec=y.vec, x.vec=x.vec, c=c, h=h.ska0)}, silent=TRUE)
  
  #extract values and give names to estimates
  
  #for the rdrobust data
  colnames(rob.ik1)[which(colnames(rob.ik1) == "tau.bc")] <- "tau.rb"
  colnames(rob.cct1)[which(colnames(rob.cct1) == "tau.bc")] <- "tau.rb"
  colnames(rob.ak1)[which(colnames(rob.ak1) == "tau.bc")] <- "tau.rb"
  colnames(rob.ska0)[which(colnames(rob.ska0) == "tau.bc")] <- "tau.rb"
  
  rob.ik1 <- as.data.frame(rob.ik1)
  rob.ik1$hl.us <- 1.959963986*rob.ik1$se.us
  rob.ik1$hl.rb <- 1.959963986*rob.ik1$se.rb
  rob.cct1 <- as.data.frame(rob.cct1)
  rob.cct1$hl.us <- 1.959963986*rob.cct1$se.us
  rob.cct1$hl.rb <- 1.959963986*rob.cct1$se.rb
  rob.ak1 <- as.data.frame(rob.ak1)
  rob.ak1$hl.us <- 1.959963986*rob.ak1$se.us
  rob.ak1$hl.rb <- 1.959963986*rob.ak1$se.rb
  rob.ska0 <- as.data.frame(rob.ska0)
  rob.ska0$hl.us <- 1.959963986*rob.ska0$se.us
  rob.ska0$hl.rb <- 1.959963986*rob.ska0$se.rb
  
  colnames(rob.ik1) <- paste(colnames(rob.ik1), "ik1", sep=".") 
  colnames(rob.cct1) <- paste(colnames(rob.cct1), "cct1", sep=".")
  colnames(rob.ak1) <- paste(colnames(rob.ak1), "ak1", sep=".")
  colnames(rob.ska0) <- paste(colnames(rob.ska0), "ska0", sep=".")
  
  #for the rdhonest data
  hon.ik1 <- t(as.matrix(c(hmod.ik1$coefficients$estimate, hmod.ik1$coefficients$std.error, 
                           hmod.ik1$coefficients$conf.high-hmod.ik1$coefficients$estimate)))
  colnames(hon.ik1) <- c("tau.flci.ik1", "se.flci.ik1", "hl.flci.ik1")
  
  hon.cct1 <- t(as.matrix(c(hmod.cct1$coefficients$estimate, hmod.cct1$coefficients$std.error, 
                            hmod.cct1$coefficients$conf.high-hmod.cct1$coefficients$estimate)))
  colnames(hon.cct1) <- c("tau.flci.cct1", "se.flci.cct1", "hl.flci.cct1")
  
  hon.ak1 <- t(as.matrix(c(hmod.ak1$coefficients$estimate, hmod.ak1$coefficients$std.error, 
                           hmod.ak1$coefficients$conf.high-hmod.ak1$coefficients$estimate)))
  colnames(hon.ak1) <- c("tau.flci.ak1", "se.flci.ak1", "hl.flci.ak1")
  
  hon.ska0 <- t(as.matrix(c(hmod.ska0$coefficients$estimate, hmod.ska0$coefficients$std.error, 
                            hmod.ska0$coefficients$conf.high-hmod.ska0$coefficients$estimate)))
  colnames(hon.ska0) <- c("tau.flci.ska0", "se.flci.ska0", "hl.flci.ska0")
  
  
  #for the plejack data 
  names(plejack.ik1) <- paste(names(plejack.ik1), "ik1", sep=".") 
  names(plejack.cct1) <- paste(names(plejack.cct1), "cct1", sep=".")
  names(plejack.ak1) <- paste(names(plejack.ak1), "ak1", sep=".")
  names(plejack.ska0) <- paste(names(plejack.ska0), "ska0", sep=".")
  
  plejack.ik1 <- t(as.matrix(plejack.ik1))
  plejack.cct1 <- t(as.matrix(plejack.cct1))
  plejack.ak1 <- t(as.matrix(plejack.ak1))
  plejack.ska0 <- t(as.matrix(plejack.ska0))
  
  #the local randomization stuff starts here
  
  #in case it doesn't work
  locrand5 <- t(as.matrix(c(tau.locrand.est5=NA, lb.locrand.est5=NA, ub.locrand.est5=NA)))
  locrand10 <- t(as.matrix(c(tau.locrand.est10=NA, lb.locrand.est10=NA, ub.locrand.est10=NA)))
  locrand20 <- t(as.matrix(c(tau.locrand.est20=NA, lb.locrand.est20=NA, ub.locrand.est20=NA)))
  
  #calculating the windows
  C=10
  x.below <- x.vec[x.vec<c]
  x.above <- x.vec[x.vec>=c]
  wsize.below <- ifelse(length(x.below)<=C, length(x.below), C) 
  length.below <- c-sort(x.below, decreasing=TRUE)[wsize.below]
  wsize.above <- ifelse(length(x.above)<=C, length(x.above), C) 
  length.above <- sort(x.above)[wsize.above]-c
  wlength10 <- max(length.above, length.below)
  C=5
  x.below <- x.vec[x.vec<c]
  x.above <- x.vec[x.vec>=c]
  wsize.below <- ifelse(length(x.below)<=C, length(x.below), C) 
  length.below <- c-sort(x.below, decreasing=TRUE)[wsize.below]
  wsize.above <- ifelse(length(x.above)<=C, length(x.above), C) 
  length.above <- sort(x.above)[wsize.above]-c
  wlength5 <- max(length.above, length.below)
  C=20
  x.below <- x.vec[x.vec<c]
  x.above <- x.vec[x.vec>=c]
  wsize.below <- ifelse(length(x.below)<=C, length(x.below), C) 
  length.below <- c-sort(x.below, decreasing=TRUE)[wsize.below]
  wsize.above <- ifelse(length(x.above)<=C, length(x.above), C) 
  length.above <- sort(x.above)[wsize.above]-c
  wlength20 <- max(length.above, length.below)
  
  est5 <- try({rdrandinf(Y=y.vec, R=x.vec, cutoff=c, wl=c-wlength5, wr=c+wlength5, p=0, ci=c(0.05), 
                         quietly=TRUE, seed=-1)}, silent=TRUE)
  est10 <- try({rdrandinf(Y=y.vec, R=x.vec, cutoff=c, wl=c-wlength10, wr=c+wlength10, p=0, ci=c(0.05), 
                          quietly=TRUE, seed=-1)}, silent=TRUE)
  est20 <- try({rdrandinf(Y=y.vec, R=x.vec, cutoff=c, wl=c-wlength20, wr=c+wlength20, p=0, ci=c(0.05), 
                          quietly=TRUE, seed=-1)}, silent=TRUE)
  locrand5 <- try({t(as.matrix(c(tau.locrand.est5=est5$obs.stat, lb.locrand.est5=est5$ci[1], 
                                 ub.locrand.est5=est5$ci[2])))}, silent=TRUE)
  locrand10 <- try({t(as.matrix(c(tau.locrand.est10=est10$obs.stat, lb.locrand.est10=est10$ci[1], 
                                  ub.locrand.est10=est10$ci[2])))}, silent=TRUE)
  locrand20 <- try({t(as.matrix(c(tau.locrand.est20=est20$obs.stat, lb.locrand.est20=est20$ci[1], 
                                  ub.locrand.est20=est20$ci[2])))}, silent=TRUE)
  
  #combining bandwdiths
  bw <- t(as.matrix(c(bw.h.ik1=h.ik1, bw.h.cct1=h.cct1, bw.h.ak1=h.ak1, bw.h.ska0=h.ska0,
                      bw.b.cct1=b.cct1, bw.h.est5=wlength5,
                      bw.h.est10=wlength10, bw.h.est20=wlength20,
                      mhat=mhat)))
  
  return(estimates=cbind(rob.ik1, rob.cct1, rob.ak1, rob.ska0, hon.ik1, hon.cct1, hon.ak1, hon.ska0,
                         plejack.ik1, plejack.cct1, plejack.ak1, plejack.ska0, 
                         locrand5, locrand10, locrand20, bw))
  
}



# setup for RV1 (AK) ====
S=50000 # number of iterations of the simulation
c=0 #cutoff for purposes of eda below

# Data generation parameters
alpha <- 1
beta <- 1

#setting up the random number generator
RNGkind("L'Ecuyer-CMRG")

set.seed(8145) #initial seed
current.seed <- .Random.seed

n <- 40 #m=10
for (i in 1:S){

  current.seed <- nextRNGStream(current.seed)
  .GlobalEnv$.Random.seed <- current.seed

  #Data generation

  #Generate values of the running variable, with half below the cutoff and half above

  z <- rbeta(n, shape1=alpha, shape2=beta)
  x <- 2*z-1

  #Generate error values
  e <- rnorm(n,0,.1295)

  #examining the distribution of x
  count.below <- sum(x<c)
  count.above <- sum(x>=c)
  h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
  rot.below <- sum(x<c & x>(c-h.rot))
  rot.above <- sum(x>c & x<(c+h.rot))

  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                       rot.below=rot.below, rot.above=rot.above)))

  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn",
                   "design4.fn", "design5.fn", "design6.fn", "design7.fn"),
                 function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)

  #Get estimates
  tauhat.design1 <- tau.hat.one.y.fn(y.vec=ymat[,1], x.vec=x, c=0)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")

  tauhat.design2 <- tau.hat.one.y.fn(y.vec=ymat[,2], x.vec=x, c=0)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")

  tauhat.design3 <- tau.hat.one.y.fn(y.vec=ymat[,3], x.vec=x, c=0)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")

  tauhat.design4 <- tau.hat.one.y.fn(y.vec=ymat[,4], x.vec=x, c=0)
  colnames(tauhat.design4) <- paste(colnames(tauhat.design4), "design4", sep=".")

  tauhat.design5 <- tau.hat.one.y.fn(y.vec=ymat[,5], x.vec=x, c=0)
  colnames(tauhat.design5) <- paste(colnames(tauhat.design5), "design5", sep=".")

  tauhat.design6 <- tau.hat.one.y.fn(y.vec=ymat[,6], x.vec=x, c=0)
  colnames(tauhat.design6) <- paste(colnames(tauhat.design6), "design6", sep=".")

  tauhat.design7 <- tau.hat.one.y.fn(y.vec=ymat[,7], x.vec=x, c=0)
  colnames(tauhat.design7) <- paste(colnames(tauhat.design7), "design7", sep=".")

  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3,
                  tauhat.design4, tauhat.design5, tauhat.design6, tauhat.design7, eda)


  #Output data
  if (i==1) write.table(output, file="paper2comp031424rv1m10estimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="paper2comp031424rv1m10estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}

n <- 101 #m=21
set.seed(50)
current.seed <- .Random.seed
for (i in 1:S){
current.seed <- nextRNGStream(current.seed)
.GlobalEnv$.Random.seed <- current.seed

#Data generation

#Generate values of the running variable, with half below the cutoff and half above

z <- rbeta(n, shape1=alpha, shape2=beta)
x <- 2*z-1

#Generate error values
e <- rnorm(n,0,.1295)

#examining the distribution of x
count.below <- sum(x<c)
count.above <- sum(x>=c)
h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
rot.below <- sum(x<c & x>(c-h.rot))
rot.above <- sum(x>c & x<(c+h.rot))

eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                     rot.below=rot.below, rot.above=rot.above)))

#Generate y matrix, one column for each function
ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn",
                 "design4.fn", "design5.fn", "design6.fn", "design7.fn"),
               function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)

#Get estimates
tauhat.design1 <- tau.hat.one.y.fn(y.vec=ymat[,1], x.vec=x, c=0)
colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")

tauhat.design2 <- tau.hat.one.y.fn(y.vec=ymat[,2], x.vec=x, c=0)
colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")

tauhat.design3 <- tau.hat.one.y.fn(y.vec=ymat[,3], x.vec=x, c=0)
colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")

tauhat.design4 <- tau.hat.one.y.fn(y.vec=ymat[,4], x.vec=x, c=0)
colnames(tauhat.design4) <- paste(colnames(tauhat.design4), "design4", sep=".")

tauhat.design5 <- tau.hat.one.y.fn(y.vec=ymat[,5], x.vec=x, c=0)
colnames(tauhat.design5) <- paste(colnames(tauhat.design5), "design5", sep=".")

tauhat.design6 <- tau.hat.one.y.fn(y.vec=ymat[,6], x.vec=x, c=0)
colnames(tauhat.design6) <- paste(colnames(tauhat.design6), "design6", sep=".")

tauhat.design7 <- tau.hat.one.y.fn(y.vec=ymat[,7], x.vec=x, c=0)
colnames(tauhat.design7) <- paste(colnames(tauhat.design7), "design7", sep=".")

iteration <- i
output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3,
                tauhat.design4, tauhat.design5, tauhat.design6, tauhat.design7, eda)


#Output data
if (i==1) write.table(output, file="paper2comp031424rv1m21estimates.csv", sep=",", row.names=FALSE)
else write.table(output, file="paper2comp031424rv1m21estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}


n <- 140 #m=27
set.seed(9869)
current.seed <- .Random.seed
for (i in 1:S){

  current.seed <- nextRNGStream(current.seed)
  .GlobalEnv$.Random.seed <- current.seed

  #Data generation

  #Generate values of the running variable, with half below the cutoff and half above

  z <- rbeta(n, shape1=alpha, shape2=beta)
  x <- 2*z-1

  #Generate error values
  e <- rnorm(n,0,.1295)

  #examining the distribution of x
  count.below <- sum(x<c)
  count.above <- sum(x>=c)
  h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
  rot.below <- sum(x<c & x>(c-h.rot))
  rot.above <- sum(x>c & x<(c+h.rot))

  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                       rot.below=rot.below, rot.above=rot.above)))

  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn",
                   "design4.fn", "design5.fn", "design6.fn", "design7.fn"),
                 function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)

  #Get estimates
  tauhat.design1 <- tau.hat.one.y.fn(y.vec=ymat[,1], x.vec=x, c=0)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")

  tauhat.design2 <- tau.hat.one.y.fn(y.vec=ymat[,2], x.vec=x, c=0)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")

  tauhat.design3 <- tau.hat.one.y.fn(y.vec=ymat[,3], x.vec=x, c=0)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")

  tauhat.design4 <- tau.hat.one.y.fn(y.vec=ymat[,4], x.vec=x, c=0)
  colnames(tauhat.design4) <- paste(colnames(tauhat.design4), "design4", sep=".")

  tauhat.design5 <- tau.hat.one.y.fn(y.vec=ymat[,5], x.vec=x, c=0)
  colnames(tauhat.design5) <- paste(colnames(tauhat.design5), "design5", sep=".")

  tauhat.design6 <- tau.hat.one.y.fn(y.vec=ymat[,6], x.vec=x, c=0)
  colnames(tauhat.design6) <- paste(colnames(tauhat.design6), "design6", sep=".")

  tauhat.design7 <- tau.hat.one.y.fn(y.vec=ymat[,7], x.vec=x, c=0)
  colnames(tauhat.design7) <- paste(colnames(tauhat.design7), "design7", sep=".")

  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3,
                  tauhat.design4, tauhat.design5, tauhat.design6, tauhat.design7, eda)


  #Output data
  if (i==1) write.table(output, file="paper2comp031424rv1m27estimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="paper2comp031424rv1m27estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}

n <- 256 #m=44
set.seed(2694)
current.seed <- .Random.seed
for (i in 1:S){

  current.seed <- nextRNGStream(current.seed)
  .GlobalEnv$.Random.seed <- current.seed

  #Data generation

  #Generate values of the running variable, with half below the cutoff and half above

  z <- rbeta(n, shape1=alpha, shape2=beta)
  x <- 2*z-1

  #Generate error values
  e <- rnorm(n,0,.1295)

  #examining the distribution of x
  count.below <- sum(x<c)
  count.above <- sum(x>=c)
  h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
  rot.below <- sum(x<c & x>(c-h.rot))
  rot.above <- sum(x>c & x<(c+h.rot))

  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                       rot.below=rot.below, rot.above=rot.above)))

  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn",
                   "design4.fn", "design5.fn", "design6.fn", "design7.fn"),
                 function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)

  #Get estimates
  tauhat.design1 <- tau.hat.one.y.fn(y.vec=ymat[,1], x.vec=x, c=0)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")

  tauhat.design2 <- tau.hat.one.y.fn(y.vec=ymat[,2], x.vec=x, c=0)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")

  tauhat.design3 <- tau.hat.one.y.fn(y.vec=ymat[,3], x.vec=x, c=0)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")

  tauhat.design4 <- tau.hat.one.y.fn(y.vec=ymat[,4], x.vec=x, c=0)
  colnames(tauhat.design4) <- paste(colnames(tauhat.design4), "design4", sep=".")

  tauhat.design5 <- tau.hat.one.y.fn(y.vec=ymat[,5], x.vec=x, c=0)
  colnames(tauhat.design5) <- paste(colnames(tauhat.design5), "design5", sep=".")

  tauhat.design6 <- tau.hat.one.y.fn(y.vec=ymat[,6], x.vec=x, c=0)
  colnames(tauhat.design6) <- paste(colnames(tauhat.design6), "design6", sep=".")

  tauhat.design7 <- tau.hat.one.y.fn(y.vec=ymat[,7], x.vec=x, c=0)
  colnames(tauhat.design7) <- paste(colnames(tauhat.design7), "design7", sep=".")

  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3,
                  tauhat.design4, tauhat.design5, tauhat.design6, tauhat.design7, eda)


  #Output data
  if (i==1) write.table(output, file="paper2comp031424rv1m44estimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="paper2comp031424rv1m44estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}

n <- 354 #m=57
set.seed(6927)
current.seed <- .Random.seed
for (i in 1:S){
  

  current.seed <- nextRNGStream(current.seed)
  .GlobalEnv$.Random.seed <- current.seed
  
  #Data generation
  
  #Generate values of the running variable, with half below the cutoff and half above
  
  z <- rbeta(n, shape1=alpha, shape2=beta)
  x <- 2*z-1
  
  #Generate error values
  e <- rnorm(n,0,.1295)
  
  #examining the distribution of x
  count.below <- sum(x<c)
  count.above <- sum(x>=c)
  h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
  rot.below <- sum(x<c & x>(c-h.rot))
  rot.above <- sum(x>c & x<(c+h.rot))
  
  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                       rot.below=rot.below, rot.above=rot.above)))
  
  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn",
                   "design4.fn", "design5.fn", "design6.fn", "design7.fn"),
                 function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)
  
  #Get estimates
  tauhat.design1 <- tau.hat.one.y.fn(y.vec=ymat[,1], x.vec=x, c=0)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")
  
  tauhat.design2 <- tau.hat.one.y.fn(y.vec=ymat[,2], x.vec=x, c=0)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")
  
  tauhat.design3 <- tau.hat.one.y.fn(y.vec=ymat[,3], x.vec=x, c=0)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")
  
  tauhat.design4 <- tau.hat.one.y.fn(y.vec=ymat[,4], x.vec=x, c=0)
  colnames(tauhat.design4) <- paste(colnames(tauhat.design4), "design4", sep=".")
  
  tauhat.design5 <- tau.hat.one.y.fn(y.vec=ymat[,5], x.vec=x, c=0)
  colnames(tauhat.design5) <- paste(colnames(tauhat.design5), "design5", sep=".")
  
  tauhat.design6 <- tau.hat.one.y.fn(y.vec=ymat[,6], x.vec=x, c=0)
  colnames(tauhat.design6) <- paste(colnames(tauhat.design6), "design6", sep=".")
  
  tauhat.design7 <- tau.hat.one.y.fn(y.vec=ymat[,7], x.vec=x, c=0)
  colnames(tauhat.design7) <- paste(colnames(tauhat.design7), "design7", sep=".")
  
  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3,
                  tauhat.design4, tauhat.design5, tauhat.design6, tauhat.design7, eda)
  
  
  #Output data
  if (i==1) write.table(output, file="paper2comp031424rv1m57estimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="paper2comp031424rv1m57estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}



#for the sensitivity analysis in the appendix
tau.hat.sens <- function(y.vec, x.vec, c, dbound, cctkernel = "triangular"){

  n <- length(x.vec)

  #getting the value of mhat
  mhat <- NA
  try({
    df <- data.frame(y.vec, x.vec)
    df2 <- RDData(df, c=0)
    mhat <- NPR_MROT.fit(df2)
  }, silent=TRUE)

  #getting values for all bandwidths
  #first setting initial NAs for each bandwidth
  h.ska0 <- NA
  h.ik1 <- NA
  h.cct1 <- NA
  h.ak1 <- NA
  b.cct1 <- NA
  h.akm1 <- NA

  #now attempting the bandwidth values from rdrobust, rdhonest, and skaband.fn

  try({
    invisible(capture.output(fit.ik1 <- rdbwselect_2014(y.vec, x.vec, c=c, p=1, kernel=cctkernel, bwselect="IK")))
    h.ik1 <- unname(fit.ik1$bws[1,1])
  }, silent=TRUE)


  try({
    invisible(capture.output(fit.cct1 <- rdbwselect(y.vec, x.vec, c=c, p=1, kernel=cctkernel)))
    h.cct1 <- fit.cct1$bws[1,1]
    b.cct1 <- fit.cct1$bws[1,3]
  }, silent=TRUE)

  try({
    suppressMessages(fit.akm1 <- RDHonest(y.vec~x.vec, cutoff=c, M=dbound, opt.criterion = "MSE"))
    h.akm1 <- unname(fit.akm1$coefficients$bandwidth)
  }, silent=TRUE)

  try({
    suppressMessages(fit.akh1 <- RDHonest(y.vec~x.vec, cutoff=c, M=mhat, opt.criterion = "MSE"))
    h.ak1 <- unname(fit.akh1$coefficients$bandwidth)
  }, silent=TRUE)

  try({
    h.ska0 <- skaband.fn(y=y.vec, x=x.vec, c=c)
  }, silent=TRUE)

  #getting estimators for each of the bandwidth values above, with default NAs if the try block has an error

  hmod.ik1 <- list(coefficients=data.frame(estimates=NA_real_, std.error=NA_real_, conf.high=NA_real_))
  hmod.cct1 <-list(coefficients=data.frame(estimates=NA_real_, std.error=NA_real_, conf.high=NA_real_))
  hmod.ak1 <- list(coefficients=data.frame(estimates=NA_real_, std.error=NA_real_, conf.high=NA_real_))
  hmod.akm1 <- list(coefficients=data.frame(estimates=NA_real_, std.error=NA_real_, conf.high=NA_real_))
  hmod.ska0 <- list(coefficients=data.frame(estimates=NA_real_, std.error=NA_real_, conf.high=NA_real_))

  hmod2.akm1 <- list(coefficients=data.frame(estimates=NA_real_, std.error=NA_real_, conf.high=NA_real_))


  #then for the rdhonest estimators
  try({hmod.ik1 <- suppressMessages(RDHonest(y.vec~x.vec, cutoff=c, M=dbound, kern=cctkernel, h=h.ik1))}, silent=TRUE)
  try({hmod.cct1 <- suppressMessages(RDHonest(y.vec~x.vec, cutoff=c, M=dbound, kern=cctkernel, h=h.cct1))}, silent=TRUE)
  try({hmod.ak1 <- suppressMessages(RDHonest(y.vec~x.vec, cutoff=c, M=dbound, kern=cctkernel, h=h.ak1))}, silent=TRUE)
  try({hmod.akm1 <- suppressMessages(RDHonest(y.vec~x.vec, cutoff=c, M=dbound, kern=cctkernel, h=h.akm1))}, silent=TRUE)
  try({hmod.ska0 <- suppressMessages(RDHonest(y.vec~x.vec, cutoff=c, M=dbound, kern=cctkernel, h=h.ska0))}, silent=TRUE)
  try({hmod2.akm1 <- suppressMessages(RDHonest(y.vec~x.vec, cutoff=c, M=mhat, kern=cctkernel, h=h.akm1))}, silent=TRUE)

  #extract values and give names to estimates

  #for the rdhonest data
  hon.ik1 <- t(as.matrix(c(hmod.ik1$coefficients$estimate, hmod.ik1$coefficients$std.error,
                           hmod.ik1$coefficients$conf.high-hmod.ik1$coefficients$estimate)))
  colnames(hon.ik1) <- c("tau.flcim.ik1", "se.flcim.ik1", "hl.flcim.ik1")

  hon.cct1 <- t(as.matrix(c(hmod.cct1$coefficients$estimate, hmod.cct1$coefficients$std.error,
                            hmod.cct1$coefficients$conf.high-hmod.cct1$coefficients$estimate)))
  colnames(hon.cct1) <- c("tau.flcim.cct1", "se.flcim.cct1", "hl.flcim.cct1")

  hon.ak1 <- t(as.matrix(c(hmod.ak1$coefficients$estimate, hmod.ak1$coefficients$std.error,
                           hmod.ak1$coefficients$conf.high-hmod.ak1$coefficients$estimate)))
  colnames(hon.ak1) <- c("tau.flcim.ak1", "se.flcim.ak1", "hl.flcim.ak1")

  hon.akm1 <- t(as.matrix(c(hmod.akm1$coefficients$estimate, hmod.akm1$coefficients$std.error,
                            hmod.akm1$coefficients$conf.high-hmod.akm1$coefficients$estimate)))
  colnames(hon.akm1) <- c("tau.flcim.akm1", "se.flcim.akm1", "hl.flcim.akm1")

  hon.ska0 <- t(as.matrix(c(hmod.ska0$coefficients$estimate, hmod.ska0$coefficients$std.error,
                            hmod.ska0$coefficients$conf.high-hmod.ska0$coefficients$estimate)))
  colnames(hon.ska0) <- c("tau.flcim.ska0", "se.flcim.ska0", "hl.flcim.ska0")

  hon2.akm1 <- t(as.matrix(c(hmod2.akm1$coefficients$estimate, hmod2.akm1$coefficients$std.error,
                             hmod2.akm1$coefficients$conf.high-hmod2.akm1$coefficients$estimate)))
  colnames(hon2.akm1) <- c("tau.flci.akm1", "se.flci.akm1", "hl.flci.akm1")

  #combining bandwdiths
  bw <- t(as.matrix(c(bw.h.ik1=h.ik1, bw.h.cct1=h.cct1, bw.h.ak1=h.ak1, bw.h.ska0=h.ska0,
                      bw.b.cct1=b.cct1, bw.h.akm1=h.akm1,
                      mhat=mhat)))

  return(estimates=cbind(hon.ik1, hon.cct1, hon.ak1, hon.akm1, hon.ska0, hon2.akm1, bw))

}

# setup for RV1 (AK) ====
S=50000 # number of iterations of the simulation
c=0 #cutoff for purposes of eda below

# Data generation parameters
alpha <- 1
beta <- 1

#setting up the random number generator
RNGkind("L'Ecuyer-CMRG")

set.seed(8145) #initial seed
current.seed <- .Random.seed

n <- 40 #m=10
for (i in 1:S){

  current.seed <- nextRNGStream(current.seed)
  .GlobalEnv$.Random.seed <- current.seed

  #Data generation

  #Generate values of the running variable, with half below the cutoff and half above

  z <- rbeta(n, shape1=alpha, shape2=beta)
  x <- 2*z-1

  #Generate error values
  e <- rnorm(n,0,.1295)

  #examining the distribution of x
  count.below <- sum(x<c)
  count.above <- sum(x>=c)
  h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
  rot.below <- sum(x<c & x>(c-h.rot))
  rot.above <- sum(x>c & x<(c+h.rot))

  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                       rot.below=rot.below, rot.above=rot.above)))

  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn"),
                 function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)

  #Get estimates
  tauhat.design1 <- tau.hat.sens(y.vec=ymat[,1], x.vec=x, c=0, dbound=2)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")

  tauhat.design2 <- tau.hat.sens(y.vec=ymat[,2], x.vec=x, c=0, dbound=233.26)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")

  tauhat.design3 <- tau.hat.sens(y.vec=ymat[,3], x.vec=x, c=0, dbound=16.2)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")

  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3, eda)


  #Output data
  if (i==1) write.table(output, file="paper2comp031424rv1m10sestimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="paper2comp031424rv1m10sestimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}

n <- 101 #m=21
set.seed(50)
current.seed <- .Random.seed
for (i in 1:S){
  current.seed <- nextRNGStream(current.seed)
  .GlobalEnv$.Random.seed <- current.seed

  #Data generation

  #Generate values of the running variable, with half below the cutoff and half above

  z <- rbeta(n, shape1=alpha, shape2=beta)
  x <- 2*z-1

  #Generate error values
  e <- rnorm(n,0,.1295)

  #examining the distribution of x
  count.below <- sum(x<c)
  count.above <- sum(x>=c)
  h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
  rot.below <- sum(x<c & x>(c-h.rot))
  rot.above <- sum(x>c & x<(c+h.rot))

  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                       rot.below=rot.below, rot.above=rot.above)))

  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn"),
                 function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)

  #Get estimates
  tauhat.design1 <- tau.hat.sens(y.vec=ymat[,1], x.vec=x, c=0, dbound=2)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")

  tauhat.design2 <- tau.hat.sens(y.vec=ymat[,2], x.vec=x, c=0, dbound=233.26)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")

  tauhat.design3 <- tau.hat.sens(y.vec=ymat[,3], x.vec=x, c=0, dbound=16.2)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")

  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3, eda)
  #Output data
  if (i==1) write.table(output, file="paper2comp031424rv1m21sestimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="paper2comp031424rv1m21sestimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}


n <- 140 #m=27
set.seed(9869)
current.seed <- .Random.seed
for (i in 1:S){

  current.seed <- nextRNGStream(current.seed)
  .GlobalEnv$.Random.seed <- current.seed

  #Data generation

  #Generate values of the running variable, with half below the cutoff and half above

  z <- rbeta(n, shape1=alpha, shape2=beta)
  x <- 2*z-1

  #Generate error values
  e <- rnorm(n,0,.1295)

  #examining the distribution of x
  count.below <- sum(x<c)
  count.above <- sum(x>=c)
  h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
  rot.below <- sum(x<c & x>(c-h.rot))
  rot.above <- sum(x>c & x<(c+h.rot))

  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                       rot.below=rot.below, rot.above=rot.above)))

  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn"),
                 function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)

  #Get estimates
  tauhat.design1 <- tau.hat.sens(y.vec=ymat[,1], x.vec=x, c=0, dbound=2)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")

  tauhat.design2 <- tau.hat.sens(y.vec=ymat[,2], x.vec=x, c=0, dbound=233.26)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")

  tauhat.design3 <- tau.hat.sens(y.vec=ymat[,3], x.vec=x, c=0, dbound=16.2)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")

  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3, eda)

  #Output data
  if (i==1) write.table(output, file="paper2comp031424rv1m27sestimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="paper2comp031424rv1m27sestimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}

n <- 256 #m=44
set.seed(2694)
current.seed <- .Random.seed
for (i in 1:S){

  current.seed <- nextRNGStream(current.seed)
  .GlobalEnv$.Random.seed <- current.seed

  #Data generation

  #Generate values of the running variable, with half below the cutoff and half above

  z <- rbeta(n, shape1=alpha, shape2=beta)
  x <- 2*z-1

  #Generate error values
  e <- rnorm(n,0,.1295)

  #examining the distribution of x
  count.below <- sum(x<c)
  count.above <- sum(x>=c)
  h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
  rot.below <- sum(x<c & x>(c-h.rot))
  rot.above <- sum(x>c & x<(c+h.rot))

  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                       rot.below=rot.below, rot.above=rot.above)))

  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn"),
                 function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)

  #Get estimates
  tauhat.design1 <- tau.hat.sens(y.vec=ymat[,1], x.vec=x, c=0, dbound=2)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")

  tauhat.design2 <- tau.hat.sens(y.vec=ymat[,2], x.vec=x, c=0, dbound=233.26)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")

  tauhat.design3 <- tau.hat.sens(y.vec=ymat[,3], x.vec=x, c=0, dbound=16.2)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")

  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3, eda)

  #Output data
  if (i==1) write.table(output, file="paper2comp031424rv1m44sestimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="paper2comp031424rv1sm44estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}

n <- 354 #m=57
set.seed(6927)
current.seed <- .Random.seed

for (i in 1:S){

  current.seed <- nextRNGStream(current.seed)
  .GlobalEnv$.Random.seed <- current.seed

  #Data generation

  #Generate values of the running variable, with half below the cutoff and half above

  z <- rbeta(n, shape1=alpha, shape2=beta)
  x <- 2*z-1

  #Generate error values
  e <- rnorm(n,0,.1295)

  #examining the distribution of x
  count.below <- sum(x<c)
  count.above <- sum(x>=c)
  h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
  rot.below <- sum(x<c & x>(c-h.rot))
  rot.above <- sum(x>c & x<(c+h.rot))

  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                       rot.below=rot.below, rot.above=rot.above)))

  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn"),
                 function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)

  #Get estimates
  tauhat.design1 <- tau.hat.sens(y.vec=ymat[,1], x.vec=x, c=0, dbound=2)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")

  tauhat.design2 <- tau.hat.sens(y.vec=ymat[,2], x.vec=x, c=0, dbound=233.26)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")

  tauhat.design3 <- tau.hat.sens(y.vec=ymat[,3], x.vec=x, c=0, dbound=16.2)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")

  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3, eda)

  #Output data
  if (i==1) write.table(output, file="paper2comp031424rv1m57sestimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="paper2comp031424rv1m57sestimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}

# setup for RV2 ====
S=50000 # number of iterations of the simulation
c=0 #cutoff for purposes of eda below

# Data generation parameters
alpha <- 2
beta <- 4

#setting up the random number generator
RNGkind("L'Ecuyer-CMRG")

set.seed(7980) #initial seed
current.seed <- .Random.seed

n <- 56 #m=10
for (i in 1:S){

  current.seed <- nextRNGStream(current.seed)
  .GlobalEnv$.Random.seed <- current.seed

  #Data generation

  #Generate values of the running variable, with half below the cutoff and half above

  z <- rbeta(n, shape1=alpha, shape2=beta)
  x <- 2*z-1

  #Generate error values
  e <- rnorm(n,0,.1295)

  #examining the distribution of x
  count.below <- sum(x<c)
  count.above <- sum(x>=c)
  h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
  rot.below <- sum(x<c & x>(c-h.rot))
  rot.above <- sum(x>c & x<(c+h.rot))

  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                       rot.below=rot.below, rot.above=rot.above)))

  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn",
                   "design4.fn", "design5.fn", "design6.fn", "design7.fn"),
                 function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)

  #Get estimates
  tauhat.design1 <- tau.hat.one.y.fn(y.vec=ymat[,1], x.vec=x, c=0)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")

  tauhat.design2 <- tau.hat.one.y.fn(y.vec=ymat[,2], x.vec=x, c=0)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")

  tauhat.design3 <- tau.hat.one.y.fn(y.vec=ymat[,3], x.vec=x, c=0)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")

  tauhat.design4 <- tau.hat.one.y.fn(y.vec=ymat[,4], x.vec=x, c=0)
  colnames(tauhat.design4) <- paste(colnames(tauhat.design4), "design4", sep=".")

  tauhat.design5 <- tau.hat.one.y.fn(y.vec=ymat[,5], x.vec=x, c=0)
  colnames(tauhat.design5) <- paste(colnames(tauhat.design5), "design5", sep=".")

  tauhat.design6 <- tau.hat.one.y.fn(y.vec=ymat[,6], x.vec=x, c=0)
  colnames(tauhat.design6) <- paste(colnames(tauhat.design6), "design6", sep=".")

  tauhat.design7 <- tau.hat.one.y.fn(y.vec=ymat[,7], x.vec=x, c=0)
  colnames(tauhat.design7) <- paste(colnames(tauhat.design7), "design7", sep=".")

  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3,
                  tauhat.design4, tauhat.design5, tauhat.design6, tauhat.design7, eda)


  #Output data
  if (i==1) write.table(output, file="paper2comp031424rv2m10estimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="paper2comp031424rv2m10estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}

n <- 140 #m=21
set.seed(4620)
current.seed <- .Random.seed
for (i in 1:S){
     current.seed <- nextRNGStream(current.seed)
      .GlobalEnv$.Random.seed <- current.seed

#Data generation

#Generate values of the running variable, with half below the cutoff and half above

z <- rbeta(n, shape1=alpha, shape2=beta)
x <- 2*z-1

#Generate error values
e <- rnorm(n,0,.1295)

#examining the distribution of x
count.below <- sum(x<c)
count.above <- sum(x>=c)
h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
rot.below <- sum(x<c & x>(c-h.rot))
rot.above <- sum(x>c & x<(c+h.rot))

eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                     rot.below=rot.below, rot.above=rot.above)))

#Generate y matrix, one column for each function
ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn",
                 "design4.fn", "design5.fn", "design6.fn", "design7.fn"),
               function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)

#Get estimates
tauhat.design1 <- tau.hat.one.y.fn(y.vec=ymat[,1], x.vec=x, c=0)
colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")

tauhat.design2 <- tau.hat.one.y.fn(y.vec=ymat[,2], x.vec=x, c=0)
colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")

tauhat.design3 <- tau.hat.one.y.fn(y.vec=ymat[,3], x.vec=x, c=0)
colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")

tauhat.design4 <- tau.hat.one.y.fn(y.vec=ymat[,4], x.vec=x, c=0)
colnames(tauhat.design4) <- paste(colnames(tauhat.design4), "design4", sep=".")

tauhat.design5 <- tau.hat.one.y.fn(y.vec=ymat[,5], x.vec=x, c=0)
colnames(tauhat.design5) <- paste(colnames(tauhat.design5), "design5", sep=".")

tauhat.design6 <- tau.hat.one.y.fn(y.vec=ymat[,6], x.vec=x, c=0)
colnames(tauhat.design6) <- paste(colnames(tauhat.design6), "design6", sep=".")

tauhat.design7 <- tau.hat.one.y.fn(y.vec=ymat[,7], x.vec=x, c=0)
colnames(tauhat.design7) <- paste(colnames(tauhat.design7), "design7", sep=".")

iteration <- i
output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3,
                tauhat.design4, tauhat.design5, tauhat.design6, tauhat.design7, eda)


#Output data
if (i==1) write.table(output, file="paper2comp031424rv2m21estimates.csv", sep=",", row.names=FALSE)
else write.table(output, file="paper2comp031424rv2m21estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}


n <- 194 #m=27
set.seed(8766)
current.seed <- .Random.seed
for (i in 1:S){
  
  current.seed <- nextRNGStream(current.seed)
  .GlobalEnv$.Random.seed <- current.seed
  
  #Data generation

  #Generate values of the running variable, with half below the cutoff and half above

  z <- rbeta(n, shape1=alpha, shape2=beta)
  x <- 2*z-1

  #Generate error values
  e <- rnorm(n,0,.1295)

  #examining the distribution of x
  count.below <- sum(x<c)
  count.above <- sum(x>=c)
  h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
  rot.below <- sum(x<c & x>(c-h.rot))
  rot.above <- sum(x>c & x<(c+h.rot))

  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                       rot.below=rot.below, rot.above=rot.above)))

  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn",
                   "design4.fn", "design5.fn", "design6.fn", "design7.fn"),
                 function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)

  #Get estimates
  tauhat.design1 <- tau.hat.one.y.fn(y.vec=ymat[,1], x.vec=x, c=0)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")

  tauhat.design2 <- tau.hat.one.y.fn(y.vec=ymat[,2], x.vec=x, c=0)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")

  tauhat.design3 <- tau.hat.one.y.fn(y.vec=ymat[,3], x.vec=x, c=0)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")

  tauhat.design4 <- tau.hat.one.y.fn(y.vec=ymat[,4], x.vec=x, c=0)
  colnames(tauhat.design4) <- paste(colnames(tauhat.design4), "design4", sep=".")

  tauhat.design5 <- tau.hat.one.y.fn(y.vec=ymat[,5], x.vec=x, c=0)
  colnames(tauhat.design5) <- paste(colnames(tauhat.design5), "design5", sep=".")

  tauhat.design6 <- tau.hat.one.y.fn(y.vec=ymat[,6], x.vec=x, c=0)
  colnames(tauhat.design6) <- paste(colnames(tauhat.design6), "design6", sep=".")

  tauhat.design7 <- tau.hat.one.y.fn(y.vec=ymat[,7], x.vec=x, c=0)
  colnames(tauhat.design7) <- paste(colnames(tauhat.design7), "design7", sep=".")

  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3,
                  tauhat.design4, tauhat.design5, tauhat.design6, tauhat.design7, eda)


  #Output data
  if (i==1) write.table(output, file="paper2comp031424rv2m27estimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="paper2comp031424rv2m27estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}

n <- 354 #m=44
set.seed(2485)
current.seed <- .Random.seed
for (i in 1:S){

  current.seed <- nextRNGStream(current.seed)
  .GlobalEnv$.Random.seed <- current.seed

  #Data generation

  #Generate values of the running variable, with half below the cutoff and half above

  z <- rbeta(n, shape1=alpha, shape2=beta)
  x <- 2*z-1

  #Generate error values
  e <- rnorm(n,0,.1295)

  #examining the distribution of x
  count.below <- sum(x<c)
  count.above <- sum(x>=c)
  h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
  rot.below <- sum(x<c & x>(c-h.rot))
  rot.above <- sum(x>c & x<(c+h.rot))

  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                       rot.below=rot.below, rot.above=rot.above)))

  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn",
                   "design4.fn", "design5.fn", "design6.fn", "design7.fn"),
                 function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)

  #Get estimates
  tauhat.design1 <- tau.hat.one.y.fn(y.vec=ymat[,1], x.vec=x, c=0)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")

  tauhat.design2 <- tau.hat.one.y.fn(y.vec=ymat[,2], x.vec=x, c=0)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")

  tauhat.design3 <- tau.hat.one.y.fn(y.vec=ymat[,3], x.vec=x, c=0)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")

  tauhat.design4 <- tau.hat.one.y.fn(y.vec=ymat[,4], x.vec=x, c=0)
  colnames(tauhat.design4) <- paste(colnames(tauhat.design4), "design4", sep=".")

  tauhat.design5 <- tau.hat.one.y.fn(y.vec=ymat[,5], x.vec=x, c=0)
  colnames(tauhat.design5) <- paste(colnames(tauhat.design5), "design5", sep=".")

  tauhat.design6 <- tau.hat.one.y.fn(y.vec=ymat[,6], x.vec=x, c=0)
  colnames(tauhat.design6) <- paste(colnames(tauhat.design6), "design6", sep=".")

  tauhat.design7 <- tau.hat.one.y.fn(y.vec=ymat[,7], x.vec=x, c=0)
  colnames(tauhat.design7) <- paste(colnames(tauhat.design7), "design7", sep=".")

  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3,
                  tauhat.design4, tauhat.design5, tauhat.design6, tauhat.design7, eda)


  #Output data
  if (i==1) write.table(output, file="paper2comp031424rv2m44estimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="paper2comp031424rv2m44estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}

n <- 490 #m=57
set.seed(6358)
current.seed <- .Random.seed

for (i in 1:S){

  current.seed <- nextRNGStream(current.seed)
  .GlobalEnv$.Random.seed <- current.seed
  
  #Data generation
  
  #Generate values of the running variable, with half below the cutoff and half above
  
  z <- rbeta(n, shape1=alpha, shape2=beta)
  x <- 2*z-1
  
  #Generate error values
  e <- rnorm(n,0,.1295)
  
  #examining the distribution of x
  count.below <- sum(x<c)
  count.above <- sum(x>=c)
  h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
  rot.below <- sum(x<c & x>(c-h.rot))
  rot.above <- sum(x>c & x<(c+h.rot))
  
  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                       rot.below=rot.below, rot.above=rot.above)))
  
  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn",
                   "design4.fn", "design5.fn", "design6.fn", "design7.fn"),
                 function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)
  
  #Get estimates
  tauhat.design1 <- tau.hat.one.y.fn(y.vec=ymat[,1], x.vec=x, c=0)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")
  
  tauhat.design2 <- tau.hat.one.y.fn(y.vec=ymat[,2], x.vec=x, c=0)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")
  
  tauhat.design3 <- tau.hat.one.y.fn(y.vec=ymat[,3], x.vec=x, c=0)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")
  
  tauhat.design4 <- tau.hat.one.y.fn(y.vec=ymat[,4], x.vec=x, c=0)
  colnames(tauhat.design4) <- paste(colnames(tauhat.design4), "design4", sep=".")
  
  tauhat.design5 <- tau.hat.one.y.fn(y.vec=ymat[,5], x.vec=x, c=0)
  colnames(tauhat.design5) <- paste(colnames(tauhat.design5), "design5", sep=".")
  
  tauhat.design6 <- tau.hat.one.y.fn(y.vec=ymat[,6], x.vec=x, c=0)
  colnames(tauhat.design6) <- paste(colnames(tauhat.design6), "design6", sep=".")
  
  tauhat.design7 <- tau.hat.one.y.fn(y.vec=ymat[,7], x.vec=x, c=0)
  colnames(tauhat.design7) <- paste(colnames(tauhat.design7), "design7", sep=".")
  
  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3,
                  tauhat.design4, tauhat.design5, tauhat.design6, tauhat.design7, eda)
  
  
  #Output data
  if (i==1) write.table(output, file="paper2comp031424rv2m57estimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="paper2comp031424rv2m57estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}

set.seed(7980) #initial seed
current.seed <- .Random.seed

n <- 56 #m=10
for (i in 1:S){

  current.seed <- nextRNGStream(current.seed)
  .GlobalEnv$.Random.seed <- current.seed

  #Data generation

  #Generate values of the running variable, with half below the cutoff and half above

  z <- rbeta(n, shape1=alpha, shape2=beta)
  x <- 2*z-1

  #Generate error values
  e <- rnorm(n,0,.1295)

  #examining the distribution of x
  count.below <- sum(x<c)
  count.above <- sum(x>=c)
  h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
  rot.below <- sum(x<c & x>(c-h.rot))
  rot.above <- sum(x>c & x<(c+h.rot))

  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                       rot.below=rot.below, rot.above=rot.above)))

  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn"),
                 function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)

  #Get estimates
  tauhat.design1 <- tau.hat.sens(y.vec=ymat[,1], x.vec=x, c=0, dbound=2)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")

  tauhat.design2 <- tau.hat.sens(y.vec=ymat[,2], x.vec=x, c=0, dbound=233.26)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")

  tauhat.design3 <- tau.hat.sens(y.vec=ymat[,3], x.vec=x, c=0, dbound=16.2)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")

  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3, eda)


  #Output data
  if (i==1) write.table(output, file="paper2comp031424rv2m10sestimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="paper2comp031424rv2m10sestimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}

n <- 140 #m=21
set.seed(4620)
current.seed <- .Random.seed
for (i in 1:S){
  current.seed <- nextRNGStream(current.seed)
  .GlobalEnv$.Random.seed <- current.seed

  #Data generation

  #Generate values of the running variable, with half below the cutoff and half above

  z <- rbeta(n, shape1=alpha, shape2=beta)
  x <- 2*z-1

  #Generate error values
  e <- rnorm(n,0,.1295)

  #examining the distribution of x
  count.below <- sum(x<c)
  count.above <- sum(x>=c)
  h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
  rot.below <- sum(x<c & x>(c-h.rot))
  rot.above <- sum(x>c & x<(c+h.rot))

  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                       rot.below=rot.below, rot.above=rot.above)))

  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn"),
                 function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)

  #Get estimates
  tauhat.design1 <- tau.hat.sens(y.vec=ymat[,1], x.vec=x, c=0, dbound=2)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")

  tauhat.design2 <- tau.hat.sens(y.vec=ymat[,2], x.vec=x, c=0, dbound=233.26)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")

  tauhat.design3 <- tau.hat.sens(y.vec=ymat[,3], x.vec=x, c=0, dbound=16.2)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")

  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3, eda)


  #Output data
  if (i==1) write.table(output, file="paper2comp031424rv2m21sestimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="paper2comp031424rv2m21sestimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}


n <- 194 #m=27
set.seed(8766)
current.seed <- .Random.seed
for (i in 1:S){

  current.seed <- nextRNGStream(current.seed)
  .GlobalEnv$.Random.seed <- current.seed

  #Data generation

  #Generate values of the running variable, with half below the cutoff and half above

  z <- rbeta(n, shape1=alpha, shape2=beta)
  x <- 2*z-1

  #Generate error values
  e <- rnorm(n,0,.1295)

  #examining the distribution of x
  count.below <- sum(x<c)
  count.above <- sum(x>=c)
  h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
  rot.below <- sum(x<c & x>(c-h.rot))
  rot.above <- sum(x>c & x<(c+h.rot))

  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                       rot.below=rot.below, rot.above=rot.above)))

  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn"),
                 function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)

  #Get estimates
  tauhat.design1 <- tau.hat.sens(y.vec=ymat[,1], x.vec=x, c=0, dbound=2)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")

  tauhat.design2 <- tau.hat.sens(y.vec=ymat[,2], x.vec=x, c=0, dbound=233.26)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")

  tauhat.design3 <- tau.hat.sens(y.vec=ymat[,3], x.vec=x, c=0, dbound=16.2)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")

  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3, eda)
  #Output data
  if (i==1) write.table(output, file="paper2comp031424rv2m27sestimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="paper2comp031424rv2m27sestimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}

n <- 354 #m=44
set.seed(2485)
current.seed <- .Random.seed
for (i in 1:S){

  current.seed <- nextRNGStream(current.seed)
  .GlobalEnv$.Random.seed <- current.seed

  #Data generation

  #Generate values of the running variable, with half below the cutoff and half above

  z <- rbeta(n, shape1=alpha, shape2=beta)
  x <- 2*z-1

  #Generate error values
  e <- rnorm(n,0,.1295)

  #examining the distribution of x
  count.below <- sum(x<c)
  count.above <- sum(x>=c)
  h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
  rot.below <- sum(x<c & x>(c-h.rot))
  rot.above <- sum(x>c & x<(c+h.rot))

  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                       rot.below=rot.below, rot.above=rot.above)))

  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn"),
                 function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)

  #Get estimates
  tauhat.design1 <- tau.hat.sens(y.vec=ymat[,1], x.vec=x, c=0, dbound=2)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")

  tauhat.design2 <- tau.hat.sens(y.vec=ymat[,2], x.vec=x, c=0, dbound=233.26)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")

  tauhat.design3 <- tau.hat.sens(y.vec=ymat[,3], x.vec=x, c=0, dbound=16.2)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")

  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3, eda)

  #Output data
  if (i==1) write.table(output, file="paper2comp031424rv2m44sestimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="paper2comp031424rv2m44sestimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}

n <- 490 #m=57
set.seed(6358)
current.seed <- .Random.seed

for (i in 1:S){

  current.seed <- nextRNGStream(current.seed)
  .GlobalEnv$.Random.seed <- current.seed

  #Data generation

  #Generate values of the running variable, with half below the cutoff and half above

  z <- rbeta(n, shape1=alpha, shape2=beta)
  x <- 2*z-1

  #Generate error values
  e <- rnorm(n,0,.1295)

  #examining the distribution of x
  count.below <- sum(x<c)
  count.above <- sum(x>=c)
  h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
  rot.below <- sum(x<c & x>(c-h.rot))
  rot.above <- sum(x>c & x<(c+h.rot))

  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                       rot.below=rot.below, rot.above=rot.above)))

  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn"),
                 function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)

  #Get estimates
  tauhat.design1 <- tau.hat.sens(y.vec=ymat[,1], x.vec=x, c=0, dbound=2)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")

  tauhat.design2 <- tau.hat.sens(y.vec=ymat[,2], x.vec=x, c=0, dbound=233.26)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")

  tauhat.design3 <- tau.hat.sens(y.vec=ymat[,3], x.vec=x, c=0, dbound=16.2)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")

  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3, eda)

  #Output data
  if (i==1) write.table(output, file="paper2comp031424rv2m57sestimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="paper2comp031424rv2m57sestimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}


# setup for RV3 ====
S=50000 # number of iterations of the simulation
c=0 #cutoff for purposes of eda below

# Data generation parameters
alpha <- 14
beta <- 7

#setting up the random number generator
RNGkind("L'Ecuyer-CMRG")

set.seed(7300) #initial seed
current.seed <- .Random.seed

n <- 140 #m=10

for (i in 1:S){

    current.seed <- nextRNGStream(current.seed)
    .GlobalEnv$.Random.seed <- current.seed

  #Data generation

  #Generate values of the running variable, with half below the cutoff and half above

  z <- rbeta(n, shape1=alpha, shape2=beta)
  x <- 2*z-1

  #Generate error values
  e <- rnorm(n,0,.1295)

  #examining the distribution of x
  count.below <- sum(x<c)
  count.above <- sum(x>=c)
  h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
  rot.below <- sum(x<c & x>(c-h.rot))
  rot.above <- sum(x>c & x<(c+h.rot))

  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                       rot.below=rot.below, rot.above=rot.above)))

  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn",
                   "design4.fn", "design5.fn", "design6.fn", "design7.fn"),
                 function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)

  #Get estimates
  tauhat.design1 <- tau.hat.one.y.fn(y.vec=ymat[,1], x.vec=x, c=0)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")

  tauhat.design2 <- tau.hat.one.y.fn(y.vec=ymat[,2], x.vec=x, c=0)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")

  tauhat.design3 <- tau.hat.one.y.fn(y.vec=ymat[,3], x.vec=x, c=0)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")

  tauhat.design4 <- tau.hat.one.y.fn(y.vec=ymat[,4], x.vec=x, c=0)
  colnames(tauhat.design4) <- paste(colnames(tauhat.design4), "design4", sep=".")

  tauhat.design5 <- tau.hat.one.y.fn(y.vec=ymat[,5], x.vec=x, c=0)
  colnames(tauhat.design5) <- paste(colnames(tauhat.design5), "design5", sep=".")

  tauhat.design6 <- tau.hat.one.y.fn(y.vec=ymat[,6], x.vec=x, c=0)
  colnames(tauhat.design6) <- paste(colnames(tauhat.design6), "design6", sep=".")

  tauhat.design7 <- tau.hat.one.y.fn(y.vec=ymat[,7], x.vec=x, c=0)
  colnames(tauhat.design7) <- paste(colnames(tauhat.design7), "design7", sep=".")

  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3,
                  tauhat.design4, tauhat.design5, tauhat.design6, tauhat.design7, eda)


  #Output data
  if (i==1) write.table(output, file="paper2comp031424rv3m10estimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="paper2comp031424rv3m10estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}

n <- 354 #m=21
set.seed(8242)
current.seed <- .Random.seed
for (i in 1:S){

      current.seed <- nextRNGStream(current.seed)
      .GlobalEnv$.Random.seed <- current.seed

#Data generation

#Generate values of the running variable, with half below the cutoff and half above

z <- rbeta(n, shape1=alpha, shape2=beta)
x <- 2*z-1

#Generate error values
e <- rnorm(n,0,.1295)

#examining the distribution of x
count.below <- sum(x<c)
count.above <- sum(x>=c)
h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
rot.below <- sum(x<c & x>(c-h.rot))
rot.above <- sum(x>c & x<(c+h.rot))

eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                     rot.below=rot.below, rot.above=rot.above)))

#Generate y matrix, one column for each function
ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn",
                 "design4.fn", "design5.fn", "design6.fn", "design7.fn"),
               function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)

#Get estimates
tauhat.design1 <- tau.hat.one.y.fn(y.vec=ymat[,1], x.vec=x, c=0)
colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")

tauhat.design2 <- tau.hat.one.y.fn(y.vec=ymat[,2], x.vec=x, c=0)
colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")

tauhat.design3 <- tau.hat.one.y.fn(y.vec=ymat[,3], x.vec=x, c=0)
colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")

tauhat.design4 <- tau.hat.one.y.fn(y.vec=ymat[,4], x.vec=x, c=0)
colnames(tauhat.design4) <- paste(colnames(tauhat.design4), "design4", sep=".")

tauhat.design5 <- tau.hat.one.y.fn(y.vec=ymat[,5], x.vec=x, c=0)
colnames(tauhat.design5) <- paste(colnames(tauhat.design5), "design5", sep=".")

tauhat.design6 <- tau.hat.one.y.fn(y.vec=ymat[,6], x.vec=x, c=0)
colnames(tauhat.design6) <- paste(colnames(tauhat.design6), "design6", sep=".")

tauhat.design7 <- tau.hat.one.y.fn(y.vec=ymat[,7], x.vec=x, c=0)
colnames(tauhat.design7) <- paste(colnames(tauhat.design7), "design7", sep=".")

iteration <- i
output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3,
                tauhat.design4, tauhat.design5, tauhat.design6, tauhat.design7, eda)


#Output data
if (i==1) write.table(output, file="paper2comp031424rv3m21estimates.csv", sep=",", row.names=FALSE)
else write.table(output, file="paper2comp031424rv3m21estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}


n <- 494 #m=27
set.seed(501)
current.seed <- .Random.seed
for (i in 1:S){

 
      current.seed <- nextRNGStream(current.seed)
      .GlobalEnv$.Random.seed <- current.seed
  #Data generation

  #Generate values of the running variable, with half below the cutoff and half above

  z <- rbeta(n, shape1=alpha, shape2=beta)
  x <- 2*z-1

  #Generate error values
  e <- rnorm(n,0,.1295)

  #examining the distribution of x
  count.below <- sum(x<c)
  count.above <- sum(x>=c)
  h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
  rot.below <- sum(x<c & x>(c-h.rot))
  rot.above <- sum(x>c & x<(c+h.rot))

  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                       rot.below=rot.below, rot.above=rot.above)))

  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn",
                   "design4.fn", "design5.fn", "design6.fn", "design7.fn"),
                 function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)

  #Get estimates
  tauhat.design1 <- tau.hat.one.y.fn(y.vec=ymat[,1], x.vec=x, c=0)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")

  tauhat.design2 <- tau.hat.one.y.fn(y.vec=ymat[,2], x.vec=x, c=0)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")

  tauhat.design3 <- tau.hat.one.y.fn(y.vec=ymat[,3], x.vec=x, c=0)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")

  tauhat.design4 <- tau.hat.one.y.fn(y.vec=ymat[,4], x.vec=x, c=0)
  colnames(tauhat.design4) <- paste(colnames(tauhat.design4), "design4", sep=".")

  tauhat.design5 <- tau.hat.one.y.fn(y.vec=ymat[,5], x.vec=x, c=0)
  colnames(tauhat.design5) <- paste(colnames(tauhat.design5), "design5", sep=".")

  tauhat.design6 <- tau.hat.one.y.fn(y.vec=ymat[,6], x.vec=x, c=0)
  colnames(tauhat.design6) <- paste(colnames(tauhat.design6), "design6", sep=".")

  tauhat.design7 <- tau.hat.one.y.fn(y.vec=ymat[,7], x.vec=x, c=0)
  colnames(tauhat.design7) <- paste(colnames(tauhat.design7), "design7", sep=".")

  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3,
                  tauhat.design4, tauhat.design5, tauhat.design6, tauhat.design7, eda)


  #Output data
  if (i==1) write.table(output, file="paper2comp031424rv3m27estimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="paper2comp031424rv3m27estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}

n <- 905 #m=44
set.seed(4474)
current.seed <- .Random.seed
for (i in 1:S){
  
  current.seed <- nextRNGStream(current.seed)
  .GlobalEnv$.Random.seed <- current.seed
  #Data generation
  
  #Generate values of the running variable, with half below the cutoff and half above
  
  z <- rbeta(n, shape1=alpha, shape2=beta)
  x <- 2*z-1
  
  #Generate error values
  e <- rnorm(n,0,.1295)
  
  #examining the distribution of x
  count.below <- sum(x<c)
  count.above <- sum(x>=c)
  h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
  rot.below <- sum(x<c & x>(c-h.rot))
  rot.above <- sum(x>c & x<(c+h.rot))
  
  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                       rot.below=rot.below, rot.above=rot.above)))
  
  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn",
                   "design4.fn", "design5.fn", "design6.fn", "design7.fn"),
                 function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)
  
  #Get estimates
  tauhat.design1 <- tau.hat.one.y.fn(y.vec=ymat[,1], x.vec=x, c=0)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")
  
  tauhat.design2 <- tau.hat.one.y.fn(y.vec=ymat[,2], x.vec=x, c=0)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")
  
  tauhat.design3 <- tau.hat.one.y.fn(y.vec=ymat[,3], x.vec=x, c=0)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")
  
  tauhat.design4 <- tau.hat.one.y.fn(y.vec=ymat[,4], x.vec=x, c=0)
  colnames(tauhat.design4) <- paste(colnames(tauhat.design4), "design4", sep=".")
  
  tauhat.design5 <- tau.hat.one.y.fn(y.vec=ymat[,5], x.vec=x, c=0)
  colnames(tauhat.design5) <- paste(colnames(tauhat.design5), "design5", sep=".")
  
  tauhat.design6 <- tau.hat.one.y.fn(y.vec=ymat[,6], x.vec=x, c=0)
  colnames(tauhat.design6) <- paste(colnames(tauhat.design6), "design6", sep=".")
  
  tauhat.design7 <- tau.hat.one.y.fn(y.vec=ymat[,7], x.vec=x, c=0)
  colnames(tauhat.design7) <- paste(colnames(tauhat.design7), "design7", sep=".")
  
  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3,
                  tauhat.design4, tauhat.design5, tauhat.design6, tauhat.design7, eda)
  
  
  #Output data
  if (i==1) write.table(output, file="paper2comp031424rv3m44estimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="paper2comp031424rv3m44estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}

n <- 1254 #m=57
set.seed(2249)
current.seed <- .Random.seed

for (i in 1:S){

  current.seed <- nextRNGStream(current.seed)
  .GlobalEnv$.Random.seed <- current.seed

#Data generation

#Generate values of the running variable, with half below the cutoff and half above

z <- rbeta(n, shape1=alpha, shape2=beta)
x <- 2*z-1

#Generate error values
e <- rnorm(n,0,.1295)

#examining the distribution of x
count.below <- sum(x<c)
count.above <- sum(x>=c)
h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
rot.below <- sum(x<c & x>(c-h.rot))
rot.above <- sum(x>c & x<(c+h.rot))

eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                     rot.below=rot.below, rot.above=rot.above)))

#Generate y matrix, one column for each function
ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn",
                 "design4.fn", "design5.fn", "design6.fn", "design7.fn"),
               function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)

#Get estimates
tauhat.design1 <- tau.hat.one.y.fn(y.vec=ymat[,1], x.vec=x, c=0)
colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")

tauhat.design2 <- tau.hat.one.y.fn(y.vec=ymat[,2], x.vec=x, c=0)
colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")

tauhat.design3 <- tau.hat.one.y.fn(y.vec=ymat[,3], x.vec=x, c=0)
colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")

tauhat.design4 <- tau.hat.one.y.fn(y.vec=ymat[,4], x.vec=x, c=0)
colnames(tauhat.design4) <- paste(colnames(tauhat.design4), "design4", sep=".")

tauhat.design5 <- tau.hat.one.y.fn(y.vec=ymat[,5], x.vec=x, c=0)
colnames(tauhat.design5) <- paste(colnames(tauhat.design5), "design5", sep=".")

tauhat.design6 <- tau.hat.one.y.fn(y.vec=ymat[,6], x.vec=x, c=0)
colnames(tauhat.design6) <- paste(colnames(tauhat.design6), "design6", sep=".")

tauhat.design7 <- tau.hat.one.y.fn(y.vec=ymat[,7], x.vec=x, c=0)
colnames(tauhat.design7) <- paste(colnames(tauhat.design7), "design7", sep=".")

iteration <- i
output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3,
                tauhat.design4, tauhat.design5, tauhat.design6, tauhat.design7, eda)


#Output data
if (i==1) write.table(output, file="paper2comp031424rv3m57estimates.csv", sep=",", row.names=FALSE)
else write.table(output, file="paper2comp031424rv3m57estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}


n <- 140 #m=10

for (i in 1:S){

  current.seed <- nextRNGStream(current.seed)
  .GlobalEnv$.Random.seed <- current.seed

  #Data generation

  #Generate values of the running variable, with half below the cutoff and half above

  z <- rbeta(n, shape1=alpha, shape2=beta)
  x <- 2*z-1

  #Generate error values
  e <- rnorm(n,0,.1295)

  #examining the distribution of x
  count.below <- sum(x<c)
  count.above <- sum(x>=c)
  h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
  rot.below <- sum(x<c & x>(c-h.rot))
  rot.above <- sum(x>c & x<(c+h.rot))

  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                       rot.below=rot.below, rot.above=rot.above)))

  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn"),
                 function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)

  #Get estimates
  tauhat.design1 <- tau.hat.sens(y.vec=ymat[,1], x.vec=x, c=0, dbound=2)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")

  tauhat.design2 <- tau.hat.sens(y.vec=ymat[,2], x.vec=x, c=0, dbound=233.26)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")

  tauhat.design3 <- tau.hat.sens(y.vec=ymat[,3], x.vec=x, c=0, dbound=16.2)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")

  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3, eda)


  #Output data
  if (i==1) write.table(output, file="paper2comp031424rv3m10sestimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="paper2comp031424rv3m10sestimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}

n <- 354 #m=21
set.seed(8242)
current.seed <- .Random.seed
for (i in 1:S){
  current.seed <- nextRNGStream(current.seed)
  .GlobalEnv$.Random.seed <- current.seed

  #Data generation

  #Generate values of the running variable, with half below the cutoff and half above

  z <- rbeta(n, shape1=alpha, shape2=beta)
  x <- 2*z-1

  #Generate error values
  e <- rnorm(n,0,.1295)

  #examining the distribution of x
  count.below <- sum(x<c)
  count.above <- sum(x>=c)
  h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
  rot.below <- sum(x<c & x>(c-h.rot))
  rot.above <- sum(x>c & x<(c+h.rot))

  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                       rot.below=rot.below, rot.above=rot.above)))

  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn"),
                 function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)

  #Get estimates
  tauhat.design1 <- tau.hat.sens(y.vec=ymat[,1], x.vec=x, c=0, dbound=2)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")

  tauhat.design2 <- tau.hat.sens(y.vec=ymat[,2], x.vec=x, c=0, dbound=233.26)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")

  tauhat.design3 <- tau.hat.sens(y.vec=ymat[,3], x.vec=x, c=0, dbound=16.2)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")

  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3, eda)

  #Output data
  if (i==1) write.table(output, file="paper2comp031424rv3m21sestimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="paper2comp031424rv3m21sestimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}


n <- 494 #m=27
set.seed(501)
current.seed <- .Random.seed
for (i in 1:S){

  current.seed <- nextRNGStream(current.seed)
  .GlobalEnv$.Random.seed <- current.seed

  #Data generation

  #Generate values of the running variable, with half below the cutoff and half above

  z <- rbeta(n, shape1=alpha, shape2=beta)
  x <- 2*z-1

  #Generate error values
  e <- rnorm(n,0,.1295)

  #examining the distribution of x
  count.below <- sum(x<c)
  count.above <- sum(x>=c)
  h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
  rot.below <- sum(x<c & x>(c-h.rot))
  rot.above <- sum(x>c & x<(c+h.rot))

  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                       rot.below=rot.below, rot.above=rot.above)))

  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn"),
                 function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)

  #Get estimates
  tauhat.design1 <- tau.hat.sens(y.vec=ymat[,1], x.vec=x, c=0, dbound=2)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")

  tauhat.design2 <- tau.hat.sens(y.vec=ymat[,2], x.vec=x, c=0, dbound=233.26)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")

  tauhat.design3 <- tau.hat.sens(y.vec=ymat[,3], x.vec=x, c=0, dbound=16.2)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")

  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3, eda)

  #Output data
  if (i==1) write.table(output, file="paper2comp031424rv3m27sestimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="paper2comp031424rv3m27sestimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}

n <- 905 #m=44
set.seed(4474)
current.seed <- .Random.seed
for (i in 1:S){

  current.seed <- nextRNGStream(current.seed)
  .GlobalEnv$.Random.seed <- current.seed
  #Data generation

  #Generate values of the running variable, with half below the cutoff and half above

  z <- rbeta(n, shape1=alpha, shape2=beta)
  x <- 2*z-1

  #Generate error values
  e <- rnorm(n,0,.1295)

  #examining the distribution of x
  count.below <- sum(x<c)
  count.above <- sum(x>=c)
  h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
  rot.below <- sum(x<c & x>(c-h.rot))
  rot.above <- sum(x>c & x<(c+h.rot))

  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                       rot.below=rot.below, rot.above=rot.above)))

  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn"),
                 function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)

  #Get estimates
  tauhat.design1 <- tau.hat.sens(y.vec=ymat[,1], x.vec=x, c=0, dbound=2)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")

  tauhat.design2 <- tau.hat.sens(y.vec=ymat[,2], x.vec=x, c=0, dbound=233.26)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")

  tauhat.design3 <- tau.hat.sens(y.vec=ymat[,3], x.vec=x, c=0, dbound=16.2)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")

  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3, eda)
  #Output data
  if (i==1) write.table(output, file="paper2comp031424rv3m44sestimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="paper2comp031424rv3m44sestimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}

n <- 1254 #m=57
set.seed(2249)
current.seed <- .Random.seed

for (i in 1:S){

  current.seed <- nextRNGStream(current.seed)
  .GlobalEnv$.Random.seed <- current.seed

  #Data generation

  #Generate values of the running variable, with half below the cutoff and half above

  z <- rbeta(n, shape1=alpha, shape2=beta)
  x <- 2*z-1

  #Generate error values
  e <- rnorm(n,0,.1295)

  #examining the distribution of x
  count.below <- sum(x<c)
  count.above <- sum(x>=c)
  h.rot <- 0.9*pmin(sd(x), IQR(x)/1.34)*n^(-.2)
  rot.below <- sum(x<c & x>(c-h.rot))
  rot.above <- sum(x>c & x<(c+h.rot))

  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot,
                       rot.below=rot.below, rot.above=rot.above)))

  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn"),
                 function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)

  #Get estimates
  tauhat.design1 <- tau.hat.sens(y.vec=ymat[,1], x.vec=x, c=0, dbound=2)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")

  tauhat.design2 <- tau.hat.sens(y.vec=ymat[,2], x.vec=x, c=0, dbound=233.26)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")

  tauhat.design3 <- tau.hat.sens(y.vec=ymat[,3], x.vec=x, c=0, dbound=16.2)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")

  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3, eda)

  #Output data
  if (i==1) write.table(output, file="paper2comp031424rv3m57sestimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="paper2comp031424rv3m57sestimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
}

