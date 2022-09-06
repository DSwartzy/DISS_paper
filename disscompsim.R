#This file contains the code for the simulation comparing small sample performances of RD estimators using DISS


# required packages ====
library(rdrobust)
library(dbscan)
library(ks)
library(RDHonest)

# true functions ====
#Data generation functions

#this is the design based on the Indiana data
design1.fn <- function(x){
  (.15-.15*x+2.5*x^2-1.5*x^3)*(x>=0) + (.05+1.5*x+3.2*x^2+2.7*x^3)*(x<0) + e
}

#this is a modified version of desgin 3 from AK
design3.fn <- function(x){
  (x+1)^2-2*pmax(x+.2,0)^2+2*pmax(x-.2,0)^2-2*pmax(x-.4,0)^2+2*pmax(x-.7,0)^2-.92+0.1*(x>=0) + e
}

#this is design 3 from from IK
design2.fn <- function(x) {
  (0.52+0.84*x-3.0*x^2+7.99*x^3-9.01*x^4+3.56*x^5)*(x>=0)+(0.42+0.84*x-3.0*x^2+7.99*x^3-9.01*x^4+3.56*x^5)*(x<0) + e
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
  h.mseik1 <- NA
  h.msecct0 <- NA
  h.msecct1 <- NA
  h.msecct2 <- NA
  h.cerccf0 <- NA
  h.cerccf1 <- NA
  h.cerccf2 <- NA
  h.mseak1 <- NA
  h.mseak2 <- NA
  b.msecct0 <- NA
  b.msecct1 <- NA
  b.msecct2 <- NA
  b.cerccf0 <- NA
  b.cerccf1 <- NA
  b.cerccf2 <- NA
  
  #added part for the ak bandwidths using mhat
  h.mseakh1 <- NA
  h.mseakh2 <- NA
  
  #now attempting the bandwidth values from rdrobust and rdhonest
  
  try({
    invisible(capture.output(fit.ik1 <- rdbwselect_2014(y.vec, x.vec, c=c, p=1, kernel=cctkernel, bwselect="IK")))
    h.mseik1 <- unname(fit.ik1$bws[1,1]) 
  }, silent=TRUE)

  try({
    invisible(capture.output(fit.cct0 <- rdbwselect(y.vec, x.vec, c=c, p=0, kernel=cctkernel)))
    h.msecct0 <- fit.cct0$bws[1,1]
    b.msecct0 <- fit.cct0$bws[1,3]
  }, silent=TRUE)
  
  try({
    invisible(capture.output(fit.cct1 <- rdbwselect(y.vec, x.vec, c=c, p=1, kernel=cctkernel)))
    h.msecct1 <- fit.cct1$bws[1,1] 
    b.msecct1 <- fit.cct1$bws[1,3]
  }, silent=TRUE)
  
  try({
    invisible(capture.output(fit.cct2 <- rdbwselect(y.vec, x.vec, c=c, p=2, kernel=cctkernel)))
    h.msecct2 <- fit.cct2$bws[1,1] 
    b.msecct2 <- fit.cct2$bws[1,3]
  }, silent=TRUE)
  
  try({
    invisible(capture.output(fit.ccf0 <- rdbwselect(y.vec, x.vec, c=c, p=0, kernel=cctkernel, bwselect="cerrd")))
    h.cerccf0 <- fit.ccf0$bws[1,1]
    b.cerccf0 <- fit.ccf0$bws[1,3]
  }, silent=TRUE)
  
  try({
    invisible(capture.output(fit.ccf1 <- rdbwselect(y.vec, x.vec, c=c, p=1, kernel=cctkernel, bwselect="cerrd")))
    h.cerccf1 <- fit.ccf1$bws[1,1] 
    b.cerccf1 <- fit.ccf1$bws[1,3]
  }, silent=TRUE)
  
  try({
    invisible(capture.output(fit.ccf2 <- rdbwselect(y.vec, x.vec, c=c, p=2, kernel=cctkernel, bwselect="cerrd")))
    h.cerccf2 <- fit.ccf2$bws[1,1] 
    b.cerccf2 <- fit.ccf2$bws[1,3]
  }, silent=TRUE)

  try({
    invisible(capture.output(fit.ak1 <- RDOptBW(y.vec~x.vec, cutoff=c, M=dbound, order=1, opt.criterion = "MSE")))
    h.mseak1 <- unname(fit.ak1$h[1])
  }, silent=TRUE)
  
  try({
    invisible(capture.output(fit.ak2 <- RDOptBW(y.vec~x.vec, cutoff=c, M=dbound, order=2, opt.criterion = "MSE")))
    h.mseak2 <- unname(fit.ak2$h[1])
  }, silent=TRUE)
  
  try({
    invisible(capture.output(fit.akh1 <- RDOptBW(y.vec~x.vec, cutoff=c, M=mhat, order=1, opt.criterion = "MSE")))
    h.mseakh1 <- unname(fit.akh1$h[1])
  }, silent=TRUE)
  
  try({
    invisible(capture.output(fit.akh2 <- RDOptBW(y.vec~x.vec, cutoff=c, M=mhat, order=2, opt.criterion = "MSE")))
    h.mseakh2 <- unname(fit.akh2$h[1])
  }, silent=TRUE)
  
  #getting estimators for each of the bandwidth values above, with default NAs if the try block has an error
  
  rob.mseik1 <- t(as.matrix(c(tau.us=NA_real_, tau.bc=NA_real_, se.us=NA_real_, se.rb=NA_real_)))
  rob.mseik1 <- t(as.matrix(c(tau.us=NA_real_, tau.bc=NA_real_, se.us=NA_real_, se.rb=NA_real_)))
  rob.msecct0 <- t(as.matrix(c(tau.us=NA_real_, tau.bc=NA_real_, se.us=NA_real_, se.rb=NA_real_)))
  rob.msecct1 <- t(as.matrix(c(tau.us=NA_real_, tau.bc=NA_real_, se.us=NA_real_, se.rb=NA_real_)))
  rob.msecct2 <- t(as.matrix(c(tau.us=NA_real_, tau.bc=NA_real_, se.us=NA_real_, se.rb=NA_real_)))
  rob.cerccf0 <- t(as.matrix(c(tau.us=NA_real_, tau.bc=NA_real_, se.us=NA_real_, se.rb=NA_real_)))
  rob.cerccf1 <- t(as.matrix(c(tau.us=NA_real_, tau.bc=NA_real_, se.us=NA_real_, se.rb=NA_real_)))
  rob.cerccf2 <- t(as.matrix(c(tau.us=NA_real_, tau.bc=NA_real_, se.us=NA_real_, se.rb=NA_real_)))
  rob.mseak1 <- t(as.matrix(c(tau.us=NA_real_, tau.bc=NA_real_, se.us=NA_real_, se.rb=NA_real_)))
  rob.mseak2 <- t(as.matrix(c(tau.us=NA_real_, tau.bc=NA_real_, se.us=NA_real_, se.rb=NA_real_)))
  rob.mseakh1 <- t(as.matrix(c(tau.us=NA_real_, tau.bc=NA_real_, se.us=NA_real_, se.rb=NA_real_)))
  rob.mseakh2 <- t(as.matrix(c(tau.us=NA_real_, tau.bc=NA_real_, se.us=NA_real_, se.rb=NA_real_)))
  
  hmod.mseik1 <- data.frame(estimate=NA_real_, sd=NA_real_, hl=NA_real_)
  hmod.msecct1 <-data.frame(estimate=NA_real_, sd=NA_real_, hl=NA_real_)
  hmod.msecct2 <- data.frame(estimate=NA_real_, sd=NA_real_, hl=NA_real_)
  hmod.cerccf1 <- data.frame(estimate=NA_real_, sd=NA_real_, hl=NA_real_)
  hmod.cerccf2 <- data.frame(estimate=NA_real_, sd=NA_real_, hl=NA_real_)
  hmod.mseak1 <- data.frame(estimate=NA_real_, sd=NA_real_, hl=NA_real_)
  hmod.mseak2 <- data.frame(estimate=NA_real_, sd=NA_real_, hl=NA_real_)
  hmod.mseakh1 <- data.frame(estimate=NA_real_, sd=NA_real_, hl=NA_real_)
  hmod.mseakh2 <- data.frame(estimate=NA_real_, sd=NA_real_, hl=NA_real_)
  
  hmod2.mseik1 <- data.frame(estimate=NA_real_, sd=NA_real_, hl=NA_real_)
  hmod2.msecct1 <-data.frame(estimate=NA_real_, sd=NA_real_, hl=NA_real_)
  hmod2.msecct2 <- data.frame(estimate=NA_real_, sd=NA_real_, hl=NA_real_)
  hmod2.cerccf1 <- data.frame(estimate=NA_real_, sd=NA_real_, hl=NA_real_)
  hmod2.cerccf2 <- data.frame(estimate=NA_real_, sd=NA_real_, hl=NA_real_)
  hmod2.mseak1 <- data.frame(estimate=NA_real_, sd=NA_real_, hl=NA_real_)
  hmod2.mseak2 <- data.frame(estimate=NA_real_, sd=NA_real_, hl=NA_real_)
  hmod2.mseakh1 <- data.frame(estimate=NA_real_, sd=NA_real_, hl=NA_real_)
  hmod2.mseakh2 <- data.frame(estimate=NA_real_, sd=NA_real_, hl=NA_real_)
  
  #first the rdrobust estimators
  try({rob.mseik1 <- rdrobust(y=y.vec, x=x.vec, c=c, p=1, h=h.mseik1, kernel=cctkernel)$Estimate}, silent=TRUE)
  
  try({rob.msecct0 <- rdrobust(y=y.vec, x=x.vec, c=c, p=0, h=h.msecct0, b=b.msecct0, kernel=cctkernel)$Estimate}, silent=TRUE)
  try({rob.msecct1 <- rdrobust(y=y.vec, x=x.vec, c=c, p=1, h=h.msecct1, b=b.msecct1, kernel=cctkernel)$Estimate}, silent=TRUE)
  try({rob.msecct2 <- rdrobust(y=y.vec, x=x.vec, c=c, p=2, h=h.msecct2, b=b.msecct2, kernel=cctkernel)$Estimate}, silent=TRUE)
  
  try({rob.cerccf0 <- rdrobust(y=y.vec, x=x.vec, c=c, p=0, h=h.cerccf0, b=b.cerccf0, kernel=cctkernel)$Estimate}, silent=TRUE)
  try({rob.cerccf1 <- rdrobust(y=y.vec, x=x.vec, c=c, p=1, h=h.cerccf1, b=b.cerccf1, kernel=cctkernel)$Estimate}, silent=TRUE)
  try({rob.cerccf2 <- rdrobust(y=y.vec, x=x.vec, c=c, p=2, h=h.cerccf2, b=b.cerccf2, kernel=cctkernel)$Estimate}, silent=TRUE)
  
  try({rob.mseak1 <- rdrobust(y=y.vec, x=x.vec, c=c, p=1, h=h.mseak1, kernel=cctkernel)$Estimate}, silent=TRUE)
  try({rob.mseak2 <- rdrobust(y=y.vec, x=x.vec, c=c, p=2, h=h.mseak2, kernel=cctkernel)$Estimate}, silent=TRUE)
  
  try({rob.mseakh1 <- rdrobust(y=y.vec, x=x.vec, c=c, p=1, h=h.mseakh1, kernel=cctkernel)$Estimate}, silent=TRUE)
  try({rob.mseakh2 <- rdrobust(y=y.vec, x=x.vec, c=c, p=2, h=h.mseakh2, kernel=cctkernel)$Estimate}, silent=TRUE)
  
  #then for the rdhonest estimators
  try({hmod.mseik1 <- RDHonest(y.vec~x.vec, cutoff=c, M=dbound, order=1, h=h.mseik1)}, silent=TRUE)
  
  try({hmod.msecct1 <- RDHonest(y.vec~x.vec, cutoff=c, M=dbound, order=1, h=h.msecct1)}, silent=TRUE)
  try({hmod.msecct2 <- RDHonest(y.vec~x.vec, cutoff=c, M=dbound, order=2, h=h.msecct2)}, silent=TRUE)
  
  try({hmod.cerccf1 <- RDHonest(y.vec~x.vec, cutoff=c, M=dbound, order=1, h=h.cerccf1)}, silent=TRUE)
  try({hmod.cerccf2 <- RDHonest(y.vec~x.vec, cutoff=c, M=dbound, order=2, h=h.cerccf2)}, silent=TRUE)
  
  try({hmod.mseak1 <- RDHonest(y.vec~x.vec, cutoff=c, M=dbound, order=1, h=h.mseak1)}, silent=TRUE)
  try({hmod.mseak2 <- RDHonest(y.vec~x.vec, cutoff=c, M=dbound, order=2, h=h.mseak2)}, silent=TRUE)
  
  try({hmod.mseakh1 <- RDHonest(y.vec~x.vec, cutoff=c, M=dbound, order=1, h=h.mseakh1)}, silent=TRUE)
  try({hmod.mseakh2 <- RDHonest(y.vec~x.vec, cutoff=c, M=dbound, order=2, h=h.mseakh2)}, silent=TRUE)
  
  #now for the rdhonest estimators with M-hat
  try({hmod2.mseik1 <- RDHonest(y.vec~x.vec, cutoff=c, M=mhat, order=1, h=h.mseik1)}, silent=TRUE)
  
  try({hmod2.msecct1 <- RDHonest(y.vec~x.vec, cutoff=c, M=mhat, order=1, h=h.msecct1)}, silent=TRUE)
  try({hmod2.msecct2 <- RDHonest(y.vec~x.vec, cutoff=c, M=mhat, order=2, h=h.msecct2)}, silent=TRUE)
  
  try({hmod2.cerccf1 <- RDHonest(y.vec~x.vec, cutoff=c, M=mhat, order=1, h=h.cerccf1)}, silent=TRUE)
  try({hmod2.cerccf2 <- RDHonest(y.vec~x.vec, cutoff=c, M=mhat, order=2, h=h.cerccf2)}, silent=TRUE)
  
  try({hmod2.mseak1 <- RDHonest(y.vec~x.vec, cutoff=c, M=mhat, order=1, h=h.mseak1)}, silent=TRUE)
  try({hmod2.mseak2 <- RDHonest(y.vec~x.vec, cutoff=c, M=mhat, order=2, h=h.mseak2)}, silent=TRUE)
  
  try({hmod2.mseakh1 <- RDHonest(y.vec~x.vec, cutoff=c, M=mhat, order=1, h=h.mseakh1)}, silent=TRUE)
  try({hmod2.mseakh2 <- RDHonest(y.vec~x.vec, cutoff=c, M=mhat, order=2, h=h.mseakh2)}, silent=TRUE)

  

  
  #extract values and give names to estimates

  #for the rdrobust data
  colnames(rob.mseik1)[which(colnames(rob.mseik1) == "tau.bc")] <- "tau.rb"
  colnames(rob.msecct0)[which(colnames(rob.msecct0) == "tau.bc")] <- "tau.rb"
  colnames(rob.msecct1)[which(colnames(rob.msecct1) == "tau.bc")] <- "tau.rb"
  colnames(rob.msecct2)[which(colnames(rob.msecct2) == "tau.bc")] <- "tau.rb"
  colnames(rob.cerccf0)[which(colnames(rob.cerccf0) == "tau.bc")] <- "tau.rb"
  colnames(rob.cerccf1)[which(colnames(rob.cerccf1) == "tau.bc")] <- "tau.rb"
  colnames(rob.cerccf2)[which(colnames(rob.cerccf2) == "tau.bc")] <- "tau.rb"
  colnames(rob.mseak1)[which(colnames(rob.mseak1) == "tau.bc")] <- "tau.rb"
  colnames(rob.mseak2)[which(colnames(rob.mseak2) == "tau.bc")] <- "tau.rb"
  colnames(rob.mseakh1)[which(colnames(rob.mseakh1) == "tau.bc")] <- "tau.rb"
  colnames(rob.mseakh2)[which(colnames(rob.mseakh2) == "tau.bc")] <- "tau.rb"
  
  colnames(rob.mseik1) <- paste(colnames(rob.mseik1), "mseik1", sep=".") 
  colnames(rob.msecct0) <- paste(colnames(rob.msecct0), "msecct0", sep=".")
  colnames(rob.msecct1) <- paste(colnames(rob.msecct1), "msecct1", sep=".")
  colnames(rob.msecct2) <- paste(colnames(rob.msecct2), "msecct2", sep=".")
  colnames(rob.cerccf0) <- paste(colnames(rob.cerccf0), "cerccf0", sep=".")
  colnames(rob.cerccf1) <- paste(colnames(rob.cerccf1), "cerccf1", sep=".")
  colnames(rob.cerccf2) <- paste(colnames(rob.cerccf2), "cerccf2", sep=".")
  colnames(rob.mseak1) <- paste(colnames(rob.mseak1), "mseak1", sep=".")
  colnames(rob.mseak2) <- paste(colnames(rob.mseak2), "mseak2", sep=".")
  colnames(rob.mseakh1) <- paste(colnames(rob.mseakh1), "mseakh1", sep=".")
  colnames(rob.mseakh2) <- paste(colnames(rob.mseakh2), "mseakh2", sep=".")
  
  #for the rdhonest data
  hon.mseik1 <- t(as.matrix(c(hmod.mseik1$estimate, hmod.mseik1$sd, hmod.mseik1$hl)))
  colnames(hon.mseik1) <- c("tau.flci.mseik1", "se.flci.mseik1", "hl.flci.mseik1")
  
  hon.msecct1 <- t(as.matrix(c(hmod.msecct1$estimate, hmod.msecct1$sd, hmod.msecct1$hl)))
  colnames(hon.msecct1) <- c("tau.flci.msecct1", "se.flci.msecct1", "hl.flci.msecct1")
  hon.msecct2 <- t(as.matrix(c(hmod.msecct2$estimate, hmod.msecct2$sd, hmod.msecct2$hl)))
  colnames(hon.msecct2) <- c("tau.flci.msecct2", "se.flci.msecct2", "hl.flci.msecct2")
  
  hon.cerccf1 <- t(as.matrix(c(hmod.cerccf1$estimate, hmod.cerccf1$sd, hmod.cerccf1$hl)))
  colnames(hon.cerccf1) <- c("tau.flci.cerccf1", "se.flci.cerccf1", "hl.flci.cerccf1")
  hon.cerccf2 <- t(as.matrix(c(hmod.cerccf2$estimate, hmod.cerccf2$sd, hmod.cerccf2$hl)))
  colnames(hon.cerccf2) <- c("tau.flci.cerccf2", "se.flci.cerccf2", "hl.flci.cerccf2")
  
  hon.mseak1 <- t(as.matrix(c(hmod.mseak1$estimate, hmod.mseak1$sd, hmod.mseak1$hl)))
  colnames(hon.mseak1) <- c("tau.flci.mseak1", "se.flci.mseak1", "hl.flci.mseak1")
  hon.mseak2 <- t(as.matrix(c(hmod.mseak2$estimate, hmod.mseak2$sd, hmod.mseak2$hl)))
  colnames(hon.mseak2) <- c("tau.flci.mseak2", "se.flci.mseak2", "hl.flci.mseak2")
  
  hon.mseakh1 <- t(as.matrix(c(hmod.mseakh1$estimate, hmod.mseakh1$sd, hmod.mseakh1$hl)))
  colnames(hon.mseakh1) <- c("tau.flci.mseakh1", "se.flci.mseakh1", "hl.flci.mseakh1")
  hon.mseakh2 <- t(as.matrix(c(hmod.mseakh2$estimate, hmod.mseakh2$sd, hmod.mseakh2$hl)))
  colnames(hon.mseakh2) <- c("tau.flci.mseakh2", "se.flci.mseakh2", "hl.flci.mseakh2")
  
  #for the rdhonest data using mhat
  hon2.mseik1 <- t(as.matrix(c(hmod2.mseik1$estimate, hmod2.mseik1$sd, hmod2.mseik1$hl)))
  colnames(hon2.mseik1) <- c("tau.flcih.mseik1", "se.flcih.mseik1", "hl.flcih.mseik1")
  
  hon2.msecct1 <- t(as.matrix(c(hmod2.msecct1$estimate, hmod2.msecct1$sd, hmod2.msecct1$hl)))
  colnames(hon2.msecct1) <- c("tau.flcih.msecct1", "se.flcih.msecct1", "hl.flcih.msecct1")
  hon2.msecct2 <- t(as.matrix(c(hmod2.msecct2$estimate, hmod2.msecct2$sd, hmod2.msecct2$hl)))
  colnames(hon2.msecct2) <- c("tau.flcih.msecct2", "se.flcih.msecct2", "hl.flcih.msecct2")
  
  hon2.cerccf1 <- t(as.matrix(c(hmod2.cerccf1$estimate, hmod2.cerccf1$sd, hmod2.cerccf1$hl)))
  colnames(hon2.cerccf1) <- c("tau.flcih.cerccf1", "se.flcih.cerccf1", "hl.flcih.cerccf1")
  hon2.cerccf2 <- t(as.matrix(c(hmod2.cerccf2$estimate, hmod2.cerccf2$sd, hmod2.cerccf2$hl)))
  colnames(hon2.cerccf2) <- c("tau.flcih.cerccf2", "se.flcih.cerccf2", "hl.flcih.cerccf2")
  
  hon2.mseak1 <- t(as.matrix(c(hmod2.mseak1$estimate, hmod2.mseak1$sd, hmod2.mseak1$hl)))
  colnames(hon2.mseak1) <- c("tau.flcih.mseak1", "se.flcih.mseak1", "hl.flcih.mseak1")
  hon2.mseak2 <- t(as.matrix(c(hmod2.mseak2$estimate, hmod2.mseak2$sd, hmod2.mseak2$hl)))
  colnames(hon2.mseak2) <- c("tau.flcih.mseak2", "se.flcih.mseak2", "hl.flcih.mseak2")
  
  hon2.mseakh1 <- t(as.matrix(c(hmod2.mseakh1$estimate, hmod2.mseakh1$sd, hmod2.mseakh1$hl)))
  colnames(hon2.mseakh1) <- c("tau.flcih.mseakh1", "se.flcih.mseakh1", "hl.flcih.mseakh1")
  hon2.mseakh2 <- t(as.matrix(c(hmod2.mseakh2$estimate, hmod2.mseakh2$sd, hmod2.mseakh2$hl)))
  colnames(hon2.mseakh2) <- c("tau.flcih.mseakh2", "se.flcih.mseakh2", "hl.flcih.mseakh2")

  #combining bandwdiths
  bw <- t(as.matrix(c(bw.h.mseik1=h.mseik1, bw.h.msecct0=h.msecct0, bw.h.msecct1=h.msecct1, bw.h.msecct2=h.msecct2,
          bw.h.cerccf0=h.cerccf0, bw.h.cerccf1=h.cerccf1, bw.h.cerccf2=h.cerccf2, bw.h.mseak1=h.mseak1, bw.h.mseak2=h.mseak2,
          bw.b.msecct0=b.msecct0, bw.b.msecct1=b.msecct1, bw.b.msecct2=b.msecct2, 
          bw.b.cerccf0=b.cerccf0, bw.b.cerccf1=b.cerccf1, bw.b.cerccf2=b.cerccf2,
          bw.h.mseakh1=h.mseakh1, bw.h.mseakh2=h.mseakh2, mhat=mhat)))
  
  return(estimates=cbind(rob.mseik1, rob.msecct0, rob.msecct1, rob.msecct2, rob.cerccf0, rob.cerccf1, rob.cerccf2, 
                         rob.mseak1, rob.mseak2, rob.mseakh1, rob.mseakh2,
                         hon.mseik1, hon.msecct1, hon.msecct2, hon.cerccf1, hon.cerccf2, 
                         hon.mseak1, hon.mseak2, hon.mseakh1, hon.mseakh2, 
                         hon2.mseik1, hon2.msecct1, hon2.msecct2, hon2.cerccf1, hon2.cerccf2, 
                         hon2.mseak1, hon2.mseak2, hon2.mseakh1, hon2.mseakh2, bw))
  

}

# setup for beta 1 (ind) ====
S=50000 # number of iterations of the simulation
c=0 #cutoff for purposes of eda below

# Data generation parameters
alpha <- 14
beta <- 7

n <- 140 #m=10
set.seed(2227)
for (i in 1:S){
  
  #Store the state of the RNG for debugging purposes
  state <- t(.Random.seed)
  if(i==1) write.table(state, file="compstudysim0519beta1m10rstates.csv", sep=",", row.names=FALSE)
  else write.table(state, file="compstudysim0519beta1m10rstates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  
  
  #Data generation 
  

  
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
  
  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot, rot.below=rot.below, rot.above=rot.above)))
  
  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn"), function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)
  
  #Get estimates
  tauhat.design1 <- tau.hat.one.y.fn(y.vec=ymat[,1], x.vec=x, c=0, dbound=16.2)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")
  
  tauhat.design2 <- tau.hat.one.y.fn(y.vec=ymat[,2], x.vec=x, c=0, dbound=233.26)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")
  
  tauhat.design3 <- tau.hat.one.y.fn(y.vec=ymat[,3], x.vec=x, c=0, dbound=6)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")
  
  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3, eda)

  
  #Output data
  if (i==1) write.table(output, file="compstudysim0519beta1m10estimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="compstudysim0519beta1m10estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)

  if(i==S) {
    state <- t(.Random.seed)
    write.table(state, file="compstudysim0519beta1m10rstates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  }
}


n <- 337 #m=20
set.seed(2795)
for (i in 1:S){
  
  #Store the state of the RNG for debugging purposes
  state <- t(.Random.seed)
  if(i==1) write.table(state, file="compstudysim0519beta1m20rstates.csv", sep=",", row.names=FALSE)
  else write.table(state, file="compstudysim0519beta1m20rstates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  
  
  #Data generation 
  

  
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
  
  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot, rot.below=rot.below, rot.above=rot.above)))
  
  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn"), function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)
  
  #Get estimates
  tauhat.design1 <- tau.hat.one.y.fn(y.vec=ymat[,1], x.vec=x, c=0, dbound=16.2)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")
  
  tauhat.design2 <- tau.hat.one.y.fn(y.vec=ymat[,2], x.vec=x, c=0, dbound=233.26)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")
  
  tauhat.design3 <- tau.hat.one.y.fn(y.vec=ymat[,3], x.vec=x, c=0, dbound=6)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")
  
  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3, eda)
  
  
  #Output data
  if (i==1) write.table(output, file="compstudysim0519beta1m20estimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="compstudysim0519beta1m20estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  
  if(i==S) {
    state <- t(.Random.seed)
    write.table(state, file="compstudysim0519beta1m20rstates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  }
}


n <- 681 #m=35
set.seed(8973)
for (i in 1:S){
  
  #Store the state of the RNG for debugging purposes
  state <- t(.Random.seed)
  if(i==1) write.table(state, file="compstudysim0519beta1m35rstates.csv", sep=",", row.names=FALSE)
  else write.table(state, file="compstudysim0519beta1m35rstates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  
  
  #Data generation 
  

  
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
  
  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot, rot.below=rot.below, rot.above=rot.above)))
  
  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn"), function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)
  
  #Get estimates
  tauhat.design1 <- tau.hat.one.y.fn(y.vec=ymat[,1], x.vec=x, c=0, dbound=16.2)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")
  
  tauhat.design2 <- tau.hat.one.y.fn(y.vec=ymat[,2], x.vec=x, c=0, dbound=233.26)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")
  
  tauhat.design3 <- tau.hat.one.y.fn(y.vec=ymat[,3], x.vec=x, c=0, dbound=6)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")
  
  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3, eda)
  
  
  #Output data
  if (i==1) write.table(output, file="compstudysim0519beta1m35estimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="compstudysim0519beta1m35estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  
  if(i==S) {
    state <- t(.Random.seed)
    write.table(state, file="compstudysim0519beta1m35rstates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  }
}


n <- 1066 #m=50
set.seed(1902)
for (i in 1:S){
  
  #Store the state of the RNG for debugging purposes
  state <- t(.Random.seed)
  if(i==1) write.table(state, file="compstudysim0519beta1m50rstates.csv", sep=",", row.names=FALSE)
  else write.table(state, file="compstudysim0519beta1m50rstates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  
  
  #Data generation 
  

  
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
  
  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot, rot.below=rot.below, rot.above=rot.above)))
  
  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn"), function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)
  
  #Get estimates
  tauhat.design1 <- tau.hat.one.y.fn(y.vec=ymat[,1], x.vec=x, c=0, dbound=16.2)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")
  
  tauhat.design2 <- tau.hat.one.y.fn(y.vec=ymat[,2], x.vec=x, c=0, dbound=233.26)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")
  
  tauhat.design3 <- tau.hat.one.y.fn(y.vec=ymat[,3], x.vec=x, c=0, dbound=6)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")
  
  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3, eda)
  
  
  #Output data
  if (i==1) write.table(output, file="compstudysim0519beta1m50estimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="compstudysim0519beta1m50estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  
  if(i==S) {
    state <- t(.Random.seed)
    write.table(state, file="compstudysim0519beta1m50rstates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  }
}


# setup for beta 2 (ik) ====
S=50000 # number of iterations of the simulation
c=0 #cutoff for purposes of eda below

# Data generation parameters
alpha <- 2
beta <- 4

n <- 56 #m=10
set.seed(5217)
for (i in 1:S){
  
  #Store the state of the RNG for debugging purposes
  state <- t(.Random.seed)
  if(i==1) write.table(state, file="compstudysim0519beta2m10rstates.csv", sep=",", row.names=FALSE)
  else write.table(state, file="compstudysim0519beta2m10rstates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  
  
  #Data generation 
  

  
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
  
  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot, rot.below=rot.below, rot.above=rot.above)))
  
  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn"), function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)
  
  #Get estimates
  tauhat.design1 <- tau.hat.one.y.fn(y.vec=ymat[,1], x.vec=x, c=0, dbound=16.2)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")
  
  tauhat.design2 <- tau.hat.one.y.fn(y.vec=ymat[,2], x.vec=x, c=0, dbound=233.26)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")
  
  tauhat.design3 <- tau.hat.one.y.fn(y.vec=ymat[,3], x.vec=x, c=0, dbound=6)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")
  
  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3, eda)
  
  
  #Output data
  if (i==1) write.table(output, file="compstudysim0519beta2m10estimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="compstudysim0519beta2m10estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  
  if(i==S) {
    state <- t(.Random.seed)
    write.table(state, file="compstudysim0519beta2m10rstates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  }
}


n <- 133 #m=20
set.seed(3266)
for (i in 1:S){
  
  #Store the state of the RNG for debugging purposes
  state <- t(.Random.seed)
  if(i==1) write.table(state, file="compstudysim0519beta2m20rstates.csv", sep=",", row.names=FALSE)
  else write.table(state, file="compstudysim0519beta2m20rstates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  
  
  #Data generation 
  

  
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
  
  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot, rot.below=rot.below, rot.above=rot.above)))
  
  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn"), function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)
  
  #Get estimates
  tauhat.design1 <- tau.hat.one.y.fn(y.vec=ymat[,1], x.vec=x, c=0, dbound=16.2)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")
  
  tauhat.design2 <- tau.hat.one.y.fn(y.vec=ymat[,2], x.vec=x, c=0, dbound=233.26)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")
  
  tauhat.design3 <- tau.hat.one.y.fn(y.vec=ymat[,3], x.vec=x, c=0, dbound=6)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")
  
  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3, eda)
  
  
  #Output data
  if (i==1) write.table(output, file="compstudysim0519beta2m20estimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="compstudysim0519beta2m20estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  
  if(i==S) {
    state <- t(.Random.seed)
    write.table(state, file="compstudysim0519beta2m20rstates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  }
}


n <- 267 #m=35
set.seed(7000)
for (i in 1:S){
  
  #Store the state of the RNG for debugging purposes
  state <- t(.Random.seed)
  if(i==1) write.table(state, file="compstudysim0519beta2m35rstates.csv", sep=",", row.names=FALSE)
  else write.table(state, file="compstudysim0519beta2m35rstates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  
  
  #Data generation 
  

  
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
  
  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot, rot.below=rot.below, rot.above=rot.above)))
  
  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn"), function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)
  
  #Get estimates
  tauhat.design1 <- tau.hat.one.y.fn(y.vec=ymat[,1], x.vec=x, c=0, dbound=16.2)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")
  
  tauhat.design2 <- tau.hat.one.y.fn(y.vec=ymat[,2], x.vec=x, c=0, dbound=233.26)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")
  
  tauhat.design3 <- tau.hat.one.y.fn(y.vec=ymat[,3], x.vec=x, c=0, dbound=6)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")
  
  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3, eda)
  
  
  #Output data
  if (i==1) write.table(output, file="compstudysim0519beta2m35estimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="compstudysim0519beta2m35estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  
  if(i==S) {
    state <- t(.Random.seed)
    write.table(state, file="compstudysim0519beta2m35rstates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  }
}


n <- 417 #m=50
set.seed(9628)
for (i in 1:S){
  
  #Store the state of the RNG for debugging purposes
  state <- t(.Random.seed)
  if(i==1) write.table(state, file="compstudysim0519beta2m50rstates.csv", sep=",", row.names=FALSE)
  else write.table(state, file="compstudysim0519beta2m50rstates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  
  
  #Data generation 
  

  
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
  
  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot, rot.below=rot.below, rot.above=rot.above)))
  
  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn"), function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)
  
  #Get estimates
  tauhat.design1 <- tau.hat.one.y.fn(y.vec=ymat[,1], x.vec=x, c=0, dbound=16.2)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")
  
  tauhat.design2 <- tau.hat.one.y.fn(y.vec=ymat[,2], x.vec=x, c=0, dbound=233.26)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")
  
  tauhat.design3 <- tau.hat.one.y.fn(y.vec=ymat[,3], x.vec=x, c=0, dbound=6)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")
  
  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3, eda)
  
  
  #Output data
  if (i==1) write.table(output, file="compstudysim0519beta2m50estimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="compstudysim0519beta2m50estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  
  if(i==S) {
    state <- t(.Random.seed)
    write.table(state, file="compstudysim0519beta2m50rstates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  }
}


# setup for beta 3 (ak) ====
S=50000 # number of iterations of the simulation
c=0 #cutoff for purposes of eda below

# Data generation parameters
alpha <- 1
beta <- 1

n <- 40 #m=10
set.seed(9246)
for (i in 1:S){
  
  #Store the state of the RNG for debugging purposes
  state <- t(.Random.seed)
  if(i==1) write.table(state, file="compstudysim0519beta3m10rstates.csv", sep=",", row.names=FALSE)
  else write.table(state, file="compstudysim0519beta3m10rstates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  
  
  #Data generation 
  

  
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
  
  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot, rot.below=rot.below, rot.above=rot.above)))
  
  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn"), function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)
  
  #Get estimates
  tauhat.design1 <- tau.hat.one.y.fn(y.vec=ymat[,1], x.vec=x, c=0, dbound=16.2)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")
  
  tauhat.design2 <- tau.hat.one.y.fn(y.vec=ymat[,2], x.vec=x, c=0, dbound=233.26)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")
  
  tauhat.design3 <- tau.hat.one.y.fn(y.vec=ymat[,3], x.vec=x, c=0, dbound=6)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")
  
  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3, eda)
  
  
  #Output data
  if (i==1) write.table(output, file="compstudysim0519beta3m10estimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="compstudysim0519beta3m10estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  
  if(i==S) {
    state <- t(.Random.seed)
    write.table(state, file="compstudysim0519beta3m10rstates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  }
}


n <- 96 #m=20
set.seed(6485)
for (i in 1:S){
  
  #Store the state of the RNG for debugging purposes
  state <- t(.Random.seed)
  if(i==1) write.table(state, file="compstudysim0519beta3m20rstates.csv", sep=",", row.names=FALSE)
  else write.table(state, file="compstudysim0519beta3m20rstates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  
  
  #Data generation 
  

  
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
  
  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot, rot.below=rot.below, rot.above=rot.above)))
  
  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn"), function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)
  
  #Get estimates
  tauhat.design1 <- tau.hat.one.y.fn(y.vec=ymat[,1], x.vec=x, c=0, dbound=16.2)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")
  
  tauhat.design2 <- tau.hat.one.y.fn(y.vec=ymat[,2], x.vec=x, c=0, dbound=233.26)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")
  
  tauhat.design3 <- tau.hat.one.y.fn(y.vec=ymat[,3], x.vec=x, c=0, dbound=6)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")
  
  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3, eda)
  
  
  #Output data
  if (i==1) write.table(output, file="compstudysim0519beta3m20estimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="compstudysim0519beta3m20estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  
  if(i==S) {
    state <- t(.Random.seed)
    write.table(state, file="compstudysim0519beta3m20rstates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  }
}


n <- 193 #m=35
set.seed(6976)
for (i in 1:S){
  
  #Store the state of the RNG for debugging purposes
  state <- t(.Random.seed)
  if(i==1) write.table(state, file="compstudysim0519beta3m35rstates.csv", sep=",", row.names=FALSE)
  else write.table(state, file="compstudysim0519beta3m35rstates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  
  
  #Data generation 
  

  
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
  
  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot, rot.below=rot.below, rot.above=rot.above)))
  
  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn"), function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)
  
  #Get estimates
  tauhat.design1 <- tau.hat.one.y.fn(y.vec=ymat[,1], x.vec=x, c=0, dbound=16.2)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")
  
  tauhat.design2 <- tau.hat.one.y.fn(y.vec=ymat[,2], x.vec=x, c=0, dbound=233.26)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")
  
  tauhat.design3 <- tau.hat.one.y.fn(y.vec=ymat[,3], x.vec=x, c=0, dbound=6)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")
  
  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3, eda)
  
  
  #Output data
  if (i==1) write.table(output, file="compstudysim0519beta3m35estimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="compstudysim0519beta3m35estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  
  if(i==S) {
    state <- t(.Random.seed)
    write.table(state, file="compstudysim0519beta3m35rstates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  }
}


n <- 301 #m=50
set.seed(0151)
for (i in 1:S){
  
  #Store the state of the RNG for debugging purposes
  state <- t(.Random.seed)
  if(i==1) write.table(state, file="compstudysim0519beta3m50rstates.csv", sep=",", row.names=FALSE)
  else write.table(state, file="compstudysim0519beta3m50rstates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  
  
  #Data generation 
  

  
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
  
  eda <- t(as.matrix(c(count.below=count.below, count.above=count.above, h.rot=h.rot, rot.below=rot.below, rot.above=rot.above)))
  
  #Generate y matrix, one column for each function
  ymat <- sapply(c("design1.fn", "design2.fn", "design3.fn"), function(myname, myx){return(do.call(myname, list(x=myx)))}, myx=x)
  
  #Get estimates
  tauhat.design1 <- tau.hat.one.y.fn(y.vec=ymat[,1], x.vec=x, c=0, dbound=16.2)
  colnames(tauhat.design1) <- paste(colnames(tauhat.design1), "design1", sep=".")
  
  tauhat.design2 <- tau.hat.one.y.fn(y.vec=ymat[,2], x.vec=x, c=0, dbound=233.26)
  colnames(tauhat.design2) <- paste(colnames(tauhat.design2), "design2", sep=".")
  
  tauhat.design3 <- tau.hat.one.y.fn(y.vec=ymat[,3], x.vec=x, c=0, dbound=6)
  colnames(tauhat.design3) <- paste(colnames(tauhat.design3), "design3", sep=".")
  
  iteration <- i
  output <- cbind(iteration, tauhat.design1, tauhat.design2, tauhat.design3, eda)
  
  
  #Output data
  if (i==1) write.table(output, file="compstudysim0519beta3m50estimates.csv", sep=",", row.names=FALSE)
  else write.table(output, file="compstudysim0519beta3m50estimates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  
  if(i==S) {
    state <- t(.Random.seed)
    write.table(state, file="compstudysim0519beta3m50rstates.csv", sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  }
}