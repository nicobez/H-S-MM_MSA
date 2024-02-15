#############
#
# Make file 
# author : nicolas.bez@ird.fr
#
# # # # # # # # # # # # # # # # # # # # # # # # #
#
# AR1 : m*(1-mu)=eta where mu=parameters[,,,3,] is the coef of auto-correlation
# AR0 : mu, eta parameters are not used
# HMM : sojourn times are geometrical ==> rgeom(1,prob=1/m.soj.1)
# HSMM : sojourn time are shifted-Poisson (time ~ shift + Poisson(m-shift) ; variance = m-shift) 
#                      or shifted-NegBinom (time ~ shift + nbinom(size,prob) where size*(1-prob)/prob = m-shift ; variance = size*(1-prob)/prob?)
#
# state 1 = NON Fishing ; state 2 = Fishing
#
# m.V.1 : means of state dependent variables in state 1
# s2.V.1 : variances of state dependent variables in state 1 ; independence ==> vector and not a matrix 
# mu.V.1 : coefficients for AR1
# m.soj : means of sojourn time in the different states ==> vector with length the number of states
#         units are in number of time steps
#
# The following codes work for either of the two sojour time pmfs 
#    parametersP <- parameters[,,c(1,2,3),,]
#    parametersNB <- parameters[,,c(1,2,4),,]
#
# # # # # # # # # # # # # # # # # # # # # # # #
# 
# **********************************************************************
#   Relationship between notation in mathematical equations and in R code
# **********************************************************************
#   Equation          | R
#   =================================
#   tau (coef AR1)    | mu (coef AR1)
#   nu                | eta
#   mu                | m
#   ---------------------------------
#   mu=nu/(1-tau)     | m=eta/(1-mu)
#
# **********************************************************************
#   Degradation
# **********************************************************************
#   d1:d4 in R objects is reverse wrt to resolution (d1 High Def ==> d4 low def)
#
#############

rm(list=ls())

# Upload required packages
source("Script/libraries.R")

# Load specific functions for AR1 and for specific tasks
source("Script/addsOnMhsmm.R")

# Loading data, simulations and working objects
source("Script/loadingModelParameters.R")

# Setting the graphical environment
source("Script/graphicalEnvironment.R")

# MSA + qality control
source("Script/MSA.R")












