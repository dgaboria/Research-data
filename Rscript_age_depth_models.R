# UTF-8
# Project : PhD in environmental sciences 2017-2020
# Author: Dorian Gaboriau - gaboriau.dorian@gmail.com - dorian.gaboriau@uqat.ca
# Thesis: Prediction of wildfire regimes up to 2100 in Northwest territories (NWT) and evaluation of the impacts for the Tlicho First Nation
# Chapter 2 : Reconstitution and caracterization of the past fire regime, vegetation and temperatures (Holocene period) in the central NWT, Canada 
# Gaboriau Dorian
# Last update : July 20th 2020
# R code for construction of Winbacon age depth models

rm(list = ls())  # Deleting variables from the environment R

##############################################################################################################################
# WinBacon

###Packages and files

#Load the package WinBacon
wd = "C://Users//client//Documents//R//win-library//3.5//WinBacon_2.2"
setwd(wd)
source("bacon.R")

#Bacon (Blaauw and Christen, 2011) is an approach to age-depth modelling that uses Bayesian statistics to reconstruct
#Bayesian accumulation histories for deposits, through combining radiocarbon and other dates with prior information.

#In many cases Bacon isn't bothered by outlying dates, since
#the dates are modelled using a student-t distribution with wide tails (Christen and Pérez, 2009). 

#Any age-depth model produces estimates of accumulation rates, either implicitly or explicitly. For example, the most
#popular method of connecting the mid-points of dated levels using linear sections (Blaauw, 2010), assumes that a
#deposit accumulated linearly between each dated level, and that changes in accumulation rate took place abruptly and
#exactly at the dated depths.
#Bacon divides a core into many thin vertical sections (by default of res=5 cm thickness), and through millions of
#Markov Chain Monte Carlo (MCMC) iterations estimates the accumulation rate (in years/cm; so more correctly,
#sedimentation times) for each of these sections. Combined with an estimated starting date for the first section, these
#accumulation rates then form the age-depth model. The accumulation rates are constrained by prior information as
#described below. 

#The section thickness (default res=5) will dictate to some degree the flexibility of the age-depth model. 

#At certain core depths, hiatuses can be inferred. 

#Bacon  is often used to age-model 14C-dated sequences

#Radiocarbon dates should be calibrated using either IntCal13 (for terrestrial northern hemisphere material; Reimer et al., 2013)
#(cc = 1)

#Bacon's calendar scale is cal BP (calendar years before AD 1950)

#Four columns (labID, age, error, depth), separated by commas, .csv file

##############################################################################################################################
# LAKE EMILE

dir("./Cores")
#head(read.csv("./Cores/MSB2K/MSB2K.csv"))
read.csv("./Cores/Emile/Lake_Emile.csv")

b1 = Bacon("Emile", d.min = 0, runname = "", d.max = 456, d.by = 0.5, cc = 1, prob = 0.95, rev.yr = T, postbomb = 1, rev.d=F, unit = "cm", rotate.axes = T, normal = T)

#Read the table with ages
Emile_bacon_ages = read.table("C:/Users/client/Documents/R/win-library/3.5/winBacon_2.2/Cores/Emile/Emile_92_ages.txt", header = T, sep = "")
Emile_bacon_ages

##############################################################################################################################
# LAKE IZAAC

read.csv("./Cores/Izaac/Lake_Izaac")

#Read the table with ages
b1 = Bacon("Izaac", d.min = 0, runname = "", d.max = 359, d.by = 0.5, cc = 1, prob = 0.95, rev.yr = T, postbomb = 1, rev.d=F, unit = "cm", rotate.axes = T, normal = T)
Izaac_bacon_ages = read.table("C:/Users/client/Documents/R/win-library/3.5/winBacon_2.2/Cores/Izaac/Izaac_72_ages.txt", header = T, sep = "")
Izaac_bacon_ages

##############################################################################################################################
# LAKE PARADIS

read.csv("./Cores/Paradis/Lake_Paradis")

#Read the table with ages
b1 = Bacon("Paradis", d.min = 0, runname = "", d.max = 450, d.by = 0.5, cc = 1, prob = 0.95, rev.yr = T, postbomb = 1, rev.d=F, unit = "cm", rotate.axes = T, normal = T)
Paradis_bacon_ages = read.table("C:/Users/client/Documents/R/win-library/3.5/winBacon_2.2/Cores/Paradis/Paradis_91_ages.txt", header = T, sep = "")
Paradis_bacon_ages

##############################################################################################################################
# LAKE SAXON

read.csv("./Cores/Saxon/Lake_Saxon")

#Read the table with ages
b1 = Bacon("Saxon", d.min = 0, runname = "", d.max = 350, d.by = 0.5, cc = 1, prob = 0.95, rev.yr = T, postbomb = 1, rev.d=F, unit = "cm", rotate.axes = T, normal = T)
Saxon_bacon_ages = read.table("C:/Users/client/Documents/R/win-library/3.5/winBacon_2.2/Cores/Saxon/Saxon_71_ages.txt", header = T, sep = "")
Saxon_bacon_ages


