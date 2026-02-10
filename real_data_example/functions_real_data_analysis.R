# This file contains all packages and functions needed to generate the results
#  in the subsection "Real data example".

if (!require("docstring")) install.packages("docstring")
if (!require("ExtremalDep")) install.packages("ExtremalDep")
if (!require("plyr")) install.packages("plyr")
if (!require("mvtnorm")) install.packages("mvtnorm")
if (!require("pracma")) install.packages("pracma")

library(docstring)
library(ExtremalDep)
library(plyr)
library(mvtnorm)
library(pracma)

source("Q_hat.R")
source("theta_hill.R")
source("cens_bbeed_mar.R")
source("bbeed_thresh_mar.R")

assignInNamespace("bbeed.thresh.mar", bbeed.thresh.mar, pos="package:ExtremalDep")
assignInNamespace("cens.bbeed.mar", cens.bbeed.mar, pos="package:ExtremalDep")
