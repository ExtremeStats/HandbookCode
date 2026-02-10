# This file contains all packages and functions needed to generate the results
#  in the section "Modeling Extreme Value Dependence".

if (!require("docstring")) install.packages("docstring")
if (!require("ExtremalDep")) install.packages("ExtremalDep")
#if (!require("plyr")) install.packages("plyr")
if (!require("mvtnorm")) install.packages("mvtnorm")
#if (!require("pracma")) install.packages("pracma")
if (!require("ggplot2")) install.packages("ggplot2")

library(docstring)
library(ExtremalDep)
#library(plyr)
library(mvtnorm)
#library(pracma)
library(ggplot2)

source("Q_hat.R")
source("truncated_bivariate_t_distribution.R")