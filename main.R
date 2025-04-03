# Problema de permissoes, é preciso instalar packages manhosos em .libPath()[2]

install.packages('TMB', type = 'source', lib = .libPaths()[2]) 

library('RTMB', lib.loc = .libPaths()[2])
library('tidyverse')
library(CatDyn)
 