##############
# CRAN packages
##############

#CRAN_package_list <- c("Matrix",'ggplot2','ggpubr','mapdata','maps','data.table','tidyverse','mhsmm',
#                       'diagram','fmsb','wesanderson','colorspace','latex2exp','png','raster')

CRAN_package_list <- c('mhsmm','gtools','tictoc','wesanderson','colorspace')

not_in_list <- CRAN_package_list[!(CRAN_package_list %in% installed.packages())]
lapply(not_in_list,install.packages,dependencies=TRUE)

lapply(CRAN_package_list, require, character.only = TRUE)
       
# clean
rm(CRAN_package_list,not_in_list)
