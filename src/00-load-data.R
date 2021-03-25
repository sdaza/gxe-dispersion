#################################
# explore HRS data
# author: sebastian daza
#################################


library(haven)
library(data.table)


# read data
dat = data.table(read_stata("data/hrs-obesity-pgs.dta"))
names(dat)
str(dat)
