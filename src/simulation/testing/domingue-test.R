# scaling_gxe
library(scalingGxE)

##############################################################
##example wherein the scaling model is the true model
set.seed(8675309)
E<-rnorm(5000)
df<-sim.data(E=E,scaling=TRUE)
head(df)

est<-mlest(df,hess=FALSE)
xi.test(est)
##note we fail to reject the null

##############################################################
##example wherein vanilla GxE is the true model

E<-rnorm(8000)
df<-sim.data(E=E,scaling=FALSE)
est<-mlest(df,hess=TRUE)
xi.test(est)

##note we reject the null

se<-sqrt(diag(est$var))
z<-est$est/se
2*pnorm(abs(z[5]),lower.tail=FALSE) #further, note that lambda1 isn't significant given that there is no heteroscedasticity here.

# test using hrs
head(df)
summary(df)

dt = readRDS("data/hrs-obesity-pgs.dta")
