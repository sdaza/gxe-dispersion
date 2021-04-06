
library(data.table)
library(brms)
library(texreg)

n = 2000
G = rnorm(n)
I = rnorm(n)
E = rnorm(n)
DN = rnorm(n)

hg = sqrt(0.4)
hi = sqrt(0.2)
he = sqrt(0.3)
hd = sqrt(0.1)

y = hg * G + he * E + hi * G * E + hd * DN

b0 = 0.2
b1 = 0.2
a = 0.1

var(hg * G) / var(y) 
var(he * E) / var(y)
var(hi * I) / var(y) 
var(hd * DN) / var(y)

sigma = exp(b0 * hd  + b1 * hd * E)
hist(sigma)

yt = 0 + a * E + b0 * hg * G + b0 * he * E + b0 * hi * G * E + 
    b1 * E * hg * G + b1 * E * he * E + b1 * E * hi * G * E + rnorm(n, 0, sigma)


dat = data.table(yt, G, E)
screenreg(lm(yt ~ G + E + E*G, data = dat))

f = bf(yt ~ G + E + E*G, sigma ~ E)
test = brm(f, data = dat)
screenreg(test)

hyp <- "G * sigma_E = G:E * sigma_Intercept"
(hyp <- hypothesis(test, hyp, alpha = 0.05))
plot(hyp)
