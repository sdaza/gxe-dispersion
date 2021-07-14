library(tidyverse)
library(brms)

set.seed(201910052)
n <- 5000
people <- tibble(
  id = 1:n,
  days_WFH = rep(0:5, length.out = n),
  parent = id > n/2,
  sleep = rnorm(n, if_else(parent, 6, 8), if_else(parent, 5, 1)),
  RandomVariation = rnorm(n)
)

table(round(people$sleep))
xtabs(~ parent + days_WFH, people)
table(people$days_WFH)
people <- people %>% mutate(
  Productivity = 20 + sleep + if_else(parent, 0.31, 0.1) * days_WFH + 1.75 * RandomVariation
)
ggplot(data = people) + geom_histogram(aes(sleep, fill = parent), position = position_identity(), alpha = 0.4)
ggplot(data = people) + geom_histogram(aes(Productivity, fill = parent), position = position_identity(), alpha = 0.4)

model <- brm(Productivity ~ parent * days_WFH, data = people, cores = 1)

model2 <- brm(bf(Productivity ~ parent * days_WFH, 
                 sigma ~ parent), data = people, cores = 1)

summary(model) 
summary(model2)

people$zproductivity = scale(people$Productivity)
people$zdays_wfh = scale(people$days_WFH)

model3 <- brm(zproductivity ~ parent * days_WFH, data = people, cores = 1)
summary(model3)

people[,c("estimate__", "std.error__", "lower__", "upper__")] <- fitted(model3, re_formula = NA)

(cors <- people %>% group_by(parent) %>%
  summarise(cor = cor(days_WFH, zproductivity),
            slope = broom::tidy(lm(zproductivity ~ days_WFH))[2, "estimate"][[1]],
            sd_resid = sd(resid(lm(zproductivity ~ days_WFH))),
            sd_prod = sd(zproductivity),
            mean_prod = mean(zproductivity),
            max_js = max(days_WFH),
            mean_js = mean(days_WFH),
            mean_est = mean(estimate__),
            sd_js = sd(days_WFH))
)
