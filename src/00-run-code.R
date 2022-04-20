# run all code
# author: sebastian daza
# R < src/00-run-code.R > output/run-code.log  --no-save  &

# empirical analysis
# R < src/01-hrs-gxe-analysis.R > output/hrs.log  --no-save  &
source("src/01-hrs-gxe-analysis.R")
# R < src/02-elsa-gxe-analysis.R > output/elsa.log  --no-save  &
source("src/02-elsa-gxe-analysis.R")

# simulation
# R < src/03-simulations.R > output/simulation.log  --no-save  &
source("src/03-simulations.R")

# create session info file
sink("output/R-session-empirical-analysis.txt")
sessionInfo()
sink()
