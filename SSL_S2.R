
source("Models.R")
source("Simulation_new.R")
#####################
Generate=Generate.S2
Case="S2"
set.seed(1)
N.all=c(1000,2000)
p.all=c(10,50)
c=1.6
T=2


Simulation()

