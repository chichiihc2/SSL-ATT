
source("Models.R")
source("Simulation_new.R")
#####################
Generate=Generate.S3
Case="S3"
set.seed(1)
N.all=c(1000,2000)
p.all=c(10,50)
c=0.2
T=2
Simulation()

