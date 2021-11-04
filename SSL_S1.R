
source("Models.R")
source("Simulation_new2.R")
#####################
Generate=Generate.S1
Case="S1"
ratio=0.4
set.seed(1)
N.all=c(1000)
p.all=c(4,100)
name=Sys.Date()

c=1.6
T=10

Simulation()

