library(tidyverse)
library(BoolNet)
library(stringr)
library(reshape2)
library(ggplot2)
source("CreateCLM.R")
source("CLMfunctions.R") 
source("simulateMutantSync.R")
source("simulateMutantAsync.R")

# create network file
createNetwork("Activity3.5.txt", "Localization3.5.txt", "FEARnet3.3.txt", "LocSpecificActivity3.5.txt", "model3.5", directory = "networkData")

# Make models
model3.5 = loadNetwork("model3.5OE.txt") # model with OE
model3.5noOE = loadNetwork("model3.5.txt") # model with no OE
model3.5forceLoc = loadNetwork("model3.5SPB.txt") # model with forced localization

# Make initial conditions
ICs3.5 = makeICs("MetaAIC3.5.csv", "MetaLIC3.5.csv", model3.5, directory = "networkData")
ICs3.5noOE = makeICs("MetaAIC3.5.csv", "MetaLIC3.5.csv", model3.5noOE, directory = "networkData")
ICs3.5FL = makeICs("MetaAIC3.5.csv", "MetaLIC3.5.csv", model3.5forceLoc, directory = "networkData")

# Synchronous steady states
WTattractors = attractorsMultiMutant(prots = "NA", mutations = "NA", model = model3.5, mutname = "WT", unregfile = "model3.5unreg.csv", ICs3.5)

bub2attractors = attractorsMultiMutant(prots = "Bub2", mutations = "delete", model = model3.5, mutname = "bub2delete", unregfile = "model3.5unreg.csv", ICs3.5)

# Asynchronous trajectories
WTEA = asyncStats(prot = "NA", mutation = "NA", model = model3.5, unregfile = "model3.5unreg.csv", ICs = ICs3.5, CCstage = "EA", timeMax =10000, N = 100)

bub2EA = asyncStats(prot = "Bub2", mutation = "delete", model = model3.5, unregfile = "model3.5unreg.csv", ICs = ICs3.5, CCstage = "EA", timeMax =10000, N = 100)

sum(as.numeric(WTEA$exits=="ME")) # find number of cells that exit mitosis
sum(as.numeric(bub2EA$exits=="ME"))


