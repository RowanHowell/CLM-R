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
createNetwork("Activity5.txt", "Localization5.txt", "FEARnet3.txt", "LocSpecificActivity5.txt", "model5", directory = "networkData")

# Make models
model5 = loadNetwork("model5OE.txt") # model with OE
model5noOE = loadNetwork("model5.txt") # model with no OE
model5forceLoc = loadNetwork("model5SPB.txt") # model with forced localization

# Make initial conditions
ICs5 = makeICs("MetaAIC5.csv", "MetaLIC5.csv", model5, directory = "networkData")
ICs5noOE = makeICs("MetaAIC5.csv", "MetaLIC5.csv", model5noOE, directory = "networkData")
ICs5FL = makeICs("MetaAIC5.csv", "MetaLIC5.csv", model5forceLoc, directory = "networkData")

# Synchronous steady states
WTattractors = attractorsMultiMutant(prots = "NA", mutations = "NA", model = model5, mutname = "WT", unregfile = "model5unreg.csv", ICs5)

bub2attractors = attractorsMultiMutant(prots = "Bub2", mutations = "delete", model = model5, mutname = "bub2delete", unregfile = "model5unreg.csv", ICs5)

# Asynchronous trajectories
WTEA = asyncStats(prot = "NA", mutation = "NA", model = model5, unregfile = "model5unreg.csv", ICs = ICs5, CCstage = "EA", timeMax =10000, N = 100)

bub2EA = asyncStats(prot = "Bub2", mutation = "delete", model = model5, unregfile = "model5unreg.csv", ICs = ICs5, CCstage = "EA", timeMax =10000, N = 100)

sum(as.numeric(WTEA$exits=="ME")) # find number of cells that exit mitosis
sum(as.numeric(bub2EA$exits=="ME"))

