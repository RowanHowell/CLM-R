library(tidyverse)
library(igraph)
library(readxl)
library(stringr)
library(CellNOptR)
library(data.table)
source("TrainingFuncs.R")

CNO = CNOlist("FEARCNO.csv")
PKN = readSIF("FEARPKNedit.txt")



processPKN = preprocessing(model = PKN, data = CNO, expansion = T, compression = F, maxInputsPerGate = 4, verbose =FALSE)

# the "priorBitString" is ued to force inclusion of specific edges
force = as.numeric(processPKN$reacID%in%grep("_OE", PKN$reacID, value = TRUE)) # OE nodes should always activate the protein, so we force these edges to be maintained
force[which(processPKN$reacID=="!Cdc5_OE=Net1")]=NA # don't include this edge in forced edges
force[force==0] = NA

# Initial run required to set up data frame to store results
fullres = gaBinaryT1(model= processPKN, CNO = CNO, popSize = 200, maxTime = Inf, maxGens = Inf, relTol = 0.01,priorBitString = force, verbose =TRUE)
allNodes_T1 = simulateTN(CNO,processPKN,bStrings = list(fullres$bString))
indexList = indexFinder(CNO, processPKN)
measuredNodes_T1 = allNodes_T1[,indexList$signals]
score = 100*length(which(getSignals(CNO)$`1`== measuredNodes_T1))/nrow(getSignals(CNO)$`1`)

results = matrix(c(score, fullres$bString), nrow = 1, ncol = length(fullres$bString)+1)
colnames(results) = c("score", processPKN$reacID)
results = as.data.frame(results)

for(j in 2:3){
  tfullres = gaBinaryT1(model= processPKN, CNO = CNO, popSize = 200, maxTime = Inf, maxGens = Inf, relTol = 0.01,priorBitString = force, verbose =FALSE)
  allNodes_T1 = simulateTN(CNO,processPKN,bStrings = list(tfullres$bString))
  measuredNodes_T1 = allNodes_T1[,indexList$signals]
  tscore = 100*length(which(getSignals(CNO)$`1`== measuredNodes_T1))/nrow(getSignals(CNO)$`1`)
  results = rbind(results, c(tscore, tfullres$bString))
  print(j)
}

write_csv(results,"output.csv")


colnames(results) = c("score", processPKN$reacID)

allModels = results[, 2:ncol(results)]

scores = results$score
hist(scores)

optModels = allModels[which(scores==max(scores)), ]
uniqModels = unique(optModels)

hist(colSums(uniqModels))

minModels = which(rowSums(uniqModels)==min(rowSums(uniqModels)))

convertToBoolNet(processPKN$namesSpecies, processPKN$reacID[which(as.logical(uniqModels[minModels[1],]))],"FEARnet.txt")

write_csv(tibble(edges = processPKN$reacID[which(as.logical(uniqModels[minModels[1],]))]),"FEARedges.csv", col_names=FALSE)

png("trainedFEARnet.png")
plotModel(model = processPKN, CNOlist = CNO, bString = as.numeric(uniqModels[minModels[1],]))
dev.off()
