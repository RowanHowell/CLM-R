# function to find attractors from the synchronous update scheme

library(tidyverse)
library(BoolNet)
library(stringr)
library(reshape2)
library(ggplot2)
source("CLMfunctions.R") # note this dependency

attractorsMultiMutant = function(prots = "NA", mutations = "NA", model, mutname, unregfile, ICs){
  # Find attractors of the synchronous update scheme
  ## INPUTS
  # prots - list of proteins eg c("Bub2", "Cdc5")
  # mutations - list of mutations eg c("delete", "KD")
  # model - BoolNet model, created by CreateNetworks
  # mutname - name for SS file
  # unregfile - list of unregulated nodes (this is needed for efficient simulation), this is created by CreateNetwork function, called modelunreg.csv
  # ICs - list of ICs, built by makeICs function
  
  ## OUTPUTS
  # List of matrices, entries $X where X is the size of the attractor (eg X = 1 are steady states, X = 2 are cyclic attractors with 2 states)
  # Each attractor is given a number, if cyclic attractors, each state is given a sub number (eg if state 9 has 2 states these are 9.1 and 9.2)
  # Each column of the matrix is a state, with the row showing the activity of a protein in that state (1 is ON, 0 is OFF)
  # This script also creates a folder called SSPlots and puts a diagram showing the localization of Cdc14 in each of the steady states. In these plot
  # if metaphase is ON, the cell is in metaphase, otherwise it is in anaphase. If Spindle alignment is ON then the spindle has aligned with the mother bud axis.
  
  unreg = read_csv(unregfile)$node
  gON = c()
  gOFF = c()
  for(n in unreg){
    if(ICs[n]==1){
      gON = c(gON, n)
    } else if(ICs[n]==0){
      gOFF = c(gOFF, n)
    }
    else{  print(paste0("Error: no initial condition for ",n))}
  }
  gOFF = c(gOFF,model$genes[grep("OE",model$genes)]) # all OE shoudl be switched off
  locs =  c("dSPB", "mSPB", "Nucleus", "Cytoplasm","Bud")
  for(j in 1:length(prots)){
    protOG = prots[j]
    if(prots[j]%in% c("Bfa1", "Bub2", "Tem1", "Cdc15", "Dbf2", "Cdc14", "CDK")){ # multilevel
      protC = paste0(prots[j], c("low","high"))
    } else {protC = prots[j]}
    if(mutations[j] == "OE"){
      allOccurences = grep(paste0("^",protOG),model$genes, value = TRUE) # just look for nodes starting with this sequence (eg for Cdc5/PP2ACdc55 issue)
      gON = c(gON, paste0(protOG, "OE"), allOccurences)
      gOFF = setdiff(gOFF, c(paste0(protOG, "OE"), allOccurences))
    } else if(mutations[j] %in% c("delete", "deplete")){
      nodes = c(apply(expand.grid(protC, paste0("A_", locs)), 1, paste, collapse=""),apply(expand.grid(protC, paste0("L_", locs)), 1, paste, collapse=""))
      nodes = intersect(nodes, model$genes) # remove SPBs if protein can't localise there
      gOFF = c(gOFF, nodes)
      gON = setdiff(gON, nodes)
    } else if( mutations[j] == "KD"){
      nodes = apply(expand.grid(protC, paste0("A_", locs)), 1, paste, collapse="")
      nodes = intersect(nodes, model$genes)
      gOFF = c(gOFF, nodes)
      gON = setdiff(gON, nodes)
    } else if(mutations[j] == "SPB"){
      nodes = apply(expand.grid(protC, paste0("L_", c("dSPB", "mSPB"))), 1, paste, collapse="")
      gON = c(gON, nodes)
      gOFF = setdiff(gOFF, nodes)
    } # for now only basic mutants and forced SPB localization.
  }
  attract = getAttractors(model, method = "sat.exhaustive", genesON = gON, genesOFF = gOFF)#, startStates = 100)
  pdf(file=NULL)
  plot = plotAttractors(attract, allInOnePlot = TRUE)
  dev.off()
  Cdc14levels = extractCdc14_2(plot$`1`)
  Cdc14levelsM = melt(Cdc14levels)
  
  name = paste0(mutname,".png")
  if(!("SSPlots"%in%ls())){
    dir.create(file.path(getwd(), "SSPlots"), showWarnings = FALSE)
  }
  setwd("SSPlots")
  png(name, height = 400, width = 1000)
  print(ggplot(Cdc14levelsM, aes(x = Var2, y = Var1, fill = value)) + geom_tile(colour="grey",size=0.25) + coord_fixed() + theme_bw() +theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),text=element_text(size=16,  family="Helvetica"), legend.title = element_blank(),legend.key.size =unit(1.5, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_fill_manual(values = c("OFF" = "white","ON" = "black","LOW" = "#FDAE6B","HIGH" = "#E6550D")) + xlab("") + ylab("") + scale_y_discrete(labels = c("ME", "Bud", "Cytoplasm", "Nucleus", "Spindle alignment", "Metaphase")))
  dev.off()
  setwd("..")
  return(plot)
}
