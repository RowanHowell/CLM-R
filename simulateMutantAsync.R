# Functions to simulate mutants using an asynchronous update scheme.

library(tidyverse)
library(BoolNet)
library(stringr)
library(reshape2)
library(ggplot2)
source("CLMfunctions.R") # note these dependencies
source("simulateMutantSync.R")

asyncStats = function(prot = "NA", mutation = "NA", model, unregfile, ICs, CCstage, timeMax, N){ 
  # simulate multiple trajectories for a mutants
  ## INPUTS
  # prots - list of proteins eg c("Bub2", "Cdc5")
  # mutations - list of mutations eg c("delete", "KD")
  # model - BoolNet model, created by CreateNetworks
  # unregfile - list of unregulated nodes (this is needed for efficient simulation), this is created by CreateNetwork function, called modelunreg.csv
  # ICs - list of ICs, built by makeICs function
  # CCstage - one of "M" (metaphase), "EA" (early anaphase, pre-spindle alignment), "LA" (late anaphase, post-spindle alignment, "NP" (non-physiological, aligned spindle in metaphase))
  # timeMax - Time when simulation terminates, usually 10,000
  # N - number of cells to simulate, usually 100
  
  ## OUTPUTs
  # list with 3 entries
  # $exits - list of length N specifying how each simulation ended, "ME" means the mitotic exit node was activated, a number means the trajectory ended
  # in a steady state where mitotic did not occur (it will not leave this state), "time" means the simulation reached timeMax, this either means a cyclic attractor,
  # where mitotic exit is not activated, exists or the simulation does not have enough time to reach an attractor and timeMax should be increased.
  # $times - list of length N specifying the number of steps that were simulated before mitotic exit (=timeMax if it never is)
  # SSs - matrix of synchronous attractors, this is calculated to save time by checking if steady state reached, format is described in simulateMutantSync.R
  
  
  
  unreg = read_csv(unregfile)$node
  gON = c()
  gOFF = c()
  for(n in unreg){
    if(ICs[n]==1){
      gON = c(gON, n)
    } else if(ICs[n]==0){
      gOFF = c(gOFF, n)
    }
    else{ print(paste0("Error: no initial condition for ",n))}
  }
  gOFF = c(gOFF,model$genes[grep("OE",model$genes)]) # all OE shoudl be switched off
  gOFF = c(gOFF,model$genes[grep("FL",model$genes)]) # all FL shoudl be switched off
  locs =  c("dSPB", "mSPB", "Nucleus", "Cytoplasm","Bud", "Nucleolus")
  
  for(j in 1:length(prot)){ # iterate throught list of mutant proteins
    protOG = prot[j]
    if(length(grep(paste0(protOG,"low"), model$genes)>0)){ # multilevel
      prott = paste0(protOG, c("low","high"))
    } else{
      prott = protOG
    }
    mutationt = mutation[j]
    if(!identical(prott,"NA")){
      if(mutationt == "OE"){
        allOccurences = grep(paste0("^",protOG),model$genes, value = TRUE) # just look for nodes starting with this sequence (eg for Cdc5/PP2ACdc55 issue)
        gON = c(gON, paste0(protOG, "OE"), allOccurences)
        gOFF = setdiff(gOFF, c(paste0(protOG, "OE"), allOccurences))
      } else if(mutationt %in% c("delete", "deplete")){
        nodes = c(apply(expand.grid(prott, paste0("A_", locs)), 1, paste, collapse=""),apply(expand.grid(prott, paste0("L_", locs)), 1, paste, collapse=""))
        nodes = intersect(nodes, model$genes) # remove SPBs if protein can't localise there
        gOFF = c(gOFF, nodes)
        gON = setdiff(gON, nodes)
      } else if( mutationt == "KD"){
        nodes = apply(expand.grid(prott, paste0("A_", locs)), 1, paste, collapse="")
        nodes = intersect(nodes, model$genes)
        gOFF = c(gOFF, nodes)
        gON = setdiff(gON, nodes)
      } else if( mutationt == "hyperactive"){
        nodes = apply(expand.grid(prott, paste0("A_", locs)), 1, paste, collapse="")
        nodes = intersect(nodes, model$genes)
        gON = c(gON, nodes)
        gOFF = setdiff(gOFF, nodes)
      }else if(grepl("!", mutationt)){
        modmutationt = sub("!","", mutationt)
        if(modmutationt == "SPB"){
          nodes = apply(expand.grid(prott, paste0("L_", c("dSPB", "mSPB"))), 1, paste, collapse="")
        } else{
          nodes = apply(expand.grid(prott, paste0("L_", modmutationt)), 1, paste, collapse="")
        }
        gOFF = c(gOFF, nodes)
        gON = setdiff(gON, nodes)
      } else if(mutationt == "Cytoplasm+Bud"){
          nodes = apply(expand.grid(prott, paste0("L_", c("Cytoplasm", "Bud"))), 1, paste, collapse="")
          gON = c(gON, nodes)
          gOFF = setdiff(gOFF, nodes)
      } else if(mutationt == "mSPB"){
        nodes = apply(expand.grid(prott, paste0("L_", "mSPB")), 1, paste, collapse="")
        gON = c(gON, c(nodes, paste0(protOG, "FLmSPB")))
        gOFF = setdiff(gOFF, c(nodes, paste0(protOG, "FLmSPB")))
      } else if(mutationt == "dSPB"){
        nodes = apply(expand.grid(prott, paste0("L_", "dSPB")), 1, paste, collapse="")
        gON = c(gON, c(nodes, paste0(protOG, "FLdSPB")))
        gOFF = setdiff(gOFF, c(nodes, paste0(protOG, "FLdSPB")))
      } else{
        if(mutationt == "SPB"){
          nodes = apply(expand.grid(prott, paste0("L_", c("dSPB", "mSPB"))), 1, paste, collapse="")
          nodes = c(nodes, paste0(protOG, "FLmSPB"),paste0(protOG, "FLdSPB"))
        } else{
          nodes = apply(expand.grid(prott, paste0("L_", mutationt)), 1, paste, collapse="")
        }
        gON = c(gON, nodes)
        gOFF = setdiff(gOFF, nodes)
      }
    }
  }
  # find synchronous SSs, the code checks if a SS is reached to save time
  gON = intersect(gON, model$genes) # make sure all gON/gOFF are really nodes
  gOFF = intersect(gOFF, model$genes)
  attract = getAttractors(model, method = "sat.exhaustive", genesON = gON, genesOFF = gOFF)#, startStates = 100)
  
  pdf(file=NULL)
  plot = plotAttractors(attract, allInOnePlot = TRUE)
  dev.off()
  
  # Now fix to specific CCstage
  if(CCstage == "M"){
    gON = c(gON, "SACA_Nucleus")
    gOFF = c(gOFF, "SpindleAlign")
  } else if(CCstage == "EA"){
    gOFF = c(gOFF, "SpindleAlign","SACA_Nucleus")
  } else if(CCstage == "LA"){
    gON = c(gON, "SpindleAlign")
    gOFF = c(gOFF, "SACA_Nucleus")
  }else if(CCstage == "NP"){
    gON = c(gON, "SACA_Nucleus", "SpindleAlign")
  } else {print("hol'up")}
  
  gON = intersect(gON, model$genes) # make sure all gON/gOFF are really nodes
  gOFF = intersect(gOFF, model$genes)
  
  fixmodel = fixGenes(model, c(gON,gOFF), c(rep(1,length(gON)), rep(0, length(gOFF)))) # enforce mutations
  
  IClist = ICs[model$genes]
  IClist[gON] = 1 # make ICs match restrictions
  IClist[gOFF] = 0
  
  exits = c()
  times = c()
  for(i in 1:N){
    flag = 0
    time = 1
    currentState = IClist
    while(flag==0){
      currentState = stateTransition(network = fixmodel, state = currentState, type = "asynchronous") # uniform randomly selected node updated
      time = time+1
      if(1%in% names(plot)){ #if there are steady states
        CC = compareColumns(df = plot$`1`, vec = currentState)
        if(currentState[which(model$genes=="ME")] == 1){ # exit condition 1: ME node activated
          exit = "ME"
          flag = 1
        } else if(CC$result){ # exit condition 2: simulation reaches steady state
          exit = CC$col
          flag = 1
        }
        else if(time == timeMax){ # # exit condition 3: time limit reached
          exit = "time"
          flag = 1
        }
      } else{ # if no steady states exist (only cyclic attractors)
        if(currentState[which(model$genes=="ME")] == 1){
          exit = "ME"
          flag = 1
        } else if(time == timeMax){
          exit = "time"
          flag = 1
        }
      }
    }
    exits = c(exits, exit)
    times = c(times, time)
  }
  return(list(exits=exits, times=times, SSs = plot))
}

asyncStatsFEAR = function(prot = "NA", mutation = "NA", model, unregfile, ICs, CCstage, timeMax, N){ # simulate N trajectories, in given cell cycle stage (CCstage), return data on FEAR release
  
  unreg = read_csv(unregfile)$node
  gON = c()
  gOFF = c()
  for(n in unreg){
    if(ICs[n]==1){
      gON = c(gON, n)
    } else if(ICs[n]==0){
      gOFF = c(gOFF, n)
    }
    else{ print("aargh")}
  }
  gOFF = c(gOFF,model$genes[grep("OE",model$genes)]) # all OE shoudl be switched off
  locs =  c("dSPB", "mSPB", "Nucleus", "Cytoplasm","Bud")
  
  for(j in 1:length(prot)){
    protOG = prot[j]
    if(length(grep(paste0(protOG,"low"), model$genes)>0)){ # multilevel
      prott = paste0(protOG, c("low","high"))
    } else{
      prott = protOG
    }
    mutationt = mutation[j]
    if(!identical(prott,"NA")){
      if(mutationt == "OE"){
        allOccurences = grep(paste0("^",protOG),model$genes, value = TRUE) # just look for nodes starting with this sequence (eg for Cdc5/PP2ACdc55 issue)
        gON = c(gON, paste0(protOG, "OE"), allOccurences)
        gOFF = setdiff(gOFF, c(paste0(protOG, "OE"), allOccurences))
      } else if(mutationt %in% c("delete", "deplete")){
        nodes = c(apply(expand.grid(prott, paste0("A_", locs)), 1, paste, collapse=""),apply(expand.grid(prott, paste0("L_", locs)), 1, paste, collapse=""))
        nodes = intersect(nodes, model$genes) # remove SPBs if protein can't localise there
        gOFF = c(gOFF, nodes)
        gON = setdiff(gON, nodes)
      } else if( mutationt == "KD"){
        nodes = apply(expand.grid(prott, paste0("A_", locs)), 1, paste, collapse="")
        nodes = intersect(nodes, model$genes)
        gOFF = c(gOFF, nodes)
        gON = setdiff(gON, nodes)
      } else if( mutationt == "hyperactive"){
        nodes = apply(expand.grid(prott, paste0("A_", locs)), 1, paste, collapse="")
        nodes = intersect(nodes, model$genes)
        gON = c(gON, nodes)
        gOFF = setdiff(gOFF, nodes)
      }else if(grepl("!", mutationt)){
        modmutationt = sub("!","", mutationt)
        if(modmutationt == "SPB"){
          nodes = apply(expand.grid(prott, paste0("L_", c("dSPB", "mSPB"))), 1, paste, collapse="")
        } else{
          nodes = apply(expand.grid(prott, paste0("L_", modmutationt)), 1, paste, collapse="")
        }
        gOFF = c(gOFF, nodes)
        gON = setdiff(gON, nodes)
      } else if(mutationt == "Cytoplasm+Bud"){
        nodes = apply(expand.grid(prott, paste0("L_", c("Cytoplasm", "Bud"))), 1, paste, collapse="")
        gON = c(gON, nodes)
        gOFF = setdiff(gOFF, nodes)
      } else if(mutationt == "mSPB"){
        nodes = apply(expand.grid(prott, paste0("L_", "mSPB")), 1, paste, collapse="")
        gON = c(gON, nodes)
        gOFF = setdiff(gOFF, nodes)
      } else{
        if(mutationt == "SPB"){
          nodes = apply(expand.grid(prott, paste0("L_", c("dSPB", "mSPB"))), 1, paste, collapse="")
        } else{
          nodes = apply(expand.grid(prott, paste0("L_", mutationt)), 1, paste, collapse="")
        }
        gON = c(gON, nodes)
        gOFF = setdiff(gOFF, nodes)
      }
    }
  }
  # find synchronous SSs
  gON = intersect(gON, model$genes) # make sure all gON/gOFF are really nodes
  gOFF = intersect(gOFF, model$genes)
  attract = getAttractors(model, method = "sat.exhaustive", genesON = gON, genesOFF = gOFF)#, startStates = 100)
  
  pdf(file=NULL)
  plot = plotAttractors(attract, allInOnePlot = TRUE)
  dev.off()
  
  # Now fix to specific CCstage
  if(CCstage == "M"){
    gON = c(gON, "SACA_Nucleus")
    gOFF = c(gOFF, "SpindleAlign")
  } else if(CCstage == "EA"){
    gOFF = c(gOFF, "SpindleAlign","SACA_Nucleus")
  } else if(CCstage == "LA"){
    gON = c(gON, "SpindleAlign")
    gOFF = c(gOFF, "SACA_Nucleus")
  }else if(CCstage == "NP"){
    gON = c(gON, "SACA_Nucleus", "SpindleAlign")
  } else {print("hol'up")}
  
  gON = intersect(gON, model$genes) # make sure all gON/gOFF are really nodes
  gOFF = intersect(gOFF, model$genes)
  
  fixmodel = fixGenes(model, c(gON,gOFF), c(rep(1,length(gON)), rep(0, length(gOFF))))
  
  IClist = ICs[model$genes]
  IClist[gON] = 1 # make ICs match restrictions
  IClist[gOFF] = 0
  
  exits = c()
  times = c()
  for(i in 1:N){
    flag = 0
    time = 1
    currentState = IClist
    while(flag==0){
      currentState = stateTransition(network = fixmodel, state = currentState, type = "asynchronous")
      time = time+1
      if(1%in% names(plot)){ # in cases where no steady states
        CC = compareColumns(df = plot$`1`, vec = currentState)
        if(currentState[which(model$genes=="Cdc14highL_Nucleus")] == 1){
          exit = "FEAR"
          flag = 1
        } else if(CC$result){
          exit = CC$col
          flag = 1
        }
        else if(time == timeMax){
          exit = "time"
          flag = 1
        }
      } else{
        if(currentState[which(model$genes=="Cdc14highL_Nucleus")] == 1){
          exit = "FEAR"
          flag = 1
        } else if(time == timeMax){
          exit = "time"
          flag = 1
        }
      }
    }
    exits = c(exits, exit)
    times = c(times, time)
  }
  return(list(exits=exits, times=times, SSs = plot))
}
