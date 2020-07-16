# CreateCLM.R
# This script makes a CLM from activity (MEN, FEAR), localization (Loc), and localization-specific (LocSpec) networks in BoolNet format in .csv files. 
library(tidyverse)
library(BoolNet)
library(stringr)
source("CLMfunctions.R") # note this dependency


createNetwork = function(MEN, Loc, FEAR, LocSpec, netname, directory){
  # This function takes in the four input files (.csv format) and writes a network file and a unreg nodes file. 
  ## INPUTS
  # MEN - MEN activity network (BoolNet format, no localization information)
  # Loc - Localization network (BoolNet format, X.loc means the node representing localization of protein X in location loc, XA means protein X activity in the relevant compartment)
  # FEAR - FEAR activity network (BoolNet format, no localization information)
  # LocSpec - Localization specific rules (BoolNet format, XL.loc means the node representing localization of protein X in location loc, XA.loc means protein X activity in location loc)
  # netname - name to give files
  # directory - name of directory containing .csv files
  
  ## OUTPUTS
  # This script will write 4 .csv files into the directory it is run in (not "directory")
  # 1. netname.csv basic model with no OE or forced localization (BoolNet format)
  # 2. netnameOE.csv model with OE (BoolNet format)
  # 3. netnameFL.csv model with forced localization at SPB (BoolNet format)
  # 4. netnameunreg.csv list of nodes that are not present but not regulated, this is used to improve efficiency of attractor search
  
  if(getwd()!= directory){ # set directories
    oldpathway = getwd()
    setwd(directory)
  }
  #---------------------------------------------------------------------------------------------------------------------------------------------
  # Load files and prepare data frames
  
  MENNet = read_csv(MEN)
  LocNet = read_csv(Loc)
  FEARNet = read_csv(FEAR)
  LocSpecNet = read_csv(LocSpec) 
  NonPhys = c("Tem1_high","Bub2_low","Bfa1_low") # Special rules for non-physiological conditions for Tem1, Bub2, Bfa1

  ActNet = rbind(MENNet, FEARNet)

  Locs = t(data.frame(strsplit(LocNet$targets, split = "\\."))) # need "\\." as "." is a wildcard!!!
  colnames(Locs) = c("Prot", "Localization")
  rownames(Locs)  =c()
  Locs = as.tibble(Locs)

  Prots = ActNet$targets
  OEnodes = paste0(Prots, "_OE")

  Cyto = Prots
  SPB = filter(Locs, Localization == "SPB")$Prot # only proteins with an SPB localization rule will have an SPB localization node created
  Nucleolus = filter(Locs, Localization == "Nucleolus")$Prot # only proteins with a nucleolus localization rule will have a nucleolus localization node created
  Bud = Cyto
  Nucleus = Prots

  TheoLocs = list() # theoretical localizations refer to categories in localization network (cytoplasm, nucleus, SPB, nucleolus)
  ProtLocs = list() # protein localizations refer to categories in CLM (cytoplasm, bud, nucleus, mSPB, dSPB, nucleolus)
  for(j in Prots){
    tvec = c()
    tpvec = c()
    if(j %in% Cyto){
      tvec = c(tvec, "Cytoplasm")
      tpvec = c(tpvec, "Cytoplasm", "Bud")
    }
    if(j %in% SPB){
      tvec = c(tvec, "SPB")
      tpvec = c(tpvec, "mSPB", "dSPB")
    }
    if(j %in% Nucleus){
      tvec = c(tvec, "Nucleus")
      tpvec = c(tpvec, "Nucleus")
    }
    if(j %in% Nucleolus){
      tvec = c(tvec, "Nucleolus")
      tpvec = c(tpvec, "Nucleolus")
    }
    TheoLocs[[j]] = tvec
    ProtLocs[[j]] = tpvec
  }

  # Create list of nodes
  LocNodes = c(paste0(Cyto, "L.Cytoplasm"),paste0(Bud, "L.Bud"),paste0(Nucleus, "L.Nucleus"),paste0(Nucleolus, "L.Nucleolus"),paste0(SPB, "L.mSPB"),paste0(SPB, "L.dSPB"))
  ActNodes = c(paste0(Cyto, "A.Cytoplasm"),paste0(Bud, "A.Bud"),paste0(Nucleus, "A.Nucleus"),paste0(Nucleolus, "A.Nucleolus"),paste0(SPB, "A.mSPB"),paste0(SPB, "A.dSPB"))
  
  #---------------------------------------------------------------------------------------------------------------------------------------------
  # make activity networks
  
  unRegA = c()
  EdgelistA = c()
    for(j in 1:length(Prots)){ # iterate through proteins to create rule for each
      tNode = Prots[j]
      tRule = ActNet[which(ActNet$targets==tNode),2]
      tRuleSplit = str_split(tRule, "\\|") # split around | (ORs)
      tRuleSplit = gsub("&","+",tRuleSplit[[1]], fixed=TRUE) # switch & to +
      tEdges = paste0(tRuleSplit, "=", tNode) # turn into edge description
      tEdges = gsub("[ ()]","", tEdges) # get rid of whitespace
      tLEdges = c()
      for(e in tEdges){
        tEdgeSplit = str_split(e, "=")
        tEdgeIns = str_split(tEdgeSplit[[1]][1], "\\+")[[1]]
        tEdgeIns2 = gsub("[!()]","",tEdgeIns) # get rid of brackets and !s
          for(L in unlist(ProtLocs[Prots[j]])){
            tRuleL = e
            for(k in 1:length(tEdgeIns2)){
              if(str_sub(tEdgeIns2[k],-3,-1)=="_OE"){
              # dont change anything
              } else if(L%in%ProtLocs[which(names(ProtLocs)==tEdgeIns2[k])][[1]]){ # [[1]] required because of list structure
                tRuleL = sub(tEdgeIns2[k],paste0(tEdgeIns2[k],"A.", L), tRuleL)
              } else if((L=="mSPB")&("Cytoplasm"%in%ProtLocs[which(names(ProtLocs)==tEdgeIns2[k])][[1]])){
                tRuleL = sub(tEdgeIns2[k],paste0(tEdgeIns2[k],"A.", "Cytoplasm"), tRuleL)
              } else if((L=="dSPB")&("Bud"%in%ProtLocs[which(names(ProtLocs)==tEdgeIns2[k])][[1]])){ # dSPB localization with bud components
                tRuleL = sub(tEdgeIns2[k],paste0(tEdgeIns2[k],"A.", "Bud+Spindle_Align"), tRuleL)
            }
        }
        tRuleL = paste0(tRuleL, "A.",L)
        if(identical(tEdgeIns2,tNode)){ # in case where no regulators, activity depends only on localization
          tRuleL = sub(paste0(tNode, "A.",L,"="),paste0(tNode, "L.",L,"="), tRuleL)
          unRegA = c(unRegA, paste0(tNode, "A.",L))
        } else{
          tRuleL = sub("=",paste0("+",tNode, "L.",L,"="), tRuleL) # make localization necessary for activation
        }
        tLEdges = c(tLEdges, tRuleL)
      }
    }
    EdgelistA = c(EdgelistA, tLEdges)
  }

  #---------------------------------------------------------------------------------------------------------------------------------------------
  # make localization network

  EdgelistL = c()

  unregNodes = rep(1, length(LocNodes))
  names(unregNodes) = LocNodes
  for(p in 1:length(Prots)){
    for(l in TheoLocs[[Prots[p]]]){ # iterate through theoretical localization for each protein 
      if(paste0(Prots[p], ".",l)%in%LocNet$targets){ # check if rule for this localization
        j = which(LocNet$targets == paste0(Prots[p], ".",l))
        tRule = LocNet[j,2]
        tRuleSplit = str_split(tRule, "\\|") # split around | (ORs)
        tRuleSplit = gsub("&","+",tRuleSplit[[1]], fixed=TRUE) # switch & to +
        htLoc = str_split(LocNet[j,1], "\\.")[[1]][2]
        tEdges = paste0(tRuleSplit, "=", LocNet[j,1]) # turn into edge description
        tEdges = gsub(" ","", tEdges) # get rid of whitespace
        if(htLoc == "Cytoplasm"){ # work out possible interpretations
          tLoc = c("Cytoplasm", "Bud")
        } else if(htLoc == "SPB"){
          tLoc = c("dSPB", "mSPB")
        } else if(htLoc == "Nucleolus"){tLoc = c("Nucleolus")
        } else {tLoc = c("Nucleus")}
        tLEdges = c()
        for(e in tEdges){
          tEdgeSplit = str_split(e, "=")
          tEdgeIns = str_split(tEdgeSplit[[1]][1], "\\+")[[1]]
          tEdgeIns2 = gsub("[!()]","",tEdgeIns) # this version does not have !s or ()s
          for(L in tLoc){
            tRuleL = e
            for(k in 1:length(tEdgeIns2)){
              if(grepl("\\.", tEdgeIns2[k])){ # if localized version
                regLoc = substr(tEdgeIns2[k],regexpr("\\.",tEdgeIns2[k])[1]+1,100) # finds text after .
                if(regLoc == "Cytoplasm"){ # work out possible interpretations
                  if(L %in% c("Cytoplasm", "mSPB", "Nucleus")){
                    tRuleL = sub(tEdgeIns2[k], paste0(str_split(tEdgeIns2[k], "\\.")[[1]][1],"L.Cytoplasm"), tRuleL)
                  } else if(L == "dSPB"){
                    tRuleL = sub(tEdgeIns2[k], paste0(str_split(tEdgeIns2[k], "\\.")[[1]][1],"L.Bud+Spindle_Align"), tRuleL)   ## some more bud-> dSPB interactions            
                  }
                  else{
                    tRuleL = sub(tEdgeIns2[k], paste0(str_split(tEdgeIns2[k], "\\.")[[1]][1],"L.Bud"), tRuleL) 
                  }
                } else if(regLoc == "SPB"){ # SPB localization should only regulate other SPB localizations
                  if(L=="mSPB"){
                    tRuleL = sub(tEdgeIns2[k], paste0(str_split(tEdgeIns2[k], "\\.")[[1]][1],"L.mSPB"), tRuleL)
                  } else{ # only other case should be dSPB
                    tRuleL = sub(tEdgeIns2[k], paste0(str_split(tEdgeIns2[k], "\\.")[[1]][1],"L.dSPB"), tRuleL)
                  }
                } else if(regLoc == "Nucleolus"){ # Only one rule for nucleolus
                  tRuleL = sub(tEdgeIns2[k], paste0(str_split(tEdgeIns2[k], "\\.")[[1]][1],"L.Nucleolus"), tRuleL)
                } else{ # only other possibility is nuclear protein, which has only one version
                  tRuleL = sub(tEdgeIns2[k], paste0(str_split(tEdgeIns2[k], "\\.")[[1]][1],"L.Nucleus"), tRuleL)
                }
                
              } else{ # activity version
                if(L%in%ProtLocs[which(paste0(names(ProtLocs),"A")==tEdgeIns2[k])][[1]]){ # [[1]] required because of list structure
                  tRuleL = sub(tEdgeIns2[k],paste0(tEdgeIns2[k],".", L), tRuleL)
                } else if((L=="mSPB")&("Cytoplasm"%in%ProtLocs[which(paste0(names(ProtLocs),"A")==tEdgeIns2[k])][[1]])){
                  tRuleL = sub(tEdgeIns2[k],paste0(tEdgeIns2[k],".", "Cytoplasm"), tRuleL)
                } else if((L=="dSPB")&("Bud"%in%ProtLocs[which(paste0(names(ProtLocs),"A")==tEdgeIns2[k])][[1]])){
                  tRuleL = sub(tEdgeIns2[k],paste0(tEdgeIns2[k],".", "Bud+Spindle_Align"), tRuleL) #  spindle_alignment required for bud->dSPB regulation 
                }
              }
              tRuleL = sub(LocNet[j,1], paste0(str_split(LocNet[j,1], "\\.")[[1]][1],"L.",L), tRuleL) # switch out hypothetical for specific
            }
            tRuleL = gsub("[()]", "", tRuleL)
            tLEdges = c(tLEdges, tRuleL)
          }
        }
        EdgelistL = c(EdgelistL, tLEdges)
        unregNodes[paste0(Prots[p],"L.",tLoc)] = 0 # these nodes are regulated and so are removed
        # self-regulated nodes are readded to the unreg list
        if(tLEdges == paste0(Prots[p],"L.",tLoc,"=",Prots[p],"L.",tLoc)){
          unregNodes[paste0(Prots[p],"L.",tLoc)] = 1
        }
      } 
    }
  }
  
  #---------------------------------------------------------------------------------------------------------------------------------------------
  # Localization specific regulation

  EdgelistLS = c()

  for(j in 1:nrow(LocSpecNet)){ # iterate though rules on Localization-specific list
    htLoc = str_split(LocSpecNet$targets[j], "\\.")[[1]][2]
    if(htLoc == "Cytoplasm"){ # work out possible interpretations
      tLoc = c("Cytoplasm", "Bud")
    } else if(htLoc == "SPB"){
      tLoc = c("dSPB", "mSPB")
    } else if(htLoc == "Nucleolus"){ tLoc = c("Nucleolus")
    } else {tLoc = c("Nucleus")}
    tRule = LocSpecNet[j,2]
    tRuleSplit = str_split(tRule, "\\|") # split around | (ORs)
    tRuleSplit = gsub("&","+",tRuleSplit[[1]], fixed=TRUE) # switch & to +
    tEdges = paste0(tRuleSplit, "=", LocSpecNet[j,1]) # turn into edge description
    tEdges = gsub(" ","", tEdges) # get rid of whitespace
    tLEdges = c()
    for(e in tEdges){
      tEdgeSplit = str_split(e, "=")
      tEdgeIns = str_split(tEdgeSplit[[1]][1], "\\+")[[1]]
      tEdgeIns2 = gsub("[!()]","",tEdgeIns) # this version does not have !s or ()s
      for(L in tLoc){
        tRuleL = e
        for(k in 1:length(tEdgeIns2)){
          if(grepl("\\.", tEdgeIns2[k])){ # if localized version
            regLoc = substr(tEdgeIns2[k],regexpr("\\.",tEdgeIns2[k])[1]+1,100) # finds text after .
            if(regLoc == "Cytoplasm"){ # work out possible interpretations
              if(L %in% c("Cytoplasm", "mSPB", "Nucleus")){
                tRuleL = sub(tEdgeIns2[k], paste0(str_split(tEdgeIns2[k], "\\.")[[1]][1],".Cytoplasm"), tRuleL)
              } else if(L == "dSPB"){
                tRuleL = sub(tEdgeIns2[k], paste0(str_split(tEdgeIns2[k], "\\.")[[1]][1],".Bud+Spindle_Align"), tRuleL)   ## some more bud-> dSPB interactions            
              } 
              else{
                tRuleL = sub(tEdgeIns2[k], paste0(str_split(tEdgeIns2[k], "\\.")[[1]][1],".Bud"), tRuleL) 
              }
            } else if(regLoc == "SPB"){ 
              if(L%in%c("mSPB", "Cytoplasm")){
                tRuleL = sub(tEdgeIns2[k], paste0(str_split(tEdgeIns2[k], "\\.")[[1]][1],".mSPB"), tRuleL)
              } else{ # only other case should be dSPB
                tRuleL = sub(tEdgeIns2[k], paste0(str_split(tEdgeIns2[k], "\\.")[[1]][1],".dSPB"), tRuleL)
              }
            } else if(regLoc == "Nucleolus"){
              tRuleL = sub(tEdgeIns2[k], paste0(str_split(tEdgeIns2[k], "\\.")[[1]][1],".Nucleolus"), tRuleL)
            } else{ # only other possibility is nuclear protein, which has only one version
              tRuleL = sub(tEdgeIns2[k], paste0(str_split(tEdgeIns2[k], "\\.")[[1]][1],".Nucleus"), tRuleL)
            }
          
          } else{ # activity version
            if(L%in%ProtLocs[which(paste0(names(ProtLocs),"A")==tEdgeIns2[k])][[1]]){ # [[1]] required because of list structure
              tRuleL = sub(tEdgeIns2[k],paste0(tEdgeIns2[k],".", L), tRuleL)
            } else if((L=="mSPB")&("Cytoplasm"%in%ProtLocs[which(paste0(names(ProtLocs),"A")==tEdgeIns2[k])][[1]])){
              tRuleL = sub(tEdgeIns2[k],paste0(tEdgeIns2[k],".", "Cytoplasm"), tRuleL)
            } else if((L=="dSPB")&("Bud"%in%ProtLocs[which(paste0(names(ProtLocs),"A")==tEdgeIns2[k])][[1]])){
              tRuleL = sub(tEdgeIns2[k],paste0(tEdgeIns2[k],".", "Bud"), tRuleL)
            }
            tRuleL = sub("=", paste0("+",sub("A\\.","L.",sub(htLoc, tLoc[1], LocSpecNet$targets[j])),"="), tRuleL) # add localization requirement to activity rule
          }
          tRuleL = sub(LocSpecNet[j,1], paste0(str_split(LocSpecNet[j,1], "\\.")[[1]][1],".",L), tRuleL) # switch out hypothetical for specific
        }
        tRuleL = gsub("[()]", "", tRuleL)
        tLEdges = c(tLEdges, tRuleL)
        if(L=="Nucleus"){
          tRuleL2 = gsub("Cytoplasm", "Bud", tRuleL)
          if(!identical(tRuleL, tRuleL2)){
            tLEdges = c(tLEdges, tRuleL2)
          }
        }
        if(grepl("A.",LocSpecNet$targets[j])){ # if activity node described here can remove self-activation
          if(length(tLoc)>1){
             EdgelistA = setdiff(EdgelistA, c(paste0(sub("A\\.","L.",sub(htLoc, tLoc[1], LocSpecNet$targets[j])), "=", sub(htLoc, tLoc[1], LocSpecNet$targets[j])),paste0(sub("A\\.","L.",sub(htLoc, tLoc[2], LocSpecNet$targets[j])), "=", sub(htLoc, tLoc[2], LocSpecNet$targets[j]))))
          } else{
             EdgelistA = setdiff(EdgelistA, c(paste0(sub("A\\.","L.",LocSpecNet$targets[j]), "=",LocSpecNet$targets[j])))
          }
        }
      }
    }
    EdgelistLS = c(EdgelistLS, tLEdges)
    if(length(tLoc)>1){ 
      EdgelistL = setdiff(EdgelistL, c(paste0(sub(htLoc, tLoc[1], LocSpecNet$targets[j]),"=",sub(htLoc, tLoc[1], LocSpecNet$targets[j])), paste0(sub(htLoc, tLoc[2], LocSpecNet$targets[j]),"=",sub(htLoc, tLoc[2], LocSpecNet$targets[j]))))
      unregNodes[c(sub(htLoc, tLoc[1], LocSpecNet$targets[j]), sub(htLoc, tLoc[2], LocSpecNet$targets[j]))] = 0
    } else{
      EdgelistL = setdiff(EdgelistL, c(paste0(sub(htLoc, tLoc, LocSpecNet$targets[j]),"=",sub(htLoc, tLoc, LocSpecNet$targets[j]))))
      unregNodes[LocSpecNet$targets[j]] = 0
    }
  }



  unregEdges = paste0(LocNodes[which(as.logical(unregNodes))],"=",LocNodes[which(as.logical(unregNodes))])
  EdgelistL = c(EdgelistL, unregEdges)

  unRegL = names(unregNodes)[which(as.logical(unregNodes))]
  # Cdc14 and Cdc20 rules are specified in algorithm so treated differently
  unRegL = setdiff(unRegL, c("Cdc14_lowL.Cytoplasm", "Cdc14_highL.Cytoplasm", "Cdc14_lowL.Bud", "Cdc14_highL.Bud","APC_Cdc20L.Cytoplasm","APC_Cdc20L.Bud"))

  #---------------------------------------------------------------------------------------------------------------------------------------------
  # fix dSPB regulation - find mSPB rules and split into pre- and post-spindle align rules
  # Note there is some redundancy here but it is fixed at the end by selecting unique rules

  ExtraL = c()
  for(j in 1:length(EdgelistL)){
    if(tail(str_split(EdgelistL[j],"\\.")[[1]], n=1)=="mSPB"){ # find edges targeting proteins at the mSPB
      tRule = gsub("mSPB","dSPB", EdgelistL[j]) # switch mSPBs to dSPBs
      if(!(tRule%in%EdgelistL)){ # in order not to add duplicate edges for dSPB-dSPB regulation
        tRule = gsub("=", "+!Spindle_Align=", tRule)# AND in Spindle_Align
        ExtraL = c(ExtraL, tRule)
      }
    }
  }

  EdgelistL = c(EdgelistL, ExtraL)

  ExtraA = c()
  for(j in 1:length(EdgelistA)){
    if(tail(str_split(EdgelistA[j],"\\.")[[1]], n=1)=="mSPB"){
      tRule = gsub("mSPB","dSPB", EdgelistA[j]) # switch mSPBs to dSPBs
      if(!(tRule%in%EdgelistA)){ # in order not to add duplicate edges for dSPB-dSPB regulation
        tRule = gsub("=", "+!Spindle_Align=", tRule)# AND in Spindle_Align
        ExtraA = c(ExtraA, tRule)
      }
    }
  }

  EdgelistA = c(EdgelistA, ExtraA)

  ExtraLS = c()
  for(j in 1:length(EdgelistLS)){
    if(tail(str_split(EdgelistLS[j],"\\.")[[1]], n=1)=="mSPB"){ # find edges targeting proteins at the mSPB
      tRule = gsub("mSPB","dSPB", EdgelistLS[j]) # switch mSPBs to dSPBs
      if(!(tRule%in%EdgelistLS)){ # in order not to add duplicate edges for dSPB-dSPB regulation
        tRule = gsub("=", "+!Spindle_Align=", tRule)# AND in Spindle_Align
        ExtraLS = c(ExtraLS, tRule)
      }
    }
  }

  EdgelistLS = c(EdgelistLS, ExtraLS)
  
  EdgelistAll = c(EdgelistA, EdgelistL, EdgelistLS)
  EdgelistSplit = str_split(EdgelistAll,"=", simplify = TRUE)

  #---------------------------------------------------------------------------------------------------------------------------------------------
  # Add OE edges
  newedges = c()
  oldedges = c()
  activatorsAll = list()
  inhibitorsAll = list()
  for(n in ActNodes){
    reledges = EdgelistSplit[which(EdgelistSplit[,2]==n),] # all edges into n
    if(length(which(EdgelistSplit[,2]==n))==1){
      reledges = as.data.frame(matrix(reledges, ncol = 2, nrow = 1))
    }
    reledgesNS = EdgelistAll[which(EdgelistSplit[,2]==n)]
    inputs = c()
    for(j in 1:nrow(reledges)){
      inputs = c(inputs,str_split(gsub("!","",reledges[j,1]),"\\+")[[1]]) # give raw input nodes names
    }
    inputs = unique(inputs)
    inputs = setdiff(inputs, c(sub("A\\.", "L.",n), "Spindle_Align")) # remove localization node and spindle align
    inputs2 = c()
    for(i in inputs){
      gene = str_sub(i, start = 1, end = (gregexpr("\\.", i)[[1]][1]-2)) # strip away location etc
      if(grepl("_high",gene)){ # get rid of any high/low
        gene = str_sub(gene, start = 1, end = -6)
      } else if(grepl("_low",gene)){
        gene = str_sub(gene, start = 1, end = -5)
      }
      inputs2 = c(inputs2, gene)
    }
    inputs2 = unique(inputs2)
    inhibitors = c()
    for(m in inputs2){ # Note this is a condition for being an inhibitor, this works as long as no cases of mixed activation inhibition.
      if(any(grepl(paste0("!",m), reledges[,1]))){
        inhibitors = c(inhibitors,m)
      }
    }
    activators = setdiff(inputs2, inhibitors)
    inhibitorsAll[[n]] = inhibitors
    activatorsAll[[n]] = activators
    for(a in activators){
      if(!(paste0(a,"+",sub("A\\.", "L_",n))%in%reledges[,1])){ # if not just direct activation
        newedges = c(newedges, paste0(a, "_OE+",sub("A\\.", "L.",n), "=",n)) # AND in localization to OE edge
      }
    }
    if(nrow(reledges)>1){ # if multiple edges in
      newedgesI = reledgesNS
      np = grep(paste(NonPhys,collapse="|"), newedgesI) # find which edges contain NonPhys states which behave differently
      edit = setdiff(1:length(newedgesI),np)
      for(i in inhibitors){
        newedgesI[edit] = gsub("=", paste0("+!",i, "_OE="), newedgesI[edit]) 
      }
      oldedges = c(oldedges,reledgesNS[edit])
      newedges = c(newedges, newedgesI[edit])
    }
    
  }
  
  High = grep("high", names(activatorsAll), value = TRUE)
  
  for(n in High){
    if(!all(inhibitorsAll[[n]]%in%inhibitorsAll[[sub("high","low",n)]])){ #Check if all inhibitors of high also inhibit low, if not:
      extraInhibitors = setdiff(inhibitorsAll[[n]],inhibitorsAll[[sub("high","low",n)]])
      lown = sub("high","low",n)
      reledges = EdgelistSplit[which(EdgelistSplit[,2]==lown),]
      if(length(which(EdgelistSplit[,2]==lown))==1){
        reledges = as.data.frame(matrix(reledges, ncol = 2, nrow = 1))
      }
      reledgesNS = EdgelistAll[which(EdgelistSplit[,2]==lown)]
      newedgesI = reledgesNS
      for(i in extraInhibitors){
        newedgesI = gsub("=", paste0("+!",i, "_OE="), newedgesI) 
      }
      oldedges = c(oldedges,reledgesNS)
      newedges = c(newedges, newedgesI)
    }
    if(!all(activatorsAll[[sub("high","low",n)]]%in%activatorsAll[[n]])){ #Check if all activators of low also activate high, if not:
      extraActivators = setdiff(activatorsAll[sub("high","low",n)],activatorsAll[n])
      for(a in extraActivators){
        newedges = c(newedges, paste0(a, "_OE+",sub("A\\.", "L.",n), "=",n))
      }
    }
  }
  
  #---------------------------------------------------------------------------------------------------------------------------------------------
  # Add OE loc edges

  for(n in LocNodes){
    reledges = EdgelistSplit[which(EdgelistSplit[,2]==n),] # all edges into n
    if(length(which(EdgelistSplit[,2]==n))==1){
      reledges = as.data.frame(matrix(reledges, ncol = 2, nrow = 1))
    }
    reledgesNS = EdgelistAll[which(EdgelistSplit[,2]==n)]
    inputs = c()
    for(j in 1:nrow(reledges)){
      inputs = c(inputs,str_split(gsub("!","",reledges[j,1]),"\\+")[[1]]) # give raw input nodes names
    }
    inputs = unique(inputs)
    inputs = setdiff(inputs, c("Spindle_Align")) # remove spindle align
    inputs2 = c()
    for(i in inputs){
      gene = str_sub(i, start = 1, end = (gregexpr("\\.", i)[[1]][1]-2)) # strip away location etc
      if(grepl("_high",gene)){ # get rid of any high/low
        gene = str_sub(gene, start = 1, end = -6)
      } else if(grepl("_low",gene)){
        gene = str_sub(gene, start = 1, end = -5)
      }
      inputs2 = c(inputs2, gene)
    }
    inputs2 = unique(inputs2)
    inhibitors = c()
    for(m in inputs2){ # Note this is a condition for being an inhibitor, I dont think there are any cases of mixed inhibition/activation so should be fine.
      if(any(grepl(paste0("!",m), reledges[,1]))){
        inhibitors = c(inhibitors,m)
      }
    }
    activators = setdiff(inputs2, inhibitors)
    inhibitorsAll[[n]] = inhibitors
    activatorsAll[[n]] = activators
    for(a in activators){
      if(!grepl(a, n)){ # if not just direct activation
        locInputs = grep("L.", inputs, value = TRUE) # inputs from localization nodes
        if(length(inhibitors)>0){
          for(j in 1:length(inhibitors)){
            locInh = grep(inhibitors[j],locInputs, value = TRUE)
            if(length(locInh)>0){
              locInputs = setdiff(locInputs, locInh)
              locInputs = c(locInputs, paste0("!",locInh))
          }
          }
        }
        if(length(locInputs>0)){
          newedges = c(newedges, paste0(a, "_OE+",concatenate(locInputs, "+"), "=",n)) # add in localization requirements
        } else{
          newedges = c(newedges, paste0(a, "_OE=",n)) # no need to include any localization
        }
      }
    }
    if(nrow(reledges)>1){ # if multiple edges in
      newedgesI = reledgesNS
      np = grep(paste(NonPhys,collapse="|"), newedgesI) # find which edges contain NonPhys states
      edit = setdiff(1:length(newedgesI),np)
      for(i in inhibitors){
        newedgesI[edit] = gsub("=", paste0("+!",i, "_OE="), newedgesI[edit]) 
      }
      oldedges = c(oldedges,reledgesNS[edit])
      newedges = c(newedges, newedgesI[edit])
    }
  }
  
  High = grep("high", names(activatorsAll), value = TRUE)
  
  for(n in High){
    if(!all(inhibitorsAll[[n]]%in%inhibitorsAll[[sub("high","low",n)]])){ #Check if all inhibitors of high also inhibit low, if not:
      extraInhibitors = setdiff(inhibitorsAll[[n]],inhibitorsAll[[sub("high","low",n)]])
      lown = sub("high","low",n)
      reledges = EdgelistSplit[which(EdgelistSplit[,2]==lown),]
      if(length(which(EdgelistSplit[,2]==lown))==1){
        reledges = as.data.frame(matrix(reledges, ncol = 2, nrow = 1))
      }
      reledgesNS = EdgelistAll[which(EdgelistSplit[,2]==lown)]
      newedgesI = reledgesNS
      for(i in extraInhibitors){
        newedgesI = gsub("=", paste0("+!",i, "_OE="), newedgesI) 
      }
      oldedges = c(oldedges,reledgesNS)
      newedges = c(newedges, newedgesI)
    }
    if(!all(activatorsAll[[sub("high","low",n)]]%in%activatorsAll[[n]])){ #Check if all activators of low also activate high, if not:
      extraActivators = setdiff(activatorsAll[sub("high","low",n)],activatorsAll[n])
      for(a in extraActivators){
        newedges = c(newedges, paste0(a, "_OE+",sub("A\\.", "L.",n), "=",n))
      }
    }
  }
  
  # Now make SPB localization models (in these localization treated as local overexpression)
  mSPBOEedgesold = gsub("OE", "FLmSPB",oldedges[which(str_sub(oldedges, start = -6, end = 100)=="L.mSPB")])
  dSPBOEedgesold = gsub("OE", "FLdSPB",oldedges[which(str_sub(oldedges, start = -6, end = 100)=="L.dSPB")])
  
  mSPBOEedgesnew = gsub("OE", "FLmSPB",newedges[which(str_sub(newedges, start = -6, end = 100)=="L.mSPB")])
  dSPBOEedgesnew = gsub("OE", "FLdSPB",newedges[which(str_sub(newedges, start = -6, end = 100)=="L.dSPB")])
  
  #---------------------------------------------------------------------------------------------------------------------------------------------
  # Combine edgelists to make network
  
  EdgelistAll2 = c(setdiff(EdgelistAll,oldedges), newedges) # OE version of the edgelist
  EdgelistAll3 = c(setdiff(EdgelistAll,c(mSPBOEedgesold,dSPBOEedgesold)), c(mSPBOEedgesnew,dSPBOEedgesnew)) # loc version of the edgelist
  
  EdgelistAll = c(EdgelistAll, "Cdc14_highA.Cytoplasm=ME") # add in ME rule
  EdgelistAll2 = c(EdgelistAll2, "Cdc14_highA.Cytoplasm=ME") 
  EdgelistAll3 = c(EdgelistAll3, "Cdc14_highA.Cytoplasm=ME") 
  

  if(getwd()!= oldpathway){setwd(oldpathway)}
  LocNodes2 = gsub("_","",LocNodes)
  ActNodes2 = gsub("_","",ActNodes)
  EdgelistAll = gsub("_","",EdgelistAll)
  EdgelistAll2 = gsub("_","",EdgelistAll2)
  EdgelistAll3 = gsub("_","",EdgelistAll3)
  LocNodes2 = gsub("\\.","_",LocNodes2)
  ActNodes2 = gsub("\\.","_",ActNodes2)
  EdgelistAll = gsub("\\.","_",EdgelistAll)
  EdgelistAll2 = gsub("\\.","_",EdgelistAll2)
  EdgelistAll3 = gsub("\\.","_",EdgelistAll3)
  allOE2 = gsub("_","",OEnodes)
  
  convertToBoolNet(c(LocNodes2, ActNodes2, "SpindleAlign","ME"), EdgelistAll, paste0(netname,".txt")) # base model
  convertToBoolNet(c(LocNodes2, ActNodes2, "SpindleAlign","ME"), EdgelistAll2, paste0(netname,"OE.txt")) # model with overexpression nodes
  convertToBoolNet(c(LocNodes2, ActNodes2, "SpindleAlign","ME"), EdgelistAll3, paste0(netname,"SPB.txt")) # model with SPB localization
  
  
  allunReg = gsub("_","",unRegL) # unRegA not included as activity always depends on localization
  allunReg2 = gsub("\\.","_",allunReg)

  allunReg2 = setdiff(allunReg2, c("SACA_Nucleus"))

  write_csv(tibble(node = allunReg2), paste0(netname, "unreg.csv"))
}

