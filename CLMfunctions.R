library(tidyverse)
library(BoolNet)
library(stringr)


reducePath = function(path){
  redPath = path[, which((!colSums(path)%in%c(0,nrow(path))))]
}

reduceTable = function(path){
  redPath = path[which((!rowSums(path)%in%c(0,ncol(path)))),]
}

inorNA = function(x, list){ # method to check whether an element is in a list or is NA (either true NA or "NA" string)
  if(is.na(x)|x=="NA"){
    out = NA
  } else {
    out = as.numeric(x%in%list)
  }
  return(out)
}

identicalNA = function(x, y){ # method to compare strings which contain NAs
  na1 = which(is.na(x))
  na2 = which(is.na(y))
  comps = setdiff(1:length(x), c(na1, na2))
  return(identical(x[comps], y[comps]))
}

compareColumns = function(vec, df){ # check if a vector is a column of a dataframe (check if model state is a steady state)
  result = FALSE
  col = NA
  for(j in 1:ncol(df)){
    if(all(unlist(df[,j])==vec)){
      result = TRUE
      col = j
    }
  }
  return(list(result = result,col = col))
}

convertToBoolNet = function(nodes, edges, filename){ # method to convert an edgelist into a BoolNet format file
  lines = c("targets, factors")
  inOut = t(data.frame(str_split(edges, "=")))
  for(j in 1:length(nodes)){
    inputEdges = which(inOut[,2]==nodes[j])
    #print(j)
    if(length(inputEdges)==0){
      # may need to do something else here 
      lines = c(lines, paste0(nodes[j],", ",nodes[j]))
    }else{
      tline = paste0(nodes[j],", ")
      for(k in 1:length(inputEdges)){
        spl = str_split(inOut[inputEdges[k],1],"\\+")
        #print(k)
        if(length(spl[[1]])>1){ 
          tline = paste0(tline, "(", gsub("!","! ",concatenate(spl[[1]]," & ")),") | ")
        } else{ # if no AND gate
          tline = paste0(tline, gsub("!","! ",spl[[1]][1])," | ")
        }
      }
      tline = str_sub(tline, end = -3)# get rid of final | symbol
      lines = c(lines, tline)
    }
  }
  writeLines(lines, filename)
  
}

concatenate = function(vec, fill){ # concatenate a vector with a string (fill) between each entry
  N = length(vec)
  result = vec[1]
  if(N>1){
    for(i in 2:N){
      result = paste0(result, fill, vec[i])
    }
  }
  return(result)
}

makeICs = function(AIC, LIC, model, directory){ # create a vector specifying initial conditions of nodes in the CLM from Activity (AIC) and Localization (LIC) initial conditions
  if(getwd()!= directory){ # set directories
    oldpathway = getwd()
    setwd(directory)
  }
  metaAIC = read_csv(AIC)
  metaLIC = read_csv(LIC)
  metaICs = c()
  for(node in model$genes){
    if(!(node%in%c("SpindleAlign", "ME"))&!(str_sub(node,-2,-1)=="OE")&!(str_sub(node,-6,-1)=="FLdSPB")&!(str_sub(node,-6,-1)=="FLmSPB")){
      tLoc = tail(str_split(node,"_")[[1]],n=1) # find specific localisation
      tRemain = sub(paste0("_",tLoc), "", node)
      tType = str_sub(tRemain, -1, -1) # just last letter
      tName = str_sub(tRemain, 1, -2) # last letter removed
      Table = eval(parse(text = paste0("meta",tType,"IC")))
      tState = Table[[which(Table$protein == tName), which(colnames(Table)==tLoc)]]
    } else{
      tState = 0
    }
    metaICs = c(metaICs, tState)
  }
  names(metaICs) = model$genes
  if(getwd()!= oldpathway){setwd(oldpathway)}
  return(metaICs)
}


extractCdc14_2 = function(data){ # get Cdc14 localization data from steady state data
  cdc14Data = as_tibble(t(data[c("SACA_Nucleus", "SpindleAlign",paste0("Cdc14high","L_",c("Nucleus", "Cytoplasm", "Bud")),paste0("Cdc14low","L_",c("Nucleus", "Cytoplasm", "Bud")), "ME"),])) #transverse to make easier to use dplyr methods
  cdc14Data = mutate(cdc14Data, Cdc14L_Nucleus = Cdc14lowL_Nucleus+Cdc14highL_Nucleus, Cdc14L_Cytoplasm = Cdc14lowL_Cytoplasm+Cdc14highL_Cytoplasm, Cdc14L_Bud = Cdc14lowL_Bud+Cdc14highL_Bud, ME = ME*3, SpindleAlign = 3*SpindleAlign, SACA_Nucleus = 3*SACA_Nucleus)
  cdc14Data = select(cdc14Data, SACA_Nucleus, SpindleAlign,paste0("Cdc14","L_",c("Nucleus", "Cytoplasm", "Bud")), ME)
  cdc14Data[cdc14Data==0] = "OFF"
  cdc14Data[cdc14Data==1] = "LOW"
  cdc14Data[cdc14Data==2] = "HIGH"
  cdc14Data[cdc14Data==3] = "ON"
  cdc14Data = t(cdc14Data)
  cdc14Data = cdc14Data[nrow(cdc14Data):1,]
  M = which(as.character(cdc14Data[6,])=="ON")
  SA = which(as.character(cdc14Data[5,])=="ON")
  order = c(setdiff(M,SA),setdiff(1:ncol(cdc14Data), union(M,SA)), setdiff(SA,M), intersect(M,SA))
  cdc14Data = cdc14Data[,order]
  return(cdc14Data)
}


