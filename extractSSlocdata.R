
extractSSlocdata = function(data, prots, locs, model){
  if(length(data)>1){ # combine SSs of varying length
    K = length(data)
    data2 = data[[1]]
    for(k in 2:K){
      data2 = cbind(data2, data[[k]])
    }
    data = data2
  } else{
    data = data[[1]]
  }
  
  
  SSs = colnames(data) # determine number and length of SSs
  SSsnum = as.tibble(str_split_fixed(SSs, "\\.", n =2))
  colnames(SSsnum) = c("SS", "state")
  N = length(unique(SSsnum$SS))
  lengths = c()
  for(j in 1:N){
    lengths[j] = max(as.numeric(filter(SSsnum, SS == paste0("Attr", j))$state))
  }
  SSslength = tibble(SS = paste0("Attr",1:N), length = lengths)
  
  nodes =c() # wrangle node names
  for(p in prots){
    nodest = grep(p, model$genes, value = TRUE)
    for(n in nodest){
      if(str_split(n, "L_")[[1]][2] %in% locs){
        nodes = c(nodes, n)
      }
    }
  }
  lownodes = grep("low", nodes, value = TRUE) # declare various useful colleciton of nodes
  highnodes = grep("high", nodes, value = TRUE)
  multinodes = sub("low","", lownodes)
  singlenodes = setdiff(nodes,c(lownodes,highnodes))
  allnodes = c(setdiff(nodes,c(lownodes,highnodes)), multinodes)
  
  allData = as_tibble(t(data[c("SACA_Nucleus", "SpindleAlign",nodes, "ME"),])) #transverse to make easier to use dplyr methods
  
  if(length(lownodes)>length(highnodes)){ # in case eg low at SPB but not high
    extra = sub("low","high",setdiff(lownodes, sub("high", "low", highnodes)))
    for(j in 1:length(extra)){
      allData = mutate(allData, !!extra[j] := 0)
    }
  }
  
  if(length(lownodes)>0){ # combine low/high nodes where necessary
    for(j in 1:length(lownodes)){
      allData = mutate(allData, !!multinodes[j] := allData[[lownodes[j]]] + allData[[sub("low","high",lownodes[j])]])
    }
  }
  if(length(singlenodes)>0){# double single level nodes to match HIGH
    for(j in 1:length(singlenodes)){
      allData = mutate(allData, !!singlenodes[j] := 2*allData[[singlenodes[j]]])
    }
  }
  allData = mutate(allData, ME = ME*3, SpindleAlign = 3*SpindleAlign, SACA_Nucleus = 3*SACA_Nucleus) # triple input/output nodes to distinguish
  allData = select(allData,SACA_Nucleus, SpindleAlign, !!allnodes, ME) # select relevant variables
  
  summaryData = data.frame(matrix(, nrow = 0, ncol = ncol(allData)))
  colnames(summaryData) = colnames(allData)
  
  for(j in 1:N){ # iterate over SSs
    if(SSslength$length[j]>1){
      subDF = allData[which(SSsnum == paste0("Attr",j)),]
      cSums = colSums(subDF)/as.numeric(SSslength$length[j]) # take mean value over all of the states in a cyclic attractor
      summaryData[j, ] = cSums
    } else{
      summaryData[j,] = allData[which(SSsnum == paste0("Attr",j)),]
    }
  }
  
  summaryData[summaryData==0] = "OFF"
  summaryData[summaryData==1] = "LOW"
  summaryData[summaryData==2] = "HIGH"
  summaryData[summaryData==3] = "ON"
  summaryData[summaryData<1 & summaryData>0] = "LOW/OFF" # note this cannot handle a scenario where oscillations occur over OFF/LOW/HIGH
  summaryData[summaryData<2 & summaryData>1] = "HIGH/LOW"
  summaryData[summaryData<3 & summaryData>2] = "MEON/OFF"
  summaryData = distinct(summaryData) # only keeps unique variables
  summaryData = t(summaryData)
  summaryData = summaryData[nrow(summaryData):1,]
  M = which(as.character(summaryData["SACA_Nucleus",])=="ON")
  SA = which(as.character(summaryData["SpindleAlign",])=="ON")
  order = c(setdiff(M,SA),setdiff(1:ncol(summaryData), union(M,SA)), setdiff(SA,M), intersect(M,SA))
  summaryData = summaryData[,order]
  return(summaryData)
}

PlotSSData = function(summaryData, filename){
  summaryDataM = melt(summaryData)
  summaryDataM$Var2 = as.factor(summaryDataM$Var2)
  setwd("SSPlots")
  tiff(filename, height = 400, width = 2000)
  print(ggplot(summaryDataM, aes(x = Var2, y = Var1, fill = value)) + geom_tile(colour="grey",size=0.25) + coord_fixed() + theme_bw() +theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),text=element_text(size=16,  family="Helvetica"), legend.title = element_blank(),legend.key.size =unit(1.5, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_fill_manual(values = c("OFF" = "white","ON" = "black","LOW" = "#DF8646","HIGH" = "#C43012", "HIGH/LOW" = "#F28835", "LOW/OFF" = "#FEE8AF", "MEON/OFF" = "grey")) + xlab("") + ylab("") )
  dev.off()
  setwd("..")
  
}