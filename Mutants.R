MutantBehaviour = function(X){
  maxT = 10000
  N = 100
  exits = c()
  tos = c()
  for(j in 1:nrow(X)){
    tICs = eval(parse(text=paste0("ICs",X$Model[j])))
    if(is.na(X$`Protein 2`[j])){# single mutation
      if(is.na(X$`Protein 1`[j])){
        tprot = "NA"
        tmut = "NA"
        tmodel = eval(parse(text=paste0("model",X$Model[j])))
        tunreg = paste0("model",X$Model[j],"unreg.csv")
      } else if(X$`Mutation 1`[j] == "PM"){
        if(X$`Protein 1`[j] == "Cdc15"){
          
          tprot = "NA"
          tmut = "NA"
          tmodel = eval(parse(text=paste0("model",X$Model[j],"Cdc157A")))
          tunreg = paste0("model",X$Model[j],"Cdc157Aunreg.csv")
        } else if(X$`Protein 1`[j] == "Mob1"){
          tprot = "NA"
          tmut = "NA"
          tmodel = eval(parse(text=paste0("model",X$Model[j],"Mob12A")))
          tunreg = paste0("model",X$Model[j],"Mob12Aunreg.csv")
        }
      } else{
        tprot = X$`Protein 1`[j]
        tmut = X$`Mutation 1`[j]
        if(X$`Mutation 1`[j] %in% c("SPB","mSPB","dSPB")){
          tmodel = eval(parse(text=paste0("model",X$Model[j],"forceLoc")))
          tICs = eval(parse(text=paste0("ICs",X$Model[j],"FL")))
        }
        else {
          tmodel = eval(parse(text=paste0("model",X$Model[j])))
        }
        tunreg = paste0("model",X$Model[j],"unreg.csv")
      }
    } else { # double mutant
      if(X$`Mutation 1`[j] == "PM"){
        if(X$`Mutation 2`[j] == "PM"){
          tprot = "NA"
          tmut = "NA"
          tmodel = eval(parse(text=paste0("model",X$Model[j],"DPM")))
          tunreg = paste0("model",X$Model[j],"DPMunreg.csv")
        } else if(X$`Protein 1`[j] == "Cdc15"){
          tprot = X$`Protein 2`[j]
          tmut = X$`Mutation 2`[j]
          tmodel = eval(parse(text=paste0("model",X$Model[j],"Cdc157A")))
          tunreg = paste0("model",X$Model[j],"Cdc157Aunreg.csv")
        } else if(X$`Protein 1`[j] == "Mob1"){
          tprot = X$`Protein 2`[j]
          tmut = X$`Mutation 2`[j]
          tmodel = eval(parse(text=paste0("model",X$Model[j],"Mob12A")))
          tunreg = paste0("model",X$Model[j],"Mob12Aunreg.csv")
        }
      } else{
        tprot = c(X$`Protein 1`[j],X$`Protein 2`[j])
        tmut = c(X$`Mutation 1`[j],X$`Mutation 2`[j])
        if(X$`Mutation 1`[j] %in% c("SPB","mSPB","dSPB")){
          tmodel = eval(parse(text=paste0("model",X$Model[j],"forceLoc")))
          tICs = eval(parse(text=paste0("ICs",X$Model[j],"FL")))
        }
        else {
          tmodel = eval(parse(text=paste0("model",X$Model[j])))
        }
        tunreg = paste0("model",X$Model[j],"unreg.csv")
      }
      
    }
    for(CC in c("M", "EA", "LA", "NP")){
      temp1 = asyncStats(prot = tprot, mutation = tmut,model = tmodel,unregfile = tunreg, ICs = tICs, CCstage = CC, timeMax = maxT, N=N)
      exits = c(exits,sum(as.numeric(temp1$exits=="ME")))
      tos = c(tos,sum(as.numeric(temp1$exits=="time")))
    }
  }
  exits2 = matrix(exits, ncol = 4, byrow = T)
  tos2 = matrix(tos, ncol = 4, byrow = T)
  colnames(exits2) = c("M_E", "EA_E", "LA_E", "NP_E")
  colnames(tos2) = c("M_T", "EA_T", "LA_T", "NP_T")
  X2 = cbind(X, exits2) %>% cbind(tos2)
  return(X2)
}

library(readxl)
inputs = read_xlsx("SimulationInputs.xlsx")

output = MutantBehaviour(inputs)

write_csv(output,"SimulationOutputs.csv")

inputs2 = read_xlsx("SimulationInputs2.xlsx")

output2 = MutantBehaviour(inputs2)

write_csv(output2,"SimulationOutputs2.csv")

createNetwork("Activity4a.txt", "Localization4b.txt", "FEARnet0.txt", "LocSpecificActivity4a.txt", "model4b", directory = "networkData")
model4b = loadNetwork("model4bOE.txt")
ICs4b = makeICs("MetaAIC3.csv", "MetaLIC3.csv", model4b, directory = "networkData")


inputs3 = read_xlsx("SimulationInputs3.xlsx")

output3 = MutantBehaviour(inputs3)
write_csv(output3,"SimulationOutputs3.csv")
