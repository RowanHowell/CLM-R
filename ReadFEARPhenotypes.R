# This script takes a handwritten list of phenotypes used for training and turns it into a CNO list as required for CellNOptR

library(tidyverse)
library(readxl)
library(stringr)

PhenosRaw = read_excel("FEARpheno.xlsx", sheet = 1) 

#------------------------------------------------------------------------------------------------------------------------------
# Pre-processing

Phenos  = PhenosRaw

source("TrainingFuncs.R")
findContradictions(Phenos, 1:11, 12) # look for contradictions, for now dealing with these case by case

red_Phenos = Phenos[,1:12]
cat(c("The following lines in the PKN file are duplicates:",which(duplicated(red_Phenos))+1)) # look for duplications
Phenos2 = Phenos[setdiff(1:nrow(Phenos),which(duplicated(red_Phenos))),] # removes duplication

#------------------------------------------------------------------------------------------------------------------------------
# create CNO table

perts = tibble(Protein = Phenos$`Protein 1`, Perturbation = Phenos$`Perturbation 1`) # make vector of all perturbations
perts = add_row(perts, Protein = c(Phenos$`Protein 2`,Phenos$`Protein 3`,Phenos$`Protein 4`,Phenos$`Protein 5`), Perturbation = c(Phenos$`Perturbation 2`,Phenos$`Perturbation 3`,Phenos$`Perturbation 4`,Phenos$`Perturbation 5`))
perts = perts[!is.na(perts$Protein),]
perts = perts[!is.na(perts$Perturbation),]# get rid of unused parts

perts = unique(perts)
gofs = perts[perts$Perturbation%in%c("OE","hyperactive","Nucleus"),] %>% arrange(Protein) # modified for v7
lofs = perts[perts$Perturbation%in%c("KD","delete","deplete","Cytoplasm"),] %>% arrange(Protein)

PKN = read_delim("FEARPKN.txt", col_names = c("n1", "sign", "n2"), delim = "\t")

PKNnoOE = PKN[which(!str_sub(PKN$n1, -3 ,-1)=="_OE"),]
inputs = setdiff(unique(PKNnoOE$n1), unique(PKNnoOE$n2))

gofs2 = setdiff(gofs$Protein,inputs)
lofs2 = setdiff(lofs$Protein,inputs)


headings = c("TR:Cell:CellLine")
for(i in 1:length(lofs2)){ # headings for perturbations
  headings = c(headings, paste0("TR:",lofs2[i],"i")) # Note need SPB node inhbition too
}

for(i in 1:nrow(gofs)){ # Note gofs still need inputs as some inputs are OE
  headings = c(headings, paste0("TR:",gofs[i,1],"_OE")) # Note change in notation
}

for(i in 1:length(inputs)){
  headings = c(headings, paste0("TR:",inputs[i]))
}

headings = c(headings, "DV:Cdc14", "DA:ALL")

CNO = data.frame(t(numeric(length(headings)))) # Note needs a (t=0) value for some reason
CNO[1,which(headings %in% c("DV:Cdc14"))] =NA # Set t=0 values to NA (don't want to train with these)
colnames(CNO) = headings
N = length(headings)
CNO = as.tibble(CNO)
for(j in 1:nrow(Phenos)){
  newrow = numeric(length(headings)) # empty row
  newrow[which(headings == "TR:SAC")] = Phenos$SAC[j]
  newrow[which(headings == "DV:Cdc14")] = Phenos$FEAR[j]
  for(k in 1:5){
    if(!is.na(Phenos[j,2*k])){ # check its not empty
      if(Phenos[j,(2*k)+1]%in%c("KD","delete","deplete","Cytoplasm")){ # If inhibition, note v7 has descriptive terms
        if(Phenos[j,2*k]%in%inputs){
          newrow[which(headings==paste0("TR:",Phenos[j,2*k]))]=0 # switch off for inhibition
        }
        else{
          newrow[which(headings==paste0("TR:",Phenos[j,2*k],"i"))]=1 # If inhibition, activates inhibition node
        }
      }
      else{ 
        if(Phenos[j,2*k]%in%inputs){
          newrow[which(headings==paste0("TR:",Phenos[j,2*k]))]=1 # switch on for activation
        }
        newrow[which(headings==paste0("TR:",Phenos[j,2*k],"_OE"))]=1 # If activation, activates OE node directly too
      } 
    }
  }
  tempProts = unlist(Phenos[j,2*(1:5)])[which(!is.na(unlist(Phenos[j,2*(1:5)])))]
  unPertInputs = setdiff(inputs,c(tempProts, "SAC")) # SAC is dealt with seperately
  for(prot in unPertInputs){
    newrow[which(headings==paste0("TR:",prot))]=1 # default position is on
  }
  newrow[N] = 1
  CNO = rbind(CNO,newrow)
  
}
findContradictions(CNO,1:(N-2),N-1) # This will find contradictions even when the order is different

write_csv(CNO,"FEARCNO.csv")

OEedges = which(str_sub(PKN$n1,-3,-1)=="_OE")
noOEedges = OEedges[which(!PKN$n1[OEedges]%in%paste0(gofs$Protein,"_OE"))]
inputsOEedges = which(PKN$n2%in%inputs)
edgesKeep = setdiff(1:nrow(PKN), c(noOEedges,inputsOEedges))

PKNedit = PKN[edgesKeep,]

# modified version of PKN based on perturbations in CNOlist
write_tsv(PKNedit, "FEARPKNedit.txt", col_names = FALSE)
