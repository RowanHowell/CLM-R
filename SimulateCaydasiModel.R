#install.packages("remotes")
#remotes::install_github("jpahle/CoRC")
#CoRC::getCopasi()

library(CoRC)
library(tidyverse)
library(ggplot2)
library(reshape2)

Caydasi2012 = loadSBML("BIOMD0000000702.xml")

noAlign = loadSBML("BIOMD0000000702.xml") # version with no alignment event
deleteEvent("SPOC deactivation upon correct spindle alignment", model = noAlign)
SS = runSS(model = noAlign) # Steady state pre-alignment
setSpecies(key = SS$species$key, initial_concentration = SS$species$concentration, model = Caydasi2012) # set ICs to be steady state 


Bfa1SPB = loadSBML("BIOMD0000000702.xml") # version modelling Bfa1-SPB
setGlobalQuantities(key = "koffB", initial_value = 1e-6, model = Bfa1SPB) # these parameters control rate of disassociation of Bfa1 from SPB
setGlobalQuantities(key = "koffB4", initial_value = 1e-6, model = Bfa1SPB)

Bfa1SPBnoAlign = loadSBML("BIOMD0000000702.xml") # version with no alignment event (this has slightly different steady state to WT)
setGlobalQuantities(key = "koffB", initial_value = 1e-6, model = Bfa1SPBnoAlign) # these parameters control rate of disassociation of Bfa1 from SPB
setGlobalQuantities(key = "koffB4", initial_value = 1e-6, model = Bfa1SPBnoAlign) # disassociation rate decreased by 10^-3
deleteEvent("SPOC deactivation upon correct spindle alignment", model = Bfa1SPBnoAlign)
SSBfa1SPB = runSS(model = Bfa1SPBnoAlign) # Steady state pre-alignment
setSpecies(key = SSBfa1SPB$species$key, initial_concentration = SSBfa1SPB$species$concentration, model = Bfa1SPB) # set ICs to be steady state 


Bfa1KO = loadSBML("BIOMD0000000702.xml") # version modelling Bfa1 delete
setSpecies(key = SS$species$key, initial_concentration = SS$species$concentration, model = Bfa1KO) # set ICs to be steady state 

# Cycle through all species containing Bfa1
binding = c("", "B-")
phosphorylation = c("","P4","P5")
Tem1binding = c("","-Tem1GTP", "-Tem1GDP")
Tem1GTPtotal = 0
Tem1GDPtotal = 0
for(b in binding){
  for(p in phosphorylation){
    for(t in Tem1binding){
      if(t == "-Tem1GTP"){ # Note we do not wish to lose Tem1 stored in Bfa1-bound complexes in the original models IC
        Tem1GTPtotal = Tem1GTPtotal+getSpecies(paste0(b,"Bfa1",p,t))$initial_number
      } else if(t == "-Tem1GDP"){
        Tem1GDPtotal = Tem1GDPtotal+getSpecies(paste0(b,"Bfa1",p,t))$initial_number
      }
      setSpecies(paste0(b,"Bfa1",p,t), initial_number = 0, model = Bfa1KO)
    }
  }
}
setSpecies("Tem1GTP", initial_number = (getSpecies("Tem1GTP")$initial_number+Tem1GTPtotal), model = Bfa1KO) # add back Tem1 from complexes
setSpecies("Tem1GDP", initial_number = (getSpecies("Tem1GDP")$initial_number+Tem1GDPtotal), model = Bfa1KO)

tcWT = runTimeCourse(duration = 3600, intervals = 1000, model = Caydasi2012) # run time course and extract info
DataWT = tcWT$result_number%>%mutate( `Tem1GTP at the SPB` = `T-Tem1GTP`+`B-Bfa1-Tem1GTP`+`B-Bfa1P4-Tem1GTP`+`B-Bfa1P5-Tem1GTP`) %>% select(Time,`Tem1GTP at the SPB`)
DataWTBfa1 = tcWT$result_number%>%mutate( `Bfa1 at the SPB` = `B-Bfa1-Tem1GTP`+`B-Bfa1P4-Tem1GTP`+`B-Bfa1P5-Tem1GTP`+`B-Bfa1-Tem1GDP`+`B-Bfa1P4-Tem1GDP`+`B-Bfa1P5-Tem1GDP`+`B-Bfa1`+`B-Bfa1P4`+`B-Bfa1P5`) %>% select(Time,`Bfa1 at the SPB`)

tcBfa1SPB = runTimeCourse(duration = 3600, intervals = 1000, model = Bfa1SPB)
DataBfa1SPB = tcBfa1SPB$result_number%>%mutate( `Tem1GTP at the SPB` = `T-Tem1GTP`+`B-Bfa1-Tem1GTP`+`B-Bfa1P4-Tem1GTP`+`B-Bfa1P5-Tem1GTP`) %>% select(Time,`Tem1GTP at the SPB`)
DataBfa1SPBBfa1 = tcBfa1SPB$result_number%>%mutate( `Bfa1 at the SPB` = `B-Bfa1-Tem1GTP`+`B-Bfa1P4-Tem1GTP`+`B-Bfa1P5-Tem1GTP`+`B-Bfa1-Tem1GDP`+`B-Bfa1P4-Tem1GDP`+`B-Bfa1P5-Tem1GDP`+`B-Bfa1`+`B-Bfa1P4`+`B-Bfa1P5`) %>% select(Time,`Bfa1 at the SPB`)

tcBfa1KO = runTimeCourse(duration = 3600, intervals = 1000, model = Bfa1KO)
DataBfa1KO = tcBfa1KO$result_number%>%mutate( `Tem1GTP at the SPB` = `T-Tem1GTP`+`B-Bfa1-Tem1GTP`+`B-Bfa1P4-Tem1GTP`+`B-Bfa1P5-Tem1GTP`) %>% select(Time,`Tem1GTP at the SPB`)
DataBfa1KOBfa1 = tcBfa1KO$result_number%>%mutate( `Bfa1 at the SPB` = `B-Bfa1-Tem1GTP`+`B-Bfa1P4-Tem1GTP`+`B-Bfa1P5-Tem1GTP`+`B-Bfa1-Tem1GDP`+`B-Bfa1P4-Tem1GDP`+`B-Bfa1P5-Tem1GDP`+`B-Bfa1`+`B-Bfa1P4`+`B-Bfa1P5`) %>% select(Time,`Bfa1 at the SPB`)

# Plot how active Tem1 (Tem1-GTP levels at the SPB) changes over time
Tem1Data = tibble(Time = DataWT$Time, WT =DataWT$`Tem1GTP at the SPB`, Bfa1SPB =DataBfa1SPB$`Tem1GTP at the SPB`, Bfa1KO =DataBfa1KO$`Tem1GTP at the SPB`  )
mTem1Data = melt(Tem1Data, id = "Time")

ggplot(mTem1Data, aes(x = Time, y = value, color = variable)) + geom_line(lwd =2) + theme_bw() + geom_hline(yintercept = 65, lty = 2)+geom_vline(xintercept = 1800)+theme(panel.grid.minor = element_blank(),panel.border = element_blank(),axis.line = element_line(colour = "black"))
ggsave("Tem1-GTPatSPB.pdf")

# Plot how active Bfa1 loading at SPB changes over time
Bfa1Data = tibble(Time = DataWT$Time, WT =DataWTBfa1$`Bfa1 at the SPB`, Bfa1SPB =DataBfa1SPBBfa1$`Bfa1 at the SPB`, Bfa1KO =DataBfa1KOBfa1$`Bfa1 at the SPB`  )
mBfa1Data = melt(Bfa1Data, id = "Time")

ggplot(mBfa1Data, aes(x = Time, y = value, color = variable)) + geom_line(lwd =2) +geom_vline(xintercept = 1800)+ theme_bw()+theme(panel.grid.minor = element_blank(),panel.border = element_blank(),axis.line = element_line(colour = "black"))
ggsave("Bfa1atSPB.pdf")

