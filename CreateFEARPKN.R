# This script takes a handwritten PKN and turns it into the form required by CellNOptR

library(tidyverse)
library(igraph)
library(readxl)
library(stringr)

ints1 = read_excel("FEARPKN.xlsx", sheet = 1)

verts = unique(c(ints1$`Interactor 1`,ints1$`Interactor 2`))

posints = ints1 %>% filter(Sign == 1) %>% select(`Interactor 1`,`Interactor 2`)
negints = ints1 %>% filter(Sign == -1) %>% select(`Interactor 1`,`Interactor 2`)
allints = posints %>% add_row(`Interactor 1` = negints$`Interactor 1`,`Interactor 2` = negints$`Interactor 2`)
net = graph_from_data_frame(allints,vertices = verts,directed = TRUE)

write_graph(net, "FEARPKN1.gml", format = "gml") # write to .gml format if required
degree = igraph::degree
degs = degree(net, mode = "in") # find in-degrees of nodes
inputs = names(degs[degs==0]) # find nodes that have in-degree 0

# Make color vector
lengths = c(nrow(posints),nrow(negints))
edgecols = rep(c("green","red"),lengths)

# plot network to view
l = layout_as_tree(net, root = which(verts%in%inputs), circular = F) 
l = norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
plot(net, rescale=F, layout=l*1, vertex.label.dist = 2,vertex.size = 7,edge.arrow.size = 0.5, edge.curved=.2,vertex.color = "black", edge.color = edgecols)

# runs some sanity checks 
red_ints = select(ints1, 'Interactor 1', 'Interactor 2') #As no sign, this could flag duplciations or contradictions!
cat(c("The following lines in the PKN file are duplicates:",which(duplicated(red_ints))+1))
cat(c("The nodes of the FEAR PKN are:",sort(verts)))
cat(c("The roots of the network are:", sort(inputs)))

ints2 = select(ints1, `Interactor 1`, Sign, `Interactor 2`)
PKN = tibble(n1 = ints2$`Interactor 1`, sign = ints2$Sign, n2 = ints2$`Interactor 2`)
nodes = unique(c(PKN$n1,PKN$n2))

# Add all OE nodes
protVerts = setdiff(nodes,c("ME", "CDK_high", "CDK_low")) # treat CDK seperately
PKN2 = rbind(PKN, tibble(n1 = paste0(protVerts,"_OE"), sign = 1, n2 = protVerts))

# Need a CDK OE node so that both high and low CDK activated by OE
PKN2 = rbind(PKN2, tibble(n1 = "CDK_OE", sign = 1, n2 = "CDK_high")) # CDK OEs
PKN2 = rbind(PKN2, tibble(n1 = "CDK_OE", sign = 1, n2 = "CDK_low")) # CDK OEs

# Cdc5 has it own OE node which has same out-edges as Cdc5 in PKN
PKN2 = rbind(PKN2, tibble(n1 = "Cdc5_OE", sign = -1, n2 = "Net1")) # Cdc5 feed forward 

# Write as .tsv
write_tsv(PKN2, "FEARPKN.txt",col_names=FALSE)



