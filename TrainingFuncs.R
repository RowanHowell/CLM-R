convertToBoolNet = function(nodes, edges, filename){
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

concatenate = function(vec, fill){
  N = length(vec)
  result = vec[1]
  for(i in 2:N){
    result = paste0(result, fill, vec[i])
  }
  return(result)
}

findContradictions = function(Data, ReadCols, OutCol, verbose = TRUE){
  # Data = tibble to search, ReadCols is columns to compare for similarity (as vector), OutCols is column to compare for contradiction
  # findContradictions(Phenos, 1:9, 10)
  # output nos take into account headings
  
  # Note this only finds contradictions when order is the same!
  output = tibble(Row1=character(), Row2=character())
  for(i in 2:nrow(Data)){
    flag = 0
    j = 1
    while(flag==0){
      if(identical(Data[i,ReadCols],Data[j,ReadCols])&(!identical(Data[i,OutCol],Data[j,OutCol]))){
        flag = 1
        output = add_row(output, Row1 = j+1, Row2 = i+1)
      }
      j = j+1
      if(j == i){
        flag = 1
      }
    }
  }
  if(verbose){cat("The following lines are contradictory:")}
  return(output)
  
}