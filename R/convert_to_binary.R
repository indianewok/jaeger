convert_to_binary<-function(n){
  if(n > 1){
    return(paste0(convert_to_binary(as.integer64.double(n/2)), n %% 2))
  }
  return(as.character(n %% 2))
}
