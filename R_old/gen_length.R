gen_length<-function(column){
  out<-sapply(column, FUN = function(x){length(x)})
  return(out)
}
