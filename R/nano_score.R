nano_score<-function(sequence){
  nsc<-DescTools::CharToAsc(sequence)-33
  return(nsc)
}
