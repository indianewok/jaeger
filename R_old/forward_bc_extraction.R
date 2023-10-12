forward_bc_extraction<-function(sequence, position_list, barcode_length){
  start<-sapply(position_list, FUN = function(x){pluck(x,2)})+1
  stop<-start+barcode_length-1
  barcode<-substr(sequence, start = start, stop = stop)
  return(barcode)
}
