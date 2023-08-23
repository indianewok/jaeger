aligner_target<-function(query, sequences,...){
  lapply(sequences, function(sequences){
    edlibalign(query = query, target = sequences, input_mode = "HW", input_task = "locations", k = -1, cigar_format = "extended")[c(1,3)]
  })
}
