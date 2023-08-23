edlibalign<-function(query, target, input_mode, input_task, k, cigar_format, input_additionalEqualities = NULL){
  .Call('_edlibR_edlibalign', PACKAGE = 'edlibR', query, target, input_mode, input_task, k, cigar_format, input_additionalEqualities)
}
