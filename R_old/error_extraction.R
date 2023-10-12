#' @export
error_extraction<-function(alignments, adapters){
  error_frame<-as.data.table(matrix(NA, nrow = nrow(alignments), ncol = length(adapters)))
  setnames(error_frame, names(adapters))
  for(adapter_name in names(adapters)){
    if(adapter_name == "poly_a" || adapter_name == "poly_t"){
      error_frame[[adapter_name]]<-alignments[[adapter_name]] %>% lapply(function(x){unlist(x[["poly_found"]])})
    }else{
      error_frame[[adapter_name]]<-alignments[[adapter_name]] %>% lapply(function(x){unlist(x[["editDistance"]])})
    }
  }
  colnames(error_frame)<-paste0(colnames(error_frame),
                                "_dist")
  error_frame<-error_frame %>% unnest(.)
  return(error_frame)
}
