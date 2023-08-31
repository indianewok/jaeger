location_extraction<-function(alignments, adapters){
  location_frame<-as.data.table(matrix(NA, nrow = nrow(alignments), ncol = length(adapters)))
  setnames(location_frame, names(adapters))
  for(adapter_name in names(adapters)){
    if(adapter_name == "poly_a"||adapter_name =="poly_t"){
      location_frame[[adapter_name]]<-alignments[[adapter_name]] %>%
        lapply(function(x){
          x[["locations"]] %>% unlist(.) %>% as.integer(.) %>% {kit::funique(.,fromLast=TRUE)+1}
        })
    }else{
      location_frame[[adapter_name]]<-alignments[[adapter_name]] %>%
        lapply(function(x){
          x[["locations"]] %>% unlist(.) %>% {kit::funique(.,fromLast=TRUE)+1}
        })
    }
  }
  colnames(location_frame)<-paste0(colnames(location_frame), "_pos")
  return(location_frame)
}
