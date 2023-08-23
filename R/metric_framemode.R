metric_framemode<-function(mode, alignments, adapters){
  if(mode == "location"){
    location_frame<-location_extraction(alignments, adapters)
    return(location_frame)
  } else{
    if(mode == "error"){
      error_frame<-error_extraction(alignments, adapters)
      return(error_frame)
    }
    else{
      if(mode == "both"){
        metric_frame<-vector(mode = "list", length = 2)
        metric_frame[[1]]<-location_extraction(alignments, adapters)
        metric_frame[[2]]<-error_extraction(alignments, adapters)
        metric_frame<-tidytable::bind_cols(metric_frame)
        return(metric_frame)
      }
    }
  }
}
