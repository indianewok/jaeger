#' @export
combo_framecore<-function(alignments, adapters, cores = NULL){
  if(is.null(cores)||cores == 1){
    combo_frame<-metric_framemode(mode = "both", alignments = alignments, adapters = adapters)
    return(combo_frame)
  }
  else{
    if(!is.null(cores)){
      combo_frame<-list("error","location")
      combo_frame<-mclapply(X = combo_frame, FUN = function(X){
        metric_framemode(mode = X, alignments = alignments, adapters = adapters)
      }, mc.cores = cores,
      mc.preschedule = TRUE,
      mc.cleanup = 9)
      combo_frame<-tidytable::bind_cols(combo_frame)
      return(combo_frame)
    }
  }
}
