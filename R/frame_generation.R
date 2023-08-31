#' @export
frame_generation<-function(df, adapters, cores){
  combo_frame<-df_dist(adapters = adapters, cores = cores, df = df) %>% combo_framecore(alignments = ., adapters = adapters, cores = cores)
  df<-tidytable::bind_cols(list(combo_frame, df))
  df<-adapter_error_calcs(df = df, adapt_metrics = adapt_metrics)
  return(df)
}
