#' @export
adapter_error_calcs<-function(df,adapt_metrics){
  null_distances<-dist_consolidation(df = df,
                                     distances_one = adapt_metrics$forwards_dist,
                                     distances_two = adapt_metrics$reverses_dist)
  df<-distance_error_batching(df, null_distances)
  return(df)
}
