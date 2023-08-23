dist_consolidation<-function(df, distances_one, distances_two){
  if(nrow(df) > 10000){
    df<-df[1:10000,]
  }
  total_bounds<-cbind(
    distance_error(df, distances_two, distances_one),
    distance_error(df, distances_one, distances_two)
  )
  upper_bounds<-ceiling(total_bounds[1,]+total_bounds[2,])
  lower_bounds<-floor(total_bounds[1,]-total_bounds[2,])
  total_bounds<-rbind(upper_bounds,average = total_bounds[1,], lower_bounds)
  return(total_bounds)
}
