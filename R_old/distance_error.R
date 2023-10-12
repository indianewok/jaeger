distance_error<-function(df, distances_one, distances_two){
  out<-df %>% select(all_of(distances_one)) %>%
    apply(X = ., MARGIN = 1, FUN=function(X){all(X == 0)}) %>%
    which(.) %>%
    df[.,] %>%
    select(., all_of(distances_two))
  out<-out %>% apply(X = ., MARGIN = 1, FUN=function(X){all(X > 0)}) %>%
    which(.) %>%
    out[.,] %>%
    apply(X = ., MARGIN = 2,FUN = function(X){c(mean(X),sd(X))})
  bounds<-c(floor(out[1,1]-out[2,1]), ceiling(out[1,1]+out[2,1]))
  return(out)
}
