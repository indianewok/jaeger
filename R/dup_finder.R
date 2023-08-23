dup_finder<-function(sequence_positions){
  plucky(sequence_positions) %>%
    sapply(., function(x){duplicated(x) %>%
        any(.)}, simplify = TRUE)
}
