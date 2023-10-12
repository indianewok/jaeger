plucky<-function(alignments, ...){
  lapply(alignments, FUN = function(x){
    purrr::pluck(x, ...) %>% unlist(.)
  })
}
