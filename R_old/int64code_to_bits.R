int64code_to_bits<-function(int64_code){
  as.bitstring(int64_code) %>%
    strsplit(., "(?<=\\G.{2})", perl = TRUE) %>%
    sapply(., FUN = function(x){
      unlist(x)
    }) %>% t(.)
}
