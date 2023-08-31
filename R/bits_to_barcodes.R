bits_to_barcodes<-function(int64_code, barcode_length = NULL){
  index<-c((32-barcode_length+1):(32))
  as.bitstring(int64_code) %>%
    strsplit(., "(?<=\\G.{2})", perl = TRUE) %>%
    sapply(., FUN = function(x){
      unlist(x) %>%
        .[index] %>%
        str_replace_all(string = ., fixed(c("00"="A","01"="C","10"="T","11"="G"))) %>%
        reduce(.,paste0)},
      simplify = TRUE)
}
