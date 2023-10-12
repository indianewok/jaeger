#' @export
barcodes_to_bits<-function(barcode){
  stringi::stri_replace_all_fixed(str = barcode, pattern = c("A","C","T","G"), replacement = c("00","01","10","11"),
                         vectorize_all = FALSE) %>% bit64::as.integer64.bitstring(.)
}
