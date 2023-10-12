#' @export
revcomp<-function(x){
  base::chartr(old = "ACTG", new = "TGAC", x = x) %>% stringi::stri_reverse(.)
}
