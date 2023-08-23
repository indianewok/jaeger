poly_scanner<-function(query, sequences){
  poly_found<-vector("logical", length(sequences))
  locations<-vector("list", length(sequences))

  poly_found<-stri_detect_regex(sequences, query)
  locations<-stri_locate_all_regex(sequences, query)

  paired_results<-mapply(function(found,loc){list(poly_found = found, locations = loc)}, poly_found, locations, SIMPLIFY = FALSE)
  return(paired_results)
}
