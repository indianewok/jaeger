lord_fastquad<-function(df,fn,type,append){
  if(type == "fq"){
    fq_list<-df %>%
      dplyr::mutate(dummy = "+") %>%
      dplyr::select(id, seq, dummy, qual) %>%
      as.matrix() %>%
      t() %>%
      as.character() %>%
      as_tibble()
    tgutil::fwrite(fq_list, file = fn, compress = "auto", col.names = FALSE, quote = FALSE, append = append)
  }
  if(type == "fa"){
    if(any(grepl(">", df$id)) == FALSE){
      df$id<-paste0(">",df$id)
      fa_list<-df %>%
        dplyr::select(id, seq) %>%
        as.matrix() %>%
        t() %>%
        as.character() %>%
        as_tibble()
      tgutil::fwrite(fa_list, file = fn, compress = "auto", col.names = FALSE, quote = FALSE, append = append)
    }
    fa_list<-df %>%
      dplyr::select(id, seq) %>%
      as.matrix() %>%
      t() %>%
      as.character() %>%
      as_tibble()
    tgutil::fwrite(fa_list, file = fn, compress = "auto", col.names = FALSE, quote = FALSE, append = append)
  }
}
