#'@export
read_fastqas<-function(fn,type, ...){
  ext<-tools::file_ext(fn)
  if (ext == "gz"){
    cat_cmd<-"zcat"
  }
  else{
    cat_cmd<-"cat"
  }
  if(type == "fq"){
    res<-fread(cmd = glue::glue("{cat_cmd} {fn} | paste - - - - | cut -f1,2,4"), col.names = c("id", "fastq_files","qc"), sep = "\t", ...)
    res<-res[,c(1,3,2)]
    res$id<-str_split(res$id, pattern = " ") %>% sapply(., function(x){x[1]}, USE.NAMES = FALSE)
    return(res)
  }
  if(type == "fa"){
    res<-fread(cmd = glue::glue("{cat_cmd} {fn} | paste - - | cut -f1,2"), col.names = c("id", "seq"), sep = "\t", ...)
    return(res)
  }
}
