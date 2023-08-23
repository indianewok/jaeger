df_dist<-function(adapters, df, cores){
  alignments<-as.data.table(matrix(NA, nrow = nrow(df), ncol = length(adapters)))
  setnames(alignments, names(adapters))
  results<-mcmapply(adapters, FUN = function(adapter){
    adapter_name<-which(adapters == adapter) %>% names(.)
    if(adapter == adapters["poly_a"] || adapter == adapters["poly_t"]){
      return(poly_scanner(query = adapter, sequences = df$fastq_files))
    }
    else{
      return(aligner_target(query = adapter, sequences = df$fastq_files))
    }
    #return(results)
  }, mc.cores = cores, mc.preschedule = TRUE, SIMPLIFY = FALSE, mc.cleanup = 9)
  for(i in seq_along(results)){
    alignments[[i]]<-results[[i]]
  }
  return(alignments)
}
