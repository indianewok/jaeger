read_stats<-function(fn){
  seqkit<-"/dartfs/rc/lab/A/AckermanLab/CMV/PostRot/cmv/miniconda3/envs/medaka/bin/seqkit"
  seqkit_stats_args<-c("stats",fn,"-T","-a")
  stats<-processx::run(seqkit, seqkit_stats_args, echo_cmd =FALSE, echo = FALSE, spinner = TRUE)
  stats<-strsplit(x = stats$stdout, split = "\n") %>%
    lapply(., function(x){strsplit(x, split = "\t")}) %>%
    unlist(.,recursive = FALSE) %>%
    data.frame(.[[2]], row.names = .[[1]]) %>%
    .[,c(1,3)]
  colnames(stats)<-NULL
  stats[,1]<-NULL
  return(stats)
}
