prepare_to_anger<-function(kit_prime, whitelist_path,...){
  adapt_metrics<-list()
  if(kit_prime == "three"){
    print("Three prime kit requested! Using the primers CTACACGACGCTCTTCCGATCT and CCCATGTACTCTGCGTTGATACCACTGCTT...")
    kit_3v2<-c("CTACACGACGCTCTTCCGATCT","T{12,}+","CCCATGTACTCTGCGTTGATACCACTGCTT",
               revcomp("CTACACGACGCTCTTCCGATCT"),"A{12,}+",revcomp("CCCATGTACTCTGCGTTGATACCACTGCTT"))
    names(kit_3v2)<-c("forw_primer","poly_t","rev_primer","rc_forw_primer","poly_a","rc_rev_primer")
    adapt_metrics$kit<-"three"
    adapt_metrics$forwards_order<-c("forw_primer","poly_t","rev_primer")
    adapt_metrics$reverses_order<-c("rc_forw_primer","poly_a","rc_rev_primer")
    adapt_metrics$forwards_dist<-c("forw_primer_dist", "rev_primer_dist")
    adapt_metrics$reverses_dist<-c("rc_forw_primer_dist", "rc_rev_primer_dist")
    adapt_metrics$forwards_pos<-c("forw_primer_pos","poly_t_pos","rev_primer_pos")
    adapt_metrics$reverses_pos<-c("rc_forw_primer_pos","poly_a_pos","rc_rev_primer_pos")
  }
  else{
    if(kit_prime == "five"){
      print("Five prime kit requested! Using the primers CTACACGACGCTCTTCCGATCT and GTACTCTGCGTTGATACCACTGCTT, and the TSO TTTCTTATATGGG...")
      kit_5v2<-c("CTACACGACGCTCTTCCGATCT","TTTCTTATATGGG","T{12,}+","GTACTCTGCGTTGATACCACTGCTT",
                 revcomp("CTACACGACGCTCTTCCGATCT"),revcomp("TTTCTTATATGGG"),"A{12,}+",revcomp("GTACTCTGCGTTGATACCACTGCTT"))
      names(kit_5v2)<-c("forw_primer","tso","poly_t","rev_primer","rc_forw_primer","rc_tso","poly_a","rc_rev_primer")
      adapt_metrics$kit<-"five"
      adapt_metrics$forwards_order<-c("forw_primer","tso","poly_a","rev_primer")
      adapt_metrics$reverses_order<-c("rc_forw_primer","rc_tso","poly_t","rc_rev_primer")
      adapt_metrics$forwards_dist<-c("forw_primer_dist", "tso_dist", "rev_primer_dist")
      adapt_metrics$reverses_dist<-c("rc_forw_primer_dist", "rc_tso_dist", "rc_rev_primer_dist")
      adapt_metrics$forwards_pos<-c("forw_primer_pos","tso_pos","poly_a_pos","rev_primer_pos")
      adapt_metrics$reverses_pos<-c("rc_forw_primer_pos","rc_tso_pos","poly_t_pos","rc_rev_primer_pos")
    }
  }
  if(!exists("whitelist")){
    print("Now setting up whitelist!")
    whitelist<-whitelist_importer(whitelist_path = whitelist_path, ...)
  }
  if(kit_prime == "three"){
    list2env(list(adapters = kit_3v2,
                  adapt_metrics = adapt_metrics,
                  whitelist = whitelist),
             envir = .GlobalEnv)
  }
  else{
    if(kit_prime == "five"){
      list2env(list(adapters = kit_5v2,
                    adapt_metrics = adapt_metrics,
                    whitelist = whitelist),
               envir = .GlobalEnv)
    }
  }
}
