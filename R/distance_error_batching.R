distance_error_batching<-function(df, null_distances){
  df$batch<-NA
  df$batch<-tidytable::case_when(
    df$rc_forw_primer_dist >= null_distances["lower_bounds","rc_forw_primer_dist"] &
      df$forw_primer_dist <= null_distances["lower_bounds","forw_primer_dist"] ~ "forwards",
    df$forw_primer_dist >= null_distances["lower_bounds","forw_primer_dist"] &
      df$rc_forw_primer_dist <= null_distances["lower_bounds","forw_primer_dist"] ~ "reverses",
    df$forw_primer_dist <= null_distances["lower_bounds","forw_primer_dist"] &
      df$rc_forw_primer_dist <= null_distances["lower_bounds","rc_forw_primer_dist"] &
      df$rev_primer_dist <= null_distances["lower_bounds","rev_primer_dist"] &
      df$rc_rev_primer_dist <= null_distances["lower_bounds", "rc_rev_primer_dist"] ~ "rescuable_concatenates",
    TRUE ~ "undecided"
  )
  df[batch == "forwards" & gen_length(forw_primer_pos) > 2
     & !dup_finder(forw_primer_pos), batch := "parallel_f_concatenates"]
  df[batch == "reverses" & gen_length(rc_forw_primer_pos) > 2
     & !dup_finder(rc_forw_primer_pos), batch := "parallel_r_concatenates"]
  #df<-df %>% unite(., col = "error_code", c(grep(pattern = "_dist", x = colnames(df), value = TRUE)), remove = TRUE)
  return(df)
}
