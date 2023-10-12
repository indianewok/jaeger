whitelist_importer<-function(whitelist_path, save = NULL, convert = FALSE){
  print("Importing a little bit of sample data...")
  whitelist<-data.table::fread(input = whitelist_path, header = FALSE, data.table = TRUE, nrows = 10, nThread = 1,
                               col.names = "whitelist_bcs")
  if(class(whitelist$whitelist_bcs) == "character"){
    print("Character type detected! Importing the whitelist.")
    whitelist<-data.table::fread(input = whitelist_path, header = FALSE, data.table = TRUE, nThread = 1, col.names = "whitelist_bcs")
    if(convert != FALSE){
      print("Converting the whitelist into a more efficient format!")
      whitelist$whitelist_bcs<-barcodes_to_bits(whitelist$whitelist_bcs)
    }
    whitelist<-setkey(whitelist, "whitelist_bcs")
    gc()
    if(!is.null(save)){
      if(save == TRUE){
        new_path<-paste0(dirname(whitelist_path),"/",file_path_sans_ext(whitelist_path)
                         %>% basename(.),
                         "_bitlist.csv.gz")
        print(paste0("Now saving the integer whitelist to ",new_path,"!"))
        data.table::fwrite(x = whitelist, file = new_path, compress = "auto", nThread = 1, showProgress = TRUE, col.names = FALSE)
        gc()
      }
    }
    return(whitelist)
  }
  else{
    if(class(whitelist$whitelist_bcs) == "integer"){
      print("Integer64-type whitelist detected! Now importing...")
      whitelist<-data.table::fread(input = whitelist_path, header = FALSE, data.table = TRUE, nThread = 1,
                                   col.names = "whitelist_bcs", key = "whitelist_bcs", integer64 = "integer64")
      gc()
      return(whitelist)
    }
  }
}
