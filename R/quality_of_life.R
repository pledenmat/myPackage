#' Change working directory
#'
#' @description This function will change the working directory to the specified folder. It
#' assumes that the destination folder is directly reachable from either the
#' current working directory or one of its parent directories
#'
#' @param folder Destination folder
#'
#' @return
#' @export
#'
go_to <- function(folder){
  mission_accomplished <- grepl(folder,getwd())
  while (!mission_accomplished) {
    if (folder %in% list.files()) {
      setwd(folder)
      mission_accomplished <- T
    }else{
      if (getwd() == "C:/") {
        print("//!\\ Could not find target directory //!\\")
        break
      }
      setwd("..")
    }
  }
}
