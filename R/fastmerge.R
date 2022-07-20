#' Merges two data frames
#'
#' @description Custom function to merge two data frames together, even if they do
#' not have the exact same column names
#'
#' @param d1 first data frame to merge
#' @param d2 second data frame to merge
#'
#' @return Merged d1 and d2
#' @export
fastmerge <- function(d1, d2) {
  d1.names <- names(d1)
  d2.names <- names(d2)

  # columns in d1 but not in d2
  d2.add <- setdiff(d1.names, d2.names)

  # columns in d2 but not in d1
  d1.add <- setdiff(d2.names, d1.names)

  # add blank columns to d2
  if(length(d2.add) > 0) {
    for(i in 1:length(d2.add)) {
      d2[d2.add[i]] <- NA
    }
  }

  # add blank columns to d1
  if(length(d1.add) > 0) {
    for(i in 1:length(d1.add)) {
      d1[d1.add[i]] <- NA
    }
  }

  return(rbind(d1, d2))
}
