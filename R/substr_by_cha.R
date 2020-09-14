######## Auxiliary function of substring ##################
substr_by_cha <- function(cha_str, delim = "|", type = c("before", "after")){
  if (!is.character(cha_str)){
    stop("Input must be a character vector!")
  }
  cha_to_return <- character()
  if (length(cha_str)!=0){
    #cha_to_return <- character()
    for (i in 1:length(cha_str)){
      loc <- which(strsplit(cha_str[i], "")[[1]] == delim)
      if (type=="before"){
        cha_to_return[i] <- substr(cha_str[i], 1, loc - 1)
      } else if (type=="after"){
        cha_to_return[i] <- substr(cha_str[i], loc + 1, nchar(cha_str[i]))
      } else{
        stop("Substring type needs to be corrected!")
      }
    }
  }
  return(cha_to_return)
}