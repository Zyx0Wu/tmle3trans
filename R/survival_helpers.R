#' Survival Data Vector to Matrix
#' 
#' @param x data vector
#' @param id patient id
#' @param time observation time
#' 
#' @export
long_to_mat = function(x, id, time){
  dt <- data.table(id=id,time=time,x=as.vector(x))
  wide <- dcast(dt, id~time, value.var="x")
  mat <- as.matrix(wide[,-1,with=FALSE])
  return(mat)
}

#' Density Matrix to Survival Matrix
#' 
#' @param dm density matrix
#' 
#' @export
dm_to_sm = function(dm){
  # TODO: check
  sm <- t(apply(1-dm,1,cumprod))
  sm <- cbind(1,sm[,-ncol(sm)])
  return(sm)
}

#' Survival Data Frame Time Variable Expansion
#' 
#' @param df data frame
#' @param t_current current time
#' 
#' @export
df_time = function(df, t_current){
  df_t <- copy(df)
  # TODO: check
  df_t$t <- t_current
  df_t$Failed <- as.numeric(t_current == df$T.tilde & df$Delta == 1)
  df_t$Censored <- as.numeric(t_current == df$T.tilde & df$Delta == 0)
  df_t$pre_failure <- as.numeric(t_current <= df$T.tilde)
  
  return(df_t)
}

