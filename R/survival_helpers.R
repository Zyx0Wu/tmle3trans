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

#' Hazard Matrix to Survival Matrix
#' 
#' @param hm hazard matrix
#' 
#' @export
hm_to_sm = function(hm){
  # TODO: check
  sm <- t(apply(1-hm,1,cumprod))
  return(sm)
}

#' Survival Data Frame Time Variable Expansion
#' 
#' @param df data frame
#' @param t_current current time
#' 
#' @export
df_time = function(df, t_current, t_tilde="T.tilde", delta="Delta"){
  df_t <- copy(df)
  # TODO: check
  df_t$t <- t_current
  df_t$Failed <- as.numeric(t_current == df[[t_tilde]] & df[[delta]] == 1)
  df_t$Censored <- as.numeric(t_current == df[[t_tilde]] & df[[delta]] == 0)
  df_t$pre_failure <- as.numeric(t_current <= df[[t_tilde]])
  
  return(df_t)
}

#' Survival Data All Time Expansion
#' 
#' @param data data
#' @param node_list node_list
#' 
#' @export
transform_data = function(data, node_list) {
  T_tilde_name <- node_list[["T.tilde"]]
  Delta_name <- node_list[["Delta"]]
  T_tilde_data <- data[T_tilde_name]
  Delta_data <- data[Delta_name]
  k_grid <- 1:max(T_tilde_data)
  
  if (is.null(node_list$id)) {
    id <- 1:nrow(data)
    data <- cbind(ID=id, data)
    node_list$id <- "ID"
  }
  
  all_times <- lapply(k_grid, function(t_current) df_time(data, t_current, T_tilde_name, Delta_name))
  df_long <- rbindlist(all_times)
  
  long_node_list <- copy(node_list)
  long_node_list$time <- "t"
  long_node_list$pre_failure <- "pre_failure"
  long_node_list$failed <- "Failed"
  long_node_list$censored <- "Censored"
  
  return(list(long_data=df_long, long_node_list=long_node_list))
}

