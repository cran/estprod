#' @import stats
ss_op <- function(y, free, capital, proxy, controls, id, time, phi, myprobit_lag, degree, free_length, k_length, proxy_length, controls_length, beta_fs_matrix, beta_fs_free, beta_fs_controls, maxiter, ...){
  
  #-----------------------------------------------------
  #Segundo estágio (ss)
  #-----------------------------------------------------

  capital_names <- colnames(capital)
  
  
  capital_lag = panel_lag(capital, id, time, verify = FALSE)
  
  names_capital_lag <- paste0("lag_", capital_names)
  
  colnames(capital_lag) <- names_capital_lag
  
  
  df_ss <- data.frame(y_ss = c(y - free %*% beta_fs_free - controls %*% beta_fs_controls - beta_fs_matrix[1]),
                      free,
                      capital,
                      capital_lag,
                      phi = c(phi),
                      phi_lag = c(panel_lag(phi, id, time, verify = FALSE)),
                      controls)
  
  df_ss <- df_ss[complete.cases(df_ss),]

  #objective function = sum of residuals ^2
  objective <- function(start, data_ss = df_ss) {
    
    data_ss <- as.matrix(data_ss)
    
    #constante: primeira variável. beta_k: (1 + 1) até o numero de k + 1  
    #beta_poly: as demais
    beta_const <- start[1]
    beta_k <- start[2:(k_length + 1)]
    beta_poly <- start[-(1:(k_length + 1))]
    
    # if (is.null(myprobit_lag)) {
    #   
    #   data_ss_qr <- cbind(const = 1 * beta_const,
    #                       capital = data_ss[, capital_names] %*% as.matrix(beta_k),
    #                       omega = (data_ss[, "phi_lag"] - data_ss[, names_capital_lag] %*% as.matrix(beta_k)) %*% beta_poly[1],
    #                       omega_sq = ((data_ss[, "phi_lag"] - data_ss[, names_capital_lag] %*% as.matrix(beta_k))^2)  %*% beta_poly[2])
    #   
    #   
    # }
    # else {
    #   
    #   data_ss_qr <- cbind(const = 1 * beta_const,
    #                       capital = data_ss[, capital_names] %*% as.matrix(beta_k),
    #                       omega = (data_ss[, "phi_lag"] - data_ss[, names_capital_lag] %*% as.matrix(beta_k)) %*% as.matrix(beta_poly[1]),
    #                       omega_sq = (data_ss[, "phi_lag"] - data_ss[, names_capital_lag] %*% as.matrix(beta_k))^2 %*% as.matrix(beta_poly[2]),
    #                       omega_probit = ((data_ss[, "phi_lag"] - data_ss[, names_capital_lag] %*% as.matrix(beta_k)) * myprobit_lag) %*% as.matrix(beta_poly[3]),
    #                       myprobit_lag  %*% as.matrix(beta_poly[4]),
    #                       myprobit_lag^2 %*% as.matrix(beta_poly[5]))
    # }
    # 
    
    poly_df <- cbind(omega = data_ss[, "phi_lag"] - data_ss[, names_capital_lag] %*% as.matrix(beta_k), myprobit_lag)
    
    poly_df <- poly(poly_df, degree = degree[2], raw = F)
    
    poly_df <- poly_df %*% diag(beta_poly, nrow = length(beta_poly))
    

    op <- rowSums(poly_df)
    
    op <- (data_ss[, "y_ss"] - 1 * beta_const - data_ss[, capital_names] %*% as.matrix(beta_k) - op) 
    
    
    #minpack.lm
    return(op)
    
    #Optim
    #return(sum(op^2))
    
  }
  
  #n_start = n_capital + n_const + n_interactions
  with_exit <- ifelse(is.null(myprobit_lag), 0, 1)
  #n_start <- k_length + 1 + 2 + with_exit
  n_start <- k_length + 1 + poly_elements(k_length + with_exit, d = degree[2])
  
  start <- rep(0, n_start)
  
  #Optimization
  ss_reg <- minpack.lm::nls.lm(par = start, fn = objective, control = minpack.lm::nls.lm.control(maxiter = maxiter,  nprint = 0))
  
  #ss_reg <- optim(par = start, fn = objective, method = "BFGS")
  #ss_reg <- optimx::optimx(par = start, fn = objective, method = c('Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'nlm', 'nlminb', 'spg', 'ucminf', 'newuoa', 'bobyqa', 'nmkb', 'hjkb', 'Rcgmin', 'Rvmmin'), itnmax = 100, control = list(ndeps = 1e-6))
  #print(ss_reg)
  
  assign("ss_reg", ss_reg, envir = parent.frame())
  
}
