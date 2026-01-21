estimate_f2_condesp_10 <- function(y, row) {
  nrow <- nrow(y)
  ncol <- ncol(y)
  yi <- y[row,]
  yni <- y[-row,]
  ynityni <- t(yni)%*%yni
  mu_i <- 1/2/ncol/(ncol-1)*(sum(yi%o%yi) - yi%*%yi + 1/(nrow-1)*(sum(ynityni) - sum(diag(ynityni))))
  return(mu_i)
}

estimate_f2_condesp_10_vec <- function(y) {
  mu <- sapply(1:nrow(y), function(x) {estimate_f2_condesp_10(y, x)})
  return(mu)
}

estimate_f2_condesp_01_vec <- function(y) {
  nrow <- nrow(y)
  ncol <- ncol(y)
  yty <- t(y)%*%y
  mu <- 1/nrow/(ncol-1)*(colSums(yty) - diag(yty))
  return(mu)
}

estimate_l2_condesp_10_vec <- function(y) {
  nrow <- nrow(y)
  ncol <- ncol(y)
  mu <- (rowSums(y)*sum(y) -
           rowSums(y%*%t(y)) -
           sapply(1:nrow(y), function(x){ sum(y[x,]%o%y[x,]) - y[x,]%*%y[x,] } )) /
    (ncol*(ncol-1)*(nrow-1))
  return(mu)
}

estimate_l2_condesp_01_vec <- function(y) {
  nrow <- nrow(y)
  ncol <- ncol(y)
  mu <- (colSums(y)*sum(y) -
           colSums(t(y)%*%y) -
           sapply(1:ncol(y), function(x){ sum(y[,x]%o%y[,x]) - y[,x]%*%y[,x] } ) )/
    (nrow*(nrow-1)*(ncol-1))
  return(mu)
}
