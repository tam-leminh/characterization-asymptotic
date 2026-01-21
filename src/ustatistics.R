up_l2f2 <- function(suby) {
  nrow <- nrow(suby)
  ncol <- ncol(suby)
  yty <- t(suby)%*%suby
  return((sum(yty) - sum(diag(yty)))/(nrow*ncol*(ncol-1)))
}

up_l2 <- function(suby) {
  nrow <- nrow(suby)
  ncol <- ncol(suby)
  suby2 <- suby^2
  yty <- t(suby)%*%suby
  yyt <- suby%*%t(suby)
  return((sum(suby)^2 - sum(yty) + sum(diag(yty)) - sum(yyt) + sum(diag(yyt)) - sum(suby2))/(nrow*(nrow-1)*ncol*(ncol-1)))
}
