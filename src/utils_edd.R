nrow_sequence <- function(N, c) {
  return(2+floor(c*(N+1)))
}

ncol_sequence <- function(N, c) {
  return(2+floor((1-c)*(N+1)))
}

simu_edd_with_N_sequence <- function(N, c, lambda, af=0, ag=0) {
  mN <- nrow_sequence(N, c)
  nN <- ncol_sequence(N, c)
  y <- simu_edd_alpha(mN, nN, lambda, af=af, ag=ag)
  return(y)
}

simu_overdisp_edd_with_N_sequence <- function(N, c, lambda, overdisp=0, af=0, ag=0) {
  mN <- nrow_sequence(N, c)
  nN <- ncol_sequence(N, c)
  y <- simu_overdisp_edd_alpha(mN, nN, lambda, mixvar=overdisp, af=af, ag=ag)
  return(y)
}

f2_to_alpha <- function(f2) {
  return(f2 - 1 + sqrt(f2*(f2-1)))
}

alpha_to_fk <- function(alpha,k) {
  return((alpha+1)^k/(k*alpha+1))
}