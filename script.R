rm(list=objects())
source("src/utils_edd.R")
source("src/ustatistics.R")
source("src/hoeffding_for_f2_g2_l2.R")

################################################
# Power test helpers
################################################
reject_unilateral_test <- function(stat, level) {
  return(stat > qnorm(level))
}
level <- 0.95

################################################
# QQ helpers
################################################
Ksim <- 500
h <- 1000
Kmc <- 10000
t <- seq(1/(h+1), h/(h+1), 1/(h+1))
q <- replicate(Kmc, quantile(rnorm(Ksim, 0, 1), probs=t))
qq_ci <- apply(q, MARGIN=1, function(x) { quantile(x, probs=c(0.005,0.995)) })

plot_qq_with_ci <- function(data, N_val, qq_ci, t) {
  qqnorm(data, main=paste0("N=", N_val),
         xlim=c(-3,3), ylim=c(-3,3),
         cex=.7, pch=16,
         xlab="", ylab="")
  lines(qnorm(t, 0, 1), qq_ci[1,], col="red", lwd=2)
  lines(qnorm(t, 0, 1), qq_ci[2,], col="red", lwd=2)
  abline(a=0, b=1, col="blue")
  abline(a=0, b=0)
  abline(v=0)
}

################################################
# Example A
################################################
simulate_er <- function(N, c) {
  nr <- nrow_sequence(N, c)
  nc <- ncol_sequence(N, c)
  y <- matrix(rnorm(n=nr*nc, mean=0, sd=1), nr, nc)
  return(y)
}

simu_estimate_l2f2 <- function(N, c) {
  y <- simulate_er(N,c)
  theta <- up_l2f2(y)
  
  mu10_l2f2 <- estimate_f2_condesp_10_vec(y)
  mu01_l2f2 <- estimate_f2_condesp_01_vec(y)
  
  vhat <- 4/c*var(mu10_l2f2) + 4/(1-c)*var(mu01_l2f2)
  return(c("theta" = theta, "vhat" = vhat))
}

#### Variables
c <- 0.5

#### Simulate
#Ns <- floor(2^seq(3,8,1))
#z <- matrix(0, length(Ns), Ksim)
#var <- rep(0, length(Ns))
#ev <- rep(0, length(Ns))
#for (i in 1:length(Ns)) {
#  print(paste0("N=",Ns[i]))
#  t1 <- Sys.time()
#  res <- replicate(Ksim, simu_estimate_l2f2(Ns[i],c))
#  saveRDS(res, paste0("ERh1-res_", Ns[i], ".rds"))
#  z[i,] <- res[1,]
#  var[i] <- var(res[1,])
#  ev[i] <- mean(res[2,])
#  t2 <- Sys.time()
#  print(t2-t1)
#}

#### QQ plot under H_0
par(mfrow=c(2,3))
par(mar=c(3,3,3,2) + 0.1)
Ns <- floor(2^seq(3,8,1))

for (i in 1:length(Ns)) {
  res <- readRDS(paste0("data/degenerate_basic/ERh1-res_", Ns[i], ".rds"))
  z <- res[1,]
  plot_qq_with_ci(Ns[i]^(3/2)*z/sqrt(16), Ns[i], qq_ci, t)
}

################################################
# Example B
################################################
simu_estimate_test <- function(N, c, lambda, af, ag) {
  y <- simu_edd_with_N_sequence(N, c, lambda, af=af, ag=ag)
  l2f2_hat <- up_l2f2(y)
  l2_hat <- up_l2(y)
  theta <-  l2f2_hat - l2_hat
  vhat <- 2*l2_hat/(c*(1-c)^2)
  mu10 <- estimate_f2_condesp_10_vec(y) - estimate_l2_condesp_10_vec(y)
  mu01 <- estimate_f2_condesp_01_vec(y) - estimate_l2_condesp_01_vec(y)
  vnd <- 4/c*var(mu10) + 4/(1-c)*var(mu01)
  return(c("theta" = theta, "vhat" = vhat, "vnd" = vnd))
}

theoretical_variance_1 <- function(lambda, f2, g2, c) {
  alpha <- f2_to_alpha(f2)
  f3 <- alpha_to_fk(alpha, 3)
  f4 <- alpha_to_fk(alpha, 4)
  vth <- lambda^4/c*(f4 - f2^2 - 4*f3 + 8*f2 - 4)+4*lambda^4/(1-c)*(f2-1)^2*(g2-1)
  return(vth)
}

theoretical_variance_d <- function(lambda, c) {
  return(2*lambda^2/(c*(1-c)^2))
}

#### Variables
c <- 0.5
lambda <- 1
f2s <- seq(1,3,0.2)
vth <- theoretical_variance_d(lambda, c)

#### Simulate
#for (i in 1:length(Ns)) {
#  for (j in 1:length(f2s)) {
#    print(paste0("N=",Ns[i]))
#    print(paste0("F2=",f2s[j]))
#    t1 <- Sys.time()
#    res <- replicate(Ksim, simu_estimate_test(Ns[i],c, lambda, afs[j], ag))
#    saveRDS(res, paste0("F2power-res_", Ns[i], "_", j, ".rds"))
#    z[j,] <- res[1,]
#    vh[j,] <- res[2,]
#    var[j] <- var(res[1,])
#    ev[j] <- mean(res[3,])
#    t2 <- Sys.time()
#    print(t2-t1)
#  }
#}

#### Power test
par(mfrow=c(1,1))
par(mar=c(4,4,2,2) + 0.1)
Ns <- floor(2^seq(3,7,1))
cols <- rainbow(length(Ns)+1)

plot(1, type="n", xlab="F2", ylab="Power", xlim=c(1, 3), ylim=c(0, 1))
abline(h = c(0,1))
abline(v = 1)
abline(h = 1-level, lty = 5, col='grey')

for (i in 1:length(Ns)) {
  z <- matrix(0, length(f2s), Ksim)
  for (j in 1:length(f2s)) {
    res <- readRDS(paste0("data/degenerate_f2/F2power-res_", Ns[i], "_", j, ".rds"))
    z[j,] <- res[1,]
  }
  rates <- rowMeans(reject_unilateral_test(Ns[i]^(3/2)*z/sqrt(vth), level), na.rm = TRUE)
  lines(f2s, rates, col=cols[i],pch=16, lwd=1.5)
  points(f2s, rates, col=cols[i], pch=i)
}
legend("bottomright", legend=paste0("N=",Ns),
       col=cols, pch=1:length(Ns), lty=1:2, cex=0.8)

#### QQ plots under H_0/H_1
par(mfrow=c(2,3))
par(mar=c(3,3,3,2) + 0.1)
Ns <- floor(2^seq(3,8,1))

for (i in 1:length(Ns)) {
  res <- readRDS(paste0("data/degenerate_f2/F2power-res_", Ns[i], "_", 1, ".rds"))
  z <- res[1,]
  plot_qq_with_ci(Ns[i]^(3/2)*z/sqrt(2/c/(1-c)^2), Ns[i], qq_ci, t)
}

for (i in 1:length(Ns)) {
  res <- readRDS(paste0("data/degenerate_f2/F2power-res_", Ns[i], "_", 2, ".rds"))
  z <- res[1,]
  plot_qq_with_ci(Ns[i]^(3/2)*z/sqrt(2/c/(1-c)^2), Ns[i], qq_ci, t)
}

################################################
# Example C
################################################
kern1 <- function(suby22) {
  return(suby22[1,1]*suby22[1,2]*suby22[2,2])
}

simple_ukern2 <- function(y) {
  nrow <- nrow(y)
  ncol <- ncol(y)
  z <- y*(y-1)
  zty <- t(z)%*%y
  zyt <- z%*%t(y)
  lt <- (sum(z)*sum(y) - sum(zty) + sum(diag(zty)) - sum(zyt) + sum(diag(zyt)) - sum(z*y))/(nrow*(nrow-1)*ncol*(ncol-1))
  return(lt)
}

ukern <- function(y, kern) {
  nrow <- nrow(y)
  ncol <- ncol(y)
  row_pairs <- expand.grid(1:nrow, 1:nrow)
  row_pairs <- as.matrix(row_pairs[row_pairs$Var1 != row_pairs$Var2, ])
  col_pairs <- expand.grid(1:ncol, 1:ncol)
  col_pairs <- as.matrix(col_pairs[col_pairs$Var1 != col_pairs$Var2, ])
  nr2 <- nrow(row_pairs)
  nc2 <- nrow(col_pairs)
  u <- 0
  for (i in 1:nr2) {
    for (j in 1:nc2) {
      u <- u + kern(y[row_pairs[i,], col_pairs[j,]])
    }
  }
  return(u)
}

simu_estimate_test_c <- function(N, c, lambda, overdisp=0, af, ag) {
  y <- simu_overdisp_edd_with_N_sequence(N, c, lambda, overdisp, af=af, ag=ag)
  nrow <- nrow(y)
  ncol <- ncol(y)
  theta <- simple_ukern2(y) - ukern(y, kern1)/(nrow*(nrow-1)*ncol*(ncol-1))
  return(theta)
}

theoretical_variance_c <- function(lambda, f2, g2, c) {
  alphaf <- f2_to_alpha(f2)
  f3 <- alpha_to_fk(alphaf, 3)
  alphag <- f2_to_alpha(g2)
  g3 <- alpha_to_fk(alphag, 3)
  vth <- lambda^4/(c*(1-c))*(lambda*(f3-f2^2)*(g3-g2^2)+2*f2*g2)
  return(vth)
}

#### Variables
c <- 0.5
lambda <- 1
g2 <- 4/3
f2 <- 4/3
overdisps <- seq(0., 1, 0.1)
af <- 1
ag <- 1
vth <- theoretical_variance_c(lambda, f2, g2, c)

#### Simulate
#library(parallel)
#for (overdisp in overdisps) {
#  print(paste0("alpha=", overdisp))
#  
#  for (i in 1:length(Ns)) {
#    print(paste0("N=",Ns[i]))
#    t1 <- Sys.time()
#    res <- simplify2array(mclapply(1:Ksim, function(x) { simu_estimate_test_c(Ns[i],c, lambda, overdisp=overdisp, af, ag) }, mc.cores = 18))
#    saveRDS(res, paste0("data/degenerate_overdisp/overdisp-res_", Ns[i], "_alpha_0",  10*overdisp, ".rds"))
#    z <- res
#    var <- var(res)
#    t2 <- Sys.time()
#    print(t2-t1)
#  }
#  
#}

#### Power test
par(mfrow=c(1,1))
par(mar=c(4,4,2,2) + 0.1)
Ns <- floor(2^seq(3,7,1))
cols <- rainbow(length(Ns)+1)

plot(1, type="n", xlab="alpha", ylab="Power", xlim=c(0, 1), ylim=c(0, 1))
abline(h = c(0,1))
abline(v = 0)
abline(h = 1-level, lty = 5, col='grey')

z <- matrix(0, length(overdisps), Ksim)
for (i in 1:length(Ns)) {
  for (j in 1:length(overdisps)) {
    res <- readRDS(paste0("data/degenerate_overdisp/overdisp-res_", Ns[i], "_alpha_0", 10*overdisps[j], ".rds"))
    z[j,] <- res
  }
  rates <- rowMeans(reject_unilateral_test(Ns[i]*z/sqrt(vth), level), na.rm = TRUE)
  lines(overdisps, rates, col=cols[i],pch=16, lwd=1.5)
  points(overdisps, rates, col=cols[i], pch=i)
}
legend("bottomright", legend=paste0("N=",Ns),
       col=cols, pch=1:length(Ns), lty=1:2, cex=0.8)

#### QQ plots under H_0/H_1
par(mfrow=c(2,3))
par(mar=c(3,3,3,2) + 0.1)
Ns <- floor(2^seq(3,8,1))

for (i in 1:length(Ns)) {
  res <- readRDS(paste0("data/degenerate_overdisp/overdisp-res_", Ns[i], "_alpha_00.rds"))
  z <- res
  plot_qq_with_ci(Ns[i]*z/sqrt(vth), Ns[i], qq_ci, t)
}

for (i in 1:length(Ns)) {
  res <- readRDS(paste0("data/degenerate_overdisp/overdisp-res_", Ns[i], "_alpha_01.rds"))
  z <- res
  plot_qq_with_ci(Ns[i]*z/sqrt(vth), Ns[i], qq_ci, t)
}