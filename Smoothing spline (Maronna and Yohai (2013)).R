require(fda)
require(pense)
require(rrcov)

m.sm.sp <- function(x, y, t, m = 2){
  n <- length(y)
  p <- length(t)
  t <- t/max(t)
  
  b.sp <- create.bspline.basis(c(0,1), breaks = t, norder = 2*m)
  b.sp.e <- eval.basis(t, b.sp)
  poly.m <- cbind(rep(1, p), poly(t, degree = m-1, raw = TRUE))
  proj.bspline <- ginv(t(b.sp.e)%*%b.sp.e)%*%t(b.sp.e)
  p.m <- p*t(proj.bspline)%*%bsplinepen(b.sp, Lfdobj=m)%*%proj.bspline + poly.m%*%solve( t(poly.m)%*%poly.m, t(poly.m) )
  
  C.m <- chol(p.m)
  Z <- x%*%solve(C.m)
  # pca.l <- PcaLocantore(Z)
  
  # p.l <- function(x){
  #   p.l <- sum( pca.l$eigenvalues/(pca.l$eigenvalues + 3*x/(4.685^2/n)) )
  #   return(p.l)
  # }
  # lambda.min <- uniroot(function(x) p.l(x)-0.5*n+1, lower = 1e-14, upper =  1e-02, tol = 1e-20)
  # lambda.min <- lambda.min$root
  # lambda.max <- uniroot(function(x) p.l(x)-1, lower = 1e-06, upper =  2, tol = 1e-20)
  # lambda.max <- lambda.max$root
  # 
  # lambda.i <- seq(p.l(lambda.min), p.l(lambda.max), length = 20)
  # lambda.m <- sapply(lambda.i, FUN = function(s) uniroot( function(x) p.l(x)- s, lower = 1e-14, upper = 1, tol = 1e-20 )$root )
  
  fit.mm <- pensem_cv(x = Z, y = y, alpha = 0, cc = 4.685, nlambda = 50)
  alpha.c = c(p*solve(C.m)%*%fit.mm$estimates[[1]]$beta)
  beta.hat <- b.sp.e%*%proj.bspline%*%alpha.c/p
  resids <- y - fit.mm$estimates[[1]]$intercept - Z%*%fit.mm$estimates[[1]]$beta
  return(list(bh = beta.hat, resids = resids, scale = fit.mm$scale))
}
