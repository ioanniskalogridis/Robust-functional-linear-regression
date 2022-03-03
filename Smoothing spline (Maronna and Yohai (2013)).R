require(fda)
require(pense)

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
  
  fit.mm <- pensem_cv(x = Z, y = y, alpha = 0, cc = 4.685, nlambda = 20)
  alpha.c = c(p*solve(C.m)%*%fit.mm$estimates[[1]]$beta)
  beta.hat <- b.sp.e%*%proj.bspline%*%alpha.c
  
}