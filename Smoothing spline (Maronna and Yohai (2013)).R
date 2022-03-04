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
  beta.hat <- b.sp.e%*%proj.bspline%*%alpha.c/p
  resids <- y - fit.mm$estimates[[1]]$intercept - Z%*%fit.mm$estimates[[1]]$beta
  return(list(bh = beta.hat, resids = resids, scale = fit.mm$scale))
}

bsb <- create.bspline.basis(rangeval = c(0,1), breaks = c(0,unique(x),1),
                            norder = 2*m )
bsbe <-  eval.basis(bsb, x)
bsbe.d <- eval.basis(bsb, x, Lfdobj=2 )
bsbe.dd <- rbind(bsbe.d[1, ], bsbe.d[100, ])
# qr.d <- qr(bsbe.dd)
sv <- svd(bsbe.dd, nv = 104)
n.sp <- sv$v[, 3:104] #### Nullspace

bsbe.dd%*%n.sp

nat.basis <- bsbe%*%n.sp
# beta <- rnorm(200)k
Pen.matrix <- t(n.sp)%*%bsplinepen(bsb, Lfdobj = m)%*%n.sp

# beta.hat <- nat.basis%*%solve(t(nat.basis)%*%nat.basis + 1e-08*Pen.matrix, t(nat.basis)%*%y)
# plot(beta.hat, type = "l")
# abline(v = 0.005025126)
