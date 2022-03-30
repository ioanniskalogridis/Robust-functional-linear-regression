require(fda)
require(robustbase)

m.sp <- function(x, y,  k = 4.685){
  
  x <- as.matrix(x)
  n <- length(y)
  p.cand <- seq( round(max(n^{1/5}, 4)), round(8 + 2*n^{1/5}) )
  
  rob.ctrl <- lmrob.control(trace.level = 0,
                            nResample   = 5000,
                            tuning.psi  = k,
                            subsampling = 'simple',
                            rel.tol     = 1e-07,
                            refine.tol  = 1e-07,
                            k.max       = 5e3,
                            maxit.scale = 5e3,
                            max.it      = 5e3)
  BIC.values <- rep(NA, length(p.cand))
  for(j in p.cand){
    b.sp = create.bspline.basis(c(0, 1), j, norder = 4)
    b.sp.e <- eval.basis(seq(1/dim(x)[2], 1-1/dim(x)[2], len = dim(x)[2]), b.sp)
    x.p <- x%*%b.sp.e/dim(x)[2]
    fit.j <- lmrob(y~x.p, method = "MM", control = rob.ctrl)
    scale.j <- fit.j$scale
    resids.j <- fit.j$residuals
    BIC.values[j-3] <- log(scale.j^2*sum(Mpsi(resids.j, psi = "bisquare", cc = k, deriv = -1) )) + log(n)*j/n
  }
  p.opt <- p.cand[which.min(BIC.values)]
  b.sp = create.bspline.basis(c(0, 1), p.opt, norder = 4)
  b.sp.e <- eval.basis(seq(1/dim(x)[2], 1-1/dim(x)[2], len = dim(x)[2]), b.sp)
  x.p <- x%*%b.sp.e/dim(x)[2]
  fit.f <- lmrob(y~x.p, method = "MM", control = rob.ctrl)
  scale.f <- fit.f$scale
  resids.f <- fit.f$residuals
  
  beta.hat <- b.sp.e%*%fit.f$coefficients[-1]/dim(x)[2]
  return(list(bh = beta.hat, beta = fit.f$coefficients[-1], alpha = fit.f$coefficients[1],  weights = fit.f$rweights, 
              resids = resids.f))
}
