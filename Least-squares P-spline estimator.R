require(fda)

ls.pen.sp <- function(x, y, norder = 4, nbasis = NULL, q  = 2, n.se = 0){
  
  x <- as.matrix(x)
  n <- length(y)
  nbasis = ifelse(is.null(nbasis),round(min(n/4, 40)), nbasis)
  b.sp <- create.bspline.basis(c(0, 1), nbasis = nbasis, norder = norder)
  b.sp.e <- eval.basis(seq(1/dim(x)[2], 1-1/dim(x)[2], len = dim(x)[2]), b.sp)
  # p.m <- t( diff(diag( (nbasis)) , differences = q ) )%*%diff(diag( (nbasis ) ), differences = q ) # + diag(nbasis)
  p.m <- bsplinepen(b.sp, Lfdobj = q )# + bsplinepen(b.sp, Lfdobj = 0 )
  
  x.p <- x%*%b.sp.e/dim(x)[2]

  p.m.a <- rbind( rep(0, (dim(x.p)[2]+1) ), cbind(rep(0, dim(x.p)[2]), p.m) )
  x.a <- cbind(rep(1, dim(x.p)[1]), x.p)
  
  GCV <- function(lambda){
    ls.est <- solve(t(x.a)%*%x.a + lambda*p.m.a, t(x.a)%*%y)
    resids <- y - x.a%*%ls.est
    hat.m <- x.a%*%solve(t(x.a)%*%x.a + lambda*p.m.a, t(x.a))
    GCV.scores <- mean( resids^2/(1-diag(hat.m))^2 )
    return(GCV.scores)
  }
  
  lambda.cand <- c(1e-09, 2e-08, 6e-08, 9e-08,  2e-07, 6e-07, 9e-07,  2e-06, 6e-06, 9e-06,
                   2e-05, 6e-05, 9e-05,  2e-04, 6e-04, 9e-04, 2e-03, 6e-03, 9e-03,  2e-02, 6e-02, 9e-02,
                   2e-01, 6e-01, 9e-01, 3, 50)
  
  lambda.e <- sapply(lambda.cand, FUN = GCV)
  wm <- which.min(lambda.e)
  if(wm==1){
    wm <- 2
  } else if(wm==length(lambda.cand)){wm <- (length(lambda.cand)-1)}
  lambda1 <- optimize(f = GCV, interval = c(lambda.cand[wm-1], lambda.cand[wm+1]) , tol = 1e-14)$minimum

  ls.est <- solve(t(x.a)%*%x.a + lambda1*p.m.a, t(x.a)%*%y)
  hat.m <- x.a%*%solve(t(x.a)%*%x.a + lambda1*p.m.a, t(x.a))
  resids <- y - x.a%*%ls.est
  
  if(n.se > 0){
    sdw <- sd(resids^2/(1-diag(hat.m))^2-GCV(lambda1))/sqrt(n)
    thr <- GCV(lambda1) + n.se*sdw
    lambda1 <- max( lambda.cand[lambda.e <= thr]  )
    
    ls.est <- solve(t(x.a)%*%x.a + lambda1*p.m.a, t(x.a)%*%y)
    hat.m <- x.a%*%solve(t(x.a)%*%x.a + lambda1*p.m.a, t(x.a))
    resids <- y - x.a%*%ls.est
  }

  beta.hat <- b.sp.e%*%ls.est[-1]/dim(x)[2]
  
  return(list(bh = beta.hat, beta = ls.est[-1], alpha = ls.est[1], lam = lambda1, 
              resids = resids, b.sp.e = b.sp.e))
}
