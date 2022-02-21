require(robustbase)
require(fda)

m.step <- function(x, y, tol = 1e-08, maxit = 200, resids.in, scale, k, p.m, lambda){
  
  weight.function <- function(x) ifelse(x==0, 1, Mpsi(x/scale, psi = "bisquare", cc = k)/(x/scale) )
  istop = 0
  ic = 0
  
  p.m <- rbind( rep(0, (dim(x)[2]+1) ), cbind(rep(0, dim(x)[2]), p.m) )
  x <- cbind(rep(1, dim(x)[1]), x)
  
  while(istop == 0 & ic <= maxit){
    ic = ic + 1
    weights.in <- weight.function(resids.in)
    weights.in.nz <- weights.in[weights.in!=0]
    x.nz <- x[weights.in!=0,]
    y.nz <- y[weights.in!=0]
    
    x.ws <- as.matrix(scale(t( x.nz), center = FALSE, scale = 1/weights.in.nz))
    M1 <-  x.ws%*%x.nz + 2*scale^2*lambda*p.m
    M2 <- x.ws%*%as.matrix(y.nz)
    v1 = SparseM::solve( M1, M2)
    resids1 <- as.vector(y-x%*%v1)
    check = max( abs(resids1-resids.in)/scale) 
    if(check < tol){istop =1}
    resids.in <- resids1
  }
  weights1 = weight.function(resids1)
  fitted.values = as.vector(x%*%v1)
  resids = y - fitted.values
  hat.tr <- diag( x%*%solve(t(x)%*%diag(weights1)%*%x + 2*scale^2*lambda*p.m, t(x)%*%diag(weights1))  )
  return(list(beta = v1[-1], alpha = v1[1], hat.tr = hat.tr, resids = resids, weights = weights1))
}

m.pen.sp <- function(x, y, nbasis = 40, norder = 4,  k = 4.685, q  = 2){
  
  x <- as.matrix(x)
  b.sp <- create.bspline.basis(c(0, 1), nbasis = nbasis, norder = norder)
  b.sp.e <- eval.basis(seq(0, 1, len = dim(x)[2]), b.sp)
  # p.m <- t( diff(diag( (nbasis)) , differences = q ) )%*%diff(diag( (nbasis ) ), differences = q )
  p.m <- bsplinepen(b.sp, Lfdobj=q )
  
  x.p <- x%*%b.sp.e/dim(x)[2]
  n <- length(y)
  
  rob.ctrl <- lmrob.control(trace.level = 0,       
                            nResample   = 10000,      
                            tuning.psi  = k,     
                            subsampling = 'simple',  
                            rel.tol     = 1e-07,      
                            refine.tol  = 1e-07,      
                            k.max       = 5e3,       
                            maxit.scale = 5e3,       
                            max.it      = 5e3)
  
  S.est <- lmrob.S(x = x.p, y = y, control = rob.ctrl)
  resids.in <- S.est$residuals
  scale <- S.est$scale
  
  GCV <- function(lambda){
    fit.r <- m.step(x.p, y, resids.in = resids.in, k = k, p.m = p.m, scale = scale, lambda = lambda)
    # GCV.scores <- scaleTau2(fit.r$resids/(1-fit.r$hat.tr), c2 = 5)^2
    # GCV.scores <- mean( fit.r$weights*fit.r$resids^2/((1-fit.r$hat.tr)^2) )
    # GCV.scores <- weighted.mean(fit.r$resids^2/((1-fit.r$hat.tr)^2), wt = fit.r$weights/n  )
    # sdw <- sqrt(sum(fit.r$weights*(fit.r$resids^2/((1-fit.r$hat.tr)^2)-GCV.scores)^2))/sqrt(n)
    weights <- Mpsi(fit.r$resids/scale, cc = k, psi = "bisquare", deriv = 1)/scale^2
    NP <- scale(t(x.p), center = FALSE, scale = c(1/weights))%*%x.p
    hessian.m <-  NP + 2*lambda*p.m
    gr <- -Mpsi(fit.r$resids/scale, cc = k, psi = "bisquare", deriv = 0)/scale
    Q.m <-  scale(t(x.p), center = FALSE, scale = 1/(gr*c(fit.r$resids/scale^2)) )%*%x.p
    hat.tr <- sum(diag(  -solve(hessian.m)%*%Q.m   )) # GIC
    # hat.tr <- sum(diag(  solve(hessian.m)%*%NP   )) # AIC
    GCV.scores <- 2*sum(Mpsi(fit.r$resids/scale, cc = k, psi = "bisquare", deriv = -1))   + log(n)*hat.tr
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
  fit.r <- m.step(x.p, y, resids.in = resids.in, k = k, p.m = p.m, scale = scale, lambda = lambda1)
  
  # sdw <- scaleTau2(fit.r$resids^2/(1-fit.r$hat.tr)^2-GCV(lambda1), c2 = 5)/sqrt(n)
  # # M <- sum(fit.r$weights !=0)
  # # sdw <- sqrt(mean(fit.r$weights*(fit.r$resids^2/((1-fit.r$hat.tr)^2)-GCV(lambda1))^2)/((M-1)*sum(fit.r$weights/n)/M))/sqrt(n)
  # # sdw <- sqrt(sum(fit.r$weights*(fit.r$resids^2/((1-fit.r$hat.tr)^2)-GCV(lambda1))^2)/sum(fit.r$weights))/sqrt(n)
  # thr <- GCV(lambda1) + 2*sdw
  # lambda1 <- max( lambda.cand[lambda.e <= thr]  )
  # fit.r <- m.step(x.p, y, resids.in = resids.in, k = k, p.m = p.m, scale = scale, lambda = lambda1)
  
  beta.hat <- b.sp.e%*%fit.r$beta
  return(list(bh = beta.hat, beta = fit.r$beta, alpha = fit.r$alpha,  weights = fit.r$weights, lam = lambda1, 
              resids = fit.r$resids, b.sp.e = b.sp.e))
}
