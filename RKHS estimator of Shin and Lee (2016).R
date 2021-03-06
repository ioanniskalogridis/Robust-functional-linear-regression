#### RKHS estimator of Shin et al. (2016)
#### The main function is based on the Tukey bisquare loss and initial values are provided with a Huber fit
require(pracma)
require(MASS)

kernel.rkhs <- function(s,t){
  kernel.rkhs = pmin(s, t)
  return(kernel.rkhs)
}

huber.reg.rkhs <- function(X, y, k = 1.345, tol = 1e-08, lambda, A.matrix, dom, maxit = 200, scale, 
                           krm, Pred.big, resids.in){
  
  rho.huber <- function(x, k = 1.345) ifelse(abs(x)<= k, x^2/2, k*abs(x)-k^2/2)
  psi.huber <- function(x, k = 1.345) ifelse(abs(x)<= k, x, k*sign(x))
  weights.huber <- function(x, k = 1.345) ifelse( x==0, 1, psi.huber(x, k = k)/x)
  
  ic = 0
  istop = 0
  
  while(istop == 0 & ic <= maxit){
    ic = ic + 1
    weights.in <- weights.huber(resids.in/scale, k = k)
    X.w <- as.matrix(scale(t(Pred.big), center = FALSE, scale = 1/weights.in))
    M1 <-  X.w%*%Pred.big + 2*scale^2*lambda*A.matrix
    M2 <- X.w%*%as.matrix(y)
    v1 = ginv(M1)%*%M2
    resids1 <- as.vector(y-Pred.big%*%v1)
    check = max( abs(resids1-resids.in)/scale )
    if(check < tol){istop =1}
    resids.in <- resids1
  }
  weights1 = weights.huber(resids1/scale, k = k)
  fitted.values = as.vector(Pred.big%*%v1)
  resids = y - fitted.values
  hat.matrix = Pred.big%*%M1%*%t(Pred.big)%*%diag(weights1)
  est.beta = as.vector(v1)
  beta = as.vector( est.beta[2] + est.beta[3]*dom/max(dom) + t(est.beta[4:length(est.beta)])%*%X%*%krm/length(dom) )
  
  return(list(beta = beta, fitted = fitted.values, ic = ic, check = check, hat.matrix = hat.matrix,
              weights = weights1, resids = resids, interc = v1[1], est.beta = v1))
}

flm.rkhs.rob.h <- function(X, y, dom, k = 1.345, lambda= NULL){
  
  T.m = matrix(NA, n, 2)
  T.m[, 1] = apply(X, 1, FUN = function(y) trapz(dom/max(dom), y ) )
  T.m[, 2] = apply(X, 1, FUN = function(y) trapz(dom/max(dom), y*dom/max(dom) ))
  
  kernel.rkhs.m = outer(dom/max(dom), dom/max(dom), FUN = kernel.rkhs)
  Sigma = (1/length(dom)^2)*X%*%kernel.rkhs.m%*%t(X)
  Pred.big = cbind(rep(1, length(y)),T.m, Sigma )
  A.matrix = cbind(matrix(0, dim(Sigma)[2]+3, 3), rbind(matrix(0, 3, dim(Sigma)[2]), Sigma)  )
  
  est.ls <- ginv(t(Pred.big)%*%Pred.big + 1e-07*A.matrix)%*%t(Pred.big)%*%y
  resids.prelim <- y - as.vector(Pred.big%*%est.ls)
  scale = mad(resids.prelim)
  
  GCV <- function(lambda){
    fit.r <- huber.reg.rkhs(X = X,  y = y, lambda = lambda, A.matrix = A.matrix, resids.in = resids.prelim, scale = scale,
                            k =  k, dom = dom, Pred.big = Pred.big, krm = kernel.rkhs.m)
    GCV.scores <- ( mean( fit.r$weights*(y-fit.r$fitted)^2  ))/((1-mean(diag(fit.r$hat.matrix)))^2)
    return(GCV.scores)
  }
  
  lambda.cand <- c(1e-11, 1e-10, 1e-09, 2e-08, 6e-08, 9e-08,  2e-07, 6e-07, 9e-07,  2e-06, 6e-06, 9e-06,
                   2e-05, 6e-05, 9e-05,  2e-04, 6e-04, 9e-04, 2e-03, 6e-03, 9e-03,  2e-02, 6e-02, 9e-02,
                   2e-01, 6e-01, 9e-01, 3, 50)
  lambda.e <- rep(NA, length(lambda.cand))
  for(e in 1:length(lambda.cand)){
    lambda.e[e] <- try(GCV(lambda.cand[e]), silent = TRUE)
  }
  lambda.e <- as.numeric(lambda.e)
  wm <- which.min(lambda.e)
  if(wm==1){
    wm <- 2
  } else if(wm==length(lambda.cand)){wm <- (length(lambda.cand)-1)}
  lambda1 <- optimize(f = GCV, interval = c(lambda.cand[wm-1], lambda.cand[wm+1]) , tol = 1e-14)$minimum
  
  fit.rf <- huber.reg.rkhs(X = X,  y = y, lambda = lambda1, A.matrix = A.matrix, resids.in = resids.prelim, scale = scale,
                          k =  k, dom = dom, Pred.big = Pred.big, krm = kernel.rkhs.m)
  beta = fit.rf$beta
  return(list(beta = beta, lambda = lambda1, resids = fit.rf$resids))
}

tukey.reg.rkhs <- function(X, y, k = 4.685, tol = 1e-08, lambda, A.matrix, dom, maxit = 200, scale, 
                           krm, Pred.big, resids.in){
  
  rho.tukey <- function(x, k = 4.685) ifelse(abs(x)<= k, (1-(1-(x/k)^2)^3), 1)
  psi.tukey <- function(x, k = 4.685) ifelse(abs(x)<= k, x*(1-(x/k)^2)^2, 0)
  weights.tukey <- function(x, k = 4.685) ifelse( x==0, 1, psi.tukey(x, k = k)/x)
  
  ic = 0
  istop = 0
  
  while(istop == 0 & ic <= maxit){
    ic = ic + 1
    weights.in <- weights.tukey(resids.in/scale, k = k)
    
    X.w <- as.matrix(scale(t(Pred.big), center = FALSE, scale = 1/weights.in))
    M1 <-  X.w%*%Pred.big + 2*scale^2*lambda*A.matrix
    M2 <- X.w%*%as.matrix(y)
    v1 = ginv(M1)%*%M2
    resids1 <- as.vector(y-Pred.big%*%v1)
    check = max( abs(resids1-resids.in)/scale )
    if(check < tol){istop =1}
    resids.in <- resids1
  }
  weights1 = weights.tukey(resids1/scale, k = k)
  fitted.values = as.vector(Pred.big%*%v1)
  resids = y - fitted.values
  hat.matrix = Pred.big%*%M1%*%t(Pred.big)%*%diag(weights1)
  est.beta = as.vector(v1)
  beta = as.vector( est.beta[2] + est.beta[3]*dom/max(dom) + t(est.beta[4:length(est.beta)])%*%X%*%krm/length(dom) )
  
  return(list(beta = beta, fitted = fitted.values, ic = ic, check = check, hat.matrix = hat.matrix,
              weights = weights1, resids = resids, interc = v1[1], est.beta = v1))
}

flm.rkhs.rob.t <- function(X, y, dom, k = 4.6855, lambda= NULL){
  
  T.m = matrix(NA, n, 2)
  T.m[, 1] = apply(X, 1, FUN = function(y) trapz(dom/max(dom), y ) )
  T.m[, 2] = apply(X, 1, FUN = function(y) trapz(dom/max(dom), y*dom/max(dom) ))
  
  kernel.rkhs.m = outer(dom/max(dom), dom/max(dom), FUN = kernel.rkhs)
  Sigma = (1/length(dom)^2)*X%*%kernel.rkhs.m%*%t(X)
  Pred.big = cbind(rep(1, length(y)),T.m, Sigma )
  A.matrix = cbind(matrix(0, dim(Sigma)[2]+3, 3), rbind(matrix(0, 3, dim(Sigma)[2]), Sigma)  )
  
  est.h <- flm.rkhs.rob.h(X = X, y = y, dom = dom)
  resids.prelim<- as.vector(est.h$resids)
  scale <- mad(resids.prelim)
  
  GCV <- function(lambda){
    fit.r <- tukey.reg.rkhs(X = X,  y = y, lambda = lambda, A.matrix = A.matrix, resids.in = resids.prelim, scale = scale,
                            k =  k, dom = dom, Pred.big = Pred.big, krm = kernel.rkhs.m)
    GCV.scores <- ( mean( fit.r$weights*(y-fit.r$fitted)^2  ))/((1-mean(diag(fit.r$hat.matrix)))^2)
    return(GCV.scores)
  }
  lambda.cand <- c(1e-11, 1e-10, 1e-09, 2e-08, 6e-08, 9e-08,  2e-07, 6e-07, 9e-07,  2e-06, 6e-06, 9e-06,
                   2e-05, 6e-05, 9e-05,  2e-04, 6e-04, 9e-04, 2e-03, 6e-03, 9e-03,  2e-02, 6e-02, 9e-02,
                   2e-01, 6e-01, 9e-01, 3, 50)
  lambda.e <- rep(NA, length(lambda.cand))
  for(e in 1:length(lambda.cand)){
    lambda.e[e] <- try(GCV(lambda.cand[e]), silent = TRUE)
  }
  lambda.e <- as.numeric(lambda.e)
  wm <- which.min(lambda.e)
  if(wm==1){
    wm <- 2
  } else if(wm==length(lambda.cand)){wm <- (length(lambda.cand)-1)}
  lambda1 <- optimize(f = GCV, interval = c(lambda.cand[wm-1], lambda.cand[wm+1]) , tol = 1e-14)$minimum

  fit.rf <- tukey.reg.rkhs(X = X,  y = y, lambda = lambda1, A.matrix = A.matrix, resids.in = resids.prelim, scale = scale,
                           k =  k, dom = dom, Pred.big = Pred.big, krm = kernel.rkhs.m)
  beta = fit.rf$beta
  return(list(beta = beta, lambda = lambda1))
}


