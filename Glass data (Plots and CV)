require(robustbase)
require(chemometrics)
require(cellWise)
require(refund)

data(data_glass)
data(glass)
x <- data_glass
x <- as.matrix(x)
matplot(t(x), col = "gray", type = "l", lwd = 3, cex.lab = 2, cex.axis = 2) ; grid()
y1 <- glass[, 1]
y2 <- glass[, 2]
y3 <- glass[, 3]
y4 <- glass[, 4]
y5 <- glass[, 5]
y6 <- glass[, 6]
y7 <- glass[, 7]
y8 <- glass[, 8]
y9 <- glass[, 9]
y10 <- glass[, 10]
y11 <- glass[, 11]
y12 <- glass[, 12]
y13 <- glass[, 13]

par(mar = c(3.1, 3.1, 3.1, 1.1)); par(mgp = c(3.8, 1, 0))

hist(y3, col = "gray", main = "", cex.lab = 2, cex.axis = 2) ; grid()
hist(y5, col = "gray", main = "", cex.lab = 2, cex.axis = 2) ; grid()
hist(y12, col = "gray", main = "", cex.lab = 2, cex.axis = 2) ; grid()

fit4 <- m.pen.sp(x = x, y = y4, norder = 4, nbasis = 40, interval = c(2e-02, 9e-02), q =2 )
plot(fit4$bh/dim(x)[2], type = "l", lwd = 3, col = "blue", cex.lab = 2, cex.axis = 2, yaxt = "n", xaxt = "n") ; grid()
fit41 <- fpcr(y = y4, xfuncs = x)
lines(fit41$fhat, lwd = 3, col = "red", lty = 2)

fit12 <- m.pen.sp(x = x, y = y12, norder = 4, nbasis = 40, interval = c(9e-02, 5e-01))
plot(fit12$bh/dim(x)[2], type = "l", lwd = 3, col = "blue", cex.lab = 2, cex.axis = 2, yaxt = "n", xaxt = "n", xlab = "", ylab = "") ; grid()
fit121 <- fpcr(y = y12, xfuncs = x)
lines(fit121$fhat, lwd = 3, col = "red", lty = 2)

# Upper trimmed mean
u.tr <- function(x, alpha){
  n <- length(x)
  m <- round(alpha*n)
  x.s <- sort(x)
  u.tr <- mean(x.s[1:(n-m)])
  return(u.tr)
}

cv.function<- function(x, y, nfolds, int){
  x <- as.matrix(x)
  splits <- split(1:dim(x)[1], sample( rep(1:nfolds, times= rep( round(dim(x)[1]/nfolds), nfolds ) )) )
  
  rmspe.r <- rep(NA, nfolds)
  rmpse.ls <- rep(NA, nfolds)
  
  for(j in 1:nfolds){
    x.test <- x[splits[[j]], ]
    x.train <- x[-splits[[j]], ]
    y.test <- y[splits[[j]]]
    y.train <- y[-splits[[j]]]
    fit.r <- m.pen.sp(x = x.train, y = y.train, norder = 4, k = 4.685, nbasis = round(min(40, length(y.train)/4)), 
                      interval =  int)
    fit.ls <- fpcr(y = y.train, x = x.train)
    pred.values.r <- predict.mpen( fit.r, x.test)
    pred.values.ls <- as.numeric(fit.ls$undecor.coef) +  x.test%*%fit.ls$fhat
    rmspe.r[j] <- sqrt(u.tr( (y.test-pred.values.r)^2, 0.1 ))
    rmpse.ls[j] <- sqrt(u.tr( (y.test-pred.values.ls)^2,0.1 ))
  }
  rmspe.r <- mean(rmspe.r)
  rmpse.ls <- mean(rmpse.ls)
  return(list(rmspe.r, rmpse.ls))
}

# Results appearing in Table 2 of the paper
# The veracity of the selected intervals for the smoothing parameter follows from inspection of plots of the RCV function

cv.function(x, y1, 5, c(5e-05, 1e-04))
cv.function(x, y2, 5, c(8, 16))
cv.function(x, y3, 5, c(9e-05, 6e-04))
cv.function(x, y4, 5, c(2e-02, 9e-02))
cv.function(x, y5, 5, c(2e-04, 9e-04))
cv.function(x, y6, 5, c(9e-03, 6e-02))
cv.function(x, y7, 5,  c(2,7))
cv.function(x, y8, 5, c(0.1, 0.9))
cv.function(x, y9, 5, c(8e-02, 5e-01))
cv.function(x, y10, 5, c(9e-03, 6e-02))
cv.function(x, y11, 5, c(8e-05, 5e-04))
cv.function(x, y12, 5, c(9e-02, 5e-01))
cv.function(x, y13, 5, c(70, 110))
