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

### Figure 3 in the paper ####
library(reshape2)
library(ggplot2)
data <- data.frame(t(x))
data$id <- 1:nrow(data)
#reshape to long format
plot_data <- melt(data, id.var="id" )
#plot
gr <- ggplot(plot_data, aes(x=id, y=value, group=variable, colour=variable)) + geom_line(col = "gray", size = 1.2)
gr <- gr + theme_bw(base_size = 40) + theme(plot.margin = margin(t = 0,  r = 0,  b = 0, l = 0))  + labs(x = "t", y = "")
gr <- gr + stat_function(aes(x = id, y = value), fun = mu, size = 1.4, colour = "black") #+ ylim(-6, 6)
gr

data <- data.frame(y3)
gr <- ggplot(data, aes(x = y3)) + geom_histogram(fill = "gray") + theme_bw(base_size = 40) + labs(x = "", y = "Frequency")
gr <- gr +  theme(plot.margin = margin(t = 0,  r = 0,  b = 0, l = 0))
gr

fit1 <- m.pen.sp(x = x, y = y1, norder = 4, q = 2, k = 4.685)
plot(fit1$bh, type = "l", lwd = 3, col = "blue", cex.lab = 2, cex.axis = 2, yaxt = "n", xaxt = "n") ; grid()
fit11 <- fpcr(y = y1, xfuncs = x)
lines(fit11$fhat, lwd = 3, col = "red", lty = 2)

fit2 <- m.pen.sp(x = x, y = y2, norder = 4, nbasis = 40, q = 2, k = 4.685)
plot(fit2$bh, type = "l", lwd = 3, col = "blue", cex.lab = 2, cex.axis = 2, yaxt = "n", xaxt = "n") ; grid()
fit21 <- fpcr(y = y2, xfuncs = x)
lines(fit21$fhat, lwd = 3, col = "red", lty = 2)

fit3 <- m.pen.sp(x = x, y = y3, norder = 4, nbasis = 40, q = 2, k = 4.685)
plot(fit3$bh, type = "l", lwd = 3, col = "blue", cex.lab = 2, cex.axis = 2, yaxt = "n", xaxt = "n") ; grid()
fit31 <- fpcr(y = y3, xfuncs = x)
lines(fit31$fhat, lwd = 3, col = "red", lty = 2)

fit4 <- m.pen.sp(x = x, y = y4, norder = 4, q = 2, k = 4.685)
plot(fit4$bh, type = "l", lwd = 3, col = "blue", cex.lab = 2, cex.axis = 2, yaxt = "n", xaxt = "n") ; grid()
fit41 <- fpcr(y = y4, xfuncs = x)
lines(fit41$fhat, lwd = 3, col = "red", lty = 2)

fit5 <- m.pen.sp(x = x, y = y5, norder = 4, q = 2, k = 4.685)
plot(fit5$bh, type = "l", lwd = 3, col = "blue", cex.lab = 2, cex.axis = 2, yaxt = "n", xaxt = "n") ; grid()
fit51 <- fpcr(y = y6, xfuncs = x)
lines(fit51$fhat, lwd = 3, col = "red", lty = 2)

fit6 <- m.pen.sp(x = x, y = y6, norder = 4, nbasis = 40, q = 2, k = 4.685)
plot(fit6$bh, type = "l", lwd = 3, col = "blue", cex.lab = 2, cex.axis = 2, yaxt = "n", xaxt = "n") ; grid()
fit61 <- fpcr(y = y6, xfuncs = x)
lines(fit61$fhat, lwd = 3, col = "red", lty = 2)

fit7 <- m.pen.sp(x = x, y = y7, norder = 4, nbasis = 40, q = 2, k = 4.685)
plot(fit7$bh, type = "l", lwd = 3, col = "blue", cex.lab = 2, cex.axis = 2, yaxt = "n", xaxt = "n") ; grid()
fit71 <- fpcr(y = y7, xfuncs = x)
lines(fit71$fhat, lwd = 3, col = "red", lty = 2)

fit8 <- m.pen.sp(x = x, y = y8, norder = 4, nbasis = 40, q = 2, k = 4.685)
plot(fit8$bh/dim(x)[2], type = "l", lwd = 3, col = "blue", cex.lab = 2, cex.axis = 2, yaxt = "n", xaxt = "n") ; grid()
fit81 <- fpcr(y = y8, xfuncs = x)
lines(fit81$fhat, lwd = 3, col = "red", lty = 2)

fit9 <- m.pen.sp(x = x, y = y9, norder = 4, nbasis = 40, q = 2, k = 4.685)
plot(fit9$bh/dim(x)[2], type = "l", lwd = 3, col = "blue", cex.lab = 2, cex.axis = 2, yaxt = "n", xaxt = "n") ; grid()
fit91 <- fpcr(y = y9, xfuncs = x)
lines(fit91$fhat, lwd = 3, col = "red", lty = 2)

fit10 <- m.pen.sp(x = x, y = y10, norder = 4, nbasis = 40, q = 2, k = 4.685)
plot(fit10$bh, type = "l", lwd = 3, col = "blue", cex.lab = 2, cex.axis = 2, yaxt = "n", xaxt = "n") ; grid()
fit101 <- fpcr(y = y10, xfuncs = x)
lines(fit101$fhat, lwd = 3, col = "red", lty = 2)

fit11 <- m.pen.sp(x = x, y = y11, norder = 4, q = 2)
plot(fit11$bh, type = "l", lwd = 3, col = "blue", cex.lab = 2, cex.axis = 2, yaxt = "n", xaxt = "n") ; grid()
fit111 <- fpcr(y = y11, xfuncs = x)
lines(fit111$fhat, lwd = 3, col = "red", lty = 2)

fit12 <- m.pen.sp(x = x, y = y12, norder = 4)
plot(fit12$bh, type = "l", lwd = 3, col = "blue", cex.lab = 2, cex.axis = 2, yaxt = "n", xaxt = "n", xlab = "", ylab = "") ; grid()
fit121 <- fpcr(y = y12, xfuncs = x)
lines(fit121$fhat, lwd = 3, col = "red", lty = 2)

fit13 <- m.pen.sp(x = x, y = y13, norder = 4, n.se = 0)
plot(fit13$bh, type = "l", lwd = 3, col = "blue", cex.lab = 2, cex.axis = 2, yaxt = "n", xaxt = "n", xlab = "", ylab = "") ; grid()
fit131 <- fpcr(y = y13, xfuncs = x)
lines(fit131$fhat, lwd = 3, col = "red", lty = 2)






# Upper trimmed mean
u.tr <- function(x, alpha){
  n <- length(x)
  m <- round(alpha*n)
  x.s <- sort(x)
  u.tr <- mean(x.s[1:(n-m)])
  return(u.tr)
}

predict.mpen <- function(fit.r, x){
  pred.values <- fit.r$alpha + x%*%fit.r$bh
  return(pred.values)
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
    fit.r <- m.pen.sp(x = x.train, y = y.train, norder = 4, k = 4.685, q = 2, nbasis = round(min(40, length(y.train)/4)))
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

# cv.function(x, y1, 5, c(5e-05, 1e-04))
# cv.function(x, y2, 5, c(8, 16))
# cv.function(x, y3, 5, c(9e-05, 6e-04))
# cv.function(x, y4, 5, c(2e-02, 9e-02))
# cv.function(x, y5, 5, c(2e-04, 9e-04))
# cv.function(x, y6, 5, c(9e-03, 6e-02))
# cv.function(x, y7, 5,  c(2,7))
# cv.function(x, y8, 5, c(0.1, 0.9))
# cv.function(x, y9, 5, c(8e-02, 5e-01))
# cv.function(x, y10, 5, c(9e-03, 6e-02))
# cv.function(x, y11, 5, c(8e-05, 5e-04))
# cv.function(x, y12, 5, c(9e-02, 5e-01))
# cv.function(x, y13, 5, c(70, 110))

cv.function(x, y1, 5)
cv.function(x, y2, 5)
cv.function(x, y3, 5)
cv.function(x, y4, 5)
cv.function(x, y5, 5)
cv.function(x, y6, 5)
cv.function(x, y7, 5)
cv.function(x, y8, 5)
cv.function(x, y9, 5)
cv.function(x, y10, 5)
cv.function(x, y11, 5)
cv.function(x, y12, 5)
cv.function(x, y13, 5)


install.packages("fda.usc")
library(fda.usc)
data(tecator)
x <- tecator$absorp.fdata
x <- x$data

y <- tecator$y
y <- y[, 1]
matplot(t(x), type = "l", col = "gray", lty = 1, lwd = 3)
fit1 <-  m.pen.sp(x = x, y = y, norder = 4, nbasis = 40, q = 2, k = 3.44)
fit.ls <- ls.pen.sp(x, y, norder = 4, nbasis = NULL, q  = 2)
plot(fit1$bh, type = "l", col = "blue", lwd = 3)
lines(fit.ls$bh, type = "l", col = "red", lwd = 3)

y <- tecator$y
y <- y[, 2]
fit1 <-  m.pen.sp(x = x, y = y, norder = 4, nbasis = 40, q = 2, k = 3.44)
fit.ls <- ls.pen.sp(x, y, norder = 4, nbasis = NULL, q  = 2)
plot(fit1$bh, type = "l", col = "blue", lwd = 3)
lines(fit.ls$bh, type = "l", col = "red", lwd = 3)

y <- tecator$y
y <- y[, 3]
fit1 <-  m.pen.sp(x = x, y = y, norder = 4, nbasis = 40, q = 2, k = 3.44)
fit.ls <- ls.pen.sp(x, y, norder = 4, nbasis = NULL, q  = 2)
plot(fit1$bh, type = "l", col = "blue", lwd = 3)
lines(fit.ls$bh, type = "l", col = "red", lwd = 3)
