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

### Figure 3 in the paper ###
library(reshape2)
library(ggplot2)
data <- data.frame(t(x))
data$id <- 1:nrow(data)
#reshape to long format
plot_data <- melt(data, id.var="id" )
#plot
gr <- ggplot(plot_data, aes(x=id, y=value, group=variable, colour=variable)) + geom_line(col = "gray", size = 1.2)
gr <- gr + theme_bw(base_size = 40) + theme(plot.margin = margin(t = 0,  r = 0,  b = 0, l = 0))  + labs(x = "t", y = "")
gr

data <- data.frame(y3)
gr <- ggplot(data, aes(x = y3)) + geom_histogram(fill = "gray", col = "black") + theme_bw(base_size = 40) + labs(x = "", y = "Frequency")
gr <- gr +  theme(plot.margin = margin(t = 0,  r = 0,  b = 0, l = 0))
gr

### Figure 4 in the paper ###

fit4 <- m.pen.sp(x = x, y = y4, norder = 4, q = 2, k = 4.685)
plot(fit4$bh, type = "l", lwd = 3, col = "blue", cex.lab = 2, cex.axis = 2, yaxt = "n", xaxt = "n") ; grid()
fit41 <- ls.pen.sp(x = x, y = y4)
lines(fit41$bh, lwd = 3, col = "red", lty = 2)

require(ggplot2)
data <- data.frame(fit.r = fit4$bh, fit.ls = fit41$bh)
data$x <- 1:dim(x)[2]/dim(x)[2]
gr <- ggplot(data = data, aes(x = x, y = fit)) + geom_line(aes(x = x, y = fit.r), colour = "blue",  size = 1.2) + theme_bw(base_size = 40)
gr <- gr + theme(plot.margin = margin(t = 0,  r = 0,  b = 0, l = 0))  + labs(x = "", y = "")
gr <- gr +  geom_line(aes(x = x, y = fit.ls), colour = "red",  size = 1.2, linetype = "longdash")
gr

fit12 <- m.pen.sp(x = x, y = y12, norder = 4)
plot(fit12$bh, type = "l", lwd = 3, col = "blue", cex.lab = 2, cex.axis = 2, yaxt = "n", xaxt = "n", xlab = "", ylab = "") ; grid()
fit121 <- ls.pen.sp(x = x, y = y12)
lines(fit121$bh, lwd = 3, col = "red", lty = 2)

require(ggplot2)
data <- data.frame(fit.r = fit12$bh, fit.ls = fit121$fhat)
data$x <- 1:dim(x)[2]/dim(x)[2]
gr <- ggplot(data = data, aes(x = x, y = fit)) + geom_line(aes(x = x, y = fit.r), colour = "blue",  size = 1.2) + theme_bw(base_size = 40)
gr <- gr + theme(plot.margin = margin(t = 0,  r = 0,  b = 0, l = 0))  + labs(x = "", y = "")
gr <- gr +  geom_line(aes(x = x, y = fit.ls), colour = "red",  size = 1.2, linetype = "longdash")
gr

# Upper trimmed mean
u.tr <- function(x, alpha){
  n <- length(x)
  m <- round(alpha*n)
  x.s <- sort(x)
  u.tr <- mean(x.s[1:(n-m)])
  return(u.tr)
}

## Predict function for PMM
predict.mpen <- function(fit.r, x){
  pred.values <- fit.r$alpha + x%*%fit.r$bh
  return(pred.values)
}

## Predict function for PLS
predict.ls <- function(fit.ls, x){
  pred.values <- fit.ls$alpha + x%*%fit.ls$bh
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
    fit.r <- m.pen.sp(x = x.train, y = y.train, norder = 4, k = 3.44, q = 2, nbasis = round(min(40, length(y.train)/4)))
    fit.ls <- ls.pen.sp(x = x.train, y = y.train,  norder = 4, q = 2, nbasis = round(min(40, length(y.train)/4)))
    pred.values.r <- predict.mpen( fit.r, x.test)
    pred.values.ls <- predict.ls(fit.ls, x.test)
    rmspe.r[j] <- sqrt(u.tr( (y.test-pred.values.r)^2, 0.1 ))
    rmpse.ls[j] <- sqrt(u.tr( (y.test-pred.values.ls)^2,0.1 ))
  }
  rmspe.r <- mean(rmspe.r)
  rmpse.ls <- mean(rmpse.ls)
  return(list(rmspe.r, rmpse.ls))
}

cvrep <- 50
cv1 <- matrix(NA, nrow = 2, ncol = cvrep)
for(j in 1:cvrep){
  cv1[, j] <- unlist(cv.function(x, y1, 5))
}
rowMeans(cv1)

cv2 <- matrix(NA, nrow = 2, ncol = cvrep)
for(j in 1:cvrep){
  cv2[, j] <- unlist(cv.function(x, y2, 5))
}
rowMeans(cv2)

cv3 <- matrix(NA, nrow = 2, ncol = cvrep)
for(j in 1:cvrep){
  cv3[, j] <- unlist(cv.function(x, y3, 5))
}
rowMeans(cv3)

cv4 <- matrix(NA, nrow = 2, ncol = cvrep)
for(j in 1:cvrep){
  cv4[, j] <- unlist(cv.function(x, y4, 5))
}
rowMeans(cv4)

cv5 <- matrix(NA, nrow = 2, ncol = cvrep)
for(j in 1:cvrep){
  cv5[, j] <- unlist(cv.function(x, y5, 5))
}
rowMeans(cv5)

cv6 <- matrix(NA, nrow = 2, ncol = cvrep)
for(j in 1:cvrep){
  cv6[, j] <- unlist(cv.function(x, y6, 5))
}
rowMeans(cv6)

cv7 <- matrix(NA, nrow = 2, ncol = cvrep)
for(j in 1:cvrep){
  cv7[, j] <- unlist(cv.function(x, y7, 5))
}
rowMeans(cv7)

cv8 <- matrix(NA, nrow = 2, ncol = cvrep)
for(j in 1:cvrep){
  cv8[, j] <- unlist(cv.function(x, y8, 5))
}
rowMeans(cv8)

cv9 <- matrix(NA, nrow = 2, ncol = cvrep)
for(j in 1:cvrep){
  cv9[, j] <- unlist(cv.function(x, y9, 5))
}
rowMeans(cv9)

cv10 <- matrix(NA, nrow = 2, ncol = cvrep)
for(j in 1:cvrep){
  cv10[, j] <- unlist(cv.function(x, y10, 5))
}
rowMeans(cv10)

cv11 <- matrix(NA, nrow = 2, ncol = cvrep)
for(j in 1:cvrep){
  cv11[, j] <- unlist(cv.function(x, y11, 5))
}
rowMeans(cv11)

cv12 <- matrix(NA, nrow = 2, ncol = cvrep)
for(j in 1:cvrep){
  cv12[, j] <- unlist(cv.function(x, y12, 5))
}
rowMeans(cv12)

  cv13 <- matrix(NA, nrow = 2, ncol = cvrep)
for(j in 1:cvrep){
  cv13[, j] <- unlist(cv.function(x, y13, 5))
}
rowMeans(cv13)

