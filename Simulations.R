require(refund)
require(fda)
require(MASS)
require(EnvStats)

Nrep <- 1000
n <- 150
p <- 100
grid <- seq(0, 1, length = p)

alpha <- sin(2*pi*grid)
# alpha = grid^2*dnorm(grid)
# alpha = 1/(1+exp(-20*(grid-0.5)))
# alpha = -dnorm(grid, mean=.2, sd=.03) + 3*dnorm(grid, mean=.5, sd=.03) + dnorm(grid, mean=.75, sd=.04)

mse.fpcr <- rep(NA, Nrep)
mse.ls <- rep(NA, Nrep)
mse.mpen <- rep(NA, Nrep)
mse.munp <- rep(NA, Nrep)
mse.smsp <- rep(NA, Nrep)

matr.fpcr <- matrix(NA, nrow = p, ncol = Nrep)
matr.ls <- matrix(NA, nrow = p, ncol = Nrep)
matr.mpen <- matrix(NA, nrow = p, ncol = Nrep)
matr.munp <- matrix(NA, nrow = p, ncol = Nrep)
matr.smsp <- matrix(NA, nrow = p, ncol = Nrep)

for(f in 1:Nrep){
  message('Iter = ', f, ' of ', Nrep)
  X0 <- sqrt(2)*matrix(rnorm(n*p), nrow = n, ncol = p)
  # X0 <- sqrt(2)*matrix(rt(n*p, df = 3), nrow = n, ncol = p)
  for(i in 1:n){
    for(j in 2:50){
      X0[i, ] <- X0[i, ] + j^{-1}*rnorm(1, 0, 1)*sqrt(2)*sapply(grid, FUN= function(x) cos((j-1)*pi*x))
      # X0[i, ] <- X0[i, ] +  1*j^{-1}*rt(1, df = 3)*sqrt(2)*sapply(grid, FUN= function(x) cos((j-1)*pi*x)) # For leverage contamination
    }
  }
  y0 = X0%*%alpha
  # y <- y0 + rnorm(n)
  # y <- y0 + rt(n, df = 3)
  y <- y0 + rnormMix(n, mean1 = 0, sd = 1, mean2 = 14, sd2 = 1, p.mix = 0.1)
  fit.mpen <- m.pen.sp(x = X0, y = y, nbasis = round(min(n/4, 40)), n.se = 0)
  fit.fpcr <- fpcr(y, xfuncs = X0, method = "GCV.Cp", pve = 0.999999, nbasis = 38)
  fit.ls <- ls.pen.sp(x = X0, y = y, nbasis = round(min(n/4, 40)), n.se = 0)
  fit.smsp <- m.sm.sp(x = X0, y = y, t = grid)
  fit.munp <- m.sp(x = X0, y = y)
  
  # require(reshape2)
  # require(ggplot2)
  # data <- data.frame(t(X0))
  # data$id <- 1:nrow(data)/nrow(data)
  # #reshape to long format
  # plot_data <- melt(data, id.var="id" )
  # gr <- ggplot(plot_data, aes(x=id, y=value, group=variable, colour=variable)) + geom_line(col = "gray", size = 1.2)
  # gr <- gr + theme_bw(base_size = 40) + theme(plot.margin = margin(t = 0,  r = 0,  b = 0, l = 0))  + labs(x = "t", y = "")
  # gr

  mse.mpen[f] <- mean( (alpha-fit.mpen$bh)^2 )  
  mse.fpcr[f] <- mean( (alpha-fit.fpcr$fhat)^2 )
  mse.ls[f]  <- mean( (alpha-fit.ls$bh)^2 )
  mse.smsp[f] <- mean( (alpha-fit.smsp$bh)^2 )
  mse.munp[f] <- mean((alpha - fit.munp$bh)^2)
  
  plot(grid, alpha, type = "l", lwd = 3, xlab = "t", ylab = "")
  lines(grid, fit.mpen$bh, lwd = 3, col = "blue")
  lines(grid, fit.fpcr$fhat, col = "gray", lwd = 3)
  lines(grid, fit.smsp$bh, lwd = 3, col = "blue")
  
  matr.mpen[, f] <- fit.mpen$bh
  matr.fpcr[, f] <- fit.fpcr$fhat
  matr.ls[, f] <- fit.ls$bh
  matr.smsp[, f] <- fit.smsp$bh
  matr.munp[, f] <- fit.munp$bh
}
mean(mse.mpen, na.rm = TRUE)*1000 ; median(mse.mpen, na.rm = TRUE)*1000
mean(mse.ls, na.rm = TRUE)*1000 ; median(mse.ls, na.rm = TRUE)*1000
mean(mse.fpcr, na.rm = TRUE)*1000 ; median(mse.fpcr, na.rm = TRUE)*1000
mean(mse.smsp, na.rm = TRUE)*1000 ; median(mse.smsp, na.rm = TRUE)*1000
mean(mse.munp, na.rm = TRUE)*1000 ; median(mse.munp, na.rm = TRUE)*1000

matplot(grid, matr1, lwd = 3, col = "gray", type = "l", cex.axis = 2, cex.lab = 2) ; lines(grid, alpha, lwd = 3, col = "black"); grid()
# matplot(grid, matr2, lwd = 3, col = "gray", type = "l", cex.axis = 2, cex.lab = 2) ; lines(grid, alpha, lwd = 3, col = "black"); grid()
matplot(grid, matr3, lwd = 3, col = "gray", type = "l", cex.axis = 2, cex.lab = 2) ; lines(grid, alpha, lwd = 3, col = "black"); grid()

##################################################################################################################################################
# Second experiment
Nrep <- 1000
n <- 150
p <- 100
grid <- seq(0, 1, length = p)
alpha = grid^2*dnorm(grid)

mse1 <- rep(NA, Nrep)
mse2 <- rep(NA, Nrep)
mse3 <- rep(NA, Nrep)

matr1 <- matrix(NA, nrow = p, ncol = Nrep)
matr2 <- matrix(NA, nrow = p, ncol = Nrep)
matr3 <- matrix(NA, nrow = p, ncol = Nrep)

for(f in 617:Nrep){
  message('Iter = ', f, ' of ', Nrep)
  X0 <- matrix(rnorm(n*p), nrow = n, ncol = p)
  for(i in 1:n){
    for(j in 1:49){
      X0[i, ] <- X0[i, ] + j^{-1}*rnorm(1, 0, 1)*sqrt(2)*sapply(grid, FUN= function(x) cos((j-1)*pi*x))
      #X0[i, ] <- X0[i, ] +  1*j^{-1}*rt(1, df = 5)*sqrt(2)*sapply(grid, FUN= function(x) cos((j-1)*pi*x)) # For leverage contamination
    }
  }
  y0 = X0%*%alpha
  y <- y0 + rnorm(n)
  # y <- y0 + rt(n, df = 3)
  # y <- y0 + rnormMix(n, mean1 = 0, sd = 1, mean2 = 14, sd2 = 1, p.mix = 0.1)
  fit1 <- m.pen.sp(x = X0, y = y, nbasis = round(min(n/4, 40)), n.se = 0)
  # fit2 <- fpcr(y, xfuncs = X0)
  fit3 <- ls.pen.sp(x = X0, y = y, nbasis = round(min(n/4, 40)), n.se = 0)
  
  mse1[f] <- mean( (alpha-fit1$bh)^2 )  
  # mse2[f] <- mean( (alpha-fit2$fhat)^2 )
  mse3[f]  <- mean( (alpha-fit3$bh)^2 ) 
  
  plot(grid, alpha, type = "l", lwd = 3, xlab = "t", ylab = "")
  lines(grid, fit1$bh, lwd = 3, col = "red")
  # lines(grid, fit2$fhat, col = "blue", lwd = 3)
  lines(grid, fit3$bh, lwd = 3, col = "blue")
  
  matr1[, f] <- fit1$bh
  # matr2[, f] <- fit2$fhat
  matr3[, f] <- fit3$bh
}
mean(mse1, na.rm = TRUE)*1000
mean(mse2, na.rm =TRUE)*1000
mean(mse3, na.rm =TRUE)*1000

median(mse1, na.rm = TRUE)*1000
median(mse2, na.rm = TRUE)*1000
median(mse3, na.rm = TRUE)*1000

matplot(grid, matr1, lwd = 3, col = "gray", type = "l", cex.axis = 2, cex.lab = 2) ; lines(grid, alpha, lwd = 3, col = "black"); grid()
# matplot(grid, matr2, lwd = 3, col = "gray", type = "l", cex.axis = 2, cex.lab = 2) ; lines(grid, alpha, lwd = 3, col = "black"); grid()
matplot(grid, matr3, lwd = 3, col = "gray", type = "l", cex.axis = 2, cex.lab = 2) ; lines(grid, alpha, lwd = 3, col = "black"); grid()

##################################################################################################################################################
# Third experiment
Nrep <- 1000
n <- 150
p <- 100
grid <- seq(0, 1, length = p)
alpha = 1/(1+exp(-20*(grid-0.5)))

mse1 <- rep(NA, Nrep)
mse2 <- rep(NA, Nrep)

matr1 <- matrix(NA, nrow = p, ncol = Nrep)
matr2 <- matrix(NA, nrow = p, ncol = Nrep)

for(f in 1:Nrep){
  message('Iter = ', f, ' of ', Nrep)
  X0 <- matrix(rnorm(n*p), nrow = n, ncol = p)
  for(i in 1:n){
    for(j in 1:49){
      X0[i, ] <- X0[i, ] + j^{-1}*rnorm(1, 0, 1)*sqrt(2)*sapply(grid, FUN= function(x) cos((j-1)*pi*x))
      #X0[i, ] <- X0[i, ] +  1*j^{-1}*rt(1, df = 5)*sqrt(2)*sapply(grid, FUN= function(x) cos((j-1)*pi*x)) # For leverage contamination
    }
  }
  y0 = X0%*%alpha
  y <- y0 + rnorm(n)
  # y <- y0 + rt(n, df = 3)
  # y <- y0 + rnormMix(n, mean1 = 0, sd = 1, mean2 = 14, sd2 = 1, p.mix = 0.1)
  fit1 <- m.pen.sp(x = X0, y = y, nbasis = round(min(n/4, 40)))
  fit2 <- fpcr(y, xfuncs = X0)
  
  mse1[f] <- mean( (alpha-fit1$bh/p)^2 )  
  mse2[f] <- mean( (alpha-fit2$fhat)^2 )
  
  plot(grid, alpha, type = "l", lwd = 3, xlab = "t", ylab = "")
  lines(grid, fit1$bh/p, lwd = 3, col = "red")
  lines(grid, fit2$fhat, col = "blue", lwd = 3)
  
  matr1[, f] <- fit1$bh/p
  matr2[, f] <- fit2$fhat
}
mean(mse1, na.rm = TRUE)*Nrep
mean(mse2, na.rm =TRUE)*Nrep

median(mse1, na.rm = TRUE)*Nrep
median(mse2, na.rm = TRUE)*Nrep

matplot(grid, matr1, lwd = 3, col = "gray", type = "l", cex.axis = 2, cex.lab = 2) ; lines(grid, alpha, lwd = 3, col = "black"); grid()
matplot(grid, matr2, lwd = 3, col = "gray", type = "l", cex.axis = 2, cex.lab = 2) ; lines(grid, alpha, lwd = 3, col = "black"); grid()


##################################################################################################################################################
# Fourth experiment
Nrep <- 1000
n <- 150
p <- 100
grid <- seq(0, 1, length = p)
alpha = -dnorm(grid, mean=.2, sd=.03) + 3*dnorm(grid, mean=.5, sd=.03) + dnorm(grid, mean=.75, sd=.04)

mse1 <- rep(NA, Nrep)
mse2 <- rep(NA, Nrep)

matr1 <- matrix(NA, nrow = p, ncol = Nrep)
matr2 <- matrix(NA, nrow = p, ncol = Nrep)

for(f in 113:Nrep){
  message('Iter = ', f, ' of ', Nrep)
  X0 <- sqrt(2)*matrix(rnorm(n*p), nrow = n, ncol = p)
  for(i in 1:n){
    for(j in 2:50){
      X0[i, ] <- X0[i, ] + j^{-1}*rnorm(1, 0, 1)*sqrt(2)*sapply(grid, FUN= function(x) cos((j-1)*pi*x))
      # X0[i, ] <- X0[i, ] +  1*j^{-1}*rt(1, df = 3)*sqrt(2)*sapply(grid, FUN= function(x) cos((j-1)*pi*x)) # For leverage contamination
    }
  }
  y0 = X0%*%alpha
  y <- y0 + rnorm(n)
  # y <- y0 + rt(n, df = 3)
  # y <- y0 + rnormMix(n, mean1 = 0, sd = 1, mean2 = 14, sd2 = 1, p.mix = 0.1)
  fit1 <- m.pen.sp(x = X0, y = y, nbasis = round(min(n/4, 40)))
  fit2 <- m.sp(x = X0, y = y)
  # fit2 <- fpcr(y, xfuncs = X0)
  
  mse1[f] <- mean( (alpha-fit1$bh)^2 )  
  # mse2[f] <- mean( (alpha-fit2$fhat)^2 )
  mse2[f] <- mean( (alpha-fit2$bh)^2 )  
  
  plot(grid, alpha, type = "l", lwd = 3, xlab = "t", ylab = "")
  lines(grid, fit1$bh, lwd = 3, col = "red")
  # lines(grid, fit2$fhat, col = "blue", lwd = 3)
  lines(grid, fit2$bh, lwd = 3, col = "blue")
  
  matr1[, f] <- fit1$bh
  matr2[, f] <- fit2$bh
}
mean(mse1, na.rm = TRUE)*1000 ; median(mse1, na.rm = TRUE)*1000
mean(mse2, na.rm =TRUE)*1000 ;median(mse2, na.rm = TRUE)*1000

matplot(grid, matr1, lwd = 3, col = "gray", type = "l", cex.axis = 2, cex.lab = 2) ; lines(grid, alpha, lwd = 3, col = "black"); grid()
matplot(grid, matr2, lwd = 3, col = "gray", type = "l", cex.axis = 2, cex.lab = 2) ; lines(grid, alpha, lwd = 3, col = "black"); grid()

library(reshape2)
library(ggplot2)
data2 <- data.frame(id = grid, alpha)
data <- data.frame(matr1)
data$id <- 1:nrow(data)/nrow(data)
#reshape to long format
plot_data <- melt(data, id.var="id" )
#plot
gr <- ggplot(plot_data, aes(x=id, y=value, group=variable, colour=variable)) + geom_line(col = "gray", size = 2) + ylim(-15, 43)
gr <- gr + theme_bw(base_size = 40) + theme(plot.margin = margin(t = 0,  r = 0,  b = 0, l = 0))  + labs(x = "", y = "")
# f0f <- function(x) -dnorm(x, mean=.2, sd=.03) + 3*dnorm(x, mean=.5, sd=.03) + dnorm(x, mean=.75, sd=.04)
gr <- gr + geom_line(data = data2, aes(x = id, y = alpha, group = 1), size = 1.2, inherit.aes = FALSE)
gr

data2 <- data.frame(id = grid, alpha)
data <- data.frame(matr2)
data$id <- 1:nrow(data)/nrow(data)
#reshape to long format
plot_data <- melt(data, id.var="id" )
#plot
gr <- ggplot(plot_data, aes(x=id, y=value, group=variable, colour=variable)) + geom_line(col = "gray", size = 2) + ylim(-15, 43)
gr <- gr + theme_bw(base_size = 40) + theme(plot.margin = margin(t = 0,  r = 0,  b = 0, l = 0))  + labs(x = "", y = "")
# f0f <- function(x) -dnorm(x, mean=.2, sd=.03) + 3*dnorm(x, mean=.5, sd=.03) + dnorm(x, mean=.75, sd=.04)
gr <- gr + geom_line(data = data2, aes(x = id, y = alpha, group = 1), size = 1.2, inherit.aes = FALSE)
gr


