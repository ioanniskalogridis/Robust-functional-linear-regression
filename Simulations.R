require(refund)
require(fda)
require(MASS)
require(EnvStats)

Nrep <- 1000
n <- 150
p <- 100
grid <- seq(0, 1, length = p)

alpha <- sin(2*pi*grid)
# alpha = grid^2*dnorm(grid, 0, 0.1)
# alpha = 1/(1+exp(-20*(grid-0.5)))
# alpha = -dnorm(grid, mean=.2, sd=.03) + 3*dnorm(grid, mean=.5, sd=.03) + dnorm(grid, mean=.75, sd=.04)

mse.ls <- rep(NA, Nrep)
mse.mpen <- rep(NA, Nrep)
mse.munp <- rep(NA, Nrep)
mse.smsp <- rep(NA, Nrep)
mse.rkhs <- rep(NA, Nrep)

matr.ls <- matrix(NA, nrow = p, ncol = Nrep)
matr.mpen <- matrix(NA, nrow = p, ncol = Nrep)
matr.munp <- matrix(NA, nrow = p, ncol = Nrep)
matr.smsp <- matrix(NA, nrow = p, ncol = Nrep)
matr.rkhs <- matrix(NA, nrow = p, ncol = Nrep)

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
  y <- y0 + rnorm(n)
  # y <- y0 + rt(n, df = 3)
  # y <- y0 + rnormMix(n, mean1 = 0, sd = 1, mean2 = 14, sd2 = 1, p.mix = 0.1)
  
  fit.mpen <- mm.pen.sp(x = X0, y = y, nbasis = round(min(n/4, 40)), n.se = 0)
  fit.ls <- ls.pen.sp(x = X0, y = y, nbasis = round(min(n/4, 40)), n.se = 0)
  fit.smsp <- m.sm.sp(x = X0, y = y, t = grid)
  fit.munp <- m.sp(x = X0, y = y)
  fit.rkhs <- flm.rkhs.rob.t(X0, y = y, dom = grid)
  
  # The following code produces the plots in Figure 1 of the paper
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
  mse.rkhs[f]<- mean((alpha - fit.rkhs$beta/p)^2)

  # plot(grid, alpha, type = "l", lwd = 3, xlab = "t", ylab = "")
  # lines(grid, fit.mpen$bh, lwd = 3, col = "blue")
  # lines(grid, fit.ls$bh, col = "gray", lwd = 3)
  # lines(grid, fit.smsp$bh, lwd = 3, col = "blue")
  # lines(grid, fit.rkhs$beta/p, lwd = 3, col = "red")
  
  matr.mpen[, f] <- fit.mpen$bh
  matr.ls[, f] <- fit.ls$bh
  matr.smsp[, f] <- fit.smsp$bh
  matr.munp[, f] <- fit.munp$bh
  matr.rkhs[, f] <- fit.rkhs$beta/p
}
mean(mse.mpen, na.rm = TRUE)*1000 ; median(mse.mpen, na.rm = TRUE)*1000
mean(mse.ls, na.rm = TRUE)*1000 ; median(mse.ls, na.rm = TRUE)*1000
mean(mse.smsp, na.rm = TRUE)*1000 ; median(mse.smsp, na.rm = TRUE)*1000
mean(mse.munp, na.rm = TRUE)*1000 ; median(mse.munp, na.rm = TRUE)*1000
mean(mse.rkhs, na.rm = TRUE)*1000 ; median(mse.rkhs, na.rm = TRUE)*1000

matplot(grid, matr.mpen, lwd = 3, col = "gray", type = "l", cex.axis = 2, cex.lab = 2) ; lines(grid, alpha, lwd = 3, col = "black"); grid()
matplot(grid, matr.smsp, lwd = 3, col = "gray", type = "l", cex.axis = 2, cex.lab = 2) ; lines(grid, alpha, lwd = 3, col = "black"); grid()
matplot(grid, matr.munp, lwd = 3, col = "gray", type = "l", cex.axis = 2, cex.lab = 2) ; lines(grid, alpha, lwd = 3, col = "black"); grid()
matplot(grid, matr.ls, lwd = 3, col = "gray", type = "l", cex.axis = 2, cex.lab = 2) ; lines(grid, alpha, lwd = 3, col = "black"); grid()
matplot(grid, matr.rkhs, lwd = 3, col = "gray", type = "l", cex.axis = 2, cex.lab = 2) ; lines(grid, alpha, lwd = 3, col = "black"); grid()

