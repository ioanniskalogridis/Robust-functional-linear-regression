# Robust-functional-linear-regression
Robust penalized estimators for the functional linear model. 
The main function m.pen.sp is in the directory named "Main function (MM P-splines).R".
This function implements the MM P-spline estimator used in Kalogridis and Van Aelst (2021).
The implementation uses an initial unpenalized S-estimator and then the penalized variant of IRLS in order to obtain the MM-estimator.
Contact me at ioannis.kalogridis@kuleuven.be for remarks and/or suggestions.



# Tecator data

Here is an example in R. First, load the functions mm.pen.sp and ls.pen.sp from the files "Main function (MM P-splines).R" and "Least-squares P-spline estimator.R", respectively. The arguments of ```mm.pen.sp()``` are

``` x ``` : the matrix of the discretized functional predictors.

``` y ``` : the response vector.

``` norder ``` : the order of the spline basis to be used. By default equal to 4 corresponding to cubic B-splines.

``` nbasis ``` : the number of B-spline basis functions with a given order. By default, ```nbasis = [max(40, n/4)]``` with ```n``` the sample size.

``` k ``` : the tuning parameter for the bisquare loss function. By default ```k = 4.685``` corresponding to ```95%``` efficiency in the Gaussian location model.

``` q ``` : the order of the penalty. By default ``` q = 2 ``` corresponding to a penalty on the integrated squared second derivative.

The ```ls.pen.sp()``` takes similar arguments except that ```k``` in this case is not needed.

```
require(fda.usc)
data(tecator)
x <- tecator$absorp.fdata
x <- x$data
y <- tecator$y
y <- y[, 1]
matplot(t(x), type = "l", col = "gray", lty = 1, lwd = 3)

fit.mm <-  mm.pen.sp(x = x, y = y)
fit.ls <- ls.pen.sp(x, y)

plot(fit.mm$bh, type = "l", col = "blue", lwd = 3)
lines(fit.ls$bh, type = "l", col = "red", lwd = 3, lty = 3)
hist(fit.mm$resids/fit.mm$scale)

qqnorm(fit.mm$resids, pch = 20, cex = 1.5) ; qqline(fit.mm$resids, lwd  = 3, col = "red")
qqnorm(fit.ls$resids, pch = 20, cex = 1.5) ; qqline(fit.ls$resids, lwd  = 3, col = "red")

```

Notice that the outliers are much less visible from the least-squares residuals.
