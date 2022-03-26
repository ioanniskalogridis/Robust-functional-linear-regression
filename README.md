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

``` norder ``` : the order of the spline basis to be used. By default it is equal to 4 corresponding to cubic B-splines.

``` nbasis ``` : the number of B-spline basis functions with a given order. By default, ```nbasis = [max(40, n/4)]``` with ```n``` the sample size.

``` k ``` : the tuning parameter for the bisquare loss function. By default ```k = 4.685``` corresponding to ```95%``` efficiency in the Gaussian location model.

``` q ``` : the order of the penalty. By default ``` q = 2 ``` corresponding to a penalty on the integrated squared second derivative.

The function ```ls.pen.sp()``` takes similar arguments except that ```k``` in this case is not needed. For this example, we need to load the data from the fda.usc package in R.

```
require(fda.usc)
data(tecator)
```

Then, we extract the functional predictor and the response variable. We may plot the predictor curves and the responses variable in order to visualize the data.
```
x <- tecator$absorp.fdata
x <- x$data
y <- tecator$y
y <- y[, 1]
par(mfrow = c(1,2))
matplot(t(x), type = "l", col = "gray", lty = 1, lwd = 3, xlab ="t")
hist(y, main = "")
```
![Fig1](https://user-images.githubusercontent.com/77626768/160242630-9457bb8a-0a4f-48de-8133-22e1bc29875a.png)

We fit the penalized MM-estimator and the least-squares estimator on this dataset.

```
fit.mm <-  mm.pen.sp(x = x, y = y)
fit.ls <- ls.pen.sp(x, y)
```
Here are plots of the estimated coefficient functions.

```
plot(fit.mm$bh, type = "l", col = "blue", lwd = 3, xlab = "t", ylab = "Estimates")
lines(fit.ls$bh, type = "l", col = "red", lwd = 3, lty = 3); grid()
legend("bottomleft", legend = c("PMM", "PLS"), col = c("blue", "red"), lwd = c(3, 3), lty = c(1, 3), cex = 1.5)
```
![Coefs](https://user-images.githubusercontent.com/77626768/160139488-dcf70cd1-a2b0-4998-81c1-3d02f1034045.png)


The least-squares estimates are considerably more wiggly than the MM-estimates. And here are normal QQ plots of the residuals.

```

qqnorm(fit.mm$resids, pch = 20, cex = 1.5) ; qqline(fit.mm$resids, lwd  = 3, col = "red")

```
![Qmm](https://user-images.githubusercontent.com/77626768/160139893-f84b7603-111a-42aa-b832-4e766d816700.png)


```

qqnorm(fit.ls$resids, pch = 20, cex = 1.5) ; qqline(fit.ls$resids, lwd  = 3, col = "red")

```
![Qls](https://user-images.githubusercontent.com/77626768/160139987-a149d15b-58a9-43ee-bf7b-76cb8ddd3bd4.png)

Notice that the MM-estimator is successful in detecting a large number of outliers whereas the residuals of the least-squares estimator suffer from the masking effect.
