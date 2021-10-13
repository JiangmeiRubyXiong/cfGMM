# cfGMM
This is the function that produces estimator for gamma mixture model. This function uses closed form estimator, so it is faster than state-of-the-art estimation. Install by running in your R console:

```{r}
devtools::install_github("JiangmeiRubyXiong/cfGMM")
```

A quick example (also can be found in `help(cfGMM)`:)
```{r}
set.seed(1)

## construct dataset
phi <- c(0.3, 0.7) # mixing proportion
a <- c(0.5, 8) #shape
b <- c(1/2, 1/3) #scale
n <- 10000 #data size
ind <- sample(c(1,2) ,size=n, replace = TRUE, prob = phi)
data.gamma <- c(rgamma(sum(ind==1), shape=a[1], rate=b[1]),
                 rgamma(sum(ind==2), shape=a[2], rate=b[2]))
                 
## run the model
out <- cfGMM(data.gamma, k=2)
out[["gamma.pars"]] # output shape and scale paramter
out[["lambda"]] # output lambda paramter
```
