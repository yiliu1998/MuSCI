#### A simulated data example
#### training data
set.seed(10012)
N <- 500
s <- 5

sT <- sort(rep(1:s, N))

X1 <- runif(N*s, -3, 3)
X2 <- rnorm(N*s, mean=1*sT, sd=1)
X3 <- rchisq(N*s, df=10, ncp=sT-1)
X4 <- X1*X2
X5 <- X1*X3
X6 <- X2^2

beta <- c(1,1,1,.1,.1,.1)
Y <- cbind(X1,X2,X3,X4,X5,X6) %*% as.matrix(beta, ncol=1) + rnorm(N*s)

e.x <- 1/(1+exp(0.5*X1^2-0.1*X2-0.2*X3))
R <- rbinom(N*s, size=1, prob=e.x)
mean(R) # non-missingness rate
Y[R==0] <- NA

data.train <- data.frame(cbind(X1,X2,X3,X4,X5,X6), Y=Y, sT=sT)

#### testing data
N.test <- 1000
X1 <- runif(N.test, -3, 3)
X2 <- rnorm(N.test, mean=1, sd=1)
X3 <- rchisq(N.test, df=10, ncp=0)
X4 <- X1*X2
X5 <- X1*X3
X6 <- X2^2

beta <- c(1,1,1,.1,.1,.1)
Y <- cbind(X1,X2,X3,X4,X5,X6) %*% as.matrix(beta, ncol=1) + rnorm(N.test)
e.x <- 1/(1+exp(0.5*X1^2-0.1*X2-0.2*X3))
R <- rbinom(N.test, size=1, prob=e.x)
mean(R) 
Y[R==0] <- NA

data.test <- data.frame(cbind(X1,X2,X3,X4,X5,X6), Y=Y)

save(file="example.Rdata", data.test, data.train)
