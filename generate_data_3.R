
p = 4


hseq = exp(seq(log(0.2 * n ^ (-1 / 5)), log(20 * n ^ (-1 / 5)), length.out = 10))
Z1 <- runif(n, -2, 2)
Z2 <- runif(n, -2, 2)
Z3 <- runif(n, -2, 2)
Z4 <- runif(n, -2, 2)

X1 <- Z1
X2 <- Z1 ^ 2 + Z2
X3 <- exp(Z3 / 2) + Z2
X4 <- sin(2 * Z1) + Z4


X = cbind(X1, X2, X3, X4)


Y1 <- 10 + Z1 ^ 2 + 2 * Z1 * sin(2 * Z1) + Z2 ^ 2 + sin(2 * Z3) * Z4 ^ 2 + rnorm(n, sd = sqrt(variance))
Y0 <- 10 - Z1 ^ 2 - 2 * Z1 * sin(2 * Z1) + Z2 ^ 2 + sin(2 * Z3) * Z4 ^ 2 + rnorm(n, sd = sqrt(variance))


prop <- 1 / (1 + exp(-X1 + X3))



treat <- rbinom(n, 1, prop)
Y <- Y1 * treat + Y0 * (1 - treat)


Vnew <- seq(from = -1.8, to = 1.8, length.out = 300)
fctn.CATE <- function(x) { 2 * x ^ 2 + 4 * x * sin(2 * x) }
real <- fctn.CATE(Vnew)

Vnew = as.matrix(Vnew)
V = X[, 1]
V = as.matrix(V)

