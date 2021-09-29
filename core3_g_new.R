################ Utility functions #############

standard <- function(x) {
  return((x - mean(x)) / sd(x))
}



transform.sob <- function(X) {
  Xlim <- apply(X, 2, range)
  Xstd <- matrix(nr = nrow(X), nc = ncol(X))
  for (i in (1:ncol(X))) {
    Xstd[, i] <- (X[, i] - Xlim[1, i]) / diff(Xlim[, i])
  }
  return(list(Xstd = Xstd, Xlim = Xlim))
}


fctn.ksm <- function(Vnew, Yadj, V, h1) {
  K1 <- getG_grid(V, Vnew, h1)
  out = apply(K1, 2, function(x) sum(x * Yadj) / sum(x))
  return(out)
}



######################################### CATE #####################



###### Vg: a grid of V ###########
ATE.ncb.core2 <- function(V, Vg, h1, ind, K, lam1s, lam2s = 0.01 * lam1s, lower = 1,
  upper = Inf, thresh.ratio = 1e-09, traceit = TRUE, w0 = NULL, maxit = 2000, xtol_rel = 1e-10,
  full = FALSE, algorithm = "lbfgsb3c") {



  N <- length(ind)

  # construct a vector
  tind <- as.logical(ind)
  n1 <- sum(tind)
  n2 <- sum(!tind)

  # lower and upper bound
  if (length(lower) == 1)
    lower <- rep(lower, n1)
  if (length(upper) == 1)
    upper <- rep(upper, n1)

  # construct P1 and q1
  e <- eigen(K)
  thresh <- e$values[1] * thresh.ratio
  eind <- (e$values >= thresh)
  r <- sum(eind)

  # construct L2 Kernal h1 matrix in the middle construct the penalty term weight
  # GG <- getG2(V, lower = lowv, upper = upv, sigma=h1)
  ForG <- getG_grid(V, Vg, h1)
  fV <- colMeans(ForG)
  VnwG <- ForG / (matrix(rep(fV, each = nrow(ForG)), nrow = nrow(ForG), ncol = ncol(ForG)))
  GG = tcrossprod(VnwG) / length(Vg) # * (max(Vg) - min(Vg))
  Vnw_w <- diag(GG)



  if (traceit)
    cat("Number of components used: r = ", r, "with threshod.ratio =", thresh.ratio,
      "\n")

  if (is.null(w0)) {
    w0 <- lower + 1
  }

  nlams <- length(lam1s)
  lam1s <- sort(lam1s) # lams is sorted
  outlist <- list()
  ws <- matrix(0, nr = nlams, nc = N)
  SNs <- array(dim = nlams)
  unorms <- array(dim = nlams)
  fittedus <- array(dim = c(nlams, nc = N))
  rmseTo1s <- array(dim = nlams)
  outlist <- NULL

  alg <- rep(1, nlams)
  obj1 <- array(dim = nlams)


  if (length(lam2s) == 1) {
    lam2s <- rep(lam2s, nlams)
  }
  status = NULL
  for (i in (1:nlams)) {
    nP1 <- e$vectors[, eind] / sqrt(N) # P1/sqrt(N)
    nq1 <- -lam1s[i] * N / e$values[eind] # -lam1s[i]*N1*q1
    # reorder
    oo <- order(nq1)
    nq1 <- nq1[oo]
    nP1 <- nP1[, oo]

    if (traceit)
      cat("####", i, ":", "\tlam1 =", lam1s[i], "\tlam2 =", lam2s[i], "\n")
    if (algorithm == "lbfgsb3c") {
      res <- lbfgsb3c::lbfgsb3c(par = w0, fn = eval_obj_cpp, gr = eval_grad_cpp,
        lower = lower, upper = upper, N = N, ind = 1 * tind, nP1 = nP1, nq1 = nq1,
        lam2 = lam2s[i], GG = GG, Vnw_w = Vnw_w, control = list(reltol = xtol_rel,
          maxit = maxit))
      res$par[res$par < 1] = 1
      #res2 <- lbfgsb3c::lbfgsb3c(par = w0, fn = eval_obj_cpp, lower = lower,
      #  upper = upper, N = N, ind = 1 * tind, nP1 = nP1, nq1 = nq1, lam2 = lam2s[i],
      #  GG = GG, Vnw_w = Vnw_w, control = list(reltol = xtol_rel, maxit = maxit))
      #reso <- optim(par = w0, fn = eval_obj_cpp, gr = eval_grad_cpp, method = 'L-BFGS-B',
      #   lower = lower, upper = upper, N = N, ind = 1 * tind, nP1 = nP1, nq1 = nq1,
      #  lam2 = lam2s[i], GG = GG, Vnw_w = Vnw_w, control = list(factr = xtol_rel,
      #    maxit = maxit))

      # cat(res$message) cat(res$counts) cat('\n')
      if (is.null(res$convergence)) {
        status <- c(status, res$convergence)
      } else {
        status <- c(status, "ERROR")
      }
      obj1[i] <- res$value
      ws[i, tind] <- res$par
    } else {
      res <- nloptr::nloptr(x0 = w0, eval_f = eval_obj_grad2_cpp, lb = lower,
        ub = upper, N = N, ind = 1 * tind, nP1 = nP1, nq1 = nq1, lam2 = lam2s[i],
        GG = GG, Vnw_w = Vnw_w, opts = list(algorithm = "NLOPT_LD_LBFGS",
          xtol_rel = xtol_rel, print_level = 0, maxeval = maxit, check_derivatives = F))
      #cat(res$message)
      #cat(res$iterations)
      #cat("\n")
      status = c(status, res$status)
      obj1[i] = res$objective
      ws[i, tind] <- res$sol
    }



    # res <- optim(par = w0, fn = eval_obj_cpp, gr = eval_grad_cpp, lower = lower,
    # upper = upper, method = 'L-BFGS-B', N=N, ind=1*tind, nP1=nP1, nq1=nq1,
    # lam2=lam2s[i],GG=GG, Vnw_w = Vnw_w, control = list(factr = xtol_rel, maxit =
    # maxit)) res <- nloptr::nloptr(x0=w0, eval_f=eval_obj_grad2_cpp, lb=lower,
    # ub=upper, N=N, ind=1*tind, nP1=nP1, nq1=nq1, lam2=lam2s[i],GG=GG, Vnw_w =
    # Vnw_w, opts=list(algorithm='NLOPT_LD_LBFGS', xtol_rel=xtol_rel, print_level=0,
    # maxeval=maxit, check_derivatives=F))

    # print(res$iterations)

    temp <- t(nP1) %*% (ws[i,] - 1)
    ee <- eigen(temp %*% t(temp) + diag(nq1))

    # compute SN
    SNs[i] <- (sum(temp * ee$vectors[, 1])) ^ 2
    unorms[i] <- sqrt(sum(nq1 / (-lam1s[i]) * ee$vectors[, 1] ^ 2))
    fittedus[i,] <- as.vector(nP1 %*% ee$vectors[, 1] * (N))
    rmseTo1s[i] <- sqrt(mean((ws[i,] - 1) ^ 2))

    w0 <- ws[i, tind]

  }

  return(list(outlist = outlist, ws = ws, SNs = SNs, unorms = unorms, fittedus = fittedus,
    rmseTo1s = rmseTo1s, lam1s = lam1s, lam2s = lam2s, alg = alg, obj1 = obj1))
}


eval_obj_grad2 <- function(w, N, r, tind, nP1, nq1, lam2, GG, Vnw_w) {
  z <- rep(-1, N)
  z[tind] <- w - 1
  V <- sum(Vnw_w[tind] * w ^ 2) / N
  Mtemp <- t(nP1) %*% diag(z) %*% GG %*% diag(z) %*% nP1
  eigen.out <- eigen(Mtemp + diag(nq1), symmetric = TRUE)
  v <- eigen.out$vectors[, 1]
  length(v)
  tZnP1 <- t(diag(z) %*% nP1)
  return(list(objective = eigen.out$values[1] + lam2 * V, gradient = as.vector(2 *
    as.vector(t(v) %*% tZnP1[, tind]) * ((GG %*% diag(z) %*% nP1)[tind,] %*%
    v) + 2 * lam2 * Vnw_w[tind] * (w) / N))) # modified
}


eval_obj <- function(w, N, tind, nP1, nq1, lam2, GG, Vnw_w) {
  z <- rep(-1, N)
  z[tind] <- w - 1
  V <- sum(Vnw_w[tind] * w ^ 2) / N
  Mtemp <- t(nP1) %*% diag(z) %*% GG %*% diag(z) %*% nP1
  eigen.out <- eigen(Mtemp + diag(nq1), symmetric = TRUE)
  return(eigen.out$values[1] + lam2 * V)
}


eval_grad <- function(w, N, tind, nP1, nq1, lam2, GG, Vnw_w) {
  z <- rep(-1, N)
  z[tind] <- w - 1
  V <- sum(Vnw_w[tind] * w ^ 2) / N
  Mtemp <- t(nP1) %*% diag(z) %*% GG %*% diag(z) %*% nP1
  eigen.out <- eigen(Mtemp + diag(nq1), symmetric = TRUE)
  v <- eigen.out$vectors[, 1]
  tZnP1 <- t(diag(z) %*% nP1)
  return(as.vector(2 * as.vector(t(v) %*% tZnP1[, tind]) * ((GG %*% diag(z) %*%
    nP1)[tind,] %*% v) + 2 * lam2 * Vnw_w[tind] * (w) / N))
}


ATE.ncb.SN2 <- function(V, Vg, h1, ind, K, lam1s, lam2s = 0.01 * lam1s, lower = 1,
  upper = Inf, thresh.ratio = 1e-08, traceit = TRUE, w0 = NULL, maxit = 2000, xtol_rel = 1e-08,
  method = 1, check = FALSE, full = FALSE, algorithm = "lbfgsb3c") {

  ores <- ATE.ncb.core2(V = V, Vg = Vg, h1 = h1, ind = ind, K = K, lam1s, lam2s = lam2s,
    lower = lower, upper = upper, thresh.ratio = thresh.ratio, traceit = traceit,
    w0 = w0, maxit = maxit, xtol_rel = xtol_rel, full = full, algorithm = algorithm)

  # find which SN is the smallest
  if (length(lam1s) == 1) {
    ind = 1
  } else {
    if (method == 1) {
      ind <- which.min(ores$SNs)
    } else if (method == 2) {
      ind <- which((diff(ores$SNs) / diff(ores$lam1s) > (-1e-06)) & (diff(ores$SNs <= 0)))[1]
    }
  }
  if (is.na(ind)) {
    ind <- which.min(ores$SNs)
  }
  return(list(w = ores$ws[ind,], ind = ind, warns = c(ind == 1, ind == length(lam1s)),
    ores = ores))
}
