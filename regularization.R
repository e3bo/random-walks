#' Run a general purpose optimizer with lasso or elasticnet regularization
#'
#' @param x matrix of predictors for linear model
#' @param y response data that will be used to fit the model
#' @param calc_convex_nll function that will be used to calculate the likelihood of the data
#' @param param_map function that maps vector of parameters into named list for likelihood function
#' @param alpha proportion of penalty function that is L1 instead of L2
#' @param nlambda number of penalty values to on regularization path
#' @param lambda.min.ratio ratio of lowest penalty on path to highest
#' penalty. Used to determine sequence of penalties if lambda argument
#' is null
#' @param lambda sequence of penalties to include in regularization path
#' @param standardize Not currently used
#' @param thresh threshold for subgradient to determine convergence at each lambda
#' @param dfmax not currently used
#' @param pmax not currently used
#' @param exclude not currently used
#' @param penalty.factor vector equal in length to the number of
#' parameters, factor by which lambda is multiplied for each parameter to calculate penalty
#' @param lower.limits lower limits for penalized parameters
#' @param upper.limits upper limits for penalized parameters
#' @param maxit maximum number of iterations in optimization attempt for each lambda
#' @param make_log whether or not to print out information about the optimization to temporary files
#' @param winit parameter values used to initialize optimization
#'
#' @export
get_gpnet <- function(x, y, calc_convex_nll, param_map, alpha=1, nlambda=100,
                      lambda.min.ratio=0.01, lambda=NULL, standardize=TRUE,
                      thresh=1e-4, dfmax=nvars + 1,
                      pmax=min(dfmax*2 + 20,nvars), exclude,
                      penalty.factor=rep(1, nvars), lower.limits=-Inf,
                      upper.limits=Inf, maxit=100, make_log = FALSE,
                      winit){
  if (alpha > 1) {
    warning("alpha >1; set to 1")
    alpha <- 1
  }
  if (alpha < 0) {
    warning("alpha<0; set to 0")
    alpha <- 0
  }
  alpha <- as.double(alpha)
  this.call <- match.call()
  nlam <- as.integer(nlambda)
  np <- dim(x)
  if (is.null(np) | (np[2] < 1))
    stop("x should be a matrix with 1 or more columns")
  nrates <- as.integer(np[1])
  #nvars <- as.integer(np[2])
  nvars <- as.integer(sum(penalty.factor > .Machine$double.eps))
  vnames <- colnames(x)
  if (is.null(vnames))
    vnames <- paste("V", seq(nvars), sep = "")
  ne <- as.integer(dfmax)
  nx <- as.integer(pmax)
  if (!missing(exclude)) {
    jd <- match(exclude, seq(nvars), 0)
    if (!all(jd > 0))
      stop("Some excluded variables out of range")
    jd <- as.integer(c(length(jd), jd))
  }
  else jd <- as.integer(0)
  vp <- as.double(penalty.factor)
  if (any(lower.limits > 0)) {
    stop("Lower limits should be non-positive")
  }
  if (any(upper.limits < 0)) {
    stop("Upper limits should be non-negative")
  }
  if (length(lower.limits) < nvars) {
    if (length(lower.limits) == 1)
      lower.limits <- rep(lower.limits, nvars)
    else stop("Require length 1 or nvars lower.limits")
  }
  else lower.limits <- lower.limits[seq(nvars)]
  if (length(upper.limits) < nvars) {
    if (length(upper.limits) == 1)
      upper.limits = rep(upper.limits, nvars)
    else stop("Require length 1 or nvars upper.limits")
  }
  else upper.limits <- upper.limits[seq(nvars)]
  cl = rbind(lower.limits, upper.limits)
  if (any(abs(cl) < thresh)) {
    stop("Cannot enforce limits this close to zero")
  }
  storage.mode(cl) <- "double"
  isd <- as.integer(standardize)
  intr <- as.integer(sum(penalty.factor < .Machine$double.eps))
  thresh <- as.double(thresh)
  if (is.null(lambda)) {
    if (lambda.min.ratio >= 1)
      stop("lambda.min.ratio should be less than 1")
    flmin <- as.double(lambda.min.ratio)
    ulam <- double(1)
  }
  else {
    flmin <- as.double(1)
    if (any(lambda < 0))
      stop("lambdas should be non-negative")
    ulam <- as.double(rev(sort(lambda)))
    nlam <- as.integer(length(lambda))
  }
  fit <- gpnet(x, y, calc_convex_nll, param_map, alpha, nobs=NULL, nvars, jd, vp,
               cl, ne, nx, nlam, flmin, ulam, thresh, isd, intr, vnames,
               maxit, make_log = make_log, winit = winit)
  fit$call <- this.call
  fit$nrates <- nrates
  class(fit) <- c(class(fit), "gpnet")
  fit
}

gpnet <- function(x, y, calc_convex_nll, param_map, alpha, nobs, nvars, jd, vp,
                  cl, ne, nx, nlam, flmin, ulam, thresh, isd, intr, vnames,
                  maxit, a=0.1, r=0.01, relStart=0.0, mubar=1, beta=0.1,
                  make_log = FALSE, debug = TRUE, initFactor=10, winit){
  maxit <- as.integer(maxit)
  niter <- 0
  dim <- nvars + intr
  parInds <- seq(dim)
  I <- diag(nrow=dim)
  nll <- function(w){
    calc_convex_nll(w=w, x=x, y=y, pm=param_map)
  }
  nll_no_penalty <- function(w_nopen){
    w <- winit
    pen_ind <- vp > .Machine$double.eps
    w[!pen_ind] <- w_nopen
    calc_convex_nll(w=w, x=x, y=y, pm=param_map)
  }
  is_unpenalized <- vp < .Machine$double.eps
  init <- winit[is_unpenalized]
  upper <- c(11.5, 18.2, 5)
  lower <- c(2.3, -2.5, 1.4)
  ans <- optim(init, nll_no_penalty, lower = lower, upper=upper, method = "L-BFGS-B")
  if (make_log) {
    logfile <- tempfile(pattern="gpnet", fileext = ".log")
    record <- function(...) cat(..., file = logfile, append = TRUE)
  }
  par <- winit
  par[is_unpenalized] <- ans$par
  nderv <- function(f, p) numDeriv::grad(f, x = p, method = "simple", method.args=list(eps=1e-3))
  gnll <- nderv(nll, par)
  mu <- mubar
  stopifnot(beta>0, beta<1)
  G <- diag(initFactor * abs(gnll), ncol=dim)
  if(flmin<1){
    lstart <- max(abs(gnll))
    loglstart <- log10(lstart) + relStart
    log10LambdaRange <- -log10(flmin)
    loglend <- loglstart - log10LambdaRange
    logstep <- (loglend - loglstart)/nlam
    lambda <- loglstart + 0:(nlam - 1)*logstep
    lambda <- 10^lambda
  } else {
    lambda <- ulam
  }
  res <- list()
  fsg <- function(p, g, l1, l2, h) {
    if(abs(p) < .Machine$double.eps) {
      max(abs(g) - l1, 0)
    } else if (p > 0){
      (-g - l1 - l2*p)/(h + l2)
    } else {
      (-g + l1 - l2*p)/(h + l2)
    }
  }
  for (i in seq_along(lambda)){
    l1penalty <- lambda[i] * vp * alpha
    l2penalty <- lambda[i] * vp * (1 - alpha)
    penFunc <- function(par) abs(par) %*% l1penalty + (par^2 %*% l2penalty)/2
    k <- 0
    nlp <- nll(par)
    F1 <- nlp + penFunc(par)
    H <- I/(2*mu) + G
    sg <- mapply(fsg, p=par, g=gnll, l1=l1penalty, l2=l2penalty, h=diag(H))
    nsg <- max(abs(sg))
    while (nsg > thresh && k < maxit){
      H <- I/(2*mu) + G
      d <- numeric(dim)
      dlist <- list()
      fmlist <- list()
      inactive <- par != 0 | sg != 0 ## inactive means not actively fixed to zero in solution
      nInactive <- sum(inactive)
      ndesc <- (i + floor(k * a)) * nInactive
      for (nd in 1:ndesc){
        j <- sample.int(n=nInactive, size=1)
        j <- parInds[inactive][j]
        Hd <- H %*% d
        gr <- gnll[j] + Hd[j]
        if (par[j] + d[j] > 0 || (abs(par[j] + d[j]) < .Machine$double.eps & -gr > 0)){
          z <- (-gr - l1penalty[j] - l2penalty[j]*(par[j] + d[j]))/(H[j,j] + l2penalty[j])
          if (par[j] + d[j] + z < 0){
            d[j] <- -par[j]
          } else {
            d[j] <- d[j] + z
          }
        } else {
          z <- (-gr + l1penalty[j] - l2penalty[j]*(par[j] + d[j]))/(H[j,j] + l2penalty[j])
          if (par[j] + d[j] + z > 0){
            d[j] <- -par[j]
          } else {
            d[j] <- d[j] + z
          }
        }
        dlist[[nd]] <- d
        fmlist[[nd]] <- d %*% (H/2) %*% d + gnll %*% d + F1 - penFunc(par) + penFunc(par + d)
      }
      if(any(diff(c(F1, unlist(fmlist)))>1e-10)) browser()
      par2 <- par + d
      if (any(par2[!is_unpenalized] > cl[2, ] || any(par2[!is_unpenalized] < cl[1, ]))){
        if(make_log) record('backtracking: out of bounds', '\n')
        mu <- mu * beta
      } else {
        Fmod <- d %*% (H/2) %*% d + gnll %*% d + F1 - penFunc(par) + penFunc(par2)
        if(Fmod  > F1 && !isTRUE(all.equal(Fmod, F1))) {
          if (debug) browser()
          if(make_log) record('backtracking: poorly solved model', '\n')
          mu <- mu * beta
        } else {
          nlp <- nll(par2)
          F2 <- nlp + penFunc(par2)
          steplen <- sqrt(sum(d^2))
          if (F2 - F1 > r * (Fmod - F1) && steplen > 1e-3){
            if(make_log) record('backtracking: insufficient decrease', '\n')
            mu <- mu * beta
            if (mu < 1e-8) {
              if(make_log) record("problem: mu < 1e-8", "\n")
              browser()
              gnll <- nderv(nll, par2)
              G <- diag(initFactor * abs(gnll), ncol=dim)
              mu <- mubar
            }
          } else {
            if (mu < mubar) {
              if(make_log) record('increasing step size: sufficient decrease', '\n')
              mu <- mu / sqrt(beta)
            }
            gnll2 <- nderv(nll, par2)
            yvec <- gnll2 - gnll
            s <- d
            ys <- yvec %*% s
            if (ys > 0){
              ## Hessian approximation update via BFGS via 8.19 in Nocedal and Wright
              sColVec <- matrix(s, ncol=1)
              yColVec <- matrix(yvec, ncol=1)
              M <- G %*% sColVec %*% t(sColVec) %*% G
              M <- M / as.numeric(t(sColVec) %*% G %*% sColVec)
              M <- G - M
              G <- M + yColVec %*% t(yColVec) / as.numeric(t(yColVec) %*% sColVec)
              if (any(diag(G) < 0)) browser()
            } else {
              if(make_log) record('resetting Hessian: ys <= 0', '\n')
              if (isTRUE(all.equal(ys, 0))) {
                if(make_log) record("problem: ys == 0", "\n")
                browser()
              }
              G <- diag(initFactor * abs(gnll2), ncol=dim)
              mu <- mubar
            }
            ## update vars
            k <- k + 1
            par <- par2
            gnll <- gnll2
            dF <- F2 - F1
            F1 <- F2
            H <- I/(2*mu) + G
            sg <- mapply(fsg, p=par, g=gnll, l1=l1penalty, l2=l2penalty, h=diag(H))
            nsg <- max(abs(sg))
            if (make_log) {
              record('lambda: ', lambda[i], '\n')
              record('k: ', k, '\n')
              record('F: ', as.numeric(F1), '\n')
              record('par: ', signif(par, 3), '\n')
              record('grad: ', signif(gnll, 3), '\n')
              record('h: ', signif(diag(H), 3), '\n')
              record('nsg: ', signif(nsg, 3), '\n')
              record('mu: ', signif(mu, 3), '\n')
              record('\n')
            }
          }
        }
      }
    }
    convergence <- ifelse(k == maxit, 'no', 'yes')
    res[[i]] <- list(par=par, nll=nlp, k=k, lambda=lambda[i],
                     convergence=convergence, mu=mu, nsg=nsg, sg=sg)
    if (make_log) record("resetting Hessian: new penalty", "\n")
    mu <- mubar
    G <- diag(initFactor * abs(gnll), ncol=dim)
  }
  path <- sapply(res, '[[', 'par')
  ret <- list(a0=path[is_unpenalized,])
  beta <- path[!is_unpenalized, , drop=FALSE]
  colnames(beta) <- paste("s", seq(ncol(beta)) - 1, sep = "")
  rownames(beta) <- vnames
  ret$beta <- beta
  ret$lambda <- sapply(res, '[[', 'lambda')
  ret$nll <- sapply(res, '[[', 'nll')
  ret$df <- colSums(abs(beta) > 0)
  ret$dim <- dim(x)
  ret$niterations <- sum(sapply(res, '[[', 'k'))
  ret$jerr <- paste('convergence: ', sapply(res, '[[', 'convergence'))
  if(make_log) record('Completed regularization path', '\n')
  ret
}

plotCoef <- function (beta, norm, lambda, df, dev, label = FALSE, xvar = c("norm",
                                                                           "lambda", "dev"), xlab = iname, ylab = "Coefficients", ...) {
  which <- which(rowSums(abs(beta)) > 0)
  nwhich <- length(which)
  switch(nwhich + 1, `0` = {
    warning("No plot produced since all coefficients zero")
    return()
  }, `1` = warning("1 or fewer nonzero coefficients; glmnet plot is not meaningful"))
  beta = as.matrix(beta[which, , drop = FALSE])
  xvar = match.arg(xvar)
  switch(xvar, norm = {
    index = if (missing(norm)) apply(abs(beta), 2, sum) else norm
    iname = "L1 Norm"
    approx.f = 1
  }, lambda = {
    index = log(lambda)
    iname = "Log Lambda"
    approx.f = 0
  }, dev = {
    index = dev
    iname = "Fraction Deviance Explained"
    approx.f = 1
  })
  dotlist = list(...)
  type = dotlist$type
  if (is.null(type))
    graphics::matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab,
                      type = "l", ...)
  else graphics::matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab,
                         ...)
  atdf = pretty(index)
  prettydf = stats::approx(x = index, y = df, xout = atdf, rule = 2,
                           method = "constant", f = approx.f)$y
  graphics::axis(3, at = atdf, labels = prettydf, tcl = NA)
  if (label) {
    nnz = length(which)
    xpos = max(index)
    pos = 4
    if (xvar == "lambda") {
      xpos = min(index)
      pos = 2
    }
    xpos = rep(xpos, nnz)
    ypos = beta[, ncol(beta)]
    graphics::text(xpos, ypos, paste(which), cex = 0.5, pos = pos)
  }
}

#' @export
plot.gpnet <- function (x, xvar = c("norm", "lambda", "dev"),
                        label = FALSE, ...){
  xvar = match.arg(xvar)
  plotCoef(x$beta, lambda = x$lambda, df = x$df, dev = x$dev.ratio,
           label = label, xvar = xvar, ...)
}

#' @export
print.gpnet <- function(x, digits = max(3, getOption("digits") - 3), ...){
  cat("\nCall: ", deparse(x$call), "\n\n")
  print(cbind(Df = x$df, `ll` = signif(-x$nll, digits),
              Lambda = signif(x$lambda, digits)))
}
