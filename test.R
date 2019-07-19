function (y, d = NA, D = NA, max.p = 5, max.q = 5, max.P = 2, 
          max.Q = 2, max.order = 5, max.d = 2, max.D = 1, start.p = 2, 
          start.q = 2, start.P = 1, start.Q = 1, stationary = FALSE, 
          seasonal = TRUE, ic = c("aicc", "aic", "bic"), stepwise = TRUE, 
          nmodels = 94, trace = FALSE, approximation = (length(x) > 
                                                          150 | frequency(x) > 12), method = NULL, truncate = NULL, 
          xreg = NULL, test = c("kpss", "adf", "pp"), test.args = list(), 
          seasonal.test = c("seas", "ocsb", "hegy", "ch"), seasonal.test.args = list(), 
          allowdrift = TRUE, allowmean = TRUE, lambda = NULL, biasadj = FALSE, 
          parallel = FALSE, num.cores = 2, x = y, ...) 
{
  if (stepwise && parallel) {
    warning("Parallel computer is only implemented when stepwise=FALSE, the model will be fit in serial.")
    parallel <- FALSE
  }
  if (trace && parallel) {
    message("Tracing model searching in parallel is not supported.")
    trace <- FALSE
  }
  series <- deparse(substitute(y))
  x <- as.ts(x)
  if (NCOL(x) > 1) {
    stop("auto.arima can only handle univariate time series")
  }
  orig.x <- x
  missing <- is.na(x)
  firstnonmiss <- head(which(!missing), 1)
  lastnonmiss <- tail(which(!missing), 1)
  serieslength <- sum(!missing[firstnonmiss:lastnonmiss])
  x <- subset(x, start = firstnonmiss)
  if (is.constant(x)) {
    if (all(is.na(x))) 
      stop("All data are missing")
    if (allowmean) {
      fit <- Arima(x, order = c(0, 0, 0), fixed = mean(x, 
                                                       na.rm = TRUE), ...)
    }
    else {
      fit <- Arima(x, order = c(0, 0, 0), include.mean = FALSE, 
                   ...)
    }
    fit$x <- orig.x
    fit$series <- series
    fit$call <- match.call()
    fit$call$x <- data.frame(x = x)
    fit$constant <- TRUE
    return(fit)
  }
  ic <- match.arg(ic)
  test <- match.arg(test)
  seasonal.test <- match.arg(seasonal.test)
  if (seasonal) {
    m <- frequency(x)
  }
  else {
    m <- 1
  }
  if (m < 1) {
    m <- 1
  }
  else {
    m <- round(m)
  }
  max.p <- min(max.p, floor(serieslength/3))
  max.q <- min(max.q, floor(serieslength/3))
  max.P <- min(max.P, floor(serieslength/3/m))
  max.Q <- min(max.Q, floor(serieslength/3/m))
  if (serieslength <= 3L) {
    ic <- "aic"
  }
  if (!is.null(lambda)) {
    x <- BoxCox(x, lambda)
    lambda <- attr(x, "lambda")
    attr(lambda, "biasadj") <- biasadj
  }
  if (!is.null(xreg)) {
    if (!is.numeric(xreg)) 
      stop("xreg should be a numeric matrix or a numeric vector")
    xreg <- as.matrix(xreg)
    if (is.null(colnames(xreg))) {
      colnames(xreg) <- if (ncol(xreg) == 1) 
        "xreg"
      else paste("xreg", 1:ncol(xreg), sep = "")
    }
    xregg <- xreg
    xx <- x
    constant_columns <- apply(xregg, 2, is.constant)
    if (all(constant_columns)) {
      xregg <- NULL
    }
    else {
      if (any(constant_columns)) {
        xregg <- xregg[, -which(constant_columns)]
      }
      sv <- svd(na.omit(cbind(rep(1, NROW(xregg)), xregg)))$d
      if (min(sv)/sum(sv) < .Machine$double.eps) {
        stop("xreg is rank deficient")
      }
      j <- !is.na(x) & !is.na(rowSums(xregg))
      xx[j] <- residuals(lm(x ~ xregg))
    }
  }
  else {
    xx <- x
    xregg <- NULL
  }
  if (stationary) {
    d <- D <- 0
  }
  if (m == 1) {
    D <- max.P <- max.Q <- 0
  }
  else if (is.na(D)) {
    D <- do.call("nsdiffs", c(list(xx, test = seasonal.test, 
                                   max.D = max.D), seasonal.test.args))
    if (D > 0 && !is.null(xregg)) {
      diffxreg <- diff(xregg, differences = D, lag = m)
      if (any(apply(diffxreg, 2, is.constant))) {
        D <- D - 1
      }
    }
    if (D > 0) {
      dx <- diff(xx, differences = D, lag = m)
      if (all(is.na(dx))) 
        D <- D - 1
    }
  }
  if (D > 0) {
    dx <- diff(xx, differences = D, lag = m)
  }
  else {
    dx <- xx
  }
  if (!is.null(xregg)) {
    if (D > 0) {
      diffxreg <- diff(xregg, differences = D, lag = m)
    }
    else {
      diffxreg <- xregg
    }
  }
  if (is.na(d)) {
    d <- do.call("ndiffs", c(list(dx, test = test, max.d = max.d), 
                             test.args))
    if (d > 0 && !is.null(xregg)) {
      diffxreg <- diff(diffxreg, differences = d, lag = 1)
      if (any(apply(diffxreg, 2, is.constant))) {
        d <- d - 1
      }
    }
    if (d > 0) {
      diffdx <- diff(dx, differences = d, lag = 1)
      if (all(is.na(diffdx))) 
        d <- d - 1
    }
  }
  if (D >= 2) {
    warning("Having more than one seasonal differences is not recommended. Please consider using only one seasonal difference.")
  }
  else if (D + d > 2) {
    warning("Having 3 or more differencing operations is not recommended. Please consider reducing the total number of differences.")
  }
  if (d > 0) {
    dx <- diff(dx, differences = d, lag = 1)
  }
  if (is.constant(dx)) {
    if (is.null(xreg)) {
      if (D > 0 && d == 0) {
        fit <- Arima(x, order = c(0, d, 0), seasonal = list(order = c(0, 
                                                                      D, 0), period = m), include.constant = TRUE, 
                     fixed = mean(dx/m, na.rm = TRUE), method = method, 
                     ...)
      }
      else if (D > 0 && d > 0) {
        fit <- Arima(x, order = c(0, d, 0), seasonal = list(order = c(0, 
                                                                      D, 0), period = m), method = method, ...)
      }
      else if (d == 2) {
        fit <- Arima(x, order = c(0, d, 0), method = method, 
                     ...)
      }
      else if (d < 2) {
        fit <- Arima(x, order = c(0, d, 0), include.constant = TRUE, 
                     fixed = mean(dx, na.rm = TRUE), method = method, 
                     ...)
      }
      else {
        stop("Data follow a simple polynomial and are not suitable for ARIMA modelling.")
      }
    }
    else {
      if (D > 0) {
        fit <- Arima(x, order = c(0, d, 0), seasonal = list(order = c(0, 
                                                                      D, 0), period = m), xreg = xreg, method = method, 
                     ...)
      }
      else {
        fit <- Arima(x, order = c(0, d, 0), xreg = xreg, 
                     method = method, ...)
      }
    }
    fit$x <- orig.x
    fit$series <- series
    fit$call <- match.call()
    fit$call$x <- data.frame(x = x)
    return(fit)
  }
  if (m > 1) {
    if (max.P > 0) {
      max.p <- min(max.p, m - 1)
    }
    if (max.Q > 0) {
      max.q <- min(max.q, m - 1)
    }
  }
  if (approximation) {
    if (!is.null(truncate)) {
      tspx <- tsp(x)
      if (length(x) > truncate) {
        x <- ts(tail(x, truncate), end = tspx[2], frequency = tspx[3])
      }
    }
    if (D == 0) {
      fit <- try(stats::arima(x, order = c(0, d, 0), xreg = xreg, 
                              ...), silent = TRUE)
    }
    else {
      fit <- try(stats::arima(x, order = c(0, d, 0), seasonal = list(order = c(0, 
                                                                               D, 0), period = m), xreg = xreg, ...), silent = TRUE)
    }
    if (!is.element("try-error", class(fit))) {
      offset <- -2 * fit$loglik - serieslength * log(fit$sigma2)
    }
    else {
      offset <- 0
    }
  }
  else {
    offset <- 0
  }
  allowdrift <- allowdrift & (d + D) == 1
  allowmean <- allowmean & (d + D) == 0
  constant <- allowdrift | allowmean
  if (approximation && trace) {
    cat("\n Fitting models using approximations to speed things up...\n")
  }
  if (!stepwise) {
    bestfit <- search.arima(x, d, D, max.p, max.q, max.P, 
                            max.Q, max.order, stationary, ic, trace, approximation, 
                            method = method, xreg = xreg, offset = offset, allowdrift = allowdrift, 
                            allowmean = allowmean, parallel = parallel, num.cores = num.cores, 
                            ...)
    bestfit$call <- match.call()
    bestfit$call$x <- data.frame(x = x)
    bestfit$lambda <- lambda
    bestfit$x <- orig.x
    bestfit$series <- series
    bestfit$fitted <- fitted.Arima(bestfit)
    if (trace) {
      cat("\n\n Best model:", arima.string(bestfit, padding = TRUE), 
          "\n\n")
    }
    return(bestfit)
  }
  if (length(x) < 10L) {
    start.p <- min(start.p, 1L)
    start.q <- min(start.q, 1L)
    start.P <- 0L
    start.Q <- 0L
  }
  p <- start.p <- min(start.p, max.p)
  q <- start.q <- min(start.q, max.q)
  P <- start.P <- min(start.P, max.P)
  Q <- start.Q <- min(start.Q, max.Q)
  results <- matrix(NA, nrow = nmodels, ncol = 8)
  bestfit <- myarima(x, order = c(p, d, q), seasonal = c(P, 
                                                         D, Q), constant = constant, ic, trace, approximation, 
                     method = method, offset = offset, xreg = xreg, ...)
  results[1, ] <- c(p, d, q, P, D, Q, constant, bestfit$ic)
  fit <- myarima(x, order = c(0, d, 0), seasonal = c(0, D, 
                                                     0), constant = constant, ic, trace, approximation, method = method, 
                 offset = offset, xreg = xreg, ...)
  results[2, ] <- c(0, d, 0, 0, D, 0, constant, fit$ic)
  if (fit$ic < bestfit$ic) {
    bestfit <- fit
    p <- q <- P <- Q <- 0
  }
  k <- 2
  if (max.p > 0 || max.P > 0) {
    fit <- myarima(x, order = c(max.p > 0, d, 0), seasonal = c((m > 
                                                                  1) & (max.P > 0), D, 0), constant = constant, ic, 
                   trace, approximation, method = method, offset = offset, 
                   xreg = xreg, ...)
    results[k + 1, ] <- c(max.p > 0, d, 0, (m > 1) & (max.P > 
                                                        0), D, 0, constant, fit$ic)
    if (fit$ic < bestfit$ic) {
      bestfit <- fit
      p <- (max.p > 0)
      P <- (m > 1) & (max.P > 0)
      q <- Q <- 0
    }
    k <- k + 1
  }
  if (max.q > 0 || max.Q > 0) {
    fit <- myarima(x, order = c(0, d, max.q > 0), seasonal = c(0, 
                                                               D, (m > 1) & (max.Q > 0)), constant = constant, 
                   ic, trace, approximation, method = method, offset = offset, 
                   xreg = xreg, ...)
    results[k + 1, ] <- c(0, d, max.q > 0, 0, D, (m > 1) & 
                            (max.Q > 0), constant, fit$ic)
    if (fit$ic < bestfit$ic) {
      bestfit <- fit
      p <- P <- 0
      Q <- (m > 1) & (max.Q > 0)
      q <- (max.q > 0)
    }
    k <- k + 1
  }
  if (constant) {
    fit <- myarima(x, order = c(0, d, 0), seasonal = c(0, 
                                                       D, 0), constant = FALSE, ic, trace, approximation, 
                   method = method, offset = offset, xreg = xreg, ...)
    results[k + 1, ] <- c(0, d, 0, 0, D, 0, 0, fit$ic)
    if (fit$ic < bestfit$ic) {
      bestfit <- fit
      p <- q <- P <- Q <- 0
    }
    k <- k + 1
  }
  startk <- 0
  while (startk < k && k < nmodels) {
    startk <- k
    if (P > 0 && newmodel(p, d, q, P - 1, D, Q, constant, 
                          results[1:k, ])) {
      k <- k + 1
      if (k > nmodels) 
        next
      fit <- myarima(x, order = c(p, d, q), seasonal = c(P - 
                                                           1, D, Q), constant = constant, ic, trace, approximation, 
                     method = method, offset = offset, xreg = xreg, 
                     ...)
      results[k, ] <- c(p, d, q, P - 1, D, Q, constant, 
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        P <- (P - 1)
        next
      }
    }
    if (Q > 0 && newmodel(p, d, q, P, D, Q - 1, constant, 
                          results[1:k, ])) {
      k <- k + 1
      if (k > nmodels) 
        next
      fit <- myarima(x, order = c(p, d, q), seasonal = c(P, 
                                                         D, Q - 1), constant = constant, ic, trace, approximation, 
                     method = method, offset = offset, xreg = xreg, 
                     ...)
      results[k, ] <- c(p, d, q, P, D, Q - 1, constant, 
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        Q <- (Q - 1)
        next
      }
    }
    if (P < max.P && newmodel(p, d, q, P + 1, D, Q, constant, 
                              results[1:k, ])) {
      k <- k + 1
      if (k > nmodels) 
        next
      fit <- myarima(x, order = c(p, d, q), seasonal = c(P + 
                                                           1, D, Q), constant = constant, ic, trace, approximation, 
                     method = method, offset = offset, xreg = xreg, 
                     ...)
      results[k, ] <- c(p, d, q, P + 1, D, Q, constant, 
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        P <- (P + 1)
        next
      }
    }
    if (Q < max.Q && newmodel(p, d, q, P, D, Q + 1, constant, 
                              results[1:k, ])) {
      k <- k + 1
      if (k > nmodels) 
        next
      fit <- myarima(x, order = c(p, d, q), seasonal = c(P, 
                                                         D, Q + 1), constant = constant, ic, trace, approximation, 
                     method = method, offset = offset, xreg = xreg, 
                     ...)
      results[k, ] <- c(p, d, q, P, D, Q + 1, constant, 
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        Q <- (Q + 1)
        next
      }
    }
    if (Q > 0 && P > 0 && newmodel(p, d, q, P - 1, D, Q - 
                                   1, constant, results[1:k, ])) {
      k <- k + 1
      if (k > nmodels) 
        next
      fit <- myarima(x, order = c(p, d, q), seasonal = c(P - 
                                                           1, D, Q - 1), constant = constant, ic, trace, 
                     approximation, method = method, offset = offset, 
                     xreg = xreg, ...)
      results[k, ] <- c(p, d, q, P - 1, D, Q - 1, constant, 
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        Q <- (Q - 1)
        P <- (P - 1)
        next
      }
    }
    if (Q < max.Q && P > 0 && newmodel(p, d, q, P - 1, D, 
                                       Q + 1, constant, results[1:k, ])) {
      k <- k + 1
      if (k > nmodels) 
        next
      fit <- myarima(x, order = c(p, d, q), seasonal = c(P - 
                                                           1, D, Q + 1), constant = constant, ic, trace, 
                     approximation, method = method, offset = offset, 
                     xreg = xreg, ...)
      results[k, ] <- c(p, d, q, P - 1, D, Q + 1, constant, 
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        Q <- (Q + 1)
        P <- (P - 1)
        next
      }
    }
    if (Q > 0 && P < max.P && newmodel(p, d, q, P + 1, D, 
                                       Q - 1, constant, results[1:k, ])) {
      k <- k + 1
      if (k > nmodels) 
        next
      fit <- myarima(x, order = c(p, d, q), seasonal = c(P + 
                                                           1, D, Q - 1), constant = constant, ic, trace, 
                     approximation, method = method, offset = offset, 
                     xreg = xreg, ...)
      results[k, ] <- c(p, d, q, P + 1, D, Q - 1, constant, 
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        Q <- (Q - 1)
        P <- (P + 1)
        next
      }
    }
    if (Q < max.Q && P < max.P && newmodel(p, d, q, P + 
                                           1, D, Q + 1, constant, results[1:k, ])) {
      k <- k + 1
      if (k > nmodels) 
        next
      fit <- myarima(x, order = c(p, d, q), seasonal = c(P + 
                                                           1, D, Q + 1), constant = constant, ic, trace, 
                     approximation, method = method, offset = offset, 
                     xreg = xreg, ...)
      results[k, ] <- c(p, d, q, P + 1, D, Q + 1, constant, 
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        Q <- (Q + 1)
        P <- (P + 1)
        next
      }
    }
    if (p > 0 && newmodel(p - 1, d, q, P, D, Q, constant, 
                          results[1:k, ])) {
      k <- k + 1
      if (k > nmodels) 
        next
      fit <- myarima(x, order = c(p - 1, d, q), seasonal = c(P, 
                                                             D, Q), constant = constant, ic, trace, approximation, 
                     method = method, offset = offset, xreg = xreg, 
                     ...)
      results[k, ] <- c(p - 1, d, q, P, D, Q, constant, 
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        p <- (p - 1)
        next
      }
    }
    if (q > 0 && newmodel(p, d, q - 1, P, D, Q, constant, 
                          results[1:k, ])) {
      k <- k + 1
      if (k > nmodels) 
        next
      fit <- myarima(x, order = c(p, d, q - 1), seasonal = c(P, 
                                                             D, Q), constant = constant, ic, trace, approximation, 
                     method = method, offset = offset, xreg = xreg, 
                     ...)
      results[k, ] <- c(p, d, q - 1, P, D, Q, constant, 
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- (q - 1)
        next
      }
    }
    if (p < max.p && newmodel(p + 1, d, q, P, D, Q, constant, 
                              results[1:k, ])) {
      k <- k + 1
      if (k > nmodels) 
        next
      fit <- myarima(x, order = c(p + 1, d, q), seasonal = c(P, 
                                                             D, Q), constant = constant, ic, trace, approximation, 
                     method = method, offset = offset, xreg = xreg, 
                     ...)
      results[k, ] <- c(p + 1, d, q, P, D, Q, constant, 
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        p <- (p + 1)
        next
      }
    }
    if (q < max.q && newmodel(p, d, q + 1, P, D, Q, constant, 
                              results[1:k, ])) {
      k <- k + 1
      if (k > nmodels) 
        next
      fit <- myarima(x, order = c(p, d, q + 1), seasonal = c(P, 
                                                             D, Q), constant = constant, ic, trace, approximation, 
                     method = method, offset = offset, xreg = xreg, 
                     ...)
      results[k, ] <- c(p, d, q + 1, P, D, Q, constant, 
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- (q + 1)
        next
      }
    }
    if (q > 0 && p > 0 && newmodel(p - 1, d, q - 1, P, D, 
                                   Q, constant, results[1:k, ])) {
      k <- k + 1
      if (k > nmodels) 
        next
      fit <- myarima(x, order = c(p - 1, d, q - 1), seasonal = c(P, 
                                                                 D, Q), constant = constant, ic, trace, approximation, 
                     method = method, offset = offset, xreg = xreg, 
                     ...)
      results[k, ] <- c(p - 1, d, q - 1, P, D, Q, constant, 
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- (q - 1)
        p <- (p - 1)
        next
      }
    }
    if (q < max.q && p > 0 && newmodel(p - 1, d, q + 1, 
                                       P, D, Q, constant, results[1:k, ])) {
      k <- k + 1
      if (k > nmodels) 
        next
      fit <- myarima(x, order = c(p - 1, d, q + 1), seasonal = c(P, 
                                                                 D, Q), constant = constant, ic, trace, approximation, 
                     method = method, offset = offset, xreg = xreg, 
                     ...)
      results[k, ] <- c(p - 1, d, q + 1, P, D, Q, constant, 
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- (q + 1)
        p <- (p - 1)
        next
      }
    }
    if (q > 0 && p < max.p && newmodel(p + 1, d, q - 1, 
                                       P, D, Q, constant, results[1:k, ])) {
      k <- k + 1
      if (k > nmodels) 
        next
      fit <- myarima(x, order = c(p + 1, d, q - 1), seasonal = c(P, 
                                                                 D, Q), constant = constant, ic, trace, approximation, 
                     method = method, offset = offset, xreg = xreg, 
                     ...)
      results[k, ] <- c(p + 1, d, q - 1, P, D, Q, constant, 
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- (q - 1)
        p <- (p + 1)
        next
      }
    }
    if (q < max.q && p < max.p && newmodel(p + 1, d, q + 
                                           1, P, D, Q, constant, results[1:k, ])) {
      k <- k + 1
      if (k > nmodels) 
        next
      fit <- myarima(x, order = c(p + 1, d, q + 1), seasonal = c(P, 
                                                                 D, Q), constant = constant, ic, trace, approximation, 
                     method = method, offset = offset, xreg = xreg, 
                     ...)
      results[k, ] <- c(p + 1, d, q + 1, P, D, Q, constant, 
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- (q + 1)
        p <- (p + 1)
        next
      }
    }
    if (allowdrift || allowmean) {
      if (newmodel(p, d, q, P, D, Q, !constant, results[1:k, 
                                                        ])) {
        k <- k + 1
        if (k > nmodels) 
          next
        fit <- myarima(x, order = c(p, d, q), seasonal = c(P, 
                                                           D, Q), constant = !constant, ic, trace, approximation, 
                       method = method, offset = offset, xreg = xreg, 
                       ...)
        results[k, ] <- c(p, d, q, P, D, Q, !constant, 
                          fit$ic)
        if (fit$ic < bestfit$ic) {
          bestfit <- fit
          constant <- !constant
        }
      }
    }
  }
  if (k > nmodels) {
    warning(sprintf("Stepwise search was stopped early due to reaching the model number limit: `nmodels = %i`", 
                    nmodels))
  }
  if (approximation && !is.null(bestfit$arma)) {
    if (trace) {
      cat("\n\n Now re-fitting the best model(s) without approximations...\n")
    }
    icorder <- order(results[, 8])
    nmodels <- sum(!is.na(results[, 8]))
    for (i in seq(nmodels)) {
      k <- icorder[i]
      fit <- myarima(x, order = c(results[k, 1], d, results[k, 
                                                            3]), seasonal = c(results[k, 4], D, results[k, 
                                                                                                        6]), constant = results[k, 7] == 1, ic, trace, 
                     approximation = FALSE, method = method, xreg = xreg, 
                     ...)
      if (fit$ic < Inf) {
        bestfit <- fit
        break
      }
    }
  }
  if (bestfit$ic == Inf && !isTRUE(method == "CSS")) {
    if (trace) {
      cat("\n")
    }
    stop("No suitable ARIMA model found")
  }
  bestfit$x <- orig.x
  bestfit$series <- series
  bestfit$ic <- NULL
  bestfit$call <- match.call()
  bestfit$call$x <- data.frame(x = x)
  bestfit$lambda <- lambda
  bestfit$fitted <- fitted.Arima(bestfit)
  if (trace) {
    cat("\n\n Best model:", arima.string(bestfit, padding = TRUE), 
        "\n\n")
  }
  return(bestfit)
}
