#' Condition Index(CI) in R
#' @param m is a matrix or data.frame
#' @author Syed Hammad Raza
#' @references Inspired by Gold Medalist Dr Shahla Faisal who teach me in GC University Faisalabad
#' @description If the Condition Index is between 10 and 30, there is moderate to strong multicollinearity and if it exceeds 30 there is severe multicollinearity
#' @export

CI <- function(m)
{
  R <- cor(m)
  cat("The Correlation matrix is: \n")
  print(R)
  cat("\n")
  l <- eigen(R)$values
  l_max <- max(l)
  cat("The maximun eigenvalue is: \n")
  print(l_max)
  cat("\n")
  l_min <- min(l)
  cat("The minimun eigenvalue is: \n")
  print(l_min)
  cat("\n")
  ci <- sqrt(l_max / l_min)
  cat("The Value of Condition Index(CI) is: \n")
  print(ci)
}

#' Variance Inflation Method in R
#' @param mod is a Regression Model which we can fit by lm() command
#' @author Syed Hammad Raza
#' @references Inspired by Gold Medalist Dr Shahla Faisal who teach me in GC University Faisalabad
#' @description When x  is nearly linearly independent on subset of remaining regressor. If R^2 is nearly 1 it mean VIF is large. As a rule of thumb high value if VIF of variable is exceed 10, mean multicollinearity is exit.
#' @export

vif <- function (mod, ...)
{
  if (any(is.na(coef(mod))))
    stop("there are aliased coefficients in the model")
  v <- vcov(mod)
  assign <- attr(model.matrix(mod), "assign")
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  }
  else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2)
    stop("model contains fewer than 2 terms")
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/(2*Df))")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs,
                                                                       -subs]))/detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1))
    result <- result[, 1]
  else result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  result
}

#' Goldfeld-Quandt Tes
#' @description Goldfeld-Quandt test against heteroscedasticity
#' @param formula a symbolic description for the model to be tested (or a fitted "lm" object).
#' @param point numerical. If point is smaller than 1 it is interpreted as percentages of data, i.e. n*point is taken to be the (potential) breakpoint in the variances, if n is the number of observations in the model. If point is greater than 1 it is interpreted to be the index of the breakpoint.
#' @param fraction numerical. The number of central observations to be omitted. If fraction is smaller than 1, it is chosen to be fraction*n if n is the number of observations in the model.
#' @param alternative a character string specifying the alternative hypothesis. The default is to test for increasing variances.
#' @param order.by Either a vector z or a formula with a single explanatory variable like ~ z. The observations in the model are ordered by the size of z. If set to NULL (the default) the observations are assumed to be ordered (e.g., a time series).
#' @param data an optional data frame containing the variables in the model. By default the variables are taken from the environment which gqtest is called from.
#' @details The Goldfeld-Quandt test compares the variances of two submodels divided by a specified breakpoint and rejects if the variances differ.
#' @references Inspired by Gold Medalist Dr Shahla Faisal who teach me in GC University Faisalabad
#' @author Syed Hammad Raza
#' @export

gqtest <- function (formula, point = 0.5, fraction = 0, alternative = c("greater",
                                                                        "two.sided", "less"), order.by = NULL, data = list())
{
  dname <- paste(deparse(substitute(formula)))
  alternative <- match.arg(alternative)
  if (!inherits(formula, "formula")) {
    X <- if (is.matrix(formula$x))
      formula$x
    else model.matrix(terms(formula), model.frame(formula))
    y <- if (is.vector(formula$y))
      formula$y
    else model.response(model.frame(formula))
  }
  else {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
  }
  k <- ncol(X)
  n <- nrow(X)
  if (point > 1) {
    if (fraction < 1)
      fraction <- floor(fraction * n)
    point1 <- point - ceiling(fraction/2)
    point2 <- point + ceiling(fraction/2 + 0.01)
  }
  else {
    if (fraction >= 1)
      fraction <- fraction/n
    point1 <- floor((point - fraction/2) * n)
    point2 <- ceiling((point + fraction/2) * n + 0.01)
  }
  if (point2 > n - k + 1 | point1 < k)
    stop("inadmissable breakpoint/too many central observations omitted")
  if (!is.null(order.by)) {
    if (inherits(order.by, "formula")) {
      z <- model.matrix(order.by, data = data)
      z <- as.vector(z[, ncol(z)])
    }
    else {
      z <- order.by
    }
    X <- as.matrix(X[order(z), ])
    y <- y[order(z)]
  }
  rss1 <- sum(lm.fit(as.matrix(X[1:point1, ]), y[1:point1])$residuals^2)
  rss2 <- sum(lm.fit(as.matrix(X[point2:n, ]), y[point2:n])$residuals^2)
  mss <- c(rss1/(point1 - k), rss2/(n - point2 + 1 - k))
  gq <- mss[2]/mss[1]
  df <- c(n - point2 + 1 - k, point1 - k)
  names(df) <- c("df1", "df2")
  PVAL <- switch(alternative, two.sided = (2 * min(pf(gq,
                                                      df[1], df[2]), pf(gq, df[1], df[2], lower.tail = FALSE))),
                 less = pf(gq, df[1], df[2]), greater = pf(gq, df[1],
                                                           df[2], lower.tail = FALSE))
  alternative <- switch(alternative, two.sided = "variance changes from segment 1 to 2",
                        less = "variance decreases from segment 1 to 2", greater = "variance increases from segment 1 to 2")
  method <- "Goldfeld-Quandt test"
  names(gq) <- "GQ"
  RVAL <- list(statistic = gq, parameter = df, method = method,
               alternative = alternative, p.value = PVAL, data.name = dname)
  class(RVAL) <- "htest"
  return(RVAL)
}

#' Breusch-Pagan Tes
#' @description  Performs the Breusch-Pagan test against heteroskedasticity.
#' @param formula symbolic description for the model to be tested (or a fitted "lm" object)
#' @param varformula a formula describing only the potential explanatory variables for the variance (no dependent variable needed). By default the same explanatory variables are taken as in the main regression model.
#' @param data an optional data frame containing the variables in the model. By default the variables are taken from the environment which bptest is called from.
#' @param weights an optional vector of weights to be used in the model.
#' @details The Breusch-Pagan test fits a linear regression model to the residuals of a linear regression model (by default the same explanatory variables are taken as in the main regression model) and rejects if too much of the variance is explained by the additional explanatory variables
#' @references Inspired by Gold Medalist Dr Shahla Faisal who teach me in GC University Faisalabad
#' @author Syed Hammad Raza
#' @export

bptest <- function (formula, varformula = NULL, data = list(),
                    weights = NULL)
{
  dname <- paste(deparse(substitute(formula)))
  if (!inherits(formula, "formula")) {
    X <- if (is.matrix(formula$x))
      formula$x
    else model.matrix(terms(formula), model.frame(formula))
    y <- if (is.vector(formula$y))
      formula$y
    else model.response(model.frame(formula))
    Z <- if (is.null(varformula))
      X
    else model.matrix(varformula, data = data)
    wts <- weights(formula)
  }
  else {
    mf <- if (is.null(weights)) {
      model.frame(formula, data = data)
    }
    else {
      model.frame(formula, weights = weights, data = data)
    }
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
    Z <- if (is.null(varformula))
      X
    else model.matrix(varformula, data = data)
    wts <- model.weights(mf)
  }
  if (is.null(wts))
    wts <- rep.int(1, NROW(X))
  if (!(all(c(row.names(X) %in% row.names(Z), row.names(Z) %in%
              row.names(X))))) {
    allnames <- row.names(X)[row.names(X) %in% row.names(Z)]
    X <- X[allnames, ]
    Z <- Z[allnames, ]
    y <- y[allnames]
    wts <- wts[row.names(X) %in% row.names(Z)]
  }
  if (ncol(Z) < 2L)
    stop("the auxiliary variance regression requires at least an intercept and a regressor")
  k <- ncol(X)
  n <- sum(wts > 0)
  resi <- lm.wfit(X, y, wts)$residuals
  sigma2 <- sum(wts * resi^2)/n
  {
    f <- resi^2/sigma2 - 1
    aux <- lm.wfit(Z, f, wts)
    bp <- 0.5 * sum(wts * aux$fitted.values^2)
    method <- "Breusch-Pagan test"
  }
  names(bp) <- "BP"
  df <- c(df = aux$rank - 1)
  RVAL <- list(statistic = bp, parameter = df, method = method,
               p.value = pchisq(bp, df, lower.tail = FALSE), data.name = dname)
  class(RVAL) <- "htest"
  return(RVAL)
}


#' White's Test for Heteroskedasticity in a Linear Regression Mode
#' @description This function implements the popular method of White (1980) for testing for heteroskedasticity in a linear regression model.
#' @param lmobj Either an object of class "lm" (e.g., generated by lm), or a list of two objects: a response vector and a design matrix. The objects are assumed to be in that order, unless they are given the names "X" and "y" to distinguish them. The design matrix passed in a list must begin with a column of ones if an intercept is to be included in the linear model. The design matrix passed in a list should not contain factors, as all columns are treated 'as is'. For tests that use ordinary least squares residuals, one can also pass a vector of residuals in the list, which should either be the third object or be named "e".
#' @param squares.only A logical. Should two-way interactions between explanatory variables be included in the auxiliary regression? Defaults to FALSE, since when interaction terms are present the test is not a pure test of heteroskedasticity but also of model specification.
#' @details White's Test entails fitting an auxiliary regression model in which the response variable is the vector of squared residuals from the original model and the design matrix includes the original explanatory variables, their squares, and (optionally) their two-way interactions. The test statistic is the number of observations multiplied by the coefficient of determination from the auxiliary regression model
#' @references Inspired by Gold Medalist Dr Shahla Faisal who teach me in GC University Faisalabad
#' @author Syed Hammad Raza
#' @export

white.test <- function(lmobj, squares.only=FALSE)
{
  stopifnot(class(lmobj)=='lm')
  mydata <- lmobj$model
  mydata[,1] <- lmobj$residual^2
  fml <- lmobj$call$formula
  formula1 <- paste(fml[2],fml[1],fml[3])
  pvs <- attr(lmobj$terms,"term.labels")
  k <- length(pvs);
  n <- length(lmobj$fit)

  for(i in 1:k){
    tmp <- NULL;
    if(substr(pvs[i],1,2)=="I("){
      tmp2 <- substr(pvs[i],3, nchar(pvs[i])-1);
    }else{
      tmp2 <- pvs[i];
    }
    for(j in 1:nchar(tmp2)){
      tmp1 <- substr(tmp2,j,j)
      if(tmp1 == ":")
        tmp <- paste(tmp, "*", sep='')
      else
        tmp <- paste(tmp, tmp1, sep='')
    }
    pvs[i] <- tmp
  }
  formula2 <- paste(fml[2],fml[1])
  for(i in 1:k){
    if(i>1)
      formula2 <- paste(formula2, "+", sep='')
    formula2 <- paste(formula2, "I(", pvs[i],")",sep='')
    if(squares.only){
      formula2 <- paste(formula2, "+I(", pvs[i],
                        "*", pvs[i], ")", sep = "")
    }else{
      for(j in i:k)
        formula2 <- paste(formula2,"+I(",pvs[i],
                          "*",pvs[j],")", sep='')
    }
  }

  method <- ifelse(squares.only,
                   "White test for constant variance, squares only",
                   "White test for constant variance")

  out <- lm(as.formula(formula2),data=mydata)
  if(summary(out)$r.squared == 1.0){
    RVAL <- NULL;
    warning("Test failed.  Possible reasons:\n\t (1) collinearity, or (2) sample size is not big enough for the White's test.");
  }else{
    LM = summary(out)$r.squared * n
    names(LM) <- "White"
    df <- out$rank - 1
    names(df) <- "df";
    RVAL <- list(statistic = LM,
                 parameter = df,
                 method = method,
                 p.value= pchisq(LM,df,lower.tail=FALSE),
                 data.name=NULL)
    class(RVAL) <- "htest"
  }
  return(RVAL)
}

#' Durbin-Watson Test
#' @description Performs the Durbin-Watson test for autocorrelation of disturbances.
#' @param formula a symbolic description for the model to be tested (or a fitted "lm" object).
#' @param order.by Either a vector z or a formula with a single explanatory variable like ~ z. The observations in the model are ordered by the size of z. If set to NULL (the default) the observations are assumed to be ordered (e.g., a time series).
#' @param alternative a character string specifying the alternative hypothesis.
#' @param iterations an integer specifying the number of iterations when calculating the p-value with the "pan" algorithm.
#' @param exact logical. If set to FALSE a normal approximation will be used to compute the p value, if TRUE the "pan" algorithm is used. The default is to use "pan" if the sample size is < 100.
#' @param tol tolerance. Eigenvalues computed have to be greater than tol to be treated as non-zero.
#' @param data an optional data frame containing the variables in the model. By default the variables are taken from the environment which dwtest is called from.
#' @details The Durbin-Watson test has the null hypothesis that the autocorrelation of the disturbances is 0. It is possible to test against the alternative that it is greater than, not equal to, or less than 0, respectively. This can be specified by the alternative argument.
#' @references Inspired by Gold Medalist Dr Shahla Faisal who teach me in GC University Faisalabad
#' @author Syed Hammad Raza
#' @export


dwtest = function (formula, order.by = NULL, alternative = c("greater",
                                                    "two.sided", "less"), iterations = 15, exact = NULL, tol = 1e-10,
          data = list())
{
  dname <- paste(deparse(substitute(formula)))
  alternative <- match.arg(alternative)
  if (!inherits(formula, "formula")) {
    if (!is.null(w <- weights(formula))) {
      if (!isTRUE(all.equal(as.vector(w), rep(1L, length(w)))))
        stop("weighted regressions are not supported")
    }
    X <- if (is.matrix(formula$x))
      formula$x
    else model.matrix(terms(formula), model.frame(formula))
    y <- if (is.vector(formula$y))
      formula$y
    else model.response(model.frame(formula))
  }
  else {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
  }
  if (!is.null(order.by)) {
    if (inherits(order.by, "formula")) {
      z <- model.matrix(order.by, data = data)
      z <- as.vector(z[, ncol(z)])
    }
    else {
      z <- order.by
    }
    X <- as.matrix(X[order(z), ])
    y <- y[order(z)]
  }
  n <- nrow(X)
  if (is.null(exact))
    exact <- (n < 100)
  k <- ncol(X)
  res <- lm.fit(X, y)$residuals
  dw <- sum(diff(res)^2)/sum(res^2)
  Q1 <- chol2inv(qr.R(qr(X)))
  if (n < 3) {
    warning("not enough observations for computing a p value, set to 1")
    pval <- 1
  }
  else {
    if (exact) {
      A <- diag(c(1, rep(2, n - 2), 1))
      A[abs(row(A) - col(A)) == 1] <- -1
      MA <- diag(rep(1, n)) - X %*% Q1 %*% t(X)
      MA <- MA %*% A
      ev <- eigen(MA)$values[1:(n - k)]
      if (any(Im(ev) > tol))
        warning("imaginary parts of eigenvalues discarded")
      ev <- Re(ev)
      ev <- ev[ev > tol]
      pdw <- function(dw) .Fortran("pan", as.double(c(dw,
                                                      ev)), as.integer(length(ev)), as.double(0),
                                   as.integer(iterations), x = double(1), PACKAGE = "lmtest")$x
      pval <- switch(alternative, two.sided = (2 * min(pdw(dw),
                                                       1 - pdw(dw))), less = (1 - pdw(dw)), greater = pdw(dw))
      if (is.na(pval) || ((pval > 1) | (pval < 0))) {
        warning("exact p value cannot be computed (not in [0,1]), approximate p value will be used")
        exact <- FALSE
      }
    }
    if (!exact) {
      if (n < max(5, k)) {
        warning("not enough observations for computing an approximate p value, set to 1")
        pval <- 1
      }
      else {
        AX <- matrix(as.vector(filter(X, c(-1, 2, -1))),
                     ncol = k)
        AX[1, ] <- X[1, ] - X[2, ]
        AX[n, ] <- X[n, ] - X[(n - 1), ]
        XAXQ <- t(X) %*% AX %*% Q1
        P <- 2 * (n - 1) - sum(diag(XAXQ))
        Q <- 2 * (3 * n - 4) - 2 * sum(diag(crossprod(AX) %*%
                                              Q1)) + sum(diag(XAXQ %*% XAXQ))
        dmean <- P/(n - k)
        dvar <- 2/((n - k) * (n - k + 2)) * (Q - P *
                                               dmean)
        pval <- switch(alternative, two.sided = (2 *
                                                   pnorm(abs(dw - dmean), sd = sqrt(dvar), lower.tail = FALSE)),
                       less = pnorm(dw, mean = dmean, sd = sqrt(dvar),
                                    lower.tail = FALSE), greater = pnorm(dw,
                                                                         mean = dmean, sd = sqrt(dvar)))
      }
    }
  }
  alternative <- switch(alternative, two.sided = "true autocorrelation is not 0",
                        less = "true autocorrelation is less than 0", greater = "true autocorrelation is greater than 0")
  names(dw) <- "DW"
  RVAL <- list(statistic = dw, method = "Durbin-Watson test",
               alternative = alternative, p.value = pval, data.name = dname)
  class(RVAL) <- "htest"
  return(RVAL)
}

#' Breusch-Godfrey Test
#' @description bgtest performs the Breusch-Godfrey test for higher-order serial correlation
#' @param formula a symbolic description for the model to be tested (or a fitted "lm" object).
#' @param order integer. maximal order of serial correlation to be tested.
#' @param order.by Either a vector z or a formula with a single explanatory variable like ~ z. The observations in the model are ordered by the size of z. If set to NULL (the default) the observations are assumed to be ordered (e.g., a time series).
#' @param type the type of test statistic to be returned. Either "Chisq" for the Chi-squared test statistic or "F" for the F test statistic
#' @param data an optional data frame containing the variables in the model. By default the variables are taken from the environment which bgtest is called from.
#' @param fill starting values for the lagged residuals in the auxiliary regression. By default 0 but can also be set to NA.
#' @details Under H_0 the test statistic is asymptotically Chi-squared with degrees of freedom as given in parameter. If type is set to "F" the function returns a finite sample version of the test statistic, employing an F distribution with degrees of freedom as given in parameter.
#' @references Inspired by Gold Medalist Dr Shahla Faisal who teach me in GC University Faisalabad
#' @author Syed Hammad Raza
#' @export


bgtest = function (formula, order = 1, order.by = NULL, type = c("Chisq",
                                                                 "F"), data = list(), fill = 0)
{
  dname <- paste(deparse(substitute(formula)))
  if (!inherits(formula, "formula")) {
    X <- if (is.matrix(formula$x))
      formula$x
    else model.matrix(terms(formula), model.frame(formula))
    y <- if (is.vector(formula$y))
      formula$y
    else model.response(model.frame(formula))
  }
  else {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
  }
  if (!is.null(order.by)) {
    if (inherits(order.by, "formula")) {
      z <- model.matrix(order.by, data = data)
      z <- as.vector(z[, ncol(z)])
    }
    else {
      z <- order.by
    }
    X <- as.matrix(X[order(z), ])
    y <- y[order(z)]
  }
  n <- nrow(X)
  k <- ncol(X)
  order <- 1:order
  m <- length(order)
  resi <- lm.fit(X, y)$residuals
  Z <- sapply(order, function(x) c(rep(fill, length.out = x),
                                   resi[1:(n - x)]))
  if (any(na <- !complete.cases(Z))) {
    X <- X[!na, , drop = FALSE]
    Z <- Z[!na, , drop = FALSE]
    y <- y[!na]
    resi <- resi[!na]
    n <- nrow(X)
  }
  auxfit <- lm.fit(cbind(X, Z), resi)
  cf <- auxfit$coefficients
  vc <- chol2inv(auxfit$qr$qr) * sum(auxfit$residuals^2)/auxfit$df.residual
  names(cf) <- colnames(vc) <- rownames(vc) <- c(colnames(X),
                                                 paste("lag(resid)", order, sep = "_"))
  switch(match.arg(type), Chisq = {
    bg <- n * sum(auxfit$fitted.values^2)/sum(resi^2)
    p.val <- pchisq(bg, m, lower.tail = FALSE)
    df <- m
    names(df) <- "df"
  }, F = {
    uresi <- auxfit$residuals
    bg <- ((sum(resi^2) - sum(uresi^2))/m)/(sum(uresi^2)/(n -
                                                            k - m))
    df <- c(m, n - k - m)
    names(df) <- c("df1", "df2")
    p.val <- pf(bg, df1 = df[1], df2 = df[2], lower.tail = FALSE)
  })
  names(bg) <- "LM test"
  RVAL <- list(statistic = bg, parameter = df, method = paste("Breusch-Godfrey test for serial correlation of order up to",
                                                              max(order)), p.value = p.val, data.name = dname, coefficients = cf,
               vcov = vc)
  class(RVAL) <- c("bgtest", "htest")
  return(RVAL)
}



