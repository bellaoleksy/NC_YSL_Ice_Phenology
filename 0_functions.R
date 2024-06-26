#Function for calculating water-year DOY. This will help facilitate plotting and analysizing trends in ice-in since they span either side of the winter-year (e.g., 2011-2012). For example, an IceInDayofYear_fed value of 150 means Ice-In occured 150 days after the start of the water-year (Oct1)

hydro.day = function(x, start.month = 10L) {
  start.yr = year(x) - (month(x) < start.month)
  start.date = make_date(start.yr, start.month, 1L)
  as.integer(x - start.date + 1L)
}


Deriv <- function(mod,
                  n = 200,
                  eps = 1e-7,
                  newdata,
                  term) {
  if (inherits(mod, "gamm"))
    mod <- mod$gam
  m.terms <- attr(terms(mod), "term.labels")
  if (missing(newdata)) {
    newD <- sapply(model.frame(mod)[, m.terms, drop = FALSE],
                   function(x)
                     seq(min(x), max(x), length = n))
    names(newD) <- m.terms
  } else {
    newD <- newdata
  }
  X0 <- predict(mod, data.frame(newD), type = "lpmatrix")
  newD <- newD + eps
  X1 <- predict(mod, data.frame(newD), type = "lpmatrix")
  Xp <- (X1 - X0) / eps
  Xp.r <- NROW(Xp)
  Xp.c <- NCOL(Xp)
  ## dims of bs
  bs.dims <- sapply(mod$smooth, "[[", "bs.dim") - 1
  ## number of smooth terms
  t.labs <- attr(mod$terms, "term.labels")
  ## match the term with the the terms in the model
  if (!missing(term)) {
    want <- grep(term, t.labs)
    if (!identical(length(want), length(term)))
      stop("One or more 'term's not found in model!")
    t.labs <- t.labs[want]
  }
  nt <- length(t.labs)
  ## list to hold the derivatives
  lD <- vector(mode = "list", length = nt)
  names(lD) <- t.labs
  for (i in seq_len(nt)) {
    Xi <- Xp * 0
    want <- grep(t.labs[i], colnames(X1))
    Xi[, want] <- Xp[, want]
    df <- Xi %*% coef(mod)
    df.sd <- rowSums(Xi %*% mod$Vp * Xi) ^ .5
    lD[[i]] <- list(deriv = df, se.deriv = df.sd)
  }
  class(lD) <- "Deriv"
  lD$gamModel <- mod
  lD$eps <- eps
  lD$eval <- newD - eps
  lD ##return
}

confint.Deriv <- function(object, term, alpha = 0.05, ...) {
  l <- length(object) - 3
  term.labs <- names(object[seq_len(l)])
  if (missing(term)) {
    term <- term.labs
  } else {
    ## how many attempts to get this right!?!?
    ##term <- match(term, term.labs)
    ##term <- term[match(term, term.labs)]
    term <- term.labs[match(term, term.labs)]
  }
  if (any(miss <- is.na(term)))
    stop(paste("'term'", term[miss], "not a valid model term."))
  res <- vector(mode = "list", length = length(term))
  names(res) <- term
  residual.df <- df.residual(object$gamModel)
  tVal <- qt(1 - (alpha / 2), residual.df)
  ##for(i in term.labs[term]) {
  for (i in term) {
    upr <- object[[i]]$deriv + tVal * object[[i]]$se.deriv
    lwr <- object[[i]]$deriv - tVal * object[[i]]$se.deriv
    res[[i]] <- list(upper = drop(upr), lower = drop(lwr))
  }
  res$alpha = alpha
  res
}

signifD <- function(x, d, upper, lower, eval = 0) {
  miss <- upper > eval & lower < eval
  incr <- decr <- x
  want <- d > eval
  incr[!want | miss] <- NA
  want <- d < eval
  decr[!want | miss] <- NA
  list(incr = incr, decr = decr)
}

plot.Deriv <- function(x,
                       alpha = 0.05,
                       polygon = TRUE,
                       sizer = FALSE,
                       term,
                       eval = 0,
                       lwd = 3,
                       col = "lightgrey",
                       border = col,
                       ylab,
                       xlab,
                       main,
                       ...) {
  l <- length(x) - 3
  ## get terms and check specified (if any) are in model
  term.labs <- names(x[seq_len(l)])
  if (missing(term)) {
    term <- term.labs
  } else {
    term <- term.labs[match(term, term.labs)]
  }
  if (any(miss <- is.na(term)))
    stop(paste("'term'", term[miss], "not a valid model term."))
  if (all(miss))
    stop("All terms in 'term' not found in model.")
  l <- sum(!miss)
  nplt <- n2mfrow(l)
  tVal <- qt(1 - (alpha / 2), df.residual(x$gamModel))
  if (missing(ylab))
    ylab <- expression(italic(hat(f) * "'" * (x)))
  if (missing(xlab)) {
    xlab <- attr(terms(x$gamModel), "term.labels")
    names(xlab) <- xlab
  }
  if (missing(main)) {
    main <- term
    names(main) <- term
  }
  ## compute confidence interval
  CI <- confint(x, term = term)
  ## plots
  layout(matrix(seq_len(l), nrow = nplt[1], ncol = nplt[2]))
  for (i in term) {
    upr <- CI[[i]]$upper
    lwr <- CI[[i]]$lower
    ylim <- range(upr, lwr)
    plot(
      x$eval[, i],
      x[[i]]$deriv,
      type = "n",
      ylim = ylim,
      ylab = ylab,
      xlab = xlab[i],
      main = main[i],
      ...
    )
    if (isTRUE(polygon)) {
      polygon(c(x$eval[, i], rev(x$eval[, i])),
              c(upr, rev(lwr)),
              col = col,
              border = border)
    } else {
      lines(x$eval[, i], upr, lty = "dashed")
      lines(x$eval[, i], lwr, lty = "dashed")
    }
    abline(h = 0, ...)
    if (isTRUE(sizer)) {
      lines(x$eval[, i], x[[i]]$deriv, lwd = 1)
      S <- signifD(x[[i]]$deriv, x[[i]]$deriv, upr, lwr,
                   eval = eval)
      lines(x$eval[, i], S$incr, lwd = lwd, col = "blue")
      lines(x$eval[, i], S$decr, lwd = lwd, col = "red")
    } else {
      lines(x$eval[, i], x[[i]]$deriv, lwd = 2)
    }
  }
  layout(1)
  invisible(x)
}


# Gam model function ####
# Fit a GAM with specified `gam_formula`  
# Make the following objects  
# model summaries  
# draw_plot  
# appraise model  
# Return a list with the following:  
#  `lake_name`  
# `mod` = model fit  
# `draw_plot` = results of the `draw()` function  
# `appraise` = results of the `apraise()` function  

gam_mod_fun <- function(data, gam_formula, lake_name){
  mod <- gam(gam_formula,
             family=Gamma(link="log"),
             data = data %>%
               filter(lake == lake_name),
             correlation = corCAR1(
               form = ~ start_year),
             method = "REML")
  mod_summary <- summary(mod)
  draw_plot <- draw(mod)
  app_plot <- appraise(mod)
  return(list(lake_name = lake_name,
              mod = mod,
              summary = mod_summary,
              draw_plot = draw_plot,
              appraise = app_plot))}


# gam summary function ####
# This function summarizes the fits from `gam_mod_fun()` above 
# and displays output in a "tidy" format

gam_summary <- function(mydata, gam_formula, group_var) {
  group_var = enquo(group_var)
  
  mydata <- mydata %>% 
    group_by(!!group_var) %>% 
    nest() 
  
  mydata %>% 
    mutate(
      model = map(data,
                  ~(gam(
                    gam_formula,
                    family=Gamma(link="log"),
                    data = .,
                    correlation = corCAR1(
                      form = ~ start_year),
                    method = "REML")))) 
}



# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
#Function from: http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  = (cormat)[ut],
    p = pmat[ut]
  )
}