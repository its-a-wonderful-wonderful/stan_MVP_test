library(rstan) #blah
library(coda)

# like sum() but with less cancelation error
pairwise_sum <- function(x, N = 2) {
  signs <- x >= 0
  if(length(unique(signs)) == 1) return(sum(x))
  else if(all(abs(x) < 1)) return(sum(x))
  else {
    m <- floor(length(x) / 2)
    return(pairwise_sum(x[1:m]) + pairwise_sum(x[-c(1:m)]))
  }
}

# This function is copied from that in the energy package and licensed under GPL 2+
# Copyright Maria Rizzo 2014 
Astar <- function(d) {
  d <- as.matrix(d)
  n <- nrow(d)
  if (n != ncol(d)) stop("Argument d should be square")
  m <- rowMeans(d)
  M <- mean(d)
  a <- sweep(d, 1, m)
  b <- sweep(a, 2, m)
  A <- b + M
  A <- A - d/n
  diag(A) <- m - M
  (n/(n - 1)) * A  
}

# bias-corrected distance correlation based on energy:::BCDCOR
bcdcor <- function(lag, Astars, Astars2, XX) {
  constant <- nrow(Astars[[1]]) / (nrow(Astars[[1]]) - 2)
  oneplus <- 1 + .Machine$double.eps
  sapply(1:length(Astars), FUN = function(i) {
    if(lag >= i) return(NA_real_)
    AABB <- Astars[[i]] * Astars[[i - lag]]
    AAAA <- Astars2[[i]]
    BBBB <- Astars2[[i - lag]]
    XY <- pairwise_sum(AABB) - constant * pairwise_sum(diag(AABB))
    XXYY <- XX[i] * XX[i - lag]
    ifelse(XXYY > 0, XY/sqrt(XXYY), ifelse(XXYY == 0, 0, oneplus))              
  })
}

setGeneric("dependence", function(object, ...) standardGeneric("dependence"))

setMethod("dependence", "array",
          function(object, lag = 1) {
            pars <- dim(object)[3]
            if(is.na(pars)) stop("'object' must be an array with 3 dimensions")
            chains <- ncol(object)
            if(chains <= 2) stop("must have at least 3 chains")
            constant <- chains / (chains - 2)
            sims <- nrow(object)
            Astars <- lapply(1:sims, FUN = function(i) Astar(dist(object[i,,])))
            Astars2 <- lapply(Astars, FUN = "^", e2 = 2)
            XX <- sapply(Astars2, FUN = function(AAAA) {
              pairwise_sum(AAAA) - constant * pairwise_sum(diag(AAAA))
            })
            bcdcor(lag, Astars, Astars2, XX)
})

get_stan_params <- function(object) {
  stopifnot(is(object, "stanfit"))
  params <- grep("context__.vals_r", fixed = TRUE, value = TRUE,
                 x = strsplit(get_cppcode(get_stanmodel(object)), "\n")[[1]])
  params <- sapply(strsplit(params, "\""), FUN = function(x) x[[2]])
  params <- intersect(params, object@sim$pars_oi)
  stopifnot(length(params) > 0)
  return(params)
}

setMethod("dependence", "stanfit",
          function(object, lag = 1:7,
                   plot = TRUE, pars = get_stan_params(object), f = 1/ncol(object), 
                   legend_location = "topright", ...) {
            stopifnot(object@mode == 0)
            chains <- ncol(object)
            if(chains <= 2) stop("must have at least 3 chains")
            constant <- chains / (chains - 2)
            sims <- nrow(object)            
            x <- rstan::extract(object, pars, permuted = FALSE, inc_warmup = TRUE)
            Astars <- lapply(1:nrow(x), FUN = function(i) Astar(dist(x[i,,])))
            Astars2 <- lapply(Astars, FUN = "^", e2 = 2)
            XX <- sapply(Astars2, FUN = function(AAAA) {
              pairwise_sum(AAAA) - constant * pairwise_sum(diag(AAAA))
            })
            dcors <- sapply(lag, FUN = function(l) bcdcor(l, Astars, Astars2, XX))
            if(length(lag) > 1) colnames(dcors) <- paste0("lag_", lag)
            warmup <- object@stan_args[[1]]$warmup / object@stan_args[[1]]$thin
            divergence <- unlist(sapply(get_sampler_params(object), 
                                 FUN = function(y) which(y[,"n_divergent__"] > 0)))
            M <- (chains * (chains - 3)) %/% 2
            df <- M - 1L
            q <- qt(0.975, df)
            cv <- q / sqrt(q^2 + df - 1)
            if(plot && length(lag) == 1) {
              lag <- lag * object@stan_args[[1]]$thin
              x <- dcors[,1]
              x[(1:nrow(dcors) %% lag) != 0] <- NA_real_
              plot(x, pch = ".", ylim = c(-cv,1), las = 1, ##log = "x",
                   xlab = "Iteration",
                   ylab = paste("Estimated distance correlation among", 
                                ncol(object), "chains"),           
                   main = paste(object@model_name, "model, lag =", lag), ...)
              lines(lowess(na.omit(dcors), f = f), type = "l", col = 2, lty = 2)
              abline(v = warmup, col = 3, lty = 3)
              abline(h = cv, col = 4, lty = 3)
              abline(h = -cv, col = 4, lty = 3)              
              legend(legend_location, legend = c("Value", "Smoothed", "95% CI"),
                     col = c(1:2, 4), lty = c(NA_integer_, 2:3), 
                     pch = c(20L, NA_integer_, NA_integer_), cex = 0.75)
              text(x = warmup, y = 1.02, labels = "Warmup  | Retained")
              rug(divergence, ticksize = -0.01)
            }
            else if(plot) {
              plot(dcors[,1], type = "n", ylim = 0:1, las = 1, ##log = "x",
                   xlab = "Iteration", 
                   ylab = paste("Estimated distance correlation among", 
                                ncol(object), "chains"),
                   main = paste(object@model_name, "model"), ...)
              sapply(1:ncol(dcors), FUN = function(i) {
                horiz <- which(!is.na(dcors[,i]))
                lines(horiz, lowess(dcors[horiz,i], f = f)$y, col = i + 1, ...)
              })
              abline(v = warmup, col = 1, lty = 3)
              legend(legend_location, legend = lag * object@stan_args[[1]]$thin, 
                     title = "lag length", ncol = 2,
                     col = 1 + 1:length(lag), lty = 1, cex = 0.75)
              text(x = warmup, y = 1.02, labels = "Warmup  | Retained")
              rug(divergence, ticksize = -0.01)
              abline(h = 0, lty = "dotted")
            }
            return(invisible(dcors))
          }
)

setOldClass("mcmc.list")

# The setMethod function takes three arguments: 
#   1) the name of the generic function, 
#   2) the signature to match for this method and 
#   3) a function to compute the result. 

setMethod("dependence", "mcmc.list",
          function(object, lag = 1:7, 
                   plot = TRUE, pars = colnames(object[[1]]), f = 1/length(object), 
                   legend_location = "topright", ...) {
            x <- as.array(object)
            if(is.matrix(x)) x <- array(x, dim = c(dim(x), 1))
            else x <- x[,pars,,drop=FALSE]
            if(dim(x)[3] == 1) Astars <- lapply(1:nrow(x), FUN = function(i) Astar(dist(x[i,,])))
            else Astars <- lapply(1:nrow(x), FUN = function(i) Astar(dist(t(x[i,,]))))
            Astars2 <- lapply(Astars, FUN = "^", e2 = 2)
            chains <- length(object)
            constant <- chains / (chains - 2)
            XX <- sapply(Astars2, FUN = function(AAAA) {
              pairwise_sum(AAAA) - constant * pairwise_sum(diag(AAAA))
            })
            dcors <- sapply(lag, FUN = function(l) bcdcor(l, Astars, Astars2, XX))            
            burnin <- attributes(object[[1]])$mcpar[1] - 1L
            title <- attributes(object[[1]])$title
            thin <- attributes(object[[1]])$mcpar[3]
            M <- (chains * (chains - 3)) %/% 2
            df <- M - 1L
            q <- qt(0.975, df)
            cv <- q / sqrt(q^2 + df - 1)            
            if(plot && length(lag) == 1) {
              x <- dcors[,1]
              x[(1:nrow(dcors) %% lag) != 0] <- NA_real_              
              plot(x, pch = ".", ylim = c(-cv,1), las = 1, ##log = "x",
                   xlab = "Iteration",
                   ylab = paste("Estimated distance correlation among", 
                                ncol(object), "chains"),           
                   main = paste(title, "lag =", lag * thin), ...)
              lines(lowess(na.omit(dcors), f = f), type = "l", col = 2, lty = 2)
              abline(h = cv, col = 4, lty = 3)
              abline(h = -cv, col = 4, lty = 3)              
              legend(legend_location, legend = c("Value", "Smoothed", "95% CI"),
                     col = c(1:2, 4), lty = c(NA_integer_, 2:3), 
                     pch = c(20L, NA_integer_, NA_integer_), cex = 0.75)
            }
            else if(plot) {
              plot(dcors[,1], type = "n", ylim = 0:1, las = 1, #log = "x",
                   xlab = "Iteration", 
                   ylab = paste("Estimated distance correlation among", 
                                chains, "chains"),
                   main = paste(title), ...)
              sapply(1:ncol(dcors), FUN = function(i) {
                horiz <- which(!is.na(dcors[,i]))
                lines(horiz, lowess(dcors[horiz,i], f = f)$y, col = i + 1, ...)
              })
              legend(legend_location, legend = lag * thin, title = "lag length", 
                     ncol = 2, col = 1 + 1:(1 + length(lag)), lty = 1, cex = 0.75)
              abline(h = 0, lty = "dotted")
            }
            return(invisible(dcors))
          }
)
