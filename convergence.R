library(rstan) # another commenting commit test thing
library(energy)
library(boot)
library(parallel)


options(boot.parallel = "multicore")
options(boot.ncpus = detectCores())

# derived from energy:::disco() which is under the GPL2+
disco_mcmc <- function (x, chain_id, split_id, index, R, fast = TRUE) {
  fast_disco1 <- function(trt, dst, total, index) {
    trt <- factor(trt)
    k <- nlevels(trt)
    n <- tabulate(trt)
    N <- sum(n)
    M <- model.matrix(~trt - 1) == 1
    if(is.matrix(dst)) {
      G <- apply(M, 2, FUN = function(x) sum(dst[x,x])) # x %*% dst %*% x
      W <- sum(G / (2 * n))
    }
    else { # AnnoyEuclidean
      W <- 0
      for(k in 1:ncol(M)) {
        x <- which(M[,k]) - 1L
        eg <- expand.grid(i = x, j = x, KEEP.OUT.ATTRS = FALSE)
        i <- eg$i
        j <- eg$j
        keep <- i < j
        i <- i[keep]
        j <- j[keep]
        W <- W + dst$getSumDistanceByVectors(i, j, index / 2) / n[k]
      }
    }
    B <- total - W
    c(B, W, k - 1, N - k)
  }
  fast_disco2 <- function(trt, dst, first, second, total) {
    both <- fast_disco1(trt, dst, total)
    fst <- fast_disco1(first, dst, total)
    snd <- fast_disco1(second, dst, total)
    df1 <- fst[3] * snd[3]
    c(B1 = fst[1], B2 = snd[1], B12 = both[1] - fst[1] - snd[1], 
      W1 = fst[2], W2 = snd[2], W12 = both[2],
      df1_1 = fst[3], df1_2 = snd[3], df1_12 = df1,
      df2_1 = fst[4], df2_2 = snd[4], df2_12 = NROW(trt) - df1)
  }  
  fast_disco2stat <- function(data, i, dst, trt, first, second, total) {
    idx <- seq_len(NROW(i))[i]
    d <- fast_disco2(trt[idx], dst, first[idx], second[idx], total)
    B1 <- (d["B1"] / d["df1_1"])
    B2 <- (d["B2"] / d["df1_2"])
    B12 <- (d["B12"] / d["df1_12"])
    W12 <- (d["W12"] / d["df2_12"])
    return(c(F1 = B1 / W12, F2 = B2 / W12, F12 = B12 / W12))
  }
  nfactors <- 3L
  N <- nrow(x)
  if(fast) {
    dst <- if(index != 1) dist(x)^index else dist(x)
    dst <- as.matrix(dst)
    total <- sum(dst) / (2 * nrow(x))
  }
  else {
    stop("fast = FALSE is not fully implemented yet")
    stopifnot(require(RcppAnnoy))
    dst <- new(AnnoyEuclidean, ncol(x))
    on.exit(dst$unload())
    for(i in 1:nrow(x)) dst$addItem(i - 1, x[i,])
    combos <- combn(0:(nrow(x) - 1L), 2)
    total <- dst$getSumDistanceByVectors(combos[1,], combos[2,], index / 2) / nrow(x)
  }
  formals(fast_disco1)$index <- index
  stats <- matrix(NA_real_, nfactors, 6)
  colnames(stats) <- c("Trt", "Within", "df1", "df2", "Stat", "p-value")
  
  # "ANOVA" by chain and split and their interaction
  trt <- as.factor(paste0(chain_id, split_id))
  d <- fast_disco2(trt = trt, dst = dst, first = chain_id, second = split_id, total = total)
  stats[1, 1:4] <- d[c("B1", "W1", "df1_1", "df2_1")]
  stats[2, 1:4] <- d[c("B2", "W2", "df1_2", "df2_2")]
  stats[3, 1:4] <- d[c("B12", "W12", "df1_12", "df2_12")]
  if(R > 0) {
    b <- boot(data = 1:N, statistic = fast_disco2stat, sim = "permutation",
              R = R, dst = dst, trt = trt, first = chain_id, second = split_id, total = total)
    stats[1:3, 5] <- b$t0
    stats[1:3, 6] <- (colSums( sweep(b$t, 2, b$t0) > 0 ) + 1) / (R + 1)  
  }
  else stats[3,5] <- fast_disco2stat(dst, i = 1:N, trt = trt, chain_id, split_id, total)["F12"]
  
  methodname <- "DISCO"
  dataname <- "posterior"
  total <- sum(stats[1, 1:2])
  within <- total - sum(stats[, 1])
  Df.trt <- stats[, 3]
  factor.names <- c("chain", "split", "chain:split")
  factor.levels <- c(nlevels(chain_id), 2, nlevels(trt))
  sizes <- c(tabulate(chain_id), tabulate(split_id), tabulate(trt))
  e <- list(call = match.call(), method = methodname, 
            statistic = stats[,5], p.value = stats[, 6], k = nfactors, N = N, 
            between = stats[, 1], withins = stats[, 2], within = within, total = total, 
            Df.trt = Df.trt, Df.e = nrow(dst) - sum(Df.trt) - 1, 
            index = index, factor.names = factor.names, factor.levels = factor.levels, 
            sample.sizes = sizes, stats = stats)
  class(e) <- "disco_mcmc"
  e
}

print.disco_mcmc <- function(x, ...) {
  k <- x$k
  md1 <- x$between/x$Df.trt
  md2 <- x$within/x$Df.e
  f0 <- x$statistic
  print(x$call)
  cat(sprintf("\nDistance Components: index %5.2f\n", x$index))
  cat(sprintf("%-20s %4s %10s %10s %10s %10s\n", "Source", 
              "Df", "Sum Dist", "Mean Dist", "B / W", "p-value"))
  for (i in 1:k) {
    fname <- x$factor.names[i]
    cat(sprintf("%-20s %4d %10.5f %10.5f %10.5f %10s\n", 
                fname, x$Df.trt[i], x$between[i], md1[i], f0[i], 
                format.pval(x$p.value[i])))
  }
  cat(sprintf("%-20s %4d %10.5f %10.5f\n", "Within", x$Df.e, 
              x$within, md2))
  cat(sprintf("%-20s %4d %10.5f\n", "Total", x$N - 1, x$total))
  return(invisible(NULL))
}

setGeneric("convergence_test", 
           function(object, ...) standardGeneric("convergence_test"))
setMethod("convergence_test", "array",
          function(object, thin = NA_integer_, index = 1, R = NA_integer_, fast = TRUE) {

  if(length(dim(object)) != 3) stop("'object' must be a 3 dimensional array")
  
  iterations <- nrow(object)
  chains <- ncol(object)
  stopifnot(is.numeric(thin), thin == as.integer(thin), thin > 0, thin < iterations)
  if((iterations %% thin) != 0) stop("'thin' value leaves a remainder")
  stopifnot(is.numeric(R), R == as.integer(R), R >= 0)  
  object <- object[(1:iterations) %% thin == 0,,,drop = FALSE]
  posterior <- apply(object, 3, FUN = function(x) x)
  thin_iterations <- iterations /thin
  chain_id <- as.factor(rep(1:chains, each = thin_iterations))
  split_id <- as.factor(rep(c(rep(1, floor(0.5 * thin_iterations)),
                              rep(2, ceiling(0.5 * thin_iterations))), times = chains))
  test <- disco_mcmc(posterior, chain_id, split_id, index = index, R = R, fast = fast)
  return(test)
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

setMethod("convergence_test", "stanfit",
          function(object, thin = NA_integer_, index = 1, R = NA_integer_, fast = TRUE) {
  convergence_test(rstan::extract(object, pars = get_stan_params(object), permuted = FALSE, inc_warmup = FALSE), 
                   thin = thin, index = index, R = R, fast = fast)          
})

setOldClass("mcmc.list")
setMethod("convergence_test", "mcmc.list",
          function(object, thin = NA_integer_, index = 1, R = NA_integer_, fast = TRUE) {
  x <- as.array(object)
  if(is.matrix(x)) x <- array(x, c(dim(x), 1))
  else x <- aperm(x, c(1,3,2))
  convergence_test(x, thin = thin, index = index, R = R, fast = fast)
})
