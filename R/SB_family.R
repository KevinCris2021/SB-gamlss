# R/SB_family.R

#' Internal helper: extract bd (or m) from ...
#' @keywords internal
.get_bd_sb <- function(...) {
  dots <- list(...)
  if (!is.null(dots$bd)) return(dots$bd)
  if (!is.null(dots$m))  return(dots$m)
  stop("No se encontró 'bd' (o 'm'). Pasa gamlss(..., bd = <vector>)")
}

#' Simplex-Binomial family for GAMLSS
#' @export
SB <- function(mu.link = "logit", sigma.link = "log") {

  if (!exists("dSB", mode = "function")) stop("No encuentro dSB().")
  if (!exists("pSB", mode = "function")) stop("No encuentro pSB().")

  if (!requireNamespace("gamlss.dist", quietly = TRUE)) stop("Falta 'gamlss.dist'")
  if (!requireNamespace("numDeriv", quietly = TRUE))    stop("Falta 'numDeriv'")

  recycle_all <- function(y, mu, sigma, bd) {
    n <- length(y)
    if (is.null(bd)) bd <- rep(max(y, na.rm = TRUE), n)
    list(
      y     = y,
      mu    = rep(mu,    length.out = n),
      sigma = rep(sigma, length.out = n),
      bd    = rep(bd,    length.out = n)
    )
  }

  clip_mu <- function(mu) pmin(pmax(mu, 1e-12), 1 - 1e-12)
  clip_sg <- function(sg) pmax(sg, 1e-12)

  mstats <- gamlss.dist::checklink("mu.link", "SB", substitute(mu.link),
                                  c("logit","probit","cloglog","log","own"))
  dstats <- gamlss.dist::checklink("sigma.link", "SB", substitute(sigma.link),
                                  c("log","identity","inverse"))

  mu.link    <- as.character(substitute(mu.link))
  sigma.link <- as.character(substitute(sigma.link))

  fam <- list(
    family     = c("SB","Simplex Binomial"),
    parameters = list(mu = TRUE, sigma = TRUE),
    nopar      = 2,
    type       = "Discrete",

    mu.link    = mu.link,
    sigma.link = sigma.link,

    mu.linkfun = function(mu)  mstats$linkfun(mu),
    mu.linkinv = function(eta) mstats$linkinv(eta),
    mu.dr      = function(eta) mstats$mu.eta(eta),

    sigma.linkfun = function(sigma) dstats$linkfun(sigma),
    sigma.linkinv = function(eta)   dstats$linkinv(eta),
    sigma.dr      = function(eta)   dstats$mu.eta(eta),

    # --- score wrt mu
    dldm = function(y, mu, sigma, ...) {
      bd <- SBgamlss:::.get_bd_sb(...)
      rr <- recycle_all(y, mu, sigma, bd)
      y <- rr$y; mu <- rr$mu; sigma <- rr$sigma; bd <- rr$bd

      vapply(seq_along(y), function(i) {
        yi <- y[i]; bd_i <- bd[i]
        mu_i <- clip_mu(mu[i])
        sg_i <- clip_sg(sigma[i])

        numDeriv::grad(function(m) {
          m <- clip_mu(m)
          as.numeric(dSB(yi, mu = m, sigma = sg_i, bd = bd_i, log = TRUE))
        }, mu_i)
      }, 0.0)
    },

    d2ldm2 = function(y, mu, sigma, ...) {
      bd <- SBgamlss:::.get_bd_sb(...)
      rr <- recycle_all(y, mu, sigma, bd)
      y <- rr$y; mu <- rr$mu; sigma <- rr$sigma; bd <- rr$bd

      vapply(seq_along(y), function(i) {
        yi <- y[i]; bd_i <- bd[i]
        mu_i <- clip_mu(mu[i])
        sg_i <- clip_sg(sigma[i])

        as.numeric(numDeriv::hessian(function(m) {
          m <- clip_mu(m)
          as.numeric(dSB(yi, mu = m, sigma = sg_i, bd = bd_i, log = TRUE))
        }, mu_i))
      }, 0.0)
    },

    # --- score wrt sigma
    dldd = function(y, mu, sigma, ...) {
      bd <- SBgamlss:::.get_bd_sb(...)
      rr <- recycle_all(y, mu, sigma, bd)
      y <- rr$y; mu <- rr$mu; sigma <- rr$sigma; bd <- rr$bd

      vapply(seq_along(y), function(i) {
        yi <- y[i]; bd_i <- bd[i]
        mu_i <- clip_mu(mu[i])
        sg_i <- clip_sg(sigma[i])

        numDeriv::grad(function(s) {
          s <- clip_sg(s)
          as.numeric(dSB(yi, mu = mu_i, sigma = s, bd = bd_i, log = TRUE))
        }, sg_i)
      }, 0.0)
    },

    d2ldd2 = function(y, mu, sigma, ...) {
      bd <- SBgamlss:::.get_bd_sb(...)
      rr <- recycle_all(y, mu, sigma, bd)
      y <- rr$y; mu <- rr$mu; sigma <- rr$sigma; bd <- rr$bd

      vapply(seq_along(y), function(i) {
        yi <- y[i]; bd_i <- bd[i]
        mu_i <- clip_mu(mu[i])
        sg_i <- clip_sg(sigma[i])

        as.numeric(numDeriv::hessian(function(s) {
          s <- clip_sg(s)
          as.numeric(dSB(yi, mu = mu_i, sigma = s, bd = bd_i, log = TRUE))
        }, sg_i))
      }, 0.0)
    },

    d2ldmdd = function(y, mu, sigma, ...) {
      bd <- SBgamlss:::.get_bd_sb(...)
      rr <- recycle_all(y, mu, sigma, bd)
      y <- rr$y; mu <- rr$mu; sigma <- rr$sigma; bd <- rr$bd

      vapply(seq_along(y), function(i) {
        yi <- y[i]; bd_i <- bd[i]
        mu_i <- clip_mu(mu[i])
        sg_i <- clip_sg(sigma[i])

        H <- numDeriv::hessian(function(v) {
          m <- clip_mu(v[1])
          s <- clip_sg(v[2])
          as.numeric(dSB(yi, mu = m, sigma = s, bd = bd_i, log = TRUE))
        }, c(mu_i, sg_i))

        as.numeric(H[1, 2])
      }, 0.0)
    },

    G.dev.incr = function(y, mu, sigma, ...) {
      bd <- SBgamlss:::.get_bd_sb(...)
      rr <- recycle_all(y, mu, sigma, bd)
      y <- rr$y; mu <- rr$mu; sigma <- rr$sigma; bd <- rr$bd

      -2 * vapply(seq_along(y), function(i) {
        as.numeric(dSB(y[i],
                      mu    = clip_mu(mu[i]),
                      sigma = clip_sg(sigma[i]),
                      bd    = bd[i],
                      log   = TRUE))
      }, 0.0)
    },

    rqres = expression({
      u1 <- pSB(y,     mu = mu, sigma = sigma, bd = bd)
      u0 <- pSB(y - 1, mu = mu, sigma = sigma, bd = bd)
      u  <- u0 + runif(length(y)) * (u1 - u0)
      u  <- pmin(pmax(u, 1e-12), 1 - 1e-12)
      qnorm(u)
    }),

    mu.initial = expression({
      n <- length(y)
      bd0 <- if (exists("bd", inherits = TRUE)) bd else NULL
      if (is.null(bd0) && exists("m", inherits = TRUE)) bd0 <- m
      if (is.null(bd0)) bd0 <- rep(max(y, na.rm = TRUE), n)
      bd0 <- rep(bd0, length.out = n)

      mu <- (y + 0.5) / (bd0 + 1)
      mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)
    }),

    sigma.initial = expression({ sigma <- rep(0.7, length(y)) }),

    mu.valid    = function(mu)    all(mu > 0 & mu < 1),
    sigma.valid = function(sigma) all(sigma > 0),
    y.valid     = function(y)     all(y >= 0)
  )

  class(fam) <- c("gamlss.family","family")
  fam
}
