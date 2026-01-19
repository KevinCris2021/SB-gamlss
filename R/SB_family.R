.get_bd_sb <- function(...) {
  dots <- list(...)

  if (!is.null(dots$bd)) return(dots$bd)
  if (!is.null(dots$m))  return(dots$m)

  stop("No se encontró 'bd' (o 'm'). Pasa gamlss(..., bd = <vector>)")
}


SB <- function(mu.link = "logit", sigma.link = "log") {

  if (!exists("dSB")) stop("Define primero dSB() (wrapper a tu dsb_vec).")
  if (!exists("pSB")) stop("Define primero pSB().")

  mstats <- gamlss.dist::checklink("mu.link","SB",substitute(mu.link),
                                   c("logit","probit","cloglog","log","own"))
  dstats <- gamlss.dist::checklink("sigma.link","SB",substitute(sigma.link),
                                   c("log","identity","inverse"))

  mu.link    <- as.character(substitute(mu.link))
  sigma.link <- as.character(substitute(sigma.link))

  list(
    family     = c("SB","Simplex Binomial"),
    parameters = list(mu=TRUE, sigma=TRUE),
    nopar      = 2, type = "Discrete",

    mu.link    = mu.link,
    sigma.link = sigma.link,
   mu.linkfun    = function(mu)   mstats$linkfun(mu),
    mu.linkinv    = function(eta)  mstats$linkinv(eta),
    mu.dr         = function(eta)  mstats$mu.eta(eta),
    
    sigma.linkfun = function(sigma) dstats$linkfun(sigma),
    sigma.linkinv = function(eta)   dstats$linkinv(eta),    
    sigma.dr      = function(eta)   dstats$mu.eta(eta),


    ## Derivadas numéricas consistentes con dSB()
    dldm = function(y, mu, sigma, ...) {
  bd <- .get_bd_sb(...)
  vapply(seq_along(y), function(i){
    mu_i <- pmin(pmax(mu[i], 1e-12), 1-1e-12)
    sg_i <- pmax(sigma[i], 1e-12)
    bd_i <- bd[i]; yi <- y[i]
    numDeriv::grad(function(m){
      m <- pmin(pmax(m, 1e-12), 1-1e-12)
      as.numeric(dSB(yi, mu = m, sigma = sg_i, bd = bd_i, log = TRUE))
    }, mu_i)
  }, 0.0)
},

    d2ldm2 = function(y, mu, sigma, ...) {
  bd <- .get_bd_sb(...)
  vapply(seq_along(y), function(i){
    mu_i <- pmin(pmax(mu[i], 1e-12), 1-1e-12)
    sg_i <- pmax(sigma[i], 1e-12)
    bd_i <- bd[i]; yi <- y[i]
    as.numeric(numDeriv::hessian(function(m){
      m <- pmin(pmax(m, 1e-12), 1-1e-12)
      as.numeric(dSB(yi, mu = m, sigma = sg_i, bd = bd_i, log = TRUE))
    }, mu_i))
  }, 0.0)
},

    dldd = function(y, mu, sigma, ...) {
  bd <- .get_bd_sb(...)
  vapply(seq_along(y), function(i){
    mu_i <- pmin(pmax(mu[i], 1e-12), 1-1e-12)
    sg_i <- pmax(sigma[i], 1e-12)
    bd_i <- bd[i]; yi <- y[i]
    numDeriv::grad(function(s){
      s <- pmax(s, 1e-12)
      as.numeric(dSB(yi, mu = mu_i, sigma = s, bd = bd_i, log = TRUE))
    }, sg_i)
  }, 0.0)
},

    d2ldd2 = function(y, mu, sigma, ...) {
  bd <- .get_bd_sb(...)
  vapply(seq_along(y), function(i){
    mu_i <- pmin(pmax(mu[i], 1e-12), 1-1e-12)
    sg_i <- pmax(sigma[i], 1e-12)
    bd_i <- bd[i]; yi <- y[i]
    as.numeric(numDeriv::hessian(function(s){
      s <- pmax(s, 1e-12)
      as.numeric(dSB(yi, mu = mu_i, sigma = s, bd = bd_i, log = TRUE))
    }, sg_i))
  }, 0.0)
},

  d2ldmdd = function(y, mu, sigma, ...) {
  bd <- .get_bd_sb(...)
  vapply(seq_along(y), function(i){
    mu_i <- pmin(pmax(mu[i], 1e-12), 1-1e-12)
    sg_i <- pmax(sigma[i], 1e-12)
    bd_i <- bd[i]; yi <- y[i]
    H <- numDeriv::hessian(function(v){
      m <- pmin(pmax(v[1], 1e-12), 1-1e-12)
      s <- pmax(v[2], 1e-12)
      as.numeric(dSB(yi, mu = m, sigma = s, bd = bd_i, log = TRUE))
    }, c(mu_i, sg_i))
    as.numeric(H[1, 2])
  }, 0.0)
},


    ## Devianza coherente: -2 sum log f_i
      G.dev.incr = function(y, mu, sigma, ...) {
    bd <- .get_bd_sb(...)
    -2 * vapply(seq_along(y), function(i){
      m  <- pmin(pmax(mu[i], 1e-12), 1-1e-12)
      sg <- pmax(sigma[i], 1e-12)
      as.numeric(dSB(y[i], mu = m, sigma = sg, bd = bd[i], log = TRUE))
    }, 0.0)
  },


    ## RQR: usa el nombre "pSB" (rqres la resuelve por nombre)
rqres = expression({
  # Randomized quantile residuals for discrete Y
  u1 <- pSB(y,     mu = mu, sigma = sigma, bd = bd, lower.tail = TRUE, log.p = FALSE)
  u0 <- pSB(y - 1, mu = mu, sigma = sigma, bd = bd, lower.tail = TRUE, log.p = FALSE)

  # Uniformización: U ~ Unif(F(y-1), F(y))
  u  <- u0 + runif(length(y)) * (u1 - u0)

  # Guardas numéricas
  u  <- pmin(pmax(u, 1e-12), 1 - 1e-12)

  # Residuo cuantilico randomizado
  qnorm(u)
}),



    ## Iniciales y dominios
  mu.initial = expression({
    mu <- pmin(pmax((y + 0.5) / (bd + 1), 1e-6), 1 - 1e-6)
  }),


    sigma.initial = expression({ sigma <- rep(0.7, length(y)) }),
    mu.valid      = function(mu)    all(mu>0 & mu<1),
    sigma.valid   = function(sigma) all(sigma>0),
    y.valid       = function(y)     all(y>=0)
  ) -> fam

  class(fam) <- c("gamlss.family","family")
  fam
}
