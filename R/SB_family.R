SB <- function(mu.link = "logit", sigma.link = "log") {

  # --- Bindings robustos: siempre resolver desde el namespace del paquete ---
  .get_bd_sb  <- getFromNamespace(".get_bd_sb",  "SBgamlss")
  .recycle_sb <- getFromNamespace(".recycle_sb", "SBgamlss")
  .clip_mu    <- getFromNamespace(".clip_mu",    "SBgamlss")
  .clip_sg    <- getFromNamespace(".clip_sg",    "SBgamlss")

  # --- Chequeos: usa exists() en el namespace del paquete para evitar falsos negativos ---
  if (!exists("dSB", envir = asNamespace("SBgamlss"), inherits = FALSE))
    stop("No se encontró dSB() en el namespace del paquete. Reinstala SBgamlss.")
  if (!exists("pSB", envir = asNamespace("SBgamlss"), inherits = FALSE))
    stop("No se encontró pSB() en el namespace del paquete. Reinstala SBgamlss.")

  # Accede a dSB/pSB del paquete (evita tomar funciones con el mismo nombre en .GlobalEnv)
  dSB_fun <- getFromNamespace("dSB", "SBgamlss")
  pSB_fun <- getFromNamespace("pSB", "SBgamlss")

  mstats <- gamlss.dist::checklink(
    "mu.link", "SB", substitute(mu.link),
    c("logit","probit","cloglog","log","own")
  )
  dstats <- gamlss.dist::checklink(
    "sigma.link", "SB", substitute(sigma.link),
    c("log","identity","inverse")
  )

  mu.link    <- as.character(substitute(mu.link))
  sigma.link <- as.character(substitute(sigma.link))

  fam <- list(
    family     = c("SB","Simplex Binomial"),
    parameters = list(mu = TRUE, sigma = TRUE),
    nopar      = 2,
    type       = "Discrete",

    mu.link    = mu.link,
    sigma.link = sigma.link,

    mu.linkfun = function(mu)   mstats$linkfun(mu),
    mu.linkinv = function(eta)  mstats$linkinv(eta),
    mu.dr      = function(eta)  mstats$mu.eta(eta),

    sigma.linkfun = function(sigma) dstats$linkfun(sigma),
    sigma.linkinv = function(eta)   dstats$linkinv(eta),
    sigma.dr      = function(eta)   dstats$mu.eta(eta),

    # ----------------------------
    # Derivadas numéricas (robustas a longitudes)
    # ----------------------------
    dldm = function(y, mu, sigma, ...) {
      bd <- .get_bd_sb(...)
      rr <- .recycle_sb(y, mu, sigma, bd)
      y <- rr$y; mu <- rr$mu; sigma <- rr$sigma; bd <- rr$bd

      vapply(seq_along(y), function(i){
        yi <- y[i]; bd_i <- bd[i]
        mu_i <- .clip_mu(mu[i])
        sg_i <- .clip_sg(sigma[i])

        numDeriv::grad(function(m){
          m <- .clip_mu(m)
          as.numeric(dSB_fun(yi, mu = m, sigma = sg_i, bd = bd_i, log = TRUE))
        }, mu_i)
      }, 0.0)
    },

    d2ldm2 = function(y, mu, sigma, ...) {
      bd <- .get_bd_sb(...)
      rr <- .recycle_sb(y, mu, sigma, bd)
      y <- rr$y; mu <- rr$mu; sigma <- rr$sigma; bd <- rr$bd

      vapply(seq_along(y), function(i){
        yi <- y[i]; bd_i <- bd[i]
        mu_i <- .clip_mu(mu[i])
        sg_i <- .clip_sg(sigma[i])

        as.numeric(numDeriv::hessian(function(m){
          m <- .clip_mu(m)
          as.numeric(dSB_fun(yi, mu = m, sigma = sg_i, bd = bd_i, log = TRUE))
        }, mu_i))
      }, 0.0)
    },

    dldd = function(y, mu, sigma, ...) {
      bd <- .get_bd_sb(...)
      rr <- .recycle_sb(y, mu, sigma, bd)
      y <- rr$y; mu <- rr$mu; sigma <- rr$sigma; bd <- rr$bd

      vapply(seq_along(y), function(i){
        yi <- y[i]; bd_i <- bd[i]
        mu_i <- .clip_mu(mu[i])
        sg_i <- .clip_sg(sigma[i])

        numDeriv::grad(function(s){
          s <- .clip_sg(s)
          as.numeric(dSB_fun(yi, mu = mu_i, sigma = s, bd = bd_i, log = TRUE))
        }, sg_i)
      }, 0.0)
    },

    d2ldd2 = function(y, mu, sigma, ...) {
      bd <- .get_bd_sb(...)
      rr <- .recycle_sb(y, mu, sigma, bd)
      y <- rr$y; mu <- rr$mu; sigma <- rr$sigma; bd <- rr$bd

      vapply(seq_along(y), function(i){
        yi <- y[i]; bd_i <- bd[i]
        mu_i <- .clip_mu(mu[i])
        sg_i <- .clip_sg(sigma[i])

        as.numeric(numDeriv::hessian(function(s){
          s <- .clip_sg(s)
          as.numeric(dSB_fun(yi, mu = mu_i, sigma = s, bd = bd_i, log = TRUE))
        }, sg_i))
      }, 0.0)
    },

    d2ldmdd = function(y, mu, sigma, ...) {
      bd <- .get_bd_sb(...)
      rr <- .recycle_sb(y, mu, sigma, bd)
      y <- rr$y; mu <- rr$mu; sigma <- rr$sigma; bd <- rr$bd

      vapply(seq_along(y), function(i){
        yi <- y[i]; bd_i <- bd[i]
        mu_i <- .clip_mu(mu[i])
        sg_i <- .clip_sg(sigma[i])

        H <- numDeriv::hessian(function(v){
          m <- .clip_mu(v[1])
          s <- .clip_sg(v[2])
          as.numeric(dSB_fun(yi, mu = m, sigma = s, bd = bd_i, log = TRUE))
        }, c(mu_i, sg_i))

        as.numeric(H[1, 2])
      }, 0.0)
    },

    # Devianza incremental
    G.dev.incr = function(y, mu, sigma, ...) {
      bd <- .get_bd_sb(...)
      rr <- .recycle_sb(y, mu, sigma, bd)
      y <- rr$y; mu <- rr$mu; sigma <- rr$sigma; bd <- rr$bd

      -2 * vapply(seq_along(y), function(i){
        as.numeric(dSB_fun(y[i], mu = .clip_mu(mu[i]), sigma = .clip_sg(sigma[i]),
                           bd = bd[i], log = TRUE))
      }, 0.0)
    },

    # Randomized quantile residuals
    rqres = expression({
      u1 <- pSB_fun(y,     mu = mu, sigma = sigma, bd = bd, lower.tail = TRUE, log.p = FALSE)
      u0 <- pSB_fun(y - 1, mu = mu, sigma = sigma, bd = bd, lower.tail = TRUE, log.p = FALSE)

      u  <- u0 + runif(length(y)) * (u1 - u0)
      u  <- pmin(pmax(u, 1e-12), 1 - 1e-12)
      qnorm(u)
    }),

    # Iniciales
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

  class(fam) <- c("gamlss.family", "family")
  fam
}
