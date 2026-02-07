# R/sb_family.R

# bd de gamlss
.get_bd_sb <- function(...) {
  dots <- list(...)

  if (!is.null(dots$bd)) return(dots$bd)
  if (!is.null(dots$m))  return(dots$m)

  # fallback (tu caso actual suele funcionar así)
  if (exists("bd", inherits = TRUE)) return(get("bd", inherits = TRUE))
  if (exists("m",  inherits = TRUE)) return(get("m",  inherits = TRUE))

  stop("No se encontró 'bd' (o 'm'). Pasa gamlss(..., bd = <vector>) o define bd/m en el entorno.")
}


SB <- function(mu.link = "logit", sigma.link = "log") {

  if (!exists("dSB")) stop("Define primero dSB() (wrapper).")
  if (!exists("pSB")) stop("Define primero pSB().")

  mstats <- gamlss.dist::checklink("mu.link","SB",substitute(mu.link),
                                   c("logit","probit","cloglog","log","own"))
  dstats <- gamlss.dist::checklink("sigma.link","SB",substitute(sigma.link),
                                   c("log","identity","inverse"))

  mu.link    <- as.character(substitute(mu.link))
  sigma.link <- as.character(substitute(sigma.link))

  # --- Constantes numéricas coherentes con dSB() ---
  eps_mu   <- 1e-12
  eps_sig  <- 1e-12
  sig_cap  <- 50

  # --- Pasos fijos para numDeriv (mejor estabilidad) ---
  eps_grad_mu  <- 1e-5
  eps_hess_mu  <- 1e-4
  eps_grad_sig <- 1e-5
  eps_hess_sig <- 1e-4
  eps_hess_mix <- 1e-4

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

    ## --- Derivadas numéricas consistentes con dSB() ---

    dldm = function(y, mu, sigma, ...) {
      bd <- .get_bd_sb(...)
      vapply(seq_along(y), function(i){

        mu_i <- pmin(pmax(mu[i], eps_mu), 1 - eps_mu)
        sg_i <- pmin(pmax(sigma[i], eps_sig), sig_cap)   
        bd_i <- bd[i]
        yi   <- y[i]

        numDeriv::grad(
          function(m){
            m <- pmin(pmax(m, eps_mu), 1 - eps_mu)
            as.numeric(dSB(yi, mu = m, sigma = sg_i, bd = bd_i, log = TRUE))
          },
          mu_i,
          method.args = list(eps = eps_grad_mu)         
        )
      }, 0.0)
    },

    d2ldm2 = function(y, mu, sigma, ...) {
      bd <- .get_bd_sb(...)
      vapply(seq_along(y), function(i){

        mu_i <- pmin(pmax(mu[i], eps_mu), 1 - eps_mu)
        sg_i <- pmin(pmax(sigma[i], eps_sig), sig_cap)   
        bd_i <- bd[i]
        yi   <- y[i]

        as.numeric(numDeriv::hessian(
          function(m){
            m <- pmin(pmax(m, eps_mu), 1 - eps_mu)
            as.numeric(dSB(yi, mu = m, sigma = sg_i, bd = bd_i, log = TRUE))
          },
          mu_i,
          method.args = list(eps = eps_hess_mu)          
        ))
      }, 0.0)
    },

    dldd = function(y, mu, sigma, ...) {
      bd <- .get_bd_sb(...)
      vapply(seq_along(y), function(i){

        mu_i <- pmin(pmax(mu[i], eps_mu), 1 - eps_mu)
        sg_i <- pmin(pmax(sigma[i], eps_sig), sig_cap)  
        bd_i <- bd[i]
        yi   <- y[i]

        numDeriv::grad(
          function(s){
            s <- pmin(pmax(s, eps_sig), sig_cap)         # <<< CAMBIO 1 (cap)
            as.numeric(dSB(yi, mu = mu_i, sigma = s, bd = bd_i, log = TRUE))
          },
          sg_i,
          method.args = list(eps = eps_grad_sig)         
        )
      }, 0.0)
    },

    d2ldd2 = function(y, mu, sigma, ...) {
      bd <- .get_bd_sb(...)
      vapply(seq_along(y), function(i){

        mu_i <- pmin(pmax(mu[i], eps_mu), 1 - eps_mu)
        sg_i <- pmin(pmax(sigma[i], eps_sig), sig_cap)   
        bd_i <- bd[i]
        yi   <- y[i]

        as.numeric(numDeriv::hessian(
          function(s){
            s <- pmin(pmax(s, eps_sig), sig_cap)        
            as.numeric(dSB(yi, mu = mu_i, sigma = s, bd = bd_i, log = TRUE))
          },
          sg_i,
          method.args = list(eps = eps_hess_sig)        
        ))
      }, 0.0)
    },

    d2ldmdd = function(y, mu, sigma, ...) {
      bd <- .get_bd_sb(...)
      vapply(seq_along(y), function(i){

        mu_i <- pmin(pmax(mu[i], eps_mu), 1 - eps_mu)
        sg_i <- pmin(pmax(sigma[i], eps_sig), sig_cap)   
        bd_i <- bd[i]
        yi   <- y[i]

        H <- numDeriv::hessian(
          function(v){
            m <- pmin(pmax(v[1], eps_mu), 1 - eps_mu)
            s <- pmin(pmax(v[2], eps_sig), sig_cap)      
            as.numeric(dSB(yi, mu = m, sigma = s, bd = bd_i, log = TRUE))
          },
          c(mu_i, sg_i),
          method.args = list(eps = eps_hess_mix)         
        )

        as.numeric(H[1, 2])
      }, 0.0)
    },


    ## Devianza:
    G.dev.incr = function(y, mu, sigma, ...) {
      bd <- .get_bd_sb(...)
      -2 * vapply(seq_along(y), function(i){
        m  <- pmin(pmax(mu[i], eps_mu), 1 - eps_mu)
        sg <- pmin(pmax(sigma[i], eps_sig), sig_cap)     # <<< CAMBIO 1 (cap)
        as.numeric(dSB(y[i], mu = m, sigma = sg, bd = bd[i], log = TRUE))
      }, 0.0)
    },

    ## RQR:
    rqres = expression({
      rqres(pfun = "pSB", type = "Discrete",
            ymin = 0, y = y, mu = mu, sigma = sigma, bd = bd)
    }),

    ## Iniciales y dominios
    mu.initial = expression({
      mu <- pmin(pmax((y + 0.5) / (bd + 1), 1e-6), 1 - 1e-6)
    }),

    sigma.initial = expression({ sigma <- rep(0.7, length(y)) }),

    mu.valid    = function(mu)    all(mu > 0 & mu < 1),
    sigma.valid = function(sigma) all(sigma > 0),
    y.valid     = function(y)     all(y >= 0)

  ) -> fam

  class(fam) <- c("gamlss.family","family")
  fam
}
