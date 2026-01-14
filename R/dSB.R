dSB<- function(x, mu = 0.5, sigma = 1, bd = 10, log = FALSE,
                     rel.tol = 1e-8, abs.tol = 0, subdivisions = 200L) {
  
  if (any(!is.finite(mu)) || any(!is.finite(sigma)) || any(!is.finite(bd)))
    stop("mu/sigma/bd must be finite")
  if (any(mu < 0 | mu > 1))  stop("mu must be between 0 and 1")
  if (any(sigma <= 0))       stop("sigma must be greater than 0")


  eps <- 1e-10
  mu    <- pmin(pmax(mu,    eps), 1 - eps)
  sigma <- pmin(pmax(sigma, eps), 50)   

  ly <- max(length(x), length(mu), length(sigma), length(bd))
  x     <- rep(x,     length.out = ly)
  mu    <- rep(mu,    length.out = ly)
  sigma <- rep(sigma, length.out = ly)
  bd    <- rep(bd,    length.out = ly)

  tiny <- .Machine$double.xmin

  safe_int <- function(y_i, n_i, mu_i, phi_i) {
    # fuera de soporte o no entero -> prob 0 (log tiny)
    if (!is.finite(y_i) || !is.finite(n_i) ||
        y_i < 0 || y_i > n_i || y_i != floor(y_i)) {
      return(if (log) log(tiny) else 0)
    }
    # caso casi binomial
    if (phi_i < 1e-6) {
      val <- dbinom(y_i, n_i, mu_i)
      val <- max(val, tiny)
      return(if (log) log(val) else val)
    }
    # integración
    val <- try(integrate(
      function(w) {
        as.numeric(cpp_joint_sb_(w, y_i, n_i, mu_i, phi_i))
      },
      lower = 1e-8, upper = 1 - 1e-8,
      rel.tol = rel.tol, abs.tol = abs.tol,
      subdivisions = subdivisions
    ), silent = TRUE)

    if (inherits(val, "try-error") || !is.finite(val$value) || val$value <= 0) {
      return(if (log) log(tiny) else tiny)
    }
    if (log) log(val$value) else val$value
  }

  out <- vapply(seq_len(ly),
                function(i) safe_int(x[i], bd[i], mu[i], sigma[i]),
                numeric(1))
  out
}
