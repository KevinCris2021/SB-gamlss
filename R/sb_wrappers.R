dSB<- function(x, mu = 0.5, sigma = 1, bd = 10, log = FALSE,
                     rel.tol = 1e-8, abs.tol = 0, subdivisions = 200L) {
  # --- Validaciones ---
  if (any(!is.finite(mu)) || any(!is.finite(sigma)) || any(!is.finite(bd)))
    stop("mu/sigma/bd must be finite")
  if (any(mu < 0 | mu > 1))  stop("mu must be between 0 and 1")
  if (any(sigma <= 0))       stop("sigma must be greater than 0")

  # --- Clamps suaves ---
  eps <- 1e-10
  mu    <- pmin(pmax(mu,    eps), 1 - eps)
  sigma <- pmin(pmax(sigma, eps), 50)   # cap sigma para evitar gradientes explosivos

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

pSB <- function(q, mu = 0.5, sigma = 1, bd = 10, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu <= 0) | any(mu >= 1)) stop("mu must be between 0 and 1")
  if (any(sigma <= 0))             stop("sigma must be greater than 0")

  # --- Igualar longitudes ANTES de cualquier chequeo cruzado ---
  ly    <- max(length(q), length(mu), length(sigma), length(bd))
  q     <- rep(q,     length.out = ly)
  mu    <- rep(mu,    length.out = ly)
  sigma <- rep(sigma, length.out = ly)
  bd    <- rep(bd,    length.out = ly)

  tiny_phi <- 1e-6

  cdf <- mapply(function(qi, mui, si, bdi) {
    # Soporte y atajos
    if (!is.finite(qi) || !is.finite(bdi)) return(NA_real_)
    if (qi < 0)  return(0)        # P(Y <= q) = 0 si q < 0
    if (qi >= bdi) return(1)      # P(Y <= q) = 1 si q >= bd

    # Aproximación binomial cuando la dispersión es pequeña o tiende a cero
    if (si < tiny_phi) return(pbinom(qi, size = bdi, prob = mui))

    # Acumular pmf discreta 0:qi desde dSB()
    pmf <- dSB(0:qi, mu = mui, sigma = si, bd = bdi)
    pmf[!is.finite(pmf) | pmf < 0] <- 0
    out <- sum(pmf)
    out <- min(max(out, 0), 1)
    out
  }, q, mu, sigma, bd)

  if (!lower.tail) cdf <- 1 - cdf
  if (log.p)       cdf <- log(cdf)
  cdf
}

   qSB <- function(p, mu = 0.5, sigma = 1, bd = 10,
                lower.tail = TRUE, log.p = FALSE, fast = TRUE)
{
  if (any(mu <= 0) | any(mu >= 1)) stop("mu must be between 0 and 1")
  if (any(sigma <= 0))             stop("sigma must be greater than 0")

  ly    <- max(length(p), length(mu), length(sigma), length(bd))
  p     <- rep(p,     length.out = ly)
  mu    <- rep(mu,    length.out = ly)
  sigma <- rep(sigma, length.out = ly)
  bd    <- rep(bd,    length.out = ly)

  if (log.p)       p <- exp(p)
  if (!lower.tail) p <- 1 - p

  eps <- 1e-12
  if (any(p < -eps | p > 1 + eps))
    stop("p must be between 0 and 1 (after applying log.p/lower.tail)")
  p <- pmin(pmax(p, 0), 1)

  tiny_phi <- 1e-6
  ans <- integer(ly)
  tol <- 1e-12

  for (i in seq_len(ly)) {
    if (p[i] <= 0)  { ans[i] <- 0L;       next }
    if (p[i] >= 1)  { ans[i] <- bd[i];    next }

    if (sigma[i] < tiny_phi) {
      ans[i] <- stats::qbinom(p[i], size = bd[i], prob = mu[i],
                              lower.tail = TRUE, log.p = FALSE)
      next
    }

    lo <- 0L
    hi <- as.integer(bd[i])
    while (lo < hi) {
      mid  <- as.integer((lo + hi) %/% 2)
      Fmid <- pSB(mid, mu = mu[i], sigma = sigma[i], bd = bd[i],
                  lower.tail = TRUE, log.p = FALSE)
      if (!is.finite(Fmid)) Fmid <- 0 

      if (Fmid + tol >= p[i]) {
        hi <- mid
      } else {
        lo <- mid + 1L
      }
    }
    ans[i] <- lo
  }

  as.integer(ans)
}



rSB <- function(n, mu = 0.5, sigma = 1, bd = 10) {
  if (any(mu <= 0) | any(mu >= 1)) stop("mu must be between 0 and 1")
  if (any(sigma <= 0))             stop("sigma must be greater than 0")
  if (any(n <= 0))                 stop("n must be a positive integer")

  n  <- as.integer(ceiling(n))
  ly <- max(length(mu), length(sigma), length(bd))
  mu    <- rep(mu,    length.out = ly)
  sigma <- rep(sigma, length.out = ly)
  bd    <- rep(bd,    length.out = ly)

  tiny_phi <- 1e-6
  out <- integer(n)

  key <- paste(mu, sigma, bd, sep = "|")
  uniq_keys <- unique(key)
  cache_cdf <- vector("list", length(uniq_keys))
  names(cache_cdf) <- uniq_keys

  for (i in seq_len(n)) {
    j <- if (ly == 1L) 1L else ((i - 1L) %% ly) + 1L

    if (sigma[j] < tiny_phi) {
      out[i] <- stats::rbinom(1L, size = bd[j], prob = mu[j])
      next
    }

    kj <- key[j]
    cdf <- cache_cdf[[kj]]
    if (is.null(cdf)) {
      pmf <- dSB(0:bd[j], mu = mu[j], sigma = sigma[j], bd = bd[j])
      pmf[!is.finite(pmf) | pmf < 0] <- 0
      s <- sum(pmf)
      if (!is.finite(s) || s <= 0) {
        out[i] <- stats::rbinom(1L, size = bd[j], prob = mu[j])
        next
      }
      cdf <- cumsum(pmf / s)
      cache_cdf[[kj]] <- cdf
    }

    u <- stats::runif(1L)
    q <- findInterval(u, c(0, cdf), left.open = FALSE) - 1L
    if (q < 0L)       q <- 0L
    if (q > bd[j])    q <- bd[j]
    out[i] <- q
  }
  as.integer(out)
}


             
