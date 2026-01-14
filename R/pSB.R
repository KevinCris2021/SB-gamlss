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

    # Acumular pmf discreta 0:qi desde tu dSB()
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
