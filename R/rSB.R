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
    # ¡ojo! mapeo 0..bd: restar 1L
    q <- findInterval(u, c(0, cdf), left.open = FALSE) - 1L
    if (q < 0L)       q <- 0L
    if (q > bd[j])    q <- bd[j]
    out[i] <- q
  }
  as.integer(out)
}
