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
      if (!is.finite(Fmid)) Fmid <- 0  # <<< FIX CLAVE

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
