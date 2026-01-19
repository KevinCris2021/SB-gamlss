# Helpers internos para el family SB (NO exportar)

.get_bd_sb <- function(...) {
  dots <- list(...)
  if (!is.null(dots$bd)) return(dots$bd)
  if (!is.null(dots$m))  return(dots$m)
  stop("No se encontró 'bd' (o 'm'). Pasa gamlss(..., bd = <vector>)")
}

.recycle_sb <- function(y, mu, sigma, bd) {
  n <- length(y)

  # bd puede venir NULL o escalar
  if (is.null(bd)) bd <- rep(max(y, na.rm = TRUE), n)

  list(
    y     = y,
    mu    = rep(mu,    length.out = n),
    sigma = rep(sigma, length.out = n),
    bd    = rep(bd,    length.out = n)
  )
}

.clip_mu <- function(mu) pmin(pmax(mu, 1e-12), 1 - 1e-12)
.clip_sg <- function(sg) pmax(sg, 1e-12)

