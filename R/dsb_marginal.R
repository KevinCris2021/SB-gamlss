## --- Densidad marginal SB por integración (wrapper a C++) ---

# wrapper seguro: integrate() puede pasar w vector
cpp_joint_sb_wrap <- function(w, y, n, mu, phi) {
  vapply(w, function(ww) cpp_joint_sb_(ww, y = y, n = n, mu = mu, phi = phi), 0.0)
}

dsb_one <- function(y, n, mu, phi) {
  integrate(
    cpp_joint_sb_wrap,
    lower = 0, upper = 1,
    y = y, n = n, mu = mu, phi = phi,
    subdivisions = 2000L,
    rel.tol = 1e-6,
    abs.tol = 0
  )$value
}

## Vectorizada para uso interno
dsb_vec <- Vectorize(dsb_one)

