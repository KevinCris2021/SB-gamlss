## --- Densidad marginal SB por integración (wrapper a C++) ---
dsb_one <- function(y, n, mu, phi) {
  integrate(
    cpp_joint_sb_,
    lower = 0, upper = 1,
    y = y, n = n, mu = mu, phi = phi,
    subdivisions = 2000L,
    rel.tol = 1e-6,
    abs.tol = 0
  )$value
}

## Vectorizada para uso interno
dsb_vec <- Vectorize(dsb_one)
