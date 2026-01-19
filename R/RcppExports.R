# Minimal Rcpp interface (manual)

cpp_joint_sb_ <- function(w, y, n, mu, phi) {
  .Call("_SBgamlss_cpp_joint_sb_", w, y, n, mu, phi)
}
