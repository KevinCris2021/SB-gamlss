#' Wrapper C++: kernel joint SB
#' @export
cpp_joint_sb <- function(w, y, n, mu, phi) {
  cpp_joint_sb_wrap(w, y, n, mu, phi)
}
