#include <Rcpp.h>
#include <R_ext/Rdynload.h>
using namespace Rcpp;

// Declarations
Rcpp::NumericVector cpp_joint_sb_(const Rcpp::NumericVector& w,
                                  const Rcpp::NumericVector& y,
                                  const Rcpp::NumericVector& n,
                                  const Rcpp::NumericVector& mu,
                                  const Rcpp::NumericVector& phi);

// Wrapper for .Call
extern "C" SEXP _SBgamlss_cpp_joint_sb_(SEXP wSEXP, SEXP ySEXP, SEXP nSEXP, SEXP muSEXP, SEXP phiSEXP) {
  BEGIN_RCPP
  Rcpp::NumericVector w(wSEXP);
  Rcpp::NumericVector y(ySEXP);
  Rcpp::NumericVector n(nSEXP);
  Rcpp::NumericVector mu(muSEXP);
  Rcpp::NumericVector phi(phiSEXP);

  Rcpp::NumericVector res = cpp_joint_sb_(w, y, n, mu, phi);
  return res;
  END_RCPP
}

// Registration
static const R_CallMethodDef CallEntries[] = {
  {"_SBgamlss_cpp_joint_sb_", (DL_FUNC) &_SBgamlss_cpp_joint_sb_, 5},
  {NULL, NULL, 0}
};

extern "C" void R_init_SBgamlss(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
