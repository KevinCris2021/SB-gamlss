// Simplex-Binomial: joint kernel for integration (Rcpp export)
#include <Rcpp.h>
#include <algorithm>
#include <math.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//



#define GETV(x, i)      x[i % x.length()]    // wrapped indexing of vector

double pi_value = M_PI;

double cpp_joint_sb(double w, int y, int n, double mu, double phi){
  return exp(R::dbinom(y,n,w,TRUE)-
             0.5*pow(w-mu,2)/(w*(1-w)*pow(mu*(1-mu),2))/phi-
             0.5*(3*log(w)+3*log(1-w)+log(phi)+log(2*pi_value))
             );
}


// [[Rcpp::export(rng = false)]]
NumericVector cpp_joint_sb_(
    const NumericVector& w,
    const NumericVector& y,
    const NumericVector& n,
    const NumericVector& mu,
    const NumericVector& phi
) {
  
  if (std::min({w.length(), y.length(), n.length(),
               mu.length(), phi.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    w.length(),
    y.length(),
    n.length(),
    mu.length(),
    phi.length()
  });
  NumericVector r(Nmax);
  
  for (int i = 0; i < Nmax; i++){
    r[i] = cpp_joint_sb(GETV(w, i), GETV(y, i), GETV(n, i),
                        GETV(mu, i), GETV(phi, i));
  }
  
  return r;
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
