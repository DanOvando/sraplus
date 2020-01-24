#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericVector testpopmodel(double r, double k, double m,double b0,double plim,int years, double sigma_proc,NumericVector catches){
  
  double growth_mult;
  
  NumericVector b_t(years);
  
  b_t(0) = b0;
  
  
  for ( int t = 1; t < years; t++) {
    
    b_t(t) = NA_REAL;
    
    if ((b_t(t - 1) / k) < plim){
      
      growth_mult = b_t(t - 1) / (plim * k);
      
    } else {
      growth_mult = 1;
    }
    
    b_t(t) = (b_t(t - 1) +  growth_mult * ((r  / (m - 1)) * b_t(t - 1) * (1 - pow(b_t(t - 1) / k,m - 1)))
                - catches(t - 1)) *  exp(R::rnorm(-pow(sigma_proc,2)/2, sigma_proc));
    
    
    // if ((b_t(t,i) / ks(i)) > 1.2){
    //   
    //   b_t(t,i) = 1.2 * ks(i);
    //   
    // }
    
    std::cout << b_t(t) << std::endl;
    
    if (b_t(t) <= 0){
      
      break;
      
    }
  } // close population model
  
  return(b_t);
  
}
