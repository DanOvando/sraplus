#include <Rcpp.h>
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

NumericVector popmodel(double r, double k, double m,double b0,double plim,int years, double sigma_proc,NumericVector catches){
  
  double growth_mult;
  
  NumericVector b_t(years,0.0);
  
  // std::vector<double> b_t(years, 0.0);
  
  b_t(0) = b0;
  
  for ( int t = 1; t < years; t++) {
    
    if ((b_t(t - 1) / k) < plim){
      
      growth_mult = b_t(t - 1) / (plim * k);
      
    } else {
      growth_mult = 1;
    }
    
    b_t(t) = (b_t(t - 1) +  growth_mult * ((r  / (m - 1)) * b_t(t - 1) * (1 - pow(b_t(t - 1) / k,m - 1)))
                - catches(t - 1)) *  exp(R::rnorm(-pow(sigma_proc,2)/2, sigma_proc));
    
    b_t(t) = std::max(0.0,b_t(t));
    
    // if ((b_t(t,i) / ks(i)) > 1.2){
    //   
    //   b_t(t,i) = 1.2 * ks(i);
    //   
    // }
    
    // std::cout << b_t(t) << std::endl;
    
    if (b_t(t) <= 0){
      
      break;
      
    }
  } // close population model
  
  return(b_t);
  
}


// [[Rcpp::export]]

NumericVector graddesc(double r, 
                       double m,
                       double init_deps,
                       double plim,
                       int years, 
                       double sigma_procs,
                       NumericVector catches, 
                       double final_state,
                       double learn_rate){
  
double growth_mult;
  
double delta;

double last_proposal;

double prop_error;

double new_proposal;

double grad_dep;

double conv_error;

double log_like;

double k;

int counter;

int crashed;

counter = 0;

delta = 100;

NumericVector proposal_result(years);

last_proposal = log(1000000 * max(catches));

// proposal_result = popmodel(r,m,exp(last_proposal), init_dep, catch_t, proc_errors, time, eps,plim);


proposal_result = popmodel(r, exp(last_proposal), m, exp(last_proposal) * init_deps, plim, years, sigma_procs, catches);


// std::cout <<  "counter is" << proposal_result(years - 1)  << std::endl;

prop_error =  log(proposal_result(years - 1) / exp(last_proposal)) - log(final_state);



while(delta > 1e-3){
  
  new_proposal = last_proposal -  learn_rate * prop_error;
  
  // std::cout << "new proposal is" << new_proposal << std::endl;
  
  // proposal_result = popmodel(r,m,exp(new_proposal), init_dep, catch_t, proc_errors, time, eps,plim);
  
  
  
  proposal_result = popmodel(r, exp(new_proposal), m, exp(new_proposal) * init_deps, plim, years, sigma_procs, catches);
  
  // std::cout <<  "initial proposal is " << proposal_result << std::endl;
  
  
  for (int t = 0; t< years; t++){
    
    // std::cout << proposal_result(t) << std::endl;
    
    if (proposal_result(t) <= 0){
      
      crashed = 1;
      
      log_like = -1e9;
      
      break;
    }
    
  }
  
  
  
  grad_dep = proposal_result(years - 1) / exp(new_proposal);
  
  
  prop_error =  log(grad_dep + 1e-6) - log(final_state);
  
  // std::cout << "prop dep is" << proposal_result[time - 1] / exp(new_proposal) << std::endl;
  
  last_proposal = new_proposal;
  
  delta = sqrt(pow(grad_dep - final_state,2));
  
  // std::cout <<  "new proposal is  is" << new_proposal  << std::endl;
  
  conv_error = delta;
  
  // std::cout<< delta << std::endl;
  
  counter += 1;
  
  if (counter > 1e6){
    
    delta = -999;
    
    crashed = 1;
  } // close escape hatch
  
  std::cout <<  "delta is" << delta  << std::endl;
  
  
} // close while looop for gradient descnet

// std::cout <<  "counter is" << counter << std::endl;



// std::cout <<  "delta is" << delta << std::endl;

k = exp(new_proposal);



// log_like(i) += log(conv_error);  

return(proposal_result);
}
