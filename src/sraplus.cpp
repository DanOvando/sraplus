#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List sraplus(NumericVector catches,
             NumericVector rs,
             NumericVector ms,
             NumericVector init_deps,
             NumericVector ks,
             NumericVector qs,
             NumericVector sigma_procs,
             NumericVector drawdex,
             NumericVector index_t,
             NumericVector sigma_obs,
             NumericVector log_final_u,
             NumericVector log_final_u_cv,
IntegerVector index_years,
IntegerVector u_years,
int draws,
int n_keep,
int ref_type,
int fit_index,
int use_final_u,
double log_final_ref,
double sigma_dep,
double plim,
int u_prior,
NumericVector u_priors,
double sigma_u
) {

  int years = catches.size();

  NumericMatrix b_t(years, draws);

  NumericMatrix index_hat_t(years, draws);

  NumericMatrix dep_t(years, draws);

  NumericMatrix b_bmsy_t(years, draws);

  NumericMatrix u_umsy_t(years, draws);

  NumericMatrix c_msy_t(years, draws);

  NumericMatrix proc_error_t(years, draws);

  NumericVector crashed(draws);

  NumericVector log_like(draws);

  NumericVector umsy(draws);

  NumericVector bmsy(draws);

  NumericVector msy(draws);

  for (int i = 0; i < draws; i++) {

  double final_ref;

  double growth_mult = 1;

  crashed(i) = 0;

  umsy(i) = (rs(i) / (ms(i) - 1)) * (1 - 1/ms(i));

  bmsy(i) = ks(i)*pow(ms(i), -1/(ms(i)- 1));

  msy(i) = umsy(i) * bmsy(i);

  b_t(0,i) = ks(i) * init_deps(i);

  u_umsy_t(0,i) = (catches(0) / b_t(0,i)) / umsy(i);

  proc_error_t(0,i) = exp(R::rnorm(0, sigma_procs(i)) - pow( sigma_procs(i),2)/2);

  double b_to_k = pow(ms(i), (-1 / (ms(i) - 1)));

  if (plim > b_to_k){

    plim = b_to_k;

  }

  for ( int t = 1; t < years; t++) {

    proc_error_t(t,i) = exp(R::rnorm(0, sigma_procs(i)));


    if ((b_t(t - 1,i) / ks(i)) < plim){

      growth_mult = b_t(t - 1,i) / (plim * ks(i));

    } else {
      growth_mult = 1;
    }

    b_t(t,i) = (b_t(t - 1,i) +  growth_mult * ((rs(i)  / (ms(i) - 1)) * b_t(t - 1,i) * (1 - pow(b_t(t - 1,i) / ks(i),ms(i) - 1)))
 - catches(t - 1)) * proc_error_t(t,i);

 u_umsy_t(t,i) = (catches(t) / b_t(t,i)) / umsy(i);

    if ((b_t(t,i) / ks(i)) > 1.2){

      b_t(t,i) = 1.2 * ks(i);

    }

    if (b_t(t,i) < 0 || umsy(i) >= 1){

        log_like(i) = -1e9;

        crashed(i) = 1;

        break;

    }

  } // close pop thing

  index_hat_t(_,i) = b_t(_,i) * qs(i);

  dep_t(_,i) = b_t(_,i) / ks(i);

  b_bmsy_t(_,i) = b_t(_,i) / bmsy(i);

  c_msy_t(_,i) = catches / msy(i);

  //// assign penalties ////

  if (crashed(i) == 0) {

  if (ref_type == 0){

    final_ref = dep_t(years - 1,i);
  }

  if (ref_type == 1){

    final_ref = b_bmsy_t(years - 1,i);
  }

  log_like(i) += R::dnorm(log(final_ref),log_final_ref, sigma_dep, true);

  if (u_prior == 1){

    for (int t = 0; t < u_years.size(); t++) {

      log_like(i) += R::dnorm(log(u_priors(t) + 1e-6), log(u_umsy_t(u_years(t) - 1,i) + 1e-6), sigma_u, true);

    }

  }

  if (fit_index == 1) {

    for (int t = 0; t < index_years.size(); t++) {

    log_like(i) += R::dnorm(log(index_t(t)), log(index_hat_t(index_years(t) - 1,i)), sigma_obs(i), true);

    }

  }

  if (use_final_u == 1){

    for (int j = 0; j < log_final_u.size(); j++){

      log_like(i) += R::dnorm(log(u_umsy_t(years - 1, i) + 1e-6), log_final_u(j),log_final_u_cv(j), true);

    } // cluse use final u loops

  } // close use final u

  } // close if crashed is false

  log_like(i) = exp(log_like(i)) * (1 - crashed(i));

} // closd SIR loop


//// run SIR ////

NumericVector scaled_like = (log_like) / (sum(log_like) + 1e-6);

NumericVector keepers = sample(drawdex, n_keep, 1, scaled_like) + 1;


    return Rcpp::List::create(
      Rcpp::Named("b_t") = b_t,
      Rcpp::Named("crashed") = crashed,
      Rcpp::Named("b_bmsy_t") = b_bmsy_t,
      Rcpp::Named("u_umsy_t") = u_umsy_t,
      Rcpp::Named("dep_t") = dep_t,
      Rcpp::Named("keepers") = keepers,
      Rcpp::Named("index_hat_t") = index_hat_t,
      Rcpp::Named("c_msy_t") = c_msy_t,
      Rcpp::Named("r") = rs,
      Rcpp::Named("m") = ms,
      Rcpp::Named("k") = ks);

} // close popmodel
