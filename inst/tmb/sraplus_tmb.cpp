#include <TMB.hpp>

template<class Type>
Type growfoo(Type r, Type m, Type b, Type plim)
{
  // if ( b > plim ){

    Type growth = (r  / (m - 1)) * b * (1 - pow(b,m - 1));

  //   return growth;
  // } else {

    // Type growth = (b * 1/plim) * ((r  / (m - 1)) * b * (1 - pow(b,m - 1)));

    return growth;
  //}

} // close function

template<class Type>
Type posfun(Type x, Type eps, Type &pen)
{
  if ( x >= eps ){
    return x;
  } else {
    pen += Type(0.01) * pow(x-eps,2);
    return eps/(Type(2.0)-x/eps);
  }
}


template<class Type>
Type objective_function<Type>::operator() ()
{

  //// data and all that ////
  // hello there

  DATA_VECTOR(catch_t);

  DATA_VECTOR(index_t);

  DATA_INTEGER(time); // number of time steps

  DATA_INTEGER(fit_index); //

  DATA_INTEGER(use_u_prior); //

  DATA_INTEGER(ref_type); //0 means that initial and final are relative to K, 1 to Bmsy

  // DATA_INTEGER(use_init); // use initial reference point

  DATA_INTEGER(use_final); // use final reference point

  DATA_INTEGER(use_final_u); // use final U/Umsy

  DATA_VECTOR(log_final_u);

  DATA_VECTOR(log_final_u_cv);

  DATA_VECTOR(u_priors);

  DATA_IVECTOR(u_years);

  // DATA_SCALAR(u_cv);

  DATA_IVECTOR(index_years);

  DATA_SCALAR(log_k_guess);

  DATA_SCALAR(plim);

  DATA_SCALAR(f_cv);

  DATA_SCALAR(log_r_prior);

  // DATA_SCALAR(sigma_proc_prior);

  // DATA_SCALAR(sigma_proc_prior_cv);

  DATA_SCALAR(log_r_cv);

  DATA_SCALAR(sigma_u);

  DATA_SCALAR(log_init_dep_prior);

  DATA_SCALAR(log_init_dep_cv);

  DATA_SCALAR(log_final_dep_prior);

  DATA_SCALAR(log_final_dep_cv);

  //// parameters ////

  PARAMETER(log_init_dep);

  PARAMETER(log_r);

  PARAMETER(log_k);

  PARAMETER(log_q);

  PARAMETER(log_sigma_proc);

  PARAMETER_VECTOR(uc_proc_errors);

  PARAMETER(log_sigma_obs);

  PARAMETER(log_m);

  PARAMETER_VECTOR(inv_f_t);

  //// model ////

  Type crashed = 0;

  Type nll = 0.0;

  Type penalty = 0.0;

  vector<Type> b_t(time);

  vector<Type> index_hat_t(time);

  vector<Type> dep_t(time);

  vector<Type> growth_t(time);

  vector<Type> b_v_bmsy(time);

  vector<Type> u_v_umsy(time);

  vector<Type> u_t(time);

  Type sigma_obs = exp(log_sigma_obs);

  Type sigma_proc = exp(log_sigma_proc);

  // Type sigma_u = exp(log_sigma_u);


  vector<Type> catch_hat_t(time - 1);

  vector<Type> proc_errors(time - 1);

  proc_errors = exp(sigma_proc * uc_proc_errors - pow(sigma_proc,2)/2);

  Type k = exp(log_k);

  Type r = exp(log_r);

  Type q = exp(log_q);

  Type m = exp(log_m);

  Type b_to_k = pow(m, (-1 / (m - 1)));

  if (plim > b_to_k){

    plim = b_to_k;

  }

  Type init_dep = exp(log_init_dep);

  Type bmsy = k*pow(m, -1/(m - 1));

  Type umsy = (r / (m - 1)) * (1 - 1/m);

  Type msy = bmsy * umsy;

  // b_t(0) = k * init_dep;

  Type fpen = 0.;
  
  b_t(0) = init_dep;

  // index_hat_t(0) = b_t(0) * q;

  for (int t = 1; t<time; t++){

    u_t(t - 1) = 1 / (1 + exp(-inv_f_t(t - 1)));

    catch_hat_t(t - 1) = u_t(t - 1) * b_t(t - 1) * k;


    growth_t(t - 1) = growfoo(r,m,b_t(t - 1),plim);

    b_t(t) = (b_t(t - 1) +  growth_t(t - 1) - catch_hat_t(t - 1) / k) * proc_errors(t - 1);

    b_t(t) = posfun(b_t(t) ,Type(.001),fpen);
    

  } // close population model


  u_t(time - 1) = catch_t(time -1)  / (b_t(time - 1) * k);

  // 1 / (1 + exp(-inv_f_t(time - 1)));

  // catch_hat_t(time - 1) = u_t(time -1) * b_t(time - 1) * k;

  growth_t(time - 1) = (r  / (m - 1)) * b_t(time - 1) * (1 - pow(b_t(time - 1),m - 1));

  b_t = b_t * k;

  b_v_bmsy = b_t / bmsy;

  u_v_umsy = u_t / umsy;

  dep_t = b_t / k;

  index_hat_t = q * b_t;
  
  nll += fpen;
  
  // umsy penalty

  if (umsy > 1){

    nll += pow(umsy - 1, 2);

  }


  if (fit_index == 1) {

    for (int t = 0; t < index_years.size(); t++){

      nll -= dnorm(log(index_t(t)), log(index_hat_t(index_years(t) - 1)), Type(0.00001) + sigma_obs, true);

    }

  }

  if (use_u_prior == 1){


    for (int t = 0; t < u_years.size(); t++){

      nll -= dnorm(log(u_v_umsy(u_years(t) - 1) + 1e-3), log(u_priors(t) + 1e-3), sigma_u, true);

    }

  }

  for (int t = 0; t < (time - 1); t++){

    nll -= dnorm(log(catch_t(t)), log(catch_hat_t(t)), Type(.05), true);

    nll -= dnorm(uc_proc_errors(t), Type(0), Type(1), true);

  }

  for (int t = 1; t < (time - 1); t++){


    nll -= dnorm(inv_f_t(t), inv_f_t(t - 1), f_cv, true);

  }

  // nll -= dbinom(crashed, Type(time), Type(0.01), true);

  // nll += crashed;

  nll -= dnorm(log_r, log_r_prior, log_r_cv, true);

  Type final_ref = 0.5;

  Type init_ref = 1.0;


  if (ref_type == 0){

    final_ref = dep_t(time - 1);

    init_ref = dep_t(0);

  }

  if (ref_type == 1) {

    // std::cout << "hello?" << "//n";

    final_ref = b_v_bmsy(time - 1);

    init_ref = b_v_bmsy(0);

  }


  if (use_final == 1) {

    nll -= dnorm(log(final_ref), log_final_dep_prior,log_final_dep_cv, true);

    // nll -= dnorm((final_ref), exp(log_final_dep_prior),log_final_dep_cv, true);

  }

  // allow for arbitarry numbers of priors on final U/Umsy, so that you can use FMI + SAR

  if (use_final_u == 1){

    for (int i = 0; i < log_final_u.size(); i++){

      nll -= dnorm(log(u_v_umsy(time - 1)), log_final_u(i),log_final_u_cv(i), true);

    } // cluse use final u loops

  } // close use final u

  nll -= dnorm(log(init_ref), log_init_dep_prior, log_init_dep_cv, true);

  nll -= dnorm(log_sigma_obs,Type(-3),Type(0.2), true);

  nll -= dnorm(log_sigma_proc,Type(-3), Type(0.1), true);

  nll -= dnorm(log_k,log_k_guess,Type(4), true);

  // nll -= dnorm(log_q,Type(-7),Type(2), true);

  nll -= dnorm(log_m,Type(0.1),Type(.1), true);


  vector<Type> log_bt = log(b_t);

  vector<Type> log_b = log(b_v_bmsy);

  vector<Type> log_dep = log(dep_t);

  vector<Type> log_u = log(u_v_umsy);

  // recompile please
  vector<Type> log_c_div_msy = log(catch_t / msy);

  vector<Type> log_chat = log(catch_hat_t);

  vector<Type> log_ihat = log(index_hat_t);

  vector<Type> ck = catch_t / k;

  REPORT(log_chat)

    REPORT(umsy);

  REPORT(k);

  REPORT(ck);

  REPORT(log_ihat);

  ADREPORT(log_ihat);

  ADREPORT(log_b);

  ADREPORT(log_bt);

  ADREPORT(log_dep);

  ADREPORT(log_u);

  ADREPORT(log_c_div_msy);
  
  ADREPORT(log_chat);

  // one more time for true report

  REPORT(log_b);

  REPORT(log_bt);

  REPORT(log_dep);

  REPORT(log_u);

  REPORT(log_c_div_msy);

  REPORT(log_chat);

  REPORT(final_ref);

  REPORT(proc_errors);

  ADREPORT(proc_errors);

  REPORT(k);

  REPORT(plim);

  REPORT(dep_t);

  REPORT(crashed);

  ADREPORT(plim);

  ADREPORT(b_to_k);

  ADREPORT(crashed);

  REPORT(b_t);

  REPORT(growth_t);

  REPORT(index_hat_t);

  ADREPORT(index_hat_t);

  ADREPORT(r);

  ADREPORT(m);

  REPORT(msy);

  ADREPORT(umsy);

  REPORT(bmsy);

  ADREPORT(k);

  REPORT(q);

  REPORT(r);

  REPORT(m);

  ADREPORT(penalty);


  return nll;
}
