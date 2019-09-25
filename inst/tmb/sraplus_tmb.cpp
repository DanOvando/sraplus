#include <TMB.hpp>

template<class Type>
Type growfoo(Type r, Type m, Type b, Type plim)
{

  Type growth = (r  / (m - 1)) * b * (1 - pow(b,m - 1));
  
  return growth;

} // close function

template<class Type>
Type posfun(Type x, Type eps, Type &pen)
{
  pen += CppAD::CondExpLt(x, eps, Type(0.01) * pow(x-eps,2), Type(0));
  return CppAD::CondExpGe(x, eps, x, eps/(Type(2)-x/eps));
}


template<class Type>
Type objective_function<Type>::operator() ()
{
  
  //// data and all that ////

  DATA_VECTOR(catch_t);
  
  DATA_VECTOR(index_t);
  
  DATA_VECTOR(effort_t);
  
  DATA_INTEGER(time); // number of time steps
  
  DATA_INTEGER(fit_index); // fit to an externally supplied index of abundance
  
  DATA_INTEGER(calc_cpue); // fit to an index of abundance generated by a catch and effort series, estimating abundance through baranov equation and q assumption
  
  DATA_INTEGER(use_baranov); // use baranov equation in CPUE standardization? requires estimate of natural mortality
  
  DATA_INTEGER(marginalize_q); // should q be marginalized out? only use for fitting index not CPUE
  
  DATA_INTEGER(use_u_prior); //
  
  DATA_INTEGER(b_ref_type); //0 means that initial and final are relative to K, 1 to Bmsy
  
  DATA_INTEGER(f_ref_type); //0 means f, 1 f/fmsy
  
  // DATA_INTEGER(use_init); // use initial reference point
  
  DATA_INTEGER(use_final_state); // use final reference point
  
  DATA_INTEGER(use_final_u); // use final U/Umsy
  
  DATA_VECTOR(log_final_u);
  
  DATA_VECTOR(log_final_u_cv);
  
  DATA_VECTOR(u_priors);
  
  DATA_IVECTOR(u_years);
  
  DATA_IVECTOR(index_years);
  
  DATA_SCALAR(k_prior);
  
  DATA_SCALAR(k_prior_cv);
  
  DATA_SCALAR(nat_m);
  
  DATA_SCALAR(shape_prior);
  
  DATA_SCALAR(shape_cv);
  
  DATA_SCALAR(plim);
  
  DATA_SCALAR(log_r_prior);
  
  DATA_SCALAR(sigma_proc_prior);
  
  DATA_SCALAR(sigma_proc_prior_cv);
  
  DATA_SCALAR(eps);
  
  DATA_SCALAR(log_r_cv);
  
  DATA_SCALAR(sigma_u);
  
  DATA_SCALAR(log_init_dep_prior);
  
  DATA_SCALAR(log_init_dep_cv);
  
  DATA_SCALAR(log_final_dep_prior);
  
  DATA_SCALAR(log_final_dep_cv);
  
  DATA_SCALAR(q_slope_prior);
  
  DATA_SCALAR(q_slope_cv);
  
  DATA_SCALAR(sigma_obs_prior);
  
  DATA_SCALAR(sigma_obs_prior_cv);
  
  //// parameters ////
  
  PARAMETER(log_init_dep);
  
  PARAMETER(log_r);
  
  PARAMETER(log_k);
  
  PARAMETER(log_q);
  
  PARAMETER(log_q_slope);
  
  PARAMETER(log_sigma_proc);
  
  PARAMETER_VECTOR(log_proc_errors);
  
  PARAMETER(log_sigma_obs);
  
  PARAMETER(log_shape);
  
  //// model ////
  
  Type crashed = 0;
  
  Type nll = 0.0;
  
  Type pen = 0.0;
  
  vector<Type> b_t(time);
  
  vector<Type> q_t(time);
  
  vector<Type> index_hat_t(time);
  
  vector<Type> dep_t(time);
  
  vector<Type> growth_t(time);
  
  vector<Type> b_v_bmsy(time);
  
  vector<Type> u_v_umsy(time);
  
  vector<Type> u_t(time);
  
  vector<Type> f_t(time - 1);
  
  vector<Type> short_u_t(time - 1);
  
  Type q_slope = exp(log_q_slope);
  
  Type sigma_obs = exp(log_sigma_obs);
  
  Type sigma_proc = exp(log_sigma_proc);
  
  vector<Type> proc_errors(time - 1);
  
  proc_errors = exp(log_proc_errors - pow(sigma_proc,2)/2);
  
  // proc_errors = exp(log_proc_errors * sigma_proc - pow(sigma_proc,2)/2);
  
  
  Type k = exp(log_k);
  
  catch_t = catch_t + Type(1e-3);
  
  Type r = exp(log_r);
  
  Type q = exp(log_q);
  
  Type m = exp(log_shape);
  
  Type b_to_k = pow(m, (-1 / (m - 1)));
  
  if (plim > b_to_k){
    
    plim = b_to_k;
    
  }
  
  
  Type init_dep = exp(log_init_dep);
  
  Type bmsy = k*pow(m, -1/(m - 1));
  
  Type umsy = (r / (m - 1)) * (1 - 1/m);
  
  Type msy = bmsy * umsy;
  
  b_t(0) = init_dep;
  
  q_t(0) = q;
  
  for (int t = 1; t<time; t++){
    
    q_t[t] = q_t[t - 1] * (1 + q_slope);
    
    u_t(t - 1) = catch_t(t - 1) / (b_t(t - 1) * k);
    
    growth_t(t - 1) = growfoo(r,m,b_t(t - 1),plim);
    
    b_t(t) = (b_t(t - 1) +  growth_t(t - 1) - catch_t(t - 1) / k) * proc_errors(t - 1);
    
    b_t(t) = posfun(b_t(t),eps,pen);
    
    
  } // close population model
  
  
  u_t(time - 1) = catch_t(time -1)  / (b_t(time - 1) * k);
  
  growth_t(time - 1) = (r  / (m - 1)) * b_t(time - 1) * (1 - pow(b_t(time - 1),m - 1));
  
  b_t = b_t * k;
  
  b_v_bmsy = b_t / bmsy;
  
  u_v_umsy = u_t / umsy;
  
  dep_t = b_t / k;
  
  if (calc_cpue == 0){
    index_hat_t = q_t * b_t;
  } else {
    
    index_hat_t = b_t;
  }
  
  nll += pen;
  
  // umsy penalty
  
  if (umsy > 1){
    
    nll += pow(umsy - 1, 2);
    
  }
  
  
  if (fit_index == 1) {
    
    vector<Type> z_t = index_t / b_t(index_years - 1);
    
    Type q_hat = z_t.mean();
    
    for (int t = 0; t < index_years.size(); t++){
      
      if (calc_cpue == 1){
        
        Type ftemp = q_t(index_years(t) - 1) * effort_t(t);
        
        Type effective_f = (ftemp / (ftemp + nat_m)) * (1 - exp(-(ftemp + nat_m)));
        
        if (use_baranov == 1){
        
        index_t(t) = catch_t(index_years(t) - 1) / effective_f;
          
        } else {
          index_t(t) = catch_t(index_years(t) - 1) / ftemp;
        }
        
      }
      
      if (marginalize_q == 0){ // standard estimation of q
        
          // test bias correction on the sigma_obs
        nll -= dnorm(log(index_t(t) + 1e-3), log(index_hat_t(index_years(t) - 1) + 1e-3), sigma_obs, true);
        
      } else { // marginalized q 
        
        nll -= dnorm(log(z_t(t) + 1e-3), log(q_hat + 1e-3), sigma_obs, true);
        
      }
    }
    
  }
  
  
  if (use_u_prior == 1){
    
    
    for (int t = 0; t < u_years.size(); t++){
      
      nll -= dnorm(log(u_v_umsy(u_years(t) - 1) + 1e-3), log(u_priors(t) + 1e-3), sigma_u, true);
      
    }
    
  }
  
  for (int t = 0; t < (time - 1); t++){
    
    nll -= dnorm(log_proc_errors(t), -pow(sigma_proc,2)/2, sigma_proc, true);
    
    // nll -= dnorm(log_proc_errors(t), Type(0), Type(1), true);
    
    
  }
  
  
  nll -= dnorm(log_r, log_r_prior, log_r_cv, true);
  
  Type final_ref = 0.5;
  
  Type init_ref = 1.0;
  
  if (b_ref_type == 0){
    
    final_ref = dep_t(time - 1);
    
    init_ref = dep_t(0);
    
  }
  
  if (b_ref_type == 1) {
    
    final_ref = b_v_bmsy(time - 1);
    
    init_ref = b_v_bmsy(0);
    
  }
  
  
  if (use_final_state == 1) {
    
    nll -= dnorm(log(final_ref), log_final_dep_prior,log_final_dep_cv, true);
    
  }
  
  // allow for arbitarry numbers of priors on final U/Umsy, so that you can use FMI + SAR
  if (use_final_u == 1){
    
    for (int i = 0; i < log_final_u.size(); i++){
      
      if (f_ref_type == 1){
      
      nll -= dnorm(log(u_v_umsy(time - 1)), log_final_u(i),log_final_u_cv(i), true);
      
      } 
      
      if (f_ref_type == 0){
        
        // std::cout<< "hello? it's u!" << std::endl;
      
      nll -= dnorm(log(u_t(time - 1)), log_final_u(i),log_final_u_cv(i), true);
      
        
      }
    } // cluse use final u loops
    
  } // close use final u
  
  nll -= dnorm(log(init_ref), log_init_dep_prior, log_init_dep_cv, true);
  
  nll -= dnorm(log_sigma_obs,log(sigma_obs_prior),sigma_obs_prior_cv, true);
  
  nll -= dnorm(log_q_slope,log(q_slope_prior),q_slope_cv, true);
  
  nll -= dnorm(log_sigma_proc,log(sigma_proc_prior), sigma_proc_prior_cv, true);
  
  nll -= dnorm(log_k,log(k_prior),Type(k_prior_cv), true);
  
  nll -= dnorm(log_shape,log(shape_prior),shape_cv, true);
  
  vector<Type> log_bt = log(b_t);
  
  vector<Type> log_b = log(b_v_bmsy);
  
  vector<Type> log_dep = log(dep_t);
  
  vector<Type> log_u = log(u_v_umsy);
  
  vector<Type> log_c_div_msy = log(catch_t / msy);
  
  vector<Type> log_ihat = log(index_hat_t);
  
  vector<Type> ck = catch_t / k;
  
  REPORT(umsy);
  
  REPORT(ck);
  
  REPORT(log_ihat);
  
  ADREPORT(log_ihat);
  
  ADREPORT(log_b);
  
  ADREPORT(log_bt);
  
  ADREPORT(log_dep);
  
  ADREPORT(log_u);
  
  ADREPORT(log_c_div_msy);
  
  REPORT(log_b);
  
  REPORT(log_bt);
  
  REPORT(log_dep);
  
  REPORT(log_u);
  
  REPORT(log_c_div_msy);
  
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
  
  ADREPORT(q_t);
  
  REPORT(r);
  
  REPORT(m);
  
  REPORT(pen);
  
  ADREPORT(q_slope);
  
  
  return nll;
}
