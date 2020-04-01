#include <TMB.hpp>

template<class Type>
Type growfoo(Type r, Type m, Type b, Type plim)
{
  
  Type growth = (r  / (m - 1)) * b * (1 - pow(b,m - 1));
  
  Type plim_growth = (r  / (m - 1)) * b * (1 - pow(b,m - 1)) * (b / (plim));
  
  return  CppAD::CondExpGe(b, plim, growth, plim_growth);
  
} // close function



template<class Type>
Type posfun(Type x, Type eps, Type &pen)
{
  
  pen += CppAD::CondExpLt(x, eps, Type(0.01) * pow(x-eps,2), Type(0));
  return CppAD::CondExpGe(x, eps, x, eps/(Type(2)-x/eps));
}


template<class Type>
vector<Type> popmodel(Type r, Type m, Type k, Type b0, vector<Type> catches, vector<Type> proc_error, int time, Type eps, Type plim)
{
  
  vector<Type> b(time);
  
  b(0) = b0;
  
  for (int t = 1; t < time; t++){
    
    
    Type growth = growfoo(r, m, b(t - 1),plim);
    
    Type temp = ((b(t - 1) + growth - catches(t - 1) / k) * proc_error(t - 1));
    
    b(t) = CppAD::CondExpGe(temp * k, eps, temp, (eps/(Type(2)-(temp * k)/eps)) / k);
    
  }
  b = b * k;
  
  return b;
  
  
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
  
  DATA_INTEGER(analytical_q); // should q be marginalized out? only use for fitting index not CPUE
  
  DATA_INTEGER(use_u_prior); //
  
  DATA_INTEGER(b_ref_type); //0 means that initial and terminal are relative to K, 1 to Bmsy
  
  DATA_INTEGER(f_ref_type); //0 means f, 1 f/fmsy
  
  DATA_INTEGER(est_k); // 1 means estiamte k, 0 means estimate terminal depletion
  
  DATA_INTEGER(estimate_proc_error); // 1 to estimate process error, 0 to not
  
  DATA_INTEGER(estimate_shape); 
  
  DATA_INTEGER(estimate_qslope);
  
  DATA_INTEGER(estimate_f);
  
  DATA_INTEGER(f_prior_form);
  
  // DATA_INTEGER(use_init); // use initial reference point
  
  DATA_INTEGER(use_terminal_state); // use terminal reference point
  
  DATA_INTEGER(use_terminal_u); // use terminal U/Umsy
  
  DATA_VECTOR(log_terminal_u);
  
  DATA_VECTOR(log_terminal_u_cv);
  
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
  
  DATA_SCALAR(sigma_ratio_prior);
  
  DATA_SCALAR(sigma_ratio_prior_cv);
  
  DATA_SCALAR(eps);
  
  DATA_SCALAR(learn_rate);
  
  DATA_SCALAR(log_r_cv);
  
  DATA_SCALAR(sigma_u);
  
  DATA_SCALAR(log_init_dep_prior);
  
  DATA_SCALAR(log_init_dep_cv);
  
  DATA_SCALAR(log_terminal_dep_prior);
  
  DATA_SCALAR(log_terminal_dep_cv);
  
  DATA_SCALAR(q_slope_prior);
  
  DATA_SCALAR(q_slope_cv);
  
  DATA_SCALAR(sigma_obs_prior);
  
  DATA_SCALAR(sigma_obs_prior_cv);
  
  
  //// parameters ////
  
  PARAMETER(log_init_dep);
  
  PARAMETER(log_r);
  
  PARAMETER(log_anchor);
  
  PARAMETER(log_q);
  
  PARAMETER(log_q_slope);
  
  PARAMETER(log_sigma_obs);
  
  // PARAMETER(sigma_obs);
  
  
  PARAMETER(log_sigma_ratio);
  
  PARAMETER_VECTOR(log_proc_errors);
  
  PARAMETER(log_shape);
  
  PARAMETER_VECTOR(log_f_t);
  
  
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
  
  vector<Type> f_t(time);
  
  vector<Type> catch_hat_t(time);
  
  vector<Type> temp_catch_t(time);
  
  Type q_slope = 0;
  
  if (estimate_qslope == 1){
    
    q_slope = exp(log_q_slope);
    
  } 
  
  if (estimate_f == 1){
    
    f_t = exp(log_f_t);
    
  }
  
  Type sigma_obs = exp(log_sigma_obs);
  
  // Type sigma_proc = exp(log_sigma_obs) * exp(log_sigma_ratio);
  
  Type sigma_proc = sigma_obs * exp(log_sigma_ratio);
  
  vector<Type> proc_errors(time - 1);
  
  proc_errors = exp(log_proc_errors);
  
  // proc_errors = exp(log_proc_errors - pow(sigma_proc,2)/2);
  
  // proc_errors = exp(log_proc_errors * sigma_proc - pow(sigma_proc,2)/2);
  
  
  
  catch_t = catch_t;
  
  Type r = exp(log_r);
  
  Type q = exp(log_q);
  
  Type m = exp(log_shape);
  
  Type k = 0;
  
  Type init_dep = exp(log_init_dep);
  
  Type delta = 100;
  
  Type counter = 0;
  
  Type new_proposal = 0;
  
  Type conv_error = 0;
  
  Type terminal_state = 0;
  
  Type grad_dep = 0;
  
  if (est_k == 1){
    
    k = exp(log_anchor);
    
  } else {
    
    // back out k form terminal status
    
    vector<Type> proposal_result(time);
    
    terminal_state = exp(log_anchor);
    
    Type last_proposal =   log(max(catch_t) * 10);
    
    proposal_result = popmodel(r,m,exp(last_proposal), init_dep, catch_t, proc_errors, time, eps,plim);
    
    Type prop_error =  log(proposal_result[time - 1] / exp(last_proposal)) - log(terminal_state);
    
    while(delta > 1e-3){
      
      new_proposal = last_proposal -  learn_rate * prop_error;
      
      // std::cout << "new proposal is" << new_proposal << std::endl;
      
      proposal_result = popmodel(r,m,exp(new_proposal), init_dep, catch_t, proc_errors, time, eps,plim);
      
      grad_dep = proposal_result[time - 1] / exp(new_proposal);
      
      prop_error =  log(grad_dep) - log(terminal_state);
      
      // std::cout << "prop dep is" << proposal_result[time - 1] / exp(new_proposal) << std::endl;
      
      
      last_proposal = new_proposal;
      
      delta = sqrt(pow(grad_dep - terminal_state,2));
      
      conv_error = delta;
      
      
      conv_error = delta;
      
      
      // std::cout<< delta << std::endl;
      
      counter = counter + 1;
      
      if (counter > 5000){
        
        delta = -999;
      } // close escape hatch
      
      
    } // close while looop for gradient descnet
    
    k = exp(new_proposal);
    
    nll += conv_error;  
  } // close estimate_k
  
  // std::cout<< counter << std::endl;
  
  Type b_to_k = pow(m, (-1 / (m - 1)));
  
  if (plim > b_to_k){
    
    plim = b_to_k;
    
  }
  
  
  Type bmsy = k*pow(m, -1/(m - 1));
  
  Type umsy = (r / (m - 1)) * (1 - 1/m);
  
  Type msy = bmsy * umsy;
  
  b_t(0) = init_dep;
  
  q_t(0) = q;
  
  
  if (estimate_f == 1){
    
    temp_catch_t(0) = b_t(0) * (1 - exp(-f_t(0))) * k;

    
  } else {
    temp_catch_t = catch_t;
  }
  
  for (int t = 1; t<time; t++){
    
    q_t[t] = q_t[t - 1] * (1 + q_slope);
    
    u_t(t - 1) = catch_t(t - 1) / (b_t(t - 1) * k);
    
    growth_t(t - 1) = growfoo(r,m,b_t(t - 1),plim);
    
    b_t(t) = (b_t(t - 1) +  growth_t(t - 1) - temp_catch_t(t - 1) / k) * proc_errors(t - 1);
    
    // b_t(t) = posfun(b_t(t) * k,eps,pen) / k;
    
    b_t(t) = posfun(b_t(t),eps,pen);
    
    if (estimate_f == 1){
      temp_catch_t(t) = b_t(t) * (1 - exp(-f_t(t))) * k;
    }
    
  } // close population model
  
  
  u_t(time - 1) = catch_t(time -1)  / (b_t(time - 1) * k);
  
  if (estimate_f == 1){
    
    u_t = 1 - exp(-f_t);
    
  }
  
  growth_t(time - 1) = (r  / (m - 1)) * b_t(time - 1) * (1 - pow(b_t(time - 1),m - 1));
  
  b_t = b_t * k;
  
  b_v_bmsy = b_t / bmsy;
  
  u_v_umsy = u_t / umsy;
  
  dep_t = b_t / k;
  
  // std::cout << "grad dep is" << grad_dep << "pop dep is" << dep_t(time - 1) << std::endl;
  
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
  
  
  if (estimate_f == 1){
    
    for (int t = 0; t < time; t++){
      
      if (f_prior_form == 1 & use_u_prior == 0 & t > 0){
        
        nll -= dnorm(log_f_t(t), log_f_t(t - 1), sigma_u);
        
      }
       
      nll -= dnorm(log(catch_t(t) + 1e-3), log(temp_catch_t(t) + 1e-3), Type(0.01), true);
      
    }
    
  }
  
  if (fit_index == 1) {
    
    nll -= dnorm(log_sigma_obs,log(sigma_obs_prior), sigma_obs_prior_cv, true);
    
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
      
      if (analytical_q == 0){ // standard estimation of q
        
        // test bias correction on the sigma_obs
        nll -= dnorm(log(index_t(t)), log(index_hat_t(index_years(t) - 1)), sigma_obs, true);
        
      } else { // marginalized q 
        
        nll -= dnorm(log(z_t(t) + 1e-3), log(q_hat + 1e-3), sigma_obs, true);
        
      }
    }
    
  }
  
  
  if (use_u_prior == 1){
    
    
    for (int t = 0; t < u_years.size(); t++){
      
      
      if (f_ref_type == 1){
        
        nll -= dnorm(log(u_v_umsy(u_years(t) - 1) + 1e-3), log(u_priors(t) + 1e-3), sigma_u, true);
      
      } 
      
      if (f_ref_type == 0){
        
        // std::cout<< "hello? it's u!" << std::endl;
        
          nll -= dnorm(log(u_t(u_years(t) - 1) + 1e-3), log(u_priors(t) + 1e-3), sigma_u, true);
          
        
      }
    }
    
  }
  
  if (estimate_proc_error == 1){
    
    for (int t = 0; t < (time - 1); t++){
      
      nll -= dnorm(log_proc_errors(t), -pow(sigma_proc,2)/2, sigma_proc, true);
      
      // nll -= dnorm(log_proc_errors(t), Type(0), Type(1), true);
      
      
    }
  }
  
  nll -= dnorm(log_r, log_r_prior, log_r_cv, true);
  
  Type terminal_ref = 0.5;
  
  Type init_ref = 1.0;
  
  if (b_ref_type == 0){
    
    terminal_ref = dep_t(time - 1);
    
    init_ref = dep_t(0);
    
  }
  
  if (b_ref_type == 1) {
    
    terminal_ref = b_v_bmsy(time - 1);
    
    init_ref = b_v_bmsy(0);
    
  }
  
  
  if (use_terminal_state == 1) {
    
    nll -= dnorm(log(terminal_ref), log_terminal_dep_prior,log_terminal_dep_cv, true);
    
  }
  
  // allow for arbitarry numbers of priors on terminal U/Umsy, so that you can use FMI + SAR
  if (use_terminal_u == 1){
    
    for (int i = 0; i < log_terminal_u.size(); i++){
      
      if (f_ref_type == 1){
        
        nll -= dnorm(log(u_v_umsy(time - 1)), log_terminal_u(i),log_terminal_u_cv(i), true);
        
      } 
      
      if (f_ref_type == 0){
        
        // std::cout<< "hello? it's u!" << std::endl;
        
        nll -= dnorm(log(u_t(time - 1)), log_terminal_u(i),log_terminal_u_cv(i), true);
        
      }
    } // cluse use terminal u loops
    
  } // close use terminal u
  
  nll -= dnorm(log_init_dep, log_init_dep_prior, log_init_dep_cv, true);
  
  
  if (estimate_qslope == 1){
    
    nll -= dnorm(log_q_slope,log(q_slope_prior),q_slope_cv, true);
    
  }
  
  if (estimate_proc_error == 1){
    
    nll -= dnorm(log_sigma_ratio, log(sigma_ratio_prior), sigma_ratio_prior_cv, true);
    
  }
  // nll -= dnorm(log_sigma_proc,log(sigma_proc_prior), sigma_proc_prior_cv, true);
  
  if (est_k == 1){
    nll -= dnorm(log_anchor,log(k_prior),k_prior_cv, true);
  }
  
  if (estimate_shape == 1){
    
    nll -= dnorm(log_shape,log(shape_prior),shape_cv, true);
    
  }
  
  vector<Type> log_bt = log(b_t);
  
  vector<Type> log_b = log(b_v_bmsy);
  
  vector<Type> log_dep = log(dep_t);
  
  vector<Type> log_u = log(u_v_umsy);
  
  vector<Type> log_c_div_msy = log(catch_t / msy);
  
  vector<Type> log_ihat = log(index_hat_t);
  
  vector<Type> ck = catch_t / k;
  
  REPORT(sigma_proc);
  
  ADREPORT(sigma_proc);
  
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
  
  REPORT(terminal_ref);
  
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
  
  REPORT(delta);
  
  ADREPORT(q_slope);
  
  
  return nll;
}
