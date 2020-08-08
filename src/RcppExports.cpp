// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// sraplus
List sraplus(NumericVector catches, NumericVector rs, NumericVector ms, NumericVector init_deps, NumericVector anchors, NumericVector qs, NumericVector sigma_procs, NumericVector drawdex, NumericVector index_t, NumericVector sigma_obs, NumericVector log_terminal_u, NumericVector log_terminal_u_cv, IntegerVector index_years, IntegerVector u_years, int draws, int n_keep, int b_ref_type, int f_ref_type, int fit_index, int use_terminal_u, int use_terminal_state, bool estimate_k, double log_terminal_ref, double sigma_dep, double plim, int use_u_prior, NumericVector u_priors, double sigma_u, double learn_rate);
RcppExport SEXP _sraplus_sraplus(SEXP catchesSEXP, SEXP rsSEXP, SEXP msSEXP, SEXP init_depsSEXP, SEXP anchorsSEXP, SEXP qsSEXP, SEXP sigma_procsSEXP, SEXP drawdexSEXP, SEXP index_tSEXP, SEXP sigma_obsSEXP, SEXP log_terminal_uSEXP, SEXP log_terminal_u_cvSEXP, SEXP index_yearsSEXP, SEXP u_yearsSEXP, SEXP drawsSEXP, SEXP n_keepSEXP, SEXP b_ref_typeSEXP, SEXP f_ref_typeSEXP, SEXP fit_indexSEXP, SEXP use_terminal_uSEXP, SEXP use_terminal_stateSEXP, SEXP estimate_kSEXP, SEXP log_terminal_refSEXP, SEXP sigma_depSEXP, SEXP plimSEXP, SEXP use_u_priorSEXP, SEXP u_priorsSEXP, SEXP sigma_uSEXP, SEXP learn_rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type catches(catchesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rs(rsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ms(msSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type init_deps(init_depsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type anchors(anchorsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type qs(qsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma_procs(sigma_procsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type drawdex(drawdexSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type index_t(index_tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma_obs(sigma_obsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type log_terminal_u(log_terminal_uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type log_terminal_u_cv(log_terminal_u_cvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type index_years(index_yearsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type u_years(u_yearsSEXP);
    Rcpp::traits::input_parameter< int >::type draws(drawsSEXP);
    Rcpp::traits::input_parameter< int >::type n_keep(n_keepSEXP);
    Rcpp::traits::input_parameter< int >::type b_ref_type(b_ref_typeSEXP);
    Rcpp::traits::input_parameter< int >::type f_ref_type(f_ref_typeSEXP);
    Rcpp::traits::input_parameter< int >::type fit_index(fit_indexSEXP);
    Rcpp::traits::input_parameter< int >::type use_terminal_u(use_terminal_uSEXP);
    Rcpp::traits::input_parameter< int >::type use_terminal_state(use_terminal_stateSEXP);
    Rcpp::traits::input_parameter< bool >::type estimate_k(estimate_kSEXP);
    Rcpp::traits::input_parameter< double >::type log_terminal_ref(log_terminal_refSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_dep(sigma_depSEXP);
    Rcpp::traits::input_parameter< double >::type plim(plimSEXP);
    Rcpp::traits::input_parameter< int >::type use_u_prior(use_u_priorSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u_priors(u_priorsSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_u(sigma_uSEXP);
    Rcpp::traits::input_parameter< double >::type learn_rate(learn_rateSEXP);
    rcpp_result_gen = Rcpp::wrap(sraplus(catches, rs, ms, init_deps, anchors, qs, sigma_procs, drawdex, index_t, sigma_obs, log_terminal_u, log_terminal_u_cv, index_years, u_years, draws, n_keep, b_ref_type, f_ref_type, fit_index, use_terminal_u, use_terminal_state, estimate_k, log_terminal_ref, sigma_dep, plim, use_u_prior, u_priors, sigma_u, learn_rate));
    return rcpp_result_gen;
END_RCPP
}
// sraplus_reserve
List sraplus_reserve(NumericVector catches, NumericVector rs, NumericVector ms, NumericVector init_deps, NumericVector ks, NumericVector qs, NumericVector sigma_procs, NumericVector drawdex, NumericVector index_t, NumericVector sigma_obs, NumericVector log_final_u, NumericVector log_final_u_cv, IntegerVector index_years, IntegerVector u_years, int draws, int n_keep, int b_ref_type, int f_ref_type, int fit_index, int use_final_u, int use_final_state, double log_final_ref, double sigma_dep, double plim, int use_u_prior, NumericVector u_priors, double sigma_u);
RcppExport SEXP _sraplus_sraplus_reserve(SEXP catchesSEXP, SEXP rsSEXP, SEXP msSEXP, SEXP init_depsSEXP, SEXP ksSEXP, SEXP qsSEXP, SEXP sigma_procsSEXP, SEXP drawdexSEXP, SEXP index_tSEXP, SEXP sigma_obsSEXP, SEXP log_final_uSEXP, SEXP log_final_u_cvSEXP, SEXP index_yearsSEXP, SEXP u_yearsSEXP, SEXP drawsSEXP, SEXP n_keepSEXP, SEXP b_ref_typeSEXP, SEXP f_ref_typeSEXP, SEXP fit_indexSEXP, SEXP use_final_uSEXP, SEXP use_final_stateSEXP, SEXP log_final_refSEXP, SEXP sigma_depSEXP, SEXP plimSEXP, SEXP use_u_priorSEXP, SEXP u_priorsSEXP, SEXP sigma_uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type catches(catchesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rs(rsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ms(msSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type init_deps(init_depsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ks(ksSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type qs(qsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma_procs(sigma_procsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type drawdex(drawdexSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type index_t(index_tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma_obs(sigma_obsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type log_final_u(log_final_uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type log_final_u_cv(log_final_u_cvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type index_years(index_yearsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type u_years(u_yearsSEXP);
    Rcpp::traits::input_parameter< int >::type draws(drawsSEXP);
    Rcpp::traits::input_parameter< int >::type n_keep(n_keepSEXP);
    Rcpp::traits::input_parameter< int >::type b_ref_type(b_ref_typeSEXP);
    Rcpp::traits::input_parameter< int >::type f_ref_type(f_ref_typeSEXP);
    Rcpp::traits::input_parameter< int >::type fit_index(fit_indexSEXP);
    Rcpp::traits::input_parameter< int >::type use_final_u(use_final_uSEXP);
    Rcpp::traits::input_parameter< int >::type use_final_state(use_final_stateSEXP);
    Rcpp::traits::input_parameter< double >::type log_final_ref(log_final_refSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_dep(sigma_depSEXP);
    Rcpp::traits::input_parameter< double >::type plim(plimSEXP);
    Rcpp::traits::input_parameter< int >::type use_u_prior(use_u_priorSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u_priors(u_priorsSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_u(sigma_uSEXP);
    rcpp_result_gen = Rcpp::wrap(sraplus_reserve(catches, rs, ms, init_deps, ks, qs, sigma_procs, drawdex, index_t, sigma_obs, log_final_u, log_final_u_cv, index_years, u_years, draws, n_keep, b_ref_type, f_ref_type, fit_index, use_final_u, use_final_state, log_final_ref, sigma_dep, plim, use_u_prior, u_priors, sigma_u));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_sraplus_sraplus", (DL_FUNC) &_sraplus_sraplus, 29},
    {"_sraplus_sraplus_reserve", (DL_FUNC) &_sraplus_sraplus_reserve, 27},
    {NULL, NULL, 0}
};

RcppExport void R_init_sraplus(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
