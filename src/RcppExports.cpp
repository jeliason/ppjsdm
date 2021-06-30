// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// compute_papangelou_cpp
Rcpp::NumericVector compute_papangelou_cpp(SEXP configuration, Rcpp::NumericVector x, Rcpp::NumericVector y, Rcpp::CharacterVector model, Rcpp::CharacterVector medium_range_model, Rcpp::NumericMatrix alpha, Rcpp::NumericVector beta0, Rcpp::NumericMatrix beta, Rcpp::NumericMatrix gamma, Rcpp::List covariates, Rcpp::NumericMatrix short_range, Rcpp::NumericMatrix medium_range, Rcpp::NumericMatrix long_range, R_xlen_t saturation, R_xlen_t type, double mark);
RcppExport SEXP _ppjsdm_compute_papangelou_cpp(SEXP configurationSEXP, SEXP xSEXP, SEXP ySEXP, SEXP modelSEXP, SEXP medium_range_modelSEXP, SEXP alphaSEXP, SEXP beta0SEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP covariatesSEXP, SEXP short_rangeSEXP, SEXP medium_rangeSEXP, SEXP long_rangeSEXP, SEXP saturationSEXP, SEXP typeSEXP, SEXP markSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type configuration(configurationSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type model(modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type medium_range_model(medium_range_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type covariates(covariatesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type short_range(short_rangeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type medium_range(medium_rangeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type long_range(long_rangeSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type saturation(saturationSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type type(typeSEXP);
    Rcpp::traits::input_parameter< double >::type mark(markSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_papangelou_cpp(configuration, x, y, model, medium_range_model, alpha, beta0, beta, gamma, covariates, short_range, medium_range, long_range, saturation, type, mark));
    return rcpp_result_gen;
END_RCPP
}
// compute_vcov
Rcpp::List compute_vcov(SEXP configuration, SEXP window, Rcpp::List covariates, Rcpp::CharacterVector model, Rcpp::CharacterVector medium_range_model, Rcpp::NumericMatrix short_range, Rcpp::NumericMatrix medium_range, Rcpp::NumericMatrix long_range, R_xlen_t saturation, Rcpp::NumericVector rho, Rcpp::NumericVector theta, Rcpp::NumericMatrix regressors, Rcpp::List data_list, Rcpp::LogicalMatrix estimate_alpha, Rcpp::LogicalMatrix estimate_gamma, bool debug, int nthreads, int npoints);
RcppExport SEXP _ppjsdm_compute_vcov(SEXP configurationSEXP, SEXP windowSEXP, SEXP covariatesSEXP, SEXP modelSEXP, SEXP medium_range_modelSEXP, SEXP short_rangeSEXP, SEXP medium_rangeSEXP, SEXP long_rangeSEXP, SEXP saturationSEXP, SEXP rhoSEXP, SEXP thetaSEXP, SEXP regressorsSEXP, SEXP data_listSEXP, SEXP estimate_alphaSEXP, SEXP estimate_gammaSEXP, SEXP debugSEXP, SEXP nthreadsSEXP, SEXP npointsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type configuration(configurationSEXP);
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type covariates(covariatesSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type model(modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type medium_range_model(medium_range_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type short_range(short_rangeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type medium_range(medium_rangeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type long_range(long_rangeSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type saturation(saturationSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type regressors(regressorsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type data_list(data_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalMatrix >::type estimate_alpha(estimate_alphaSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalMatrix >::type estimate_gamma(estimate_gammaSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< int >::type npoints(npointsSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_vcov(configuration, window, covariates, model, medium_range_model, short_range, medium_range, long_range, saturation, rho, theta, regressors, data_list, estimate_alpha, estimate_gamma, debug, nthreads, npoints));
    return rcpp_result_gen;
END_RCPP
}
// has_duplicates
bool has_duplicates(Rcpp::List configuration);
RcppExport SEXP _ppjsdm_has_duplicates(SEXP configurationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type configuration(configurationSEXP);
    rcpp_result_gen = Rcpp::wrap(has_duplicates(configuration));
    return rcpp_result_gen;
END_RCPP
}
// make_default_model_parameters
SEXP make_default_model_parameters(SEXP alpha, SEXP beta0, SEXP covariates, SEXP beta, SEXP gamma, SEXP short_range, SEXP medium_range, SEXP long_range, SEXP types, int default_number_types);
RcppExport SEXP _ppjsdm_make_default_model_parameters(SEXP alphaSEXP, SEXP beta0SEXP, SEXP covariatesSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP short_rangeSEXP, SEXP medium_rangeSEXP, SEXP long_rangeSEXP, SEXP typesSEXP, SEXP default_number_typesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< SEXP >::type covariates(covariatesSEXP);
    Rcpp::traits::input_parameter< SEXP >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type short_range(short_rangeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type medium_range(medium_rangeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type long_range(long_rangeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type types(typesSEXP);
    Rcpp::traits::input_parameter< int >::type default_number_types(default_number_typesSEXP);
    rcpp_result_gen = Rcpp::wrap(make_default_model_parameters(alpha, beta0, covariates, beta, gamma, short_range, medium_range, long_range, types, default_number_types));
    return rcpp_result_gen;
END_RCPP
}
// prepare_gibbsm_data
Rcpp::List prepare_gibbsm_data(Rcpp::List configuration_list, SEXP window, Rcpp::List covariates, Rcpp::CharacterVector model, Rcpp::CharacterVector medium_range_model, SEXP short_range, SEXP medium_range, SEXP long_range, R_xlen_t saturation, Rcpp::NumericVector mark_range, R_xlen_t max_dummy, R_xlen_t min_dummy, double dummy_factor, Rcpp::LogicalMatrix estimate_alpha, Rcpp::LogicalMatrix estimate_gamma, int nthreads, bool debug, std::string dummy_distribution, Rcpp::CharacterVector type_names);
RcppExport SEXP _ppjsdm_prepare_gibbsm_data(SEXP configuration_listSEXP, SEXP windowSEXP, SEXP covariatesSEXP, SEXP modelSEXP, SEXP medium_range_modelSEXP, SEXP short_rangeSEXP, SEXP medium_rangeSEXP, SEXP long_rangeSEXP, SEXP saturationSEXP, SEXP mark_rangeSEXP, SEXP max_dummySEXP, SEXP min_dummySEXP, SEXP dummy_factorSEXP, SEXP estimate_alphaSEXP, SEXP estimate_gammaSEXP, SEXP nthreadsSEXP, SEXP debugSEXP, SEXP dummy_distributionSEXP, SEXP type_namesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type configuration_list(configuration_listSEXP);
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type covariates(covariatesSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type model(modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type medium_range_model(medium_range_modelSEXP);
    Rcpp::traits::input_parameter< SEXP >::type short_range(short_rangeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type medium_range(medium_rangeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type long_range(long_rangeSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type saturation(saturationSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mark_range(mark_rangeSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type max_dummy(max_dummySEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type min_dummy(min_dummySEXP);
    Rcpp::traits::input_parameter< double >::type dummy_factor(dummy_factorSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalMatrix >::type estimate_alpha(estimate_alphaSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalMatrix >::type estimate_gamma(estimate_gammaSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    Rcpp::traits::input_parameter< std::string >::type dummy_distribution(dummy_distributionSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type type_names(type_namesSEXP);
    rcpp_result_gen = Rcpp::wrap(prepare_gibbsm_data(configuration_list, window, covariates, model, medium_range_model, short_range, medium_range, long_range, saturation, mark_range, max_dummy, min_dummy, dummy_factor, estimate_alpha, estimate_gamma, nthreads, debug, dummy_distribution, type_names));
    return rcpp_result_gen;
END_RCPP
}
// rbinomialpp_cpp
SEXP rbinomialpp_cpp(SEXP window, SEXP n, R_xlen_t nsim, SEXP types, bool drop, Rcpp::NumericVector mark_range);
RcppExport SEXP _ppjsdm_rbinomialpp_cpp(SEXP windowSEXP, SEXP nSEXP, SEXP nsimSEXP, SEXP typesSEXP, SEXP dropSEXP, SEXP mark_rangeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< SEXP >::type n(nSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< SEXP >::type types(typesSEXP);
    Rcpp::traits::input_parameter< bool >::type drop(dropSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mark_range(mark_rangeSEXP);
    rcpp_result_gen = Rcpp::wrap(rbinomialpp_cpp(window, n, nsim, types, drop, mark_range));
    return rcpp_result_gen;
END_RCPP
}
// rgibbs_cpp
SEXP rgibbs_cpp(SEXP window, SEXP alpha, Rcpp::NumericVector beta0, SEXP covariates, SEXP beta, SEXP gamma, SEXP short_range, SEXP medium_range, SEXP long_range, R_xlen_t saturation, R_xlen_t steps, R_xlen_t nsim, SEXP types, Rcpp::CharacterVector model, Rcpp::CharacterVector medium_range_model, bool drop, Rcpp::NumericVector mark_range, SEXP starting_configuration);
RcppExport SEXP _ppjsdm_rgibbs_cpp(SEXP windowSEXP, SEXP alphaSEXP, SEXP beta0SEXP, SEXP covariatesSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP short_rangeSEXP, SEXP medium_rangeSEXP, SEXP long_rangeSEXP, SEXP saturationSEXP, SEXP stepsSEXP, SEXP nsimSEXP, SEXP typesSEXP, SEXP modelSEXP, SEXP medium_range_modelSEXP, SEXP dropSEXP, SEXP mark_rangeSEXP, SEXP starting_configurationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< SEXP >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< SEXP >::type covariates(covariatesSEXP);
    Rcpp::traits::input_parameter< SEXP >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type short_range(short_rangeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type medium_range(medium_rangeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type long_range(long_rangeSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type saturation(saturationSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type steps(stepsSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< SEXP >::type types(typesSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type model(modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type medium_range_model(medium_range_modelSEXP);
    Rcpp::traits::input_parameter< bool >::type drop(dropSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mark_range(mark_rangeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type starting_configuration(starting_configurationSEXP);
    rcpp_result_gen = Rcpp::wrap(rgibbs_cpp(window, alpha, beta0, covariates, beta, gamma, short_range, medium_range, long_range, saturation, steps, nsim, types, model, medium_range_model, drop, mark_range, starting_configuration));
    return rcpp_result_gen;
END_RCPP
}
// rppp_cpp
SEXP rppp_cpp(SEXP window, SEXP lambda, R_xlen_t nsim, SEXP types, bool drop, Rcpp::NumericVector mark_range);
RcppExport SEXP _ppjsdm_rppp_cpp(SEXP windowSEXP, SEXP lambdaSEXP, SEXP nsimSEXP, SEXP typesSEXP, SEXP dropSEXP, SEXP mark_rangeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< SEXP >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< SEXP >::type types(typesSEXP);
    Rcpp::traits::input_parameter< bool >::type drop(dropSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mark_range(mark_rangeSEXP);
    rcpp_result_gen = Rcpp::wrap(rppp_cpp(window, lambda, nsim, types, drop, mark_range));
    return rcpp_result_gen;
END_RCPP
}
// rstratpp_cpp
SEXP rstratpp_cpp(SEXP window, SEXP delta_x, SEXP delta_y, R_xlen_t nsim, SEXP types, bool drop, Rcpp::NumericVector mark_range);
RcppExport SEXP _ppjsdm_rstratpp_cpp(SEXP windowSEXP, SEXP delta_xSEXP, SEXP delta_ySEXP, SEXP nsimSEXP, SEXP typesSEXP, SEXP dropSEXP, SEXP mark_rangeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< SEXP >::type delta_x(delta_xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type delta_y(delta_ySEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< SEXP >::type types(typesSEXP);
    Rcpp::traits::input_parameter< bool >::type drop(dropSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mark_range(mark_rangeSEXP);
    rcpp_result_gen = Rcpp::wrap(rstratpp_cpp(window, delta_x, delta_y, nsim, types, drop, mark_range));
    return rcpp_result_gen;
END_RCPP
}
// show_short_range_models
void show_short_range_models();
RcppExport SEXP _ppjsdm_show_short_range_models() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    show_short_range_models();
    return R_NilValue;
END_RCPP
}
// show_medium_range_models
void show_medium_range_models();
RcppExport SEXP _ppjsdm_show_medium_range_models() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    show_medium_range_models();
    return R_NilValue;
END_RCPP
}

RcppExport SEXP run_testthat_tests(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_ppjsdm_compute_papangelou_cpp", (DL_FUNC) &_ppjsdm_compute_papangelou_cpp, 16},
    {"_ppjsdm_compute_vcov", (DL_FUNC) &_ppjsdm_compute_vcov, 18},
    {"_ppjsdm_has_duplicates", (DL_FUNC) &_ppjsdm_has_duplicates, 1},
    {"_ppjsdm_make_default_model_parameters", (DL_FUNC) &_ppjsdm_make_default_model_parameters, 10},
    {"_ppjsdm_prepare_gibbsm_data", (DL_FUNC) &_ppjsdm_prepare_gibbsm_data, 19},
    {"_ppjsdm_rbinomialpp_cpp", (DL_FUNC) &_ppjsdm_rbinomialpp_cpp, 6},
    {"_ppjsdm_rgibbs_cpp", (DL_FUNC) &_ppjsdm_rgibbs_cpp, 18},
    {"_ppjsdm_rppp_cpp", (DL_FUNC) &_ppjsdm_rppp_cpp, 6},
    {"_ppjsdm_rstratpp_cpp", (DL_FUNC) &_ppjsdm_rstratpp_cpp, 7},
    {"_ppjsdm_show_short_range_models", (DL_FUNC) &_ppjsdm_show_short_range_models, 0},
    {"_ppjsdm_show_medium_range_models", (DL_FUNC) &_ppjsdm_show_medium_range_models, 0},
    {"run_testthat_tests", (DL_FUNC) &run_testthat_tests, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_ppjsdm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
