// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// compute_A2_A3
Rcpp::List compute_A2_A3(SEXP configuration, Rcpp::List covariates, Rcpp::CharacterVector model, Rcpp::CharacterVector medium_range_model, SEXP short_range, SEXP medium_range, SEXP long_range, R_xlen_t saturation, Rcpp::NumericMatrix alpha, Rcpp::NumericVector beta0, Rcpp::NumericMatrix beta, Rcpp::NumericMatrix gamma, double rho, Rcpp::List t_over_papangelou);
RcppExport SEXP _ppjsdm_compute_A2_A3(SEXP configurationSEXP, SEXP covariatesSEXP, SEXP modelSEXP, SEXP medium_range_modelSEXP, SEXP short_rangeSEXP, SEXP medium_rangeSEXP, SEXP long_rangeSEXP, SEXP saturationSEXP, SEXP alphaSEXP, SEXP beta0SEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP rhoSEXP, SEXP t_over_papangelouSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type configuration(configurationSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type covariates(covariatesSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type model(modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type medium_range_model(medium_range_modelSEXP);
    Rcpp::traits::input_parameter< SEXP >::type short_range(short_rangeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type medium_range(medium_rangeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type long_range(long_rangeSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type saturation(saturationSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type t_over_papangelou(t_over_papangelouSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_A2_A3(configuration, covariates, model, medium_range_model, short_range, medium_range, long_range, saturation, alpha, beta0, beta, gamma, rho, t_over_papangelou));
    return rcpp_result_gen;
END_RCPP
}
// compute_papangelou
Rcpp::NumericVector compute_papangelou(SEXP configuration, Rcpp::NumericVector x, Rcpp::NumericVector y, R_xlen_t type, Rcpp::CharacterVector model, Rcpp::CharacterVector medium_range_model, Rcpp::NumericMatrix alpha, Rcpp::NumericVector beta0, Rcpp::NumericMatrix beta, Rcpp::NumericMatrix gamma, Rcpp::List covariates, Rcpp::NumericMatrix short_range, Rcpp::NumericMatrix medium_range, Rcpp::NumericMatrix long_range, R_xlen_t saturation, double mark);
RcppExport SEXP _ppjsdm_compute_papangelou(SEXP configurationSEXP, SEXP xSEXP, SEXP ySEXP, SEXP typeSEXP, SEXP modelSEXP, SEXP medium_range_modelSEXP, SEXP alphaSEXP, SEXP beta0SEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP covariatesSEXP, SEXP short_rangeSEXP, SEXP medium_rangeSEXP, SEXP long_rangeSEXP, SEXP saturationSEXP, SEXP markSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type configuration(configurationSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type type(typeSEXP);
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
    Rcpp::traits::input_parameter< double >::type mark(markSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_papangelou(configuration, x, y, type, model, medium_range_model, alpha, beta0, beta, gamma, covariates, short_range, medium_range, long_range, saturation, mark));
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
SEXP make_default_model_parameters(SEXP alpha, SEXP beta0, SEXP covariates, SEXP beta, SEXP gamma, SEXP short_range, SEXP medium_range, SEXP long_range, SEXP types);
RcppExport SEXP _ppjsdm_make_default_model_parameters(SEXP alphaSEXP, SEXP beta0SEXP, SEXP covariatesSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP short_rangeSEXP, SEXP medium_rangeSEXP, SEXP long_rangeSEXP, SEXP typesSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(make_default_model_parameters(alpha, beta0, covariates, beta, gamma, short_range, medium_range, long_range, types));
    return rcpp_result_gen;
END_RCPP
}
// prepare_gibbsm_data
Rcpp::List prepare_gibbsm_data(Rcpp::List configuration_list, SEXP window, Rcpp::List covariates, Rcpp::List traits, Rcpp::CharacterVector model, Rcpp::CharacterVector medium_range_model, SEXP short_range, SEXP medium_range, SEXP long_range, R_xlen_t saturation, Rcpp::NumericVector mark_range, bool approximate, R_xlen_t ndummy);
RcppExport SEXP _ppjsdm_prepare_gibbsm_data(SEXP configuration_listSEXP, SEXP windowSEXP, SEXP covariatesSEXP, SEXP traitsSEXP, SEXP modelSEXP, SEXP medium_range_modelSEXP, SEXP short_rangeSEXP, SEXP medium_rangeSEXP, SEXP long_rangeSEXP, SEXP saturationSEXP, SEXP mark_rangeSEXP, SEXP approximateSEXP, SEXP ndummySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type configuration_list(configuration_listSEXP);
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type covariates(covariatesSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type traits(traitsSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type model(modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type medium_range_model(medium_range_modelSEXP);
    Rcpp::traits::input_parameter< SEXP >::type short_range(short_rangeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type medium_range(medium_rangeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type long_range(long_rangeSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type saturation(saturationSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mark_range(mark_rangeSEXP);
    Rcpp::traits::input_parameter< bool >::type approximate(approximateSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type ndummy(ndummySEXP);
    rcpp_result_gen = Rcpp::wrap(prepare_gibbsm_data(configuration_list, window, covariates, traits, model, medium_range_model, short_range, medium_range, long_range, saturation, mark_range, approximate, ndummy));
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
SEXP rgibbs_cpp(SEXP window, SEXP alpha, Rcpp::NumericVector beta0, SEXP covariates, SEXP beta, SEXP gamma, SEXP short_range, SEXP medium_range, SEXP long_range, R_xlen_t saturation, R_xlen_t steps, R_xlen_t nsim, SEXP types, Rcpp::CharacterVector model, Rcpp::CharacterVector medium_range_model, bool drop, Rcpp::NumericVector mark_range);
RcppExport SEXP _ppjsdm_rgibbs_cpp(SEXP windowSEXP, SEXP alphaSEXP, SEXP beta0SEXP, SEXP covariatesSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP short_rangeSEXP, SEXP medium_rangeSEXP, SEXP long_rangeSEXP, SEXP saturationSEXP, SEXP stepsSEXP, SEXP nsimSEXP, SEXP typesSEXP, SEXP modelSEXP, SEXP medium_range_modelSEXP, SEXP dropSEXP, SEXP mark_rangeSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(rgibbs_cpp(window, alpha, beta0, covariates, beta, gamma, short_range, medium_range, long_range, saturation, steps, nsim, types, model, medium_range_model, drop, mark_range));
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

static const R_CallMethodDef CallEntries[] = {
    {"_ppjsdm_compute_A2_A3", (DL_FUNC) &_ppjsdm_compute_A2_A3, 14},
    {"_ppjsdm_compute_papangelou", (DL_FUNC) &_ppjsdm_compute_papangelou, 16},
    {"_ppjsdm_has_duplicates", (DL_FUNC) &_ppjsdm_has_duplicates, 1},
    {"_ppjsdm_make_default_model_parameters", (DL_FUNC) &_ppjsdm_make_default_model_parameters, 9},
    {"_ppjsdm_prepare_gibbsm_data", (DL_FUNC) &_ppjsdm_prepare_gibbsm_data, 13},
    {"_ppjsdm_rbinomialpp_cpp", (DL_FUNC) &_ppjsdm_rbinomialpp_cpp, 6},
    {"_ppjsdm_rgibbs_cpp", (DL_FUNC) &_ppjsdm_rgibbs_cpp, 17},
    {"_ppjsdm_rppp_cpp", (DL_FUNC) &_ppjsdm_rppp_cpp, 6},
    {"_ppjsdm_show_short_range_models", (DL_FUNC) &_ppjsdm_show_short_range_models, 0},
    {"_ppjsdm_show_medium_range_models", (DL_FUNC) &_ppjsdm_show_medium_range_models, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_ppjsdm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
