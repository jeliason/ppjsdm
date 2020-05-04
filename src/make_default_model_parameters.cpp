#include <Rcpp.h>
#include <Rinternals.h>

#include "utility/construct_if_missing.hpp"
#include "utility/get_number_types.hpp"
#include "utility/is_symmetric_matrix.hpp"
#include "utility/make_default_types.hpp"

// [[Rcpp::export]]
SEXP make_default_model_parameters(SEXP alpha,
                                   SEXP lambda,
                                   SEXP covariates,
                                   SEXP beta,
                                   SEXP gamma,
                                   SEXP short_range,
                                   SEXP medium_range,
                                   SEXP long_range,
                                   SEXP types) {
  const auto number_types(ppjsdm::get_number_types_and_check_conformance(alpha, gamma, lambda, short_range, medium_range, long_range, types));
  alpha = ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(alpha, 0., number_types);
  gamma = ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(gamma, 0., number_types);
  lambda = ppjsdm::construct_if_missing<Rcpp::NumericVector>(lambda, 1., number_types);
  if(!ppjsdm::is_symmetric_matrix(alpha) || !ppjsdm::is_symmetric_matrix(gamma)) {
    Rcpp::stop("Either alpha or gamma is not symmetric.");
  }

  const auto beta_nrows(number_types);
  const auto beta_ncols(Rcpp::as<Rcpp::List>(covariates).size());
  beta = ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(beta, 1., beta_nrows, beta_ncols);
  if(beta_ncols != 0 && (Rcpp::as<Rcpp::NumericMatrix>(beta).nrow() != beta_nrows
                            || Rcpp::as<Rcpp::NumericMatrix>(beta).ncol() != beta_ncols)) {
    Rcpp::stop("The parameter `beta` does not have the right dimensions.");
  }

  types = ppjsdm::make_types(types, number_types, lambda);
  const auto sh(ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(short_range, 0.1, number_types));
  const auto me(ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(medium_range, 0., number_types));
  const auto lo(ppjsdm::construct_if_missing<Rcpp::NumericMatrix>(long_range, 0., number_types));

  return Rcpp::List::create(Rcpp::Named("alpha") = alpha,
                            Rcpp::Named("lambda") = lambda,
                            Rcpp::Named("covariates") = covariates,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("gamma") =gamma,
                            Rcpp::Named("short_range") = sh,
                            Rcpp::Named("medium_range") = me,
                            Rcpp::Named("long_range") = lo,
                            Rcpp::Named("types") = types);
}
