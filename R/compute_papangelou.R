#' Compute the Papangelou conditional intensity of the model.
#'
#' IMPORTANT: Check ?compute_papangelou.default for the documentation.
#'
#' @param ... Parameters to be forwarded to the relevant method. Currently, either
#' the parameters of the Gibbs point process, or a fit object obtained from running `gibbsm`.
#' @export
compute_papangelou <- function(...) {
  UseMethod("compute_papangelou")
}

#' Compute the Papangelou conditional intensity of the model with given parameters.
#'
#' @param configuration Configuration.
#' @param x Coordinates along the x-axis of the points at which to evaluate the Papangelou conditional intensity.
#' @param y Coordinates along the y-axis of the points at which to evaluate the Papangelou conditional intensity.
#' @param type Type of the point (as an integer >= 1 or string representing the type).
#' @param mark Mark of the point to add.
#' @param alpha Matrix of short-range interaction coefficients.
#' @param gamma Matrix of medium-range interaction coefficients.
#' @param beta0 Vector representing the intercept.
#' @param covariates List of covariates. These are converted to the `im` format by applying `as.im` to all elements in the list.
#' @param beta Fitted regression coefficients with respect to covariates.
#' @param short_range Symmetric matrix of short-range interaction radii. Can also be a list of matrices, each entry representing a different potential.
#' @param medium_range Symmetric matrix of medium-range interaction radii.
#' @param long_range Symmetric matrix of long-range interaction radii.
#' @param saturation Saturation parameter of the point process.
#' @param types Character vector, with entry i representing the name of type i.
#' @param model String representing the short-range model to use. The currently authorised models are obtained with a call to `show_short_range_models()`.
#' @param medium_range_model String representing the medium-range model to use. The currently authorised models are obtained with a call to `show_medium_range_models()`.
#' @param drop_type_from_configuration Should we remove the considered type(s) from the configuration?
#' @param nthreads Maximum number of threads for parallel computing.
#' @param ... Ignored.
#' @export
#' @method compute_papangelou default
compute_papangelou.default <- function(configuration,
                                       x,
                                       y,
                                       type = rep(1, length(x)),
                                       mark = rep(1.0, length(x)),
                                       alpha,
                                       gamma,
                                       beta0,
                                       covariates,
                                       beta,
                                       short_range,
                                       medium_range,
                                       long_range,
                                       saturation,
                                       types,
                                       model,
                                       medium_range_model,
                                       drop_type_from_configuration = FALSE,
                                       nthreads = 1,
                                       ...) {
  # If user did not supply types, by default they should be those of the configuration
  if(missing(types) & !missing(configuration)) {
    if(nlevels(as.Configuration(configuration)$types) > 0) {
      types <- levels(as.Configuration(configuration)$types)
    } else {
      types <- NULL
    }

    # If the user supplied type, then add them to the list of types
    if(all(is.character(type))) {
      types <- sort(union(types, unique(type)))
    }
  }

  # Check format of type and mark parameters
  check_type_mark <- function(obj, force_numeric) {
    if(length(obj) != length(x)) {
      if(length(obj) == 1) {
        rep(obj, length(x))
      } else {
        stop("Unknown format for type or mark.")
      }
    } else {
      if(all(is.numeric(obj)) | (all(is.character(obj)) & !force_numeric)) {
        obj
      } else {
        stop("Type or mark are not of the right class.")
      }
    }
  }

  type <- check_type_mark(type, force_numeric = FALSE)
  mark <- check_type_mark(mark, force_numeric = TRUE)

  # Relevant types contains types mentioned by the user by name, we want to make sure that these are included
  relevant_types <- if(is.character(type)) {
    sort(unique(type))
  } else {
    NULL
  }

  if(!missing(configuration)) {
    configuration <- as.Configuration(configuration)
    if(nlevels(configuration$types) > 0) {
      relevant_types <- sort(union(relevant_types, levels(configuration$types)))
    }
  } else {
    configuration <- Configuration()
  }

  # Relevant indices contains types mentioned by the user by index
  relevant_indices <- if(is.numeric(type)) {
    unique(type)
  } else {
    NULL
  }

  # Construct defaults for the rest of the parameters
  parameters <- model_parameters(alpha = alpha,
                                 gamma = gamma,
                                 beta0 = beta0,
                                 covariates = covariates,
                                 beta = beta,
                                 short_range = short_range,
                                 medium_range = medium_range,
                                 long_range = long_range,
                                 saturation = saturation,
                                 model = model,
                                 types = types,
                                 medium_range_model = medium_range_model,
                                 relevant_types = relevant_types,
                                 relevant_indices = relevant_indices)

  # If the user supplies a string as the type, they want to evaluate the intensity
  # at that type.
  if(is.character(type)) {
    type <- sapply(type, function(t) which(t == names(parameters$beta0))[1])
  } else {
    type <- sapply(type, function(i) which(parameters$full_types[i] == names(parameters$beta0))[1])
  }

  # Remove type from the configuration?
  if(drop_type_from_configuration) {
    configuration <- configuration[setdiff(configuration$types, names(parameters$beta0)[unique(type)])]
  }

  # Check type format to avoid segfaults in C++ code
  if(!(is.integer(type) & all(type >= 1 & type <= length(parameters$beta0)))) {
    stop("Type does not have the right format")
  }

  # Fix type levels (configuration types might have fewer levels than the ones implied by other parameters,
  # causing issues in the C++ code)
  configuration$types <- factor(as.character(configuration$types), levels = names(parameters$beta0))

  compute_papangelou_cpp(x = x,
                         y = y,
                         type = type,
                         mark = mark,
                         configuration = configuration,
                         model = parameters$model,
                         medium_range_model = parameters$medium_range_model,
                         alpha = parameters$alpha,
                         beta0 = parameters$beta0,
                         beta = parameters$beta,
                         gamma = parameters$gamma,
                         covariates = parameters$covariates,
                         short_range = parameters$short_range,
                         medium_range = parameters$medium_range,
                         long_range = parameters$long_range,
                         saturation = parameters$saturation,
                         nthreads = nthreads)
}

#' Compute the Papangelou conditional intensity of the model.
#'
#' @param ... Parameters to be forwarded to the relevant method. Currently, either
#' the parameters of the Gibbs point process, or a fit object obtained from running `gibbsm`.
#' @export
model_matrix_papangelou <- function(...) {
  UseMethod("model_matrix_papangelou")
}

#' Compute the Papangelou conditional intensity of the model with given parameters.
#'
#' @param configuration Configuration.
#' @param x Coordinates along the x-axis of the points at which to evaluate the Papangelou conditional intensity.
#' @param y Coordinates along the y-axis of the points at which to evaluate the Papangelou conditional intensity.
#' @param type Type of the point (as an integer >= 1).
#' @param model String representing the model to use. You can check the currently authorised models with a call to `show_models()`.
#' @param medium_range_model String representing the medium range model to use. You can check the currently authorised models with a call to `show_medium_range_models()`.
#' @param alpha Short range repulsion matrix.
#' @param beta0 A vector representing the log-intensities of the point processes.
#' @param beta Coefficients related to covariates.
#' @param gamma Medium range repulsion matrix.
#' @param covariates Covariates.
#' @param short_range Short range interaction radii.
#' @param medium_range Medium range interaction radii.
#' @param long_range Long range interaction radii.
#' @param saturation Saturation parameter.
#' @param types Types of the points. Default is a vector (type1, type2, ...) of same size as n.
#' @param drop_type_from_configuration Should we remove the considered type(s) from the configuration?
#' @importFrom stats na.omit
#' @param mark Mark of the point to add.
#' @param nthreads Maximum number of threads for parallel computing.
#' @param ... Ignored.
#' @export
#' @method model_matrix_papangelou default
model_matrix_papangelou.default <- function(configuration,
                                       x,
                                       y,
                                       type = rep(1, length(x)),
                                       mark = rep(1.0, length(x)),
                                       model,
                                       medium_range_model,
                                       alpha,
                                       beta0,
                                       beta,
                                       gamma,
                                       covariates,
                                       short_range,
                                       medium_range,
                                       long_range,
                                       saturation,
                                       types,
                                       drop_type_from_configuration = FALSE,
                                       nthreads = 1,
                                       ...) {
  # If user did not supply types, by default they should be those of the configuration
  if(missing(types) & !missing(configuration)) {
    types <- levels(as.Configuration(configuration)$types)
  }

  # Construct a default configuration if not supplied
  if(missing(configuration)) {
    configuration <- Configuration()
  }
  configuration <- as.Configuration(configuration)

  # Construct defaults for the rest of the parameters
  parameters <- model_parameters(alpha = alpha,
                                 gamma = gamma,
                                 beta0 = beta0,
                                 covariates = covariates,
                                 beta = beta,
                                 short_range = short_range,
                                 medium_range = medium_range,
                                 long_range = long_range,
                                 saturation = saturation,
                                 model = model,
                                 types = types,
                                 medium_range_model = medium_range_model)

  # There are some cases where configuration might only have a subset of all types,
  # e.g., if one wants to predict at a new location.
  if(nlevels(configuration$types) != length(parameters$types)) {
    # In this case, everything is ok, we can evaluate the parameters on the subset of types
    if(all(levels(configuration$types) %in% parameters$types)) {
      parameters$alpha <- lapply(parameters$alpha, function(a) a[levels(configuration$types), levels(configuration$types)])
      parameters$gamma <- parameters$gamma[levels(configuration$types), levels(configuration$types)]
      parameters$short_range <- lapply(parameters$short_range, function(s) s[levels(configuration$types), levels(configuration$types)])
      parameters$medium_range <- parameters$medium_range[levels(configuration$types), levels(configuration$types)]
      parameters$long_range <- parameters$long_range[levels(configuration$types), levels(configuration$types)]

      parameters$beta0 <- parameters$beta0[levels(configuration$types)]
      parameters$beta <- parameters$beta[levels(configuration$types), ]
    } else {
      stop(paste0("The types of the configuration are not a subset of those given by the parameters, configuration: ",
                  paste0(levels(configuration$types), collapse = ", "),
                  " and supplied types: ",
                  paste0(parameters$types, collapse = ", ")))
    }
  } else if(!all(levels(configuration$types) == parameters$types)) {
    stop(paste0("The types of the configuration do not correspond to those given by the parameters, configuration: ",
                paste0(levels(configuration$types), collapse = ", "),
                " and supplied types: ",
                paste0(parameters$types, collapse = ", ")))
  }

  # At this point, the parameters should exactly correspond to the types of the configuration
  if(!all(names(parameters$beta0) == levels(configuration$types))) {
    stop(paste0("Unexpected setting when computing the Papangelou intensity: the supplied parameters do not refer to the same types as the configuration, configuration: ",
                paste0(levels(configuration$types), collapse = ", "),
                " and supplied parameters: ",
                paste0(names(parameters$beta0), collapse = ", ")))
  }

  check_type_mark <- function(obj) {
    if(length(obj) != length(x)) {
      if(length(obj) == 1) {
        rep(obj, length(x))
      } else {
        stop("Unknown format for type or mark.")
      }
    } else {
      obj
    }
  }

  type <- check_type_mark(type)
  mark <- check_type_mark(mark)

  # If the user supplies a string as the type, they want to evaluate the intensity
  # at that type.
  if(is.character(type)) {
    type <- sapply(type, function(t) which(t == names(parameters$beta0))[1])
  }

  # Remove type from the configuration?
  if(drop_type_from_configuration) {
    configuration <- Configuration(x = configuration$x[!(configuration$types %in% names(parameters$beta0)[unique(type)])],
                                   y = configuration$y[!(configuration$types %in% names(parameters$beta0)[unique(type)])],
                                   types = configuration$types[!(configuration$types %in% names(parameters$beta0)[unique(type)])],
                                   marks = configuration$marks[!(configuration$types %in% names(parameters$beta0)[unique(type)])])
  }

  sr_disp <- get_sr_disp_papangelou_cpp(x = x,
                         y = y,
                         type = type,
                         mark = mark,
                         configuration = configuration,
                         model = parameters$model,
                         medium_range_model = parameters$medium_range_model,
                         alpha = parameters$alpha,
                         beta0 = parameters$beta0,
                         beta = parameters$beta,
                         gamma = parameters$gamma,
                         covariates = parameters$covariates,
                         short_range = parameters$short_range,
                         medium_range = parameters$medium_range,
                         long_range = parameters$long_range,
                         saturation = parameters$saturation,
                         nthreads = nthreads)

  sr_disp  <- lapply(1:length(sr_disp), \(i) {
    disp <- sr_disp[[i]]
    disp <- as.data.frame(do.call(rbind,disp))
    disp <- disp[-(1:length(configuration)),]
    names(disp) <- paste0("sr",i,"_",types)
    disp$x <- x
    disp$y <- y
    disp
  }) %>%
  purrr::reduce(dplyr::left_join)

  mr_disp <- get_mr_disp_papangelou_cpp(x = x,
                         y = y,
                         type = type,
                         mark = mark,
                         configuration = configuration,
                         model = parameters$model,
                         medium_range_model = parameters$medium_range_model,
                         alpha = parameters$alpha,
                         beta0 = parameters$beta0,
                         beta = parameters$beta,
                         gamma = parameters$gamma,
                         covariates = parameters$covariates,
                         short_range = parameters$short_range,
                         medium_range = parameters$medium_range,
                         long_range = parameters$long_range,
                         saturation = parameters$saturation,
                         nthreads = nthreads)

  mr_disp <- as.data.frame(do.call(rbind,mr_disp))
  mr_disp <- mr_disp[-(1:length(configuration)),]
  names(mr_disp) <- paste0("mr","_",types)
  mr_disp$x <- x
  mr_disp$y <- y

  disp <- dplyr::left_join(mr_disp,sr_disp,by=c("x","y"))

  disp$beta0 <- parameters$beta0[type]

  beta_covars <- get_beta_covariates_papangelou_cpp(x = x,
                        y = y,
                        type = type,
                        mark = mark,
                        configuration = configuration,
                        model = parameters$model,
                        medium_range_model = parameters$medium_range_model,
                        alpha = parameters$alpha,
                        beta0 = parameters$beta0,
                        beta = parameters$beta,
                        gamma = parameters$gamma,
                        covariates = parameters$covariates,
                        short_range = parameters$short_range,
                        medium_range = parameters$medium_range,
                        long_range = parameters$long_range,
                        saturation = parameters$saturation,
                        nthreads = nthreads)

  beta_covars <- as.data.frame(do.call(rbind,beta_covars))
  names(beta_covars) <- names(parameters$covariates)
  beta_covars$x <- x
  beta_covars$y <- y

  # beta_covars <- lapply(1:length(parameters$covariates), \(i) {
  #   obj <- parameters$covariates[[i]]
  #   df <- as.data.frame(obj)
  #   df$value <- df$value * parameters$beta[type,i]
  #   names(df)[3] <- paste0("beta_",names(parameters$covariates)[i])
  #   df
  # }) %>%
  #   purrr::reduce(dplyr::left_join,by=c("x","y"))

  model_matrix <- dplyr::left_join(disp,beta_covars,by=c("x","y"))

  # model_matrix %>%
  #   dplyr::select(dplyr::contains("mr") | dplyr::contains("sr")) %>%
  #   tibble::rownames_to_column() %>%
  #   tidyr::pivot_longer(-1, names_to = c(".value", "set"), names_sep = "[_]") -> d

  # sum_cols <- colnames(d)[-(1:2)]

  # d <- d %>%
  #   dplyr::mutate(overall = purrr::pmap_dbl(list(!!!rlang::parse_exprs(sum_cols)),sum)) %>%
  #   tidyr::pivot_wider(names_from = set,values_from = c(!!!rlang::parse_exprs(sum_cols),overall))

  # model_matrix <- model_matrix %>%
  #   tibble::rownames_to_column() %>%
  #   dplyr::left_join(d) %>%
  #   dplyr::select(-rowname) %>%
  #   dplyr::relocate(x,y)

  return(model_matrix)
}

#' Compute the Papangelou conditional intensity of the model from a fit object.
#'
#' @param fit Fit object obtained by running gibbsm.
#' @param type Type of the point at which to evaluate the Papangelou conditional intensity.
#' @param configuration Configuration of points at which to evaluate the Papangelou conditional intensity.
#' @param alpha Matrix of short-range interaction coefficients.
#' @param gamma Matrix of medium-range interaction coefficients.
#' @param beta0 Vector representing the intercept.
#' @param covariates List of covariates. These are converted to the `im` format by applying `as.im` to all elements in the list.
#' @param beta Fitted regression coefficients with respect to covariates.
#' @param short_range Symmetric matrix of short-range interaction radii. Can also be a list of matrices, each entry representing a different potential.
#' @param medium_range Symmetric matrix of medium-range interaction radii.
#' @param long_range Symmetric matrix of long-range interaction radii.
#' @param saturation Saturation parameter of the point process.
#' @param types Character vector, with entry i representing the name of type i.
#' @param model String representing the short-range model to use. The currently authorised models are obtained with a call to `show_short_range_models()`.
#' @param medium_range_model String representing the medium-range model to use. The currently authorised models are obtained with a call to `show_medium_range_models()`.
#' @param nthreads Maximum number of threads for parallel computing.
#' @param ... Forwarded to plot_papangelou
#' @export
#' @method compute_papangelou gibbsm
compute_papangelou.gibbsm <- function(fit,
                                      type = 1,
                                      configuration = fit$configuration_list[[1]],
                                      alpha = fit$coefficients$alpha,
                                      gamma = fit$coefficients$gamma,
                                      beta0 = fit$coefficients$beta0,
                                      covariates = fit$parameters$covariates,
                                      beta = fit$coefficients$beta,
                                      short_range = fit$coefficients$short_range,
                                      medium_range = fit$coefficients$medium_range,
                                      long_range = fit$coefficients$long_range,
                                      saturation = fit$parameters$saturation,
                                      types = fit$type_names,
                                      model = fit$parameters$model,
                                      medium_range_model = fit$parameters$medium_range_model,
                                      nthreads = fit$nthreads,
                                      ...) {
  compute_papangelou(type = type,
                     configuration = configuration,
                     alpha = alpha,
                     gamma = gamma,
                     beta0 = beta0,
                     covariates = covariates,
                     beta = beta,
                     short_range = short_range,
                     medium_range = medium_range,
                     long_range = long_range,
                     saturation = saturation,
                     types = types,
                     model = model,
                     medium_range_model = medium_range_model,
                     nthreads = nthreads,
                     ...)
}
