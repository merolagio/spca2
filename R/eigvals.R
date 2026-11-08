#' Computes All Eigenvalues of a Covariance or Correlation Matrix
#'
#' Compute all eigenvalues (without eigenvectors) of the covariance or
#' correlation matrix of the data, together with the proportions of
#' variance explained. Optionally produce a scree plot and a Wachter
#' qq-plot, and print the variance explained.
#'
#' @param M A data matrix, correlation matrix, or covariance matrix.
#' @param center_data A logical value (default \code{FALSE}). If \code{TRUE},
#'   center variables to zero mean. If \code{M} is detected as a data matrix and
#'   any column mean is nonzero, centering is performed automatically.
#' @param scale_data A logical value (default \code{FALSE}). If \code{TRUE},
#'   scale variables to unit variance, so that the eigenvalues are those of
#'   the correlation matrix.
#' @param screeplot A logical value (default \code{TRUE}). If \code{TRUE},
#'   produce a scree plot with \code{\link{scree_plot}}.
#' @param qq_plot A logical value (default \code{FALSE}). If \code{TRUE},
#'   produce a Wachter qq-plot with \code{\link{mp_qqplot}}.
#' @param n_obs An integer scalar or \code{NULL} (default \code{NULL}).
#'   Number of rows in the original data set. Required when
#'   \code{qq_plot = TRUE} and \code{M} is a covariance or correlation
#'   matrix. Ignored when \code{M} is a data matrix.
#' @param neigen_toplot An integer scalar or \code{NULL} (default
#'   \code{NULL}). Number of eigenvalues to show in the plots and in the
#'   printout. If \code{NULL}, all eigenvalues are shown.
#' @param cor A logical value (default \code{TRUE}). Passed to
#'   \code{\link{mp_qqplot}}: if \code{TRUE}, the reference variance of the
#'   Marchenko--Pastur quantiles is 1.
#' @param common_var A positive numeric scalar or \code{NULL} (default
#'   \code{NULL}). Common population variance assumed for the Wachter
#'   qq-plot. Required when \code{qq_plot = TRUE} and \code{cor = FALSE}.
#' @param print_vexp A logical value (default \code{FALSE}). If \code{TRUE},
#'   print the eigenvalues with the proportions of variance explained and
#'   cumulative variance explained.
#'
#' @details For a data matrix \eqn{X} (centred, and scaled if requested),
#'   the eigenvalues are those of the covariance matrix \eqn{X'X/(n - 1)},
#'   or of the correlation matrix if the data are scaled. They are computed
#'   from the smaller of \eqn{X'X} and \eqn{XX'}, so for a data matrix with
#'   fewer rows than columns only \eqn{n} eigenvalues are returned (the
#'   others are zero). For a covariance or correlation matrix the
#'   eigenvalues are those of \code{M}.
#'
#'   Only the eigenvalues are computed, which is cheaper than a full PCA.
#'   Use the result to choose the number of components, then compute them
#'   with \code{\link{pca}} or \code{\link{spca}}.
#'
#' @return An object of class \code{eigvals}: a list with elements
#' \describe{
#'   \item{eigenvalues}{All eigenvalues, in decreasing order.}
#'   \item{vexp}{Proportion of the total variance explained by each
#'     eigenvalue.}
#'   \item{cvexp}{Cumulative proportion of the total variance explained.}
#'   \item{tot_var}{Total variance (trace of the covariance or correlation
#'     matrix).}
#'   \item{n_vars}{Number of variables.}
#'   \item{n_obs}{Number of observations, or \code{NULL} if not available.}
#'   \item{call}{The matched call.}
#' }
#'
#' @examples
#' data(holzinger)
#' ho_eig = eigvals(holzinger, n_obs = 144, qq_plot = TRUE,
#'                  print_vexp = TRUE, neigen_toplot = 10)
#' # retain 4 components
#' ho_pca = pca(holzinger, n_comps = 4, screeplot = FALSE)
#' @family pca
#' @export
eigvals = function(M, center_data = FALSE, scale_data = FALSE,
                   screeplot = TRUE, qq_plot = FALSE, n_obs = NULL,
                   neigen_toplot = NULL, cor = TRUE, common_var = NULL,
                   print_vexp = FALSE) {

  # validation ==========

  vbool = validate_booleans(center_data = center_data,
                            scale_data = scale_data,
                            screeplot = screeplot,
                            qq_plot = qq_plot,
                            cor = cor,
                            print_vexp = print_vexp)
  if (!vbool) {
    stop("wrong input for boolean argument")
  }

  if (is.data.frame(M))
    M = as.matrix(M)
  if (!is.matrix(M))
    stop("M must be a matrix")
  if (any(is.na(M)))
    stop("The data matrix cannot contain missing values")

  n = nrow(M)
  p = ncol(M)

  is_datamatrix_M = FALSE
  if ((n != p) || !is_symmetric_fast(M)) {
    is_datamatrix_M = TRUE
    n_obs = n
  }

  if (!is.null(n_obs) &&
      (!is.numeric(n_obs) || length(n_obs) != 1 || is.na(n_obs) ||
       n_obs < 2 || n_obs != floor(n_obs)))
    stop("n_obs must be NULL or an integer greater than 1")

  if (is_datamatrix_M) {
    if (n < 2)
      stop("a data matrix needs at least two rows")
    # data matrix must be column mean centered
    if (any(abs(colMeans(M)) > 1e-6))
      center_data = TRUE
    if (center_data || scale_data)
      M = standardize_data(M, center_data, scale_data)
  }

  # computation ========

  eout = eigvalsC(M, data_matrix = is_datamatrix_M)

  eigenvalues = as.vector(eout$eigenvalues)
  tot_var = eout$totvar
  vexp = eigenvalues / tot_var

  out = list(eigenvalues = eigenvalues,
             vexp = vexp,
             cvexp = cumsum(vexp),
             tot_var = tot_var,
             n_vars = p,
             n_obs = n_obs)
  out$call = match.call()
  class(out) = c("eigvals", "list")

  # output ========

  n_eig = length(eigenvalues)
  if (is.null(neigen_toplot))
    neigen_toplot = n_eig
  if (!is.numeric(neigen_toplot) || length(neigen_toplot) != 1 ||
      is.na(neigen_toplot) || neigen_toplot < 1) {
    warning(paste("neigen_toplot must be a positive integer, setting it to",
                  n_eig))
    neigen_toplot = n_eig
  }
  if (neigen_toplot > n_eig) {
    warning(paste("neigen_toplot exceeds the number of eigenvalues;
                  setting it to", n_eig))
    neigen_toplot = n_eig
  }
  neigen_toplot = as.integer(neigen_toplot)

  if (print_vexp) {
    tab = data.frame(eigenvalue = eigenvalues,
                     vexp = vexp,
                     cvexp = cumsum(vexp))[seq_len(neigen_toplot), ,
                                           drop = FALSE]
    rownames(tab) = paste0("PC", seq_len(neigen_toplot))
    print(round(tab, 3))
  }

  if (qq_plot && is.null(n_obs)) {
    warning("qq-plot cannot be produced because n_obs is not available")
    qq_plot = FALSE
  }

  if (screeplot)
    scree_plot(out, n_plot = neigen_toplot, ylab = "eigenvalues")

  if (qq_plot)
    mp_qqplot(out, cor = cor, common_var = common_var,
              n_plot = neigen_toplot, n_fitline = NULL)

  invisible(out)
}

#' @rdname scree_plot
#' @exportS3Method
scree_plot.eigvals = function(
    pca_fit, n_plot = NULL, ylab = "eigenvalues", addtitle = TRUE,
    show_plot = TRUE, return_plot = FALSE) {

  if (!inherits(pca_fit, "eigvals"))
    stop("pca_fit must be an eigvals object", call. = FALSE)

  scree_plot_values(eigenvalues = pca_fit$eigenvalues, n_plot = n_plot,
                    ylab = ylab, addtitle = addtitle, show_plot = show_plot,
                    return_plot = return_plot, what = "pca_fit")
}

#' @rdname mp_qqplot
#' @exportS3Method
mp_qqplot.eigvals = function(
    pca_fit, n_vars = NULL, n_obs = NULL, gamma = NULL, cor = TRUE,
    common_var = NULL, n_plot = NULL, n_fitline = NULL, addtitle = TRUE,
    show_plot = TRUE, return_plot = FALSE) {

  if (!inherits(pca_fit, "eigvals"))
    stop("pca_fit must be an eigvals object", call. = FALSE)

  if (is.null(n_vars))
    n_vars = pca_fit$n_vars
  if (is.null(n_obs))
    n_obs = pca_fit$n_obs

  mp_qqplot_values(
    eigenvalues = pca_fit$eigenvalues, n_vars = n_vars, n_obs = n_obs,
    gamma = gamma, cor = cor, common_var = common_var, n_plot = n_plot,
    n_fitline = n_fitline, addtitle = addtitle, show_plot = show_plot,
    return_plot = return_plot, what = "pca_fit")
}

#' Print Eigenvalues
#'
#' Print the eigenvalues stored in an \code{eigvals} object.
#'
#' @param x An object of class \code{eigvals}, as returned by
#'   \code{\link{eigvals}()}.
#' @param digits An integer scalar (default \code{2}). Number of decimal
#'   places.
#' @param ... Further arguments. These are currently unused and trigger an
#'   error if supplied.
#'
#' @return \code{x}, invisibly.
#' @examples
#' data(holzinger)
#' ho_eig = eigvals(holzinger, screeplot = FALSE)
#' print(ho_eig)
#' @family pca
#' @exportS3Method
print.eigvals = function(x, digits = 2, ...) {
  dots = list(...)
  if (length(dots) > 0)
    stop("Unused arguments: ", paste(names(dots), collapse = ", "))
  if (!inherits(x, "eigvals"))
    stop("x must be an eigvals object", call. = FALSE)
  if (!is.numeric(digits) || length(digits) != 1 || is.na(digits) ||
      digits < 0)
    stop("digits must be a nonnegative integer", call. = FALSE)

  ev = round(x$eigenvalues, digits)
  names(ev) = paste0("PC", seq_along(ev))
  print(ev)
  invisible(x)
}

#' Summarize Eigenvalues
#'
#' Print a table with the eigenvalues, the percentage of variance explained
#' (Vexp) and the cumulative percentage of variance explained (Cvexp).
#'
#' @param object An object of class \code{eigvals}, as returned by
#'   \code{\link{eigvals}()}.
#' @param cols An integer vector of eigenvalue indices or \code{NULL}
#'   (default \code{NULL}: all eigenvalues). A single integer \code{k}
#'   selects the first \code{k} eigenvalues.
#' @param digits An integer scalar (default \code{1}). Number of decimal
#'   places of the percentages. Eigenvalues are shown with one decimal.
#' @param return_table A logical value (default \code{FALSE}). If
#'   \code{TRUE}, return the numeric table (proportions, not percentages).
#' @param print_table A logical value (default \code{TRUE}). If \code{TRUE},
#'   print the formatted table.
#' @param ... Further arguments. These are currently unused and trigger an
#'   error if supplied.
#'
#' @return If \code{return_table = TRUE}, a numeric matrix with rows
#'   \code{Eigval}, \code{Vexp} and \code{Cvexp}; otherwise \code{NULL}
#'   invisibly.
#' @examples
#' data(holzinger)
#' ho_eig = eigvals(holzinger, screeplot = FALSE)
#' summary(ho_eig, cols = 6)
#' @family pca
#' @exportS3Method
summary.eigvals = function(object, cols = NULL, digits = 1,
                           return_table = FALSE, print_table = TRUE, ...) {
  dots = list(...)
  if (length(dots) > 0)
    stop("Unused arguments: ", paste(names(dots), collapse = ", "))
  if (!inherits(object, "eigvals"))
    stop("object must be an eigvals object", call. = FALSE)
  validate_booleans(return_table = return_table, print_table = print_table)
  if (!is.numeric(digits) || length(digits) != 1 || is.na(digits) ||
      digits < 0)
    stop("digits must be a nonnegative integer", call. = FALSE)

  n_eig = length(object$eigenvalues)
  if (is.null(cols)) {
    cols = seq_len(n_eig)
  } else if (length(cols) == 1L) {
    cols = seq_len(cols)
  }
  if (!is.numeric(cols) || anyNA(cols) || any(cols < 1) ||
      any(cols > n_eig) || any(cols != floor(cols)))
    stop("cols must contain integers between 1 and the number of eigenvalues",
         call. = FALSE)

  out = rbind(Eigval = object$eigenvalues,
              Vexp = object$vexp,
              Cvexp = object$cvexp)[, cols, drop = FALSE]
  colnames(out) = paste0("PC", cols)

  if (print_table) {
    outp = rbind(
      Eigval = sprintf("%.1f", round(out["Eigval", ], 1)),
      Vexp = paste0(sprintf(paste0("%.", digits, "f"),
                            round(100 * out["Vexp", ], digits)), "%"),
      Cvexp = paste0(sprintf(paste0("%.", digits, "f"),
                             round(100 * out["Cvexp", ], digits)), "%"))
    colnames(outp) = colnames(out)
    print(outp, quote = FALSE, right = TRUE)
  }

  if (return_table)
    return(out)
  invisible()
}

#' Plot Eigenvalues
#'
#' Produce the scree plot and, optionally, the Wachter qq-plot of an
#' \code{eigvals} object, using \code{\link{scree_plot}} and
#' \code{\link{mp_qqplot}}.
#'
#' @param x An object of class \code{eigvals}, as returned by
#'   \code{\link{eigvals}()}.
#' @param screeplot A logical value (default \code{TRUE}). Produce the scree
#'   plot.
#' @param qq_plot A logical value (default \code{FALSE}). Produce the
#'   Wachter qq-plot.
#' @param n_plot An integer scalar or \code{NULL} (default \code{NULL}).
#'   Number of leading eigenvalues to plot. If \code{NULL}, all are plotted.
#' @param n_fitline An integer scalar or \code{NULL} (default \code{NULL}:
#'   no line). Passed to \code{\link{mp_qqplot}}.
#' @param n_obs,cor,common_var Passed to \code{\link{mp_qqplot}}.
#' @param ... Further arguments. These are currently unused and trigger an
#'   error if supplied.
#'
#' @return \code{x}, invisibly.
#' @examples
#' data(holzinger)
#' ho_eig = eigvals(holzinger, screeplot = FALSE)
#' plot(ho_eig, qq_plot = TRUE, n_fitline = -4)
#' @family pca
#' @exportS3Method
plot.eigvals = function(x, screeplot = TRUE, qq_plot = FALSE, n_plot = NULL,
                        n_fitline = NULL, n_obs = NULL, cor = TRUE,
                        common_var = NULL, ...) {
  dots = list(...)
  if (length(dots) > 0)
    stop("Unused arguments: ", paste(names(dots), collapse = ", "))
  if (!inherits(x, "eigvals"))
    stop("x must be an eigvals object", call. = FALSE)
  validate_booleans(screeplot = screeplot, qq_plot = qq_plot)

  if (screeplot)
    scree_plot(x, n_plot = n_plot)
  if (qq_plot)
    mp_qqplot(x, n_obs = n_obs, cor = cor, common_var = common_var,
              n_plot = n_plot, n_fitline = n_fitline)
  invisible(x)
}
