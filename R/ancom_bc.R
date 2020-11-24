#' @title Differential abundance (DA) analysis for
#' microbial absolute abundance data.
#'
#' @aliases ancom
#'
#' @description Determine taxa whose absolute abundances, per unit volume, of
#' the ecosystem (e.g. gut) are significantly different with changes in the
#' covariate of interest (e.g. the group effect). The current version of
#' \code{ancombc} function implements Analysis of Compositions of Microbiomes
#' with Bias Correction (ANCOM-BC) in cross-sectional data while allowing
#' the adjustment of covariates.
#'
#' @details The definition of structural zero can be found at
#' \href{https://doi.org/10.3389/fmicb.2017.02114}{ANCOM-II}.
#' Setting \code{neg_lb = TRUE} indicates that you are using both criteria
#' stated in section 3.2 of
#' \href{https://doi.org/10.3389/fmicb.2017.02114}{ANCOM-II}
#' to detect structural zeros; otherwise, the algorithm will only use the
#' equation 1 in section 3.2 for declaring structural zeros. Generally, it is
#' recommended to set \code{neg_lb = TRUE} when the sample size per group is
#' relatively large (e.g. > 30).
#'
#' @param phyloseq a phyloseq-class object, which consists of a feature table
#' (microbial observed abundance table), a sample metadata, a taxonomy table
#' (optional), and a phylogenetic tree (optional). The row names of the
#' metadata must match the sample names of the feature table, and the row names
#' of the taxonomy table must match the taxon (feature) names of the feature
#' table. See \code{\link[phyloseq]{phyloseq}} for more details.
#' @param formula the character string expresses how the microbial absolute
#' abundances for each taxon depend on the variables in metadata.
#' @param p_adj_method method to adjust p-values by. Default is "holm".
#' Options include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#' "fdr", "none". See \code{\link[stats]{p.adjust}} for more details.
#' @param zero_cut a numerical fraction between 0 and 1. Taxa with proportion of
#' zeroes greater than \code{zero_cut} will be excluded in the analysis. Default
#' is 0.90.
#' @param lib_cut a numerical threshold for filtering samples based on library
#' sizes. Samples with library sizes less than \code{lib_cut} will be
#' excluded in the analysis.
#' @param group the name of the group variable in metadata. Specifying
#' \code{group} is required for detecting structural zeros and
#' performing global test.
#' @param struc_zero whether to detect structural zeros. Default is FALSE.
#' @param neg_lb whether to classify a taxon as a structural zero in the
#' corresponding study group using its asymptotic lower bound.
#' Default is FALSE.
#' @param tol the iteration convergence tolerance for the E-M algorithm.
#' Default is 1e-05.
#' @param max_iter the maximum number of iterations for the E-M algorithm.
#' Default is 100.
#' @param conserve whether to use a conservative variance estimate of
#' the test statistic. It is recommended if the sample size is small and/or
#' the number of differentially abundant taxa is believed to be large.
#' Default is FALSE.
#' @param alpha level of significance. Default is 0.05.
#' @param global whether to perform global test. Default is FALSE.
#'
#' @return a \code{list} with components:
#'         \itemize{
#'         \item{ \code{feature_table}, a \code{data.frame} of pre-processed
#'         (based on \code{zero_cut} and \code{lib_cut}) microbial observed
#'         abundance table. }
#'         \item{ \code{zero_ind}, a logical \code{matrix} with TRUE indicating
#'         the taxon is identified as a structural zero for the specified
#'         \code{group} variable.}
#'         \item{ \code{samp_frac}, a numeric vector of estimated sampling
#'         fractions in log scale (natural log). }
#'         \item{ \code{resid}, a \code{matrix} of residuals from the ANCOM-BC
#'         log-linear (natural log) model.
#'         Rows are taxa and columns are samples.}
#'         \item{ \code{delta_em}, estimated bias terms through E-M algorithm. }
#'         \item{ \code{delta_wls}, estimated bias terms through weighted
#'         least squares (WLS) algorithm.}
#'         \item{ \code{res},  a \code{list} containing ANCOM-BC primary result,
#'         which consists of:}
#'         \itemize{
#'         \item{ \code{beta}, a \code{data.frame} of coefficients obtained
#'         from the ANCOM-BC log-linear (natural log) model. }
#'         \item{ \code{se}, a \code{data.frame} of standard errors (SEs) of
#'         \code{beta}. }
#'         \item{ \code{W}, a \code{data.frame} of test statistics.
#'         \code{W = beta/se}. }
#'         \item{ \code{p_val}, a \code{data.frame} of p-values. P-values are
#'         obtained from two-sided Z-test using the test statistic \code{W}. }
#'         \item{ \code{q_val}, a \code{data.frame} of adjusted p-values.
#'         Adjusted p-values are obtained by applying \code{p_adj_method}
#'         to \code{p_val}.}
#'         \item{ \code{diff_abn}, a logical \code{data.frame}. TRUE if the
#'         taxon has \code{q_val} less than \code{alpha}.}
#'         }
#'         \item{ \code{res_global},  a \code{data.frame} containing ANCOM-BC
#'         global test result for the variable specified in \code{group},
#'         each column is:}
#'         \itemize{
#'         \item{ \code{W}, test statistics.}
#'         \item{ \code{p_val}, p-values, which are obtained from two-sided
#'         Chi-square test using \code{W}.}
#'         \item{ \code{q_val}, adjusted p-values. Adjusted p-values are
#'         obtained by applying \code{p_adj_method} to \code{p_val}.}
#'         \item{ \code{diff_abn}, A logical vector. TRUE if the taxon has
#'         \code{q_val} less than \code{alpha}.}
#'         }
#'         }
#'
#' @examples
#' #================Build a Phyloseq-Class Object from Scratch==================
#' library(phyloseq)
#'
#' otu_mat = matrix(sample(1:100, 100, replace = TRUE), nrow = 10, ncol = 10)
#' rownames(otu_mat) = paste0("taxon", 1:nrow(otu_mat))
#' colnames(otu_mat) = paste0("sample", 1:ncol(otu_mat))
#'
#'
#' meta = data.frame(group = sample(LETTERS[1:4], size = 10, replace = TRUE),
#'                   row.names = paste0("sample", 1:ncol(otu_mat)),
#'                   stringsAsFactors = FALSE)
#'
#' tax_mat = matrix(sample(letters, 70, replace = TRUE),
#'                  nrow = nrow(otu_mat), ncol = 7)
#' rownames(tax_mat) = rownames(otu_mat)
#' colnames(tax_mat) = c("Kingdom", "Phylum", "Class", "Order",
#'                       "Family", "Genus", "Species")
#'
#' OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
#' META = sample_data(meta)
#' TAX = tax_table(tax_mat)
#' physeq = phyloseq(OTU, META, TAX)
#'
#' #========================Run ANCOMBC Using a Real Data=======================
#'
#' library(phyloseq)
#' library(tidyverse)
#' data(GlobalPatterns)
#'
#' # Aggregate to phylum level
#' phylum_data = tax_glom(GlobalPatterns, "Phylum")
#' # The taxonomy table
#' tax_mat = as(tax_table(phylum_data), "matrix")
#'
#' # Run ancombc function
#' out = ancombc(phyloseq = phylum_data, formula = "SampleType",
#'               p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
#'               group = "SampleType", struc_zero = TRUE, neg_lb = FALSE,
#'               tol = 1e-5, max_iter = 100, conserve = TRUE,
#'               alpha = 0.05, global = TRUE)
#'
#' res = out$res
#' res_global = out$res_global
#'
#' @author Huang Lin
#'
#
#'
#' @import stats
#' @import phyloseq
#' @importFrom MASS ginv
#' @importFrom nloptr neldermead
#' @importFrom Rdpack reprompt
#'
#' @export
ancombc = function(phyloseq, formula, p_adj_method = "holm", zero_cut = 0.90,
                   lib_cut, group = NULL, struc_zero = FALSE, neg_lb = FALSE,
                   tol = 1e-05, max_iter = 100, conserve = FALSE, alpha = 0.05,
                   global = FALSE){
  # 1. Data pre-processing
  fiuo_prep = data_prep(phyloseq, group, zero_cut, lib_cut, global = global)
  feature_table = fiuo_prep$feature_table
  meta_data = fiuo_prep$meta_data
  global = fiuo_prep$global
  # samp_id = colnames(feature_table)
  # taxa_id = rownames(feature_table)
  # n_samp = ncol(feature_table)
  n_taxa = nrow(feature_table)
  # Add pseudocount (1) and take logarithm.
  y = log(feature_table + 1)
  x = get_x(formula, meta_data)
  covariates = colnames(x)
  n_covariates = length(covariates)

  # 2. Identify taxa with structural zeros
  if (struc_zero) {
    if (is.null(group)) {
      stop("Please specify the group variable for detecting structural zeros.")
    }
    zero_ind = get_struc_zero(feature_table, meta_data, group, neg_lb)
  }else{ zero_ind = NULL }

  # 3. Estimation of parameters
  fiuo_para = para_est(y, meta_data, formula, tol, max_iter)
  beta = fiuo_para$beta
  d = fiuo_para$d
  e = fiuo_para$e
  var_cov_hat = fiuo_para$var_cov_hat
  var_hat = fiuo_para$var_hat

  # 4. Estimation of the between-sample bias
  fiuo_bias = bias_est(beta, var_hat, tol, max_iter, n_taxa)
  delta_em = fiuo_bias$delta_em
  delta_wls = fiuo_bias$delta_wls
  var_delta = fiuo_bias$var_delta

  # 5. Coefficients, standard error, and sampling fractions
  fiuo_fit = fit_summary(y, x, beta, var_hat, delta_em, var_delta, conserve)
  beta_hat = fiuo_fit$beta_hat
  se_hat = fiuo_fit$se_hat
  d_hat = fiuo_fit$d_hat

  # 6. Primary results
  W = beta_hat/se_hat
  p = 2 * pnorm(abs(W), mean = 0, sd = 1, lower.tail = FALSE)
  q = apply(p, 2, function(x) p.adjust(x, method = p_adj_method))
  diff_abn = q < alpha & !is.na(q)
  res = list(beta = data.frame(beta_hat, check.names = FALSE),
             se = data.frame(se_hat, check.names = FALSE),
             W = data.frame(W, check.names = FALSE),
             p_val = data.frame(p, check.names = FALSE),
             q_val = data.frame(q, check.names = FALSE),
             diff_abn = data.frame(diff_abn, check.names = FALSE))

  # 7. Global test results
  if (global) {
    res_global = global_test(y, x, group, beta_hat, var_cov_hat,
                             p_adj_method, alpha)
  } else { res_global = NULL }

  # 8. Combine the information of structural zeros
  fiuo_out = res_combine_zero(x, group, struc_zero, zero_ind, alpha,
                              global, res, res_global)
  res = fiuo_out$res
  res_global = fiuo_out$res_global

  # 9. Outputs
  out = list(feature_table = feature_table, zero_ind = zero_ind,
             samp_frac = d_hat, resid = e,
             delta_em = delta_em, delta_wls = delta_wls,
             res = res, res_global = res_global)
  return(out)
}



# E-M algorithm for estimating the bias term
bias_est = function(beta, var_hat, tol, max_iter, n_taxa) {
  delta_em = rep(NA, ncol(beta) - 1)
  delta_wls = rep(NA, ncol(beta) - 1)
  var_delta = rep(NA, ncol(beta) - 1)
  for (i in seq_along(delta_em)) {
    # Ignore the intercept
    Delta = beta[, i + 1]
    Delta = Delta[!is.na(Delta)]
    nu0 = var_hat[, i + 1]
    nu0 = nu0[!is.na(nu0)]

    # Initials
    pi0_0 = 0.75
    pi1_0 = 0.125
    pi2_0 = 0.125
    delta_0 = mean(Delta[Delta >= quantile(Delta, 0.25, na.rm = TRUE)&
                           Delta <= quantile(Delta, 0.75, na.rm = TRUE)],
                   na.rm = TRUE)
    l1_0 = mean(Delta[Delta < quantile(Delta, 0.125, na.rm = TRUE)],
                na.rm = TRUE)
    l2_0 = mean(Delta[Delta > quantile(Delta, 0.875, na.rm = TRUE)],
                na.rm = TRUE)
    kappa1_0 = var(Delta[Delta < quantile(Delta, 0.125, na.rm = TRUE)],
                   na.rm = TRUE)
    if(is.na(kappa1_0)|kappa1_0 == 0) kappa1_0 = 1
    kappa2_0 = var(Delta[Delta > quantile(Delta, 0.875, na.rm = TRUE)],
                   na.rm = TRUE)
    if(is.na(kappa2_0)|kappa2_0 == 0) kappa2_0 = 1

    # Apply E-M algorithm
    fiuo_em = em_iter(Delta, nu0, pi0_0, pi1_0, pi2_0, delta_0,
                      l1_0, l2_0, kappa1_0, kappa2_0, tol, max_iter)

    # The EM estimator of bias
    delta_em[i] = fiuo_em$delta

    # The WLS estimator of bias
    pi1 = fiuo_em$pi1
    pi2 = fiuo_em$pi2
    l1 = fiuo_em$l1
    l2 = fiuo_em$l2
    kappa1 = fiuo_em$kappa1
    kappa2 = fiuo_em$kappa2
    # Cluster 0
    C0 = which(Delta >= quantile(Delta, pi1, na.rm = TRUE) &
                 Delta < quantile(Delta, 1 - pi2, na.rm = TRUE))
    # Cluster 1
    C1 = which(Delta < quantile(Delta, pi1, na.rm = TRUE))
    # Cluster 2
    C2 = which(Delta >= quantile(Delta, 1 - pi2, na.rm = TRUE))
    # Numerator of the WLS estimator
    nu = nu0
    nu[C1] = nu[C1] + kappa1
    nu[C2] = nu[C2] + kappa2
    wls_deno = sum(1 / nu)
    # Denominator of the WLS estimator
    wls_nume = 1 / nu
    wls_nume[C0] = (wls_nume * Delta)[C0]
    wls_nume[C1] = (wls_nume * (Delta - l1))[C1]
    wls_nume[C2] = (wls_nume * (Delta - l2))[C2]
    wls_nume = sum(wls_nume)

    delta_wls[i] = wls_nume / wls_deno

    # Estimate the variance of bias
    var_delta[i] = 1 / wls_deno
    if (is.na(var_delta[i])) var_delta[i] = 0
  }

  fiuo_bias = list(delta_em = delta_em, delta_wls = delta_wls,
                   var_delta = var_delta)
}

# Data pre-processing
data_prep = function(phyloseq, group, zero_cut, lib_cut, global = global) {
  feature_table = as(otu_table(phyloseq), "matrix")
  feature_table = data.frame(feature_table, check.names = FALSE)
  meta_data = as(sample_data(phyloseq), "data.frame")
  # Drop unused levels
  meta_data[] = lapply(meta_data, function(x)
    if(is.factor(x)) factor(x) else x)
  # Check the group variable
  if (is.null(group)) {
    if (global) {
      stop("Please specify the group variable for the global test.")
    }
  } else {
    n_level = length(unique(meta_data[, group]))
    if (n_level < 2) {
      stop("The group variable should have >= 2 categories.")
    } else if (n_level < 3) {
      global = FALSE
      warning("The multi-group comparison will be deactivated as the group variable has < 3 categories.")
    }
  }

  # Discard taxa with zeros >= zero_cut
  zero_prop = apply(feature_table, 1, function(x)
    sum(x == 0, na.rm = TRUE)/length(x[!is.na(x)]))
  tax_del = which(zero_prop >= zero_cut)
  if (length(tax_del) > 0) {
    feature_table = feature_table[- tax_del, ]
  }

  # Discard samples with library size < lib_cut
  lib_size = colSums(feature_table, na.rm = TRUE)
  if(any(lib_size < lib_cut)){
    subj_del = which(lib_size < lib_cut)
    feature_table = feature_table[, - subj_del]
    meta_data = meta_data[- subj_del, ]
  }
  fiuo_prep = list(feature_table = feature_table,
                   meta_data = meta_data,
                   global = global)
  return(fiuo_prep)
}

em_iter = function(Delta, nu0, pi0_0, pi1_0, pi2_0, delta_0,
                   l1_0, l2_0, kappa1_0, kappa2_0, tol, max_iter) {
  # Store all paras in vectors/matrices
  pi0_vec = pi0_0
  pi1_vec = pi1_0
  pi2_vec = pi2_0
  delta_vec = delta_0
  l1_vec = l1_0
  l2_vec = l2_0
  kappa1_vec = kappa1_0
  kappa2_vec = kappa2_0
  n_taxa = length(Delta)

  # E-M iteration
  iterNum = 0
  epsilon = 100
  while (epsilon > tol & iterNum < max_iter) {
    # Current value of paras
    pi0 = pi0_vec[length(pi0_vec)]
    pi1 = pi1_vec[length(pi1_vec)]
    pi2 = pi2_vec[length(pi2_vec)]
    delta = delta_vec[length(delta_vec)]
    l1 = l1_vec[length(l1_vec)]
    l2 = l2_vec[length(l2_vec)]
    kappa1 = kappa1_vec[length(kappa1_vec)]
    kappa2 = kappa2_vec[length(kappa2_vec)]

    # E-step
    pdf0 = vapply(seq(n_taxa), function(i)
      dnorm(Delta[i], delta, sqrt(nu0[i])), FUN.VALUE = double(1))
    pdf1 = vapply(seq(n_taxa), function(i)
      dnorm(Delta[i], delta + l1, sqrt(nu0[i] + kappa1)),
      FUN.VALUE = double(1))
    pdf2 = vapply(seq(n_taxa), function(i)
      dnorm(Delta[i], delta + l2, sqrt(nu0[i] + kappa2)),
      FUN.VALUE = double(1))
    r0i = pi0*pdf0/(pi0*pdf0 + pi1*pdf1 + pi2*pdf2)
    r0i[is.na(r0i)] = 0
    r1i = pi1*pdf1/(pi0*pdf0 + pi1*pdf1 + pi2*pdf2)
    r1i[is.na(r1i)] = 0
    r2i = pi2*pdf2/(pi0*pdf0 + pi1*pdf1 + pi2*pdf2)
    r2i[is.na(r2i)] = 0

    # M-step
    pi0_new = mean(r0i, na.rm = TRUE)
    pi1_new = mean(r1i, na.rm = TRUE)
    pi2_new = mean(r2i, na.rm = TRUE)
    delta_new = sum(r0i*Delta/nu0 + r1i*(Delta-l1)/(nu0+kappa1) +
                      r2i*(Delta-l2)/(nu0+kappa2), na.rm = TRUE)/
      sum(r0i/nu0 + r1i/(nu0+kappa1) + r2i/(nu0+kappa2), na.rm = TRUE)
    l1_new = min(sum(r1i*(Delta-delta)/(nu0+kappa1), na.rm = TRUE)/
                   sum(r1i/(nu0+kappa1), na.rm = TRUE), 0)
    l2_new = max(sum(r2i*(Delta-delta)/(nu0+kappa2), na.rm = TRUE)/
                   sum(r2i/(nu0+kappa2), na.rm = TRUE), 0)

    # Nelder-Mead simplex algorithm for kappa1 and kappa2
    obj_kappa1 = function(x){
      log_pdf = log(vapply(seq(n_taxa), function(i)
        dnorm(Delta[i], delta+l1, sqrt(nu0[i]+x)),
        FUN.VALUE = double(1)))
      log_pdf[is.infinite(log_pdf)] = 0
      -sum(r1i*log_pdf, na.rm = TRUE)
    }
    kappa1_new = nloptr::neldermead(x0 = kappa1,
                                    fn = obj_kappa1, lower = 0)$par

    obj_kappa2 = function(x){
      log_pdf = log(vapply(seq(n_taxa), function(i)
        dnorm(Delta[i], delta+l2, sqrt(nu0[i]+x)),
        FUN.VALUE = double(1)))
      log_pdf[is.infinite(log_pdf)] = 0
      -sum(r2i*log_pdf, na.rm = TRUE)
    }
    kappa2_new = nloptr::neldermead(x0 = kappa2,
                                    fn = obj_kappa2, lower = 0)$par

    # Merge to the paras vectors/matrices
    pi0_vec = c(pi0_vec, pi0_new)
    pi1_vec = c(pi1_vec, pi1_new)
    pi2_vec = c(pi2_vec, pi2_new)
    delta_vec = c(delta_vec, delta_new)
    l1_vec = c(l1_vec, l1_new)
    l2_vec = c(l2_vec, l2_new)
    kappa1_vec = c(kappa1_vec, kappa1_new)
    kappa2_vec = c(kappa2_vec, kappa2_new)

    # Calculate the new epsilon
    epsilon = sqrt((pi0_new-pi0)^2 + (pi1_new-pi1)^2 + (pi2_new-pi2)^2 +
                     (delta_new-delta)^2 + (l1_new-l1)^2 + (l2_new-l2)^2 +
                     (kappa1_new-kappa1)^2 + (kappa2_new-kappa2)^2)
    iterNum = iterNum + 1
  }
  fiuo_em = list(pi0 = pi0_new, pi1 = pi1_new, pi2 = pi2_new,
                 delta = delta_new, l1 = l1_new, l2 = l2_new,
                 kappa1 = kappa1_new, kappa2 = kappa2_new)
  return(fiuo_em)
}
fit_summary = function(y, x, beta, var_hat, delta_em, var_delta, conserve) {
  n_taxa = nrow(y)

  beta_hat = beta
  beta_hat[, -1] = t(t(beta_hat[, -1]) - delta_em)

  if (conserve) {
    # Account for the variance of delta_hat
    se_hat = sqrt(sweep(var_hat, 2, c(0, var_delta), "+") +
                    2 * sqrt(sweep(var_hat, 2, c(0, var_delta), "*")))
  }else{ se_hat = sqrt(var_hat) }

  d_hat = vector()
  for (i in seq_len(n_taxa)) {
    d_hat_i = y[i, ] - x %*% beta_hat[i, ]
    d_hat = rbind(d_hat, d_hat_i)
  }
  d_hat = colMeans(d_hat, na.rm = TRUE)

  fiuo_fit = list(beta_hat = beta_hat, se_hat = se_hat, d_hat = d_hat)
  return(fiuo_fit)
}
# Identify structural zeros
get_struc_zero = function(feature_table, meta_data, group, neg_lb) {
  group_data = factor(meta_data[, group])
  present_table = as.matrix(feature_table)
  present_table[is.na(present_table)] = 0
  present_table[present_table != 0] = 1
  n_taxa = nrow(feature_table)
  n_group = nlevels(group_data)

  p_hat = matrix(NA, nrow = n_taxa, ncol = n_group)
  rownames(p_hat) = rownames(feature_table)
  colnames(p_hat) = levels(group_data)
  for (i in seq_len(n_taxa)) {
    p_hat[i, ] = tapply(present_table[i, ], group_data,
                        function(x) mean(x, na.rm = TRUE))
  }

  samp_size = matrix(NA, nrow = n_taxa, ncol = n_group)
  rownames(samp_size) = rownames(feature_table)
  colnames(samp_size) = levels(group_data)
  for (i in seq_len(n_taxa)) {
    samp_size[i, ] = tapply(as.matrix(feature_table)[i, ], group_data,
                            function(x) length(x[!is.na(x)]))
  }

  p_hat_lo = p_hat - 1.96 * sqrt(p_hat * (1 - p_hat)/samp_size)

  zero_ind = (p_hat == 0)
  # Do we classify a taxon as a structural zero by its negative lower bound?
  if (neg_lb) zero_ind[p_hat_lo <= 0] = TRUE

  colnames(zero_ind) = paste0("structural_zero (", group,
                              " = ", colnames(zero_ind), ")")
  return(zero_ind)
}
get_x = function(formula, meta_data) {
  opt = options(na.action = "na.pass") # Keep NA's in rows of x
  on.exit(options(opt)) # Switch it back
  x = model.matrix(formula(paste0("~", formula)), data = meta_data)
  return(x)
}
# Global test
global_test = function(y, x, group, beta_hat, var_cov_hat, p_adj_method, alpha){
  taxa_id = rownames(y)
  n_taxa = nrow(y)
  covariates = colnames(x)

  res_global = data.frame(matrix(NA, nrow = n_taxa, ncol = 4))
  rownames(res_global) = taxa_id
  colnames(res_global) = c("W", "p_val", "q_val", "diff_abn")

  group_ind = grepl(group, covariates)
  # Loop over the parameter of interest
  beta_hat_sub = beta_hat[, group_ind]
  var_cov_hat_sub = lapply(var_cov_hat, function(x)
    x = x[group_ind, group_ind])

  for (i in seq_len(n_taxa)) {
    # Loop over taxa
    beta_hat_sub_i = beta_hat_sub[i, ]
    var_cov_hat_sub_i = var_cov_hat_sub[[i]]
    A = diag(x = 1, nrow = length(beta_hat_sub_i))
    W = t(A %*% beta_hat_sub_i) %*%
      MASS::ginv(A %*% var_cov_hat_sub_i %*% t(A)) %*%
      (A %*% beta_hat_sub_i)
    p = 2 * min(pchisq(W, df = length(beta_hat_sub_i), lower.tail = TRUE),
                pchisq(W, df = length(beta_hat_sub_i), lower.tail = FALSE))
    res_global[i, "W"] = W
    res_global[i, "p_val"] = p
  }
  # Model summary
  q_global = p.adjust(res_global[, "p_val"], method = p_adj_method)
  q_global[is.na(q_global)] = 1
  diff_global = q_global < alpha & !is.na(q_global)

  res_global$q_val = q_global
  res_global$diff_abn = diff_global
  return(res_global)
}
para_est = function(y, meta_data, formula, tol, max_iter) {
  x = get_x(formula, meta_data)
  taxa_id = rownames(y)
  n_taxa = nrow(y)
  samp_id = colnames(y)
  n_samp = ncol(y)
  covariates = colnames(x)

  # Sampling fractions
  d = rep(0, n_samp)
  tformula = formula(paste0("y ~ ", formula))
  fits = lapply(seq_len(n_taxa), function(i) {
    df = data.frame(y = unlist(y[i, ]) - d, meta_data)
    return(lm(tformula, data = df))
  })
  # Regression coefficients
  beta = lapply(fits, function(i) {
    beta_i = rep(NA, length(covariates)) # prevent errors of missing values
    coef_i = coef(i)
    beta_i[match(names(coef_i), covariates)] = coef_i
    return(beta_i)
  })
  beta = Reduce('rbind', beta)

  # Iterative least square
  iterNum = 0
  epsilon = 100
  while (epsilon > tol & iterNum < max_iter) {
    # Updating beta
    fits = lapply(seq_len(n_taxa), function(i) {
      df = data.frame(y = unlist(y[i, ]) - d, meta_data)
      return(lm(tformula, data = df))
    })
    beta_new = lapply(fits, function(i) {
      beta_i = rep(NA, length(covariates))
      coef_i = coef(i)
      beta_i[match(names(coef_i), covariates)] = coef_i
      return(beta_i)
    })
    beta_new = Reduce('rbind', beta_new)

    # Updating d
    y_hat = lapply(fits, function(i) {
      y_hat_i = rep(NA, n_samp)
      fit_i = fitted(i)
      y_hat_i[match(names(fit_i), samp_id)] = fit_i
      return(y_hat_i)

    })
    y_hat = Reduce('rbind', y_hat)
    d_new = colMeans(y - y_hat, na.rm = TRUE)

    # Iteration
    epsilon = sqrt(sum((beta_new - beta)^2, na.rm = TRUE) +
                     sum((d_new - d)^2, na.rm = TRUE))
    iterNum = iterNum + 1
    beta = beta_new
    d = d_new
  }

  # Regression residuals
  y_hat = lapply(fits, function(i) {
    y_hat_i = rep(NA, n_samp)
    fit_i = fitted(i)
    y_hat_i[match(names(fit_i), samp_id)] = fit_i
    return(y_hat_i)

  })
  y_hat = Reduce('rbind', y_hat)
  e = t(t(y - y_hat) - d)

  # Variance-covariance matrices of coefficients
  fiuo_var_cov = var_cov_est(x, e, n_taxa)
  var_cov_hat = fiuo_var_cov$var_cov_hat
  var_hat = fiuo_var_cov$var_hat

  colnames(beta) = covariates
  rownames(beta) = taxa_id
  names(d) = samp_id
  names(var_cov_hat) = taxa_id
  colnames(var_hat) = covariates
  rownames(var_hat) = taxa_id

  fiuo_para = list(beta = beta, d = d, e = e,
                   var_cov_hat = var_cov_hat, var_hat = var_hat)
  return(fiuo_para)
}
res_combine_zero = function(x, group, struc_zero, zero_ind, alpha,
                            global, res, res_global) {
  covariates = colnames(x)

  # Set p/q-values of structural zeros to be 0s.
  if (struc_zero) {
    group_ind = grepl(group, covariates)
    zero_mask = 1 - apply(zero_ind, 1, function(x) any(x == 1))
    res$p_val[, group_ind] = res$p_val[, group_ind] * zero_mask
    res$q_val[, group_ind] = res$q_val[, group_ind] * zero_mask
    res$diff_abn = res$q_val < alpha & !is.na(res$q_val)

    # Global test
    if (global) {
      res_global[, "p_val"] = res_global[, "p_val"] * zero_mask
      res_global[, "q_val"] = res_global[, "q_val"] * zero_mask
      res_global[, "diff_abn"] = res_global[, "q_val"] < alpha &
        !is.na(res_global[, "q_val"])
    }
  }
  fiuo_out = list(res = res, res_global = res_global)
  return(fiuo_out)
}
var_cov_est = function(x, e, n_taxa) {
  covariates = colnames(x)
  n_covariates = length(covariates)
  n_samp = nrow(x)
  XTX_inv = MASS::ginv(t(x[complete.cases(x), ]) %*% x[complete.cases(x), ])
  var_cov_hat = vector(mode = "list", length = n_taxa) # Covariances
  var_hat = matrix(NA, nrow = n_taxa, ncol = n_covariates) # Variances
  for (i in seq_len(n_taxa)) {
    sigma2_xxT = matrix(0, ncol = n_covariates, nrow = n_covariates)
    for (j in seq_len(n_samp)) {
      sigma2_xxT_j = e[i, j]^2 * x[j, ] %*% t(x[j, ])
      sigma2_xxT_j[is.na(sigma2_xxT_j)] = 0
      sigma2_xxT = sigma2_xxT + sigma2_xxT_j
    }
    var_cov_hat[[i]] = XTX_inv %*% sigma2_xxT %*% XTX_inv
    rownames(var_cov_hat[[i]]) = covariates
    colnames(var_cov_hat[[i]]) = covariates
    var_hat[i, ] = diag(var_cov_hat[[i]])
  }
  fiuo_var_cov = list(var_cov_hat = var_cov_hat, var_hat = var_hat)
  return(fiuo_var_cov)
}
