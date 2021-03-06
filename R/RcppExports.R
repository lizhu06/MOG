# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

set_seed_mb <- function(seed) {
    invisible(.Call('_MOG_set_seed_mb', PACKAGE = 'MOG', seed))
}

update_gamma1_mb <- function(pi1, sigma2, C, P, c_index, gamma_1, gamma_1_rep_per_feature, gamma_2_rep_per_feature, gamma_3, b, Y, X) {
    invisible(.Call('_MOG_update_gamma1_mb', PACKAGE = 'MOG', pi1, sigma2, C, P, c_index, gamma_1, gamma_1_rep_per_feature, gamma_2_rep_per_feature, gamma_3, b, Y, X))
}

update_gamma2_mb <- function(pi2_rep_per_group, sigma2, G, P, g_index, c_index_group, m_gc, gamma_1_rep_per_feature, gamma_2, gamma_2_rep_per_feature, gamma_3, b, Y, X, weight_Dk_int, weight_for_Dk) {
    invisible(.Call('_MOG_update_gamma2_mb', PACKAGE = 'MOG', pi2_rep_per_group, sigma2, G, P, g_index, c_index_group, m_gc, gamma_1_rep_per_feature, gamma_2, gamma_2_rep_per_feature, gamma_3, b, Y, X, weight_Dk_int, weight_for_Dk))
}

update_gamma3_mb <- function(pi3_rep_per_feature, sigma2, P, c_index, g_index, m_fg, gamma_1_rep_per_feature, gamma_2_rep_per_feature, gamma_3, b, beta, Y, X, weight_Tj_int, weight_for_Tj, weight_Rj_int, weight_for_Rj) {
    invisible(.Call('_MOG_update_gamma3_mb', PACKAGE = 'MOG', pi3_rep_per_feature, sigma2, P, c_index, g_index, m_fg, gamma_1_rep_per_feature, gamma_2_rep_per_feature, gamma_3, b, beta, Y, X, weight_Tj_int, weight_for_Tj, weight_Rj_int, weight_for_Rj))
}

update_b_mb <- function(pi3_rep_per_feature, sigma2, s2_vec, P, c_index, gamma_1_rep_per_feature, gamma_2_rep_per_feature, gamma_3, b, beta, Y, X, weight_s2_int, weight_for_s2) {
    invisible(.Call('_MOG_update_b_mb', PACKAGE = 'MOG', pi3_rep_per_feature, sigma2, s2_vec, P, c_index, gamma_1_rep_per_feature, gamma_2_rep_per_feature, gamma_3, b, beta, Y, X, weight_s2_int, weight_for_s2))
}

update_pi1_mb <- function(gamma_1, alpha1, beta1, C, G) {
    .Call('_MOG_update_pi1_mb', PACKAGE = 'MOG', gamma_1, alpha1, beta1, C, G)
}

update_pi2_mb <- function(gamma_2, pi2, c_index_group, m_gc, alpha2, beta2, C, G, weight_Dk_int, weight_for_Dk, pi2_loglikeli, pi2_prop_n, MH_ind) {
    invisible(.Call('_MOG_update_pi2_mb', PACKAGE = 'MOG', gamma_2, pi2, c_index_group, m_gc, alpha2, beta2, C, G, weight_Dk_int, weight_for_Dk, pi2_loglikeli, pi2_prop_n, MH_ind))
}

update_pi3_mb <- function(gamma_3, pi3, g_index, m_fg, alpha3, beta3, G, P, weight_Tj_int, weight_for_Tj, weight_Rj_int, weight_for_Rj, pi3_loglikeli, pi3_prop_n, MH_ind) {
    invisible(.Call('_MOG_update_pi3_mb', PACKAGE = 'MOG', gamma_3, pi3, g_index, m_fg, alpha3, beta3, G, P, weight_Tj_int, weight_for_Tj, weight_Rj_int, weight_for_Rj, pi3_loglikeli, pi3_prop_n, MH_ind))
}

update_sigma2_mb <- function(N, alpha_sigma, beta_sigma, gamma_1_rep_per_feature, gamma_2_rep_per_feature, gamma_3, b, Y, X) {
    .Call('_MOG_update_sigma2_mb', PACKAGE = 'MOG', N, alpha_sigma, beta_sigma, gamma_1_rep_per_feature, gamma_2_rep_per_feature, gamma_3, b, Y, X)
}

update_s2_mb <- function(P, alpha_s, beta_s, uni_types, types, b, weight_s2_int, weight_for_s2, s2) {
    invisible(.Call('_MOG_update_s2_mb', PACKAGE = 'MOG', P, alpha_s, beta_s, uni_types, types, b, weight_s2_int, weight_for_s2, s2))
}

MCMC_mb <- function(seed, burnInIter, keepIter, C, G, P, N, alpha1, alpha2, alpha3, alpha_sigma, alpha_s, beta1, beta2, beta3, beta_sigma, beta_s, pi1, sigma2, uni_types, types, m_gc, m_fg, c_index, c_index_group, g_index, gamma_1, gamma_1_rep_per_feature, gamma_2, gamma_2_rep_per_feature, gamma_3, b, pi2, pi2_rep_per_group, pi3, pi3_rep_per_feature, s2, Yb, X, printInt, weight_for_s2, weight_for_Rj, weight_for_Tj, weight_for_Dk, weight_s2_int, weight_Rj_int, weight_Tj_int, weight_Dk_int, pi2_loglikeli, pi3_loglikeli, pi2_prop_n, pi3_prop_n, MH_ind) {
    .Call('_MOG_MCMC_mb', PACKAGE = 'MOG', seed, burnInIter, keepIter, C, G, P, N, alpha1, alpha2, alpha3, alpha_sigma, alpha_s, beta1, beta2, beta3, beta_sigma, beta_s, pi1, sigma2, uni_types, types, m_gc, m_fg, c_index, c_index_group, g_index, gamma_1, gamma_1_rep_per_feature, gamma_2, gamma_2_rep_per_feature, gamma_3, b, pi2, pi2_rep_per_group, pi3, pi3_rep_per_feature, s2, Yb, X, printInt, weight_for_s2, weight_for_Rj, weight_for_Tj, weight_for_Dk, weight_s2_int, weight_Rj_int, weight_Tj_int, weight_Dk_int, pi2_loglikeli, pi3_loglikeli, pi2_prop_n, pi3_prop_n, MH_ind)
}

set_seed_mc <- function(seed) {
    invisible(.Call('_MOG_set_seed_mc', PACKAGE = 'MOG', seed))
}

update_gamma1_mc <- function(pi1, sigma2, C, P, c_index, gamma_1, gamma_1_rep_per_feature, gamma_2_rep_per_feature, gamma_3, b, Y, X) {
    invisible(.Call('_MOG_update_gamma1_mc', PACKAGE = 'MOG', pi1, sigma2, C, P, c_index, gamma_1, gamma_1_rep_per_feature, gamma_2_rep_per_feature, gamma_3, b, Y, X))
}

update_gamma2_mc <- function(pi2_rep_per_group, sigma2, G, P, g_index, c_index_group, m_gc, gamma_1_rep_per_feature, gamma_2, gamma_2_rep_per_feature, gamma_3, b, Y, X, weight_Dk_int, weight_for_Dk) {
    invisible(.Call('_MOG_update_gamma2_mc', PACKAGE = 'MOG', pi2_rep_per_group, sigma2, G, P, g_index, c_index_group, m_gc, gamma_1_rep_per_feature, gamma_2, gamma_2_rep_per_feature, gamma_3, b, Y, X, weight_Dk_int, weight_for_Dk))
}

update_gamma3_mc <- function(pi3_rep_per_feature, sigma2, P, c_index, g_index, m_fg, gamma_1_rep_per_feature, gamma_2_rep_per_feature, gamma_3, b, beta, Y, X, weight_Tj_int, weight_for_Tj, weight_Rj_int, weight_for_Rj) {
    invisible(.Call('_MOG_update_gamma3_mc', PACKAGE = 'MOG', pi3_rep_per_feature, sigma2, P, c_index, g_index, m_fg, gamma_1_rep_per_feature, gamma_2_rep_per_feature, gamma_3, b, beta, Y, X, weight_Tj_int, weight_for_Tj, weight_Rj_int, weight_for_Rj))
}

update_b_mc <- function(pi3_rep_per_feature, sigma2, s2_vec, P, c_index, gamma_1_rep_per_feature, gamma_2_rep_per_feature, gamma_3, b, beta, Y, X, weight_s2_int, weight_for_s2) {
    invisible(.Call('_MOG_update_b_mc', PACKAGE = 'MOG', pi3_rep_per_feature, sigma2, s2_vec, P, c_index, gamma_1_rep_per_feature, gamma_2_rep_per_feature, gamma_3, b, beta, Y, X, weight_s2_int, weight_for_s2))
}

update_pi1_mc <- function(gamma_1, alpha1, beta1, C, G) {
    .Call('_MOG_update_pi1_mc', PACKAGE = 'MOG', gamma_1, alpha1, beta1, C, G)
}

update_pi2_mc <- function(gamma_2, pi2, c_index_group, m_gc, alpha2, beta2, C, G, weight_Dk_int, weight_for_Dk, pi2_loglikeli, pi2_prop_n, MH_ind) {
    invisible(.Call('_MOG_update_pi2_mc', PACKAGE = 'MOG', gamma_2, pi2, c_index_group, m_gc, alpha2, beta2, C, G, weight_Dk_int, weight_for_Dk, pi2_loglikeli, pi2_prop_n, MH_ind))
}

update_pi3_mc <- function(gamma_3, pi3, g_index, m_fg, alpha3, beta3, G, P, weight_Tj_int, weight_for_Tj, weight_Rj_int, weight_for_Rj, pi3_loglikeli, pi3_prop_n, MH_ind) {
    invisible(.Call('_MOG_update_pi3_mc', PACKAGE = 'MOG', gamma_3, pi3, g_index, m_fg, alpha3, beta3, G, P, weight_Tj_int, weight_for_Tj, weight_Rj_int, weight_for_Rj, pi3_loglikeli, pi3_prop_n, MH_ind))
}

update_sigma2_mc <- function(N, alpha_sigma, beta_sigma, gamma_1_rep_per_feature, gamma_2_rep_per_feature, gamma_3, beta, Y, X) {
    .Call('_MOG_update_sigma2_mc', PACKAGE = 'MOG', N, alpha_sigma, beta_sigma, gamma_1_rep_per_feature, gamma_2_rep_per_feature, gamma_3, beta, Y, X)
}

update_s2_mc <- function(P, alpha_s, beta_s, uni_types, types, b, weight_s2_int, weight_for_s2, s2) {
    invisible(.Call('_MOG_update_s2_mc', PACKAGE = 'MOG', P, alpha_s, beta_s, uni_types, types, b, weight_s2_int, weight_for_s2, s2))
}

MCMC_mc <- function(seed, burnInIter, keepIter, C, G, P, N, alpha1, alpha2, alpha3, alpha_sigma, alpha_s, beta1, beta2, beta3, beta_sigma, beta_s, pi1, sigma2, uni_types, types, m_gc, m_fg, c_index, c_index_group, g_index, gamma_1, gamma_1_rep_per_feature, gamma_2, gamma_2_rep_per_feature, gamma_3, b, pi2, pi2_rep_per_group, pi3, pi3_rep_per_feature, s2, Y, X, print_int, weight_for_s2, weight_for_Rj, weight_for_Tj, weight_for_Dk, weight_s2_int, weight_Rj_int, weight_Tj_int, weight_Dk_int, pi2_loglikeli, pi3_loglikeli, pi2_prop_n, pi3_prop_n, MH_ind) {
    .Call('_MOG_MCMC_mc', PACKAGE = 'MOG', seed, burnInIter, keepIter, C, G, P, N, alpha1, alpha2, alpha3, alpha_sigma, alpha_s, beta1, beta2, beta3, beta_sigma, beta_s, pi1, sigma2, uni_types, types, m_gc, m_fg, c_index, c_index_group, g_index, gamma_1, gamma_1_rep_per_feature, gamma_2, gamma_2_rep_per_feature, gamma_3, b, pi2, pi2_rep_per_group, pi3, pi3_rep_per_feature, s2, Y, X, print_int, weight_for_s2, weight_for_Rj, weight_for_Tj, weight_for_Dk, weight_s2_int, weight_Rj_int, weight_Tj_int, weight_Dk_int, pi2_loglikeli, pi3_loglikeli, pi2_prop_n, pi3_prop_n, MH_ind)
}

set_seed_sb <- function(seed) {
    invisible(.Call('_MOG_set_seed_sb', PACKAGE = 'MOG', seed))
}

update_gamma1_sb <- function(pi1, sigma2, G, P, g_index, gamma_1, gamma_1_rep_per_feature, gamma_2, b, Y, X) {
    invisible(.Call('_MOG_update_gamma1_sb', PACKAGE = 'MOG', pi1, sigma2, G, P, g_index, gamma_1, gamma_1_rep_per_feature, gamma_2, b, Y, X))
}

update_gamma2_sb <- function(pi2_rep_per_feature, sigma2, P, g_index, m_fg, gamma_1_rep_per_feature, gamma_2, b, beta, Y, X, v_weight, pi2_weighted_int) {
    invisible(.Call('_MOG_update_gamma2_sb', PACKAGE = 'MOG', pi2_rep_per_feature, sigma2, P, g_index, m_fg, gamma_1_rep_per_feature, gamma_2, b, beta, Y, X, v_weight, pi2_weighted_int))
}

update_b_sb <- function(pi2_rep_per_feature, sigma2, s2_vec, P, gamma_1_rep_per_feature, gamma_2, b, Y, X, w_weight) {
    invisible(.Call('_MOG_update_b_sb', PACKAGE = 'MOG', pi2_rep_per_feature, sigma2, s2_vec, P, gamma_1_rep_per_feature, gamma_2, b, Y, X, w_weight))
}

update_pi1_sb <- function(gamma_1, alpha1, beta1, G) {
    .Call('_MOG_update_pi1_sb', PACKAGE = 'MOG', gamma_1, alpha1, beta1, G)
}

update_pi2_sb <- function(gamma_2, pi2, g_index, m_fg, alpha2, beta2, G, P, v_weight, pi2_loglikeli, pi2_prop_n, MH_ind) {
    invisible(.Call('_MOG_update_pi2_sb', PACKAGE = 'MOG', gamma_2, pi2, g_index, m_fg, alpha2, beta2, G, P, v_weight, pi2_loglikeli, pi2_prop_n, MH_ind))
}

update_sigma2_sb <- function(N, alpha_sigma, beta_sigma, gamma_1_rep_per_feature, gamma_2, b, Y, X) {
    .Call('_MOG_update_sigma2_sb', PACKAGE = 'MOG', N, alpha_sigma, beta_sigma, gamma_1_rep_per_feature, gamma_2, b, Y, X)
}

update_s2_sb <- function(P, alpha_s, beta_s, uni_types, types, b, s2, w_weight) {
    invisible(.Call('_MOG_update_s2_sb', PACKAGE = 'MOG', P, alpha_s, beta_s, uni_types, types, b, s2, w_weight))
}

MCMC_sb <- function(seed, burnInIter, keepIter, G, P, printInt, N, alpha1, alpha2, alpha_sigma, alpha_s, beta1, beta2, beta_sigma, beta_s, pi1, sigma2, uni_types, types, m_fg, g_index, gamma_1, gamma_1_rep_per_feature, gamma_2, b, pi2, pi2_rep_per_feature, s2, Yb, X, w_weight, v_weight, pi2_weighted_int, pi2_loglikeli, pi2_prop_n, MH_ind) {
    .Call('_MOG_MCMC_sb', PACKAGE = 'MOG', seed, burnInIter, keepIter, G, P, printInt, N, alpha1, alpha2, alpha_sigma, alpha_s, beta1, beta2, beta_sigma, beta_s, pi1, sigma2, uni_types, types, m_fg, g_index, gamma_1, gamma_1_rep_per_feature, gamma_2, b, pi2, pi2_rep_per_feature, s2, Yb, X, w_weight, v_weight, pi2_weighted_int, pi2_loglikeli, pi2_prop_n, MH_ind)
}

set_seed_sc <- function(seed) {
    invisible(.Call('_MOG_set_seed_sc', PACKAGE = 'MOG', seed))
}

update_gamma1_sc <- function(pi1, sigma2, G, P, g_index, gamma_1, gamma_1_rep_per_feature, gamma_2, b, Y, X) {
    invisible(.Call('_MOG_update_gamma1_sc', PACKAGE = 'MOG', pi1, sigma2, G, P, g_index, gamma_1, gamma_1_rep_per_feature, gamma_2, b, Y, X))
}

update_gamma2_sc <- function(pi2_rep_per_feature, sigma2, P, g_index, m_fg, gamma_1_rep_per_feature, gamma_2, b, beta, Y, X, v_weight, BernoulliWeighted_int) {
    invisible(.Call('_MOG_update_gamma2_sc', PACKAGE = 'MOG', pi2_rep_per_feature, sigma2, P, g_index, m_fg, gamma_1_rep_per_feature, gamma_2, b, beta, Y, X, v_weight, BernoulliWeighted_int))
}

update_b_sc <- function(pi2_rep_per_feature, sigma2, s2_vec, P, gamma_1_rep_per_feature, gamma_2, b, beta, Y, X, w_weight) {
    invisible(.Call('_MOG_update_b_sc', PACKAGE = 'MOG', pi2_rep_per_feature, sigma2, s2_vec, P, gamma_1_rep_per_feature, gamma_2, b, beta, Y, X, w_weight))
}

update_pi1_sc <- function(gamma_1, alpha1, beta1, G) {
    .Call('_MOG_update_pi1_sc', PACKAGE = 'MOG', gamma_1, alpha1, beta1, G)
}

update_pi2_sc <- function(gamma_2, pi2, g_index, m_fg, alpha2, beta2, G, P, v_weight, pi2_loglikeli, pi2_prop_n, MH_ind) {
    invisible(.Call('_MOG_update_pi2_sc', PACKAGE = 'MOG', gamma_2, pi2, g_index, m_fg, alpha2, beta2, G, P, v_weight, pi2_loglikeli, pi2_prop_n, MH_ind))
}

update_sigma2_sc <- function(N, alpha_sigma, beta_sigma, gamma_1_rep_per_feature, gamma_2, beta, Y, X) {
    .Call('_MOG_update_sigma2_sc', PACKAGE = 'MOG', N, alpha_sigma, beta_sigma, gamma_1_rep_per_feature, gamma_2, beta, Y, X)
}

update_s2_sc <- function(P, alpha_s, beta_s, uni_types, types, b, w_weight, s2) {
    invisible(.Call('_MOG_update_s2_sc', PACKAGE = 'MOG', P, alpha_s, beta_s, uni_types, types, b, w_weight, s2))
}

MCMC_sc <- function(seed, burnInIter, keepIter, G, P, printInt, N, alpha1, alpha2, alpha_sigma, alpha_s, beta1, beta2, beta_sigma, beta_s, pi1, sigma2, uni_types, types, m_fg, g_index, gamma_1, gamma_1_rep_per_feature, gamma_2, b, pi2, pi2_rep_per_feature, s2, Y, X, w_weight, v_weight, s2Weighted_int, BernoulliWeighted_int, pi2_loglikeli, pi2_prop_n, MH_ind) {
    .Call('_MOG_MCMC_sc', PACKAGE = 'MOG', seed, burnInIter, keepIter, G, P, printInt, N, alpha1, alpha2, alpha_sigma, alpha_s, beta1, beta2, beta_sigma, beta_s, pi1, sigma2, uni_types, types, m_fg, g_index, gamma_1, gamma_1_rep_per_feature, gamma_2, b, pi2, pi2_rep_per_feature, s2, Y, X, w_weight, v_weight, s2Weighted_int, BernoulliWeighted_int, pi2_loglikeli, pi2_prop_n, MH_ind)
}

rcppeigen_hello_world <- function() {
    .Call('_MOG_rcppeigen_hello_world', PACKAGE = 'MOG')
}

rcppeigen_outerproduct <- function(x) {
    .Call('_MOG_rcppeigen_outerproduct', PACKAGE = 'MOG', x)
}

rcppeigen_innerproduct <- function(x) {
    .Call('_MOG_rcppeigen_innerproduct', PACKAGE = 'MOG', x)
}

rcppeigen_bothproducts <- function(x) {
    .Call('_MOG_rcppeigen_bothproducts', PACKAGE = 'MOG', x)
}

