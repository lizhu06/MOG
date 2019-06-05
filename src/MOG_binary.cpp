#include <RcppEigen.h>
#include <R.h>
#include <Rmath.h>

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace Rcpp;
using namespace std;
using Eigen::Map;                 
using Eigen::MatrixXd;                  
using Eigen::VectorXd;                  
using Eigen::MatrixXi;                  
using Eigen::VectorXi;  

// [[Rcpp::export]]
void set_seed_mb(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);  
}

// [[Rcpp::export]]
void update_gamma1_mb(double pi1, double sigma2, int C, int P,  
  Eigen::VectorXi& c_index, Eigen::VectorXd& gamma_1, 
  Eigen::VectorXd& gamma_1_rep_per_feature,
  Eigen::VectorXd& gamma_2_rep_per_feature, 
  Eigen::VectorXd& gamma_3, Eigen::VectorXd& b, 
  Eigen::VectorXd& Y, Eigen::MatrixXd& X) {
  
  Eigen::VectorXd tailvec = gamma_2_rep_per_feature.array() * gamma_3.array() * 
    b.array();
  Eigen::ArrayXd beta_old = gamma_1_rep_per_feature.array() * tailvec.array();
  Eigen::VectorXd res_old = Y - X * beta_old.matrix();
  double mul = 1.0 / (2.0 * sigma2);
  double L0_L1;

  for(int c = 0; c < C; c++) {
    int fc = c+1;
    Eigen::VectorXd res_change = res_old;
    if(gamma_1[c] == 0){
      for(int p = 1; p < (P+1); p++){
        if(c_index[p] == fc){
          res_change = res_change - X.col(p) * tailvec[p];
        }
      }
      L0_L1 = exp(mul * (res_change.dot(res_change) - 
        res_old.dot(res_old)));
    }else{
      for(int p = 1; p < (P+1); p++){
        if(c_index[p] == fc){
          res_change = res_change + X.col(p) * tailvec[p];
        }
      }
      L0_L1 = exp(mul * (res_old.dot(res_old) - 
        res_change.dot(res_change)));
    }
    double p_gamma_1 = pi1 / (pi1 + (1.0 - pi1) * L0_L1);
    double new_gamma_1 = R::rbinom(1, p_gamma_1);
    if(new_gamma_1 != gamma_1[c]){
      res_old = res_change;
      gamma_1[c] = new_gamma_1;
      for(int p = 1; p < (P+1); p++){
        if(c_index[p] == fc){
          gamma_1_rep_per_feature[p] = new_gamma_1;
        }
      }  
    }else{
      res_change = res_old;
    }
  }
  //return (gamma_1);
}

// [[Rcpp::export]]
void update_gamma2_mb(Eigen::VectorXd& pi2_rep_per_group, 
  double sigma2, int G, 
  int P, Eigen::VectorXi& g_index, Eigen::VectorXi& c_index_group, 
  Eigen::VectorXd& m_gc, 
  Eigen::VectorXd& gamma_1_rep_per_feature, 
  Eigen::VectorXd& gamma_2, Eigen::VectorXd& gamma_2_rep_per_feature, 
  Eigen::VectorXd& gamma_3, Eigen::VectorXd& b, Eigen::VectorXd& Y, 
  Eigen::MatrixXd& X, int weight_Dk_int, Eigen::VectorXd& weight_for_Dk) {

  Eigen::VectorXd tailvec = gamma_1_rep_per_feature.array() * gamma_3.array() * 
    b.array();
  Eigen::ArrayXd beta_old = gamma_2_rep_per_feature.array() * tailvec.array();
  Eigen::VectorXd res_old = Y - X * beta_old.matrix();

  double mul = 1.0 / (2 * sigma2);
  double L0_L1;
  for(int g = 0; g < G; g++) {
    int fg = g + 1;
    Eigen::VectorXd res_change = res_old;
    if(gamma_2[g] == 0){
      for(int p = 1; p < (P+1); p++){
        if(g_index[p] == fg){
          res_change = res_change - X.col(p) * tailvec[p];
        }
      }
      L0_L1 = exp(mul * (res_change.dot(res_change) - 
        res_old.dot(res_old)));
    }else{
      for(int p = 1; p < (P+1); p++){
        if(g_index[p] == fg){
          res_change = res_change + X.col(p) * tailvec[p];
        }
      }
      L0_L1 = exp(mul * (res_old.dot(res_old) - 
        res_change.dot(res_change)));
    }
    double p_gamma_2 = 0.0; 
    if(weight_Dk_int == 1){
      double pi2_weighted = pi2_rep_per_group[g]*weight_for_Dk[g];
      p_gamma_2 = pi2_weighted / (pi2_weighted 
        + (1 - pi2_weighted) * L0_L1);
    }else{
      p_gamma_2 = pi2_rep_per_group[g] / (pi2_rep_per_group[g] 
        + (1 - pi2_rep_per_group[g]) * L0_L1);
    }

    double new_gamma_2 = R::rbinom(1, p_gamma_2);
    if(new_gamma_2 != gamma_2[g]){
      res_old = res_change;
      gamma_2[g] = new_gamma_2;
      for(int p = 1; p < (1+P); p++){
        if(g_index[p] == fg){
          gamma_2_rep_per_feature[p] = new_gamma_2;
        }
      }  
    }else{
      res_change = res_old;
    }
  }
  //return (gamma_2);
}

// [[Rcpp::export]]
void update_gamma3_mb(Eigen::VectorXd& pi3_rep_per_feature, double sigma2, int P, 
  Eigen::VectorXi& c_index, Eigen::VectorXi& g_index, Eigen::VectorXd& m_fg,
  Eigen::VectorXd& gamma_1_rep_per_feature, 
  Eigen::VectorXd& gamma_2_rep_per_feature, Eigen::VectorXd& gamma_3, 
  Eigen::VectorXd& b, Eigen::VectorXd& beta, Eigen::VectorXd& Y, Eigen::MatrixXd& X,
  int weight_Tj_int, Eigen::VectorXd& weight_for_Tj, 
  int weight_Rj_int, Eigen::VectorXd& weight_for_Rj) {

  Eigen::VectorXd tailvec = gamma_1_rep_per_feature.array() * 
    gamma_2_rep_per_feature.array() * b.array();
  
  Eigen::VectorXd res_old = Y - X * beta.matrix();
  double mul = 1.0 / (2.0 * sigma2);
  double L0_L1;

  for(int p = 1; p < (P+1); p++) {
    Eigen::VectorXd res_change = res_old;
    if(gamma_3[p] == 0){
      res_change = res_change - X.col(p) * tailvec[p];
      L0_L1 = exp(mul * (res_change.dot(res_change) - 
        res_old.dot(res_old)));
    }else{
      res_change = res_change + X.col(p) * tailvec[p];
      L0_L1 = exp(mul * (res_old.dot(res_old) - 
        res_change.dot(res_change)));
    }

    double p_gamma_3 = 0.0; 
    if(weight_Tj_int == 1){
      double pi3_weighted = pi3_rep_per_feature[p]*weight_for_Tj[p];
      p_gamma_3 = pi3_weighted / (pi3_weighted 
        + (1 - pi3_weighted) * L0_L1);
    }else if(weight_Rj_int == 1){
      double pi3_weighted = pi3_rep_per_feature[p]*weight_for_Rj[p];
      p_gamma_3 = pi3_weighted / (pi3_weighted 
        + (1 - pi3_weighted) * L0_L1);
    }else{
      p_gamma_3 = pi3_rep_per_feature[p] / (pi3_rep_per_feature[p] 
        + (1 - pi3_rep_per_feature[p]) * L0_L1);
    }
    double new_gamma_3 = R::rbinom(1, p_gamma_3);

    if(new_gamma_3 != gamma_3[p]){
      res_old = res_change;
      gamma_3[p] = new_gamma_3;
      beta[p] = gamma_1_rep_per_feature[p] * 
        gamma_2_rep_per_feature[p] * gamma_3[p] * b[p];
    }else{
      res_change = res_old;
    }
  }
  //return (gamma_3);
}

// [[Rcpp::export]]
void update_b_mb(Eigen::VectorXd& pi3_rep_per_feature, 
  double sigma2, 
  Eigen::VectorXd& s2_vec, int P, Eigen::VectorXi& c_index, 
  Eigen::VectorXd& gamma_1_rep_per_feature, 
  Eigen::VectorXd& gamma_2_rep_per_feature, 
  Eigen::VectorXd& gamma_3, 
  Eigen::VectorXd& b, Eigen::VectorXd& beta, Eigen::VectorXd& Y, Eigen::MatrixXd& X,
  int weight_s2_int, Eigen::VectorXd& weight_for_s2) {

  Eigen::VectorXd res = Y - X * beta.matrix();

  for(int p = 0; p < (P+1); p++) {
    if(gamma_1_rep_per_feature[p] * gamma_2_rep_per_feature[p] * gamma_3[p] == 0){
      if(weight_s2_int == 1){
        b[p] = R::rnorm(0.0, sqrt(s2_vec[p] * weight_for_s2[p]));
      }else{
        b[p] = R::rnorm(0.0, sqrt(s2_vec[p]));
      }
    }else{
      double sigma2_b = 0.0;
      double mu_b = 0.0;
      if(weight_s2_int == 1){
        sigma2_b = 1.0 / ((X.col(p).dot(X.col(p))) / sigma2 + 
          1.0 / (weight_for_s2[p]*s2_vec[p]));
        mu_b = X.col(p).dot((res + X.col(p) * beta[p])) / sigma2 * sigma2_b;
      }else{
        sigma2_b = 1.0 / ((X.col(p).dot(X.col(p))) / sigma2 + 
          1.0 / s2_vec[p]);
        mu_b = X.col(p).dot((res + X.col(p) * beta[p])) / sigma2 * sigma2_b;
      }
      b[p] = R::rnorm(mu_b, sqrt(sigma2_b));
    }
    double betap_new = gamma_1_rep_per_feature[p] * 
      gamma_2_rep_per_feature[p] * gamma_3[p] * b[p];
    res = res - X.col(p) * (betap_new - beta[p]);
    beta[p] = betap_new;
  }
  //return (b);
}

// [[Rcpp::export]]
double update_pi1_mb(Eigen::VectorXd& gamma_1, 
  double alpha1, double beta1, int C, int G) {
  double sum_c = gamma_1.sum();
  double pi1 = R::rbeta(sum_c + alpha1, C - sum_c + beta1);
  return (pi1);
}

// [[Rcpp::export]]
void update_pi2_mb(Eigen::VectorXd& gamma_2, Eigen::VectorXd& pi2, 
  Eigen::VectorXi& c_index_group, 
  Eigen::VectorXd& m_gc, double alpha2, double beta2, int C, int G,
  int weight_Dk_int ,Eigen::VectorXd& weight_for_Dk, 
  Eigen::VectorXd& pi2_loglikeli, double pi2_prop_n, int MH_ind) {
  //Eigen::VectorXd pi2(C);
  for(int c = 0; c < C; c++){
    double fc = c+1;
    if(MH_ind == 1){
      double loglikeli_new = 0.0;
      double pi2_new = R::rbeta(pi2_prop_n * pi2[c], 
        pi2_prop_n * (1.0 - pi2[c]));
      if(pi2_new > 0.9) {pi2_new = 0.9;}
      if(pi2_new < 0.1) {pi2_new = 0.1;}
      for(int g = 0; g < G; g++){
        if(c_index_group[g] == fc){
          if(gamma_2[g]==1){
            if(weight_Dk_int==1){
              loglikeli_new = loglikeli_new + 
                log(pi2_new * weight_for_Dk[g]);
            }else{
              loglikeli_new = loglikeli_new + 
                log(pi2_new);
            }
          }else{
            if(weight_Dk_int==1){
              loglikeli_new = loglikeli_new + 
                log(1-(pi2_new * weight_for_Dk[g]));
            }else{
              loglikeli_new = loglikeli_new + log(1-(pi2_new));
            }
          }
        }
      }
      double acc_p = min(0.0, loglikeli_new - pi2_loglikeli[c] + 
        R::dbeta(pi2[c], pi2_prop_n*pi2_new, pi2_prop_n*(1-pi2_new), 1) -
        R::dbeta(pi2_new, pi2_prop_n*pi2[c], pi2_prop_n*(1-pi2[c]), 1));
      double u_uni = R::runif(0, 1);
      if(log(u_uni) < acc_p){
        pi2[c] = pi2_new;
        pi2_loglikeli[c] = loglikeli_new;
      }
    }else{
      double sum_g = 0;
      for(int g = 0; g < G; g++){
        if(c_index_group[g] == fc){
          sum_g = sum_g + gamma_2[g];
        }
      }
      pi2[c] = R::rbeta(sum_g + alpha2, m_gc[c] - sum_g + beta2);
    }
  }
  //return (pi2);
}

// [[Rcpp::export]]
void update_pi3_mb(Eigen::VectorXd& gamma_3, Eigen::VectorXd& pi3, 
  Eigen::VectorXi& g_index, 
  Eigen::VectorXd& m_fg, double alpha3, double beta3, int G, int P,
  int weight_Tj_int ,Eigen::VectorXd& weight_for_Tj, 
  int weight_Rj_int ,Eigen::VectorXd& weight_for_Rj, 
  Eigen::VectorXd& pi3_loglikeli, double pi3_prop_n, int MH_ind) {
  //Eigen::VectorXd pi3(G);
  for(int g = 0; g < G; g++){
    double fg = g + 1;
    if(MH_ind == 1){
      double loglikeli_new = 0.0;
      double pi3_new = R::rbeta(pi3_prop_n * pi3[g], 
        pi3_prop_n * (1.0 - pi3[g]));
      if(pi3_new > 0.999999) {pi3_new = 0.999999;}
      if(pi3_new < 0.000001) {pi3_new = 0.000001;}
      for(int p = 0; p < P; p++){
        if(g_index[p] == fg){
          if(gamma_3[p]==1){
            if(weight_Tj_int==1){
              loglikeli_new = loglikeli_new + 
                log(pi3_new * weight_for_Tj[p]);
            }else if(weight_Rj_int==1){
              loglikeli_new = loglikeli_new + 
                log(pi3_new * weight_for_Rj[p]);
            }else{
              loglikeli_new = loglikeli_new + 
                log(pi3_new);
            }
          }else{
            if(weight_Tj_int==1){
              loglikeli_new = loglikeli_new + 
                log(1-(pi3_new * weight_for_Tj[p]));
            }else if(weight_Rj_int==1){
              loglikeli_new = loglikeli_new + 
                log(1-(pi3_new * weight_for_Rj[p]));
            }else{
              loglikeli_new = loglikeli_new + log(1-(pi3_new));
            }
          }
        }
      }
      double acc_p = min(0.0, loglikeli_new - pi3_loglikeli[g] + 
        R::dbeta(pi3[g], pi3_prop_n*pi3_new, pi3_prop_n*(1-pi3_new), 1) -
        R::dbeta(pi3_new, pi3_prop_n*pi3[g], pi3_prop_n*(1-pi3[g]), 1));
      double u_uni = R::runif(0, 1);
      if(log(u_uni) < acc_p){
        pi3[g] = pi3_new;
        pi3_loglikeli[g] = loglikeli_new;
      }
    }else{
      double sum_p = 0.0;
      for(int p = 0; p < P; p++){
        if(g_index[p] == fg){
          sum_p = sum_p + gamma_3[p];
        }
      }
      pi3[g] = R::rbeta(sum_p + alpha3, m_fg[g] - sum_p + beta3);
    }
  }
  //return (pi3);
}

// [[Rcpp::export]]
double update_sigma2_mb(double N, double alpha_sigma, double beta_sigma,
  Eigen::VectorXd& gamma_1_rep_per_feature, Eigen::VectorXd& gamma_2_rep_per_feature,
  Eigen::VectorXd& gamma_3, Eigen::VectorXd& b, Eigen::VectorXd& Y, Eigen::MatrixXd& X) {
    
    Eigen::VectorXd beta = gamma_1_rep_per_feature.array() * 
      gamma_2_rep_per_feature.array() * gamma_3.array() * 
      b.array();
    Eigen::VectorXd res = Y - X * beta.matrix();
    double arg1 = (N*1.0) / 2.0 + alpha_sigma;
    double arg2 = res.dot(res) / 2.0 + beta_sigma;
    double sigma2 = 1.0 / R::rgamma(arg1, 1.0 /arg2);
    return (sigma2);
}

// [[Rcpp::export]]
void update_s2_mb(int P, double alpha_s, double beta_s,
  Eigen::VectorXi& uni_types, Eigen::VectorXi& types, Eigen::VectorXd& b,
  int weight_s2_int, Eigen::VectorXd& weight_for_s2,
  Eigen::VectorXd& s2) {

  int T = uni_types.size(); 
  //Eigen::VectorXd s2(T);
  //s2[T-1] = 0.0001;
  for(int t = 0; t < (T-1); t++){
    double sum_t = 0.0;
    double sum_b2 = 0.0;
    for(int p = 1; p < (P+1); p++){
      if(types[p] == uni_types[t]){
        sum_t = sum_t + 1;
        if(weight_s2_int==1){
          sum_b2 = sum_b2 + pow(b[p], 2)/weight_for_s2[p];
        }else{
          sum_b2 = sum_b2 + pow(b[p], 2);
        }
      }
    }
    double arg1 = sum_t / 2.0 + alpha_s;
    double arg2 = sum_b2 / 2.0 + beta_s;
    s2[t] = 1.0 / R::rgamma(arg1, 1.0/arg2);
  }
  //return (s2);
}

// [[Rcpp::export]]
Rcpp::List MCMC_mb(int seed, int burnInIter, int keepIter, 
  int C, int G, int P, double N, 
  double alpha1, double alpha2, double alpha3, double alpha_sigma,
  double alpha_s, double beta1, double beta2, double beta3,
  double beta_sigma, double beta_s,
  double pi1, double sigma2, 
  Eigen::VectorXi& uni_types, Eigen::VectorXi& types,
  Eigen::VectorXd& m_gc, Eigen::VectorXd& m_fg,
  Eigen::VectorXi& c_index, Eigen::VectorXi& c_index_group, Eigen::VectorXi& g_index, 
  Eigen::VectorXd& gamma_1, Eigen::VectorXd& gamma_1_rep_per_feature, Eigen::VectorXd& gamma_2,
  Eigen::VectorXd& gamma_2_rep_per_feature, Eigen::VectorXd& gamma_3, Eigen::VectorXd& b, 
  Eigen::VectorXd& pi2, Eigen::VectorXd& pi2_rep_per_group, Eigen::VectorXd& pi3, 
  Eigen::VectorXd& pi3_rep_per_feature,
  Eigen::VectorXd& s2,
  Eigen::VectorXd& Yb, Eigen::MatrixXd& X, int printInt, 
  Eigen::VectorXd& weight_for_s2, Eigen::VectorXd& weight_for_Rj, 
  Eigen::VectorXd& weight_for_Tj, Eigen::VectorXd& weight_for_Dk,
  int weight_s2_int, int weight_Rj_int, int weight_Tj_int, 
  int weight_Dk_int, 
  Eigen::VectorXd& pi2_loglikeli, Eigen::VectorXd& pi3_loglikeli,
  double pi2_prop_n, double pi3_prop_n, int MH_ind){

  int num_t = uni_types.size();

  Eigen::MatrixXd GAMMA_1(C, keepIter);
  Eigen::MatrixXd GAMMA_2(G, keepIter);
  Eigen::MatrixXd GAMMA_3(P+1, keepIter);
  Eigen::MatrixXd B(P+1, keepIter);
  Eigen::MatrixXd PI2(C, keepIter);
  Eigen::MatrixXd PI3(G, keepIter);
  Eigen::MatrixXd S2(num_t, keepIter);
  Eigen::VectorXd SIGMA2(keepIter);
  Eigen::VectorXd PI1(keepIter);
  
  set_seed_mb(seed);
  int totalIter = burnInIter + keepIter;
  Eigen::VectorXd s2_vec(P+1);
  for(int p = 0; p < (P+1); p++){
    int index = types[p] - 1;
    s2_vec[p] = s2[index];
  }
  Eigen::VectorXd beta = gamma_1_rep_per_feature.array() * 
    gamma_2_rep_per_feature.array() * gamma_3.array() * 
    b.array();

  for(int i = 0; i < totalIter; i++){
    Eigen::VectorXd mu = X * beta.matrix();
    Eigen::VectorXd Y = Yb; 
    for(int n = 0; n < N; n++){
      double u = R::runif(0, 1);
      if(Yb[n] == 1){
        double phi_a = R::pnorm(0.0 - mu[n], 0, 1, 1, 0);
        if(phi_a > 0.9){
          phi_a = 0.9;
        }
        double phi_b = 1.0;
        Y[n] = R::qnorm(phi_a + (phi_b - phi_a) * u, 0, 1, 1, 0) + mu[n];
      }else{
        double phi_a = 0.0;
        double phi_b = R::pnorm(0.0 - mu[n], 0, 1, 1, 0);
        if(phi_b < 0.1){
          phi_b = 0.1;
        }
        Y[n] = R::qnorm(phi_a + (phi_b - phi_a) * u, 0, 1, 1, 0) + mu[n];
      }
    }
    
    update_gamma1_mb(pi1, sigma2, C, P, c_index, 
      gamma_1, gamma_1_rep_per_feature, gamma_2_rep_per_feature, 
      gamma_3, b, Y, X);
    for(int p = 1; p < (P+1); p++){
      int index = c_index[p] - 1;
      gamma_1_rep_per_feature[p] = gamma_1[index];
    }
  
    update_gamma2_mb(pi2_rep_per_group, sigma2, G, 
      P, g_index, c_index_group, m_gc, gamma_1_rep_per_feature, gamma_2, 
      gamma_2_rep_per_feature, gamma_3, b, Y, X, 
      weight_Dk_int, weight_for_Dk);
    for(int p = 1; p < (P+1); p++){
      int index = g_index[p] - 1;
      gamma_2_rep_per_feature[p] = gamma_2[index];
    }
   
    update_gamma3_mb(pi3_rep_per_feature, sigma2, P, 
      c_index, g_index, m_fg, gamma_1_rep_per_feature, gamma_2_rep_per_feature, 
      gamma_3, b, beta, Y, X, weight_Tj_int, weight_for_Tj, 
      weight_Rj_int, weight_for_Rj);
    

    update_b_mb(pi3_rep_per_feature, sigma2, s2_vec, P, 
      c_index, gamma_1_rep_per_feature,
       gamma_2_rep_per_feature, 
      gamma_3, b, beta, Y, X, weight_s2_int, weight_for_s2);
    

    pi1 = update_pi1_mb(gamma_1, alpha1, beta1, C, G);
    

    update_pi2_mb(gamma_2, pi2, c_index_group, m_gc, alpha2, beta2, C, G,
      weight_Dk_int, weight_for_Dk, pi2_loglikeli, pi2_prop_n, MH_ind);
    
    for(int g = 0; g < G; g++){
      int index = c_index_group[g] - 1;
      pi2_rep_per_group[g] = pi2[index];
    }

    update_pi3_mb(gamma_3, pi3, g_index, m_fg, alpha3, beta3, G, P,
      weight_Tj_int, weight_for_Tj, weight_Rj_int, weight_for_Rj, 
      pi3_loglikeli,pi3_prop_n, MH_ind);
    
    for(int p = 0; p < P; p++){
      int index = g_index[p] - 1;
      pi3_rep_per_feature[p] = pi3[index];
    }

    update_s2_mb(P, alpha_s, beta_s, uni_types, types, b, weight_s2_int, 
      weight_for_s2, s2);
    
    for(int p = 0; p < (P+1); p++){
      int index = types[p] - 1;
      s2_vec[p] = s2[index];
    }
    if(i > burnInIter - 1){
      int store_index = i - burnInIter;
      GAMMA_1.col(store_index) = gamma_1;
      GAMMA_2.col(store_index) = gamma_2;
      GAMMA_3.col(store_index) = gamma_3;
      B.col(store_index) = b;
      PI1[store_index] = pi1;
      PI2.col(store_index) = pi2;
      PI3.col(store_index) = pi3;
      SIGMA2[store_index] = sigma2;
      S2.col(store_index) = s2;
    }

    if(i%printInt == 0){
      cout<< "MCMC iteration: " << i << ".\n";
    }
  }

  return Rcpp::List::create(Rcpp::Named("GAMMA_1") = GAMMA_1,
                            Rcpp::Named("GAMMA_2") = GAMMA_2,
                            Rcpp::Named("GAMMA_3") = GAMMA_3,
                            Rcpp::Named("B") = B,
                            Rcpp::Named("PI1") = PI1,
                            Rcpp::Named("PI2") = PI2,
                            Rcpp::Named("PI3") = PI3,
                            Rcpp::Named("S2") = S2
                            ); 
}

































